module mod_initial_data
!=============================================================================
   use mod_prec
   use mod_field,  only: field
   use mod_cheb,   only: Rvec
   use mod_swal,   only: Yvec, swal
   use mod_params, only: &
      dt, nx, ny, &
      min_m, max_m, &
      lin_m, &
      spin=>psi_spin, &
      cl=>compactification_length, &
      bhm=>black_hole_mass, &
      bhs=>black_hole_spin, &
      R_max, &
! initial data parameters 
      initial_data_direction_v=>initial_data_direction, &
      amp_re_v=>amp_re, &
      amp_im_v=>amp_im, &
      rl_v=>rl, &
      ru_v=>ru, &
      l_ang_v=>l_ang
!=============================================================================
   implicit none
   private

   public :: set_initial_data

   complex(rp), parameter :: ZI = (0.0_rp, 1.0_rp)
!=============================================================================
contains
!=============================================================================
   subroutine set_initial_data(m_i,p,q,f)
      integer(ip), intent(in) :: m_i
      type(field), intent(inout) :: p, q, f

      integer(ip) :: i, j
      integer(ip) :: m_ang, l_ang
      real(rp) :: max_val, bump, r, y
      real(rp) :: width, rl, ru 

      complex(rp) :: amp

      character(:), allocatable :: initial_data_direction

      m_ang = lin_m(m_i)
      amp   = amp_re_v(m_i) + ZI*amp_im_v(m_i)
      ru    = ru_v(m_i)
      rl    = rl_v(m_i)
      l_ang = l_ang_v(m_i)
      initial_data_direction = initial_data_direction_v(m_i:m_i)

      max_val = 0.0_rp
      width = maxval([ru-rl,1e-16_rp]) ! need width>0 so no NaNs
!-----------------------------------------------------------------------------
      y_loop: do j=1,ny
      x_loop: do i=1,nx-1 ! index 'i=nx' is at 'r=infinity'
         r = (cl**2) / Rvec(i)

         if ((r<ru).and.(r>rl)) then
            bump = exp(-1.0_rp*width/(r-rl))*exp(-2.0_rp*width/(ru-r))
         else
            bump = 0.0_rp
         end if

         f%n(i,j,m_ang) = (((r-rl)/width)**2) * (((ru-r)/width)**2) * bump
         
         q%n(i,j,m_ang) = ( &
            (2.0_rp*(((r-rl)/width)   )*(((ru-r)/width)**2)) &
         -  (2.0_rp*(((r-rl)/width)**2)*( (ru-r)/width    )) &
         +  (1.0_rp*(1.0_rp           )*(((ru-r)/width)**2)) &
         -  (2.0_rp*(((r-rl)/width)**2)*(1.0_rp              )) &
         )*bump/width
         !--------------------------------------------------------------------
         ! rescale q as q = \partial_R f = -(r/cl)^2 partial_r f
         !--------------------------------------------------------------------
         q%n(i,j,m_ang) = q%n(i,j,m_ang)*(-(r**2)/(cl**2))

         p%n(i,j,m_ang) = 0.0_rp
                             
         f%n(i,j,m_ang) = f%n(i,j,m_ang) * swal(j,l_ang,m_ang,spin)
         q%n(i,j,m_ang) = q%n(i,j,m_ang) * swal(j,l_ang,m_ang,spin)

         max_val = max(abs(f%n(i,j,m_ang)),max_val)

      end do x_loop
      end do y_loop
!-----------------------------------------------------------------------------
! rescale to make max val = 'amp'
!-----------------------------------------------------------------------------
      f%n(:,:,m_ang) = f%n(:,:,m_ang) * (amp / max_val)
      q%n(:,:,m_ang) = q%n(:,:,m_ang) * (amp / max_val)
!-----------------------------------------------------------------------------
! p = (...)\phi_{,t} + (...), so use p to set \phi_{,t}
!-----------------------------------------------------------------------------
   select case (initial_data_direction)
      !-----------------------------------------------------------------------
      ! ingoing
      case ("i")
         do j=1,ny
         do i=1,nx
            R = Rvec(i)
            Y = Yvec(j)

            p%n(i,j,m_ang) = &
               (f%n(i,j,m_ang)*( &
               -  2*bhs**2*R*(cl**2 + 6*bhm*R) &
               +  4*cl**2*bhm*(-(cl**2*spin) + 2*bhm*R*(2 + spin)) &
               +  (0,2)*bhs*(4*cl**2*bhm*m_ang*R + cl**4*(m_ang - spin*Y)) &
               ))/cl**4 &
            +  q%n(i,j,m_ang)*( &
                  (-2*( &
                     cl**6 &
                  +  cl**2*(bhs**2 - 8*bhm**2)*R**2 &
                  +  4*bhs**2*bhm*R**3 &
                  ))/cl**4 &
               -  (R**3*( &
                     8*bhm*(2*bhm - (bhs**2*R)/cl**2)*(1 + (2*bhm*R)/cl**2) &
                  +  bhs**2*(-1 + Y**2) &
                  ))/(cl**2*(2 + (4*bhm*R)/cl**2)) &
               )
         end do
         end do
      !-----------------------------------------------------------------------
      ! outgoing
      case ("o")
         do j=1,ny
         do i=1,nx
            R = Rvec(i)
            Y = Yvec(j)

            p%n(i,j,m_ang) = &
               (q%n(i,j,m_ang)*( &
                  64*cl**4*bhm**4*R**2 &
               +  bhs**4*R**2*(16*bhm**2*R**2 + cl**4*(-1 + Y**2)) &
               +  bhs**2*(-64*cl**2*bhm**3*R**3 + cl**8*(-1 + Y**2) - 2*cl**6*bhm*R*(-1 + Y**2)) &
               ))/( &
                  4.*cl**4*bhm*(2*cl**2*bhm - bhs**2*R) &
               ) &
            +  f%n(i,j,m_ang)*( &
                  (-2*bhs**2*R*(cl**2 + 6*bhm*R))/cl**4 &
               -  4*bhm*spin &
               +  (8*bhm**2*R*(2 + spin))/cl**2 &
               -  (0,2)*bhs*spin*Y &
               -  (0,0.5_rp)*m_ang*( &
                     8 &
                  +  (16*bhm*R)/cl**2 &
                  +  bhs*(-4 - (16*bhm*R)/cl**2) &
                  +  (bhs**2*cl**2*(-1 + Y**2))/(bhm*(2*cl**2*bhm - bhs**2*R)) &
                  ) &
               )
         end do
         end do
      !-----------------------------------------------------------------------
      ! time symmetric 
      case ("t")
         do j=1,ny
         do i=1,nx
            R = Rvec(i)
            Y = Yvec(j)
            
            p%n(i,j,m_ang) = &
               (-2*q%n(i,j,m_ang)*( &
                  cl**6 &
               +  cl**2*(bhs**2 - 8*bhm**2)*R**2 &
               +  4*bhs**2*bhm*R**3 &
               ))/cl**4 &
            +  (f%n(i,j,m_ang)*( &
               -  2*bhs**2*R*(cl**2 + 6*bhm*R) &
               +  4*cl**2*bhm*(-(cl**2*spin) + 2*bhm*R*(2 + spin)) &
               +  (0,2)*bhs*(4*cl**2*bhm*m_ang*R + cl**4*(m_ang - spin*Y)) &
               ))/cl**4 
         end do
         end do
      case default
         ! do nothing
   end select 
!-----------------------------------------------------------------------------
! copy to np1 so can be saved
!-----------------------------------------------------------------------------
      p%np1(:,:,m_ang) = p%n(:,:,m_ang) 
      q%np1(:,:,m_ang) = q%n(:,:,m_ang) 
      f%np1(:,:,m_ang) = f%n(:,:,m_ang) 

   end subroutine set_initial_data
!=============================================================================
end module mod_initial_data

