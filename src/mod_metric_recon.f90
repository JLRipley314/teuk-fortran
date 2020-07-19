!
! Metric reconstruction evolution equations
!
!=============================================================================
module mod_metric_recon
!=============================================================================
   use mod_prec
   use mod_cheb,     only: Rvec=>R, compute_DR
   use mod_field,    only: field, set_field, set_level
   use mod_ghp,      only: set_edth, set_edth_prime, set_thorn, set_thorn_prime
   use mod_bkgrd_np, only: mu_0, ta_0, pi_0, rh_0, thorn_prime_ta_0, psi2_0
   use mod_teuk,     only: psi4_f
   use mod_params,   only: &
      dt, nx, ny, min_m, max_m, &
      cl=>compactification_length, &
      bhm=>black_hole_mass, &
      bhs=>black_hole_spin

!=============================================================================
   implicit none
   private

   public :: metric_recon_time_step, metric_recon_indep_res

   type(field), public :: &
      psi3, psi2,         &
      la, pi,             &
      muhll, hlmb, hmbmb, &
      bianci3_res, bianci2_res, hll_res 
!=============================================================================
   contains
!=============================================================================
   pure subroutine set_k(m_ang, fname, falloff, level, level_DR, kl)
      integer(ip),  intent(in) :: m_ang
      character(*), intent(in) :: fname
      integer(ip),  intent(in) :: falloff
      complex(rp), dimension(nx,ny,min_m:max_m), intent(in)    :: level
      complex(rp), dimension(nx,ny,min_m:max_m), intent(inout) :: level_DR
      complex(rp), dimension(nx,ny,min_m:max_m), intent(inout) :: kl

      real(rp) :: R

      integer(ip) :: i,j

      select_field: select case (fname)
         !-------------------------------------------------------------------
         case ("psi3")
            do j=1,ny
            do i=1,nx
               R = Rvec(i)

               kl(i,j,m_ang) = &
               -  R*4.0_rp*mu_0(i,j)*level(i,j,m_ang) &
               -  R*ta_0(i,j)*(psi4_f%level(i,j,m_ang)) &
               +    (psi4_f%edth(i,j,m_ang))
            end do 
            end do
         !--------------------------------------------------------------------
         case ("la")
            do j=1,ny
            do i=1,nx
               R = Rvec(i)

               kl(i,j,m_ang) = &
               -  R*(mu_0(i,j) + conjg(mu_0(i,j)))*level(i,j,m_ang) &
               -    (psi4_f%level(i,j,m_ang))
            end do
            end do
         !--------------------------------------------------------------------
         case ("psi2")
            do j=1,ny
            do i=1,nx
               R = Rvec(i)

               kl(i,j,m_ang) = &
               -  R*(3.0_rp*mu_0(i,j))*level(i,j,m_ang) &
               -  R*2.0_rp*ta_0(i,j)*(psi3%level(i,j,m_ang)) &
               +    (psi3%edth(i,j,m_ang))
            end do
            end do
         !--------------------------------------------------------------------
         case ("hmbmb")
            do j=1,ny
            do i=1,nx
               R = Rvec(i)

               kl(i,j,m_ang) = & 
                  R*(mu_0(i,j) - conjg(mu_0(i,j)))*level(i,j,m_ang) &
               -    2.0_rp*(la%level(i,j,m_ang))
            end do
            end do
         !--------------------------------------------------------------------
         case ("pi")
            do j=1,ny
            do i=1,nx
               R = Rvec(i)

               kl(i,j,m_ang) = &
               -       R*(conjg(pi_0(i,j)) + ta_0(i,j))*(la%level(i,j,m_ang)) &
               +  (R**2)*0.5_rp*mu_0(i,j)*(conjg(pi_0(i,j))+ta_0(i,j))*(hmbmb%level(i,j,m_ang)) &
               -         (psi3%level(i,j,m_ang))
            end do
            end do
         !--------------------------------------------------------------------
         case ("hlmb")
            do j=1,ny
            do i=1,nx
               R = Rvec(i)

               kl(i,j,m_ang) = &
               -  R*conjg(mu_0(i,j))*level(i,j,m_ang) &
               -    2.0_rp*(pi%level(i,j,m_ang)) &
               -  R*ta_0(i,j)*(hmbmb%level(i,j,m_ang)) 
            end do
            end do
         !--------------------------------------------------------------------
         case ("muhll")
            do j=1,ny
            do i=1,nx
               R = Rvec(i)

               kl(i,j,m_ang) = &
               -       R*conjg(mu_0(i,j))*level(i,j,m_ang) &
               -       R*mu_0(i,j)*(hlmb%edth(i,j,m_ang)) &
               -  (R**2)*mu_0(i,j)*(conjg(pi_0(i,j))+2.0_rp*ta_0(i,j))*(hlmb%level(i,j,m_ang)) &
               -         2.0_rp*(pi%edth(i,j,m_ang)) &
               -       R*2.0_rp*conjg(pi_0(i,j))*(pi%level(i,j,m_ang)) &
               -         2.0_rp*(psi2%level(i,j,m_ang)) &
               -  (R**2)*pi_0(i,j)*conjg(hmbmb%edth(i,j,-m_ang)) &
               +  (R**2)*(pi_0(i,j)**2)*conjg(hmbmb%level(i,j,-m_ang)) &
               +       R*mu_0(i,j)*conjg(hlmb%edth(i,j,-m_ang))  &
               -  (R**2)*3.0_rp*mu_0(i,j)*pi_0(i,j)*conjg(hlmb%level(i,j,-m_ang)) &
               +  (R**2)*2.0_rp*conjg(mu_0(i,j))*pi_0(i,j)*conjg(hlmb%level(i,j,-m_ang)) &
               -  (R**2)*2.0_rp*mu_0(i,j)*conjg(ta_0(i,j))*conjg(hlmb%level(i,j,-m_ang)) &
               -       R*2.0_rp*pi_0(i,j)*conjg(pi%level(i,j,-m_ang))
            end do
            end do
         !--------------------------------------------------------------------
         case default
            kl = -1.0_rp

      end select select_field
      !-----------------------------------------------------------------------
      call compute_DR(m_ang, level, level_DR)

      do j=1,ny
      do i=1,nx
         R = Rvec(i)

         kl(i,j,m_ang) = kl(i,j,m_ang) &
         -  ((R/cl)**2)*level_DR(i,j,m_ang) &
         -  (falloff*R/(cl**2))*level(i,j,m_ang)
      
         kl(i,j,m_ang) = kl(i,j,m_ang) / (2.0_rp*(1.0_rp+(2.0_rp*bhm*R/(cl**2))))
      end do
      end do

   end subroutine set_k
!=============================================================================
! RK4 evolution
!=============================================================================
   pure subroutine take_step(step, m_ang, f)
      integer(ip), intent(in)    :: step, m_ang
      type(field), intent(inout) :: f
      
      integer(ip) :: i, j

      !-----------------------------------------------------------------------
      select_step: select case (step)
         !--------------------------------------------------------------------
         case (1)
            call set_k(m_ang, f%name, f%falloff, f%n, f%DR, f%k1)

            do j=1,ny
            do i=1,nx
               f%l2(i,j,m_ang)= f%n(i,j,m_ang)+0.5_rp*dt*f%k1(i,j,m_ang)
            end do
            end do
         !--------------------------------------------------------------------
         case (2)
            call set_k(m_ang, f%name, f%falloff, f%l2, f%DR, f%k2)

            do j=1,ny
            do i=1,nx
               f%l3(i,j,m_ang)= f%l2(i,j,m_ang)+0.5_rp*dt*f%k2(i,j,m_ang)
            end do
            end do
         !--------------------------------------------------------------------
         case (3)
            call set_k(m_ang, f%name, f%falloff, f%l3, f%DR, f%k3)

            do j=1,ny
            do i=1,nx
               f%l4(i,j,m_ang)= f%l3(i,j,m_ang)+dt*f%k3(i,j,m_ang)
            end do
            end do
         !--------------------------------------------------------------------
         case (4)
            call set_k(m_ang, f%name, f%falloff, f%l4, f%DR, f%k4)

            do j=1,ny
            do i=1,nx
               f%np1(i,j,m_ang)= f%n(i,j,m_ang) &
               +  (dt/6.0_rp)*(f%k1(i,j,m_ang)+2.0_rp*f%k2(i,j,m_ang)+2.0_rp*f%k3(i,j,m_ang)+f%k4(i,j,m_ang))
            end do
            end do
         !--------------------------------------------------------------------
         case (5)
            call set_k(m_ang, f%name, f%falloff, f%np1, f%DR, f%k5)
         !--------------------------------------------------------------------
         case default
            f%error = -1.0_ip
      end select select_step
      !-----------------------------------------------------------------------
   end subroutine take_step
!=============================================================================
   subroutine set_indep_res(m_ang, fname)
      integer(ip),  intent(in) :: m_ang
      character(*), intent(in) :: fname

      integer(ip) :: i, j
      real(rp) :: R
      !-----------------------------------------------------------------------
      select_field: select case (fname)
         !--------------------------------------------------------------------
         case ("bianci3_res")
            do j=1,ny
            do i=1,nx
               R = Rvec(i)

               bianci3_res%np1(i,j,m_ang) = &
                       R*psi3%edth_prime(i,j,m_ang) &
               +  (R**2)*4.0_rp*pi_0(i,j)*psi3%level(i,j,m_ang) &
               -       R*psi4_f%thorn(i,j,m_ang) &
               +         rh_0(i,j)*psi4_f%level(i,j,m_ang) &
               -  (R**2)*3.0_rp*psi2_0(i,j)*la%level(i,j,m_ang) 
            end do
            end do
         !--------------------------------------------------------------------
         case ("bianci2_res")
            do j=1,ny
            do i=1,nx
               R = Rvec(i)

               bianci2_res%np1(i,j,m_ang) = 0.0_rp 
            end do
            end do

         !--------------------------------------------------------------------
         case ("hll_res")
         
         case default
            continue
      end select select_field
   end subroutine set_indep_res
!=============================================================================
   subroutine step_all_fields(step, m_ang)
      integer(ip), intent(in) :: step, m_ang
      !-----------------------------------------------------------------------
      call set_level(step,m_ang,psi4_f)
      call set_edth( step,m_ang,psi4_f)

      call take_step(step,m_ang,psi3)

      call take_step(step,m_ang,la)
      !-----------------------------------------------------------------------
      call set_level(step,m_ang,psi3)
      call set_edth( step,m_ang,psi3)

      call take_step(step,m_ang,psi2)
      !-----------------------------------------------------------------------
      call set_level(step,m_ang,la)

      call take_step(step,m_ang,hmbmb)
      !-----------------------------------------------------------------------
      call set_level(step,m_ang,hmbmb)

      call take_step(step,m_ang,pi)
      !-----------------------------------------------------------------------
      call set_level(step,m_ang,pi)

      call take_step(step,m_ang,hlmb)
      !-----------------------------------------------------------------------
      call set_level(step,m_ang,psi2)
      call set_level(step,m_ang,hlmb)

      call set_edth( step,m_ang,pi)
      call set_edth( step,m_ang,hlmb)

      call take_step(step,m_ang,muhll)

   end subroutine step_all_fields
!=============================================================================
   subroutine metric_recon_time_step(m_ang)
      integer(ip), intent(in) :: m_ang

      call step_all_fields(1_ip,m_ang)
      call step_all_fields(2_ip,m_ang)
      call step_all_fields(3_ip,m_ang)
      call step_all_fields(4_ip,m_ang)
      call step_all_fields(5_ip,m_ang)

   end subroutine metric_recon_time_step
!=============================================================================
   subroutine metric_recon_indep_res(m_ang)
      integer(ip), intent(in) :: m_ang
      !-----------------------------------------------------------------------
      call set_level(5_ip,m_ang,psi4_f)
      call set_level(5_ip,m_ang,psi3)

      call set_thorn(     5_ip,m_ang,psi4_f)
      call set_edth_prime(5_ip,m_ang,psi3)

      call set_indep_res(m_ang,"bianci3_res")
      !-----------------------------------------------------------------------
      call set_indep_res(m_ang,"bianci2_res")

   end subroutine metric_recon_indep_res
!=============================================================================
end module mod_metric_recon 