module mod_cheb 
!=============================================================================
   use mod_prec
   use mod_field,  only: field, set_level
   use mod_io,     only: set_arr
   use mod_params, only: tables_dir, R_max, nx, ny, min_m, max_m 

!=============================================================================
! for use in fftw
   use, intrinsic :: iso_c_binding, only: c_int, c_double

   implicit none
!=============================================================================
   include 'fftw3.f03'
!=============================================================================
   private

   real(rp) :: dx_over_dR

   integer(ip) :: N

   ! fftw3 objects
   type(c_ptr) :: plan_dct

   ! radial points
   real(rp), allocatable, protected, public :: Rvec(:)
   real(rp), allocatable, protected, public :: Rarr(:,:)

   ! filter array
   real(rp), allocatable :: filter_arr(:,:)

   ! subroutines
   public :: cheb_init, compute_DR, cheb_filter, cheb_test

   public :: cheb_real_to_coef

   interface compute_DR
      module procedure compute_DR_arr, compute_DR_field
   end interface
!=============================================================================
contains
!=============================================================================
   subroutine cheb_init()
      integer :: j
      complex(rp), allocatable :: test_in(:), test_out(:) 

      dx_over_dR = 2.0_rp / R_max
      N = nx - 1

      allocate(Rvec(nx))
      allocate(Rarr(nx,ny))
      allocate(filter_arr(nx,ny))

      call set_arr('cheb_pts.txt', nx, Rvec)

      Rvec = (R_max/2.0_rp) * (Rvec + 1.0_rp)

      do j=1,ny
         Rarr(:,j) = Rvec
      end do

      ! setup dct fftw plan
      allocate(test_in( nx))
      allocate(test_out(nx))

      call dfftw_plan_r2r_1d( &
         plan_dct, &
         nx, &
         test_in,test_out, &
         FFTW_REDFT00,FFTW_PATIENT)
      
      do j=1,nx
         filter_arr(j,:) = exp(-40.0_rp*(real(j-1,rp)/real(nx-1,rp))**16)
      end do

   end subroutine cheb_init
!=============================================================================
   subroutine cheb_real_to_coef(m_ang,vals,coefs,re_v,im_v,re_c,im_c)
      integer(ip), intent(in) :: m_ang
      complex(rp), dimension(nx,ny,min_m:max_m), intent(in)  :: vals
      complex(rp), dimension(nx,ny,min_m:max_m), intent(out) :: coefs
      real(rp),    dimension(nx,ny,min_m:max_m), intent(inout) :: re_v, im_v, re_c, im_c

      integer(ip) :: i,j 

      do j=1,ny
         re_v(:,j,m_ang) = real( vals(:,j,m_ang),kind=rp)
         im_v(:,j,m_ang) = aimag(vals(:,j,m_ang))

         call dfftw_execute_r2r(plan_dct,re_v(:,j,m_ang),re_c(:,j,m_ang))
         call dfftw_execute_r2r(plan_dct,im_v(:,j,m_ang),im_c(:,j,m_ang))

         re_c(1,j,m_ang) = re_c(1,j,m_ang)/(2.0_rp*N)
         im_c(1,j,m_ang) = im_c(1,j,m_ang)/(2.0_rp*N)
         
         re_c(nx,j,m_ang) = re_c(nx,j,m_ang)/(2.0_rp*N)
         im_c(nx,j,m_ang) = im_c(nx,j,m_ang)/(2.0_rp*N)

         do i=2,N
            re_c(i,j,m_ang) = re_c(i,j,m_ang)/N
            im_c(i,j,m_ang) = im_c(i,j,m_ang)/N
         end do

         coefs(:,j,m_ang) = cmplx(re_c(:,j,m_ang),im_c(:,j,m_ang),kind=rp)
      end do

   end subroutine cheb_real_to_coef
!=============================================================================
   subroutine cheb_coef_to_real(m_ang,coefs,vals,re_v,im_v,re_c,im_c)
      integer(ip), intent(in) :: m_ang
      complex(rp), dimension(nx,ny,min_m:max_m), intent(in)  :: coefs
      complex(rp), dimension(nx,ny,min_m:max_m), intent(out) :: vals
      real(rp),    dimension(nx,ny,min_m:max_m), intent(inout) :: re_v, im_v, re_c, im_c

      integer(ip) :: i, j 

      do j=1,ny
         re_c(:,j,m_ang) = real( coefs(:,j,m_ang),kind=rp)
         im_c(:,j,m_ang) = aimag(coefs(:,j,m_ang))

         do i=2,N
            re_c(i,j,m_ang) = re_c(i,j,m_ang) / 2.0_rp
            im_c(i,j,m_ang) = im_c(i,j,m_ang) / 2.0_rp
         end do

         call dfftw_execute_r2r(plan_dct,re_c(:,j,m_ang),re_v(:,j,m_ang))
         call dfftw_execute_r2r(plan_dct,im_c(:,j,m_ang),im_v(:,j,m_ang))

         vals(:,j,m_ang) = cmplx(re_v(:,j,m_ang),im_v(:,j,m_ang),kind=rp)
      end do

   end subroutine cheb_coef_to_real
!=============================================================================
   subroutine compute_DR_arr(m_ang,vals,coefs,DR,re_v,im_v,re_c,im_c)
      integer(ip), intent(in) :: m_ang
      complex(rp), dimension(nx,ny,min_m:max_m), intent(in)  :: vals
      complex(rp), dimension(nx,ny,min_m:max_m), intent(out) :: coefs
      complex(rp), dimension(nx,ny,min_m:max_m), intent(out) :: DR
      real(rp),    dimension(nx,ny,min_m:max_m), intent(inout) :: re_v, im_v, re_c, im_c

      integer(ip) :: i
   
      call cheb_real_to_coef(m_ang,vals,coefs,re_v,im_v,re_c,im_c)

      ! Any error incurred by setting the last two cheb coefs to zero 
      ! should converge away with higher resolution.
      coefs(nx,  :,m_ang) = 0
      coefs(nx-1,:,m_ang) = 0

      ! use D_vals as a temporary array
      DR(:,:,m_ang) = coefs(:,:,m_ang)

      ! note fortran indexing 1..nx (and not 0..nx-1)
      do i=nx-1, 2, -1
         coefs(i-1,:,m_ang) = 2.0_rp*(i-1)*DR(i,:,m_ang) + coefs(i+1,:,m_ang)
      end do

      coefs(1,:,m_ang) = coefs(1,:,m_ang) / 2.0_rp

      call cheb_coef_to_real(m_ang,coefs,DR,re_v,im_v,re_c,im_c)

      DR(:,:,m_ang) = dx_over_dR * DR(:,:,m_ang)

   end subroutine compute_DR_arr
!=============================================================================
   subroutine compute_DR_field(step,m_ang,f)
      integer(ip), intent(in) :: step, m_ang
      type(field), intent(inout) :: f

      call set_level(step,m_ang,f)
      call compute_DR_arr(m_ang,f%level,f%coefs_cheb,f%DR, &
         f%re,f%im, &
         f%coefs_cheb_re,f%coefs_cheb_im &
      )
   end subroutine compute_DR_field
!=============================================================================
! Low pass filter. A smooth filter appears to help prevent Gibbs-like ringing
!=============================================================================
   subroutine cheb_filter(m_ang, f)
      integer(ip), intent(in)    :: m_ang
      type(field), intent(inout) :: f

      call cheb_real_to_coef(m_ang,f%np1,f%coefs_cheb, &
         f%re,f%im, &
         f%coefs_cheb_re,f%coefs_cheb_im &
      ) 
      f%coefs_cheb(:,:,m_ang) = filter_arr*f%coefs_cheb(:,:,m_ang)

      call cheb_coef_to_real(m_ang,f%coefs_cheb,f%np1, &
         f%re,f%im, &
         f%coefs_cheb_re,f%coefs_cheb_im &
      ) 
   end subroutine cheb_filter
!=============================================================================
   subroutine cheb_test()
      complex(rp), dimension(nx,ny,min_m:max_m) :: &
         vals, coefs, DR_vals, computed_DR_vals
      real(rp), dimension(nx,ny,min_m:max_m) :: re, im, re_c, im_c
      integer(ip) :: i
      integer(ip) :: m_ang = 0_ip

      do i=1,nx
         vals(i,:,:)    = sin(Rvec(i))**2
         DR_vals(i,:,:) = 2*sin(Rvec(i))*cos(Rvec(i))
      end do

      call compute_DR(m_ang, vals, coefs, computed_DR_vals, re, im, re_c, im_c)

      do i=1,nx
         write (*,*) computed_DR_vals(i,:,m_ang) - DR_vals(i,:,m_ang)
      end do


   end subroutine cheb_test
!=============================================================================
end module mod_cheb
