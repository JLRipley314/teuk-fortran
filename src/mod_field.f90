!
! Provides the field type, and some basic utility functions 
!
module mod_field
!=============================================================================
   use mod_prec
   use mod_params, only: nx, ny, max_l, min_m, max_m, dt
   implicit none
!=============================================================================
   private

   public :: field, &
      set_field, set_level, set_DT, &
      shift_time_step, &
      time_integrate_field
!=============================================================================
! field type:
! always have two levels: n, np1 and intermediate levels for RK4 time
! evolution. k1, etc. are the time derivatives (note 1==n, 5==np1).
! These functions are complex as we work in NP formalism
!-----------------------------------------------------------------------------
   type :: field

   character(:), allocatable :: fname

   logical :: first_time

   integer(ip) :: spin, boost, falloff

   complex(rp), allocatable :: &
      n(:,:,:),  l2(:,:,:), l3(:,:,:), l4(:,:,:), np1(:,:,:), &
      k1(:,:,:), k2(:,:,:), k3(:,:,:), k4(:,:,:), k5(:,:,:),  &

      level(:,:,:), &
      DT(:,:,:), &
      DR(:,:,:), &
      raised(:,:,:), &
      lowered(:,:,:), & 
      swal_lap(:,:,:), &

      coefs_cheb(:,:,:), &

      edth(:,:,:), &
      edth_prime(:,:,:), &
      thorn(:,:,:), &
      thorn_prime(:,:,:), & 

      coefs_swal(:,:,:), coefs_both(:,:,:) 

   real(rp), allocatable :: &
      re(:,:,:), im(:,:,:), coefs_cheb_re(:,:,:), coefs_cheb_im(:,:,:) 

   end type field
!=============================================================================
contains
!=============================================================================
   subroutine set_field(fname, spin, boost, falloff, f)
      character(*), intent(in)    :: fname ! field fname
      integer(ip),  intent(in)    :: spin, boost, falloff
      type(field),  intent(inout) :: f

      f % fname = fname 

      f % first_time = .true.

      f % spin    = spin
      f % boost   = boost
      f % falloff = falloff

      allocate(f%n(  nx,ny,min_m:max_m))
      allocate(f%l2( nx,ny,min_m:max_m))
      allocate(f%l3( nx,ny,min_m:max_m))
      allocate(f%l4( nx,ny,min_m:max_m))
      allocate(f%np1(nx,ny,min_m:max_m))
      f % n   = 0.0_rp
      f % l2  = 0.0_rp
      f % l3  = 0.0_rp
      f % l4  = 0.0_rp
      f % np1 = 0.0_rp

      allocate(f%k1(nx,ny,min_m:max_m))
      allocate(f%k2(nx,ny,min_m:max_m))
      allocate(f%k3(nx,ny,min_m:max_m))
      allocate(f%k4(nx,ny,min_m:max_m))
      allocate(f%k5(nx,ny,min_m:max_m))
      f % k1 = 0.0_rp
      f % k2 = 0.0_rp
      f % k3 = 0.0_rp
      f % k4 = 0.0_rp
      f % k5 = 0.0_rp

      allocate(f%level(   nx,ny,min_m:max_m))
      allocate(f%DT(      nx,ny,min_m:max_m))
      allocate(f%DR(      nx,ny,min_m:max_m))
      allocate(f%raised(  nx,ny,min_m:max_m))
      allocate(f%lowered( nx,ny,min_m:max_m))
      allocate(f%swal_lap(nx,ny,min_m:max_m))
      f % level    = 0.0_rp
      f % DT       = 0.0_rp
      f % DR       = 0.0_rp
      f % raised   = 0.0_rp
      f % lowered  = 0.0_rp
      f % swal_lap = 0.0_rp

      allocate(f%coefs_swal(nx,0:max_l,min_m:max_m))
      allocate(f%coefs_cheb(nx,     ny,min_m:max_m))
      allocate(f%coefs_both(nx,0:max_l,min_m:max_m))
      f % coefs_swal = 0.0_rp 
      f % coefs_cheb = 0.0_rp 
      f % coefs_both = 0.0_rp 

      allocate(f%edth(       nx,ny,min_m:max_m))
      allocate(f%edth_prime( nx,ny,min_m:max_m))
      allocate(f%thorn(      nx,ny,min_m:max_m))
      allocate(f%thorn_prime(nx,ny,min_m:max_m))
      f % edth        = 0.0_rp
      f % edth_prime  = 0.0_rp
      f % thorn       = 0.0_rp
      f % thorn_prime = 0.0_rp

      allocate(f%re(           nx,ny,min_m:max_m))
      allocate(f%im(           nx,ny,min_m:max_m))
      allocate(f%coefs_cheb_re(nx,ny,min_m:max_m))
      allocate(f%coefs_cheb_im(nx,ny,min_m:max_m))
      f % re = 0.0_rp
      f % im = 0.0_rp
      f % coefs_cheb_re = 0.0_rp
      f % coefs_cheb_im = 0.0_rp

   end subroutine set_field
!=============================================================================
   subroutine set_level(step, m_ang, f)
      integer(ip), intent(in)    :: step, m_ang
      type(field), intent(inout) :: f

      select case (step)
         case (1)
            f % level(:,:,m_ang) = f % n(:,:,m_ang) 
         case (2)
            f % level(:,:,m_ang) = f % l2(:,:,m_ang) 
         case (3)
            f % level(:,:,m_ang) = f % l3(:,:,m_ang) 
         case (4)
            f % level(:,:,m_ang) = f % l4(:,:,m_ang) 
         case (5)
            f % level(:,:,m_ang) = f % np1(:,:,m_ang)  
         case default
            write(*,*) "ERROR(set_field_level), step=", step 
      end select
   end subroutine set_level
!=============================================================================
   subroutine set_DT(step, m_ang, f)
      integer(ip), intent(in)    :: step, m_ang
      type(field), intent(inout) :: f

      select case (step)
         case (1)
            f % DT(:,:,m_ang) = f % k1(:,:,m_ang) 
         case (2)
            f % DT(:,:,m_ang) = f % k2(:,:,m_ang) 
         case (3)
            f % DT(:,:,m_ang) = f % k3(:,:,m_ang) 
         case (4)
            f % DT(:,:,m_ang) = f % k4(:,:,m_ang) 
         case (5)
            f % DT(:,:,m_ang) = f % k5(:,:,m_ang) 
         case default
            write(*,*) "ERROR(set_field_DT), step=", step 
      end select
   end subroutine set_DT
!=============================================================================
   subroutine shift_time_step(m_ang, f)
      integer(ip), intent(in)    :: m_ang
      type(field), intent(inout) :: f 

      f % n( :,:,m_ang) = f % np1(:,:,m_ang) 
      f % k1(:,:,m_ang) = f % k5( :,:,m_ang) 
   end subroutine shift_time_step
!=============================================================================
! RK4 time integrator: integrate field in time 
!=============================================================================
   subroutine time_integrate_field(m_ang, f, int_f) 
      integer(ip), intent(in)    :: m_ang
      type(field), intent(in)    :: f
      type(field), intent(inout) ::int_f ! int_f :time integral of f 
   !--------------------------------------------------------
      int_f%k1(:,:,m_ang)= f%n(:,:,m_ang)
      int_f%l2(:,:,m_ang)= int_f%n(:,:,m_ang)+0.5_rp*dt*int_f%k1(:,:,m_ang)
   !--------------------------------------------------------
      int_f%k2(:,:,m_ang)= f%l2(:,:,m_ang)
      int_f%l3(:,:,m_ang)= int_f%n(:,:,m_ang)+0.5_rp*dt*int_f%k2(:,:,m_ang)
   !--------------------------------------------------------
      int_f%k3(:,:,m_ang)= f%l3(:,:,m_ang)
      int_f%l4(:,:,m_ang)= int_f%n(:,:,m_ang)+dt*int_f%k3(:,:,m_ang)
   !--------------------------------------------------------
      int_f%k4(:,:,m_ang) = f%l4(:,:,m_ang)
      int_f%np1(:,:,m_ang)= int_f%n(:,:,m_ang) &
      +  (dt/6.0_rp)*( &
            int_f%k1(:,:,m_ang) &
         +  2.0_rp*int_f%k2(:,:,m_ang) &
         +  2.0_rp*int_f%k3(:,:,m_ang) &
         +  int_f%k4(:,:,m_ang) &
         )
   !------------------------------------------------------------
   ! want k5 for computing source term and independent residuals
   !------------------------------------------------------------
      int_f%k5(:,:,m_ang) = f%np1(:,:,m_ang)

   end subroutine time_integrate_field
!=============================================================================
end module mod_field
