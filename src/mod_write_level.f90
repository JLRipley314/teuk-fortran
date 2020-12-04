!
! module for writing fields to file
!
!=============================================================================
module mod_write_level
!=============================================================================
   use, intrinsic :: iso_fortran_env, only: &
      stdout=>output_unit, stdin=>input_unit, stderr=>error_unit

   use mod_prec

   use mod_params, only: nx, ny 
   use mod_field,  only: field

   use mod_io, only: write_csv

   use mod_cheb, only: cheb_real_to_coef
   use mod_swal, only: swal_real_to_coef

   use mod_teuk, only: compute_res_q

   use mod_metric_recon, only: metric_recon_indep_res

   use mod_params, only: &
      write_lin_m, len_write_lin_m, &
      write_scd_m, len_write_scd_m, &
      metric_recon, scd_order, &
      write_metric_recon_fields, &
      write_indep_res, write_scd_order_source, &
      write_coefs, write_coefs_swal, constrained_evo

   use mod_fields_list, only: &
      psi4_lin_p, psi4_lin_q, psi4_lin_f, &
      res_lin_q, & 

      psi4_scd_p, psi4_scd_q, psi4_scd_f, &
      res_scd_q, & 

      psi4_integrated_lin_f, &
      psi4_integrated_scd_f, &

      psi4_twice_integrated_lin_f, &
      psi4_twice_integrated_scd_f, &

      psi3, psi2, la, pi, muhll, hlmb, hmbmb, &
      res_bianchi3, res_bianchi2, res_hll

   use mod_scd_order_source, only: source
!=============================================================================
   implicit none
   private

   public :: write_out, write_diagnostics, write_level 

   interface write_horizon_or_scriplus
      module procedure &
         write_field_horizon_or_scriplus, &
         write_source_horizon_or_scriplus
   end interface
!=============================================================================
contains
!=============================================================================
   subroutine write_field_horizon_or_scriplus(time,location,m_ang,f)
      real(rp),     intent(in) :: time
      character(*), intent(in) :: location
      integer(ip),  intent(in) :: m_ang
      type(field),  intent(inout) :: f
      !----------------------------------------------
      character(:), allocatable :: fn
      complex(rp),  allocatable :: vals(:)

      fn = location//"_"//f%fname 

      !------------------------------------------------
      !------------------------------------------------
      select case (location)
      !------------------------------------------------
      case ("horizon")
         allocate(vals(lbound(f%np1,2):ubound(f%np1,2)))
         vals = f%np1(nx,:,m_ang) 
         call write_csv(fn, time, m_ang, vals)
      !------------------------------------------------
      case ("scriplus")
         allocate(vals(lbound(f%np1,2):ubound(f%np1,2)))
         vals = f%np1(1,:,m_ang) 
         call write_csv(fn, time, m_ang, vals)
      !------------------------------------------------
      case ("coef_horizon")
         allocate(vals(lbound(f%coefs_swal,2):ubound(f%coefs_swal,2)))
         call swal_real_to_coef( &
            f%spin, &
            m_ang, &
            f%np1, &
            f%coefs_swal &
         )
         vals = f%coefs_swal(nx,:,m_ang) 
         call write_csv(fn, time, m_ang, vals)
      !------------------------------------------------
      case ("coef_scriplus")
         allocate(vals(lbound(f%coefs_swal,2):ubound(f%coefs_swal,2)))
         call swal_real_to_coef( &
            f%spin, &
            m_ang, &
            f%np1, &
            f%coefs_swal &
         )
         vals = f%coefs_swal(1,:,m_ang) 
         call write_csv(fn, time, m_ang, vals)
      !------------------------------------------------
      case default
         ! do nothing
      !------------------------------------------------
      end select 
      !------------------------------------------------
      !----------------------------------------------
   end subroutine write_field_horizon_or_scriplus
!=============================================================================
   subroutine write_source_horizon_or_scriplus(time,location,m_ang)
      real(rp),     intent(in) :: time
      character(*), intent(in) :: location
      integer(ip),  intent(in) :: m_ang
      !----------------------------------------------
      character(:), allocatable :: fn
      complex(rp),  allocatable :: vals(:)

      fn = location//"_"//source%fname 

      !------------------------------------------------
      !------------------------------------------------
      select case (location)
      !------------------------------------------------
      case ("horizon")
         allocate(vals(lbound(source%np1,2):ubound(source%np1,2)))
         vals = source%np1(nx,:,m_ang) 
         call write_csv(fn, time, m_ang, vals)
      !------------------------------------------------
      case ("scriplus")
         allocate(vals(lbound(source%np1,2):ubound(source%np1,2)))
         vals = source%np1(1,:,m_ang) 
         call write_csv(fn, time, m_ang, vals)
      !------------------------------------------------
      case ("coef_horizon")
         allocate(vals(lbound(source%coefs_swal,2):ubound(source%coefs_swal,2)))
         call swal_real_to_coef( &
            source%spin, &
            m_ang, &
            source%np1, &
            source%coefs_swal &
         )
         vals = source%coefs_swal(nx,:,m_ang) 
         call write_csv(fn, time, m_ang, vals)
      !------------------------------------------------
      case ("coef_scriplus")
         allocate(vals(lbound(source%coefs_swal,2):ubound(source%coefs_swal,2)))
         call swal_real_to_coef( &
            source%spin, &
            m_ang, &
            source%np1, &
            source%coefs_swal &
         )
         vals = source%coefs_swal(1,:,m_ang) 
         call write_csv(fn, time, m_ang, vals)
      !------------------------------------------------
      case default
         ! do nothing
      end select 
      !------------------------------------------------
      !------------------------------------------------
   end subroutine write_source_horizon_or_scriplus
!=============================================================================
! computes two norm
!-----------------------------------------------------------------------------
   subroutine write_norm(time,m_ang,f)
      real(rp),     intent(in) :: time
      integer(ip),  intent(in) :: m_ang
      type(field),  intent(in) :: f
      !----------------------------------------------
      character(:), allocatable :: fn
      real(rp)                  :: norm

      fn = "norm_"//f%fname 

      norm = sum(abs(conjg(f%np1(:,:,m_ang))*f%np1(:,:,m_ang)))

      norm = sqrt(norm / real(size(f%np1(:,:,m_ang)),kind=rp))

      call write_csv(fn, time, m_ang, norm)

   end subroutine write_norm
!=============================================================================
   subroutine write_out(time)
      real(rp), intent(in) :: time

      write(stdout,'(e14.6)') time

   end subroutine write_out
!=============================================================================
   subroutine write_diagnostics(time)
      real(rp), intent(in) :: time

      integer(ip) :: i 
      !-----------------------------------------------------------------------
      write(stdout,'(e14.6)') time
      !-----------------------------------------------------------------------
      ! field values at future null infinity and horizon
      !-----------------------------------------------------------------------
      do i=1,len_write_lin_m
         call write_horizon_or_scriplus(time,"horizon", write_lin_m(i),psi4_lin_f)
         call write_horizon_or_scriplus(time,"scriplus",write_lin_m(i),psi4_lin_f)
         call write_horizon_or_scriplus(time,"scriplus",write_lin_m(i),psi4_integrated_lin_f)
         call write_horizon_or_scriplus(time,"scriplus",write_lin_m(i),psi4_twice_integrated_lin_f)
      end do
      !-----------------------------------------------------------------------
      if (write_metric_recon_fields) then
         do i=1,len_write_lin_m
            call write_horizon_or_scriplus(time,"horizon", write_lin_m(i),psi3)
            call write_horizon_or_scriplus(time,"scriplus",write_lin_m(i),psi3)

            call write_horizon_or_scriplus(time,"horizon", write_lin_m(i),psi2)
            call write_horizon_or_scriplus(time,"scriplus",write_lin_m(i),psi2)

            call write_horizon_or_scriplus(time,"horizon", write_lin_m(i),hmbmb)
            call write_horizon_or_scriplus(time,"scriplus",write_lin_m(i),hmbmb)

            call write_horizon_or_scriplus(time,"horizon", write_lin_m(i),hlmb)
            call write_horizon_or_scriplus(time,"scriplus",write_lin_m(i),hlmb)

            call write_horizon_or_scriplus(time,"horizon", write_lin_m(i),muhll)
            call write_horizon_or_scriplus(time,"scriplus",write_lin_m(i),muhll)
         end do
      end if
      !--------------------------------------------------------------------
      if (write_indep_res) then
         do i=1,len_write_lin_m
            call metric_recon_indep_res(write_lin_m(i))

            call write_norm(time,write_lin_m(i),res_bianchi3)
            call write_norm(time,write_lin_m(i),res_bianchi2)
            call write_norm(time,write_lin_m(i),res_hll)
         end do
      end if
      !-----------------------------------------------------------------------
      if (  (write_indep_res) &
      .and. (.not. constrained_evo) &
      ) then
         do i=1,len_write_lin_m
            call compute_res_q( write_lin_m(i),psi4_lin_q,psi4_lin_f,res_lin_q)
            call write_norm(time,write_lin_m(i),res_lin_q)
         end do
      end if
      !-----------------------------------------------------------------------
      if (scd_order) then
         !--------------------------------------------------------------------
         if (  (write_indep_res) &
         .and. (.not. constrained_evo) &
         ) then
            do i=1,len_write_scd_m
               call compute_res_q( write_scd_m(i),psi4_scd_q,psi4_scd_f,res_scd_q)
               call write_norm(time,write_scd_m(i),res_scd_q)
            end do
         end if
         !--------------------------------------------------------------------
         do i=1,len_write_scd_m
            call write_horizon_or_scriplus(time,"horizon", write_scd_m(i),psi4_scd_f)
            call write_horizon_or_scriplus(time,"scriplus",write_scd_m(i),psi4_scd_f)
            call write_horizon_or_scriplus(time,"scriplus",write_scd_m(i),psi4_integrated_scd_f)
            call write_horizon_or_scriplus(time,"scriplus",write_scd_m(i),psi4_twice_integrated_scd_f)
         end do

         if (write_scd_order_source) then 
            do i=1,len_write_scd_m
               call write_horizon_or_scriplus(time,"horizon",write_scd_m(i))
            end do
         end if
      end if
      !-----------------------------------------------------------------------
      if (write_coefs_swal) then
         !--------------------------------------------------------------------
         do i=1,len_write_lin_m
            call write_horizon_or_scriplus(time,"coef_horizon", write_lin_m(i),psi4_lin_f)
            call write_horizon_or_scriplus(time,"coef_scriplus",write_lin_m(i),psi4_lin_f)
            call write_horizon_or_scriplus(time,"coef_scriplus",write_lin_m(i),psi4_integrated_lin_f)
            call write_horizon_or_scriplus(time,"coef_scriplus",write_lin_m(i),psi4_twice_integrated_lin_f)
         end do 
         !--------------------------------------------------------------------
         if (scd_order) then
            do i=1,len_write_scd_m
               call write_horizon_or_scriplus(time,"coef_horizon", write_scd_m(i),psi4_scd_f)
               call write_horizon_or_scriplus(time,"coef_scriplus",write_scd_m(i),psi4_scd_f)
               call write_horizon_or_scriplus(time,"coef_scriplus",write_scd_m(i),psi4_integrated_scd_f)
               call write_horizon_or_scriplus(time,"coef_scriplus",write_scd_m(i),psi4_twice_integrated_scd_f)
            end do 
         end if
         !--------------------------------------------------------------------
         if (write_scd_order_source) then 
            do i=1,len_write_scd_m
               call write_horizon_or_scriplus(time,"coef_horizon",write_scd_m(i))
            end do
         end if
         !--------------------------------------------------------------------
      end if

   end subroutine write_diagnostics
!=============================================================================
   subroutine write_level(time)
      real(rp), intent(in) :: time

      integer(ip) :: i
      !-----------------------------------------------------------------------
      ! \Psi_4^{(1)} and linear metric reconstruction 
      !-----------------------------------------------------------------------
      do i=1,len_write_lin_m
         call write_csv(time,write_lin_m(i),psi4_lin_f)
      end do 
      !-----------------------------------------------------------------------
      if (write_metric_recon_fields) then
         do i=1,len_write_lin_m
            call write_csv(time,write_lin_m(i),psi3)
            call write_csv(time,write_lin_m(i),psi2)
            call write_csv(time,write_lin_m(i),hmbmb)
            call write_csv(time,write_lin_m(i),hlmb)
            call write_csv(time,write_lin_m(i),muhll)
         end do
      end if
      !-----------------------------------------------------------------------
      if (write_indep_res) then
         if (.not. constrained_evo) then
            do i=1,len_write_lin_m
               call compute_res_q( write_lin_m(i),psi4_lin_q,psi4_lin_f,res_lin_q)
               call write_csv(time,write_lin_m(i),res_lin_q)
            end do
         end if
         if (metric_recon) then
            do i=1,len_write_lin_m
               call metric_recon_indep_res(write_lin_m(i))

               call write_csv(time,write_lin_m(i),res_bianchi3)
               call write_csv(time,write_lin_m(i),res_bianchi2)
               call write_csv(time,write_lin_m(i),res_hll)
            end do
         end if
      end if
      !-----------------------------------------------------------------------
      ! \Psi_4^{(2)} and 2nd order source term 
      !-----------------------------------------------------------------------
      if (scd_order) then

         do i=1,len_write_scd_m
            call write_csv(time,write_scd_m(i),psi4_scd_f)
         end do

         if (write_indep_res) then
            if (.not. constrained_evo) then
               do i=1,len_write_scd_m
                  call compute_res_q( write_scd_m(i),psi4_scd_q,psi4_scd_f,res_scd_q)
                  call write_csv(time,write_scd_m(i),res_scd_q)
               end do
            end if
         end if

         if (write_scd_order_source) then 
            do i=1,len_write_scd_m
               call write_csv( &
                  source%fname, &
                  time, &
                  write_scd_m(i), &
                  source%np1(:,:,write_scd_m(i)) &
               )
            end do
         end if

      end if
      !-----------------------------------------------------------------------
      if (write_coefs) then
         !--------------------------------------------------------------------
         do i=1,len_write_lin_m
            call cheb_real_to_coef( &
               write_lin_m(i), &
               psi4_lin_f%np1, &
               psi4_lin_f%coefs_cheb, &
               psi4_lin_f%re, &
               psi4_lin_f%im, &
               psi4_lin_f%coefs_cheb_re, &
               psi4_lin_f%coefs_cheb_im  &
            )
            call swal_real_to_coef( &
               psi4_lin_f%spin, &
               write_lin_m(i), &
               psi4_lin_f%coefs_cheb, &
               psi4_lin_f%coefs_both &
            )
            call write_csv( &
               "coefs_"//psi4_lin_f%fname, &
               time, &
               write_lin_m(i), &
               psi4_lin_f%coefs_both(:,:,write_lin_m(i)) &
            )
         end do 
         !--------------------------------------------------------------------
         if (scd_order) then
            do i=1,len_write_scd_m
               call cheb_real_to_coef( &
                  write_scd_m(i), &
                  psi4_scd_f%np1, &
                  psi4_scd_f%coefs_cheb, &
                  psi4_scd_f%re, &
                  psi4_scd_f%im, &
                  psi4_scd_f%coefs_cheb_re, &
                  psi4_scd_f%coefs_cheb_im  &
               )
               call swal_real_to_coef( &
                  psi4_scd_f%spin, &
                  write_scd_m(i), &
                  psi4_scd_f%coefs_cheb, &
                  psi4_scd_f%coefs_both &
               )
               call write_csv( &
                  "coefs_"//psi4_scd_f%fname, &
                  time, &
                  write_scd_m(i), &
                  psi4_scd_f%coefs_both(:,:,write_scd_m(i)) &
               )
            end do 
         end if
         !--------------------------------------------------------------------
      end if
      !-----------------------------------------------------------------------
      if (write_coefs_swal) then
         !--------------------------------------------------------------------
         do i=1,len_write_lin_m
            call swal_real_to_coef( &
               psi4_lin_f%spin, &
               write_lin_m(i), &
               psi4_lin_f%np1, &
               psi4_lin_f%coefs_swal &
            )
            call write_csv( &
               "swal_coefs_"//psi4_lin_f%fname, &
               time, &
               write_lin_m(i), &
               psi4_lin_f%coefs_swal(:,:,write_lin_m(i)) &
            )
         end do 
         !--------------------------------------------------------------------
         if (scd_order) then
            do i=1,len_write_scd_m
               call swal_real_to_coef( &
                  psi4_scd_f%spin, &
                  write_scd_m(i), &
                  psi4_scd_f%np1, &
                  psi4_scd_f%coefs_swal &
               )
               call write_csv( &
                  "swal_coefs_"//psi4_scd_f%fname, &
                  time, &
                  write_scd_m(i), &
                  psi4_scd_f%coefs_swal(:,:,write_scd_m(i)) &
               )
            end do 
         end if
         !--------------------------------------------------------------------
      end if
      !--------------------------------------------------------------------
   end subroutine write_level
!=============================================================================
end module mod_write_level
