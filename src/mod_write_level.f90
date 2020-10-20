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

   use mod_io, only: write_csv, write_horizon_or_scriplus_csv

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

      psi3, psi2, la, pi, muhll, hlmb, hmbmb, &
      res_bianchi3, res_bianchi2, res_hll

   use mod_scd_order_source, only: source
!=============================================================================
   implicit none
   private

   public :: write_level, write_diagnostics
!=============================================================================
contains
!=============================================================================
   subroutine write_diagnostics(time)
      real(rp), intent(in) :: time

      integer(ip) :: i 
      !-----------------------------------------------------------------------
      write(stdout,'(e14.6)',advance='no') time
      !-----------------------------------------------------------------------
      ! field values at future null infinity and horizon
      !-----------------------------------------------------------------------
      do i=1,len_write_lin_m
         call write_horizon_or_scriplus_csv(time,"horizon", write_lin_m(i),psi4_lin_f)
         call write_horizon_or_scriplus_csv(time,"scriplus",write_lin_m(i),psi4_lin_f)

         call write_horizon_or_scriplus_csv(time,"horizon", write_lin_m(i),psi3)
         call write_horizon_or_scriplus_csv(time,"scriplus",write_lin_m(i),psi3)

         call write_horizon_or_scriplus_csv(time,"horizon", write_lin_m(i),psi2)
         call write_horizon_or_scriplus_csv(time,"scriplus",write_lin_m(i),psi2)

         call write_horizon_or_scriplus_csv(time,"horizon", write_lin_m(i),hmbmb)
         call write_horizon_or_scriplus_csv(time,"scriplus",write_lin_m(i),hmbmb)

         call write_horizon_or_scriplus_csv(time,"horizon", write_lin_m(i),hlmb)
         call write_horizon_or_scriplus_csv(time,"scriplus",write_lin_m(i),hlmb)

         call write_horizon_or_scriplus_csv(time,"horizon", write_lin_m(i),muhll)
         call write_horizon_or_scriplus_csv(time,"scriplus",write_lin_m(i),muhll)
      end do

      do i=1,len_write_scd_m
!         call write_horizon_or_scriplus_csv(time,"horizon",write_scd_m(i),source)

         call write_horizon_or_scriplus_csv(time,"horizon", write_scd_m(i),psi4_scd_f)
         call write_horizon_or_scriplus_csv(time,"scriplus",write_scd_m(i),psi4_scd_f)
      end do

   end subroutine write_diagnostics
!=============================================================================
   subroutine write_level(time)
      real(rp), intent(in) :: time

      integer(ip) :: i
      !-----------------------------------------------------------------------
      ! \Psi_4^{(1)} and linear metric reconstruction 
      !--------------------------------------------------------------------------
      do i=1,len_write_lin_m
         call write_csv(time,write_lin_m(i),psi4_lin_f)
      end do 

      if (write_metric_recon_fields) then
         do i=1,len_write_lin_m
            call write_csv(time,write_lin_m(i),psi3)
            call write_csv(time,write_lin_m(i),psi2)
            call write_csv(time,write_lin_m(i),hmbmb)
            call write_csv(time,write_lin_m(i),hlmb)
            call write_csv(time,write_lin_m(i),muhll)
         end do
      end if

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
