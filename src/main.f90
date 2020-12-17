!
! Evolve linear Teukolsky field, reconstruct metric, and evolve
! second order Teukolsky field
!
!=============================================================================
program main
!=============================================================================
   use, intrinsic :: iso_fortran_env, only: stdout=>output_unit

   use mod_prec
   use mod_params, only: &
      read_params, &
      sparse_save, &
      nt, dt, t_step_save, black_hole_mass, &
      lin_m,     scd_m, &
      lin_pos_m, &
      len_lin_m, len_scd_m, &
      len_lin_pos_m, &
      psi_spin, psi_boost, &
      metric_recon, scd_order, &
      integrate_psi4_start_time, &
      scd_order_start_time

   use mod_field,        only: set_field, shift_time_step, time_integrate_field
   use mod_cheb,         only: cheb_init, cheb_filter, cheb_test
   use mod_swal,         only: swal_init, swal_filter, swal_test_orthonormal
   use mod_ghp,          only: ghp_init
   use mod_teuk,         only: teuk_init, teuk_time_step
   use mod_initial_data, only: set_initial_data
   use mod_bkgrd_np,     only: bkgrd_np_init
   use mod_metric_recon, only: &
      metric_recon_time_step_preserve_m_ang, &
      metric_recon_time_step_mix_m_ang
   use mod_write_level,  only: write_level, write_diagnostics

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

   use mod_scd_order_source, only: &
      scd_order_source, &
      source, &
      scd_order_source_init, &
      set_scd_order_source_fields, &
      scd_order_source_compute, &
      scd_order_source_shift_time_step

   implicit none
!=============================================================================
! Put everything in a block so valgrind doesn't get confused about 
! automatically deallocated memory
!=============================================================================
clean_memory: block
!=============================================================================
! declare and initialize variables, fields, etc.
!=============================================================================
   integer(ip) :: m_i, t_step
   real(rp)    :: time
!=============================================================================
   write (*,*) "Reading in params"
   call read_params()
!=============================================================================
   write (*,*) "Initializing fields"   
!-----------------------------------------------------------------------------
! first order metric field
!-----------------------------------------------------------------------------
   call set_field(fname="lin_p",spin=psi_spin,boost=psi_boost,falloff=1_ip,f=psi4_lin_p)
   call set_field(fname="lin_q",spin=psi_spin,boost=psi_boost,falloff=2_ip,f=psi4_lin_q)
   call set_field(fname="lin_f",spin=psi_spin,boost=psi_boost,falloff=1_ip,f=psi4_lin_f)

   call set_field(fname="integrated_lin_f", &
      spin=psi_spin,boost=psi_boost,falloff=1_ip,f=psi4_integrated_lin_f &
   )
   call set_field(fname="twice_integrated_lin_f", &
      spin=psi_spin,boost=psi_boost,falloff=1_ip,f=psi4_twice_integrated_lin_f &
   )
   if (scd_order) then
      call set_field(fname="scd_p",spin=psi_spin,boost=psi_boost,falloff=1_ip,f=psi4_scd_p)
      call set_field(fname="scd_q",spin=psi_spin,boost=psi_boost,falloff=2_ip,f=psi4_scd_q)
      call set_field(fname="scd_f",spin=psi_spin,boost=psi_boost,falloff=1_ip,f=psi4_scd_f)

      call set_field(fname="integrated_scd_f", &
         spin=psi_spin,boost=psi_boost,falloff=1_ip,f=psi4_integrated_scd_f &
      )
      call set_field(fname="twice_integrated_scd_f", &
         spin=psi_spin,boost=psi_boost,falloff=1_ip,f=psi4_twice_integrated_scd_f &
      )
   end if
!-----------------------------------------------------------------------------
! metric reconstructed fields
!-----------------------------------------------------------------------------
   call set_field(fname="psi3",spin=-1_ip,boost=-1_ip,falloff=2_ip,f=psi3)
   call set_field(fname="psi2",spin= 0_ip,boost= 0_ip,falloff=3_ip,f=psi2)

   call set_field(fname="la",spin=-2_ip,boost=-1_ip,falloff=1_ip,f=la)
   call set_field(fname="pi",spin=-1_ip,boost= 0_ip,falloff=2_ip,f=pi)

   call set_field(fname="muhll",spin= 0_ip,boost=1_ip,falloff=3_ip,f=muhll)
   call set_field(fname="hlmb" ,spin=-1_ip,boost=1_ip,falloff=2_ip,f= hlmb)
   call set_field(fname="hmbmb",spin=-2_ip,boost=0_ip,falloff=1_ip,f=hmbmb)
!-----------------------------------------------------------------------------
! independent residual fields
!-----------------------------------------------------------------------------
   call set_field(fname="res_lin_q",spin=-2_ip,boost=-2_ip,falloff=2_ip,f=res_lin_q)

   if (scd_order) then
      call set_field(fname="res_scd_q",spin=-2_ip,boost=-2_ip,falloff=2_ip,f=res_scd_q)
   end if

   call set_field(fname="res_bianchi3",spin=-2_ip,boost=-1_ip,falloff=2_ip,f=res_bianchi3)
   call set_field(fname="res_bianchi2",spin=-1_ip,boost= 0_ip,falloff=2_ip,f=res_bianchi2)
   call set_field(fname="res_hll",     spin= 0_ip,boost= 2_ip,falloff=2_ip,f=res_hll)
!-----------------------------------------------------------------------------
! source term for \psi_4^{(2)}
!-----------------------------------------------------------------------------
   if (scd_order) then
      call scd_order_source_init(fname="scd_order_source",sf=source)
   end if
!-----------------------------------------------------------------------------
! initialize chebyshev diff matrices, swal matrices, etc.
!-----------------------------------------------------------------------------
   call cheb_init()
   call swal_init()
   call ghp_init()
   call bkgrd_np_init()
   call teuk_init()
!=============================================================================
! initial data 
!=============================================================================
   write (stdout,*) "Setting up initial data"
!-----------------------------------------------------------------------------
   time = 0.0_rp

   do m_i=1,len_lin_m
      call set_initial_data(m_i, psi4_lin_p, psi4_lin_q, psi4_lin_f)
   end do

   if (sparse_save) then
      call write_diagnostics(time / black_hole_mass)
   else
      call write_level(      time / black_hole_mass)
   end if
!=============================================================================
! integrate in time 
!=============================================================================
   write (stdout,*) "Beginning time evolution"
!-----------------------------------------------------------------------------
   time_evolve: do t_step=1,nt
      time = t_step*dt
      !--------------------------------------------------------------------
      ! \Psi_4^{(1)} evolution 
      !--------------------------------------------------------------------
      !$OMP PARALLEL DO NUM_THREADS(len_lin_m) IF(len_lin_m>1)
      preserve_m_evo: do m_i=1,len_lin_m
         call teuk_time_step(lin_m(m_i),psi4_lin_p,psi4_lin_q,psi4_lin_f)
         !------------------------------------
         ! low pass filter (in spectral space)
         !------------------------------------
         call cheb_filter(lin_m(m_i),psi4_lin_p)
         call cheb_filter(lin_m(m_i),psi4_lin_q)
         call cheb_filter(lin_m(m_i),psi4_lin_f)

         call swal_filter(lin_m(m_i),psi4_lin_p)
         call swal_filter(lin_m(m_i),psi4_lin_q)
         call swal_filter(lin_m(m_i),psi4_lin_f)
         !------------------------------------
         if (time > integrate_psi4_start_time) then
            call time_integrate_field(lin_m(m_i),psi4_lin_f,           psi4_integrated_lin_f)
            call time_integrate_field(lin_m(m_i),psi4_integrated_lin_f,psi4_twice_integrated_lin_f)
         end if
      !--------------------------------------------------------------------
      ! metric recon evolves +/- m_ang so only evolve m_ang>=0
      !--------------------------------------------------------------------
         if (metric_recon) then
            call metric_recon_time_step_preserve_m_ang(lin_m(m_i))
            !------------------------------------
            ! low pass filter (in spectral space)
            !------------------------------------
            call cheb_filter(lin_m(m_i),psi3)
            call cheb_filter(lin_m(m_i),psi2)
            call cheb_filter(lin_m(m_i),la)
            call cheb_filter(lin_m(m_i),pi)
            call cheb_filter(lin_m(m_i),hmbmb)
            call cheb_filter(lin_m(m_i),hlmb)

            call swal_filter(lin_m(m_i),psi3)
            call swal_filter(lin_m(m_i),psi2)
            call swal_filter(lin_m(m_i),la)
            call swal_filter(lin_m(m_i),pi)
            call swal_filter(lin_m(m_i),hmbmb)
            call swal_filter(lin_m(m_i),hlmb)
         end if
      end do preserve_m_evo
      !$OMP END PARALLEL DO

      if (metric_recon) then
         !$OMP PARALLEL DO NUM_THREADS(len_lin_pos_m) IF(len_lin_pos_m>1)
         mix_m_evo: do m_i=1,len_lin_pos_m
            call metric_recon_time_step_mix_m_ang(lin_pos_m(m_i))

            call cheb_filter( lin_pos_m(m_i),muhll)
            call swal_filter( lin_pos_m(m_i),muhll)

            if (lin_pos_m(m_i)>0) then
               call cheb_filter(-lin_pos_m(m_i),muhll)
               call swal_filter(-lin_pos_m(m_i),muhll)
            end if
         end do mix_m_evo 
         !$OMP END PARALLEL DO
      end if
      !-----------------------------------------------------------------------
      ! \Psi_4^{(2)} evolution 
      !-----------------------------------------------------------------------
      if (scd_order) then

         !$OMP PARALLEL DO NUM_THREADS(len_lin_m) IF(len_lin_m>1)
         do m_i=1,len_lin_m
            call set_scd_order_source_fields(lin_m(m_i))
         end do
         !$OMP END PARALLEL DO

         !$OMP PARALLEL DO NUM_THREADS(len_scd_m) IF(len_scd_m>1)
         do m_i=1,len_scd_m
            call scd_order_source_compute(scd_m(m_i),source) 
            if (time>=scd_order_start_time) then
               call teuk_time_step(scd_m(m_i),source,psi4_scd_p,psi4_scd_q,psi4_scd_f)
               !------------------------------------
               ! low pass filter (in spectral space)
               !------------------------------------
               call cheb_filter(scd_m(m_i),psi4_scd_p)
               call cheb_filter(scd_m(m_i),psi4_scd_q)
               call cheb_filter(scd_m(m_i),psi4_scd_f)

               call swal_filter(scd_m(m_i),psi4_scd_p)
               call swal_filter(scd_m(m_i),psi4_scd_q)
               call swal_filter(scd_m(m_i),psi4_scd_f)
               !------------------------------------
               if (time > integrate_psi4_start_time) then
                  call time_integrate_field(scd_m(m_i),psi4_scd_f,           psi4_integrated_scd_f)
                  call time_integrate_field(scd_m(m_i),psi4_integrated_scd_f,psi4_twice_integrated_scd_f)
               end if 
            end if
         end do
         !$OMP END PARALLEL DO

      end if
      !-----------------------------------------------------------------------
      ! save to file 
      !-----------------------------------------------------------------------
      if (mod(t_step,t_step_save)==0) then
         if (sparse_save) then
            call write_diagnostics(time / black_hole_mass)
         else
            call write_level(      time / black_hole_mass)
         end if
      end if
      !-----------------------------------------------------------------------
      ! shift time steps
      !-----------------------------------------------------------------------
      do m_i=1,len_lin_m
         call shift_time_step(lin_m(m_i),psi4_lin_p)
         call shift_time_step(lin_m(m_i),psi4_lin_q)
         call shift_time_step(lin_m(m_i),psi4_lin_f)

         call shift_time_step(lin_m(m_i),psi4_integrated_lin_f)
         call shift_time_step(lin_m(m_i),psi4_twice_integrated_lin_f)
      end do

      if (metric_recon) then
         do m_i=1,len_lin_m
            call shift_time_step(lin_m(m_i),psi3)
            call shift_time_step(lin_m(m_i),psi2)

            call shift_time_step(lin_m(m_i),la)
            call shift_time_step(lin_m(m_i),pi)

            call shift_time_step(lin_m(m_i),hmbmb)
            call shift_time_step(lin_m(m_i),hlmb)
            call shift_time_step(lin_m(m_i),muhll) 
         end do 
      end if

      if (scd_order) then
         do m_i=1,len_scd_m
            call shift_time_step(scd_m(m_i),psi4_scd_p)
            call shift_time_step(scd_m(m_i),psi4_scd_q)
            call shift_time_step(scd_m(m_i),psi4_scd_f)

            call shift_time_step(lin_m(m_i),psi4_integrated_scd_f)
            call shift_time_step(lin_m(m_i),psi4_twice_integrated_scd_f)

            call scd_order_source_shift_time_step(scd_m(m_i),source)
         end do
      end if
   end do time_evolve
!=============================================================================
   write (*,*) "Finished evolution"
!=============================================================================
end block clean_memory
!=============================================================================
end program main
