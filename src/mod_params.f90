!
! automatically generated from sim_class.py
!
module mod_params
use mod_prec
implicit none
   character(*), parameter :: home_dir = '/home/jripley/teuk-fortran'
   character(*), parameter :: run_type = 'basic_run'
   logical, parameter :: debug = .false.
   character(*), parameter :: computer = 'home'
   real(rp), parameter :: black_hole_mass = 0.5_rp
   real(rp), parameter :: black_hole_spin = 0.35_rp
   real(rp), parameter :: compactification_length = 1.0_rp
   real(rp), parameter :: evolve_time = 30.0_rp
   integer(ip), parameter :: num_saved_times = 100_ip
   character(*), parameter :: initial_data_direction = 'time_symmetric'
   character(*), parameter :: initial_data_type = 'compact_bump'
   real(rp), parameter :: amp_pm = 0.1_rp
   real(rp), parameter :: rl_pm = -1.5_rp
   real(rp), parameter :: ru_pm = 1.5_rp
   integer(ip), parameter :: l_ang_pm = 2_ip
   real(rp), parameter :: amp_nm = 10.0_rp
   real(rp), parameter :: rl_nm = -1.5_rp
   real(rp), parameter :: ru_nm = 1.5_rp
   integer(ip), parameter :: l_ang_nm = 2_ip
   integer(ip), parameter :: pm_ang = 2_ip
   integer(ip), parameter :: psi_spin = -2_ip
   integer(ip), parameter :: psi_boost = -2_ip
   integer(ip), parameter :: nx = 48_ip
   integer(ip), parameter :: nl = 16_ip
   logical, parameter :: metric_recon = .true.
   logical, parameter :: scd_order = .true.
   logical, parameter :: write_indep_res = .true.
   logical, parameter :: write_metric_recon_fields = .false.
   logical, parameter :: write_scd_order_source = .true.
   real(rp), parameter :: start_multiple = 1.0_rp
   character(*), parameter :: walltime = '12:00:00'
   character(*), parameter :: memory = '512'
   character(*), parameter :: num_nodes = '1'
   character(*), parameter :: num_tasks_per_node = '1'
   integer(ip), parameter :: max_l = 15_ip
   real(rp), parameter :: horizon = 0.8570714214271424_rp
   real(rp), parameter :: R_max = 1.1667639067172042_rp
   real(rp), parameter :: scd_order_start_time = 4.810510846598895_rp
   integer(ip), parameter :: ny = 28_ip
   real(rp), parameter :: dt = 0.00390625_rp
   integer(ip), parameter :: nt = 3840_ip
   integer(ip), parameter :: t_step_save = 38_ip
   integer(ip), parameter :: max_m = 4_ip
   integer(ip), parameter :: min_m = -4_ip
   integer(ip), parameter :: max_s = 3_ip
   integer(ip), parameter :: min_s = -3_ip
   character(*), parameter :: output_stem = 'Sat_18_57_bhm0.5_bhs0.35_nx48_ny28_nl16_s-2_lpm2_lnm2_pm2'
   character(*), parameter :: output_dir = 'output/Sat_18_57_bhm0.5_bhs0.35_nx48_ny28_nl16_s-2_lpm2_lnm2_pm2'
   character(*), parameter :: bin_name = 'Sat_18_57_bhm0.5_bhs0.35_nx48_ny28_nl16_s-2_lpm2_lnm2_pm2.run'
   character(*), parameter :: tables_dir = 'output/Sat_18_57_bhm0.5_bhs0.35_nx48_ny28_nl16_s-2_lpm2_lnm2_pm2/tables'
end module mod_params
