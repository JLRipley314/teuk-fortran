!
! automatically generated from sim_class.py
!
module mod_params
use mod_prec
implicit none
   character(*), parameter :: home_dir = '/home/jripley/teuk-fortran'
   character(*), parameter :: run_type = 'basic_run'
   logical, parameter :: debug = .false.
   real(rp), parameter :: black_hole_mass = 0.5_rp
   real(rp), parameter :: black_hole_spin = 0.35_rp
   real(rp), parameter :: compactification_length = 1.0_rp
   real(rp), parameter :: evolve_time = 150.0_rp
   integer(ip), parameter :: num_saved_times = 500_ip
   integer(ip), parameter :: nx = 160_ip
   integer(ip), parameter :: nl = 30_ip
   logical, parameter :: metric_recon = .true.
   logical, parameter :: scd_order = .true.
   logical, parameter :: constrained_evo = .true.
   logical, parameter :: write_indep_res = .true.
   logical, parameter :: write_metric_recon_fields = .false.
   logical, parameter :: write_scd_order_source = .true.
   logical, parameter :: write_coefs = .false.
   character(*), parameter :: computer = 'della'
   character(*), parameter :: della_out_stem = '/tigress/jripley/tf-out/'
   character(*), parameter :: walltime = '24:00:00'
   character(*), parameter :: memory = '2048'
   character(*), parameter :: email = 'lloydripley@gmail.com'
   integer(ip), parameter :: psi_spin = -2_ip
   integer(ip), parameter :: psi_boost = -2_ip
   real(rp), parameter :: start_multiple = 1.0_rp
   integer(ip), dimension(4), parameter :: lin_m = [-2, 2, -3, 3]
   integer(ip), dimension(4), parameter :: l_ang = [2, 2, 3, 3]
   real(rp), dimension(4), parameter :: amp_re = [0.0, 0.4, 0.0, 0.1]
   real(rp), dimension(4), parameter :: amp_im = [0.0, 0.0, 0.0, 0.0]
   real(rp), dimension(4), parameter :: rl_0 = [1.1, 1.1, 1.1, 1.1]
   real(rp), dimension(4), parameter :: ru_0 = [2.5, 2.5, 2.5, 2.5]
   character(*), parameter :: initial_data_direction = 'iiii'
   integer(ip), dimension(3), parameter :: scd_m = [0, -1, 1]
   integer(ip), dimension(4), parameter :: write_lin_m = [-2, 2, -3, 3]
   integer(ip), dimension(3), parameter :: write_scd_m = [0, -1, 1]
   integer(ip), parameter :: max_l = 29_ip
   real(rp), parameter :: horizon = 0.8570714214271424_rp
   real(rp), parameter :: R_max = 1.1667639067172042_rp
   real(rp), dimension(4), parameter :: rl = [0.9427785635698568, 0.9427785635698568, 0.9427785635698568, 0.9427785635698568]
   real(rp), dimension(4), parameter :: ru = [2.142678553567856, 2.142678553567856, 2.142678553567856, 2.142678553567856]
   real(rp), parameter :: constraint_damping = 36.51483710615301_rp
   real(rp), parameter :: scd_order_start_time = 4.4037957280297375_rp
   integer(ip), dimension(2), parameter :: lin_pos_m = [2, 3]
   integer(ip), parameter :: len_lin_pos_m = 2_ip
   integer(ip), parameter :: len_lin_m = 4_ip
   integer(ip), parameter :: len_scd_m = 3_ip
   integer(ip), parameter :: len_write_lin_m = 4_ip
   integer(ip), parameter :: len_write_scd_m = 3_ip
   integer(ip), parameter :: num_threads = 2_ip
   integer(ip), parameter :: ny = 48_ip
   real(rp), parameter :: dt = 0.0003515625_rp
   integer(ip), parameter :: nt = 213333_ip
   integer(ip), parameter :: t_step_save = 426_ip
   integer(ip), parameter :: max_m = 3_ip
   integer(ip), parameter :: min_m = -3_ip
   integer(ip), parameter :: max_s = 3_ip
   integer(ip), parameter :: min_s = -3_ip
   character(*), parameter :: output_stem = 'Oct_7_14_14_26_m0.5_s0.35_nx160_nl30'
   character(*), parameter :: output_dir = '/tigress/jripley/tf-out/Oct_7_14_14_26_m0.5_s0.35_nx160_nl30'
   character(*), parameter :: bin_name = 'Oct_7_14_14_26_m0.5_s0.35_nx160_nl30.run'
   character(*), parameter :: tables_dir = '/tigress/jripley/tf-out/Oct_7_14_14_26_m0.5_s0.35_nx160_nl30/tables'
end module mod_params
