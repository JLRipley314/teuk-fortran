#!/usr/bin/env python

import sys, time
from sim_class import Sim
#=============================================================================
args= sys.argv
#=============================================================================
## input paramters: set these by hand 
#=============================================================================
sim= Sim(args)
#-----------------------------------------------------------------------------
#sim.base_output_dir= '/mnt/grtheory/second_order_pert/'
sim.base_output_dir= '/mnt/grtheory/test_second_order_pert/'
#-----------------------------------------------------------------------------
sim.black_hole_mass= float(0.5)	
sim.black_hole_spin= float(0.35)
sim.compactification_length= float(1)
#-----------------------------------------------------------------------------
## evolve time: in units of black hole mass
#-----------------------------------------------------------------------------
sim.evolve_time= float(20)
sim.num_saved_times= int(50)
#-----------------------------------------------------------------------------
## initial data
#-----------------------------------------------------------------------------
sim.initial_data_direction= str("ingoing")#str("time_symmetric")#
sim.initial_data_type= str("compact_bump")
sim.amp_pm= float(10.0)  ## amplitude of the initial perturbation
sim.rl_pm=  float(-1.5)  ## compact support: lower r value
sim.ru_pm=  float( 1.5)  ## compact support: upper r value 
sim.l_ang_pm= int(2)     ## support over single spin weighted spherical harmonic

sim.amp_nm= float( 0.0)  ## amplitude of the initial perturbation
sim.rl_nm=  float(-1.5)  ## compact support: lower r value
sim.ru_nm=  float( 1.5)  ## compact support: upper r value 
sim.l_ang_nm= int(2)     ## support over single spin weighted spherical harmonic
#-----------------------------------------------------------------------------
##  Teukolsky equation preserves m 
sim.pm_ang= int(2)
assert(sim.pm_ang>=0)
#-----------------------------------------------------------------------------
## psi_4 is spin -2, psi_0 is spin +2 (code only reconstructs for psi_4) 
sim.spin= int(-2)
#-----------------------------------------------------------------------------
sim.nx= int(pow(2,6)*pow(3,0)*pow(5,0)*pow(7,0)) ## num radial pts 
sim.nl= int(pow(2,4)*pow(3,0)*pow(5,0)*pow(7,0)) ## num swaL angular pts 
#-----------------------------------------------------------------------------
## further diagnostics
sim.save_indep_res= "true"#"false"#
sim.save_coefs= "true"#"false"#
sim.save_metric= "true"#"false"#
sim.save_source= "true"#"false"#
#-----------------------------------------------------------------------------
## change start time
sim.start_multiple= float(1.0)
#-----------------------------------------------------------------------------
sim.set_derived_params()
#=============================================================================
sim.walltime= '168:00:00' ### (hh:mm:ss)
sim.memory=  '512' ### MB 
#=============================================================================
if (sim.run_type == "basic_run"):
	sim.set_derived_params()
	sim.make_output_dir()
	sim.make_swaL_dir()
	sim.make_Legendre_pts()
	sim.make_Gauss_pts()
	sim.write_sim_params()
	sim.launch_run()
#=============================================================================
elif (sim.run_type == "convergence_test"):
	nx_vals= [64, 80, 96]
	nl_vals= [16, 20, 24]

#	nx_vals= [80, 96,112]
#	nl_vals= [16, 20, 24]

#	nx_vals= [112, 128, 144]
#	nl_vals= [24,   28,  32]

	for i in range(len(nx_vals)):
		sim.nx= nx_vals[i]
		sim.nl= nl_vals[i]

		sim.set_derived_params()
		sim.make_output_dir()
		sim.make_swaL_dir()
		sim.make_Legendre_pts()
		sim.make_Gauss_pts()
		sim.write_sim_params()
		sim.launch_run()

		time.sleep(1)
#=============================================================================
else:
	raise ValueError("run_type = "+str(sim.run_type)) 
