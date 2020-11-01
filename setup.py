#!/usr/bin/env python
#=============================================================================
## parameter file for evolution
## usage:
## ./setup.py [run_type] [debug]
## example:
## ./setup.py basic_run
#=============================================================================
import sys, time
from sim_class import Sim
#=============================================================================
args= sys.argv
sim= Sim(args)
#=============================================================================
sim.bin_name= 'default.run'
sim.recompile= False
#=============================================================================
sim.black_hole_mass= float(0.5)	
sim.black_hole_spin= round(1.0*sim.black_hole_mass,16)
sim.compactification_length= float(1)
#=============================================================================
sim.evolve_time= float(300) ## units of black hole mass
sim.num_saved_times= int(2000)
#=============================================================================
sim.nx= 160 ## num radial pts 
sim.nl= 30 ## num angular values
#=============================================================================
## whether or not only to save horizon/scriplus/norm or full field vals 
sim.sparse_save= True
#=============================================================================
sim.metric_recon=     True
sim.scd_order=        True
sim.constrained_evo = True

sim.write_indep_res=           True
sim.write_metric_recon_fields= False
sim.write_scd_order_source=    True
sim.write_coefs=               False
sim.write_coefs_swal=          False
#=============================================================================
## details of computer setup 

sim.computer= 'della'#'home'#
sim.della_out_stem= '/tigress/jripley/tf-out/'

## for della cluster/slurm script

sim.walltime= '71:00:00' ## (hh:mm:ss)
sim.memory=   '2048' ## MB 
sim.email=    'lloydripley@gmail.com' ## for slurm notification
#=============================================================================
## we can only do metric reconstruction starting from psi4 for now.
## For pure first order Teukolsky evolution we can consider other
## spin weighted fields though.
## psi4 is spin -2, boost -2
## psi3 is spin -1, boost -1
## psi2 is spin  0, boost  0 

sim.psi_spin=  int(-2)
sim.psi_boost= int(-2)
#=============================================================================
## Initial data
#=============================================================================
## l_ang:                  initial data is a particular swal function
## initial_data_direction: which way pulse is approximately "heading"
## amp_re(im):             initial amplitude of real/imaginary parts of psi4
## rl(ru)_0:               lower(upper) bounds of initial data as a multiple
##                         of the black hole horizon
#=============================================================================
## initial data for all linear modes 
#=============================================================================
sim.lin_m=  [-2,    2,   -3,    3,   -4,    4,   -5,    5   ]
sim.l_ang=  [ 2,    2,    3,    3,    4,    4,    5,    5   ]
sim.amp_re= [ 0.1,  0.1,  0.1,  0.1,  0.1,  0.1,  0.1,  0.1 ]
sim.amp_im= [ 0.1,  0.1,  0.1,  0.1,  0.1,  0.1,  0.1,  0.1 ]
sim.rl_0=   [ 1.01, 1.01, 1.01, 1.01, 1.01, 1.01, 1.01, 1.01]
sim.ru_0=   [ 3.0 , 3.0,  3.0,  3.0,  3.0,  3.0,  3.0,  3.0 ]
sim.initial_data_direction= "iiiiiiii"

## rescale for smaller amps
sim.amp_re= [4.0*x for x in sim.amp_re]
sim.amp_im= [4.0*x for x in sim.amp_im]

#sim.lin_m=  [-2,    2,  -5,    5   ]
#sim.l_ang=  [ 2,    2,   5,    5   ]
#sim.amp_re= [ 0.0,  0.1, 0.0,  0.1 ]
#sim.amp_im= [ 0.0,  0.0, 0.0,  0.0 ]
#sim.rl_0=   [ 1.01, 1.01,1.01, 1.01]
#sim.ru_0=   [ 3.0 , 3.0, 3.0 , 3.0 ]
#sim.initial_data_direction= "iiii"#
#=============================================================================
## second order modes to evolve

sim.scd_m= [-7,-6,-5,-4,-3,-2,-1,0, 1, 2, 3, 4, 5, 6, 7]
#sim.scd_m= [-7, -3, 3, 7]
#=============================================================================
## which m angular values to write to file

sim.write_lin_m= sim.lin_m 

sim.write_scd_m= sim.scd_m
#=============================================================================
## start multiple for second order metric evolution 

sim.start_multiple= float(6.0)
#=============================================================================
if (sim.run_type == "basic_run"):
   sim.launch_run()
#=============================================================================
elif (sim.run_type == "res_study"):
   all_res= [
      [160,30],
      [170,34],
      [180,38]
   ]
  
   for i in range(len(all_res)):
      sim.nx= all_res[i][0]
      sim.nl= all_res[i][1]
      sim.launch_run()
#=============================================================================
elif (sim.run_type == "amp_study"):
  
   sim.launch_run()

   for i in range(3):
      sim.amp_re= [2.0*x for x in sim.amp_re]
      sim.amp_im= [2.0*x for x in sim.amp_im]
      sim.launch_run()
#=============================================================================
elif (sim.run_type == "spin_ramp"):
   for bhs in [0,0.01,0.02,0.04,0.08,0.12,0.16,0.2,0.24,0.28,0.32]:
      sim.black_hole_spin= bhs
      sim.launch_run()
#=============================================================================
else:
   raise ValueError("run_type = "+str(sim.run_type)) 
