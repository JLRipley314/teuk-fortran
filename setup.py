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
sim.black_hole_spin= round(0.999*sim.black_hole_mass,6)
sim.compactification_length= float(1)
#=============================================================================
sim.evolve_time= float(200) ## units of black hole mass
sim.num_saved_times= int(600)
#=============================================================================
sim.nx= 180 ## num radial pts 
sim.nl= 34 ## num angular values
#=============================================================================
sim.metric_recon=     True
sim.scd_order=        True
sim.constrained_evo = True

sim.write_indep_res=           True
sim.write_metric_recon_fields= False
sim.write_scd_order_source=    True
sim.write_coefs=               False
#=============================================================================
sim.computer= 'della'#'home'#
sim.della_out_stem= '/tigress/jripley/tf-out/'

## for della cluster/slurm script

sim.walltime= '36:00:00' ## (hh:mm:ss)
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
sim.lin_m= [-2, 2, -3, 3]
#-----------------------------------------------------------------------------
sim.l_ang= [ 2, 2,  3, 3]

sim.initial_data_direction= "iiii"

sim.amp_re= [0.0, 4.0, 0.0, 2.0]
sim.amp_im= [0.0, 0.0, 0.0, 0.0]

sim.rl_0= [1.1, 1.1, 1.1, 1.1]
sim.ru_0= [3.0, 3.0, 2.5, 2.5]
#=============================================================================
## second order modes to evolve

sim.scd_m= [0, -1, 1]
#=============================================================================
## which m angular values to write to file

sim.write_lin_m= sim.lin_m 

sim.write_scd_m= sim.scd_m
#=============================================================================
## start multiple for second order metric evolution 

sim.start_multiple= float(8.0)
#=============================================================================
if (sim.run_type == "basic_run"):
   sim.launch_run()
#=============================================================================
elif (sim.run_type == "multiple_runs"):
#-----------------------------------------------------------------------------
   default_sm  = 1

   default_nx_07 = 176
   default_nl_07 = 32
#   default_nx_07 = 144
#   default_nl_07 = 24

   default_nx_0998 = 214
   default_nl_0998 = 44
#   default_nx_0998 = 160
#   default_nl_0998 = 28
#-----------------------------------------------------------------------------
#-----------------------------------------------------------------------------
   sim.black_hole_spin= round(0.7*sim.black_hole_mass,6)
   sim.start_multiple = default_sm
   sim.nx = default_nx_07
   sim.nl = default_nl_07

   nxs = [160, 176, 192]
   nls = [ 28,  32,  36]
#   nxs = [128, 144, 160]
#   nls = [ 20,  24,  28]
   for i in range(len(nxs)):
      sim.nx = nxs[i]
      sim.nl = nls[i]
      sim.launch_run() 
#-----------------------------------------------------------------------------
   sim.black_hole_spin= round(0.7*sim.black_hole_mass,6)
   sim.start_multiple = default_sm
   sim.nx = default_nx_07
   sim.nl = default_nl_07

   sms = [2,3]
   for sm in sms:
      sim.start_multiple= sm
      sim.launch_run()
#-----------------------------------------------------------------------------
   sys.exit()
#-----------------------------------------------------------------------------
   sim.black_hole_spin= round(0.99998*sim.black_hole_mass,6)
   sim.start_multiple = default_sm
   sim.nx = default_nx_0998
   sim.nl = default_nl_0998

   nxs = [208, 214, 230]
   nls = [ 40,  44,  48]
#   nxs = [144, 160, 176]
#   nls = [ 24,  28,  32]
   for i in range(len(nxs)):
      sim.nx = nxs[i]
      sim.nl = nls[i]
      sim.launch_run() 
#-----------------------------------------------------------------------------
   sim.black_hole_spin= round(0.99998*sim.black_hole_mass,6)
   sim.start_multiple = default_sm
   sim.nx = default_nx_0998
   sim.nl = default_nl_0998

   sms = [2,3]
   for sm in sms:
      sim.start_multiple= sm
      sim.launch_run()

#=============================================================================
elif (sim.run_type == "spin_ramp"):
   for bhs in [0,0.01,0.02,0.04,0.08,0.12,0.16,0.2,0.24,0.28,0.32]:
      sim.black_hole_spin= bhs
      sim.launch_run()
#=============================================================================
else:
   raise ValueError("run_type = "+str(sim.run_type)) 
