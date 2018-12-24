#import Atrium_new as AC

import Atrium_new_Jack as AC

"""Defaults: L=200,v=0.2,d=0.05,e=0.05,rp=50,tot_time = 10**6
                 ,pace_rate = 220,seed1 = 1,seed2=2,seed3=3"""
   
from line_profiler import LineProfiler            

A = AC.SourceSinkModel(hexagonal = True, tot_time = 10**3)
#A = AC.DysfuncModel(hexagonal = True, tot_time = 10**3)
lp = LineProfiler()
lp_wrapper = lp(A.cmp_full)
lp.add_function(A.cmp_timestep)
lp.add_function(A.sinus_rhythm)
lp.add_function(A.conduct)
lp.add_function(A.relaxing)
#lp.add_function(A.get_inward_current)
#lp.add_function(A.find_resting_neighbours)

lp_wrapper()
lp.print_stats()


#A = AC.Atrium(tot_time = 10**3)
#
#lp = LineProfiler()
#lp_wrapper = lp(A.CMP2D)
#lp.add_function(A.Conduct2)
#lp.add_function(A.CMP2D_timestep2)
#lp.add_function(A.Relaxing)
#lp.add_function(A.SinusRhythm2)
#
#
#lp_wrapper()
#lp.print_stats()