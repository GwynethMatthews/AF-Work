
import Atrium_new_Jack as AC

"""Defaults: L=200,v=0.2,d=0.05,e=0.05,rp=50,tot_time = 10**6
                 ,pace_rate = 220,seed1 = 1,seed2=2,seed3=3"""
   
from line_profiler import LineProfiler            

A = AC.SourceSinkModel(tot_time = 10**3)

lp = LineProfiler()
lp_wrapper = lp(A.cmp_full)
lp.add_function(A.cmp_timestep)
lp.add_function(A.sinus_rhythm)
lp.add_function(A.conduct)
lp.add_function(A.relaxing)
lp.add_function(A.get_inward_current)

lp_wrapper()
lp.print_stats()
