#import Atrium_class as AC
import Atrium2 as AC
#import Risk_Curve as RC
#import cProfile
#import pstats
#pr = cProfile.Profile()
"""Defaults: L=200,v=0.2,d=0.05,e=0.05,rp=50,tot_time = 10**6
                 ,pace_rate = 220,seed1 = 1,seed2=2,seed3=3"""
   
from line_profiler import LineProfiler            
A= AC.Atrium()        
lp = LineProfiler()
lp_wrapper = lp(A.CMP2D)
lp.add_function(A.CMP2D_timestep)
lp.add_function(A.SinusRhythm)
lp.add_function(A.Conduct)
lp.add_function(A.Relaxing)
lp_wrapper()
lp.print_stats()
#from line_profiler import LineProfiler            
#A= AC.Atrium()        
#lp = LineProfiler()
#lp_wrapper = lp(RC.RiskCurve(2))
#lp.add_function(A.CMP2D_timestep)
#lp.add_function(A.SinusRhythm)
#lp.add_function(A.Conduct)
#lp.add_function(A.Relaxing)
#lp_wrapper()
#lp.print_stats()
#pr = cProfile.Profile()
#pr.enable()
#A.CMP2D()
#print('Atrium')
#pr.disable()
#sortby = 'cumulative'
#ps = pstats.Stats(pr).sort_stats(sortby)
#ps.print_stats()