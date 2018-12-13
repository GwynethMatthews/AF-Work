<<<<<<< HEAD
#import Atrium_class as AC
import Atrium as AC
#import Risk_Curve as RC
#import cProfile
#import pstats
#pr = cProfile.Profile()
=======

import Atrium_new as AC

>>>>>>> Excitation time AF definition
"""Defaults: L=200,v=0.2,d=0.05,e=0.05,rp=50,tot_time = 10**6
                 ,pace_rate = 220,seed1 = 1,seed2=2,seed3=3"""
   
from line_profiler import LineProfiler            
<<<<<<< HEAD
A= AC.Atrium(v_tran_1 = 0.145, tot_time = 10**4, s2 = 49)        
lp = LineProfiler()
lp_wrapper = lp(A.CMP2D)
lp.add_function(A.CMP2D_timestep)
lp.add_function(A.SinusRhythm)
lp.add_function(A.Conduct)
lp.add_function(A.Relaxing)
lp.add_function(A.TimeInAF)
#lp.add_function(mode)
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
=======

A = AC.Atrium(hexagonal = True, model = 2, L = 100, v_para = 0.49,
                 v_tran_1 = 0.49, v_tran_2 = 0.49, d = 0.05,
                 threshold = 0.5, p = 0.75, e = 0.05, rp = 30, tot_time = 10**5, 
                 pace_rate = 220, s1 = 100, s2 = 4760, s3 = 3306, s4 = 476)

lp = LineProfiler()
lp_wrapper = lp(A.CMP2D)
lp.add_function(A.CMP2D_timestep)
#lp.add_function(A.SinusRhythm)
#lp.add_function(A.Conduct)
lp.add_function(A.Relaxing)
lp.add_function(A.TimeInAF)

lp_wrapper()
lp.print_stats()
>>>>>>> Excitation time AF definition
