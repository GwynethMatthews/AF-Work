import numpy as np
import pickle
import matplotlib.pyplot as plt

results1 = pickle.load(open('RiskCurveTrialforPoster1.p',"rb"))
results2 = pickle.load(open('RiskCurveTrialforPoster2.p',"rb"))
results3 = pickle.load(open('RiskCurveTrialforPoster3.p',"rb"))
results4 = pickle.load(open('RiskCurveTrialforPoster4.p',"rb"))

plt.scatter(results1[0],results1[1])
plt.scatter(results2[0],results2[1])
plt.scatter(results3[0],results3[1])
