import matplotlib.pyplot as plt
import pickle 

results = pickle.load(open("riskphasedata.p", "rb"))

levels = [0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1]
#levels = [0,0.1,0.9,1]

plt.contourf(results[1],results[2],results[3],levels)
#plt.pcolormesh(results[1],results[2],results[3])

plt.colorbar()
plt.xlabel('v_para')
plt.ylabel('v_tran')    
