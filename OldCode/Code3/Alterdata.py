import matplotlib.pyplot as plt
import pickle 
import numpy as np
results = pickle.load(open("data_risk_curve_compare.p","rb"))
results2 = pickle.load(open("data_risk_curve_extra_data_0.16.p","rb"))
print(results[0])
print(len(results[1]))
print(len(results[2]))
#print(results2[0])
#results[0] = np.insert(results[0],5,results2[0][0])
#results[0] = np.insert(results[0],7,results2[0][1])
#results[0] = np.insert(results[0],32,results2[0][3])
#results[0] = np.insert(results[0],33,results2[0][4])
#results[0] = np.insert(results[0],34,results2[0][5])
#results[0] = np.insert(results[0],38,results2[0][1])
#results[1] = np.insert(results[1],38,results2[1][1])
#results[1] = np.delete(results[1],11)
#results[1] = np.insert(results[1],11,results2[1][0])
#results[1] = np.insert(results[1],24,results2[1][1])
#results[0] = np.insert(results[0],24,results2[0][1])
#results[1] = np.insert(results[1],33,results2[1][4])
#results[1] = np.insert(results[1],34,results2[1][5])
#results[1] = np.insert(results[1],35,results2[1][6])
#results[2] = np.insert(results[2],30,results2[2][1])
#results[2] = np.insert(results[2],31,results2[2][2])
#results[2] = np.insert(results[2],32,results2[2][3])
#results[2] = np.insert(results[2],33,results2[2][4])
#results[2] = np.insert(results[2],34,results2[2][5])
#results[2] = np.insert(results[2],35,results2[2][6])
#print(len(results[2]))
#results[2] = list(np.delete(results[1],37))
#a = results[2][:21]
#b = results[2][22:]
#a.extend([results2[2][0]])
#a.extend(b)
#results[2] = a
#a = list(results[2][:24])
#b = results[2][24:]
#a.extend([results2[2][1]])
#a.extend(b)
#results[2] = a
#
##a = results[2][:7]
#b = results[2][7:]
#a.extend([results2[2][1]])
#a.extend(b)
#results[2] = a
#print(results[2][30:])
#results[2] = np.delete(results[2],30)
#results[2] = np.delete(results[2],30)
#results[2] = np.delete(results[2],30)
#print(results[2][30:])
#print(results[2])
#print(len(results[2]))
#b = list(results[2])
#b = b.append(results2[2][1])
#results[2] = list(results[2]).extend(a)
#print(len(results[2]))
#print(len(results[1]))
#print(len(results[0]))
#data = [results[0],results[1],results[2]]
#pickle.dump(data,open( "data_risk_curve_compare.p", "wb" ) )
#plt.scatter(results[0],results[1])
