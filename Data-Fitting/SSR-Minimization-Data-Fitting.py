import numpy as np
import matplotlib.pyplot as plt
import os
import gc
import pathlib
from scipy.optimize import minimize

SCTime = [2, 3, 4, 5, 6, 7, 8, 9, 10, 12, 15, 18]
SC = [ 0.00007109709432312438, 0.0002586416205275976, 0.0019073187599725216, 0.039137456019803916, 0.048733846139805366, 0.16885472875899837, 0.2618119214401843, 0.31815994481219106, 0.3961718357595989, 0.2493592004984161, 0.4262158829015338, 1.2152232910634688]

MCTime = [12, 18, 24, 30, 36, 48, 60, 72, 84, 96, 108, 120]
MC = np.array([819.4679979024306, 7759.590246392491, 2905620.46821282, 73765857.98521078, 999437728.9462978, 6640163884.645269, 6709615082.419942, 5663087196.46287, 1040447862.7422994, 1051330177.0262353, 54967688.493010186, 1403746.5929611742])

#MCError = np.empty([len(MC)])
#for i in range(len(MC)):   
#    power = np.around(MC[i],-1*int(np.floor(np.log10(MC[i]))+1))
#    if power == 0.0:
#        MCError[i] = np.power(10.0,int(np.floor(np.log10(MC[i]))))
#    else:
#        MCError[i] = np.around(MC[i],-1*int(np.floor(np.log10(MC[i]))+1))
MCError = MC*0.5

#plt.plot(MCTime, MC, "o")
#plt.yscale("log")
#plt.show()

def SSR(x):

    beta = x[0]
    rho = x[1]
    TauI = x[2]
    TauE = x[3]
    c = x[4]
    
    print(beta, rho, TauI, TauE, c)
        
    MedianData = []
    MedianCD = []

    LongestTime = 0

    os.system("nvcc DataFitting.cu -o program.out && ./program.out"+" "+str(beta)+" "+str(rho)+" "+str(TauI)+" "+str(TauE)+" "+str(c))
    gc.collect()

    with open("Run:"+str(beta)+","+str(rho)+","+str(TauI)+","+str(TauE)+","+str(c)+".txt") as infile:
        index = 0
        for line in infile:
            Data = line
            Data = Data.split(",")
            del Data[-1]
            VirusData[index] = np.asarray(Data, dtype=np.float)
            if(LongestTime < len(VirusData[index])):
                LongestTime = len(VirusData[index])
            index = index + 1
            if index == NumberOfRuns:
                break

    for j in range(LongestTime):
        InbetweenArray = []
        for i in range(NumberOfRuns):
            try:
                InbetweenArray.append(VirusData[i][j])
            except IndexError:
                InbetweenArray.append(0.0)
            
        MedianData.append(np.median(InbetweenArray))
        
    for i in MCTime:
        try:
            MedianCD.append(MedianData[i])
        except IndexError:
            MedianCD.append(0.0)
        
    SSR = 0
    for i in range(len(MCTime)):
        SSR = SSR + (np.log10(MC[i])-np.log10(MedianCD[i]))**2

    gc.collect()

    return SSR

NumberOfRuns = 10

VirusData = np.empty([NumberOfRuns],dtype=tuple )

beta = 54.21
rho = 3000.0
TauI = 26.24
TauE = 15.16
c = 0.25

#beta = 50.0
#rho = 200.0
#TauI = 50.0
#TauE = 15.0

x0 = [beta, rho, TauI, TauE, c]

res = minimize(SSR, x0, method="Nelder-Mead", options={'xatol':1e-4})

print(res)
print("\n")
print("DONE")
