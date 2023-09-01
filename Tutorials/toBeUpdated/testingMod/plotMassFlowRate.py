### IMPORTS

import sys
import matplotlib.pyplot as plt

### MAIN





for q in range(len(sys.argv)-1) :

    filename = sys.argv[q+1]
    file = open(filename, "r")
    lines = file.readlines()
    
    times = []
    mFlowRates = [0]
    i = 0
   
    time = 0
    
    for line in lines :
        
        if "Time =" in line and not "ExecutionTime" in line :
            time = float(line.split()[2])
        if "faceZone massFlowSurface_z1 massFlow" in line :
            mFlowRate = float(line.split()[4])
            
            mFlowRates.append(mFlowRate)
            times.append(time)
           
    
    
    t0 = times[0]
    for i in range(len(times)) :
    	times[i] = times[i] - t0  
    del mFlowRates[-1]
 
 

plt.figure()
plt.plot(times, mFlowRates)
plt.xlabel("iteration")
plt.ylabel("Mass flow rate")
plt.grid()
plt.show()
