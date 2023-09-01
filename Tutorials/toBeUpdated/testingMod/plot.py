import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import sys


for q in range(len(sys.argv)-1) :

    filename = sys.argv[q+1]
    file = open(filename, "r")
    lines = file.readlines()
    
    timesLog = [0]
    mFlowRates = [0]
    i = 0
   
    time = 0
    
    for line in lines :
        
        if "Time =" in line and not "ExecutionTime" in line :
            time = float(line.split()[2])
        if "faceZone massFlowSurface_z1 massFlow" in line :
            mFlowRate = float(line.split()[4])
            
            mFlowRates.append(mFlowRate)
            timesLog.append(time)
           

    
    
    for i in range(len(timesLog)) :
    	timesLog[i] = timesLog[i] - timesLog[0]
    del mFlowRates[0]
    del timesLog[0]


mFlowRates = np.array(mFlowRates)


mFlowRatesScaled = (mFlowRates - mFlowRates[0])/mFlowRates[0]

data = pd.read_csv("momentumSourceTest.csv")


time = data["time"]
momentumSource = data["model.root.system1.modelicaMomentumSource"]
power = data["model.root.system1.power"]
desiredPower = data["model.root.system1.pi.SP"]

desiredPowerScaled = data["model.root.system1.pi.SPs"]
powerScaled = data["model.root.system1.pi.PVs"]

momentumSourceScaled = data["model.root.system1.pi.CSs"]





plt.plot(time, powerScaled, label='normalized power variation')
plt.plot(time, desiredPowerScaled, '--', label='normalized desired power')
#plt.plot(time, momentumSourceScaled, label='normalized momentumSource variation')
plt.plot(timesLog, mFlowRatesScaled, label='normalized momentumSource variation')
plt.xlabel("Time (s)")
plt.ylabel(" normalized variation")
plt.legend()
plt.show()
