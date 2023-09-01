### IMPORTS

import copy
import sys
import math
import matplotlib.pyplot as plt

### SETTINGS

beta = 1e5*0.00313126

SMALL_SIZE = 10
MEDIUM_SIZE = 12
BIGGER_SIZE = 14

plt.rc('font', size=SMALL_SIZE)          # controls default text sizes
plt.rc('axes', titlesize=SMALL_SIZE)     # fontsize of the axes title
plt.rc('axes', labelsize=SMALL_SIZE)    # fontsize of the x and y labels
plt.rc('xtick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
plt.rc('ytick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
plt.rc('legend', fontsize=SMALL_SIZE-1)    # legend fontsize
plt.rc('figure', titlesize=BIGGER_SIZE)  # fontsize of the figure title

### FUNCTIONS

def readTimes(loglines) :

    times = []
    totalRunTime = 0
    for line in loglines :
        if "Time =" in line and not "ExecutionTime" in line :
            times.append(float(line.split()[2]))
        if "ExecutionTime = " in line :
            totalRunTime = float(line.split()[2])
    print("    total run time = "+str(totalRunTime)+" s")
    return times

def readValues(times, loglines, keyword, pos, scale=1) :

    values = []
    for line in loglines :
        if keyword in line:
            values.append(scale*float(line.split()[pos]))
            if len(values) == len(times) :
                return values
    return values

def readExpTimesAndValues(loglines, timeOffset, valueOffset, scaleFactor) :
    times = []
    values = []
    for line in loglines :
        try :
            times.append(float(line.split()[0])+timeOffset)
            values.append(scaleFactor*float(line.split()[1])+valueOffset)
        except :
            pass
    return times, values

def matchSize(list1, list2):
    
    while len(list1) > len(list2) :
        del list1[-1]
    while len(list2) > len(list1) :
        del list2[-1]

###

fig, axes = plt.subplots(4, 1, sharex=True)

axP, axF, axT, axR = axes.flatten()

for ax in axes.flatten(): 
    ax.grid(True)

axP.set_ylabel(r'$P\;(W)$')
axF.set_ylabel(r'$\dot{m}\;(kg/s)$')
axT.set_ylabel(r'$T\;(K)$')
axR.set_ylabel(r'$\rho\;(\$)$')
axR.set_xlabel(r'$t\;(s)$')

for q in range(len(sys.argv)-1) :

    logname = sys.argv[q+1]
    print("Processing "+logname)
    log = open(logname, "r")
    loglines = log.readlines()

    times = readTimes(loglines)

    for i in range(10) :
        del times[-1]

    flowPrimary = readValues(times, loglines, "faceZone pumpMiddleCutPrimary massFlow", 4)
    flowCore = readValues(times, loglines, "faceZone coreMiddleCut massFlow", 4)
    flowInnerCore = readValues(times, loglines, "faceZone innerCoreMiddleCut massFlow", 4)
    flowOuterCore = readValues(times, loglines, "faceZone outerCoreMiddleCut massFlow", 4)
    flowColdLegSecondaryInlet = readValues(times, loglines, "patch coldLegSecondaryInlet massFlow", 4)


    axF.plot(times, flowPrimary, label="primary")
    axF.plot(times, flowCore, label="core")
    axF.plot(times, flowInnerCore, label="innerCore")
    axF.plot(times, flowOuterCore, label="outerCore")

    TCoreIn = readValues(times, loglines, "faceZone coreInlet TBulk", 4)
    TCoreOut = readValues(times, loglines, "faceZone coreOutlet TBulk", 4)
    TInnerCoreOut = readValues(times, loglines, "faceZone innerCoreOutlet TBulk", 4)
    TOuterCoreOut = readValues(times, loglines, "faceZone outerCoreOutlet TBulk", 4)
    TIHXInPrimary = readValues(times, loglines, "faceZone primaryIHXInlet TBulk", 4)
    TIHXOutPrimary = readValues(times, loglines, "faceZone primaryIHXOutlet TBulk", 4)
    TIHXInSecondary = readValues(times, loglines, "faceZone secondaryIHXInlet TBulk", 4)
    TIHXOutSecondary = readValues(times, loglines, "faceZone secondaryIHXOutlet TBulk", 4)
    THotLegPrimary = readValues(times, loglines, "faceZone hotLegMiddleCut TBulk", 4)
    TColdLegPrimary = readValues(times, loglines, "faceZone coldLegMiddleCut TBulk", 4)


    axT.plot(times, TCoreIn, label="coreInlet")
    axT.plot(times, TCoreOut, label="coreOutlet")
    axT.plot(times, TInnerCoreOut, label="innerCoreOutlet")
    axT.plot(times, TOuterCoreOut, label="outerCoreOutlet")
    axT.plot(times, THotLegPrimary, label="THotLegPrimary")
    axT.plot(times, TColdLegPrimary, label="TColdLegPrimary")
    #axT.plot(times, TIHXInPrimary, label="IHXPrimaryInlet")
    #axT.plot(times, TIHXOutPrimary, label="IHXPrimaryOutlet")
    #axT.plot(times, TIHXInSecondary, label="IHXSecondaryInlet")
    #axT.plot(times, TIHXOutSecondary, label="IHXSecondaryOutlet")

    try:
        totalPower = readValues(times, loglines, "totalPower =", 2, 180)
        fissionPower = readValues(times, loglines, "-> fission", 3, 180)
        decayPower = readValues(times, loglines, "-> decay", 3, 180)
        axP.plot(times, totalPower, label="totalPower")
        axP.plot(times, fissionPower, label="fissionPower")
        axP.plot(times, decayPower, label="decayPower")
    except :
        pass

    try:
        RTot = readValues(times, loglines, "totalReactivity", 2, 1.0/beta)
        RDoppler = readValues(times, loglines, "Doppler (fast)", 4, 1.0/beta)
        RTFuel = readValues(times, loglines, "-> TFuel", 3, 1.0/beta)
        RTClad = readValues(times, loglines, "-> TClad", 3, 1.0/beta)
        RRhoCool = readValues(times, loglines, "-> rhoCool", 3, 1.0/beta)
        RTStruct = readValues(times, loglines, "-> TStruct", 3, 1.0/beta)
        RDriveline = readValues(times, loglines, "-> driveline", 3, 1.0/beta)
        RGEM = readValues(times, loglines, "-> GEM", 3, 1.0/beta)
        axR.plot(times, RTot, label="total")
        axR.plot(times, RDoppler, label="Doppler")
        axR.plot(times, RTFuel, label="Fuel axial expansion")
        axR.plot(times, RTClad, label="Cladding temperature")
        axR.plot(times, RRhoCool, label="Coolant density")
        axR.plot(times, RTStruct, label="Core radial expansion")
        axR.plot(times, RDriveline, label="Driveline")
        axR.plot(times, RGEM, label="GEM")
    except :
        pass

###

timeOffset = 910

log = open("expFlow", "r")
loglines = log.readlines()
expFlowTimes, expFlowValues = readExpTimesAndValues(loglines, timeOffset, 0, 1)
axF.plot(expFlowTimes, expFlowValues, label="expTotal")

log = open("expPIOTA2", "r")
loglines = log.readlines()
expPIOTA2Times, expPIOTA2Values = readExpTimesAndValues(loglines, timeOffset, 273.15, 1)
axT.plot(expPIOTA2Times, expPIOTA2Values, label="expPIOTA2")

'''
log = open("expPIOTA6", "r")
loglines = log.readlines()
expPIOTA6Times, expPIOTA6Values = readExpTimesAndValues(loglines, timeOffset, 273.15, 1)
axT.plot(expPIOTA6Times, expPIOTA6Values, label="expPIOTA6")
'''

log = open("expPow", "r")
loglines = log.readlines()
expPowTimes, expPowValues = readExpTimesAndValues(loglines, timeOffset, 0, 1000000)
axP.plot(expPowTimes, expPowValues, label="expPrimary")

log = open("expReactivity", "r")
loglines = log.readlines()
expRTimes, expRValues = readExpTimesAndValues(loglines, timeOffset, 0, 1)
axR.plot(expRTimes, expRValues, label="expTotal")

for ax in axes.flatten() :
    ax.legend()

plt.show()
