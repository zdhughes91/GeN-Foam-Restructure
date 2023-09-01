### IMPORTS

import sys
import matplotlib.pyplot as plt

### FUNCTIONS

def nPIMPLE(lines) :
    
    n = 0;
    for line in lines :
        if "PIMPLE: iteration" in line :
            nn = int(line.split()[2])
            if nn > n :
                n = nn
            else :
                break
    return n

def plotEntryToAx(lines, key, pos, xAxis, plotType, lab, normalize=False, offset=0, col="") :

    entryValues = []
    if offset == 0 :
        for line in lines :
            if key in line :
                entry = line.split()[pos]
                entry = entry.split(")")[0]
                entry = entry.split("(")[0]
                entryValues.append(float(entry))
    else :
        for i in range(len(lines)) :
            if key in lines[i] :
                entry = lines[i+offset].split()[pos]
                entry = entry.split(")")[0]
                entry = entry.split("(")[0]
                entryValues.append(float(entry))

    if normalize :
        maxEntryValue = max(entryValues)
        for i in range(len(entryValues)) :
            entryValues[i]/=maxEntryValue
    lenX = len(xAxis)
    lenE = len(entryValues)
    newXAxis = []
    newEntryValues = []
    if lenX > lenE :
        for i in range(len(entryValues)) :
            newXAxis.append(xAxis[i])
        newEntryValues = entryValues
    elif lenX < lenE :
        for i in range(len(xAxis)) :
            newEntryValues.append(entryValues[i])
        newXAxis = xAxis
    else :
        newXAxis = xAxis
        newEntryValues = entryValues
    if col != "" :
        plotType(newXAxis, newEntryValues, label=lab, color=col)
    else :
        plotType(newXAxis, newEntryValues, label=lab)  

def axesGridLegend(axes) :
    for ax in axes :
        ax.grid(True)
        ax.legend(fontsize=6)
    axes[-1].set_xlabel("time (s)")

### MAIN

#fig, (ax1, ax2, ax3, ax4, ax5) = plt.subplots(5, sharex=True)
#axes = [ax1, ax2, ax3, ax4, ax5]
fig, (ax1, ax2, ax3, ax5) = plt.subplots(4, sharex=True)
axes = [ax1, ax2, ax3, ax5]
ax1.set_ylabel("mDot/mDot0 (-)")
ax1.set_ylim([-1, 4])
ax2.set_ylabel("p (Pa)")
ax2.set_ylim([8e4, 2e5])
ax3.set_ylabel("alpha (-)")
#ax4.set_ylabel("T (K)")
#ax4.set_ylim([600, 1500])
ax5.set_ylabel("Field extent (m)")

for q in range(len(sys.argv)-1) :

    # Read log lines
    filename = sys.argv[q+1]
    file = open(filename, "r")
    lines = file.readlines()

    # Read time, boiling inception time, power off time
    time = 0
    times = []
    tStartOfBoiling = 0
    tPowerOff = 0
    tEndOfBoiling= 0
    searchForStartOfBoiling = True
    searchForPowerOff = True
    searchForEndOfBoiling = True
    totalRunTime = 0
    for line in lines :
        if "Time =" in line and not "ExecutionTime" in line :
            time = float(line.split()[2])
            times.append(time)
        elif searchForStartOfBoiling :
            if "volAverage(fluidRegion) of alpha.vapour " in line :
                alphaAvg = float(line.split()[4])
                if alphaAvg > 0 :
                    tStartOfBoiling = time
                    searchForStartOfBoiling = False
        elif searchForPowerOff : 
            if "Power off" in line :
                tPowerOff = time
                searchForPowerOff = False
        elif searchForEndOfBoiling :
            if "dmdt liquid->vapour" in line :
                dmdtMax = float(line.split()[8])
                if dmdtMax == 0 :
                    tEndOfBoiling = time
                    searchForEndOfBoiling = False
        if "ExecutionTime = " in line :
            totalRunTime = float(line.split()[2])

    '''                
    t0 = times[0]
    for i in range(len(times)) :
        times[i] = times[i] - t0  
    tPowerOff -= t0
    tStartOfBoiling -= t0
    tEndOfBoiling -= t0
    '''

    # Plot stuff
    print(filename)
    print(" -> tStartOfBoiling = "+str(tStartOfBoiling))
    print(" -> tPowerOff = "+str(tPowerOff))
    print(" -> tEndOfBoiling = "+str(tEndOfBoiling))
    print(" -> totalRunTime = "+str(totalRunTime/3600.0)+" h")
    plotEntryToAx(lines, "average(inlet) of U.liquid ", 6, times, ax1.plot, filename)
    #plotEntryToAx(lines, "volAverage(Z545) of p ", 4, times, ax2.plot, "p.Z545."+filename)
    plotEntryToAx(lines, "volAverage(Z779) of p ", 4, times, ax2.plot, "p.Z779."+filename)
    plotEntryToAx(lines, "volAverage(Z545) of alpha.vapour ", 4, times, ax3.plot, "alpha.Z545."+filename) 
    plotEntryToAx(lines, "volAverage(Z779) of alpha.vapour ", 4, times, ax3.plot, "alpha.Z779."+filename)
    plotEntryToAx(lines, "volAverage(fluidRegion) of alpha.vapour ", 4, times, ax3.plot, "alpha.avg."+filename)
    #plotEntryToAx(lines, "max(fluidRegion) of T.activeStructure ", 4, times, ax4.plot, "T.structure.max"+filename)
    #plotEntryToAx(lines, "max(fluidRegion) of T.vapour ", 4, times, ax4.plot, "T.vapour.max."+filename)
    #plotEntryToAx(lines, "max(fluidRegion) of T.liquid ", 4, times, ax4.plot, "T.liquid.max."+filename)
    #plotEntryToAx(lines, "max(fluidRegion) of T.interface ", 4, times, ax4.plot, "T.interface.max."+filename)
    plotEntryToAx(lines, "field: dmdt", 4, times, ax5.plot, "dmdt."+filename, False, 1)
    plotEntryToAx(lines, "field: dmdt", 7, times, ax5.plot, "", False, 1)
    plotEntryToAx(lines, "field: alpha.vapour", 4, times, ax5.plot, "alpha."+filename, False, 1)
    plotEntryToAx(lines, "field: alpha.vapour", 7, times, ax5.plot, "", False, 1)

# Read experimental data
eFilename = "experimentalResults/inletMassFlow.txt"
eFile = open(eFilename, "r")
eLines = eFile.readlines()
eTime = []
eFlow = []
u0 = 3.38
for line in eLines :
    eTime.append(float(line.split(",")[0])+3)
    eFlow.append(u0*float(line.split(",")[1]))
ax1.plot(eTime, eFlow, label="exp")     

axesGridLegend(axes)
plt.show()
