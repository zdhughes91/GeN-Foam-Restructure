###############################################################################
### IMPORTS
###############################################################################

import os
import re
import sys
import math
import time
import string

###############################################################################
### GLOBAL
###############################################################################

global root
global gPointList
global gFaceList
global gCellList
debug = False

root = os.getcwd()
gPointList = []
gFaceList = []
gCellList = []

###############################################################################
### CLASSES
###############################################################################

class pointClass :

    def __init__(self, cmpts, pointLabel) :

        global gPointList
        
        if len(cmpts) != 3 :
            sys.exit("Invalid point size")
        for i in range(0, len(cmpts)) :
            cmpts[i] = float(cmpts[i])
        self.cmpts = cmpts
        self.label = pointLabel
        gPointList.append(self)

    def __eq__(self, rhs) :
        if self.label == rhs.label :
            return True
        else :
            return False

class faceClass :

    def __init__(self, pointLabelList, faceLabel) :
        
        global gPointList
        global gFaceList
        self.pointList = []

        for pointLabel in pointLabelList :
            self.pointList.append(gPointList[pointLabel])
        self.label = faceLabel
        gFaceList.append(self)

    def center(self) :

        x = 0
        y = 0
        z = 0
        N = len(self.pointList)

        for point in self.pointList :
            x += point.cmpts[0]/N
            y += point.cmpts[1]/N
            z += point.cmpts[2]/N
        self.C = [x, y, z]
        return [x, y, z]

class cellClass :

    def __init__(self, cellLabel) :

        global gCellList
        self.faceList = []
        self.label = cellLabel
        self.zone = ""

        gCellList.append(self)

    def __eq__(self, rhs) :
        if self.label == rhs.label :
            return True
        else :
            return False

    def center(self) :
        
        x = 0
        y = 0
        z = 0
        N = len(self.faceList)

        for face in self.faceList :
            x += face.center()[0]/N
            y += face.center()[1]/N
            z += face.center()[2]/N
        self.C = [x, y, z]
        return [x, y, z]

    def setDelta(self) :

        points = []
        for face in self.faceList :
            for point in face.pointList :
                points.append(point)
        cmptsx = [point.cmpts[0] for point in points]
        cmptsy = [point.cmpts[1] for point in points]
        cmptsz = [point.cmpts[2] for point in points]
        dx = max(cmptsx)-min(cmptsx)
        dy = max(cmptsy)-min(cmptsy)
        dz = max(cmptsz)-min(cmptsz)
        self.delta = [dz, dy, dz]

    def volume(self) : ### SPECIFIC FOR CURRENT MESH TYPE (i.e. trapzeio)
        
        H1 = self.delta[0]
        B = self.delta[1]
        H2 = self.delta[2]
        points = []

        for face in self.faceList :
            for point in face.pointList :
                if point not in points :
                    points.append(point)
        cmptsx = [point.cmpts[0] for point in points]
        cmptsz = [point.cmpts[2] for point in points]
        minxPoints = [point for point in points if ((point.cmpts[0] == min(cmptsx)) and (point.cmpts[2] == min(cmptsz)))]
        if len(minxPoints) == 1 :
            b = 0
        else :
            b = abs(minxPoints[0].cmpts[1]-minxPoints[1].cmpts[1])
        V = (B + b)*H1*H2/2
        self.V = V
        return V

class scalarFieldClass :

    def __init__(self, field, name, units, time) :
        self.field = field
        self.name = name
        self.units = units
        self.time = time

class vectorFieldClass :

    def __init__(self, field, name, units, time) :
        self.field = field
        self.name = name
        self.units = units
        self.time = time

    def mag(self) :
        field = []
        for i in range(0, len(self.field)) :
            field.append((self.field[i][0]**2 + self.field[i][1]**2 + self.field[i][2]**2)**0.5)
        return scalarFieldClass(field, "mag"+self.name, self.units, time)

    def cmpt(self, cmpti) :
        field = []
        if len(self.field) != 1 :
            for i in range(0, len(self.field)) :
                field.append(self.field[i][cmpti])
        else:
            field.append(self.field[0][cmpti])
        return scalarFieldClass(field, self.name+"["+str(cmpti)+"]", self.units, time)

###############################################################################
### FUNCTIONS
###############################################################################

### UTILITY

def valuesInBrakets(strList, mode) :

    valList = []

    inBrakets = False
    for strg in strList :
        tmpStrg = ""
        for char in strg :
            if char == "(" :
                inBrakets = True
            elif char == ")" :
                inBrakets = False
            elif inBrakets :
                tmpStrg+=char
        if mode == "int" and tmpStrg != "" :
            valList.append(int(float(tmpStrg)))
        if mode == "float" and tmpStrg != "" :
            valList.append(float(tmpStrg))
    return valList

def readFile(filename) :
    
    lines = []

    tmpLines = open(filename, "r").readlines()
    for line in tmpLines :
        line = re.split("\s", line)
        while '' in line:
            line.remove('')
        if len(line) > 0 :
            lines.append(line)
    return lines
    tmpLines.close()

def getTimes(root) :
    times = []
    subDirs = [name for name in os.listdir(root) if os.path.isdir(os.path.join(root, name))]
    for i in range(0, len(subDirs)):
        try:
            float(subDirs[i])
            times.append(subDirs[i])
        except:
            pass
    return times

def printTable(axis, axisLabel, values, valuesLabel, valuesSI) :
    if isinstance(axis, list) :
        print(axisLabel+" [m]"+"    "+valuesLabel+" ["+valuesSI+"]")
        for i in range(0, len(axis)) :
            print("{:.3F}".format(axis[i]), "   ", "{:.3E}".format(values[i]))
    print("\n")

### PLOYMESH READING FUNCTIONS

def readPoints() :

    global gPointList
    global gFaceList
    global gCellList
    
    file = readFile("points")
    for i in range(0, len(file)) :
        line = file[i]
        if len(line) == 1 and line[0] == "(" :
            while True :
                i += 1
                line = file[i]
                if len(line) != 1 :
                    point = pointClass(valuesInBrakets(line, "float"), len(gPointList))
                else :
                    break

def readFaces() :

    global gPointList
    global gFaceList
    global gCellList

    file = readFile("faces")
    for i in range(0, len(file)) :
        line = file[i]
        if len(line) == 1 and line[0] == "(" :
            while True :
                i += 1
                line = file[i]
                if len(line) != 1 :
                    face = faceClass(valuesInBrakets(line, "int"), len(gFaceList))
                else :
                    break
    for face in gFaceList :
        face.center()

def readCells() :

    global gPointList
    global gFaceList
    global gCellList

    file = readFile("owner")
    for i in range(0, len(file)) :
        line = file[i]
        if len(line) == 1 and line[0] == "(" :
            faceLabel = 0
            while True :
                i += 1
                if file[i][0] != ")" :
                    cellLabel = int(float(file[i][0]))
                    if cellLabel not in [cell.label for cell in gCellList] :
                        cell = cellClass(cellLabel)
                    else :
                        cell = [cell for cell in gCellList if cell.label == cellLabel][0]
                    cell.faceList.append(gFaceList[faceLabel])
                    faceLabel += 1
                else :
                    break
    file = readFile("neighbour")
    for i in range(0, len(file)) :
        line = file[i]
        if len(line) == 1 and line[0] == "(" :
            faceLabel = 0
            while True :
                i += 1
                if file[i][0] != ")" :
                    cellLabel = int(float(file[i][0]))
                    cell = [cell for cell in gCellList if cell.label == cellLabel][0]
                    cell.faceList.append(gFaceList[faceLabel])
                    faceLabel += 1
                else :
                    break
    file = readFile("cellZones")
    startOfData = False
    cellZone = ""
    startOfCellZoneData = False
    for i in range(0, len(file)) :
        line = file[i]
        if line[0] == "(" and not startOfData :
            startOfData = True
        elif startOfData :
            if line[0] == "{" :
                cellZone = file[i-1][0]
            elif cellZone != "" and line[0] == "(" :
                startOfCellZoneData = True
            elif startOfCellZoneData : 
                if line[0] != ")":
                    cellLabel = int(float(line[0]))
                    cell = [cell for cell in gCellList if cell.label == cellLabel][0]
                    cell.zone = cellZone
                else :
                    startOfCellZoneData = False
                    cellZone = ""
    for cell in gCellList :
        cell.center()
        cell.setDelta()
        cell.volume()

def domainSize() :
    global gCellList
    Cx = []
    Cy = []
    Cz = []
    dx = min([cell.delta[0] for cell in gCellList])/2
    dy = min([cell.delta[1] for cell in gCellList])/2
    dz = min([cell.delta[2] for cell in gCellList])/2
    for cell in gCellList :
        Cx.append(cell.C[0])
        Cy.append(cell.C[1])
        Cz.append(cell.C[2])
    return min(Cx)-dx, min(Cy)-dy, min(Cz)-dz, max(Cx)+dx, max(Cy)+dy, max(Cz)+dz

### FIELD LOADER

# Load and return a single field for a given time
def loadField(name, units, time, region, mode) :
    path = root+"/"+time+"/"+region
    os.chdir(path)
    file = readFile(name)
    field = []
    stopFlag = False
    for i in range(0, len(file)) :
        if stopFlag != True :
            line = file[i]
            if line[0] == "internalField" and line[1] == "uniform" :
                if mode == "scalar" :
                    field.append(float(line[2].replace(";", "")))
                    stopFlag = True
                elif mode == "vector" :
                    field.append(valuesInBrakets(line, "float"))
                    stopFlag = True
            elif line[0] == "(" :
                while True :
                    i += 1
                    line = file[i]
                    if line[0] == ")" :
                        stopFlag = True # be sure not to read anything further
                        break
                    if mode == "scalar" :
                        field.append(float(line[0]))
                    elif mode == "vector" :
                        field.append(valuesInBrakets(line, "float"))
    if mode == "scalar" :
        return scalarFieldClass(field, name, units, time)
    if mode == "vector" :
        return vectorFieldClass(field, name, units, time)

### DATA SAMPLER

def cellFromCoordinates(x0, y0, z0) :

    global gCellList

    for cell in gCellList :
        dx = cell.delta[0]/2
        dy = cell.delta[1]/2
        dz = cell.delta[2]/2
        withinX = ( (x0-dx) <= cell.C[0] <= (x0+dx) )
        withinY = ( (y0-dy) <= cell.C[1] <= (y0+dy) )
        withinZ = ( (z0-dz) <= cell.C[2] <= (z0+dz) ) 
        if withinX and withinY and withinZ :
            return cell

def sampleCellsAlongAxis(x0, y0, z0, axis, zones) :
    
    global gCellList
    dx = min([cell.delta[0] for cell in gCellList])
    dy = min([cell.delta[1] for cell in gCellList])
    dz = min([cell.delta[2] for cell in gCellList])
    cells = []
    x = x0
    y = y0
    z = z0
    cx = 0
    cy = 0
    cz = 0
    minx, miny, minz, maxx, maxy, maxz = domainSize()

    if axis == "x" :
        cx = 1
    if axis == "y" :
        cy = 1
    if axis == "z" :
        cz = 1
    while True :
        cell = cellFromCoordinates(x, y, z)
        if (minx <= x and x <= maxx) and (miny <= y and y <= maxy) and (minz <= z and z <= maxz) :
            if cell is not None :
                if cell not in cells :
                    if zones[0] == "All" :
                        cells.append(cell)
                    else :
                        if cell.zone in zones :
                            cells.append(cell)
        else :
            break
        x += cx*dx
        y += cy*dy
        z += cz*dz

    return cells

# intgAx = integration axis, profAx = profileAxis
def integrationParameters(zones, intgAx, profAx) :

    # List cell coordinates in more useful manner
    global gCellList
    Cs = [] # vector of cell centers
    if zones[0] != "All" :
        cellsInZone = [cell for cell in gCellList if cell.zone in zones]
    else :
        cellsInZone = gCellList
    for cell in cellsInZone :
        Cs.append(cell.C)
    x0 = min([C[0] for C in Cs])
    y0 = min([C[1] for C in Cs])
    z0 = min([C[2] for C in Cs])
    dx = min([cell.delta[0] for cell in gCellList])
    dy = min([cell.delta[1] for cell in gCellList])
    dz = min([cell.delta[2] for cell in gCellList])

    # Define max of profile range
    if profAx == "x" :
        profI = 0
        delta = dx
        profCoord = x0
    elif profAx == "y" :
        profI = 1
        delta = dy 
        profCoord = y0
    elif profAx == "z" :
        profI = 2
        delta = dz
        profCoord = z0
    maxProfCoord = max([C[profI] for C in Cs]) + delta

    # Define direction index of integration range
    if intgAx == "x" :
        intgI = 0
    elif intgAx == "y" :
        intgI = 1
    elif intgAx == "z" :
        intgI = 2
    
    profRows = [] # profile
    x = x0
    y = y0
    z = z0
    while profCoord < maxProfCoord :
        if profAx == "x" :
            intgCells = sampleCellsAlongAxis(profCoord, y, z, intgAx, zones)
            #print("x = ", profCoord)
        elif profAx == "y" :
            intgCells = sampleCellsAlongAxis(x, profCoord, z, intgAx, zones)
            #print("y = ", profCoord)
        elif profAx == "z" :
            intgCells = sampleCellsAlongAxis(x, y, profCoord, intgAx, zones)
            #print("z = ", profCoord)

        if len(intgCells) > 0 and intgCells not in profRows:
            profRows.append(intgCells)
        profCoord += delta

    return intgI, profI, profRows

def singleFieldZoneIntegral(fieldObj, zones, intgAx, profAx) :
    if len(fieldObj.field) == 1 :
        print(" | Field is uniform, not integrating")
        return [0], [0]
    else :
        intgI, profI, profRows = integrationParameters(zones, intgAx, profAx)
        integratedField = []
        axis = []
        for row in profRows :
            #print("Integrating field ", fieldObj.name, " at time ", fieldObj.time, " over ", intgAx, " at ", profAx, " : ", row[0].C[profI])
            integratedValue = 0
            totRowV = 0
            axis.append(row[0].C[profI])
            for cell in row :
                if debug :
                    print(" | Cell "+str(cell.label)+" has a "+fieldObj.name+" of "+str(fieldObj.field[cell.label])+" at time "+str(fieldObj.time))
                totRowV += cell.V
                integratedValue += fieldObj.field[cell.label]*cell.V
            integratedValue /= totRowV
            integratedField.append(integratedValue)
        return axis, integratedField

def fieldsProductZoneIntegral(fieldObj1, fieldObj2, zones, intgAx, profAx) :
    f1Len = len(fieldObj1.field)
    f2Len = len(fieldObj2.field)
    if f1Len == 1 and f2Len == 1 :
        print(" | Fields are both uniform, not integrating")
        return [0], [0]
    else :
        intgI, profI, profRows = integrationParameters(zones, intgAx, profAx)
        integratedField = []
        axis = []
        for row in profRows :
            integratedValue = 0
            totRowV = 0
            axis.append(row[0].C[profI])
            for cell in row :
                totRowV += cell.V
                i1 = cell.label
                i2 = i1
                if f1Len == 1 :
                    i1 = 0
                elif f2Len == 1 :
                    i2 = 0
                integratedValue += fieldObj1.field[i1]*fieldObj2.field[i2]*cell.V
            integratedValue /= totRowV
            integratedField.append(integratedValue)
        return axis, integratedField

def fieldMinMax(fieldObj, zones, f=None) :
    if f is None:
        f = 1
    global gCellList
    if len(fieldObj.field) == 1 :
        print(" | Uniform field, min = max = "+str(fieldObj.field[0]))
        minC = "nan"
        minC = "nan"
    else :
        cells = [cell for cell in gCellList if cell.zone in zones]
        field = []
        volField = []
        totVol = 0
        avg = 0
        for cell in cells :
            field.append(fieldObj.field[cell.label])
            volField.append(cell.V)
        for i in range(0, len(field)) :
            totVol += volField[i]
            avg += (volField[i]*field[i])
        avg /= totVol
        minLabel = fieldObj.field.index(min(field))
        maxLabel = fieldObj.field.index(max(field))
        minC = [cell for cell in gCellList if cell.label == minLabel][0].C
        maxC = [cell for cell in gCellList if cell.label == maxLabel][0].C
        print(" | Zones: "+" ".join(str(zone) for zone in zones))
        print(" | Max of "+fieldObj.name+" is "+str(max(field)*f)+" "+fieldObj.units+" at"+str(maxC)+"m")
        print(" | Min of "+fieldObj.name+" is "+str(min(field)*f)+" "+fieldObj.units+" at"+str(minC)+"m")
        print(" | Avg of "+fieldObj.name+" is "+str(avg*f)+" "+fieldObj.units+"\n")

### CASE SPECIFIC FUNCTIONS :

def massFlow(zones, surf, time) :
    print("MASS FLOW AT TIME "+str(time)+"s IN ZONES : "+" ".join(str(zone) for zone in zones))
    Uz = loadField("U", "m/s", str(time), "fluidRegion", "vector").cmpt(2)
    rho = loadField("rho", "kg/m3", str(time), "fluidRegion", "scalar")
    axis, flux = fieldsProductZoneIntegral(Uz, rho, zones, "x", "z")

    flow = []
    for fluxi in flux :
        fluxi *= surf
        flow.append(fluxi)
    #printTable(axis, "z", flow, "flow", "kg/s")
    avgFlow = 0
    for flowi in flow :
        avgFlow += flowi
    avgFlow /= len(flow)
    print(" | Total core flow rate = "+str(avgFlow)+" kg/s\n\n")
    Umag = loadField("U", "m/s", str(time), "fluidRegion", "vector").mag()
    fieldMinMax(Umag, zones, (1/0.32593527))

def coolantTemperature(zones, time) :
    print("TEMPERATURE AT TIME "+str(time)+"s IN ZONES : "+" ".join(str(zone) for zone in zones))
    T = loadField("T", "K", str(time), "fluidRegion", "scalar")
    print(" | In core")
    fieldMinMax(T, zones)
    #axis, Tax = singleFieldZoneIntegral(T, zones, "x", "z")
    #printTable(axis, "z", Tax, T.name, T.units)
    print(" | In diagrid")
    fieldMinMax(T, ["diagrid"])
    print(" | In upper plenum")
    fieldMinMax(T, ["sodium"])
    #axis, Tax = singleFieldZoneIntegral(T, ["diagrid"], "x", "z")
    #printTable(axis, "z", Tax, T.name, T.units)

def fuelCladTemperature(zones, time) :
    print("FUEL TEMPERATURE AT TIME "+str(time)+"s IN ZONES : "+" ".join(str(zone) for zone in zones))
    TiF = loadField("TinnerFuel", "K", str(time), "fluidRegion", "scalar")
    print(" | Inner fuel")
    fieldMinMax(TiF, zones)
    TaF = loadField("TavFuel", "K", str(time), "fluidRegion", "scalar")
    print(" | Average fuel")
    fieldMinMax(TaF, zones)
    ToF = loadField("TouterFuel", "K", str(time), "fluidRegion", "scalar")
    print(" | Outer fuel")
    fieldMinMax(ToF, zones)
    print("CLAD TEMPERATURE AT TIME "+str(time)+"s IN ZONES : "+" ".join(str(zone) for zone in zones))
    TiC = loadField("TinnerClad", "K", str(time), "fluidRegion", "scalar")
    print(" | Inner clad")
    fieldMinMax(TiC, zones)
    TaC = loadField("TavClad", "K", str(time), "fluidRegion", "scalar")
    print(" | Average clad")
    fieldMinMax(TaC, zones)
    ToC = loadField("TouterClad", "K", str(time), "fluidRegion", "scalar")
    print(" | Outer clad")
    fieldMinMax(ToC, zones)

def volFuelPower(zones, time) :
    print("VOLFUELPOWER AT TIME "+str(time)+"s IN ZONES : "+" ".join(str(zone) for zone in zones))
    VFP = loadField("volFuelPower", "W/m3", str(time), "fluidRegion", "scalar")
    fieldMinMax(VFP, zones)
    #axis, VFPax = singleFieldZoneIntegral(VFP, zones, "x", "z")
    #printTable(axis, "z", VFPax, VFP.name, VFP.units)
    #axis, VFPrad = singleFieldZoneIntegral(VFP, zones, "z", "x")
    #printTable(axis, "r", VFPrad, VFP.name, VFP.units)

def volFuelPowerNeutro(zones, time) :
    print("VOLFUELPOWER AT TIME "+str(time)+"s IN ZONES : "+" ".join(str(zone) for zone in zones))
    VFP = loadField("volFuelPowerNeutro", "W/m3", str(time), "neutroRegion", "scalar")
    fieldMinMax(VFP, zones)
    #axis, VFPax = singleFieldZoneIntegral(VFP, zones, "x", "z")
    #printTable(axis, "z", VFPax, VFP.name, VFP.units)
    #axis, VFPrad = singleFieldZoneIntegral(VFP, zones, "z", "x")
    #printTable(axis, "r", VFPrad, VFP.name, VFP.units)

def pressureDrop(zones, time) :
    print("PRESSURE DROP AT TIME "+str(time)+"s IN ZONES : "+" ".join(str(zone) for zone in zones))
    p = loadField("p_rgh", "Pa", str(time), "fluidRegion", "scalar")
    axis, p = singleFieldZoneIntegral(p, zones, "x", "z")
    print(" | Pressure drop = "+str(max(p)-min(p))+" Pa\n\n")

def coreRadDisplacement(zones, time) :
    print("CORE RADIAL DISPLACEMENT AT TIME "+str(time)+"s IN ZONES : "+" ".join(str(zone) for zone in zones))
    Dx = loadField("Disp", "%", str(time), "thermoMechanicalRegion", "vector").cmpt(0)
    fieldMinMax(Dx, zones, (100/2.44))
    #axis, Dxrad = singleFieldZoneIntegral(Dx, zones, "z", "x")
    #printTable(axis, "z", Dxrad, Dx.name, Dx.units)

def coreAxialDisplacement(zones, time) :
    print("CORE AXIAL DISPLACEMENT AT TIME "+str(time)+"s IN ZONES : "+" ".join(str(zone) for zone in zones))
    Dz = loadField("fuelDisp", "%", str(time), "thermoMechanicalRegion", "scalar")
    fieldMinMax(Dz, zones, 100)
    #axis, Dzax = singleFieldZoneIntegral(Dz, zones, "x", "z")
    #printTable(axis, "z", Dzax, Dz.name, Dz.units)    

def flux(zones, time) :
    print("NEUTRON FLUX AT TIME "+str(time)+"s IN ZONES : "+" ".join(str(zone) for zone in zones))
    phi = loadField("flux0", "1/m2/s", str(time), "neutroRegion", "scalar")
    fieldMinMax(phi, zones)
    #axis, phiRad = singleFieldZoneIntegral(phi, zones, "z", "x")
    #printTable(axis, "r", phiRad, phi.name, phi.units)
    #print("AXIAL FLUX DISTRIBUTION AT TIME ", str(time), "s IN ZONES : ", zones)
    #axis, phiAx = singleFieldZoneIntegral(phi, zones, "x", "z")
    #printTable(axis, "z", phiAx, phi.name, phi.units)

def rhok(zones, time) :
    print("RHOK AT TIME "+str(time)+"s IN ZONES : "+" ".join(str(zone) for zone in zones))
    rhok = loadField("rhok", "kg/m3", str(time), "fluidRegion", "scalar")
    fieldMinMax(rhok, zones, 860)
    #axis, rhokax = singleFieldZoneIntegral(rhok, zones, "x", "z")
    #printTable(axis, "z", rhokax, rhok.name, rhok.units)
    #axis, rhokrad = singleFieldZoneIntegral(rhok, zones, "z", "x")
    #printTable(axis, "r", rhokrad, rhok.name, rhok.units)

def postProcess(parametersList) :

    global gPointList
    global gFaceList
    global gCellList

    root = os.getcwd()

    endTime = "endTime_"+str(parametersList[10])
    processTH = bool(parametersList[4])
    processEn = bool(parametersList[5])
    processNt = bool(parametersList[6])
    processTM = bool(parametersList[8])

    if processTH :
        polyMeshDir = root+"/constant/fluidRegion/polyMesh"
        os.chdir(polyMeshDir)
        print("========================================")
        print("\nReading fluidRegion polyMesh\n")
        gPointList = []
        gFaceList = []
        gCellList = []
        readPoints()
        readFaces()
        readCells()
        massFlow(["innerCore", "outerCore"], 18.703786, endTime) #18.703786
        #volFuelPower(["innerCore", "outerCore"], endTime)
        pressureDrop(["innerCore", "outerCore"], endTime)
        os.chdir(root)

    if processEn :
        polyMeshDir = root+"/constant/fluidRegion/polyMesh"
        os.chdir(polyMeshDir)
        print("========================================")
        print("\nReading fluidRegion polyMesh\n")
        gPointList = []
        gFaceList = []
        gCellList = []
        readPoints()
        readFaces()
        readCells()
        coolantTemperature(["innerCore", "outerCore"], endTime)
        fuelCladTemperature(["innerCore", "outerCore"], endTime)
        rhok(["innerCore", "outerCore"], endTime)
        os.chdir(root)

    if processNt :
        polyMeshDir = root+"/constant/neutroRegion/polyMesh"
        os.chdir(polyMeshDir)
        print("========================================")
        print("Reading neutroRegion polyMesh\n")
        gPointList = []
        gFaceList = []
        gCellList = []
        readPoints()
        readFaces()
        readCells()
        volFuelPowerNeutro(["innerCore", "outerCore"], endTime)
        flux(["innerCore", "outerCore"], endTime)
        os.chdir(root)

    if processTM :
        polyMeshDir = root+"/constant/thermoMechanicalRegion/polyMesh"
        os.chdir(polyMeshDir)
        print("========================================")
        print("Reading thermoMechanicalRegion polyMesh\n")
        gPointList = []
        gFaceList = []
        gCellList = []
        readPoints()
        readFaces()
        readCells()
        coreRadDisplacement(["outerCore"], endTime)
        coreAxialDisplacement(["innerCore", "outerCore"], endTime)
        os.chdir(root)