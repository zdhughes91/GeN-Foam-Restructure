### IMPORTS ###

from dependencies.fileManipulation import *
from dependencies.folderManipulation import *
from dependencies.postProcessing import *
import sys

### FUNCTIONS ###

def setupCase(regions, parameterList) :

    editKeywordInFile("system/controlDict", "endTime", parameterList[0])
    editKeywordInFile("system/controlDict", "deltaT", parameterList[1])
    editKeywordInFile("system/controlDict", "writeInterval", parameterList[2])
    editKeywordInFile("system/fvSolution",  "tightlyCoupled", parameterList[3])
    editKeywordInFile("system/controlDict", "solveFluidMechanics", parameterList[4])
    editKeywordInFile("system/controlDict", "solveEnergy", parameterList[5])
    editKeywordInFile("system/controlDict", "solveNeutronics", parameterList[6])
    editKeywordInFile("constant/"+regions[1]+"/neutronicsProperties", "eigenvalueNeutronics", parameterList[7])
    editKeywordInFile("system/controlDict", "solveThermalMechanics", parameterList[8])
    editKeywordInFile("constant/"+regions[0]+"/fuelProperties", "fuelVolPower", parameterList[11])

def processLog(process, printLog, saveLog, saveName = None) :
    
    logLines = []
    for line in iter(process.stdout.readline, b'') :
        logLine = line.decode('ascii')
        if printLog :
            sys.stdout.write(logLine)
        if saveLog :
            logLines.append(logLine)
    if saveLog :
        save(saveName, logLines)

def runCase(solverName, regions, parallelOptions, parameterList, i) :

    preprocessFolders(parameterList, i)
    setupCase(regions, parameterList)

    # Extra variables just to make it more readable
    
    index = str(parameterList[10])
    parallel = bool(parallelOptions[0])
    scheme = str(parallelOptions[1])
    nX = int(parallelOptions[2])
    nY = int(parallelOptions[3])
    nZ = int(parallelOptions[4])
    

    if parallel :

        # Decompose domain
        editDecomposeParDict(regions, scheme, [nX, nY, nZ])
        for region in regions :
            print("Running decomposePar for case "+index+" | region : "+region+" | method : "+scheme+" | coeffs (X Y Z): "+str(nX)+" "+str(nY)+" "+str(nZ))
            decomposePar = subprocess.Popen(["decomposePar", "-region", region], stdout=subprocess.PIPE)
            processLog(decomposePar, False, True, "decomposeParLog_"+region+"_"+index)
            decomposePar.wait()
        
        # Run in parallel
        nProc = nX*nY*nZ
        print("Running "+solverName+" in parallel for case "+index+" ...")
        solver = subprocess.Popen(["mpirun", "-np", str(nProc), solverName, "-parallel"], stdout=subprocess.PIPE)    
        processLog(solver, False, True, solverName+"Log_"+index)
        solver.wait()
        
        # Reconstruct domain
        for region in regions :
            print("Running reconstructPar for case "+index+" | region : "+region)
            reconstructPar = subprocess.Popen(["reconstructPar", "-region", region], stdout=subprocess.PIPE)
            processLog(reconstructPar, True, False)
            reconstructPar.wait()
    
    else :

        # Run single core
        print("Running "+solverName+" for case "+index+" ...")
        solver = subprocess.Popen([solverName], stdout=subprocess.PIPE)    
        processLog(solver, False, True, solverName+"Log_"+index)
        solver.wait()

    postprocessFolders(parameterList)
    
    with open("postProcess_"+index, "w") as out :
        sys.stdout = out
        postProcess(parameterList)
    sys.stdout = sys.__stdout__

    i += 1
    return i