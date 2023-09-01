### IMPORTS ###

import subprocess
import shutil
import os

### GLOBAL

root = os.getcwd()

### FUNCTIONS ###

def preprocessFolders(parameterList, i) :

    # Variables to make the code more readable
    startFromCase = parameterList[9]
    index = parameterList[10]

    # Backup the 0 folder at the very first case
    if i == 0 :
        try :
            shutil.copytree("0", "0.orig")
        except :
            pass

    # Backup 0 folder
    try : 
        print("Backing up time directory 0 for case "+str(index)+" ...")
        shutil.copytree("0", "startTime_"+str(index))
    except :
        pass

    # Set the 0 folder
    startFromCase = parameterList[9]
    index = parameterList[10]
    if startFromCase != -1 :
        shutil.rmtree("0")
        shutil.copytree("endTime_"+str(startFromCase), "0")

    # Remove all other eventual time directories from previous cases
    timeDirs = getTimeDirs(root)
    for timeDir in timeDirs :
        if timeDir != "0" :
            print("Removing time directory "+timeDir+" ...")
            subprocess.call(["rm", "-r", timeDir])

    # Remove eventual processor directories from previous cases
    procDirs = [procDir for procDir in os.listdir(root) if "processor" in procDir]
    for procDir in procDirs :
        print("Removing processor directory "+procDir)
        subprocess.call(["rm", "-r", procDir])

def postprocessFolders(parameterList) :

    # Extra variables just to make the code more readables
    index = parameterList[10]
    endTime = parameterList[0]

    if int(endTime) - float(endTime) == 0 :
        endTime = str(int(endTime))
    else :
        endTime = str(float(endTime))

    print("Backing up end time directory "+endTime+ " of case "+str(index))
    os.chdir(endTime)
    if "uniform" in os.listdir(os.getcwd()) :
        subprocess.call(["rm", "-r", "uniform"])
    os.chdir(root)
    try :
        shutil.copytree(endTime, "endTime_"+str(index))
    except :
        pass

def getTimeDirs(rootDir) :
    timeDirs = []
    subDirs = [name for name in os.listdir(rootDir) if os.path.isdir(os.path.join(rootDir, name))]
    for i in range(0, len(subDirs)):
        try:
            float(subDirs[i])
            timeDirs.append(subDirs[i])
        except:
            pass
    return timeDirs