### IMPORTS ###

import re

### FUNCTIONS ###

def readDict() :

    lines = []
    file = open("runDict", "r")
    tmp = file.readlines()
    for line in tmp :
        line = re.split("\s", line)
        while '' in line:
            line.remove('')
        if len(line) > 0 :
            lines.append(line)
    read = True
    regions = []
    parallelOptions = []
    parameterTable = []
    for line in lines :
        if line[0] == "'''" :
            if read == True :
                read = False
            else :
                read = True
        elif read and line[0] != "#" :
            if line[0] == "regions" :
                for region in line :
                    if region != "regions" :
                        regions.append(region)
            elif line[0] == "solverName" :
                solverName = line[1]
            elif line[0] == "parallel" :
                for option in line :
                    if option != "parallel" :
                        parallelOptions.append(option)
            else :
                parameterList = []
                for i, parameter in enumerate(line) :
                    if i <= 2 :
                        parameterList.append(float(parameter))
                    else :
                        try :
                            parameterList.append(int(parameter))
                        except :
                            parameterList.append(parameter)
                parameterTable.append(parameterList)
    return solverName, regions, parallelOptions, parameterTable