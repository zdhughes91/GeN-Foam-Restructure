### IMPORTS ###

import shutil
import math

### FUNCTIONS ###

### Generic read/save

def read(filename) :
    lines = [];
    with open(filename) as file :
        lines = file.readlines()
    file.close()
    return lines;

def save(filename, lines) :
    with open(filename, "w") as file :
        for line in lines :
            file.writelines(line)
    file.close()

### OpenFOAM dict manipulation

def editKeywordInFile(file, key, value) :
    lines = read(file)
    for linei, line in enumerate(lines) :
        tmp = line.split()
        if len(tmp) > 0 :
            if tmp[0] == key :
                lines[linei] = key+" "+str(value)+";\n"
    save(file, lines)

def editDecomposeParDict(regions, method, coeffs) :
    for region in regions :
        path = "system/"+region+"/decomposeParDict"
        lines = read(path)
        coeffsSubDict = False
        nSubDomains = coeffs[0]*coeffs[1]*coeffs[2]
        for linei, line in enumerate(lines) :
            tmp = line.split()
            if len(tmp) > 0 :
                if tmp[0] == "numberOfSubdomains" :
                    lines[linei] = "numberOfSubdomains   "+str(nSubDomains)+";\n"
                elif tmp[0] == "method" :
                    lines[linei] = "method   "+str(method)+";\n"
                elif tmp[0] == method+"Coeffs" :
                    coeffsSubDict = True
                elif coeffsSubDict == True :
                    if tmp[0] == "n" :
                        lines[linei] = "n      ("+str(coeffs[0])+" "+str(coeffs[1])+" "+str(coeffs[2])+");\n"
                    elif tmp[0] == "}" :
                        coeffsSubDict = False
        save(path, lines)
    shutil.copy(path, "system/decomposeParDict")