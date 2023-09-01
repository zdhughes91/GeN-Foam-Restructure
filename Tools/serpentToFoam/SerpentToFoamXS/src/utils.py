"""
SerpentToFoamXS

Generic usefull functions

Author: Thomas Guilbaud, EPFL/Transmutex SA
Last Update: 24/03/2022
"""

import tkinter
import os.path


#------------------------------------------------------------------------------*
# Boolean

def isInputFileNameDefined(filename: str, isCreate: bool=False) -> bool:
    if (os.path.isfile(filename) or isCreate):
        return(filename not in ["", "Open"])
    return(False)

def isNameVector(name: str="") -> bool:
    if (len(name) >= 11):
        return(name[len(name)-11:] == "Orientation")
    else:
        return(False)

def isBinary(var: int=0) -> bool:
    return(var in [0, 1])

def isStringVar(var) -> bool:
    return(type(var) == type(tkinter.StringVar()))

def isIntVar(var) -> bool:
    return(type(var) == type(tkinter.IntVar()))

#------------------------------------------------------------------------------*
# Printer

def printFileName(filename: str, prefixMessage: str="") -> None:
    if (isInputFileNameDefined(filename)):
        if (prefixMessage == ""):
            print("File selected: {}\n".format(filename))
        else:
            print("{} {}\n".format(prefixMessage, filename))

#------------------------------------------------------------------------------*
# Translate

def binaryToStrBool(bin: int) -> str:
    return("true" if bin == 1 else "false")

def binaryToBool(bin: int) -> bool:
    return(bin == 1)

#------------------------------------------------------------------------------*
# Writer

def writeScalarLine(
    file,
    parameter: str,
    value,
    unit: str="",
    comment: str="",
    isBoolValue: bool=False,
    isWriteForNuclearData: bool=False,
    isCommented: bool=False
) -> None:
    if (isCommented):
        file.write("// ")

    if (isBoolValue):
        file.write("{:25} {}{}  // 1: true; 0: false\n\n".format(
            parameter, value, ";" if isWriteForNuclearData else ""
        ))
    else:
        file.write("{:25} {}".format(parameter, value))
        if (isWriteForNuclearData):
            file.write(";")
        if (unit != ""):
            file.write(" // [{}]".format(unit))
        if (comment != ""):
            file.write("{} {}".format("," if (unit != "") else " //", comment))
        file.write("\n\n")

def writeSimpleList(outputFile, label: str, liste: list) -> None:
    outputFile.write("\t\t{:18} nonuniform List<scalar> {} ( ".format(label, len(liste)))
    for value in liste:
        outputFile.write("{:.6e} ".format(value))
    outputFile.write(");\n")

def writeMatrix(outputFile, label: str, matrix: list) -> None:
    n = int(len(matrix) ** 0.5)
    outputFile.write("\t\t{0:18} {1} {1} ( ".format(label, n))
    idx = 0
    for i in range(n):
        outputFile.write("( ")
        for j in range(n):
            outputFile.write("{:.6e} ".format(matrix[idx]))
            idx += 1
        outputFile.write(") ")
    outputFile.write(");\n")
