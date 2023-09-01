"""
SerpentToFoamXS

Author: Thomas Guilbaud, EPFL/Transmutex SA
Last Update: 18/03/2022
"""

#------------------------------------------------------------------------------*
# Imports
import sys

#------------------------------------------------------------------------------*
# Key words
def allWords() -> list:
    return(["--all", "-a"])

def bashWords() -> list:
    return(["--bash", "-b"])

def helpWords() -> list:
    return(["--help", "-h"])

#------------------------------------------------------------------------------*
# Boolean
def isAllKey(key: str) -> bool:
    return(key in allWords())

def isBashKey(key: str) -> bool:
    return(key in bashWords())

def isHelpKey(key: str) -> bool:
    return(key in helpWords())

def isExtractAllMode() -> bool:
    for key in sys.argv:
        if (isAllKey(key)):
            return(True)
    return(False)

def isBashMode() -> bool:
    for key in sys.argv:
        if (isBashKey(key)):
            return(True)
    return(False)

#------------------------------------------------------------------------------*
# Printer
def printUsage() -> None:
    print("\nSerpentToFoamXS Application for cross-section extraction\n")
    print("Options:")
    print("\t-b, --bash\tBash mode")
    print("\t-a, --all \tExtract all the nuclearData files\n")
    print("Usage in GUI : ./SerpentToFoamXS")
    print("  or         : ./SerpentToFoamXS userInputFile")
    print("  or         : ./SerpentToFoamXS userInputFile serpentFile_res.m")
    print("Usage in Bash: ./SerpentToFoamXS --bash userInputFile")
    print("  or         : ./SerpentToFoamXS --bash userInputFile serpentFile_res.m")
    print("Usage in Bash: ./SerpentToFoamXS --bash --all userInputFile")
    print("  or         : ./SerpentToFoamXS --bash --all userInputFile serpentFile_res.m\n")

#------------------------------------------------------------------------------*
# Main function
def entries():
    # User Input File Path
    userInputFileName, serpentFileName = "", ""

    isBash = isBashMode()
    isExtractAll = isExtractAllMode()
    isUserCmdCorrect = True

    # Graphical User Interface Mode without input
    if (len(sys.argv) == 1):
        pass
    # Graphical User Interface Mode with input
    elif (len(sys.argv) == 2 and not isBash):
        if (isHelpKey(sys.argv[1])):
            printUsage()
            isUserCmdCorrect = False
        else:
            userInputFileName = sys.argv[1]
    elif (len(sys.argv) == 3 and not isBash):
        userInputFileName, serpentFileName = sys.argv[1], sys.argv[2]
    # Bash Mode
    elif (len(sys.argv) >= 3 and isBash):
        for key in sys.argv[1:]:
            if (not isBashKey(key) and not isAllKey(key) and not isHelpKey(key)):
                if (userInputFileName == ""):
                    userInputFileName = key
                elif (serpentFileName == ""):
                    serpentFileName = key
                else:
                    print("Too many parameters")
    # Wrong usage
    else:
        printUsage()
        isUserCmdCorrect = False

    return(userInputFileName, serpentFileName, isUserCmdCorrect, isBash, isExtractAll)
