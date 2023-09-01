"""
SerpentToFoamXS

Application for Easy Extraction of Serpent output files for GeN-Foam.
Can be use in Bash mode or using a Graphical User Interface.

Author: Thomas Guilbaud, EPFL/Transmutex SA
Last Update: 18/03/2022
"""

#-------------------------------------------------------------------------------*
# Imports
import tkinter as tk
import sys

import src

#-------------------------------------------------------------------------------*
# Main
if __name__ == "__main__":
    # User Input File Path and System
    userInputFileName, serpentFileName, isUserCmdCorrect, isBashMode, isExtractAll = src.system.entries()

    #---------------------------------------------------------------------------*
    # Build Application
    if (isUserCmdCorrect):
        root = tk.Tk()

        # Call the Application
        app = src.App(
            root, userInputFileName, serpentFileName, isBashMode, isExtractAll
        )

        # Main loop
        if (not isBashMode):
            root.mainloop()

# End of Main
#-------------------------------------------------------------------------------*
