"""
SerpentToFoamXS

Author: Thomas Guilbaud, EPFL/Transmutex SA
Last Update: 24/03/2022
"""

from . import widget
from . import utils


class Parameter(object):
    def __init__(self,
        name: str,
        value,
        default="",
        unit: str="",
        comment: str="",
        isWrite: bool=True,
        isCommented: bool=False
    ):
        self.name     = name
        self.value    = value
        self.default  = default
        self.unit     = unit
        self.comment  = comment
        self.isVector = name[len(name)-11:] == "Orientation"
        self.isWrite  = isWrite
        self.isCommented = isCommented

        # Widget
        self.labelText = None
        self.entry     = None
        self.labelUnit = None
        self.row       = 0

    #--------------------------------------------------------------------------*
    # Setter
    def setWidget(self, frame, row: int) -> None:
        self.labelText, self.entry, self.labelUnit = widget.line(
            frame,
            self.name,
            self.value,
            default=self.default,
            row=row,
            unit=self.unit
        )
        self.row = row

    def showWidget(self) -> None:
        self.labelText.grid(row=self.row, column=0, sticky="W")
        self.entry.grid(row=self.row, column=2, sticky="W")
        if (self.labelUnit != None):
            self.labelUnit.grid(row=self.row, column=0, sticky="E")

    def hideWidget(self) -> None:
        self.labelText.grid_forget()
        self.entry.grid_forget()
        if (self.labelUnit != None):
            self.labelUnit.grid_forget()

    #--------------------------------------------------------------------------*
    # Writer
    def write(self, outputFile, isWriteForNuclearData: bool=False) -> None:
        # If parameter is meant to not be written
        if (not self.isWrite): # and isWriteForNuclearData):
            return

        if (self.isVector):
            if (isWriteForNuclearData):
                outputFile.write("{:25} ( {} );\n\n".format(self.name, " ".join(self.value.get().split())))
            else:
                outputFile.write("{:25} {}\n\n".format(self.name, " ".join(self.value.get().split())))

        # General scalar/boolean case
        elif (utils.isStringVar(self.value)):
            # special as expansionFromNominal can be used for radial and axial
            # expansion. For userInputSFXS, differentiation between expansionFromNominalR
            # and expansionFromNominalA. For nuclearData, same key word but
            # in different file. Simply remove the last letter (R or A) in case
            # "expansionFromNominal"
            utils.writeScalarLine(outputFile,
                # remove last letter (R or A), same key word in output nuclearData files
                self.name[0:-1] if ("expansionFromNominal" in self.name) else self.name,
                self.value.get(),
                self.unit,
                comment=self.comment,
                isWriteForNuclearData=isWriteForNuclearData,
                isCommented=self.isCommented
            )
        elif (utils.isIntVar(self.value)):
            if (isWriteForNuclearData):
                utils.writeScalarLine(outputFile,
                    self.name,
                    utils.binaryToStrBool(self.value.get()),
                    comment=self.comment,
                    isWriteForNuclearData=isWriteForNuclearData,
                    isCommented=self.isCommented
                )
            else:
                utils.writeScalarLine(outputFile,
                    self.name,
                    self.value.get(),
                    comment=self.comment,
                    isWriteForNuclearData=isWriteForNuclearData,
                    isCommented=self.isCommented
                )
