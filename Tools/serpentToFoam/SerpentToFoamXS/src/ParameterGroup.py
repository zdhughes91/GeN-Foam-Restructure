"""
SerpentToFoamXS

Author: Thomas Guilbaud, EPFL/Transmutex SA
Last Update: 18/03/2022
"""

from . import Parameter


class ParameterGroup(object):
    def __init__(self,
        name: str,
        accronym: str,
        fileName: str=""
    ):
        self.name     = name
        self.accronym = accronym
        self.fileName = fileName
        self.widgetList = []

        self.parameterList = []

    #--------------------------------------------------------------------------*
    # Adder
    def addParameter(self,
        name,
        value,
        default="",
        unit: str="",
        comment: str="",
        isWrite: bool=True,
        isCommented: bool=False
    ) -> None:
        self.parameterList.append(
            Parameter.Parameter(name, value, default=default, unit=unit, isWrite=isWrite
        ))

    #--------------------------------------------------------------------------*
    # Setter
    def setParameterFromLine(self, line: list):
        for parameter in self.parameterList:
            # General scalar/boolean case
            if (line[0] == parameter.name):
                # Special case for vector, end with "Orientation"
                if (parameter.isVector):
                    if (len(line) >= 4): # name + 3 components vector
                        parameter.value.set(" ".join(line[1:4]))
                    else:
                        parameter.value.set(line[1]+" 0 0")
                # General scalar/boolean cases
                else:
                    parameter.value.set(line[1])

    def setWidgets(self, frame, row: int) -> int:
        for parameter in self.parameterList:
            parameter.setWidget(frame, row)
            row += 1
        return(row)

    #--------------------------------------------------------------------------*
    # Displayer
    def displayParameters(self, isDisplay: bool):
        if (isDisplay):
            self.showParameters()
        else:
            self.hideParameters()

    def showParameters(self) -> None:
        for parameter in self.parameterList:
            parameter.showWidget()

    def hideParameters(self) -> None:
        for parameter in self.parameterList:
            parameter.hideWidget()

    #--------------------------------------------------------------------------*
    # Writer
    def writeParameters(self, outputFile, isWriteForNuclearData: bool=False) -> None:
        for parameter in self.parameterList:
            parameter.write(outputFile, isWriteForNuclearData)
