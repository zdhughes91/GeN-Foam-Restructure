"""
SerpentToFoamXS

Application for Easy Extraction of Serpent output files for GeN-Foam

Author: Thomas Guilbaud, EPFL/Transmutex SA
Last Update: 25/03/2022
"""

#------------------------------------------------------------------------------*
# Imports
import tkinter as tk
import tkinter.ttk as ttk
import datetime

from . import ParameterGroup
from . import utils
from . import extraction
from . import widget


#------------------------------------------------------------------------------*
# App
class App(object):
    def __init__(
        self,
        master,
        userInputFileName: str="",
        serpentFileName: str="",
        isBashMode: bool=False,
        isExtractAll: bool=False
    ) -> None:
        #----------------------------------------------------------------------*
        # Initial prints
        self.version      = "1.0 - 25/03/2022"
        self.isBashMode   = isBashMode
        self.isExtractAll = isExtractAll

        print("\n"+
            "SerpentToFoamXS - Version: {}\n\n".format(self.version)+
            "This application converts a SERPENT(2) output into the crossSections dictionary\n"+
            "for neutronDiffusionFoam.\n\n"
        )

        # Mode used
        if (self.isBashMode):
            print("Bash mode\n")
        else:
            print("Graphical User Inteface (GUI)\n")

        # Check user input file name
        if (utils.isInputFileNameDefined(userInputFileName)):
            print("User input file '{}'\n".format(userInputFileName))
        else:
            print("No user input file provided in cmd\n")

        # Check Serpent output file name
        if (utils.isInputFileNameDefined(serpentFileName)):
            print("Serpent output file {}\n".format(serpentFileName))
        else:
            print("No Serpent output file provided in cmd\n")

        #----------------------------------------------------------------------*
        # Set Variables

        # Master for GUI
        self.master = master
        self.master.title("Serpent to Foam XS")

        # Set style
        style = ttk.Style(self.master)
        style.theme_use('clam')

        # User variables
        self.userInputFileName = tk.StringVar(self.master)
        self.userInputFileName.set(userInputFileName)
        self.serpentFileName   = tk.StringVar(self.master)
        self.serpentFileName.set(serpentFileName)
        self.coreStateSelected = tk.StringVar(self.master)

        # Diffusion Neutronics
        self.isDiffusionNeutronicsExtract = tk.IntVar(self.master) # 1: true; 0: false
        self.diffusionNeutronicsKeyWord   = "Diffusion-Neutronics"
        self.diffusionNeutronicsParameterList = ParameterGroup.ParameterGroup(self.diffusionNeutronicsKeyWord, "DiffNeutron")
        self.diffusionNeutronicsParameterList.addParameter("adjustDiscFactors", tk.IntVar(self.master), default="0")
        self.diffusionNeutronicsParameterList.addParameter("useGivenDiscFactors", tk.IntVar(self.master), default="0")

        # Point-Kinetics variables, only written if self.isPointKineticsExtract == 1: true
        self.isPointKineticsExtract = tk.IntVar(self.master) # 1: true; 0: false
        self.pointKineticsKeyWord   = "Point-Kinetics"
        self.zeroeffSelected        = tk.StringVar(self.master) # physical or effective

        self.pointKineticsParameterList = ParameterGroup.ParameterGroup(self.pointKineticsKeyWord, "PK")
        # self.pointKineticsParameterList.addParameter("promptGenerationTime", tk.StringVar(self.master), default="1e-7", unit="s")
        self.pointKineticsParameterList.addParameter("feedbackCoeffFastDoppler", tk.StringVar(self.master), default="0", unit="1/K")
        # feedbackCoeffTFuel = tk.StringVar(self.master) # in case the code had to display multiple variables
        self.pointKineticsParameterList.addParameter("feedbackCoeffTFuel", tk.StringVar(self.master), default="0", unit="1/K")
        self.pointKineticsParameterList.addParameter("feedbackCoeffTClad", tk.StringVar(self.master), default="0", unit="1/K")
        self.pointKineticsParameterList.addParameter("feedbackCoeffTCool", tk.StringVar(self.master), default="0", unit="1/K")
        self.pointKineticsParameterList.addParameter("feedbackCoeffRhoCool", tk.StringVar(self.master), default="0", unit="1/(kg/m3)")
        self.pointKineticsParameterList.addParameter("feedbackCoeffTStruct", tk.StringVar(self.master), default="0", unit="1/K")

        # Perturbation variables, Nominal and non-Nominal cases
        self.reactorState = ParameterGroup.ParameterGroup("Reactor State", "S", "reactorState")
        self.reactorState.addParameter("pTarget", tk.StringVar(self.master), default="1", unit="W")
        # self.reactorState.addParameter("keff", tk.StringVar(self.master), default="1")
        self.reactorState.addParameter("externalReactivity", tk.StringVar(self.master), default="0")

        self.nominalKeyWord = "Nominal"
        nominalState = ParameterGroup.ParameterGroup(self.nominalKeyWord, "N", "nuclearData")
        nominalState.addParameter("fastNeutrons", tk.IntVar(self.master), default=1) # 1: true; 0: false
        nominalState.addParameter(self.pointKineticsKeyWord, self.isPointKineticsExtract, default=1, isWrite=False)
        nominalState.addParameter(self.diffusionNeutronicsKeyWord, self.isDiffusionNeutronicsExtract, default=0, isWrite=False)

        axialExpansionState = ParameterGroup.ParameterGroup("Axially expansion", "A", "nuclearDataAxialExp")
        axialExpansionState.addParameter("expansionFromNominalA", tk.StringVar(self.master), default="5.0e-2") # Last letter is removed during nuclearData writing

        radialExpansionState = ParameterGroup.ParameterGroup("Radially expansion", "R", "nuclearDataRadialExp")
        radialExpansionState.addParameter("expansionFromNominalR", tk.StringVar(self.master), default="5.0e-2") # Last letter is removed during nuclearData writing
        radialExpansionState.addParameter("radialOrientation", tk.StringVar(self.master), default="1 0 0")
        radialExpansionState.addParameter("axialOrientation", tk.StringVar(self.master), default="0 0 1")

        fuelTemperatureState = ParameterGroup.ParameterGroup("Fuel temperature", "T", "nuclearDataFuelTemp")
        # fuelTemperatureState.addParameter("feedbackCoeffTFuel", feedbackCoeffTFuel, default="0", unit="1/K", isWrite=False)
        fuelTemperatureState.addParameter("TFuelRef", tk.StringVar(self.master), default="600", unit="K")
        fuelTemperatureState.addParameter("TFuelPerturbed", tk.StringVar(self.master), default="650", unit="K")

        coolantDensityState = ParameterGroup.ParameterGroup("Coolant density", "C", "nuclearDataRhoCool")
        coolantDensityState.addParameter("rhoCoolRef", tk.StringVar(self.master), default="1000", unit="kg/m3")
        coolantDensityState.addParameter("rhoCoolPerturbed", tk.StringVar(self.master), default="1100", unit="kg/m3")

        coolantTemperatureState = ParameterGroup.ParameterGroup("Coolant temperature", "CT", "nuclearDataTCool")
        coolantTemperatureState.addParameter("TCoolRef", tk.StringVar(self.master), default="500", unit="K")
        coolantTemperatureState.addParameter("TCoolPerturbed", tk.StringVar(self.master), default="600", unit="K")

        claddingTemperatureState = ParameterGroup.ParameterGroup("Cladding temperature", "CL", "nuclearDataCladExp")
        claddingTemperatureState.addParameter("TCladRef", tk.StringVar(self.master), default="500", unit="K")
        claddingTemperatureState.addParameter("TCladPerturbed", tk.StringVar(self.master), default="600", unit="K")

        self.coreStateList = [ # do the same for point kinetics and diffusion variables
            self.reactorState,
            nominalState,
            axialExpansionState,
            radialExpansionState,
            fuelTemperatureState,
            coolantDensityState,
            coolantTemperatureState,
            claddingTemperatureState
        ]

        # Serpent name to Foam name (Tkinter.Text -> Text Area)
        self.serpentUnivToFoamText = None

        # App variables
        self.isExtracted     = tk.IntVar(self.master)
        self.outputTextLabel = tk.StringVar(self.master)


        #----------------------------------------------------------------------*
        # Build frame with default parameters
        self.buildDefaultFrame()

        #----------------------------------------------------------------------*
        # User input file and update default parameters
        if (utils.isInputFileNameDefined(self.userInputFileName.get())):
            self.setDataFromUserInputFile(utils.isInputFileNameDefined(self.serpentFileName.get()))
        else:
            # Update the coreState dropdown
            self.coreStateSelected.set(self.coreStateSelected.get())

        #----------------------------------------------------------------------*
        # Bash mode
        if (self.isBashMode):
            self.writeResults()


    #--------------------------------------------------------------------------*
    # Builder

    def buildDefaultFrame(self) -> None:
        # Add a grid for components
        masterFrame = ttk.Frame(self.master, padding=5, relief=tk.GROOVE)
        masterFrame.grid()
        row = 0

        # Vertical seperator between the parameter name and parameter value
        widget.seperator(masterFrame, row+1, column=1, orientation="v")

        #----------------------------------------------------------------------*
        # Header
        # ttk.Label(masterFrame, text="Serpent to Foam XS", font=35, width=20).grid(row=row, column=0, sticky="W")
        logo, image = widget.logo(masterFrame)
        logo.grid(row=row, column=0, sticky="W")
        ttk.Button(masterFrame, text="Quit", command=self.master.destroy, width=10).grid(row=row, column=2, sticky="E")
        row = widget.seperator(masterFrame, row+1)

        # Add the logo to the application
        self.master.iconphoto(False, image)

        #----------------------------------------------------------------------*
        # Inputs

        # User input file
        row = widget.fileSelection(masterFrame, self.userInputFileName, "User input file", row=row)
        self.userInputFileName.trace('w', lambda _, __, ___: self.setDataFromUserInputFile())
        row = widget.seperator(masterFrame, row)

        # Serpent File Selection
        row = widget.fileSelection(masterFrame, self.serpentFileName, "Serpent file", row=row)
        row = widget.seperator(masterFrame, row)

        #----------------------------------------------------------------------*
        # General Paremeters
        # Delayed neutron fraction (zeroeff)
        zeroeffList = ["physical", "effective"] # zero/eff
        row = widget.dropdownSelection(masterFrame, self.zeroeffSelected, "Delayed neutron fraction", zeroeffList, row=row)
        self.zeroeffSelected.trace('w', lambda _, __, ___: self.getDelayedNeutronFractionMessage())

        #----------------------------------------------------------------------*
        # Core State Selection
        row = widget.dropdownSelection(masterFrame, self.coreStateSelected, "Core state variation", [cs.name for cs in self.coreStateList], row=row)
        # Update widget when change core state
        self.coreStateSelected.trace('w', self.widgetActionUpdateDropdownCoreState)
        row = widget.seperator(masterFrame, row)

        # Reference and Perturbed values for core state variables
        for state in self.coreStateList:
            row = state.setWidgets(masterFrame, row)


        #----------------------------------------------------------------------*
        # Point-Kinetics variables
        row = self.pointKineticsParameterList.setWidgets(masterFrame, row)
        self.isPointKineticsExtract.set(0)
        self.isPointKineticsExtract.trace('w', self.widgetActionUpdateDropdownPointKinetics)

        #----------------------------------------------------------------------*
        # Diffusion Neutronics variables
        row = self.diffusionNeutronicsParameterList.setWidgets(masterFrame, row)
        self.isDiffusionNeutronicsExtract.set(0)
        self.isDiffusionNeutronicsExtract.trace('w', self.widgetActionUpdateDropdownDiffusionNeutronics)

        # Horizontal Separation
        row = widget.seperator(masterFrame, row+1)

        #----------------------------------------------------------------------*
        # Universe conversion Serpent/Foam

        # Create textarea widget for universe conversion Serpent/Foam
        row += 30
        ttk.Label(
            masterFrame,
            text="Serpent/Foam cell conversion\n<serpent> <foam> <fuelFrac>"
        ).grid(row=row, column=0, sticky="nw")
        self.serpentUnivToFoamText = tk.Text(masterFrame, height=6, width=26)
        self.serpentUnivToFoamText.grid(row=row, column=2, sticky="w")
        row = widget.seperator(masterFrame, row+1)

        #----------------------------------------------------------------------*
        # Footer: Status message + Save and Extract buttons

        # Status message
        self.outputTextLabel.set("")
        ttk.Label(masterFrame, textvariable=self.outputTextLabel).grid(row=row, column=0, sticky="W")

        # Save and Extract buttons frame
        buttonOptionFrame = ttk.Frame(masterFrame)
        buttonOptionFrame.grid(row=row, column=2)

        # Save the input file for the user
        ttk.Button(buttonOptionFrame, text="Save", command=self.writeUserInputFile, width=11).grid(row=0, column=0)

        # Create output file for GeN-Foam
        self.isExtracted.set(0) # Default at false=0
        ttk.Button(buttonOptionFrame, text="Extract", command=self.writeResults, width=12).grid(row=0, column=1)


    #--------------------------------------------------------------------------*
    # Getter

    def getDelayedNeutronFractionMessage(self) -> str: # zero/eff
        if (self.zeroeffSelected.get() == "physical"):
            return("Physical delayed neutron fraction and spectrum\n")
        elif (self.zeroeffSelected.get() == "effective"):
            return("Effective delayed neutron fraction, prompt neutron spectrum\n")
        return("")


    #--------------------------------------------------------------------------*
    # Printer

    def printOutputTextLabel(self, text: str="") -> None:
        self.outputTextLabel.set(text)
        print(text+"\n")


    #--------------------------------------------------------------------------*
    # Setter

    def setDataFromUserInputFile(self, isSerpentFileNameSet: bool=False) -> None:
        userInputFileName = self.userInputFileName.get()
        if (utils.isInputFileNameDefined(userInputFileName)):
            userInputFile = open(userInputFileName, 'r')

            # Clear the text area to avoid Serpent/Foam corresponce accumulation
            self.serpentUnivToFoamText.delete(1.0, 'end')

            # Read the user input file
            for line in userInputFile.readlines():
                line = line.split()
                if (len(line) > 1):
                    # Variables from core state list
                    for state in self.coreStateList:
                        state.setParameterFromLine(line)

                    # Variables for Point-Kinetics
                    self.pointKineticsParameterList.setParameterFromLine(line)

                    # Variables for Diffusion Neutronics
                    self.diffusionNeutronicsParameterList.setParameterFromLine(line)

                    # Special variables
                    if (line[0] == "serpentFilePath" and not isSerpentFileNameSet):
                        self.serpentFileName.set(line[1])
                    elif (line[0] == "coreState"):
                        self.coreStateSelected.set(" ".join(line[1:3]))
                    elif (line[0] == "zeroeff"):
                        self.zeroeffSelected.set(line[1])
                    elif (line[0] == "isPointKineticsExtract"):
                        if (utils.isBinary(int(line[1]))):
                            self.isPointKineticsExtract.set(int(line[1]))
                    elif (line[0] == "isDiffusionNeutronicsExtract"):
                        if (utils.isBinary(int(line[1]))):
                            self.isDiffusionNeutronicsExtract.set(int(line[1]))
                    elif (line[0] == "serpentUnivToFoam"):
                        self.serpentUnivToFoamText.insert('end', " ".join(line[1:4])+"\n")
        else:
            self.printOutputTextLabel("User input file name not defined")


    #--------------------------------------------------------------------------*
    # Widget Action

    # On change dropdown value
    def widgetActionUpdateDropdownCoreState(self, *args) -> None:
        # Update Core State
        for state in self.coreStateList:
            state.displayParameters(state.name == self.coreStateSelected.get())

        # Update Point-Kinetics
        self.widgetActionUpdateDropdownPointKinetics()

        # Update Diffusion Neutronics
        self.widgetActionUpdateDropdownDiffusionNeutronics()

    # On change of Point-Kinetics checkbox
    def widgetActionUpdateDropdownPointKinetics(self, *args) -> None:
        isPKShow = self.coreStateSelected.get() == self.nominalKeyWord
        isPKShow = isPKShow and utils.binaryToBool(self.isPointKineticsExtract.get())
        self.pointKineticsParameterList.displayParameters(isPKShow)

    # On change of Diffusion Neutronics checkbox
    def widgetActionUpdateDropdownDiffusionNeutronics(self, *args) -> None:
        isDNShow = self.coreStateSelected.get() == self.nominalKeyWord
        isDNShow = isDNShow and utils.binaryToBool(self.isDiffusionNeutronicsExtract.get())
        self.diffusionNeutronicsParameterList.displayParameters(isDNShow)


    #--------------------------------------------------------------------------*
    # Writer

    def writeOutputFile(self, outputFilename: str="", coreState: str="") -> bool:
        """
        return success of Serpent file extraction and writing
        """
        if (utils.isInputFileNameDefined(outputFilename, isCreate=True)):
            print("Opening output file: {}\n".format(outputFilename))
            outputFile = open(outputFilename, 'w')

            # Header
            outputFile.write(
                "/*--------------------------------*- C++ -*----------------------------------*\ \n"+
                "|   {:73} |\n".format("Cross-section dictionary")+
                "|   {:73} |\n".format("Generated for GeN-Foam by SerpentToFoamXS Version {}".format(self.version))+
                "|   {:73} |\n".format("Date: {}".format(datetime.date.today().strftime("%d/%m/%Y")))+
                "|   {:73} |\n".format("From SERPENT results file: {}".format(self.serpentFileName.get().split("/")[-1]))+
                "\*--------------------------------------------------------------------------*/\n"
            )
            # Foam header
            outputFile.write(
                "FoamFile\n"+
                "{\n"+
                "    version     2.0;\n"+
                "    format      ascii;\n"+
                "    class       dictionary;\n"+
                "    location    \"constant/neutroRegion\";\n"+
                "    object      {};\n".format(outputFilename)+
                "}\n"
            )
            outputFile.write("// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //\n\n")

            # Special case for reactorState file
            if (coreState == self.reactorState.name):
                self.reactorState.writeParameters(outputFile, isWriteForNuclearData=True)
                # Extract keff values
                extraction.keff(outputFile, self.serpentFileName.get())
                outputFile.write("// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //\n")
                return(True)

            # Delayed neutron fraction (zero/eff)
            print(self.getDelayedNeutronFractionMessage())
            outputFile.write("//- "+self.getDelayedNeutronFractionMessage()+"\n")

            # Fast Neutrons
            # outputFile.write("\n//- Relevant for both Point-Kinetics and Diffusion Neutronics\n\n")

            #------------------------------------------------------------------*
            # Reference and Perturbed values (only core state selected by the user)
            coreStateSelected = coreState if (coreState != "") else self.coreStateSelected.get()
            for state in self.coreStateList:
                if (state.name == coreStateSelected):
                    state.writeParameters(outputFile, isWriteForNuclearData=True)

            #------------------------------------------------------------------*
            # Print Serpent file path used for data extraction
            utils.printFileName(self.serpentFileName.get(), "Serpent file selected:")

            #------------------------------------------------------------------*
            # Diffusion Neutronics
            outputFile.write("\n//- Parameters relevant only for Diffusion Neutronics\n\n")
            # Energy and Precursors Groups
            nEnergyGroups, nPrecursorGroups = extraction.groups(outputFile, self.serpentFileName.get())

            #------------------------------------------------------------------*
            # Extraction from Serpent File if it exist AND core state nominal
            if (utils.isInputFileNameDefined(self.serpentFileName.get())):
                if (coreStateSelected == self.nominalKeyWord):
                    # Diffusion Neutronics
                    if (self.isDiffusionNeutronicsExtract.get() == 1):
                        self.diffusionNeutronicsParameterList.writeParameters(outputFile, isWriteForNuclearData=True)

                    # Point-Kinetics
                    if (self.isPointKineticsExtract.get() == 1):
                        outputFile.write("\n//- Parameters relevant only for Point-Kinetics\n\n")

                        # Point-Kinetics Beta and Lambda
                        isPKExtractedCorrectly = extraction.pointKinetics(
                            outputFile,
                            self.serpentFileName.get(),
                            zeroeff=self.zeroeffSelected.get()
                        )

                        # Point-Kinetics feedback coefficients
                        self.pointKineticsParameterList.writeParameters(outputFile, isWriteForNuclearData=True)

                        if (not isPKExtractedCorrectly):
                            self.printOutputTextLabel("Point-Kinetics not extracted")
                            return(False)
            else:
                self.printOutputTextLabel("No Serpent file name provided")
                return(False)

            #------------------------------------------------------------------*
            # Universes
            outputFile.write("\n//- Foam zones from Serpent universes\n\n")
            serpentUnivToFoam = []
            for line in self.serpentUnivToFoamText.get(1.0, 'end').split("\n"):
                serpentUnivToFoam.append(line.split())
            isXSExtractedCorrectly = extraction.crossSections(
                outputFile=outputFile,
                serpentFileName=self.serpentFileName.get(),
                serpentUnivToFoam=serpentUnivToFoam,
                zeroeff=self.zeroeffSelected.get()
            )

            if (not isXSExtractedCorrectly):
                self.printOutputTextLabel("XS not extracted")

            # End line
            outputFile.write("// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //\n")

            # Close file
            outputFile.close()
            return(True)
        else:
            self.printOutputTextLabel("No output file name provided")
            return(False)


    def writeResults(self) -> None:
        isExtractionSuccessful = False

        # Writing based on core state selected by the user
        coreStateSelected = self.coreStateSelected.get()
        if (coreStateSelected in [cs.name for cs in self.coreStateList]):
            for state in self.coreStateList:
                if (coreStateSelected == state.name or (self.isBashMode and self.isExtractAll)):
                    self.printOutputTextLabel("\nExtraction for core state: "+state.name)
                    isExtractionSuccessful = self.writeOutputFile(state.fileName, state.name)
        else:
            self.printOutputTextLabel("Select a core state")

        # Extraction success message
        if (isExtractionSuccessful):
            self.writeUserInputFile()
            self.printOutputTextLabel("Extraction completed")


    def writeUserInputFile(self) -> None:
        # If no user input file is provided, the code create a new one with the
        # Serpent output file name
        defaultUserInputFileName = self.userInputFileName.get()
        if (not utils.isInputFileNameDefined(self.userInputFileName.get())):
            if (utils.isInputFileNameDefined(self.serpentFileName.get())):
                defaultUserInputFileName = self.serpentFileName.get().split("/")[-1]+".in"
            else:
                defaultUserInputFileName = "userInputSFXS.in"

        # Open
        userInputFile = open(defaultUserInputFileName, "w")

        utils.writeScalarLine(userInputFile, "serpentFilePath", self.serpentFileName.get())
        utils.writeScalarLine(userInputFile, "coreState", self.coreStateSelected.get())
        utils.writeScalarLine(userInputFile, "zeroeff", self.zeroeffSelected.get())

        utils.writeScalarLine(userInputFile, "isPointKineticsExtract", self.isPointKineticsExtract.get(), isBoolValue=True)
        utils.writeScalarLine(userInputFile, "isDiffusionNeutronicsExtract", self.isDiffusionNeutronicsExtract.get(), isBoolValue=True)

        # Loop over the Point-Kinetics variables
        userInputFile.write("\n// Point-Kinetics variables\n")
        self.pointKineticsParameterList.writeParameters(userInputFile, isWriteForNuclearData=False)

        # Loop over the Diffusion-Neutronics variables
        userInputFile.write("\n// Diffusion-Neutronics variables\n")
        self.diffusionNeutronicsParameterList.writeParameters(userInputFile, isWriteForNuclearData=False)

        # Loop over all the core state variables
        for state in self.coreStateList:
            userInputFile.write("\n// {}\n".format(state.name))
            state.writeParameters(userInputFile, isWriteForNuclearData=False)

        # Write the content of the universe textbox
        userInputFile.write("\n// Universes <serpentName> <foamName> <fuelFrac>\n")
        for line in self.serpentUnivToFoamText.get(1.0, 'end').split("\n"):
            if (len(line) > 1):
                userInputFile.write("serpentUnivToFoam {}\n".format(line))

        # Close
        userInputFile.close()

        # Refresh the data for the GUI
        self.userInputFileName.set(defaultUserInputFileName)

        # Print successful save
        self.printOutputTextLabel("User input file saved")
