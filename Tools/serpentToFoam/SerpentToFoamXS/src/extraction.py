"""
SerpentToFoamXS

Generic functions for test and data extraction from the Serpent output file
_res.m

Author: Thomas Guilbaud, EPFL/Transmutex SA
Last Update: 25/03/2022
"""

#------------------------------------------------------------------------------*
# Imports

from . import utils

#------------------------------------------------------------------------------*
# Keff

def keff(outputFile, serpentFileName: str="") -> None:
    class keffParameter:
        def __init__(self, serpentName: str):
            self.serpentName = serpentName
            self.mean        = 0
            self.error       = 1
            self.isExtracted = False

    if (utils.isInputFileNameDefined(serpentFileName)):
        inputFile = open(serpentFileName, 'r')

        outputFile.write("//- keff values using different estimators (relative uncertainty in unit)\n")
        outputFile.write("//  keff has no effect on pointKinetics\n\n")

        # List of keff estimators implemented in Serpent
        keffList = [
            keffParameter("ANA_KEFF"), # ANA_KEFF 	6(2) 	Analog estimate of keff: total, prompt and delayed neutron contribution.
            keffParameter("IMP_KEFF"), # IMP_KEFF 	2 	Implicit estimate of keff.
            keffParameter("COL_KEFF"), # COL_KEFF 	2 	Collision estimate of keff.
            keffParameter("ABS_KEFF")  # ABS_KEFF 	2 	Absorption estimate of keff.
        ]
        for line in inputFile.readlines():
            line = line.split()
            if (len(line) > 0):
                for k in keffList:
                    if (line[0] == k.serpentName):
                        k.mean, k.error = float(line[6]), float(line[7])
                        k.isExtracted = True

            if (all([k.isExtracted for k in keffList])):
                # Write all keff values find and uncomment the estimator with
                # the minimum uncertainty
                minError = min([k.error for k in keffList])
                isFirstToBeWritten = True
                for k in keffList:
                    isFirstMinValue = k.error == minError and isFirstToBeWritten
                    utils.writeScalarLine(
                        outputFile, "keff", k.mean,
                        comment="({}), from {}".format(k.error, k.serpentName),
                        isWriteForNuclearData=True, isCommented=not isFirstMinValue
                    )
                    if (isFirstMinValue):
                        isFirstToBeWritten = False
                break

#------------------------------------------------------------------------------*
# Groups

def groups(outputFile, serpentFileName: str=""):
    """
    Extract the energy group and the precursor group
    """
    numberOfEnergyGroups, numberOfPrecursorGroups = 0, 0

    if (utils.isInputFileNameDefined(serpentFileName)):
        inputFile = open(serpentFileName, 'r')

        isNEnergyGroupExtracted, isNPrecursorsGroupExtracted = False, False
        for line in inputFile.readlines():
            line = line.split()
            if (len(line) > 0):
                # Number of energy groups
                if (line[0] == "MACRO_NG"):
                    numberOfEnergyGroups = int(line[4])
                    print("Number of energy groups: {}\n".format(numberOfEnergyGroups))
                    utils.writeScalarLine(
                        outputFile, "energyGroups", numberOfEnergyGroups,
                        isWriteForNuclearData=True, comment="from MACRO_NG"
                    )
                    isNEnergyGroupExtracted = True
                # Number of precursor groups
                elif (line[0] == "FWD_ANA_BETA_ZERO"): # FWD_ANA_BETA_ZERO
                    numberOfPrecursorGroups = len(line[6:-2:2])-1 # -1 to remove the total value
                    print("Number of delayed neutron precursor group: {}\n".format(numberOfPrecursorGroups))
                    utils.writeScalarLine(
                        outputFile, "precGroups", numberOfPrecursorGroups,
                        isWriteForNuclearData=True, comment="from FWD_ANA_BETA_ZERO"
                    )
                    isNPrecursorsGroupExtracted = True
            # Stop when data extracted once
            if (isNEnergyGroupExtracted and isNPrecursorsGroupExtracted):
                break

        inputFile.close()

    else:
        print("No Serpent file name provided")

    return(numberOfEnergyGroups, numberOfPrecursorGroups)

#------------------------------------------------------------------------------*
# Point-Kinetics Parameters Extraction

def getPKKeyWord(zeroeff: str="",
    isBeta: bool=False, isLambda: bool=False, isPromptGenerationTime: bool=False
) -> str:
    if (isBeta):
        if (zeroeff == "effective"):
            return("ADJ_IFP_IMP_BETA_EFF")
        elif (zeroeff == "physical"):
            return("FWD_ANA_BETA_ZERO")
    elif (isLambda):
        return("FWD_ANA_LAMBDA")
    elif (isPromptGenerationTime):
        return("ADJ_IFP_GEN_TIME")

def pointKinetics(
    outputFile,
    serpentFileName: str="",
    zeroeff: str=""
) -> bool:
    """
    Extract the Point-Kinetics parameters.
    return the status of the extraction. True: success; False: fail
    """
    print("Point-Kinetics parameters extraction ...")

    if (utils.isInputFileNameDefined(serpentFileName)):
        inputFile = open(serpentFileName, 'r')

        isBetaExtracted, isLambdaExtracted, isPertGenTime = False, False, False
        textBeta, textLambda = "", ""
        for line in inputFile.readlines():
            line = line.split()
            if (len(line) > 0):
                # Beta delayed neutron fraction
                if (line[0] == getPKKeyWord(zeroeff, isBeta=True)
                    and not isBetaExtracted
                ):
                    textBeta += "Beta ( // from {}\n".format(getPKKeyWord(zeroeff, isBeta=True))
                    for beta_i in line[6+2:-2:2]: # 6+2 to remove the total (first value)
                        textBeta += "\t{:.6e}\n".format(float(beta_i))
                    textBeta += ");\n\n"
                    isBetaExtracted = True
                # Lambda decay
                elif (line[0] == getPKKeyWord(isLambda=True) and not isLambdaExtracted):
                    textLambda += "lambda ( // from {}\n".format(getPKKeyWord(isLambda=True))
                    for lambda_i in line[6+2:-2:2]: # 6+2 to remove the total (first value)
                        textLambda += "\t{:.6e}\n".format(float(lambda_i))
                    textLambda += ");\n\n"
                    isLambdaExtracted = True
                # Neutron Generation time
                elif (line[0] == getPKKeyWord(isPromptGenerationTime=True) and not isPertGenTime):
                    utils.writeScalarLine(
                        outputFile, "promptGenerationTime", line[8], "s",
                        isWriteForNuclearData=True, comment="from "+getPKKeyWord(isPromptGenerationTime=True)
                    )
                    isPertGenTime = True
            # Stop when data extracted once
            if (isBetaExtracted and isLambdaExtracted and isPertGenTime):
                outputFile.write(textBeta)
                outputFile.write(textLambda)
                break

        inputFile.close()

        print("OK\n")
        return(True)
    else:
        print("No Serpent file selected\n")
        return(False)

#------------------------------------------------------------------------------*
# Extract Cross Section

def crossSections(
    outputFile,
    serpentFileName: str="",
    serpentUnivToFoam: list=[],
    zeroeff: str=""
) -> bool:
    """
    Extract the cross-sections and the multi-group variables.
    return the status of the extraction. True: success; False: fail
    """

    # Unit conversions
    cm, MeV = 1, 1
    m = 1e2*cm
    j = MeV/1.602176487e-13

    print("Cross-sections extraction ...")

    if (utils.isInputFileNameDefined(serpentFileName)):
        inputFile = open(serpentFileName, 'r')

        # List of the Serpent universe's name used in the simulation
        serpentUniverseList = []
        # Storing variables of the Serpent output file
        numberOfEnergyGroups = 0
        INF_INVV, INF_DIFFCOEF, INF_NSF, INF_FISS, FISSE = {}, {}, {}, {}, {}
        INF_SP0, INF_SP1, INF_SP2, INF_SP3, INF_SP4, INF_SP5, INF_TOT = {}, {}, {}, {}, {}, {}, {}
        INF_CHIP, INF_CHID, FWD_ANA_BETA_ZERO, ADJ_IFP_IMP_BETA_EFF = {}, {}, {}, {}
        FWD_ANA_LAMBDA, INF_FLX = {}, {}
        # Values that come before the entry "GC_UNIVERSE_NAME"
        FISSE_temp, FWD_ANA_BETA_ZERO_temp, FWD_ANA_LAMBDA_temp, ADJ_IFP_IMP_BETA_EFF_temp = 0, [], [], []

        # Serpent values extraction
        for line in inputFile.readlines():
            line = line.split()
            if (len(line) > 0):
                if (line[0] == "GC_UNIVERSE_NAME"):
                    serpentUniverseList.append(line[5].split("'")[1])
                    FISSE[serpentUniverseList[-1]] = FISSE_temp
                    FWD_ANA_BETA_ZERO[serpentUniverseList[-1]] = FWD_ANA_BETA_ZERO_temp
                    FWD_ANA_LAMBDA[serpentUniverseList[-1]] = FWD_ANA_LAMBDA_temp
                    ADJ_IFP_IMP_BETA_EFF[serpentUniverseList[-1]] = ADJ_IFP_IMP_BETA_EFF_temp
                elif (line[0] == "MACRO_NG"):
                    numberOfEnergyGroups = int(line[4])
                elif (line[0] == "INF_INVV"):
                    INF_INVV[serpentUniverseList[-1]] = [float(v)*m for v in line[6:-2:2]] # Conversion from cm to m
                elif (line[0] == "INF_DIFFCOEF"):
                    INF_DIFFCOEF[serpentUniverseList[-1]] = [float(v)/m for v in line[6:-2:2]] # Conversion from cm to m
                elif (line[0] == "INF_NSF"):
                    INF_NSF[serpentUniverseList[-1]] = [float(v)*m for v in line[6:-2:2]] # Conversion from cm to m
                elif (line[0] == "INF_FISS"):
                    INF_FISS[serpentUniverseList[-1]] = [float(v)*m for v in line[6:-2:2]] # Conversion from cm to m
                elif (line[0] == "FISSE"):
                    FISSE_temp = float(line[6])/j # Conversion from MeV to J
                elif (line[0] == "INF_SP0"):
                    INF_SP0[serpentUniverseList[-1]] = [float(v)*m for v in line[6:-2:2]] # Conversion from cm to m
                elif (line[0] == "INF_SP1"):
                    INF_SP1[serpentUniverseList[-1]] = [float(v)*m for v in line[6:-2:2]] # Conversion from cm to m
                elif (line[0] == "INF_SP2"):
                    INF_SP2[serpentUniverseList[-1]] = [float(v)*m for v in line[6:-2:2]] # Conversion from cm to m
                elif (line[0] == "INF_SP3"):
                    INF_SP3[serpentUniverseList[-1]] = [float(v)*m for v in line[6:-2:2]] # Conversion from cm to m
                elif (line[0] == "INF_SP4"):
                    INF_SP4[serpentUniverseList[-1]] = [float(v)*m for v in line[6:-2:2]] # Conversion from cm to m
                elif (line[0] == "INF_SP5"):
                    INF_SP5[serpentUniverseList[-1]] = [float(v)*m for v in line[6:-2:2]] # Conversion from cm to m
                elif (line[0] == "INF_TOT"):
                    INF_TOT[serpentUniverseList[-1]] = [float(v)*m for v in line[6:-2:2]] # Conversion from cm to m

                elif (line[0] == "INF_CHIP"):
                    INF_CHIP[serpentUniverseList[-1]] = [float(v) for v in line[6:-2:2]]
                elif (line[0] == "INF_CHID"):
                    INF_CHID[serpentUniverseList[-1]] = [float(v) for v in line[6:-2:2]]
                elif (line[0] == "FWD_ANA_BETA_ZERO"):
                    FWD_ANA_BETA_ZERO_temp = [float(v) for v in line[6+2:-2:2]] # Avoid first because total
                elif (line[0] == "ADJ_IFP_IMP_BETA_EFF"):
                    ADJ_IFP_IMP_BETA_EFF_temp = [float(v) for v in line[6+2:-2:2]] # Avoid first because total
                elif (line[0] == "FWD_ANA_LAMBDA"):
                    FWD_ANA_LAMBDA_temp = [float(v) for v in line[6+2:-2:2]] # Avoid first because total
                elif (line[0] == "INF_FLX"):
                    INF_FLX[serpentUniverseList[-1]] = [float(v) for v in line[6:-2:2]]

        # All the necessary values have been extracted and the Serpent output
        # file is no more needed
        inputFile.close()

        # Writing in the output file for Foam
        outputFile.write("zones\n(\n")

        # Loop over the conversion asked by the user
        for sentence in serpentUnivToFoam:
            if (len(sentence) > 0):
                serpentUniv, foamCell, fuelFraction = sentence[0], sentence[1], sentence[2]

                # Check is the universe provided by the user exist
                isUniverseFoundInSerpent = False
                for idx, univFromSerpentFile in enumerate(serpentUniverseList):
                    if (serpentUniv == univFromSerpentFile):
                        isUniverseFoundInSerpent = True

                if (isUniverseFoundInSerpent):
                    outputFile.write("\t"+foamCell+"  //- Serpent Universe: "+serpentUniv+"\n\t{\n")

                    # Fuel Fraction
                    outputFile.write("\t\t{:18} {};\n".format("fuelFraction", fuelFraction))

                    # 1/V -> IV
                    utils.writeSimpleList(outputFile, "IV", INF_INVV[serpentUniv])

                    # Diffusion coefficient -> D
                    utils.writeSimpleList(outputFile, "D", INF_DIFFCOEF[serpentUniv])

                    # nuSigmaEff -> A
                    utils.writeSimpleList(outputFile, "nuSigmaEff", INF_NSF[serpentUniv])

                    # sigmaPow -> A
                    utils.writeSimpleList(
                        outputFile, "sigmaPow",
                        [value*FISSE[serpentUniv] for value in INF_FISS[serpentUniv]]
                    )

                    # scattering matrix P0 -> MS
                    utils.writeMatrix(outputFile, "scatteringMatrixP0", INF_SP0[serpentUniv])

                    # scattering matrix P1 -> MS1
                    utils.writeMatrix(outputFile, "scatteringMatrixP1", INF_SP1[serpentUniv])

                    # scattering matrix P2 -> MS2
                    utils.writeMatrix(outputFile, "scatteringMatrixP2", INF_SP2[serpentUniv])

                    # scattering matrix P3 -> MS3
                    utils.writeMatrix(outputFile, "scatteringMatrixP3", INF_SP3[serpentUniv])

                    # scattering matrix P4 -> MS4
                    utils.writeMatrix(outputFile, "scatteringMatrixP4", INF_SP4[serpentUniv])

                    # scattering matrix P5 -> MS5
                    utils.writeMatrix(outputFile, "scatteringMatrixP5", INF_SP5[serpentUniv])
                    """

                    fprintf(fid,"\n  nonuniform List<scalar> %i (",ng)

                    DISAPP = zeros(ng,1)
                    for i = 1:ng
                            DISAPP(i) = INF_TOT(idx,(i*2)-1) - MS(i,i)
                            fprintf(fid,"%.6e ",DISAPP(i)/cm2m)
                    end
                    fprintf(fid," );")

                    """
                    # sigma disappearence (abs+ capture + group transfer below) -> DISAPP
                    # sigmaDisapp[i] = INF_TOT[i] - INF_SP0[i,i]
                    utils.writeSimpleList(
                        outputFile, "sigmaDisapp",
                        [v - INF_SP0[serpentUniv][i*(numberOfEnergyGroups+1)] for i, v in enumerate(INF_TOT[serpentUniv])]
                    )

                    # Prompt Neutron Spectrum -> XP
                    utils.writeSimpleList(outputFile, "chiPrompt", INF_CHIP[serpentUniv])

                    if (zeroeff == "effective"):
                        # Delayed Neutron Spectrum -> XD
                        utils.writeSimpleList(outputFile, "chiDelayed", INF_CHIP[serpentUniv])
                        # Beta (delayed neutron fraction) -> BE(eff), BZ(zero)
                        utils.writeSimpleList(outputFile, "Beta", ADJ_IFP_IMP_BETA_EFF[serpentUniv])
                    elif (zeroeff == "physical"):
                        # Delayed Neutron Spectrum -> XD
                        utils.writeSimpleList(outputFile, "chiDelayed", INF_CHID[serpentUniv])
                        # Beta (delayed neutron fraction) -> BE(eff), BZ(zero)
                        utils.writeSimpleList(outputFile, "Beta", FWD_ANA_BETA_ZERO[serpentUniv])

                    # Decay constant (lambda) -> LAM
                    utils.writeSimpleList(outputFile, "lambda", FWD_ANA_LAMBDA[serpentUniv])

                    # Discontinuity factors
                    utils.writeSimpleList(
                        outputFile, "discFactor",
                        [1]*numberOfEnergyGroups
                    )

                    # Integral Fluxes -> integralFlux
                    utils.writeSimpleList(
                        outputFile, "integralFlux",
                        [
                            value/INF_FLX["0"][i] if (INF_FLX["0"][i] != 0) else 1
                            for i, value in enumerate(INF_FLX[serpentUniv])
                        ]
                    )

                    # End material data extraction
                    outputFile.write("\t}\n\n")

                else:
                    print("The Serpent universe {} doesn't exist".format(serpentUniv))

        outputFile.write(");\n\n")

        print("OK\n")
        return(True)
    else:
        print("No Serpent file selected\n")
        return(False)
