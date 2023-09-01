"""
Author: Thomas Guilbaud, 2022/10/13

Script to generate time dependant core parameters graphs.

Comment: The power are well computed by the integral of the powerDensity in the
    neutroRegion and fluidRegion and using int(A h (T_clad - T_fluid) dV).
    For some reason, there is a problen with P=m cp dT + m (vo²-vi²)/2.
    Even with changing to another form there is the same result
    (m gamma/(gamma-1)(Ps/rhos-Pe/rhoe)).
"""


#==============================================================================*
#                                   Imports
#==============================================================================*

import sys
import numpy as np
import matplotlib.pyplot as plt
from hydrogen import specificHeatCapacity
# from unitsOpenMC import *
from unitsSI import *


#==============================================================================*
#                              Usefull Functions
#==============================================================================*

def addVectors(vec1: list, vec2: list) -> list:
    return([v1+v2 for v1, v2 in zip(vec1, vec2)])

def divideVectors(vec1: list, vec2: list) -> list:
    return([v1/v2 for v1, v2 in zip(vec1, vec2)])

def appendByLine(list: list, pattern: str, line: str, pos: int=4) -> list:
    if (list == None):
        list = []
    if (pattern in line):
        list.append(float(line.split()[pos].split(')')[0]))
    return(list)

def powerFromKinetic(mdot: float, vi: float, vo: float) -> float:
    """
        return the power from the kinetic energy
    """
    return(0.5 * mdot * (vo**2 - vi**2))

def power(mdot, cp, Ti, To, vi=0, vo=0) -> float:
    """
        Power from the first law of thermodynamics of open systems
        return(mdot * cp * (To - Ti) + powerFromKinetic(mdot, vi, vo))

        Power for compressible fluids: Q = mdot * int_T cp(T) dT
    """
    deltaT = To-Ti
    Tavg = (To+Ti)/2
    Z = 1.32
    dz = 1.32/1000
    T = lambda z: Tavg - deltaT/2 * np.cos(np.pi/Z * z)
    power = 0
    for z in np.arange(0, Z, dz):
        power += cp(T(z+dz/2)) * (T(z+dz)-T(z))
    return(power * mdot)


#==============================================================================*

cp = lambda T: specificHeatCapacity(T, isFromFit=True) # J/(kg.K)

nFuelAssembly = 1500/6 # Number of fuel assembly in core
nominalPower = 937*MW/nFuelAssembly
fuelElementStructureFraction = 0.696911
SfuelElement = 0.00190587*m2 * (1-fuelElementStructureFraction) # m²
gConst = 9.8 * m/(s**2)
Rconst = 8.314 * J/mol/K# SI units
gamma = 1.4
M_H2 = 2.014 * g/mol


fig1, (axTemperature, axPower) = plt.subplots(2, 1, sharex=True)
fig2, axMassFlowrate = plt.subplots(1)

# --- Loop over files
for arg in sys.argv[1:]:
    file = open(arg, 'r')

    # Set the lists
    times = []
    inletFuelElemTemp, outletFuelElemTemp, inletFuelElemMFlow, outletFuelElemMFlow = [], [], [], []
    inletCentralUTemp, outletCentralUTemp, inletCentralUMFlow, outletCentralUMFlow = [], [], [], []
    inletFuelElemVelo, outletFuelElemVelo, inletCentralUVelo, outletCentralUVelo = [], [], [], []
    inletFuelElemPres, outletFuelElemPres, inletCentralUPres, outletCentralUPres = [], [], [], []
    inletFuelElemAlphaRhoPhi, outletFuelElemAlphaRhoPhi, inletCentralUAlphaRhoPhi, outletCentralUAlphaRhoPhi = [], [], [], []
    inletFuelElemAlphaPhi, outletFuelElemAlphaPhi, inletCentralUAlphaPhi, outletCentralUAlphaPhi = [], [], [], []
    powerIntegratedFluidRegion, powerIntegratedNeutroRegion, powerAhdeltaT = [], [], []
    timeTemp = 0

    # Loop on the lines to recover the parameters
    for line in file.readlines():
        if (line[:6] == "Time ="):
            timeTemp = float(line.split()[2])
        if ("patch inletCentralUnloaded massFlow =" in line):
            times.append(timeTemp)
            inletCentralUMFlow = appendByLine(inletCentralUMFlow, "patch inletCentralUnloaded massFlow =", line)
        outletCentralUMFlow = appendByLine(outletCentralUMFlow, "patch outletCentralUnloaded massFlow =", line)
        inletFuelElemMFlow  = appendByLine(inletFuelElemMFlow , "patch inletFuelElement massFlow =", line)
        outletFuelElemMFlow = appendByLine(outletFuelElemMFlow, "patch outletFuelElement massFlow =", line)

        inletCentralUTemp   = appendByLine(inletCentralUTemp  , "patch inletCentralUnloaded TBulk =", line)
        outletCentralUTemp  = appendByLine(outletCentralUTemp , "patch outletCentralUnloaded TBulk =", line)
        inletFuelElemTemp   = appendByLine(inletFuelElemTemp  , "patch inletFuelElement TBulk =", line)
        outletFuelElemTemp  = appendByLine(outletFuelElemTemp , "patch outletFuelElement TBulk =", line)

        inletFuelElemVelo  = appendByLine(inletFuelElemVelo , "areaAverage(inletFuelElement) of magU =", line)
        outletFuelElemVelo = appendByLine(outletFuelElemVelo, "areaAverage(outletFuelElement) of magU =", line)
        inletCentralUVelo  = appendByLine(inletCentralUVelo , "areaAverage(inletCentralUnloaded) of magU =", line)
        outletCentralUVelo = appendByLine(outletCentralUVelo, "areaAverage(outletCentralUnloaded) of magU =", line)

        inletFuelElemPres  = appendByLine(inletFuelElemPres , "areaAverage(inletFuelElement) of p =", line)
        outletFuelElemPres = appendByLine(outletFuelElemPres, "areaAverage(outletFuelElement) of p =", line)
        inletCentralUPres  = appendByLine(inletCentralUPres , "areaAverage(inletCentralUnloaded) of p =", line)
        outletCentralUPres = appendByLine(outletCentralUPres, "areaAverage(outletCentralUnloaded) of p =", line)

        inletFuelElemAlphaRhoPhi  = appendByLine(inletFuelElemAlphaRhoPhi , "areaAverage(inletFuelElement) of alphaRhoPhi =", line)
        outletFuelElemAlphaRhoPhi = appendByLine(outletFuelElemAlphaRhoPhi, "areaAverage(outletFuelElement) of alphaRhoPhi =", line)
        inletCentralUAlphaRhoPhi  = appendByLine(inletCentralUAlphaRhoPhi , "areaAverage(inletCentralUnloaded) of alphaRhoPhi =", line)
        outletCentralUAlphaRhoPhi = appendByLine(outletCentralUAlphaRhoPhi, "areaAverage(outletCentralUnloaded) of alphaRhoPhi =", line)

        inletFuelElemAlphaPhi  = appendByLine(inletFuelElemAlphaPhi , "areaAverage(inletFuelElement) of alphaPhi =", line)
        outletFuelElemAlphaPhi = appendByLine(outletFuelElemAlphaPhi, "areaAverage(outletFuelElement) of alphaPhi =", line)
        inletCentralUAlphaPhi  = appendByLine(inletCentralUAlphaPhi , "areaAverage(inletCentralUnloaded) of alphaPhi =", line)
        outletCentralUAlphaPhi = appendByLine(outletCentralUAlphaPhi, "areaAverage(outletCentralUnloaded) of alphaPhi =", line)

        powerIntegratedFluidRegion  = appendByLine(powerIntegratedFluidRegion , "volIntegrate(fluidRegion) of powerDensity.nuclearFuelPin =", line)
        powerIntegratedNeutroRegion = appendByLine(powerIntegratedNeutroRegion, "volIntegrate(neutroRegion) of powerDensity =", line)
        powerAhdeltaT = appendByLine(powerAhdeltaT, "volIntegrate(fluidRegion) of hdeltaT =", line)


    # Compute the total mass flowrate
    inletMassFlowrate = addVectors(inletFuelElemMFlow, inletCentralUMFlow)
    outletMassFlowrate = addVectors(outletFuelElemMFlow, outletCentralUMFlow)

    # Compute density
    inletFuelElemRho = divideVectors(inletFuelElemAlphaRhoPhi, inletFuelElemAlphaPhi)
    outletFuelElemRho = divideVectors(outletFuelElemAlphaRhoPhi, outletFuelElemAlphaPhi)
    inletCentralURho = divideVectors(inletCentralUAlphaRhoPhi, inletCentralUAlphaPhi)
    outletCentralURho = divideVectors(outletCentralUAlphaRhoPhi, outletCentralUAlphaPhi)

    # Compute mass flowrate from density and velocity
    fuelElementMDot = [rho*SfuelElement*u for rho, u in zip(inletFuelElemRho, inletFuelElemVelo)]

    # Compute the power in each cellZones
    powerFuelElem = [power((mdotIn+mdotOut)/2, cp, Ti, To, vi, vo)
        for mdotIn, mdotOut, Ti, To, vi, vo in zip(
            fuelElementMDot, fuelElementMDot,
            inletFuelElemTemp, outletFuelElemTemp,
            inletFuelElemVelo, outletFuelElemVelo
    )]
    powerCentralU = [power((mdotIn+mdotOut)/2, cp, Ti, To, vi, vo)
        for mdotIn, mdotOut, Ti, To, vi, vo in zip(
            inletCentralUMFlow, outletCentralUMFlow,
            inletCentralUTemp, outletCentralUTemp,
            inletCentralUVelo, outletCentralUVelo
    )]

    powerMcpdT = addVectors(powerFuelElem, powerCentralU)

    powerFromKineticLast = powerFromKinetic((inletFuelElemMFlow[-1]+outletFuelElemMFlow[-1])/2, inletFuelElemVelo[-1], outletFuelElemVelo[-1])


    # --- Plots
    # Temperature
    axTemperature.plot(times, outletFuelElemTemp, label='Fuel outlet', color='tab:blue')
    axTemperature.plot(times, outletCentralUTemp, label='Central outlet', color='tab:orange')
    axTemperature.hlines([83.3], xmin=times[0], xmax=times[-1], linestyle='--', label='Inlet', color='gray')

    # Power
    # axPower.plot(times, powerFuelElem, label='Fuel')
    # axPower.plot(times, powerCentralU, label='Central')
    axPower.plot(times, powerIntegratedFluidRegion, label='Fluid integrated', marker="o")
    axPower.plot(times, powerIntegratedNeutroRegion, label='Neutro integrated', marker="|")
    axPower.plot(times, powerAhdeltaT, label=r"$\int A h (T_{fluid}-T_{clad,outer}) dV$", marker="s")
    axPower.plot(times, powerMcpdT, label=r"$\dot{m} c_p \Delta T + \frac{\dot{m}}{2} (v_o^2 - v_i^2)$", marker="D")

    # Mass flowrate
    axMassFlowrate.plot(times, outletFuelElemMFlow, label='Fuel outlet', color='tab:blue')
    axMassFlowrate.plot(times, inletFuelElemMFlow , label='Fuel inlet', color='tab:blue', linestyle='--')
    axMassFlowrate.plot(times, outletCentralUMFlow, label='Central outlet', color='tab:orange')
    axMassFlowrate.plot(times, inletCentralUMFlow , label='Central inlet', color='tab:orange', linestyle='--')
    axMassFlowrate.plot(times, outletMassFlowrate , label='Total outlet', color='tab:green')
    axMassFlowrate.plot(times, inletMassFlowrate  , label='Total inlet', color='tab:green', linestyle='--')


    # Print
    print("Time: {} s".format(times[-1]))
    print("")
    print("Power:")
    print("    Expected          = {:8.6e} W".format(nominalPower))
    print("    Neutro integrated = {:8.6e} W ({:7.3f} %)".format(powerIntegratedNeutroRegion[-1], powerIntegratedNeutroRegion[-1]/nominalPower*100))
    print("    Fluid  integrated = {:8.6e} W ({:7.3f} %)".format(powerIntegratedFluidRegion[-1], powerIntegratedFluidRegion[-1]/nominalPower*100))
    print("    int Ah(Tf-Tco) dV = {:8.6e} W ({:7.3f} %)".format(powerAhdeltaT[-1], powerAhdeltaT[-1]/nominalPower*100))
    # print("    mcpDT+m/2(vo²-vi²)= {:8.6e} W ({:7.3f} %)".format(powerMcpdT[-1], powerMcpdT[-1]/nominalPower*100))
    print("    m int_T cp(T) dT  = {:8.6e} W ({:7.3f} %)".format(powerMcpdT[-1], powerMcpdT[-1]/nominalPower*100))
    print("    m/2 (vo² - vi²)   = {:8.6e} W ({:7.3f} %)".format(powerFromKineticLast, powerFromKineticLast/nominalPower*100))
    print("")
    print("Temperature:")
    print("    Expected         : 83.3 -> 1833.3 K")
    print("    Fuel element     : {:.1f} -> {:6.1f} K (last : {:.3f} K/s)".format(
        inletFuelElemTemp[-1], outletFuelElemTemp[-1], (outletFuelElemTemp[-1]-outletFuelElemTemp[-2])/(times[-1]-times[-2])
    ))
    print("    Central unloaded : {:.1f} -> {:6.1f} K (last : {:.3f} K/s)".format(
        inletCentralUTemp[-1], outletCentralUTemp[-1], (outletCentralUTemp[-1]-outletCentralUTemp[-2])/(times[-1]-times[-2])
    ))
    print("")
    print("Velocity:")
    print("    Fuel element     : {:.2f} -> {:6.2f} m/s".format(inletFuelElemVelo[-1], outletFuelElemVelo[-1]))
    print("    Central unloaded : {:.2f} -> {:6.2f} m/s".format(inletCentralUVelo[-1], outletCentralUVelo[-1]))
    print("")
    print("Density:")
    print("    Fuel element     : {:.2f} -> {:6.2f} kg/m3".format(inletFuelElemRho[-1], outletFuelElemRho[-1]))
    print("    Central unloaded : {:.2f} -> {:6.2f} kg/m3".format(inletCentralURho[-1], outletCentralURho[-1]))
    print("")
    print("Pressure:")
    print("    Expected         : 4.1370 -> 3.4470 MPa")
    print("    Fuel element     : {:.4f} -> {:.4f} MPa".format(inletFuelElemPres[-1]/MPa, outletFuelElemPres[-1]/MPa))
    print("    Central unloaded : {:.4f} -> {:.4f} MPa".format(inletCentralUPres[-1]/MPa, outletCentralUPres[-1]/MPa))
    print("")
    print("Mass flowrate:")
    print("    Fuel element     : {:6.3f} -> {:6.3f} kg/s".format(inletFuelElemMFlow[-1], outletFuelElemMFlow[-1]))
    print("    Central unloaded : {:6.3f} -> {:6.3f} kg/s".format(inletCentralUMFlow[-1], outletCentralUMFlow[-1]))
    print("    Total            : {:6.3f} -> {:6.3f} kg/s".format(inletMassFlowrate[-1], outletMassFlowrate[-1]))
    print("    Total (core)     : {:6.3f} -> {:6.3f} kg/s".format(nFuelAssembly*inletMassFlowrate[-1], nFuelAssembly*outletMassFlowrate[-1]))
    print("    Expected         : 31.750 kg/s")
    print("")
    print("From mass conservation:")
    print("    rho_i S_i u_i = rho_o S_o u_o => {:.3f} = {:.3f} kg/s".format(
        inletFuelElemRho[-1] * SfuelElement * inletFuelElemVelo[-1],
        outletFuelElemRho[-1] * SfuelElement * outletFuelElemVelo[-1] # Here S_i = S_o
    ))
    print("")
    print("From Perfect Gas law:")
    print("    r = P_i/(rho_i T_i) = P_o/(rho_o T_o) => {:.3f} = {:.3f} m²/s²/K".format(
        inletFuelElemPres[-1]/(inletFuelElemRho[-1]*inletFuelElemTemp[-1]),
        outletFuelElemPres[-1]/(outletFuelElemRho[-1]*outletFuelElemTemp[-1])
    ))
    print("")

    print("Rocketery (to-do: need to further understand, not sure on value):")
    print("    From Rocket propulsion elements 8th edition")

    outletNozzlePres = 0# 0.1*MPa
    nozzleMinArea = np.pi * (8.720 * inch/m /2)**2
    nozzleOutletArea = np.pi * (30.370 * inch/m /2)**2

    gm1og = (gamma-1)/gamma
    nozzleExhaustVelocity = np.sqrt(
        2/gm1og * Rconst*outletFuelElemTemp[-1]/M_H2 * (
            1 - (outletNozzlePres/outletFuelElemPres[-1])**gm1og
        )
        # + outletFuelElemVelo[-1]**2
    )

    # thrustCore = nFuelAssembly * SfuelElement * outletFuelElemRho[-1] * outletFuelElemPres[-1]
    # thrustCore = nFuelAssembly * outletMassFlowrate[-1] * outletFuelElemVelo[-1]
    # IspCore = thrustCore/(nFuelAssembly * inletMassFlowrate[-1] * gConst)

    IspNozzle = nozzleExhaustVelocity/gConst
    thrustNozzle = nFuelAssembly * inletMassFlowrate[-1] * nozzleExhaustVelocity

    # effectiveExhaustVelocity = IspCore * gConst
    effectiveExhaustVelocity = nozzleExhaustVelocity + outletNozzlePres * nozzleOutletArea/(nFuelAssembly*inletMassFlowrate[-1])

    # thrustCore2 = effectiveExhaustVelocity2*inletMassFlowrate[-1]*nFuelAssembly

    characteristicVelocity = outletFuelElemPres[-1] * nozzleMinArea / (nFuelAssembly*inletMassFlowrate[-1])

    # print("    Fuel element")
    # print("        thrust : {:.1f} N (core (S*rho*P1))".format(thrustCore/N))
    # print("        thrust : {:.1f} N (core (m*c))".format(nFuelAssembly*inletMassFlowrate[-1]*effectiveExhaustVelocity /N))
    # print("        Isp    : {:.1f} s (core (F/(m*g)))".format(IspCore/s))
    # print("        c      : {:.2f} m/s (Effective exhaust velocity (Isp*g))".format(
    #     effectiveExhaustVelocity/(m/s)
    # ))
    # print("        c      : {:.2f} m/s (Effective exhaust velocity (v2 + p2*A2/m))".format(
    #     effectiveExhaustVelocity2/(m/s)
    # ))
    # print("    Nozzle")
    # print("        p1/p2  : {:.3f}".format(outletFuelElemPres[-1]/outletNozzlePres))
    print("    v2     : {:.2f} m/s (Exhaust velocity)".format(
        nozzleExhaustVelocity/(m/s)
    ))
    print("    c      : {:.2f} m/s (Effective exhaust velocity (Isp*g))".format(
        effectiveExhaustVelocity/(m/s)
    ))
    print("    c*     : {:.2f} m/s (Characteristic velocity)".format(
        characteristicVelocity/(m/s)
    ))
    print("    thrust : {:.1f} N (m*c)".format(thrustNozzle/N))
    print("    Isp    : {:.1f} s (v2/g)".format(IspNozzle/s))
    print("")

    # plt.figure()
    # plt.plot(times, [1-mcpdt/pow for mcpdt, pow in zip(powerMcpdT, powerIntegratedFluidRegion)])
    # plt.grid(True)

# --- Plotting settings
# Temperature
axTemperature.set_ylabel('Temperature [K]')
axTemperature.grid(True)
axTemperature.legend()

# Power
axPower.hlines(nominalPower, xmin=times[0], xmax=times[-1], linestyle='--', color='gray', label='Nominal')
axPower.set_xlabel('Time [s]')
axPower.set_ylabel('Power [W]')
axPower.grid(True)
axPower.legend()


# Mass flowrate
axMassFlowrate.set_ylabel('Mass flowrate [kg/s]')
axMassFlowrate.grid(True)
axMassFlowrate.legend()

fig1.savefig('steadyState_temperatureAndPower.png')
fig2.savefig('steadyState_massFlowrate.png')

# plt.show()

#==============================================================================*
