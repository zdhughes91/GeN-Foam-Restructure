"""
Script to analyse the last time in a log.GeN-Foam file.

Author: Thomas Guilbaud, 03/12/2022
"""

from numpy import pi, sqrt
import os
import sys

rootPath = os.path.abspath(__file__)
for i in range(1, 4):
    sys.path.append("/".join(rootPath.split("/")[:-i]))


from unitsSI import *


class SuperDict:
    def __init__(self) -> None:
        self.dict = {}

    def __getitem__(self, key):
        return(self.dict[key])

    def __setitem__(self, key, value) -> None:
        self.dict[key] = value

    def any(self) -> bool:
        return(any(self.dict))

    def computeDensity(self) -> None:
        keys = self.dict.keys()
        if ('alphaRhoPhi' in keys and 'alphaPhi' in keys):
            self['rho'] = self['alphaRhoPhi'] / self['alphaPhi']

    def findInLine(self, line: str, pattern: str, key: str=None, pos: int=4) -> None:
        if (pattern in line):
            if (not key):
                key = pattern.split()[-1]
            self[key] = float(line.split()[pos])

# User information on usage
if (len(sys.argv) != 2):
    print("Usage: python analysis.py path/to/log.GeN-Foam")
    sys.exit(0)

# Open the log file
filename = sys.argv[1]
logfile = open(filename, 'r')

# Data
time, executionTime, clockTime = 0, 0, 0
powers = SuperDict()
inlet,        outlet        = SuperDict(), SuperDict()
inletFuel,    outletFuel    = SuperDict(), SuperDict()
inletCentral, outletCentral = SuperDict(), SuperDict()

# Read the file in reverse
for line in logfile.readlines()[::-1]:
    powers.findInLine(line, "volIntegrate(fluidRegion) of powerDensity.nuclearFuelPin", "fluid")
    powers.findInLine(line, "volIntegrate(neutroRegion) of powerDensity", "neutro")
    powers.findInLine(line, "volIntegrate(fluidRegion) of hdeltaT", "AhdeltaT")


    inlet.findInLine(line, "patch inlet massFlow")
    inlet.findInLine(line, "patch inlet massFlow", 'S', pos=7)
    inlet.findInLine(line, "areaAverage(inlet) of magU", 'U')
    inlet.findInLine(line, "areaAverage(inlet) of p")
    inlet.findInLine(line, "areaAverage(inlet) of T")
    inlet.findInLine(line, "areaAverage(inlet) of alphaRhoPhi")
    inlet.findInLine(line, "areaAverage(inlet) of alphaPhi")
    inlet.findInLine(line, "patch inlet TBulk")
    
    outlet.findInLine(line, "patch outlet massFlow")
    outlet.findInLine(line, "patch outlet massFlow", 'S', pos=7)
    outlet.findInLine(line, "areaAverage(outlet) of magU", 'U')
    outlet.findInLine(line, "areaAverage(outlet) of p")
    outlet.findInLine(line, "areaAverage(outlet) of T")
    outlet.findInLine(line, "areaAverage(outlet) of alphaRhoPhi")
    outlet.findInLine(line, "areaAverage(outlet) of alphaPhi")
    outlet.findInLine(line, "patch outlet TBulk")

    
    inletFuel.findInLine(line, "faceZone inletFuelElement massFlow")
    inletFuel.findInLine(line, "faceZone inletFuelElement massFlow", 'S', pos=7)
    inletFuel.findInLine(line, "areaAverage(inletFuelElement) of magU", 'U')
    inletFuel.findInLine(line, "areaAverage(inletFuelElement) of p")
    inletFuel.findInLine(line, "areaAverage(inletFuelElement) of T")
    inletFuel.findInLine(line, "areaAverage(inletFuelElement) of alphaRhoPhi")
    inletFuel.findInLine(line, "areaAverage(inletFuelElement) of alphaPhi")
    inletFuel.findInLine(line, "faceZone inletFuelElement TBulk")
    
    outletFuel.findInLine(line, "faceZone outletFuelElement massFlow")
    outletFuel.findInLine(line, "faceZone outletFuelElement massFlow", 'S', pos=7)
    outletFuel.findInLine(line, "areaAverage(outletFuelElement) of magU", 'U')
    outletFuel.findInLine(line, "areaAverage(outletFuelElement) of p")
    outletFuel.findInLine(line, "areaAverage(outletFuelElement) of T")
    outletFuel.findInLine(line, "areaAverage(outletFuelElement) of alphaRhoPhi")
    outletFuel.findInLine(line, "areaAverage(outletFuelElement) of alphaPhi")
    outletFuel.findInLine(line, "faceZone outletFuelElement TBulk")

    
    inletCentral.findInLine(line, "faceZone inletCentralUnloaded massFlow")
    inletCentral.findInLine(line, "faceZone inletCentralUnloaded massFlow", 'S', pos=7)
    inletCentral.findInLine(line, "areaAverage(inletCentralUnloaded) of magU", 'U')
    inletCentral.findInLine(line, "areaAverage(inletCentralUnloaded) of p")
    inletCentral.findInLine(line, "areaAverage(inletCentralUnloaded) of T")
    inletCentral.findInLine(line, "areaAverage(inletCentralUnloaded) of alphaRhoPhi")
    inletCentral.findInLine(line, "areaAverage(inletCentralUnloaded) of alphaPhi")
    inletCentral.findInLine(line, "faceZone inletCentralUnloaded TBulk")
    
    outletCentral.findInLine(line, "faceZone outletCentralUnloaded massFlow")
    outletCentral.findInLine(line, "faceZone outletCentralUnloaded massFlow", 'S', pos=7)
    outletCentral.findInLine(line, "areaAverage(outletCentralUnloaded) of magU", 'U')
    outletCentral.findInLine(line, "areaAverage(outletCentralUnloaded) of p")
    outletCentral.findInLine(line, "areaAverage(outletCentralUnloaded) of T")
    outletCentral.findInLine(line, "areaAverage(outletCentralUnloaded) of alphaRhoPhi")
    outletCentral.findInLine(line, "areaAverage(outletCentralUnloaded) of alphaPhi")
    outletCentral.findInLine(line, "faceZone outletCentralUnloaded TBulk")
    

    if (powers.any() and inlet.any() and outlet.any() and "Time =" == line[:6]):
        time = float(line.split()[2])
        break
    if (powers.any() and inlet.any() and outlet.any() and "ExecutionTime" == line[:13]):
        executionTime = float(line.split()[2])
        clockTime = float(line.split()[6])

# Compute the densities
inlet.computeDensity()
outlet.computeDensity()
inletFuel.computeDensity()
outletFuel.computeDensity()
inletCentral.computeDensity()
outletCentral.computeDensity()

# Reactor parameters
nFuelAssembly = 1500/6 # Number of fuel assembly in core
nominalPower = 937*MW/nFuelAssembly
fuelElementStructureFraction = 0.696911
SfuelElement = 0.00190587*m2 * (1-fuelElementStructureFraction)
SPlenum = 0.00222351*m2
gConst = 9.8 * m/(s**2)
Rconst = 8.314 * J/mol/K
gamma = 1.4
M_H2 = 2.014 * g/mol

outletNozzlePres = 0 * MPa
nozzleMinArea = pi * (8.720 * inch/m /2)**2
nozzleOutletArea = pi * (30.370 * inch/m /2)**2

gm1og = (gamma-1)/gamma
nozzleExhaustVelocity = sqrt(
    2/gm1og * Rconst*outlet['TBulk']/M_H2 * (
        1 - (outletNozzlePres/outlet['p'])**gm1og
    )
    # + outletFuelElemVelo[-1]**2
)

# thrustCore = nFuelAssembly * SfuelElement * outletFuelElemRho[-1] * outletFuelElemPres[-1]
# thrustCore = nFuelAssembly * outletMassFlowrate[-1] * outletFuelElemVelo[-1]
# IspCore = thrustCore/(nFuelAssembly * inletMassFlowrate[-1] * gConst)

IspNozzle = nozzleExhaustVelocity/gConst
thrustNozzle = nFuelAssembly * inlet['massFlow'] * nozzleExhaustVelocity

# effectiveExhaustVelocity = IspCore * gConst
effectiveExhaustVelocity = nozzleExhaustVelocity + outletNozzlePres * nozzleOutletArea/(nFuelAssembly*inlet['massFlow'])

# thrustCore2 = effectiveExhaustVelocity2*inletMFlow[-1]*nFuelAssembly

characteristicVelocity = outlet['p'] * nozzleMinArea / (nFuelAssembly*inlet['massFlow'])

print(f"""
Time           = {time} s
Execution Time = {executionTime} s
Clock Time     = {clockTime} s

Power:
    Expected          = {nominalPower      :8.6e} W
    Neutro integrated = {powers['neutro']  :8.6e} W ({powers['neutro']  /nominalPower*100:7.3f} %)
    Fluid  integrated = {powers['fluid']   :8.6e} W ({powers['fluid']   /nominalPower*100:7.3f} %)
    int Ah(Tf-Tco) dV = {powers['AhdeltaT']:8.6e} W ({powers['AhdeltaT']/nominalPower*100:7.3f} %)

Temperature:
    Expected         : 83.3 -> 1833.3 K
    Plenums          : {inlet['TBulk']       :.1f} -> {outlet['TBulk']       :6.1f} K
    Fuel element     : {inletFuel['TBulk']   :.1f} -> {outletFuel['TBulk']   :6.1f} K
    Central unloaded : {inletCentral['TBulk']:.1f} -> {outletCentral['TBulk']:6.1f} K

Velocity:
    Plenums          : {inlet['U']       :5.2f} -> {outlet['U']       :6.2f} m/s
    Fuel element     : {inletFuel['U']   :5.2f} -> {outletFuel['U']   :6.2f} m/s
    Central unloaded : {inletCentral['U']:5.2f} -> {outletCentral['U']:6.2f} m/s

Density:
    Plenums          : {inlet['rho']       :5.2f} -> {outlet['rho']       :6.2f} kg/m3
    Fuel element     : {inletFuel['rho']   :5.2f} -> {outletFuel['rho']   :6.2f} kg/m3
    Central unloaded : {inletCentral['rho']:5.2f} -> {outletCentral['rho']:6.2f} kg/m3

Pressure:
    Expected         : 4.1370 -> 3.4470 MPa
    Plenums          : {inlet['p']/MPa       :.4f} -> {outlet['p']/MPa       :.4f} MPa
    Fuel element     : {inletFuel['p']/MPa   :.4f} -> {outletFuel['p']/MPa   :.4f} MPa
    Central unloaded : {inletCentral['p']/MPa:.4f} -> {outletCentral['p']/MPa:.4f} MPa

Mass flowrate:
    Fuel element     : {inletFuel['massFlow']   :6.3f} -> {outletFuel['massFlow']   :6.3f} kg/s
    Central unloaded : {inletCentral['massFlow']:6.3f} -> {outletCentral['massFlow']:6.3f} kg/s
    Plenums          : {inlet['massFlow']       :6.3f} -> {outlet['massFlow']       :6.3f} kg/s
    Plenums (core)   : {nFuelAssembly*inlet['massFlow']:6.3f} -> {nFuelAssembly*outlet['massFlow']:6.3f} kg/s
    Expected         : 31.750 kg/s

From mass conservation:
    rho_i S_i u_i = rho_o S_o u_o => {
        inlet['rho']  * inlet['S']  * inlet['U']:.3f} = {
        outlet['rho'] * outlet['S'] * outlet['U']:.3f} kg/s

From Perfect Gas law:
    r = P_i/(rho_i T_i) = P_o/(rho_o T_o) => {
        inlet['p'] /(inlet['rho'] *inlet['TBulk']) :.3f} = {
        outlet['p']/(outlet['rho']*outlet['TBulk']):.3f} m²/s²/K

Rocketery (to-do: need to further understand, not sure on value):
    From Rocket propulsion elements 8th edition
    Exhaust velocity: {nozzleExhaustVelocity/(m/s):.2f} m/s (v_2)
    c               : {effectiveExhaustVelocity/(m/s):.2f} m/s (Effective exhaust velocity (Isp*g))
    c*              : {characteristicVelocity/(m/s):.2f} m/s (Characteristic velocity)
    thrust          : {thrustNozzle/N:.1f} N (m*c)
    Isp             : {IspNozzle/s:.1f} s (v2/g)
""")


# c      : {effectiveExhaustVelocity2/(m/s):.2f} m/s (Effective exhaust velocity (v2 + p2*A2/m))
# Fuel element
#             thrust : {:.1f} N (core (S*rho*P1)thrustCore/N))
#             thrust : {:.1f} N (core (m*c)nFuelAssembly*inletMassFlowrate[-1]*effectiveExhaustVelocity /N))
#             Isp    : {IspCore/s:.1f} s (core (F/(m*g))
#             c      : {effectiveExhaustVelocity/(m/s):.2f} m/s (Effective exhaust velocity (Isp*g))
#         Nozzle
#             p1/p2  : {:.3f}".format(outletFuelElemPres[-1]/outletNozzlePres))