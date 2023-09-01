"""
Author: Thomas Guilbaud, 2022/10/21

Hydrogen thermophysical properties.
"""

#==============================================================================*
#                                   Imports
#==============================================================================*

import numpy as np
import matplotlib.pyplot as plt
import os
import sys

rootPath = os.path.abspath(__file__)
for i in range(1, 4):
    sys.path.append("/".join(rootPath.split("/")[:-i]))

from unitsOpenMC import *
from physicsFunctions import gasDensity


#==============================================================================*

interpolation = lambda x, ya, yb, xa, xb: (yb-ya)/(xb-xa) * (x-xa) + ya

def specificHeatCapacity(T: float, A: float=2.016, isFromFit: bool=False) -> float:
    # temperatureList = [175.0, 200.0, 225.0, 250.0, 275.0, 300.0, 325.0, 350.0, 375.0, 400.0, 450.0, 500.0, 550.0, 600.0, 650.0, 700.0, 750.0, 800.0, 850.0, 900.0, 950.0, 1000.0, 1050.0, 1100.0, 1150.0, 1200.0, 1250.0, 1300.0, 1350.0, 1400.0, 1500.0, 1600.0, 1700.0, 1800.0, 1900.0, 2000.0, 2100.0, 2200.0, 2300.0, 2400.0, 2500.0, 2600.0, 2700.0, 2800.0, 2900.0, 3000.0, 3500.0, 4000.0, 4500.0, 5000.0, 5500.0, 6000.0]
    # specificHeatCapacityList = [13.12, 13.53, 13.83, 14.05, 14.2, 14.31, 14.38, 14.43, 14.46, 14.48, 14.5, 14.51, 14.53, 14.55, 14.57, 14.6, 14.65, 14.71, 14.77, 14.83, 14.9, 14.98, 15.06, 15.15, 15.25, 15.34, 15.44, 15.54, 15.65, 15.77, 16.02, 16.23, 16.44, 16.64, 16.83, 17.01, 17.18, 17.35, 17.5, 17.65, 17.8, 17.93, 18.06, 18.17, 18.28, 18.39, 18.91, 19.39, 19.83, 20.23, 20.61, 20.96]
    #
    # for i in range(len(temperatureList)-1):
    #     Ta, Tb = temperatureList[i], temperatureList[i+1]
    #     if (Ta <= T and T <= Tb):
    #         cpa, cpb = specificHeatCapacityList[i], specificHeatCapacityList[i+1]
    #         cp = interpolation(T, cpa, cpb, Ta, Tb)
    #         return(cp * 1e3) # J/kg/K

    if (isFromFit):
        z = [1.39051574e-19, -1.77480852e-15, 9.28566400e-12, -2.54178633e-08, 3.81902598e-05, -2.98179826e-02, 1.17094171e+01, 1.26644563e+04]
        p = np.poly1d(z)
        return(p(T))

    cp = 0
    t = T/1000
    if (298 <= T and T <= 1000):
        cp = 33.066178 -11.363417*t + 11.432816*t**2 -2.772874*t**3 -0.158558/(t**2)
    elif (1000 <= T and T <= 2500):
        cp = 18.563083 +12.257357*t -2.859786*t**2 +0.268238*t**3 +1.977990/(t**2)
    elif (2500 <= T and T <= 6000):
        cp = 43.413560 -4.293079*t +1.272428*t**2 -0.096876*t**3 -20.533862/(t**2)
    return(cp/(A*1e-3))

def thermalConductivity(T: float) -> float:
    temperatureList = [25, 100, 175, 250, 325, 400]
    thermalConductivityList = [20.8, 68.3, 117, 161, 198, 234]

    for i in range(len(temperatureList)-1):
        Ta, Tb = temperatureList[i], temperatureList[i+1]
        if (Ta <= T and T <= Tb):
            ka, kb = thermalConductivityList[i], thermalConductivityList[i+1]
            k = interpolation(T, ka, kb, Ta, Tb)
            return(k * 1e-3) # W/m/K
    return(0)


#==============================================================================*

if __name__ == '__main__':
    temperatureRange = np.linspace(300, 3000, 100)
    cpList = [specificHeatCapacity(T) for T in temperatureRange]

    z = np.polyfit(temperatureRange, cpList, 7)
    print("Cp :", z)
    p = np.poly1d(z)

    plt.figure()
    plt.plot(temperatureRange, cpList, label="Table")
    plt.plot(temperatureRange, p(temperatureRange), label="Fit")
    plt.xlabel("Temperature [K]")
    plt.ylabel("Specific heat capacity [J/kg/K]")
    plt.grid(True)
    plt.legend()


    temperatureRange = np.linspace(25, 400, 20)
    kList = [thermalConductivity(T) for T in temperatureRange]

    z = np.polyfit(temperatureRange, kList, 1)
    print("k :", z)
    p = np.poly1d(z)

    plt.figure()
    plt.plot(temperatureRange, kList, label="Table")
    plt.plot(temperatureRange, p(temperatureRange), label="Fit")
    plt.xlabel("Temperature [K]")
    plt.ylabel("Thermal conductivity [W/m/K]")
    plt.grid(True)
    plt.legend()


    temperatureRange = np.linspace(150, 3000, 100)
    rhoList = [gasDensity(1e5, T)/(kg/m3) for T in temperatureRange]

    z = np.polyfit(temperatureRange, rhoList, 7)
    print("rho :", z)
    p = np.poly1d(z)

    plt.figure()
    plt.plot(temperatureRange, rhoList, label="Table")
    plt.plot(temperatureRange, p(temperatureRange), label="Fit")
    plt.xlabel("Temperature [K]")
    plt.ylabel(r"Density [kg/m$^3$]")
    plt.grid(True)
    plt.legend()

    plt.show()

#==============================================================================*
