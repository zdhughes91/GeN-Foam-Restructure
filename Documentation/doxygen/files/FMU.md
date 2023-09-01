# FMU coupling {#FMU}

GeN-Foam provides several interface points to communicate with [Functional Mock-up Units](https://fmi-standard.org/) (FMUs). FMUs are containers of software and data that are based on a widely employed communication standard called Functional Mockup Interface (FMI). The FMI is developed by an industrial consortium led by the Modelica Association.


## Compiling

To use the FMI coupling interface in GeN-Foam, the user have to install the [FMU4FOAM](https://github.com/DLR-RY/FMU4FOAM) project developed by the DLR using the following commands:
```
cd GeN-Foam
git clone https://github.com/DLR-RY/FMU4FOAM.git
cd FMU4FOAM/
./build-ECI4FOAM.sh && ./Allwmake
```

Then the GeN-Foam project can be build as usual.


## Potential Issues

Make sure that: 
- `-I./../../../FMU4FOAM/ECI4FOAM/src/externalComm/lnInclude` has been added to `EXE_INC` in [GeN-Foam/classes/thermalHydraulics/src/Make/options](GeN-Foam/classes/thermalHydraulics/src/Make/options).
- `-lexternalComm` has been added to `LIB_LIBS` in [GeN-Foam/classes/thermalHydraulics/src/Make/options](GeN-Foam/classes/thermalHydraulics/src/Make/options).