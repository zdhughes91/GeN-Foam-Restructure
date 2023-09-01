# Serpent to Foam XS

Author: Thomas Guilbaud, EPFL/Transmutex SA  
Last update: 13/09/2022  

---

This application aims to post-process the output file from the Monte Carlo code
Serpent2 to OpenFOAM/GeN-Foam. It has been greatly inspired by the Octave
script.

## Build

Execute the following command:
``` bash
./script/buildApp.sh
```
The executable file is generated in `build/bin`.

## How to use it

### Graphical User Interface (GUI)

Run:
``` bash
python SerpentToFoamXS.py
# or
./SerpentToFoamXS
```
Or if the user has already an input file:
``` bash
python SerpentToFoamXS.py path/to/userInputFile
# or
./SerpentToFoamXS path/to/userInputFile
# or
./SerpentToFoamXS path/to/userInputFile path/to/serpentFile_res.m
```

`path/to/serpentFile_res.m` can be relative. But the user output file of Serpent
(`serpentFile_res.m`) must be given using an absolute path only in the
`userInputFile`.

The user can edit the `userInputFile` using the GUI view and then click `Save`
to keep a copy of the input. A copy is always made when the user press
`Extract` (The name of the saved user input file will never be changed, which
will erase the previous user input file).

### Bash

Simply execute the previous command using `-b` or `--bash`:
``` bash
python SerpentToFoamXS.py --bash path/to/userInputFile
# or
./SerpentToFoamXS --bash path/to/userInputFile
# or
./SerpentToFoamXS --bash path/to/userInputFile path/to/serpentFile_res.m
```

To extract all the `nuclearData` files, add `-a` or `--all`:
``` bash
python SerpentToFoamXS.py --bash --all path/to/userInputFile
# or
./SerpentToFoamXS --bash --all path/to/userInputFile
# or
./SerpentToFoamXS --bash --all path/to/userInputFile path/to/serpentFile_res.m
```
This parameter only works in bash mode.

--------------------------------------------------------------------------------

## To check
- [x] Which Beta and Lambda values when Point-Kinetics is checked
- [ ] `INF_FLX=0` issue (see "Integral Fluxes" in `extraction.py`)
- [x] Why `radialOrientation` and `axialOrientation` are only in
    `nuclearDataRadialExp` and not in `nuclearDataAxialExp` ?

## To do
See in `nuclearData` from `2D_onePhaseAndPointKineticsCoupling`
- Add the other parameters for **point-kinetics** which are not in the Octave file
    - [x] promptGenerationTime              float; (user or extracted ?) (use ADJ_PERT_GEN_TIME ?)
    - [x] feedbackCoeffFastDoppler          float; (user)
    - [x] feedbackCoeffTFuel                float; (user)
    - [x] feedbackCoeffTClad                float; (user)
    - [x] feedbackCoeffTCool                float; (user)
    - [x] feedbackCoeffRhoCool              float; (user)
    - [x] feedbackCoeffTStruct              float; (user)
    - [ ] absoluteDrivelineExpansionCoeff   float; (extract or user ?)
    - [ ] controlRodReactivityMap           table; (extract or user ?)
    - [ ] fuelFeedbackZones                 table; (user)
    - [ ] coolFeedbackZones                 table; (user)
    - [ ] structFeedbackZones               table; (user)
    - [ ] fuelFeedbackZones                 table; (user)
- Add the other parameters for **diffusion** which are not in the Octave file
    - [x] adjustDiscFactors                 bool; (user)
    - [x] useGivenDiscFactors               bool; (user)
- Add the other parameters for **both** which are not in the Octave file
    - [x] fastNeutrons                      bool; (user)

- Add meaning of Serpent/Foam cell conversion (how to use it)
