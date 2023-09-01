# Tips and tricks {#TIPS}

* When unsure about the functionalities of a class, take a look at the corresponding header file in the source code. In most case you'll find useful information.

* When unsure about the meaning of an input parameters, take a look at the tutorials. Many of them are commented. Another option is the famous "banana method". Type in anything (for instance "banana"). In many cases, the solver will tell you that "banana" is not a valid option, and it will suggest the valid options.

* defaultPrec has 1/m3 units except for the adjoint solver that needs 1/m2/s.

* In point kinetics, pTarget (in reactorState) MUST be the same as the one used for reaching the steady-state. As power, GeN-Foam uses what it finds under powerDensity, or under the powerDensity of the fluidRegion if it does not find a powerDensity in the neutroRegion. pTarget does not enter the calculation, it is used simply to plot the results 

