# 2D Void Motion No Phase Change

Very simple tutorial displaying a two-phase case without mass transfer between
phases (obtained by not inserting the phaseChangeModel subdictionary in the
phaseProperties dictionary). Please refer to the 1D_boiling tutorial for
details on how to use the two-phase flow solver. The main additional feature
employed in this tutorial compared to 1D_boiling is the use of the
initialAlphas subdict in the vapourProperties subdict of the phaseProperties
dictionary. It is used to provide potentially different initital phase
fractions for different cellZones (as an alternative to the use of setFields).
