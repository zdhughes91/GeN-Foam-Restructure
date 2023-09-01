### IMPORTS ###

from dependencies.readDict import readDict
from dependencies.runCase  import runCase
import sys

### MAIN ###

i = 0
solverName, regions, parallelOptions, parameterTable = readDict()
for parameterList in parameterTable :
	i = runCase(solverName, regions, parallelOptions, parameterList, i)