rad = "minicore_3D"; %name of serpent input file
coreState = "N"; % N nominal; A axially expanded; R radially expanded; T different fuel temp; C different coolant density; CT different coolant temp; CL different clad temp
zeroeff = "zero";%zero or eff
pTarget = 1e6;%core power
keff = 1;%initial guess

% do not touch from here
if (strcmp("R",coreState))
	expansionFromNominalR = 1;
	radialOrientationX = 1;
	radialOrientationY = 1;
	radialOrientationZ = 0;
	AxialOrientationX = 0 ;
	AxialOrientationY = 0 ;
	AxialOrientationZ = 1 ;
end

if (strcmp("A",coreState))
	expansionFromNominalA = 1;
end

if (strcmp("T",coreState))
	TfuelRef = 560;
	TfuelPerturbed = 860;
end
if (strcmp("C",coreState))
	rhoCoolRef = 752.06;
	rhoCoolPerturbed = 552.06;
end
if (strcmp("CT",coreState))
	TCoolRef = 560;
	TCoolPerturbed = 633;
end
if (strcmp("CL",coreState))
	Tcladref = 560;
	TcladPerturbed = 560;
end
% do not touch before here

% Serpent name of the univereses you want to include
SERPENT_NAME(1, [1:  2])  = '11' ;
SERPENT_NAME(2, [1:  3])  = '100' ;
SERPENT_NAME(3, [1:  3])  = '200' ;
SERPENT_NAME(4, [1:  3])  = '300' ;
SERPENT_NAME(5, [1:  3])  = '150' ;
SERPENT_NAME(6, [1:  3])  = '250' ;
SERPENT_NAME(7, [1:  3])  = '350' ;

% corresponding cellZones names in the neutronics mesh
OF_NAME = [
'radialReflector';
'UO2bundle';
'MOXbundleCenter';
'MOXbundleCorner';
'UO2bundleCR';
'MOXbundleCRCenter';
'MOXbundleCRCorner';
];


% corresponding volumetric fuel fractions //useful only if there is fuel. Otherwise not used (but still necessary to put a number)
fuelFraction = [
1;
1;
1;
1;
1;
1;
1;
];

