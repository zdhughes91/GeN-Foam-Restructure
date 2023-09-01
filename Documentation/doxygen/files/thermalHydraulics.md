# Thermal-hydraulics {#TH}

## Introduction

Both single- and two-phase simulations can be performed using GeN-Foam. All sub-solvers were developed for a coarse-mesh porous-medium treatment of complex structures such as core and heat exchanger, and for a standard RANS treatment of clear-fluid regions. The sub-solvers automatically switch from a porous-medium (coarse-mesh) treatment to a standard CFD (fine-mesh) treatment when the colume fraction of the sub-scale structures is set to zero. This allows for an implict coupling of porous-medium (sub-channel-like in 2D and 3D, or system-code-like) treatment of compelex structures (e.g., core and heat exchnagers) with a standard CFD treatment of clear-fluid regions (e.g., plena and pools).

A coarse-mesh porous-medium treatment of the core implies that the core is modeled without resolving the sub-scale structure (e.g., the fuel rods or the heat exchanger tubes). As a matter of fact, in principle and for consistency, the finest radial mesh chosen by a user should not finer than one cell per pin cell. A porous-medium formulation derives from a volume averaging of the Navier-Stokes equations. The volume averaging results in source terms that describe the interaction (drag and heat transfer) of the fluid with the sub-scale structure. In GeN-Foam, these source terms are modeled using user-selectable correlations for drag (e.g., correlations for the Darcy friction factor) and heat transfer (e.g., correlations for the Nusselt number). In this sense, a porous-medium model can be associated with a 3-D version of a system code.

With regards to the modelling of the sub-scale structures, GeN-Foam allows to model simultaneously in the same region both a "power model" and a "passive structure". Power models are used to model fo instance the nuclear fuel (based on a 1-D approximation), electrically heated rods, or a fixed temperature body (which can be used to approximate a heat exchanger). Passive structures are structures that passively heats up or cool down based on their own heat capacity, volumetric area, and heat tranfer with the coolant. This can be used to model structures like the assembly wrappers or the reflectors.

All thermal-hydraulics functionalities are handled by the class *thermalHydraulicsModel.H*, the derived classes for the various sub-solvers (see below), and a thermal-hydraulic library that can be found under */GeN-Foam/classes/thermalHydraulics/src*.

<div class="border-box" style='padding:0.1em; margin-left: 4em;  margin-right: 8em;  border: 1px solid gray; background-color:#f2f3fa; color:#05134a'>
<b>The porous-medium approach in GeN-Foam</b>

GeN-Foam was born for safety analyses and, to reduce computational footprint, its base approach is to model for instance the core as a porous medium. In a porous-medium approach, the fuel and other structures (e.g., the assembly wrappers) are modelled using sub-scale models. This means that, in each cell, we have the fluid and one or two lumped models for the sub-scale structures. The simplest structures in GeN-Foam are the passive structures. These passive structures are simply modelled as a heat capacity and they can be used for instance to model assembly wrappers, or reflector structures. In essence, the fluid will interact cell-by-cell with these passive structures: they will take energy from the fluid if the temperature of the fluid is higher than that of the surface of the structure, and viceversa. A more complex example of sub-scale structure is given by the powerModels. An example of a power model is the nuclearFuelPin, which can be used to model a standard pin-type fuel. This power model is capable of getting the power density from neutronics, solving a 1-D model for heat transfer in the fuel, and give back to the fluid the temperature at the surface of the cladding. The fluid will then be capable of calculating the heat transfer with the fuel based on the cladding surface temperature and the Nusselt number. Each cellZone can host one passive structure and one powerModel.
<br>
</div>

## Sub-solvers

Thermal-hydraulics calculations are performed by classes derived from *thermalHydraulicsModel.H* that contain specific sub-solvers:
* *onePhase* for single-phase calculations, using the formulation proposed in Refs. \cite Radman2019ADesign \cite RADMAN2021111178 \cite RADMAN2021111422  (see *onePhase.H*)
* *onePhaseLegacy* for single-phase calculations, using the formulation proposed in Ref. \cite FIORINA201524 (see *onePhaseLegacy.H*)
* *twoPhase* for adjoint diffusion calculations, using the formulation proposed in Refs. \cite Radman2019ADesign \cite RADMAN2021111178 \cite RADMAN2021111422 (see *twoPhase.H*)
For the user, the derived classes translate into runtime selectable models. The specific sub-solver  to be used in a simulation can be selected at runtime in the *phaseProperties* dictionary in *constant/fluidRegion/*, using the keyword *thermalHydraulicsType*. 


## Porous-medium properties

The various parameters to be used in a porous-medium simulation can be set using the *phaseProperties* dictionary. 
<div class="border-box" style='padding:0.1em; margin-left: 4em;  margin-right: 8em;  border: 1px solid gray; background-color:#f2f3fa; color:#05134a'>
<b>The *phaseProperties* dictionary</b>

The *phaseProperties* dictionary can be found in *constant/fluidRegion/*. It is a large dictionary that can be used to: chose the sub-solver to be used (one-phase, legacy one-phase or two-phase); set various properties of the phases (beside basic thermo-physical properties defined in the *thermophysicalProperties* dictionary); set the properties of the sub-scale structures (fuel pins, heat exchangers, etc) in the porous zones, including the possibility to assign a *powerModel* for power production (e.g., nuclear fuel, or constant power) and the *passiveProperties* of another sub-structure that interacts thermally with the fluid (for instance the wrappers in sodium fast reactors). In addition models are available to model pumps and heat exchangers.  The name of the porous zones must coincide with that of the cellZones of the fluidRegion mesh.  
<br>
<br>
One can find detailed, commented examples in the tutorials 
[3D_SmallESFR](https://gitlab.com/foam-for-nuclear/GeN-Foam/-/tree/master/Tutorials/3D_SmallESFR/rootCase/constant/fluidRegion/phaseProperties) (single phase) and
[1D_boiling](https://gitlab.com/foam-for-nuclear/GeN-Foam/-/tree/master/Tutorials/1D_boiling/constant/fluidRegion/phaseProperties) (two phases). In addition, an example on how to use a two-dimensional flow-regime map can be found in [1D_PSBT_SC](https://gitlab.com/foam-for-nuclear/GeN-Foam/-/tree/master/Tutorials/1D_PSBT_SC/constant/fluidRegion/phaseProperties).
<br>
<p>
**Drag models**. Currently available models to describe presure drops induced by the sub-scale structure or by a second phase include:
<UL>
<LI> Fluid-fluid drag models (see *FFDragCoefficientModel.H*)	
	<UL>
	<LI> Autruffe (see *AutruffeFFDragCoefficient.H*)
	<LI> Bestion (see *BestionFFDragCoefficient.H*)
	<LI> Bestion as in TRACE (see *BestionTRACEFFDragCoefficient.H*)
	<LI> No Kazimi (see *NoKazimiFFDragCoefficient.H*)
	<LI> Schiller Naumann (see *SchillerNaumannFFDragCoefficient.H*)
	<LI> Wallis (see *WallisFFDragCoefficient.H*)
	</UL>	
<LI> Fluid-structure drag models (see *FSDragCoefficientModel.H*)
	<UL>
	<LI> Baxi Dalle Donne (see *BaxiDalleDonneFSDragCoefficient.H*)
	<LI> Churchill (see *ChurchillFSDragCoefficient.H*)
	<LI> Engel as in TRACE (see *EngelFSDragCoefficient.H*)
	<LI> Modified Engel (see *modifiedEngelFSDragCoefficient.H*)
	<LI> No Kazimi (see *NoKazimiFSDragCoefficient.H*)
	<LI> Rehme (see *RehmeFSDragCoefficient.H*)
	<LI> Drag coefficient as a A*Re^B+C (see *ReynoldsPowerFSDragCoefficient.H*)
	<LI> Drag coefficient as a (A*log10(Re)+B)^C (see *ColebrookFSDragCoefficient.H*)
	</UL>	
<LI> Two-phase drag multipliers	(see *twoPhaseDragMultiplierModel.H*)
	<UL>
	<LI> Chen Kalish (see *ChenKalishTwoPhaseDragMultiplier.H*)
	<LI> Constant (see *constantTwoPhaseDragMultiplier.H*)
	<LI> Kaiser 74 (see *Kaiser74TwoPhaseDragMultiplier.H*)
	<LI> Kaiser 88 (see *Kaiser88TwoPhaseDragMultiplier.H*)
	<LI> Kottowski Savatteri (see *KottowskiSavatteriTwoPhaseDragMultiplier.H*)
	<LI> Lockhart Martinelli (see *LockhartMartinelli.H*)
	<LI> Lottes Flinn(see *LottesFlinnTwoPhaseDragMultiplier.H*)
	<LI> Lottes Flinn Nguyen(see *LottesFlinnNguyenTwoPhaseDragMultiplier.H*)
	</UL>	
</UL>
<br>
<p>
**Heat transfer models**. Currently available models to describe the enrgy exchange with a sub-scale structure or with a second phase include:
<UL>
<LI> Fluid-fluid heat-tranfer models (see *FFHeatTransferCoefficientModel.H*)	
	<UL>
	<LI> No Kazimi (see *NoKazimiFFHeatTransferCoefficient.H*)
	<LI> Nusselt number correlation in the form Nu = A_+B_*(Re^C_)*(Pr^D_) (see *NusseltFFHeatTransferCoefficient.H*)
	</UL>	
<LI> Fluid-structure heat-tranfer models (see FSHeatTransferCoefficientModel.H)
	<UL>
	<LI>  Nusselt number correlation in the form Nu = A_+B_*(Re^C_)*(Pr^D_) (see *NusseltFSHeatTransferCoefficient.H*)
	<LI>  Nusselt number correlation in the form Nu = A_+B_*(Re^C_)*(Pr^D_), plus and additional heat transfer coefficient to take into account the resistance of a wall, such that  H = Nu * kappa / Dh + H_wall (see *NusseltAndWallFSHeatTransferCoefficient.H*)
	<LI>  Shah (see *ShahFSHeatTransferCoefficient.H*)
	<LI>  Gorenflo (see *GorenfloFSHeatTransferCoefficient.H*)
	<LI>  multiRegimeBoilingTRACE - multi-regime heat transfer coefficient that replicates what TRACE does below CHF (see *multiRegimeBoilingTRACEFSHeatTransferCoefficient.H*) ADD REF GAUTHIER
	<LI>  multiRegimeBoilingTRACE - multi-regime heat transfer coefficient that replicates what TRACE does, including CHF and post-CHF. Not verified! It requires specifiyng the *multiRegimeBoilingTRACECHF* model for the water.structure heat tranfer, and the *multiRegimeBoilingVapourTRACE* model for the vapour.structure heat tranfer (see *multiRegimeBoilingTRACECHFFSHeatTransferCoefficient.H* and *multiRegimeBoilingVapourTRACEFSHeatTransferCoefficient.H*). Please notice that a lookup table for CHF is still missing. 
	<LI>  multiRegimeBoiling - multi-regime heat transfer coefficient that replicates  the same logic as TRACE, but with more flexibility for user-selectable sub-models (see *multiRegimeBoiling.H*). Can be used below CHF.
	<LI>  Sub-models employed by the multi-regime models:
		<UL>
		<LI> Critical heat flux related models (see *CHFModel.H*):
			<UL>
			<LI> Critical heat flux models
				<UL> 
				<LI> Constant, user-selctable value (see *constantCHF.H*)
				<LI> Lookup table, not yet implemented (empty class at *lookUpTableCHF.H*)
				</UL>
			<LI> Leidenfrost models (see *TLFModel.H*)
				<UL> 
				<LI> Groeneveld Stewart (see *GroeneveldStewartTLF.H*)
				</UL>
			</UL>
		<LI> Flow Enhancement Factor Models (see *flowEnhancementFactorModel.H*)
			<UL>
			<LI> Chen (see *ChenFlowEnhancementFactor.H*)
			<LI> COBRA-TF (see *COBRA-TFFlowEnhancementFactor.H*)
			<LI> Rezkallah Sims (see *RezkallahSimsFlowEnhancementFactor.H*)
			</UL>
		<LI> Post-CHF models
			<UL>
			<LI> Cachard (for liquid) (see *CachardLiquidFSHeatTransferCoefficient.H*)
			<LI> Cachard (for vapour) (see *CachardVapourFSHeatTransferCoefficient.H*)
			</UL>
		<LI> Sub-Cooled Boiling Fraction Models (see *subCooledBoilingFractionModel.H*)
			<UL>
			<LI> Constant (see *constantSubCooledBoilingFraction.H*)
			<LI> Saha Zuber (see *SahaZuberSubCooledBoilingFraction.H*)
			</UL>
		<LI> Superposition Nucleate Boiling (see *superpositionNucleateBoilingFSHeatTransferCoefficient.H*)
		<LI> Suppression factor models (see *suppressionFactorModel.H*)
			<UL>
			<LI> Chen (see *ChenSuppressionFactor.H*)
			<LI> COBRA-TF (see *COBRA-TFSuppressionFactor.H*)
			</UL>
		<LI> Temperature of the onset of nucleate boiling (see *TONBModel.H*)
			<UL>
			<LI> Basu (see *BasuTONB.H*)
			</UL>
		</UL>
	</UL>	
</UL>
<br>
<p>
**Special models**. Dedicated models for specific sub-scale structures inlcude:
<UL>
<LI> Power models, i.e., active media that can provide and substract energy, including:
	<UL>
	<LI> Fixed (possibly time-dependent) power (see *fixedPower.H* and the tutorial *1D_CHF/imposedPower*)
	<LI> Fixed (possibly time-dependent) temperature (see *fixedTemperature.H* and the tutorial *1D_CHF/imposedTemperature*)
	<LI> Heated pin, typically used for electrically heated pins (see *heatedPin.H* and the tutorial *2D_KNS37-L22*)
	<LI> Nuclear fuel pin (see *nuclearFuelPin.H* and the tutorials *3D_SmallESFR* and *2D_FFTF*)
	<LI> Lumped-parameter nuclear structure (see *lumpedNuclearStructure.H*) and the tutorial *1D_thermalMSR_pointKinetics*
	</UL>
<LI> A heat exchanger model that is used to model the heat transfer between two disconnected regions, for instance representing the primary and secondary circuit (see *heatExchanger.H* and the tutorials *1D_HX* and *2D_FFTF*)
<LI> A pump model used to set a (possibly time-dependent) momentum source (see *pump.H* and tutorials *2D_FFTF* and *2D_MSFR*).
</UL>
<br>
<p>
**Models for two-phase flows**. Currently available models for two-phase flow simulations include:
<UL>
<LI> Contact partition models (see *contactPartitionModel.H*)
	<UL>
	<LI> Linear (see *linearContactPartition.H*)
	<LI> Complementary (see *complementaryContactPartition.H*)
	</UL>
<LI> Disperions models (see *dispersionModel.H*)
	<UL>
	<LI> Constant (see *constantDispersion.H*)
	</UL>	
<LI> Fluid diameter models (see *fluidDiameterModel.H*)
	<UL>
	<LI> Iso-molar bubble (see *isomolarBubbleFluidDiameter.H*)
	<LI> Iso-thermal bubble  (see *isothermalBubbleFluidDiameter.H*)
	<LI> Pipe film (see *pipeFilmFluidDiameter.H*)
	<LI> Wallis film (see *WallisFilmFluidDiameter.H*)
	</UL>
<LI> Virtual mass models (see *virtualMass.H*)
	<UL>
	<LI> Virtual mass coefficient (see *virtualMassCoefficientModel.H*)
	</UL>	
<LI> Interfacial area models (see *interfacialAreaModel.H*)
	<UL>
	<LI> Annular (see *annularInterfacialArea.H*)
	<LI> No Kazimi  (see *NoKazimiInterfacialArea.H*)
	<LI> Schor (see *SchorInterfacialArea.H*)
	<LI> Spherical (see *sphericalInterfacialArea.H*)
	</UL>	
<LI> Phase change models
	<UL>
	<LI> Forced constant (see *forcedConstantPhaseChange.H*)
	<LI> Heat driven  (see *heatDrivenPhaseChange.H*)
	<LI> Latent heat (see *latentHeatModel.H*)
		<UL>
		<LI> Fink Leibowitz for sodium (see *FinkLeibowitzLatentHeat.H*)
		<LI> NIST interpolation for water (see *waterLatentHeat.H*)
		<LI> Use value for thermoPhysicalProerties dictionary (see *fromThermophysicalPropertiesLatentHeat.H*)
		</UL>
	<LI> Saturation temperature/pressure (see *saturationModel.H*)
		<UL>
		<LI> Browning Potter for sodium (see *BrowningPotterSaturation.H*)
		<LI> NIST interpolation for water (see *waterSaturation.H*)
		<LI> TRACE model interpolation for water - YET TO BE VEIFIED (see *waterTRACESaturation.H*)
		<LI> Constant temperature (see *constantTemperatureSaturation.H*)
		</UL>
	</UL>	
</UL>
<br>
<p> 
**Regime maps**. In GeN-Foam the possibility exists to employ 1- and 2-dimensinal regime maps in order to employ different models for different flow conditions. This can be used for instance in one-phae simulation to provide different correlations for turbulent and laminar flow  (see [3D_SmallESFR](https://gitlab.com/foam-for-nuclear/GeN-Foam/-/blob/master/Tutorials/3D_SmallESFR/rootCase/constant/fluidRegion/phaseProperties) for a commented example), or in 2-phase flow simulations to provide full regime maps (see [1D_boiling](https://gitlab.com/foam-for-nuclear/GeN-Foam/-/blob/master/Tutorials/1D_boiling/constant/fluidRegion/phaseProperties) for a commented example of 1-dimensinal map, [1D_PSBT_SC](https://gitlab.com/foam-for-nuclear/GeN-Foam/-/blob/master/Tutorials/1D_PSBT_SC/constant/fluidRegion/phaseProperties) and for a commented example of 2-dimensinal map). Multiple maps can be used in the same simulation. In addition, regime-maps models can be mixed with multi-regime models, as in [1D_PSBT_SC](https://gitlab.com/foam-for-nuclear/GeN-Foam/-/blob/master/Tutorials/1D_PSBT_SC/constant/fluidRegion/phaseProperties), where a *preCHFTraceRegimeMap* is employed to assign models for phase dispersion, interfacial area and bubble diameter, while a single multi-regime model is employed to describe heat transfer between liquid and structure throughout the various regimes.
<p>
N.B.: Anisotropic pressure drops can be set by using the keywords *transverseDragModel* (Blasius, GunterShaw, same) and *principalAxis*(localX, localY, localZ) in the sub-dictionary *dragModels.(nameOfPhase).structure.(nameOfCellZones)*. *principalAxis* sets the axis on which the nominal dragModel is used. *transverseDragModel* sets the model to be used on the two directions that are perpendicular to *principalAxis*. If *same* is chosen as *transverseDragModel*, the code will use the nominal model in all directions, but with the possibility of an anisotropic hydraulic diameter. The anisotropy of the hydraulic diameter can be set using the keyword *localDhAnisotrpy* and assign to it a vector of 3 scaling factors (one for each local directions).
<p>
N.B.2: The thermal hydraulic class can make use of a local coordinate system, which can be used by setting the keywords *localX* and *localY*  in the sub-dictionary *dragModels.(nameOfPhase).structure.(nameOfCellZones)* of the dictionary *constant/fluidRegion/phaseProperties*. A local coordinate system can be used for instance when one knows the pressure drop correlation in a direction that is different from the x, y, and z directions of the global coordinate system. Besides drag models, the local coordinate system can be used also for defining a tortuosity (keyword *localTortuosity*, to be defined as a vector in the local coordinate system).
</div>
<br>



## Physical properties

<div class="border-box" style='padding:0.1em; margin-left: 4em;  margin-right: 8em;  border: 1px solid gray; background-color:#f2f3fa; color:#05134a'>
<b>The *g* dictionary</b>

The *g* dictionary can be found under *constant/fluidRegion/*. It is a standard OpenFOAM dictionary that allows specifying the gravitational acceleration.
</div>
<br>

<div class="border-box" style='padding:0.1em; margin-left: 4em;  margin-right: 8em;  border: 1px solid gray; background-color:#f2f3fa; color:#05134a'>
<b>The *thermophysicalProperties* dictionary</b>

The *thermophysicalProperties* dictionary can be found under *constant/fluidRegion/*. It is a standard OpenFOAM dictionary that allows defining the thermo-physical properties of the coolant. When performing two-phase flow analyses, two dictionaries must be employed named*thermophysicalProperties.(name of fluid)*. The name of the two fluids are defined in the *phaseProperties* dictionary.
<br><br>One can find a detailed, commented example in hte tutorials 
[3D_SmallESFR](https://gitlab.com/foam-for-nuclear/GeN-Foam/-/tree/master/Tutorials/3D_SmallESFR/rootCase/constant/fluidRegion/thermophysicalProperties) (one-phase), 
[1D_boiling (liquid)](https://gitlab.com/foam-for-nuclear/GeN-Foam/-/tree/master/Tutorials/1D_boiling/constant/fluidRegion/thermophysicalProperties.liquid) (two-phase, liquid)
[1D_boiling (vapour)](https://gitlab.com/foam-for-nuclear/GeN-Foam/-/tree/master/Tutorials/1D_boiling/constant/fluidRegion/thermophysicalProperties.vapour) (two-phase, vapour)
</div>
<br>


## Turbulence properties

<div class="border-box" style='padding:0.1em; margin-left: 4em;  margin-right: 8em;  border: 1px solid gray; background-color:#f2f3fa; color:#05134a'>
<b>The *turbulenceProperties* dictionary</b>

The *turbulenceProperties* dictionary can be found under *constant/fluidRegion/*. It is a standard OpenFOAM dictionary that allows defining the turbulence model to be used. 
<br><br>When clear-fluid simulations (i.e., without porous zones) are performed, on can used the standard kEpsilon model of OpenFOAM.
<br><br>When porous zones are present in the simulation, it is reccomended to use *porousKEpsilon* (see *porousKEpsilon.H*). The only difference w.r.t. the standard k-epsilon model is that it forces k and epsilon to equilibrium values inside the porous zones. These equilibium values can be set in the *porousKepsionProperties* sub-dictionary. Please notice that a porous medium simulation using the equilibium values of k and epsilon for the sub-scale structure (viz., the values inside a fuel sub-channel) would entail the risk of an unstable solution. This is due to the fact that the turbulent viscosity will be that of the sub-scale structure, and  thus potentially not enough to stabilize a solution on the length scale of the coarse mesh. To address this problem, one can define the keyword DhStruct in *constant/fluidRegion/phaseProperties/dragModels.(nameOfPhase).structure.(nameOfCellZones)*. This keyword defines the hydraulic diameter of the whole porous structure (viz., the dimension of the assembly, if using baffles to model wrappers, or of the entire core). The code uses it to make sure the turbulent viscosity results in a laminar Reynolds number (defaulted to 500).
<br><br>While some approaches to model k and epsilon for two-phase flow simulations are presently included in the code. In particular, the Lahey model (see *LaheyKEpsilon.H*) and a mixture model (see *mixtureKEpsilon.H*) can be uses for clear-fluids, or for mixed clear-fluid and porous-medium simulations when in case of strongly advective two-phase flow scenarios where turbulent mixing mat be neglected. In addition, as simple extension of the *porousKEpsilon* model has been implemented that allows to correct the turbulent intensity using a term that is proportional to the fraction of the other phase (see *porousKEpsilon2PhaseCorrected.H*). 
<br><br>One can find a detailed, commented example for a porous one-phase simulation in the tutorial 
[3D_SmallESFR](https://gitlab.com/foam-for-nuclear/GeN-Foam/-/tree/master/Tutorials/3D_SmallESFR/rootCase/constant/fluidRegion/turbulenceProperties).
</div>
<br>

## Initial and boundary conditions

Initial and boundary conditions adopt the usual OpenFOAM logic for one- and two-phase solvers. A couple of things to be kept in mind:
<UL>
<LI> The pressure field we solve for is p_rgh (pressure minus the gravitational head)
<LI> When performing turbulent analyses, one need to add the fields *k*, *epsilon*, nut and *alphat*
</UL>

One thing that instead specific to GeN-Foam (except for the one-phase legacy sub-solver) and that one needs to keep in mind is that U (or u.(name of fluid)) are the real velocities, not the Darcy velocities. In a porous structure, they represent the actual velocity of the fluid, and not the velocity multiplied by the fluid fraction. For instance, U will increase when transiting from a high- to a low-porosity region.


OpenFOAM provides most of the boundary conditions one may need for thermal-hydraulics models. In addition, a few boundary conditions have been included in GeN-Foam and can be found in *GeN-Foam/classes/thermalHydraulics/src/boundaryConditions*. Information on the use of each boundary condition can be found in the header files (.H). 


<div class="border-box" style='padding:0.1em; margin-left: 4em;  margin-right: 8em;  border: 1px solid gray; background-color:#f2f3fa; color:#05134a'>
<b>Setting the initial power</b>

There are several ways to set the power in GeN-Foam. For the power generated in subscale structures:
<UL>
<LI> The thermal-hydraulics solver will normally use the *powerDensity.* fields (for instance, *powerDensity.nuclearFuelPin* for pin-based reactors) that it finds in the 0 (or *startTime*) folder.
<LI> As an alternative, one can provide cellZone-by-cellZone values via the keyword *powerDensity* in the various power models in the *phase* properties (see for instance [1D_boiling](https://gitlab.com/foam-for-nuclear/GeN-Foam/-/blob/master/Tutorials/1D_boiling/constant/fluidRegion/phaseProperties)). However, if GeN-Foam finds the corresponding field in  the 0 (or *startTime*) folder, this will take priority.
</UL>
For the power generated in the fluid itself:
<UL>
<LI> The thermal-hydraulics solver will normally use the *powerDensity* field that it finds in the 0 (or *startTime*) folder.
<LI> One can override this behavior by using the *initialPowerDensity* keyword in the *phaseProperties* (see [1D_MSR_pointKinetics](https://gitlab.com/foam-for-nuclear/GeN-Foam/-/blob/master/Tutorials/1D_MSR_pointKinetics/rootCase/constant/fluidRegion/phaseProperties)) for an example. Also in this case the field in the 0 (or *startTime*) folder, this will take priority.
</UL>
If neutronics if activated, the neutronic sub-solver will overwrite everything with the power it calculates. For eigenvalue calculations, this power is set in the *pTarget* keyword in the *reactorState* dictionary. For transients, the power is a result of calculations. There is one important exception to this behavior: the point kinetics solver will only rescale the power densities (see below about why plural) that it finds in *neutroRegion*, or, if it does not find it, the one(s) that it finds in *fluidRegion*. The rescaled power density will be written to both *neutroRegion* and *fluidRegion*. For point kinetics, the *pTarget* keyword in *reactorState* is not used by the solver itslef. However, to correctly plot point kinetics results, pTarget must be consistent with the mentioned power densities.
<p>
NB: In two-phase simulations with liquid fuel, the powerDensity in neutronics goes to anything that is liquid in thermal-hydraulics. You are supposed to have one liquid and one gas. Otherwise, power will be counted twice.
</div>

<br>

<div class="border-box" style='padding:0.1em; margin-left: 4em;  margin-right: 8em;  border: 1px solid gray; background-color:#f2f3fa; color:#05134a'>
<b>Power densities and secondary power densities, and liquid fuel</b>

The spatial neutronics solvers always create a *powerDensity* and a *secondaryPowerDensity* fields. By default, *secondaryPowerDensity* is set to zero and the *fuelFraction* keyword in *nuclearData* is used to translate the volume-average power density that is normally calculated by multiplying cross-sections and fluxes into the fuel-averaged power density that is needed by the thermal-hydraulic sub-solver
<br><br> However, a *secondaryPowerDensity* might sometimes be needed. It might be used to provide some power to the coolant in a solid-fuel reactor and, more important, to provide some power to the graphite in a liquid-fuel reactor. In order to calculate a *secondaryPowerDensity*, GeN-Foam needs to know how much of the total power goes into the *secondaryPowerDensity*, and what is the volume fraction of the secondary power-producing structure or liquid. This can be done by using the *fractionToSecondaryPower* and *secondaryPowerVolumeFraction* keywords in each cellZone in nuclearData (the same place as *fuelFraction*). If these keywords are present, GeN-Foam will calculate power densities as follows:
<UL>
<LI> secondaryPowerDenisty_ = powerDensity_  / max(secondaryPowerVolumeFraction, SMALL) * fractionToSecondaryPower  ;
<LI> powerDensity_ /= max(fuelFraction, SMALL) * (1.0 - fractionToSecondaryPower);
</UL>
When the *liquidFuel* flag is set to false, the thermal-hydraulic sub-solver will:
<UL>
<LI> take the *powerDensity* field from neutronics and project it to its own *powerDensityNeutronics* field;
<LI> take the *secondaryPowerDensity* field from neutronics and project it to its own *powerDensityNeutronicsToLiquid* field.
</UL>
When the *liquidFuel* flag is set to true, the thermal-hydraulic sub-solver will:
<UL>
<LI> take the *powerDensity* field from neutronics and project it to its own *powerDensityNeutronicsToLiquid* field;
<LI> take the *secondaryPowerDensity* field from neutronics and project it to its own *powerDensityNeutronics* field.
</UL>
When point kinetics is used, the solver will simply rescale the *powerDensity* and *secondaryPowerDensity* it finds, and the thermal-hydraulic solver will take them depending on the *liquidFuel* flag as described above. The only exception is when *liquidFuel* is true and the *initialPowerDensity* keyword is used. In this case, *initialPowerDensity* will take priority and this is the value that GeN-Foam will rescale and print to the powerDensityToLqiuid


NB1: The power densities in the thermal-hydraulic sub-solver are ALWAYS the physical ones: for instance, when the *nuclearFuelPin* model is used for pin-based reactors, *powerDensity* refers to the power density inside the fuel matrix. For liquid fuel, the *powerDensity* is the power density in the liquid. They are not the power densities smeared over the whole volume.

</div>


## Discretization and solution

Details for discretization and solution of equations are handled in a standard OpenFOAM way, i.e., through the *fvSolution* and *fvSchemes* dictionaries in *constant/fluidRegion*. 

© All rights reserved. ECOLE POLYTECHNIQUE FEDERALE DE LAUSANNE, Switzerland, 2021

