/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright held by original author
     \\/     M anipulation  |
-------------------------------------------------------------------------------

-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Class
    Foam::XS

Description
    Class for XS parametrization and other nuclear data.

Author of this file (and associated .C or included .H files)
    Carlo Fiorina <carlo.fiorina@outlook.com; carlo.fiorina@epfl.ch;>, EPFL (Switzerland)

SourceFiles
    XS.C

\*---------------------------------------------------------------------------*/

#ifndef XS_H
#define XS_H

#include "volFields.H"
#include "rhoThermo.H"
#include "fvMatrices.H"
#include "autoPtr.H"
#include "runTimeSelectionTables.H"
#include "coordinateSystem.H"
#include "fvCFD.H"
#include "turbulentFluidThermoModel.H"
#include "volPointInterpolation.H"
#include "meshToMesh.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{


/*---------------------------------------------------------------------------*\
                         Class thermoMechanicsSolver Declaration
\*---------------------------------------------------------------------------*/

class XS
{
    // Private data

    // Private Member Functions
        
        //- Disallow default bitwise copy construct
        XS(const XS&);

        //- Disallow default bitwise assignment
        void operator=(const XS&);
        
protected:
    
    // Protected data

        //- Reference to mesh
            const fvMesh& mesh_;

        //Data
            IOdictionary nuclearData_;
            IOdictionary nuclearDataRadialExp_;
            IOdictionary nuclearDataAxialExp_;
            IOdictionary nuclearDataFuelTemp_;
            IOdictionary nuclearDataRhoCool_;
            IOdictionary nuclearDataTCool_;
            IOdictionary nuclearDataCladExp_;

            label energyGroups_;
            label precGroups_;
            label legendreMoments_;
        
        //Data fields
            PtrList<volScalarField> IV_;
            PtrList<volScalarField> D_;
            PtrList<volScalarField> nuSigmaEff_;
            PtrList<volScalarField> sigmaPow_;
            PtrList<volScalarField> sigmaDisapp_;
            List<List<PtrList<volScalarField > > > sigmaFromTo_;
            PtrList<volScalarField> chiPrompt_;
            PtrList<volScalarField> chiDelayed_;
            PtrList<volScalarField> Beta_;
            volScalarField BetaTot_;
            PtrList<volScalarField> lambda_;
            volScalarField fuelFraction_;
            volScalarField secondaryPowerVolumeFraction_; // Volume fraction of secondary power-producing structure
                                                         // e.g., graphite in MSRs
            volScalarField fractionToSecondaryPower_;  // Fraction of total power that goes to secondary
                                                       // power-producing structure
            PtrList<volScalarField > discFactor_;
            scalar ScNo_;
        
        //XS
            PtrList<entry> entries_;
            label zoneNumber_;

        //Disc factors
            // fuel fraction used to determine vol power IN the fuel starting from vol power from sigmaPow
            PtrList<scalar > fuelFractionList_;
            PtrList<scalar > secondaryPowerVolumeFractionList_;
            PtrList<scalar > fractionToSecondaryPowerList_;
            // adjust or not DF
            PtrList<bool > dfAdjustList_;
            //discontinuity factors
            PtrList<scalarField > discFactorList_;
            //integral fluxes for adapting disc factors
            PtrList<scalarField > integralFluxList_;

        //Feedbacks Coefficients

            bool fastNeutrons_;
            bool adjustDiscFactors_;
            bool useGivenDiscFactors_;
            List<int> groupsWoDF_;
            List<int> doNotParametrize_;
            //Constants
            PtrList<scalarField > IVList_;
            PtrList<scalarField > chiPromptList_;
            PtrList<scalarField > chiDelayedList_;
            PtrList<scalarField > BetaList_;
            PtrList<scalar > BetaTotList_;
            PtrList<scalarField > lambdaList_;
            //Nominal 
            PtrList<scalarField > DList_;
            PtrList<scalarField > nuSigmaEffList_;
            PtrList<scalarField > sigmaPowList_;
            PtrList<scalarField > sigmaDisappList_;
            List<PtrList<scalarSquareMatrix > > sigmaFromToList_; 
            // fuel temperature coefficient
            scalar TfuelRef_;
            scalar TfuelPerturbed_;
            PtrList<scalarField > fuelTempDList_;
            PtrList<scalarField > fuelTempNuSigmaEffList_;
            PtrList<scalarField > fuelTempSigmaPowList_;
            PtrList<scalarField > fuelTempSigmaDisappList_;
            List<PtrList<scalarSquareMatrix > > fuelTempSigmaFromToList_;
            // axial expansion coefficient
            scalar AxExp_;
            PtrList<scalarField > axialExpDList_;
            PtrList<scalarField > axialExpNuSigmaEffList_;
            PtrList<scalarField > axialExpSigmaPowList_;
            PtrList<scalarField > axialExpSigmaDisappList_;
            List<PtrList<scalarSquareMatrix > > axialExpSigmaFromToList_;
            // radial expansion coefficient
            scalar RadExp_;
            vector axialOrientation_;
            PtrList<scalarField > radialExpDList_;
            PtrList<scalarField > radialExpNuSigmaEffList_;
            PtrList<scalarField > radialExpSigmaPowList_;
            PtrList<scalarField > radialExpSigmaDisappList_;
            List<PtrList<scalarSquareMatrix > > radialExpSigmaFromToList_;
            // coolant expansion coefficient
            scalar rhoCoolRef_;
            scalar rhoCoolPerturbed_;
            PtrList<scalarField > rhoCoolDList_;
            PtrList<scalarField > rhoCoolNuSigmaEffList_;
            PtrList<scalarField > rhoCoolSigmaPowList_;
            PtrList<scalarField > rhoCoolSigmaDisappList_;
            List<PtrList<scalarSquareMatrix > > rhoCoolSigmaFromToList_;
            // coolant temperature coefficient
            scalar TCoolRef_;
            scalar TCoolPerturbed_;
            PtrList<scalarField > TCoolDList_;
            PtrList<scalarField > TCoolNuSigmaEffList_;
            PtrList<scalarField > TCoolSigmaPowList_;
            PtrList<scalarField > TCoolSigmaDisappList_;
            List<PtrList<scalarSquareMatrix > > TCoolSigmaFromToList_;
            // cladding expansion coefficient
            scalar TcladRef_;
            scalar TcladPerturbed_;
            PtrList<scalarField > cladExpDList_;
            PtrList<scalarField > cladExpNuSigmaEffList_;
            PtrList<scalarField > cladExpSigmaPowList_;
            PtrList<scalarField > cladExpSigmaDisappList_;
            List<PtrList<scalarSquareMatrix > > cladExpSigmaFromToList_;

        //For control rods
            IOdictionary CRmove_;
            PtrList<entry> CRentries_;
            label CRNumber_;
            PtrList<scalar > CRstart_;
            PtrList<scalar > CRfinish_;
            PtrList<scalar > CRspeed_;
            PtrList<word >   CRFollowerName_;
            PtrList<scalar > CRinitialPosition_;
            PtrList<scalar > CRposition_;
            PtrList<scalar > initialDistanceFromMeshCR_;

public:

    //- Runtime type information
    TypeName("XS");
    
    // Declare run-time constructor selection table

    declareRunTimeSelectionTable
    (
        autoPtr,
        XS,
        dictionary,
        (const fvMesh& mesh),
        (mesh)
    );    
    // Access

        label energyGroups() const
        {
            return energyGroups_;
        }
        label precGroups() const
        {
            return precGroups_;
        }
        label legendreMoments() const
        {
            return legendreMoments_;
        }
        const PtrList<volScalarField> IV() const
        {
            return IV_;
        }
        const PtrList<volScalarField> D() const
        {
            return D_;
        }
        const PtrList<volScalarField> nuSigmaEff() const
        {
            return nuSigmaEff_;
        }
        const PtrList<volScalarField> sigmaPow() const
        {
            return sigmaPow_;
        }
        const PtrList<volScalarField> sigmaDisapp() const
        {
            return sigmaDisapp_;
        }
        const List<List<PtrList<volScalarField > > > sigmaFromTo() const
        {
            return sigmaFromTo_;
        }
        const PtrList<volScalarField> chiPrompt() const
        {
            return chiPrompt_;
        }
        const PtrList<volScalarField> chiDelayed() const
        {
            return chiDelayed_;
        }
        const PtrList<volScalarField> Beta() const
        {
            return Beta_;
        }
        const volScalarField BetaTot() const
        {
            return BetaTot_;
        }
        const PtrList<volScalarField> lambda() const
        {
            return lambda_;
        }
        const volScalarField fuelFraction() const
        {
            return fuelFraction_;
        }
        const volScalarField secondaryPowerVolumeFraction() const
        {
            return secondaryPowerVolumeFraction_;
        }
        const volScalarField fractionToSecondaryPower() const
        {
            return fractionToSecondaryPower_;
        }
        const PtrList<volScalarField> discFactor() const
        {
            return discFactor_;
        }
        scalar ScNo() const
        {
            return ScNo_;
        }

    // Constructors

        //- Construct from mesh and thermo-mechanical properties
        XS
        (
            const fvMesh& mesh
        );
        
        //- Return a pointer to a new XS
        static autoPtr<XS> New
        (
            const fvMesh& mesh
        );           

    //- Destructor
        virtual ~XS();

    // Member Functions
    
        //- Update XS      

        void correct
        (
            const volScalarField& Tfuel, 
            const volScalarField& Tclad, 
            const volScalarField& rhoCool, 
            const volScalarField& TCool,
            const volVectorField& Disp             
        );

        void init();

        void adjustDiscFactors(const PtrList<volScalarField>& fluxStar);


};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
