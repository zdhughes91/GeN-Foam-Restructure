/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2015-2018 OpenFOAM Foundation
     \\/     M anipulation  |
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
    Foam::saturationModels::water

Description
    Saturation model for water in the 0.01-350 *C range based on an 
    interpolation of NIST data in the form Psat = A*(T-B)^C. The NIST data
    can be found at:
    https://www.nist.gov/system/files/documents/srd/NISTIR5078-Tab1.pdf

Author
    Stefan Radman (sradman@pm.me / stefan.radman@epfl.ch)

SourceFiles
    waterSaturation.C

\*---------------------------------------------------------------------------*/

#ifndef waterSaturation_H
#define waterSaturation_H

#include "saturationModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace saturationModels
{

/*---------------------------------------------------------------------------*\
                         Class water Declaration
\*---------------------------------------------------------------------------*/

class water
:
    public saturationModel
{
private :

    const volScalarField& iT_;

    const volScalarField& p_;

    //- Coeffs
    const scalar A_;
    const scalar B_;
    const scalar C_;

public:

    TypeName("water");

    //- Constructors

        //- Construct from a dictionary
        water
        (
            const phaseChangeModel& pcm,
            const dictionary& dict, 
            const objectRegistry& objReg
        );

        //- Destructor
        virtual ~water(){}

    //- Member Functions

        virtual scalar valueTSat(const label& celli) const;

        virtual scalar valuePSat(const label& celli) const;

        virtual scalar valueLnPSat(const label& celli) const;

        virtual scalar valuePSatPrime(const label& celli) const;
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace saturationModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
