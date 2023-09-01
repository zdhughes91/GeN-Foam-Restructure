/*---------------------------------------------------------------------------*\
|       ______          _   __           ______                               |
|      / ____/  ___    / | / /          / ____/  ____   ____ _   ____ ___     |
|     / / __   / _ \  /  |/ /  ______  / /_     / __ \ / __ `/  / __ `__ \    |
|    / /_/ /  /  __/ / /|  /  /_____/ / __/    / /_/ // /_/ /  / / / / / /    |
|    \____/   \___/ /_/ |_/          /_/       \____/ \__,_/  /_/ /_/ /_/     |
|    Copyright (C) 2015 - 2022 EPFL                                           |
|                                                                             |
|    Built on OpenFOAM v2212                                                  |
|    Copyright 2011-2016 OpenFOAM Foundation, 2017-2022 OpenCFD Ltd.         |
-------------------------------------------------------------------------------
License
    This file is part of GeN-Foam.

    GeN-Foam is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    GeN-Foam is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    This offering is not approved or endorsed by the OpenFOAM Foundation nor
    OpenCFD Limited, producer and distributor of the OpenFOAM(R)software via
    www.openfoam.com, and owner of the OPENFOAM(R) and OpenCFD(R) trademarks.

    This particular snippet of code is developed according to the developer's
    knowledge and experience in OpenFOAM. The users should be aware that
    there is a chance of bugs in the code, though we've thoroughly test it.
    The source code may not be in the OpenFOAM coding style, and it might not
    be making use of inheritance of classes to full extent.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Class
    Foam::regimeMapModels::twoParameters

Description
    Regime map model that depends on two parameters.

    Dictionary example
    \verbatim
    regimeMaps
    {
        "mapName"
        {
            type twoParameters;
            
            // These correspond to the desired parameters for the 2-D map.
            // These can only correspond to the IOobject names of 
            // tvolScalarField ype objects in the mesh object registry. No 
            // other volTypeFields are supported
            parameter1      "name1";
            parameter2      "name2";

            // This sets the size (width) of the interpolation region that
            // sits across the line that separates two different regimes. This
            // is used to avoid sharp changes in the transfer coefficients that
            // might depend on this flow regime map to avoid numerical
            // instabilities. In order to deal with the fact that the two axes
            // (or dimensions) of the regime map might refer to parameters with
            // wildy different domains (e.g. one axis corresponds to a volume
            // fraction, which spans 0-1, while the other corresponds to a 
            // mass flux, which can spans up to e.g. 10000 kg/m2/s), it is
            // difficult to define an absolute meaningful "distance" in this
            // 2-D domain. Thus, this width refers to the relative width of
            // this interpolation region compared to a normalized map wherein
            // both axes span the 0-1 range only. Thus, this width needs to be
            // comprised between 0 and 1, wherein 0 = no interpolation (not
            // recommended) and 1 = everything is interpolated (then the map
            // loses meaning). A generic value of 0.1 would serve a good 
            // starting guess.
            interpolationWidth  value;

            // This specifies how the transfer coefficients are interpolated
            // in the interpolation regions described above. The only possile
            // modes are linear or quadratic.
            interpolationMode linear; //quadratic;

            // Dictionary that holds pairs of point names and point coordinates
            // (par1_X, par2_X) for a generic point pX in the parameters space.
            // These points are used to construct the regime domains in the 
            // regime map. The names of the points can be set to anything, as
            // long as they are referenced correctly in the reigmes dictionary
            // when building the regime domains
            regimePoints
            {
                p0 (par1_0, par2_0)
                p1 (par1_1, par2_1)
                //...
                pN (par1_N, par2_N)
            }

            // Dictionary of all regime domains and the points that define
            // their boundaries. For any domain, the point list can start
            // from any point as long as the resulting domain boundaries
            // (i.e. the lines that connect one point to the next one, and
            // last back to the first) never cross each other. Both convex
            // and concave shapes are possible. Please note, however, that no
            // consistency checks are made to enforce the no-crossing
            // condition (it will be implemented some time in the future).
            // The last condition for this approach to work is with regards
            // to scenarios like the one show below, wherein the boundary
            // of a domain (AB of ABEF) intersects the straight boundary (CE)
            // of another domain (EDC) at a given point (B).
            /*
                ^ parameter2
                |
                |   
                |                     C
                |                     |\
                |                     | \
                |          A----------B  \
                |          |          |   \
                |          |          |    \
                |          |          |     \
                |          F----------E------D
                |
                i-------------------------------------------------> parameter1 
            *//*
            // In this scenario, while the domain EDC could be constructed
            // from only those three points (as CE is straight), the code 
            // requires it to be constructed as EDCB (or DCBE, CBED, EBCD, you 
            // get the idea), as if the B point is not included, the "internal
            // boundary" BE between the two domains will never be recognized
            // as such, thus leading to outright interpolation errors.
            regimes
            {
                "regime0" (p0 p1 ... pI);
                "regime1" (p7 p9 ... pJ);
                //...
                "regimeM" (p10 p27 ... pK);

            }
        }
    }
    \endverbatim

    Author: 
    Stefan Radman (sradman@protonmail.com / stefan.radman@alumni.epfl.ch)

SourceFiles
    twoParameters.C

\*---------------------------------------------------------------------------*/

#ifndef twoParameters_H
#define twoParameters_H

#include "regimeMapModel.H"
#include "regimeDomain2D.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

namespace regimeMapModels
{

/*---------------------------------------------------------------------------*\
                           Class twoParameters Declaration
\*---------------------------------------------------------------------------*/

class twoParameters
:
    public regimeMapModel
{
public:

    enum class interpolationMode
    {
        linear, 
        quadratic
    };

    static const Enum<interpolationMode> interpolationModeNames_;
    
private:

    //- Parameter names
    const word parameter1Name_;
    const word parameter2Name_;

    //- Parameter refs
    const scalarField* parameter1Ptr_;
    const scalarField* parameter2Ptr_;

    //- Table of regime points used to define regimeDomains
    typedef Vector2D<scalar> scalarVector2D;
    typedef HashTable
    <
        scalarVector2D,
        word,
        word::hash
    >   pointTable; 
    pointTable points_;

    //- Table of regimeDomains
    typedef HashTable
    <
        regimeDomain2D,
        word,
        word::hash
    >   regimeDomainTable;
    regimeDomainTable regimeDomains_; 

    //- Scale factors for normalization
    scalar DX_;
    scalar DY_;

    //- Interpolation related
    const interpolationMode interpolationMode_;
    scalar interpolationWidth_;
    bool interpolate_;
    

public:

    //- Runtime type information
    TypeName("twoParameters");

    //- Constructors

        twoParameters
        (
            const fvMesh& mesh,
            const dictionary& dict
        );

        //- Destructor
        virtual ~twoParameters(){}

    //- Member Functions

        virtual void correct();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace regimeMapModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //