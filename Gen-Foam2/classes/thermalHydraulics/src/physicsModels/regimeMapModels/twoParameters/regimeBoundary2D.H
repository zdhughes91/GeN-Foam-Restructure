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
    Foam::regimeMapModels::regimeBoundary2D

Description
    2D line segment representing the boundary between two different
    regimeDomains.

    Author: 
    Stefan Radman (stefan.radman@epfl.ch / stefanradman92@gmail.com)

SourceFiles
    twoParameters.C

\*---------------------------------------------------------------------------*/

#ifndef regimeBoundary2D_H
#define regimeBoundary2D_H

#include "dictionary.H"
#include "dimensionedScalar.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "HashPtrTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

namespace regimeMapModels
{

//- Forward declaration
class regimeDomain2D;

/*---------------------------------------------------------------------------*\
                           Class regimeBoundary2D Declaration
\*---------------------------------------------------------------------------*/

class regimeBoundary2D
{
private:

    //- Start and end points of the segment in parameter space (x is
    //  parameter1, y is parameter2)
    const Vector2D<scalar>* p0Ptr_;
    const Vector2D<scalar>* p1Ptr_;

    //- Oriented segment vector
    Vector2D<scalar> v_;

    //- Oriented segment versor
    Vector2D<scalar> vNorm_;

    //- Oriented segment normal (according to right-hand-rule)
    Vector2D<scalar> vNormPerp_;

    //- Segment length
    scalar len_;

    //- regimeDomain that owns the segment (i.e. pointed by vNormPerp_)
    const regimeDomain2D* ownerPtr_;

    //- regimeDomain on the other side of the segment (i.e. opposite to 
    //  vNormPerp_)
    const regimeDomain2D* nbrPtr_;

public:

    //- Runtime type information
    TypeName("regimeBoundary2D");

    //- Constructors

        //- Default constructor required if using this class to template a List
        //  (which regimeDomain2D does)
        regimeBoundary2D();

        //- Construct from points
        regimeBoundary2D
        (
            const Vector2D<scalar>* p0Ptr,
            const Vector2D<scalar>* p1Ptr
        );

        //- Destructor
        virtual ~regimeBoundary2D(){}//- Equivalence operators

    //- Member Functions

        //- Re-calculate non-ptr members, useful after point normalization
        void correct();

        //-
        void setOwner(const regimeDomain2D& owner)
        {
            ownerPtr_ = &owner;
        }

        //-
        void setNbr(const regimeDomain2D& nbr)
        {
            nbrPtr_ = &nbr;
        }

        //- Returns label of owner regimeDomain, if set, otherwise returns -1
        label ownerId() const;

        //- Returns label of nbr regimeDomain, if set, otherwise returns -1
        label nbrId() const;

        //- Distance of the line segment p0-p1 to a point p in parameter space
        scalar distanceTo(const Vector2D<scalar>& p) const;

        //- Access

            const Vector2D<scalar>& p0() const
            {
                return *p0Ptr_;
            }

            const Vector2D<scalar>& p1() const
            {
                return *p1Ptr_;
            }

            const Vector2D<scalar>& v() const
            {
                return v_;
            }

            const Vector2D<scalar>& vNorm() const
            {
                return vNorm_;
            }

            const Vector2D<scalar>& vNormPerp() const
            {
                return vNormPerp_;
            }

            const scalar& len() const
            {
                return len_;
            }

            //- Careful when you call these, if the ptr is not set (via an 
            //  earlier call of setOnwer, setNbr)
            
            const regimeDomain2D& owner() const
            {
                return *ownerPtr_;
            }

            const regimeDomain2D& nbr() const
            {
                return *nbrPtr_;
            }

    //- Operators

        bool operator==(const regimeBoundary2D& rhs) const;
        bool operator!=(const regimeBoundary2D& rhs) const
        {
            return !((*this) == rhs);
        }
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace regimeMapModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //