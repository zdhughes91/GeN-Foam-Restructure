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
    Foam::regimeMapModels::regimeDomain2D

Description
    Regime domain representation in 2D parameter-space.

    Autor: 
    Stefan Radman (stefan.radman@epfl.ch / stefanradman92@gmail.com)

SourceFiles
    twoParameters.C

\*---------------------------------------------------------------------------*/

#ifndef regimeDomain2D_H
#define regimeDomain2D_H

#include "regimeBoundary2D.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

namespace regimeMapModels
{

/*---------------------------------------------------------------------------*\
                           Class regimeDomain2D Declaration
\*---------------------------------------------------------------------------*/

class regimeDomain2D
{
private:

    //- Domain label
    label id_;

    //- List of points that make up the domain (relatively useless)
    List<const Vector2D<scalar>*> pointPtrs_;

    //- List of boundaries that make up the domain
    List<regimeBoundary2D> boundaries_;

    //- List of boundaries in boundaries_ that are shared by this and another
    //  regimeDomain. Used for interpolation purposes
    List<const regimeBoundary2D*> internalBoundaryPtrs_;

    //- If positive, the polygon is convex with counter-clockwise ordered 
    //  points. If negaative, it is convex with clockwise ordered points.
    //  if 0 then the polygon is concave
    int sign_;

    //-
    Vector2D<scalar> rayNorm_;

public:

    //- Runtime type information
    TypeName("regimeDomain2D");

    //- Constructors

        //- Default
        regimeDomain2D();

        //- Constructor
        regimeDomain2D
        (
            label id,
            List<const Vector2D<scalar>*> pointPtrs
        );

        //- Destructor
        virtual ~regimeDomain2D(){}

    //- Member Functions

        //-
        Tuple2<Vector2D<scalar>, Vector2D<scalar>> boundingBox() const;

        //- 
        void setInternalBoundaries();

        //- Check whether point p is contained within the domain
        bool containsPoint(const Vector2D<scalar>& p) const;

        //- Check whether point p is contained within the domain, but works
        //  only if the domain is convex (regardless of point orientation)
        bool convexContainsPoint(const Vector2D<scalar>& p) const;

        //- Check whether point p is contained within the domain. Also works
        //  for convex polygons though more expensive than the other
        //  formulation and so it should really be used only for concave 
        //  polygons
        bool concaveContainsPoint(const Vector2D<scalar>& p) const;

        //- Access

            const label& id() const
            {
                return id_;
            }

            //- Non const-access as well to allow the use of the setNbr 
            //  function of regimeBoundary2D
            List<regimeBoundary2D>& boundaries()
            {
                return boundaries_;
            }

            const List<regimeBoundary2D>& boundaries() const
            {
                return boundaries_;
            }

            const List<const regimeBoundary2D*>& internalBoundaryPtrs() const
            {
                return internalBoundaryPtrs_;
            }
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace regimeMapModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //