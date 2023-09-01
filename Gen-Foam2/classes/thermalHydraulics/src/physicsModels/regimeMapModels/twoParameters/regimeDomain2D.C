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

\*---------------------------------------------------------------------------*/

#include "regimeDomain2D.H"
#include "Random.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace regimeMapModels
{
    defineTypeNameAndDebug(regimeDomain2D, 0);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::regimeMapModels::regimeDomain2D::regimeDomain2D()
:
    id_(-1)
{}


Foam::regimeMapModels::regimeDomain2D::regimeDomain2D
(
    label id,
    List<const Vector2D<scalar>*> pointPtrs
)
:
    id_(id),
    pointPtrs_(pointPtrs),
    boundaries_(0),
    internalBoundaryPtrs_(0)
{
    //Info << "    constructing regimeDomain2D ID " << id_ << endl;
    int N(pointPtrs_.size());
    forAll(pointPtrs_, i)
    {
        const Vector2D<scalar>* p0Ptr(nullptr);
        const Vector2D<scalar>* p1Ptr(nullptr);
        if (i < N-1)
        {
            p0Ptr = pointPtrs_[i];
            p1Ptr = pointPtrs_[i+1];
        }
        else
        {
            p0Ptr = pointPtrs_[N-1];
            p1Ptr = pointPtrs_[0];
        }
        boundaries_.append
        (
            regimeBoundary2D
            (
                p0Ptr,
                p1Ptr
            )
        );
        //Info<< "        added boundary " << boundaries_[i].p0() << " " 
        //    << boundaries_[i].p1() << endl;
    }

    //- Determine polygon sign (1 = convex counter-clockwise, -1 = convex 
    //  clockwise, 0 = concave (ordering is inconsequential in that scenario))
    forAll(boundaries_, i)
    {
        const regimeBoundary2D& b0(boundaries_[i]);
        const regimeBoundary2D* b1Ptr(nullptr);
        if (i == boundaries_.size()-1)
            b1Ptr = &(boundaries_[0]);
        else
            b1Ptr = &(boundaries_[i+1]); 
        const regimeBoundary2D& b1(*b1Ptr);
        const Vector2D<scalar>& v0(b0.vNorm());
        const Vector2D<scalar>& v1(b1.vNorm());
        scalar x(v0[0]*v1[1] - v0[1]*v1[0]);
        int sign(int(x/max(mag(x), 1e-9)));
        if (i == 0)
        {
            if (sign >= 0)
                sign_ = 1;
            else
                sign_ = -1;
        }
        else 
        {
            if (sign_*sign < 0)
            {
                //- Then domain is concave
                sign_ = 0;
                break;
            }
        }
    }

    //- Construct the ray to be used for intercept testing in a manner so
    //  that is not parallel to any boundary in the domain
    Random rng(-(1922-1991)*(id+1917));
    bool flag = true;
    while (flag)
    {
        forAll(boundaries_, i)
        {
            rayNorm_ = Vector2D<scalar>
            (
                rng.sample01<scalar>(), 
                rng.sample01<scalar>()
            ).normalise();
            scalar dot(mag(rayNorm_&boundaries_[i].vNorm()));
            Info << dot << endl;
            scalar tol(1e-2);
            if ((dot < tol) or (dot > 1.0-tol))
            {
                flag = true;
                break;
            }
            else
                flag = false;
        }
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::Tuple2<Foam::Vector2D<Foam::scalar>, Foam::Vector2D<Foam::scalar>> 
Foam::regimeMapModels::regimeDomain2D::boundingBox() const
{
    const Vector2D<scalar>& p0(*pointPtrs_[0]);
    scalar minX=p0[0], minY=p0[1], maxX=p0[0], maxY=p0[1];

    forAll(pointPtrs_, i)
    {
        const Vector2D<scalar>& p(*pointPtrs_[i]);
        minX = min(minX, p[0]);
        maxX = max(maxX, p[0]);
        minY = min(minY, p[1]);
        maxY = max(maxY, p[1]);
    }

    return 
        Tuple2<Vector2D<scalar>, Vector2D<scalar>>
        (
            Vector2D<scalar>(minX, minY),
            Vector2D<scalar>(maxX, maxY)
        );
}


void Foam::regimeMapModels::regimeDomain2D::setInternalBoundaries()
{
    forAll(boundaries_, i)
    {
        const regimeBoundary2D& b(boundaries_[i]);
        if (b.nbrId() != -1)
            internalBoundaryPtrs_.append(&b);
    }
}


bool Foam::regimeMapModels::regimeDomain2D::containsPoint
(
    const Vector2D<scalar>& p
) const
{
    //- If the polygon is not concave
    if (sign_ != 0)
        return convexContainsPoint(p);
    //- Else 
    return concaveContainsPoint(p);
}


bool Foam::regimeMapModels::regimeDomain2D::convexContainsPoint
(
    const Vector2D<scalar>& p
) const
{
    //- The logic here is to perform a cross-product (wherein the thrid cmpt of
    //  the involved vectors is always 0, as these objects are 2-D) between
    //  each of the vectors that can be constructed between each boundary
    //  starting poing p0 and the point under exam p with the corresponding
    //  p0-p1 vector for each boundary, for all boundaries. If the sign of 
    //  these products is the same for all boundaries, then the point lies in
    //  the domain. This only works if the domain is convex (i.e. sign_ != 0)
    forAll(boundaries_, i)
    {
        const regimeBoundary2D& b(boundaries_[i]);
        const Vector2D<scalar> v01(b.v());
        Vector2D<scalar> v0p(p-b.p0());
        scalar x(v01[0]*v0p[1] - v01[1]*v0p[0]);
        if (sign_*x < 0)
            return false;
    }
    return true;
}


bool Foam::regimeMapModels::regimeDomain2D::concaveContainsPoint
(
    const Vector2D<scalar>& p
) const
{
    //- The idea is to count the number of times a test ray that starts at
    //  p and extends towards infinity intersects domain boundaries. If said
    //  number is odd, the point is inside the domain, otherwise it is outside.
    //  The test ray is rayNorm_. This approax works for both concave and
    //  convex polygons, but it is farily more expensive, so it is used only
    //  for concave shapes in containsPoint(p)
    int numberOfIntersections= 0;
    forAll(boundaries_, i)
    {
        const regimeBoundary2D& b(boundaries_[i]);
        const Vector2D<scalar> v01(b.v());
        Vector2D<scalar> vp0(b.p0()-p);
        Vector2D<scalar> vp1(b.p1()-p);
        scalar vxvp0(rayNorm_[0]*vp0[1]-rayNorm_[1]*vp0[0]);
        scalar vxvp1(rayNorm_[0]*vp1[1]-rayNorm_[1]*vp1[0]);
        scalar vp0xv01(vp0[0]*v01[1]-vp0[1]*v01[0]);

        //- The condition on vxvp1 is strictly greater or strictly lesser so to
        //  not consider the intersections with end-points in the count. Only
        //  intersections with start points are considered. The condition on 
        //  vp0xv010 is striclty greater than 0 so that points on the edge are
        //  considered outtside of the domain. In my experience with typical
        //  regime maps, this is the best option
        if ((vp0xv01 > 0.0 and vxvp0 <= 0.0 and vxvp1 > 0.0) or
            (vp0xv01 < 0.0 and vxvp0 >= 0.0 and vxvp1 < 0.0))
            numberOfIntersections += 1.0;
    }
    return (numberOfIntersections%2 == 1);
}

// ************************************************************************* //
