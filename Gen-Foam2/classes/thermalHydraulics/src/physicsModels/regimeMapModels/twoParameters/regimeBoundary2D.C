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

#include "regimeBoundary2D.H"
#include "regimeDomain2D.H"


// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace regimeMapModels
{
    defineTypeNameAndDebug(regimeBoundary2D, 0);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::regimeMapModels::regimeBoundary2D::regimeBoundary2D()
:
    p0Ptr_(nullptr),
    p1Ptr_(nullptr),
    v_(Vector2D<scalar>(0.0, 0.0)),
    vNorm_(Vector2D<scalar>(0.0, 0.0)),
    vNormPerp_(Vector2D<scalar>(0.0, 0.0)),
    len_(0.0),
    ownerPtr_(nullptr),
    nbrPtr_(nullptr)
{
}


Foam::regimeMapModels::regimeBoundary2D::regimeBoundary2D
(
    const Vector2D<scalar>* p0Ptr,
    const Vector2D<scalar>* p1Ptr
)
:
    p0Ptr_(p0Ptr),
    p1Ptr_(p1Ptr),
    v_(p1()-p0()),
    vNorm_((p1()-p0()).normalise()),
    vNormPerp_(vNorm_[1], -vNorm_[0]),
    len_(Foam::mag(v_)),
    ownerPtr_(nullptr),
    nbrPtr_(nullptr)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::regimeMapModels::regimeBoundary2D::correct()
{
    v_ = p1()-p0();
    vNorm_ = (p1()-p0()).normalise();
    vNormPerp_ = Vector2D<scalar>(vNorm_[1], -vNorm_[0]);
    len_ = Foam::mag(v_);
}


Foam::label Foam::regimeMapModels::regimeBoundary2D::ownerId() const
{
    if (ownerPtr_ == nullptr)
        return -1;
    return ownerPtr_->id();
}


Foam::label Foam::regimeMapModels::regimeBoundary2D::nbrId() const
{
    if (nbrPtr_ == nullptr)
        return -1;
    return nbrPtr_->id();
}


Foam::scalar Foam::regimeMapModels::regimeBoundary2D::distanceTo
(
    const Vector2D<scalar>& p
) const
{    
    Vector2D<scalar> v0p(p-p0());
    scalar v0pParLen(v0p & vNorm_);
    if (v0pParLen < 0)
        return Foam::mag(v0p);
    Vector2D<scalar> v1p(p-p1());
    scalar v1pParLen(v1p & vNorm_);
    if (v1pParLen > 0)
        return Foam::mag(v1p);
    return Foam::mag(v0p-v0pParLen*vNorm_);
}


bool Foam::regimeMapModels::regimeBoundary2D::operator==
(
    const regimeBoundary2D& rhs
) const
{
    if 
    (
        (p0() == rhs.p0() and p1() == rhs.p1())
    or  (p0() == rhs.p1() and p1() == rhs.p0())
    )
        return true;
    return false;
}


// ************************************************************************* //
