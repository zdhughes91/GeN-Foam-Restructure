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

#include "twoParameters.H"
#include "addToRunTimeSelectionTable.H"
#include "myOps.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace regimeMapModels
{
    defineTypeNameAndDebug(twoParameters, 0);
    addToRunTimeSelectionTable
    (
        regimeMapModel, 
        twoParameters, 
        regimeMapModels
    );

    typedef Vector2D<scalar> scalarVector2D;
    typedef HashTable
    <
        scalarVector2D,
        word,
        word::hash
    > pointTable; 
    typedef HashTable
    <
        regimeDomain2D,
        word,
        word::hash
    >   regimeDomainTable;
}
}

const Foam::Enum
<
    Foam::regimeMapModels::twoParameters::interpolationMode
>
Foam::regimeMapModels::twoParameters::interpolationModeNames_
(
    {
        { 
            interpolationMode::linear, 
            "linear" 
        },
        { 
            interpolationMode::quadratic, 
            "quadratic" 
        }
    }
);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::regimeMapModels::twoParameters::twoParameters
(
    const fvMesh& mesh,
    const dictionary& dict
)
:
    regimeMapModel
    (
        mesh,
        dict
    ),
    parameter1Name_(word(this->lookup("parameter1"))),
    parameter2Name_(word(this->lookup("parameter2"))),
    parameter1Ptr_(nullptr),
    parameter2Ptr_(nullptr),
    interpolationMode_
    (
        interpolationModeNames_.get
        (
            this->lookupOrDefault<word>
            (
                "interpolationMode", 
                "linear"
            )
        )
    ),
    interpolationWidth_(this->get<scalar>("interpolationWidth")),
    interpolate_(interpolationWidth_ != 0.0)
{
    //- Read regimePoints and insert them in the points_ hashTable
    forAllIter
    (
        dictionary,
        this->subDict("regimePoints"),
        iter
    )
    {   
        points_.insert
        (
            iter->keyword(),
            Vector2D<scalar>
            (
                List<scalar>(iter->stream())[0],
                List<scalar>(iter->stream())[1]
            )
        );
    }

    //- Construct regimeDomains, regimeNameToLabel, regimeLabelToName
    label id(0);
    forAllIter
    (
        dictionary,
        this->subDict("regimes"),
        iter
    )
    {
        word regimeName(iter->keyword());
        wordList pointNames(wordList(iter->stream()));
        List<const Vector2D<scalar>*> domainPoints(0);
        forAll(pointNames, i)
        {
            domainPoints.append(&points_[pointNames[i]]);
        }
        regimeDomains_.insert
        (
            regimeName,
            regimeDomain2D
            (
                id,
                domainPoints
            )
        );
        regimeLabelToName_.append(regimeName);
        regimeNameToLabel_.insert(regimeName, id);
        id += 1;
    }

    //- Set regime boundaries owners and neighbours, flag regime boundaries 
    //  that are internal, for each domain (i.e. shared by two domains). The
    //  latter are only used for interpolation purposes
    forAllIter
    (
        regimeDomainTable,
        regimeDomains_,
        iter0
    )
    {
        regimeDomain2D& d0(iter0());
        forAllIter
        (
            regimeDomainTable,
            regimeDomains_,
            iter1
        )
        {
            regimeDomain2D& d1(iter1());
            if (d0.id() != d1.id())
            {
                forAll(d0.boundaries(), i)
                {
                    regimeBoundary2D& b0(d0.boundaries()[i]);
                    b0.setOwner(d0);
                    forAll(d1.boundaries(), j)
                    {
                        regimeBoundary2D& b1(d1.boundaries()[j]);
                        if (b0 == b1)
                        {
                            b0.setNbr(d1);
                            b1.setNbr(d0);
                            break;
                        }
                    }
                }
            }
        }
        d0.setInternalBoundaries();
    }

    //- Debug & development Infos
    forAllIter
    (
        regimeDomainTable,
        regimeDomains_,
        iter
    )
    {
        regimeDomain2D& d(iter());
        Info << "    regime " << iter.key() << endl;
        Info << "        ID " << d.id() << endl;
        Info << "        boundingBox " << d.boundingBox() << endl;
        Info << "        all boundaries" << endl;
        forAll(d.boundaries(), i)
        {
            regimeBoundary2D& b(d.boundaries()[i]);
            Info<< "            " << b.p0() << " " << b.p1() << " " 
                << b.ownerId() << " " << b.nbrId() << endl;
        }
        Info << "        internal boundaries" << endl;
        forAll(d.internalBoundaryPtrs(), i)
        {
            const regimeBoundary2D& b(*d.internalBoundaryPtrs()[i]);
            Info<< "            " << b.p0() << " " << b.p1() << " " 
                << b.ownerId() << " " << b.nbrId() << endl;
        }
    }

    //- Calculate scale factor to normalize regime extent in parameter space.
    //  This gives meaning to the concept of "distance" of a point in parameter
    //  space to a line segment (i.e. regime boundary) in the same parameter
    //  space, which is used for interpolation purposes
    Tuple2<Vector2D<scalar>, Vector2D<scalar>> bb0
    (
        regimeDomains_[regimeLabelToName_[0]].boundingBox()
    );
    scalar  minX=bb0.first()[0], minY=bb0.first()[1], maxX=bb0.second()[0], 
            maxY=bb0.second()[1];
    forAllConstIter
    (
        regimeDomainTable, 
        regimeDomains_, 
        iter
    )
    {
        const regimeDomain2D& d(iter());
        Tuple2<Vector2D<scalar>, Vector2D<scalar>> bb(d.boundingBox());
        minX = min(minX, bb.first()[0]);
        minY = min(minY, bb.first()[1]);
        maxX = max(maxX, bb.second()[0]);
        maxY = max(maxY, bb.second()[1]);
    }
    DX_ = maxX-minX;
    DY_ = maxY-minY;

    //- Normalize all points with respect to the width of the domains
    //  in the space of each parameter (x is parameter1, y is parameter2)
    forAllIter
    (
        pointTable,
        points_,
        iter
    )
    {
        Vector2D<scalar>& p(iter());
        p[0] /= DX_;
        p[1] /= DY_;
    }

    //- Re-calculate regime boundary members
    forAllIter
    (
        regimeDomainTable,
        regimeDomains_,
        iter
    )
    {
        regimeDomain2D& d(iter());
        forAll(d.boundaries(), i)
        {
            d.boundaries()[i].correct();
        }
    }

    //- Init regimeLabelCoeffs
    regimeLabelCoeffs_ = DynamicList<DynamicList<Tuple2<label,scalar>>>
    (
        mesh_.cells().size()
    );
    
    Info << endl;
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::regimeMapModels::twoParameters::correct()
{
    //- Lookup parameter field ptrs if not set
    if (parameter1Ptr_ == nullptr)
    {
        parameter1Ptr_ = 
            &mesh_.lookupObjectRef<volScalarField>(parameter1Name_).field();
    } 
    if (parameter2Ptr_ == nullptr)
    {
        parameter2Ptr_ = 
            &mesh_.lookupObjectRef<volScalarField>(parameter2Name_).field();
    } 
    const scalarField& x(*parameter1Ptr_);
    const scalarField& y(*parameter2Ptr_);

    int N(mesh_.cells().size());

    //- Loop over all mesh cell labels
    for (int i=0; i < N; i++)
    {
        //- Reset regimeLabelCoeffs for the i-th mesh cell
        DynamicList<Tuple2<label,scalar>>& regimeLabelCoeffsi
        (
            regimeLabelCoeffs_[i]
        );
        regimeLabelCoeffsi = DynamicList<Tuple2<label,scalar>>(0);

        //- Value of the tuple parameter1-parameter2 which is used to determine
        //  the existing regime in the i-th mesh cell (normalized so distances
        //  make sense for the purposes of interpolation)
        Vector2D<scalar> p(x[i]/DX_, y[i]/DY_);

        //- Find the regime the point belongs to
        const regimeDomain2D* regimePtr = nullptr;
        forAllConstIter
        (
            regimeDomainTable,
            regimeDomains_,
            iter
        )
        {
            //- I recall that regimeLabelCoeffs_ is a list of tuples of each 
            //  mesh cell, wherein for each cell, the tuple consists of the ID
            //  of the regime that exists in said cell (first tuple element)
            //  and the coefficient weighting the contribution of said regime
            //  (which is 1.0 if only one regime exists). To check how this
            //  is factually used to interpolate the value of models of 
            //  interest, go see the interpolateValue, interpolateValueCmpt
            //  template functions in regimeMapModelTemplates.C
            if (iter().containsPoint(p))
            {
                regimePtr = &(iter());
                regimeLabelCoeffsi.append
                (
                    Tuple2<label,scalar>(regimePtr->id(), 1.0)
                );
                break;
            }
        }
        if (regimePtr != nullptr) //- Set regimeLabelCoeffs
        {
            if (interpolate_)
            {
                const regimeDomain2D& regime(*regimePtr);
                const List<const regimeBoundary2D*>& internalBoundaryPtrs
                (
                    regime.internalBoundaryPtrs()
                );
                if (internalBoundaryPtrs.size() != 0)
                {
                    //- In theory, a starting max distance of sqrt(2) should
                    //  be sufficient as these distances are in the normalized
                    //  domain, whose range for both parameters is [0, 1] (so
                    //  the max possible distance is between (0,0) and (1,1)). 
                    //  Nonetheless, I want to play it safe so 2 it is
                    scalar distance(2);
                    scalar nbrId(-1);
                    forAll(internalBoundaryPtrs, j)
                    {   
                        const regimeBoundary2D& boundaryj
                        (
                            *internalBoundaryPtrs[j]
                        );
                        scalar distancej(boundaryj.distanceTo(p));
                        if (distancej < distance)
                        {
                            distance = distancej;
                            nbrId = boundaryj.nbrId();
                        }
                    }

                    //- The coefficient is calculated so that, at the boundary,
                    //  distance = 0 and c = 0.5, regardless of the
                    //  interpolation type
                    scalar d(distance/interpolationWidth_);
                    if (d < 1.0)
                    {
                        scalar c(0.0);
                        switch (interpolationMode_)
                        {   
                            case interpolationMode::linear :
                            {
                                c = 0.5*(d + 1.0);
                                break;
                            }
                            case interpolationMode::quadratic :
                            {
                                c = 1.0-0.5*sqr(d-1.0);
                                break;
                            }                
                        }
                        regimeLabelCoeffsi[0].second() = c;
                        regimeLabelCoeffsi.append
                        (
                            Tuple2<label,scalar>(nbrId, 1.0-c)
                        );
                    }
                }
            }
            //- Debug & development infos
            
            //Info<< i << ", (" << x[i] << " " << y[i] << ") => " << p << ", "
            //    << regimePtr->id() << ", " << regimeLabelCoeffsi << endl;
        }
        else //- Point is outside regime map bounds, throw an error
        {
            FatalErrorInFunction
                << "Point (" << x[i] << ", " << y[i] << ") in parameter space "
                << "(" << parameter1Name_ << ", " << parameter2Name_ << ") "
                << "is outside the bounds of " << this->name() << "!"
                << exit(FatalError);
        }
    }
}


// ************************************************************************* //
