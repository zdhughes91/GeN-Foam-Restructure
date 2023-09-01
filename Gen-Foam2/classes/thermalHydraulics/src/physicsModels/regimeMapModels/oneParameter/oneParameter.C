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

#include "oneParameter.H"
#include "addToRunTimeSelectionTable.H"
#include "myOps.H"


// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace regimeMapModels
{
    defineTypeNameAndDebug(oneParameter, 0);
    addToRunTimeSelectionTable
    (
        regimeMapModel, 
        oneParameter, 
        regimeMapModels
    );
}
}

const Foam::Enum
<
    Foam::regimeMapModels::oneParameter::interpolationMode
>
Foam::regimeMapModels::oneParameter::interpolationModeNames_
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

Foam::regimeMapModels::oneParameter::oneParameter
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
    parameterPtr_(nullptr),
    names_(0),
    thresholds_(0),
    interpolationFlags_(0),
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
    )
{
    //- Read the regime names and the bounds from the dict
    List<List<scalar>> bounds(0);
    forAllIter
    (
        dictionary,
        this->subDict("regimeBounds"),
        iter
    )
    {   
        names_.append(iter->keyword());  
        bounds.append
        (
            List<scalar>(iter->stream())
        );
    }

    //- Order bounds and names into ascending threshold values, and
    //  re-arrange the data so that something in a format e.g.:
    /*
            regimeBounds
            {
                "regime0" (0.5  1);
                "regime1" (1.2  1.5);
                "regime2" (0    0.37);
                "regime4" (1    1.1);
            }

        is converted into two ordered lists, ordered by threshold value, i.e.:
        
            thresholds_ = [0, 0.37, 0.5, 1, 1.1, 1.2, 1.5]
            names_ = 
            [
                "regime2", 
                "",
                "regime0",
                "regime4",
                "",
                "regime1"
            ]
            interpolationFlags = 
            [
                false,
                true,
                false,
                false,
                true,
                false
            ]
        so that names_[i] is bounded by thresholds_[i] and
        thresholds_[i+1]. If the upper threshold of regime i and the lower
        threshold of regime i+1 mismatch, an nameless interpolation regime is
        created. 
    */

    //- Re-arrange individual bounds so that leftmost value is smaller
    //  than the rightmost one
    for (int i = 0; i < bounds.size(); i++)
    {
        if (bounds[i][0] > bounds[i][1])
        {
            scalar tmp = bounds[i][0];
            bounds[i][0] = bounds[i][1];
            bounds[i][1] = tmp;
        }
    }
    
    //- Then do the re-arranging as described before
    List<scalar> tmpBound;
    word tmpName;
    bool swap(true);
    while (swap)
    {
        swap = false;
        for (int i = 0; i < bounds.size()-1; i++)
        {
            if (bounds[i][0] > bounds[i+1][0])
            {
                tmpBound = bounds[i];
                bounds[i] = bounds[i+1];
                bounds[i+1] = tmpBound;
                tmpName = names_[i];
                names_[i] = names_[i+1];
                names_[i+1] = tmpName;
                swap = true;
                break;
            }
        }
    }

    //- Reset regimeLabelToName and regimeLabelCoeffs
    forAll(names_, i)
    {
        word name(names_[i]);
        regimeNameToLabel_.set
        (
            name,
            i
        );
        regimeLabelToName_.append(name);
    }
    regimeLabelCoeffs_ = DynamicList<DynamicList<Tuple2<label,scalar>>>
    (
        mesh_.cells().size()
    );

    //- Assemble all. Note that thresholds is longer than names_ and 
    //  interpolationFlags by 1
    wordList tmpNames(names_);
    names_ = wordList(0);
    names_.append(tmpNames[0]);
    thresholds_.append(-1e69);
    interpolationFlags_.append(false);
    for (int i = 1; i < bounds.size(); i++)
    {
        if (bounds[i-1][1] != bounds[i][0])
        {
            names_.append("");
            thresholds_.append(bounds[i-1][1]);
            interpolationFlags_.append(true);
        }
        names_.append(tmpNames[i]);
        thresholds_.append(bounds[i][0]);
        interpolationFlags_.append(false);
    }
    thresholds_.append(1e69);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::regimeMapModels::oneParameter::correct()
{
    if (parameterPtr_ == nullptr)
    {
        parameterPtr_ = 
            &mesh_.lookupObjectRef<volScalarField>
            (
                word(this->lookup("parameter"))
            ).field();
    } 
    const scalarField& p(*parameterPtr_);

    int N(mesh_.cells().size());

    //- Reset (faster like this when compared to using a mesh-size constructor)
    for (int i=0; i < N; i++)
    {
        regimeLabelCoeffs_[i] = DynamicList<Tuple2<label,scalar>>(0);
    }

    //-
    forAll(names_, i)
    {
        bool interpolated(interpolationFlags_[i]);
        
        const scalar& t0(thresholds_[i]);
        const scalar& t1(thresholds_[i+1]);
        scalar dt(t1-t0);
        scalarField t1_pByDt((t1-p)/dt);
        scalarField p_t0ByDt((p-t0)/dt);
        
        if (!interpolated)
        {
            label regimeLabel(regimeNameToLabel_[names_[i]]);
            
            for (int j = 0; j < N; j++)
            {
                if (t1_pByDt[j] > 0 and p_t0ByDt[j] >= 0)
                {
                    regimeLabelCoeffs_[j].append
                    (
                        Tuple2<label,scalar>(regimeLabel, 1.0)
                    );
                }
            }
        }
        else
        {
            label regime0Label(regimeNameToLabel_[names_[i-1]]);
            label regime1Label(regimeNameToLabel_[names_[i+1]]);

            switch (interpolationMode_)
            {
                case interpolationMode::linear :
                {
                    for (int j = 0; j < N; j++)
                    {
                        const scalar& c0(t1_pByDt[j]);
                        const scalar& c1(p_t0ByDt[j]);
                        if (c0 >= 0 and c1 > 0)
                        {
                            DynamicList<Tuple2<label,scalar>>& regimeLabelCoeffsj
                            (
                                regimeLabelCoeffs_[j]
                            );
                            regimeLabelCoeffsj.append
                            (
                                Tuple2<label,scalar>(regime0Label, c0)
                            );
                            regimeLabelCoeffsj.append
                            (
                                Tuple2<label,scalar>(regime1Label, c1)
                            );
                        }
                    }
                    break;
                }
                case interpolationMode::quadratic :
                {
                    scalar mid(t0+dt/2.0);
                    for (int j = 0; j < N; j++)
                    {
                        scalar c0(p_t0ByDt[j]);
                        scalar c1(t1_pByDt[j]);
                        if (c0 >= 0 and c1 > 0)
                        {
                            c0 = 
                            (
                                (p[j] <= mid) ?
                                1.0-2.0*sqr(c0) :
                                2.0*sqr(c1)
                            );
                            c1 = 1.0-c0;
                            DynamicList<Tuple2<label,scalar>>& regimeLabelCoeffsj
                            (
                                regimeLabelCoeffs_[j]
                            );
                            regimeLabelCoeffsj.append
                            (
                                Tuple2<label,scalar>(regime0Label, c0)
                            );
                            regimeLabelCoeffsj.append
                            (
                                Tuple2<label,scalar>(regime1Label, c1)
                            );
                        }
                    }
                    break;
                }
            }
        }
    }
    /*
    forAll(orderedRegimeNames_, i)
    {
        regime& regime = regimes_[orderedRegimeNames_[i]]();
        labelList& cellList(regime.cellList());
        cellList = labelList(0);
        const scalar& t0(thresholds_[i]);
        const scalar& t1(thresholds_[i+1]);
        scalar dt(t1-t0);
        scalarField t1_pByDt((t1-p)/dt);
        scalarField p_t0ByDt((p-t0)/dt);
        bool isCurrentlyPresent = false;//- I can't use the omonymous method
                                        //  of the regime class as the len of
                                        //  its cellList is still 0 at this
                                        //  point
        //- This is equivalent to 
        //  isPresent = ((max(t1_pByDt) > 0) and max(p_t0ByDt[j]) >= 0)); 
        //  but faster
        forAll(mesh_.cells(), j)
        {
            if (t1_pByDt[j] > 0  and p_t0ByDt[j] >= 0)
            {
                isCurrentlyPresent = true;
                break;
            }
        }

        scalar f(myOps::relaxationFactor(mesh_, "regime"));
        scalarField& cellField(regime.cellField());
        if (isCurrentlyPresent)
        {
            cellField = pos0(p_t0ByDt)*pos(t1_pByDt);
            if (regime.isInterpolated())
            {
                switch (interpolationMode_)
                {
                    case interpolationMode::linear :
                    {
                        if (f == 1.0)
                        {
                            forAll(mesh_.cells(), j)
                            {
                                if (cellField[j] == 1)
                                {
                                    cellList.append(j);
                                    regime.coeffs1()[j] = t1_pByDt[j];
                                    regime.coeffs2()[j] = p_t0ByDt[j];
                                }
                            }
                        }
                        else
                        {
                            scalarField& C1(regime.coeffs1());
                            scalarField& C2(regime.coeffs2());
                            forAll(mesh_.cells(), j)
                            {
                                if (cellField[j] == 1)
                                {
                                    cellList.append(j);
                                    scalar& c1(C1[j]);
                                    scalar& c2(C2[j]);
                                    c1 = f*t1_pByDt[j] + (1.0-f)*c1;
                                    c2 = 1.0-c1;
                                }
                            }
                        }
                        break;
                    }
                    case interpolationMode::quadratic :
                    {
                        scalar mid(t0+dt/2.0);
                        if (f == 1.0)
                        {
                            forAll(mesh_.cells(), j)
                            {
                                if (cellField[j] == 1)
                                {
                                    cellList.append(j);
                                    scalar& c1(regime.coeffs1()[j]);
                                    scalar& c2(regime.coeffs2()[j]);
                                    c1 = 
                                    (
                                        (p[j] <= mid) ?
                                        1.0-2.0*sqr(p_t0ByDt[j]) :
                                        2.0*sqr(t1_pByDt[j])
                                    );
                                    c2 = 1.0-c1;
                                }
                            }
                        }
                        else
                        {
                            forAll(mesh_.cells(), j)
                            {
                                if (cellField[j] == 1)
                                {
                                    cellList.append(j);
                                    scalar& c1(regime.coeffs1()[j]);
                                    scalar& c2(regime.coeffs2()[j]);
                                    c1 = 
                                    (
                                        f*
                                        (
                                            (p[j] <= mid) ?
                                            1.0-2.0*sqr(p_t0ByDt[j]) :
                                            2.0*sqr(t1_pByDt[j])
                                        ) 
                                    +   (1.0-f)*c1
                                    );
                                    c2 = 1.0-c1;
                                }
                            }
                        }
                        break;
                    }
                }
            }
            else
            {
                forAll(mesh_.cells(), j)
                {
                    if (cellField[j] == 1)
                    {
                        cellList.append(j);
                    }
                }
            }
        }
        else
        {
            cellField *= 0;
        }
    }

    //- Update the requiresModelCorrection flag. The way this is done is
    //  not related to the specific run-time-selectable regimeMapModel,
    //  so its wrapped in a function in the base class
    this->setRequiresModelCorrection();
    */
}


// ************************************************************************* //
