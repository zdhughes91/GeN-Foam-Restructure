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

#include "regimeMapModel.H"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class modelType, class firstArgType>
void Foam::regimeMapModel::constructModels
(
    List<autoPtr<modelType>>& models,
    const firstArgType& firstArg,
    const dictionary& dict,
    const objectRegistry& objReg
) const
{
    models = List<autoPtr<modelType>>(0);
    forAll(regimeLabelToName_, i)
    {
        word regimeName(regimeLabelToName_[i]);
        const dictionary& regimeDict(dict.subDict(regimeName));
        models.append
        (
            modelType::New
            (
                firstArg,
                regimeDict,
                objReg
            )
        );
    }
}

template<class valueType, class modelType>
valueType Foam::regimeMapModel::interpolateValue
(
    const List<autoPtr<modelType>>& models,
    const label& celli
) const
{
    //- The If-else appears to be useless, but it is to force a Return Value
    //  Optimization for the case n==1. I should check the actual performance 
    //  gain...
    const DynamicList<Tuple2<label,scalar>>& rlci
    (
        regimeLabelCoeffs_[celli]
    );
    int n(rlci.size());
    if (n == 1)
        return valueType(models[rlci[0].first()]->value(celli));
    //- else
    const Tuple2<label,scalar>& rlci0(rlci[0]);
    valueType value
    (
        rlci0.second()*models[rlci0.first()]->value(celli)
    );
    for(int j=1; j<n; j++)
    {
        //- The first element of the tuple is the regime label (which 
        //  corresponds to the label of the sub-model in models
        //  while the second element of the tuple corresponds to 
        //  the coefficient of the associated regime, used to weight the
        //  contribution of the correpsonding regime model value valuej in 
        //  the total value value
        const Tuple2<label,scalar>& rlcij(rlci[j]);
        value += rlcij.second()*models[rlcij.first()]->value(celli);
    }
    return value;
}
