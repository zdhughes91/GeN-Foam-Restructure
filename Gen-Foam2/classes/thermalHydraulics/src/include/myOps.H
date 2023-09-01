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
    N/A

Description
    A collection of functions that are currently only used for under-relaxation
    and word-porcessing purposes

    Author: 
    Stefan Radman (stefan.radman@epfl.ch / stefanradman92@gmail.com)

SourceFiles
    myOps.C


\*---------------------------------------------------------------------------*/

#ifndef MYOPS_H
#define MYOPS_H

#include "fvMesh.H"

namespace Foam
{
namespace myOps
{
    /*
    Given a word name, the function returns:
    -   if NOT on the last PIMPLE iteration:
        --  true if the keyword fvSolution.relaxationFactors.fields.name
            exists;
        --  false otherwise;
    -   if on the last PIMPLE iteration:
        --  true if the keyword 
            fvSolution.relaxationFactors.fields.(name+"Final") exists;
        --  false otherwise;
    */
    bool relax(const fvMesh& mesh, word name);

    /*
    Given a word name, the function returns:
    -   if NOT on the last PIMPLE iteration:
        --  the scalar found under 
            fvSolution.relaxationFactors.fields.name if name exists;
        --  1.0 otherwise;
    -   if on the last PIMPLE iteration:
        --  the scalar found under 
            fvSolution.relaxationFactors.fields.(name+"Final") if name 
            exists;
        --  1.0 otherwise;

    This is the standard behaviour of the relax function of a geometricField,
    however, I need to be able to read relaxation factors in this way even
    without them being tied to a geometricField object
    */
    scalar relaxationFactor(const fvMesh& mesh, word name);

    //- Store prev iter only if the field is to be relaxed. To be used before
    //  a field update if the field update is follwed by a standard OpenFOAM
    //  field.relax().
    template<class Type>
    void storePrevIterIfRelax(Type& field);

    //- Split a string or word according to a delimiter char d, and return a
    //  list of split elements
    template<class Type>
    List<Type> split(const Type& compoundItem, char d);

    //- Split a string or word according to any delimiter char in ds, and 
    //  return a list of split elements
    template<class Type>
    List<Type> split(const Type& compoundItem, List<char> ds);

    template<class Type>
    void fieldInfo(const Type& field, word units = "");
}
}

//- Templated functions must be defined in the header file too, not in myOps.C,
//  so their definitions are in a source file that is included in the header
#ifdef NoRepository
    #include "myOpsTemplates.C"
#endif

#endif

/*---------------------------------------------------------------------------*/
