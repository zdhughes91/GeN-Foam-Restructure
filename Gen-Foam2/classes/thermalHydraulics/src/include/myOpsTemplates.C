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

namespace Foam
{
namespace myOps
{
    template<class Type>
    void storePrevIterIfRelax(Type& field)
    {
        if (myOps::relax(field.mesh(), field.name()))
            field.storePrevIter();
    }

    template<class Type>
    List<Type> split(const Type& compoundItem, char d)
    {
        List<Type> items(1, "");
        int i = 0;
        for (char const &c : compoundItem)
        {
            if (c == d)
            {
                i += 1;
                items.append("");
            }
            else items[i] += c;
        }
        return items;
    }

    template<class Type>
    List<Type> split(const Type& compoundItem, List<char> ds)
    {
        List<Type> items(1, "");
        int i = 0;
        int j = 0;
        int compoundItemSize(0);
        for (char const &c : compoundItem)
        {
            (void) c;
            compoundItemSize += 1;
        }
        for (char const &c : compoundItem)
        {
            bool addChar = true;
            for (char const &d : ds)
            {
                if (c == d)
                {
                    if 
                    (
                        items[items.size()-1] != ""
                    and j != compoundItemSize-1)
                    {
                        i += 1;
                        items.append("");
                    }
                    addChar = false;
                    break;
                }
            }
            j += 1;
            if (addChar)
                items[i] += c;
        }
        return items;
    }

    template<class Type>
    void fieldInfo(const Type& field, word units)
    {
        Info<< field.name() << " (avg min max) = "
        << field.weightedAverage(field.mesh().V()).value()
        << " " << min(field).value()
        << " " << max(field).value()
        << " " << units << endl;
    }
}
}
