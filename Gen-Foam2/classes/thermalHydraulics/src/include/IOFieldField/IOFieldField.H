/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Class
    Foam::IOFieldField

Description
    A primitive FieldField of type \<T\> with automated input and output. 
    Includes functionalities for time dependent simulations.

SourceFiles
    IOFieldField.C

Author
    Carlo Fiorina (carlo.fiorina@outlook.com; carlo.fiorina@epfl.ch)

\*---------------------------------------------------------------------------*/

#ifndef IOFieldField_H
#define IOFieldField_H

#include "regIOobject.H"
#include "Field.H"
#include "FieldField.H"
#include "Time.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class IOFieldField Declaration
\*---------------------------------------------------------------------------*/

template<template<class> class Field, class Type>
class IOFieldField
:
    public regIOobject,
    public FieldField<Field,Type>
{

//- Next few lines needed to overcome the problem described in line 54 of 
//  IOFieldField.C
public:

    typedef Field<Type> typeField;

    //- Class INew<T> is the existing one but does not work with Fields. Create
    //  class INew (nested in IOFieldField) that works (only) with Field<Type>
    class INew 
    {
    public:
        //- Construct from Istream
        tmp<typeField> operator()(Foam::Istream& is) const 
        {
            return tmp<typeField>
            (
                //- The original one gives: T::new(is) where T is a type 
                // (Field<Type> for me) and the problem is that "new" is not 
                //  defined for Field<Type>. "Field" has a contructor from 
                //  istream
                new typeField(is) 
                                   
            );
        }
    };

private:

    //- Current time index.
    //  Used to trigger the storing of the old-time value
    mutable label timeIndex_;

    //- Pointer to old time field
    mutable IOFieldField<Field,Type>* field0Ptr_;

public:

    TypeName("FieldField");


    // Constructors

        //- Construct from IOobject
        IOFieldField(const IOobject&);

        //- Construct from IOobject and size (does not set values)
        IOFieldField(const IOobject&, const label size);

        //- Construct from components
        IOFieldField(const IOobject&, const FieldField<Field,Type>&);

        //- Construct by transferring the FieldField contents
        IOFieldField(const IOobject&, FieldField<Field,Type> &);


    //- Destructor
    virtual ~IOFieldField();


    // Member functions

        //- Return non-const old time field
        //  (Not a good idea but it is used for sub-cycling)
        IOFieldField<Field,Type>& oldTime();

        //- Return old time field
        const IOFieldField<Field,Type>& oldTime() const;

        //- Store the old-time field
        void storeOldTime() const;

        bool writeData(Ostream&) const;


    // Member operators

        void operator=(const IOFieldField<Field,Type>&);

        void operator=(const FieldField<Field,Type>&);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

//- When you create a template, the compiler reads and discards everything 
//  unless you use the template.
#include "IOFieldField.C"  


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
