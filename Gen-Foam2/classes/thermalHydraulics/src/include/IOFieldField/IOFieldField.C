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


\*---------------------------------------------------------------------------*/

#include "IOFieldField.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<template<class> class Field, class Type>
Foam::IOFieldField<Field,Type>::IOFieldField(const IOobject& io)
:
    regIOobject(io),
    timeIndex_(this->time().timeIndex()),
    field0Ptr_(NULL)
{
    // Temporary warning
    if (io.readOpt() == IOobject::MUST_READ_IF_MODIFIED)
    {
        WarningIn("IOFieldField::IOFieldField(const IOobject&)")
            << "IOFieldField " << name()
            << " constructed with IOobject::MUST_READ_IF_MODIFIED"
            " but IOFieldField does not support automatic rereading."
            << endl;
    }

    if
    (
        (
            io.readOpt() == IOobject::MUST_READ
         || io.readOpt() == IOobject::MUST_READ_IF_MODIFIED
        )
     || (io.readOpt() == IOobject::READ_IF_PRESENT && headerOk())
    )
    {
        //- readStream(typeName) >> *this; // there was a problem with the 
        //  operator ">>" for Field<Type> (the function "new" is missing for
        //  Field<Type>) problem work around: skip that operator, directly use 
        //  the read function and provide it with the right pointer.
        FieldField<Field, Type>::readIstream
        (
            readStream(typeName), 
            IOFieldField<Field, Type>::INew()
        );
        close();
    }
}


template<template<class> class Field, class Type>
Foam::IOFieldField<Field,Type>::IOFieldField
(
    const IOobject& io, 
    const label size
)
:
    regIOobject(io),
    timeIndex_(this->time().timeIndex()),
    field0Ptr_(NULL)
{
    // Temporary warning
    if (io.readOpt() == IOobject::MUST_READ_IF_MODIFIED)
    {
        WarningIn("IOFieldField::IOFieldField(const IOobject&, const label)")
            << "IOFieldField " << name()
            << " constructed with IOobject::MUST_READ_IF_MODIFIED"
            " but IOFieldField does not support automatic rereading."
            << endl;
    }

    if
    (
        (
            io.readOpt() == IOobject::MUST_READ
         || io.readOpt() == IOobject::MUST_READ_IF_MODIFIED
        )
     || (io.readOpt() == IOobject::READ_IF_PRESENT && headerOk())
    )
    {
        //- readStream(typeName) >> *this; // there was a problem with the 
        //  operator ">>" for Field<Type> (the function "new" is missing for
        //  Field<Type>) problem work around: skip that operator, directly use 
        //  the read function and provide it with the right pointer.
        FieldField<Field, Type>::readIstream
        (
            readStream(typeName), 
            IOFieldField<Field, Type>::INew()
        );
        close();
    }
    else
    {
        FieldField<Field,Type>::setSize(size);
    }
}


template<template<class> class Field, class Type>
Foam::IOFieldField<Field,Type>::IOFieldField
(
    const IOobject& io, 
    const FieldField<Field,Type>& f
)
:
    regIOobject(io),
    //FieldField<Field,Type>(f),
    timeIndex_(this->time().timeIndex()),
    field0Ptr_(NULL)
{
    // Temporary warning
    if (io.readOpt() == IOobject::MUST_READ_IF_MODIFIED)
    {
        WarningIn
        (
"IOFieldField::IOFieldField(const IOobject&, const FieldField<Field,Type>&)"
        )
            << "IOFieldField " << name()
            << " constructed with IOobject::MUST_READ_IF_MODIFIED"
            " but IOFieldField does not support automatic rereading."
            << endl;
    }

    if
    (
        (
                io.readOpt() == IOobject::MUST_READ
            ||  io.readOpt() == IOobject::MUST_READ_IF_MODIFIED
        )
        ||  (io.readOpt() == IOobject::READ_IF_PRESENT && headerOk())
    )
    {
        //- readStream(typeName) >> *this; // there was a problem with the 
        //  operator ">>" for Field<Type> (the function "new" is missing for
        //  Field<Type>) problem work around: skip that operator, directly use 
        //  the read function and provide it with the right pointer.
        FieldField<Field, Type>::readIstream
        (
            readStream(typeName), 
            IOFieldField<Field, Type>::INew()
        );
        close();
    }
    else
    {
        PtrList<Field<Type> >::operator=(f);
        //FieldField<Field,Type>(f);
        //FieldField<Field,Type>::operator=(f);
    }
}


template<template<class> class Field, class Type>
Foam::IOFieldField<Field,Type>::IOFieldField
(
    const IOobject& io, 
    FieldField<Field,Type>& f
)
:
    regIOobject(io),
    timeIndex_(this->time().timeIndex()),
    field0Ptr_(NULL)
{
    // Temporary warning
    if (io.readOpt() == IOobject::MUST_READ_IF_MODIFIED)
    {
        WarningIn
        (
        "IOFieldField::IOFieldField(const IOobject&, FieldField<Field,Type>&)"
        )   << "IOFieldField " << name()
            << " constructed with IOobject::MUST_READ_IF_MODIFIED"
            " but IOFieldField does not support automatic rereading."
            << endl;
    }

    //- FieldField<Field,Type>::transfer(f());//careful. This might not work.
    //  In such case, try the PtrList<Field<Type> > verions of transfer
    const tmp<FieldField<Field,Type>>& tf(f);
    FieldField<Field,Type>* fieldPtr = tf.ptr();
    PtrList<Field<Type>>::transfer(*fieldPtr);
    delete fieldPtr;

    if
    (
        (
                io.readOpt() == IOobject::MUST_READ
            ||  io.readOpt() == IOobject::MUST_READ_IF_MODIFIED
        )
        ||  (io.readOpt() == IOobject::READ_IF_PRESENT && headerOk())
    )
    {
        //- readStream(typeName) >> *this; // there was a problem with the 
        //  operator ">>" for Field<Type> (the function "new" is missing for
        //  Field<Type>) problem work around: skip that operator, directly use 
        //  the read function and provide it with the right pointer.
        FieldField<Field, Type>::readIstream
        (
            readStream(typeName), 
            IOFieldField<Field, Type>::INew()
        );
        close();
    }
}


// * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * * //

template<template<class> class Field, class Type>
Foam::IOFieldField<Field,Type>::~IOFieldField()
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

// Return old time internal field
template<template<class> class Field, class Type>
const Foam::IOFieldField<Field,Type>&
Foam::IOFieldField<Field,Type>::oldTime() const
{
    if (!field0Ptr_)
    {
        field0Ptr_ = new IOFieldField<Field, Type>(*this);
    }
    else
    {
        storeOldTime();
    }

    return *field0Ptr_;
}

// Return old time internal field
template<template<class> class Field, class Type>
Foam::IOFieldField<Field,Type>&
Foam::IOFieldField<Field,Type>::oldTime()
{
    //- static_cast<const subscaleFuel&>(*this) converts *this into 
    //  const subscaleFuel&. Applying the static_cast operator to a null 
    //  pointer will convert it to a null pointer value of the target 
    //  type
    static_cast<const IOFieldField<Field,Type>&>(*this).oldTime();
                                                        
    return *field0Ptr_;
}



// Store old-time field
template<template<class> class Field, class Type>
void Foam::IOFieldField<Field,Type>::storeOldTime() const
{
    if (timeIndex_ != this->time().timeIndex())
    {
        if (field0Ptr_)
        {
            *field0Ptr_ = *this;
        }
        timeIndex_ = this->time().timeIndex();
    }
}

template<template<class> class Field, class Type>
bool Foam::IOFieldField<Field,Type>::writeData(Ostream& os) const
{
    return (os << static_cast<const FieldField<Field,Type>&>(*this)).good();
}

// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

template<template<class> class Field, class Type>
void Foam::IOFieldField<Field,Type>::operator=
(
    const IOFieldField<Field,Type>& rhs
)
{
    FieldField<Field,Type>::operator=(rhs);
}


template<template<class> class Field, class Type>
void Foam::IOFieldField<Field,Type>::operator=
(
    const FieldField<Field,Type>& rhs
)
{
    FieldField<Field,Type>::operator=(rhs);
}


// ************************************************************************* //
