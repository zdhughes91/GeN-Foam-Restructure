/*---------------------------------------------------------------------------*\
|       ______          _   __           ______                               |
|      / ____/  ___    / | / /          / ____/  ____   ____ _   ____ ___     |
|     / / __   / _ \  /  |/ /  ______  / /_     / __ \ / __ `/  / __ `__ \    |
|    / /_/ /  /  __/ / /|  /  /_____/ / __/    / /_/ // /_/ /  / / / / / /    |
|    \____/   \___/ /_/ |_/          /_/       \____/ \__,_/  /_/ /_/ /_/     |
|    Copyright (C) 2015 - 2022 EPFL                                           |
|                                                                             |
|    Built on OpenFOAM v2112                                                  |
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
    Foam::fvMatrixExt

Description
    Extended fvMatrix to have matrix multiplication and source with effect 
    of boundaries and including already dividion by cell volume. Works only 
    for scalars!

    Author: 
    Carlo Fiorina <carlo.fiorina@outlook.com; carlo.fiorina@epfl.ch;>


\*---------------------------------------------------------------------------*/

#ifndef fvMatrixExt_H
#define fvMatrixExt_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                            Class fvMatrixExt Declaration
\*---------------------------------------------------------------------------*/
template<class Type>
class fvMatrixExt
:
    public fvMatrix<Type>
{

public:

// Constructors

// Construct from mesh and dictionary
fvMatrixExt(const tmp<fvMatrix<Type> >& tfvm):
    fvMatrix<Type>(tfvm)
{}

//- Destructor
~fvMatrixExt()
{}

// Member Functions

Foam::tmp< GeometricField<Type, fvPatchField, volMesh> > Source() const
{

    tmp<GeometricField<Type, fvPatchField, volMesh> > tSource
    (
        new GeometricField<Type, fvPatchField, volMesh>
        (
            IOobject
            (
                "Source("+this->psi().name()+')',
                this->psi().instance(),
                this->psi().mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            this->psi().mesh(),
            this->dimensions()/dimVol,
            zeroGradientFvPatchScalarField::typeName
        )
    );
    GeometricField<Type, fvPatchField, volMesh>& Source = tSource();

    Source.primitiveFieldRef() = this->source();
    Source.primitiveFieldRef() /= this->psi().mesh().V();
    Source.correctBoundaryConditions();

    return tSource;
        
}

Foam::tmp< GeometricField<Type, fvPatchField, volMesh> > SourceWithBC() const //beta version, may give problems
{
    tmp<GeometricField<Type, fvPatchField, volMesh> > tSource
    (
        new GeometricField<Type, fvPatchField, volMesh>
        (
            IOobject
            (
                "Source("+this->psi().name()+')',
                this->psi().instance(),
                this->psi().mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            this->psi().mesh(),
            this->dimensions()/dimVol,
            zeroGradientFvPatchScalarField::typeName
        )
    );
    GeometricField<Type, fvPatchField, volMesh>& Source = tSource.ref();

    Source.primitiveFieldRef() = this->source();
    Foam::fvMatrix<Type>::addBoundarySource(Source.primitiveFieldRef(),false);//false exclude processor boundaries
    Source.primitiveFieldRef() /= this->psi().mesh().V();
    Source.correctBoundaryConditions(); 
    return tSource;
}

Foam::tmp< GeometricField<Type, fvPatchField, volMesh> > AmulWithBC() 
{
    tmp<GeometricField<Type, fvPatchField, volMesh> > tAmulWithBC
    (
        new GeometricField<Type, fvPatchField, volMesh>
        (
            IOobject
            (
                "AmulWithBC("+this->psi().name()+')',
                this->psi().instance(),
                this->psi().mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            this->psi().mesh(),
            this->dimensions()/dimVol,
            zeroGradientFvPatchScalarField::typeName
        )
    );
    GeometricField<Type, fvPatchField, volMesh>& AmulWithBC = tAmulWithBC.ref();

    this->Foam::lduMatrix::Amul(AmulWithBC.primitiveFieldRef(), this->psi().internalField(), this->boundaryCoeffs(),  this->psi().boundaryField().scalarInterfaces() ,  direction(0));

    Field<Type> boundaryDiagCmpt(this->psi().size(), 0.0);
    Foam::fvMatrix<Type>::addBoundaryDiag(boundaryDiagCmpt, false);
    AmulWithBC.primitiveFieldRef() += boundaryDiagCmpt*this->psi().internalField();
    AmulWithBC.primitiveFieldRef() /= this->psi().mesh().V();
    AmulWithBC.correctBoundaryConditions();
    return tAmulWithBC;     
}



Foam::tmp< GeometricField<Type, fvPatchField, volMesh> > AmulWithBC(GeometricField<Type, fvPatchField, volMesh>& phi) 
{
    tmp<GeometricField<Type, fvPatchField, volMesh> > tAmulWithBC
    (
        new GeometricField<Type, fvPatchField, volMesh>
        (
            IOobject
            (
                "AmulWithBC("+phi.name()+')',
                phi.instance(),
                phi.mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            phi.mesh(),
            this->dimensions()/dimVol,
            zeroGradientFvPatchScalarField::typeName
        )
    );
    GeometricField<Type, fvPatchField, volMesh>& AmulWithBC = tAmulWithBC.ref();


    this->Foam::lduMatrix::Amul(AmulWithBC.primitiveFieldRef(), phi.internalField(), this->boundaryCoeffs(),  phi.boundaryField().scalarInterfaces() ,  direction(0));

    Field<Type> boundaryDiagCmpt(phi.size(), 0.0);
    Foam::fvMatrix<Type>::addBoundaryDiag(boundaryDiagCmpt, false);
    AmulWithBC.primitiveFieldRef() += boundaryDiagCmpt*phi.internalField();
    AmulWithBC.primitiveFieldRef() /= phi.mesh().V();
    AmulWithBC.correctBoundaryConditions();
    return tAmulWithBC;     
}

Foam::tmp< GeometricField<Type, fvPatchField, volMesh> > AmulWithoutBC(GeometricField<Type, fvPatchField, volMesh>& phi) 
{
    tmp<GeometricField<Type, fvPatchField, volMesh> > tAmulWithBC
    (
        new GeometricField<Type, fvPatchField, volMesh>
        (
            IOobject
            (
                "AmulWithBC("+phi.name()+')',
                phi.instance(),
                phi.mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            phi.mesh(),
            this->dimensions()/dimVol,
            zeroGradientFvPatchScalarField::typeName
        )
    );
    GeometricField<Type, fvPatchField, volMesh>& AmulWithBC = tAmulWithBC.ref();


    this->Foam::lduMatrix::Amul(AmulWithBC.primitiveFieldRef(), phi.internalField(), this->boundaryCoeffs(),  phi.boundaryField().scalarInterfaces() ,  direction(0));

    AmulWithBC.primitiveFieldRef() /= phi.mesh().V();
    AmulWithBC.correctBoundaryConditions();
    return tAmulWithBC;     
}

Foam::tmp< GeometricField<Type, fvPatchField, volMesh> > SourceWithBCnoV() const //beta version, may give problems
{
    tmp<GeometricField<Type, fvPatchField, volMesh> > tSource
    (
        new GeometricField<Type, fvPatchField, volMesh>
        (
            IOobject
            (
                "Source("+this->psi().name()+')',
                this->psi().instance(),
                this->psi().mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            this->psi().mesh(),
            this->dimensions(),
            zeroGradientFvPatchScalarField::typeName
        )
    );
    GeometricField<Type, fvPatchField, volMesh>& Source = tSource.ref();

    Source.primitiveFieldRef() = this->source();
    Foam::fvMatrix<Type>::addBoundarySource(Source.primitiveFieldRef(),false);//false exclude processor boundaries
    Source.correctBoundaryConditions();
    return tSource;
}

Foam::tmp< GeometricField<Type, fvPatchField, volMesh> > AmulWithBCnoV() 
{
    tmp<GeometricField<Type, fvPatchField, volMesh> > tAmulWithBC
    (
        new GeometricField<Type, fvPatchField, volMesh>
        (
            IOobject
            (
                "AmulWithBC("+this->psi().name()+')',
                this->psi().instance(),
                this->psi().mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            this->psi().mesh(),
            this->dimensions(),
            zeroGradientFvPatchScalarField::typeName
        )
    );
    GeometricField<Type, fvPatchField, volMesh>& AmulWithBC = tAmulWithBC.ref();


    this->Foam::lduMatrix::Amul(AmulWithBC.primitiveFieldRef(), this->psi().internalField(), this->boundaryCoeffs(),  this->psi().boundaryField().scalarInterfaces() ,  direction(0));

    Field<Type> boundaryDiagCmpt(this->psi().size(), 0.0);
    Foam::fvMatrix<Type>::addBoundaryDiag(boundaryDiagCmpt, false);
    AmulWithBC.primitiveFieldRef() += boundaryDiagCmpt*this->psi().internalField();
    AmulWithBC.correctBoundaryConditions();
    return tAmulWithBC;
}

Foam::tmp< GeometricField<Type, fvPatchField, volMesh> > AmulWithBCnoV(GeometricField<Type, fvPatchField, volMesh>& phi) 
{
    tmp<GeometricField<Type, fvPatchField, volMesh> > tAmulWithBC
    (
        new GeometricField<Type, fvPatchField, volMesh>
        (
            IOobject
            (
                "AmulWithBC("+phi.name()+')',
                phi.instance(),
                phi.mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            phi.mesh(),
            this->dimensions(),
            zeroGradientFvPatchScalarField::typeName
        )
    );
    GeometricField<Type, fvPatchField, volMesh>& AmulWithBC = tAmulWithBC.ref();


    this->Foam::lduMatrix::Amul(AmulWithBC.primitiveFieldRef(), phi.internalField(), this->boundaryCoeffs(),  phi.boundaryField().scalarInterfaces() ,  direction(0));

    Field<Type> boundaryDiagCmpt(phi.size(), 0.0);
    Foam::fvMatrix<Type>::addBoundaryDiag(boundaryDiagCmpt, false);
    AmulWithBC.primitiveFieldRef() += boundaryDiagCmpt*phi.internalField();
    AmulWithBC.correctBoundaryConditions();
    return tAmulWithBC;         
}

Type initialResidual() 
{
    scalarField saveDiag(this->diag());

    Field<Type> Source(this->source());

    Foam::fvMatrix<Type>::addBoundarySource(Source,true);//false exclude processor boundaries

     scalarField psi = this->psi().primitiveField();
    Foam::fvMatrix<Type>::addBoundaryDiag(this->diag(), 0);

    FieldField<Field, scalar> bouCoeffs
    (
        this->boundaryCoeffs()
    );

    FieldField<Field, scalar> intCoeffs
    (
        this->internalCoeffs()
    );

    lduInterfaceFieldPtrsList interfaces =
        this->psi().boundaryField().scalarInterfaces();

    const label startRequest = Pstream::nRequests();

    this->initMatrixInterfaces
    (
        true,
        bouCoeffs,
        interfaces,
        psi,
        Source,
        0
    );

    this->updateMatrixInterfaces
    (
        true,
        bouCoeffs,
        interfaces,
        psi,
        Source,
        0,
        startRequest
    );

    label nCells = psi.size();

    scalarField wA(nCells);
    scalarField pA(nCells);

    this->Foam::fvMatrix<Type>::Amul(wA, psi, this->boundaryCoeffs(),  this->psi().boundaryField().scalarInterfaces() ,  direction(0));

    // --- Calculate initial residual field
    scalarField rA(Source - wA);

    // --- Calculate normalisation factor

    this->Foam::fvMatrix<Type>::sumA(pA, bouCoeffs, interfaces);

    pA *= gAverage(psi,this->Foam::lduMatrix::mesh().comm());

    //this->Foam::lduMatrix::mesh().comm()

    scalar normFactor = 
    gSum
    (
        (mag(wA - pA) + mag(Source - pA))(),
        this->Foam::lduMatrix::mesh().comm()
    )
    + SMALL;


        this->diag() = saveDiag;

        return gSumMag(rA,this->Foam::lduMatrix::mesh().comm())  / normFactor;
            
    }



};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
