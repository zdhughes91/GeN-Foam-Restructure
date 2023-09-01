//- Correct fixed-flux phi BCs to be consistent with the velocity BCs
FieldField<fvsPatchField, scalar> Uf
(
    mesh_.Sf().boundaryField() & U.boundaryField()
);
surfaceScalarField::Boundary& phibf = fluid_.phi().boundaryFieldRef();
forAll(mesh_.boundary(), patchi)
{
    if (isA<fixedValueFvsPatchScalarField>(phibf[patchi]))
        phibf[patchi] == Uf[patchi];
}
