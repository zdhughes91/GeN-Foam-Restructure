//- Correct fixed-flux phi BCs to be consistent with the velocity BCs
FieldField<fvsPatchField, scalar> U1f
(
    mesh_.Sf().boundaryField() & U1.boundaryField()
);
surfaceScalarField::Boundary& phi1bf = phi1.boundaryFieldRef();
forAll(mesh_.boundary(), patchi)
{
    if (isA<fixedValueFvsPatchScalarField>(phi1bf[patchi]))
        phi1bf[patchi] == U1f[patchi];
}

FieldField<fvsPatchField, scalar> U2f
(
    mesh_.Sf().boundaryField() & U2.boundaryField()
);
surfaceScalarField::Boundary& phi2bf = phi2.boundaryFieldRef();
forAll(mesh_.boundary(), patchi)
{
    if (isA<fixedValueFvsPatchScalarField>(phi2bf[patchi]))
        phi2bf[patchi] == U2f[patchi];
}
