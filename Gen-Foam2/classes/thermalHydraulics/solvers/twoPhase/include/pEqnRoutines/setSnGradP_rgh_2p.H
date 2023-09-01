//- Update the fixedFluxPressure BCs to ensure flux consistency
surfaceScalarField::Boundary phib(phi_.boundaryField());
phib = 
    alphaf1.boundaryField()*phi1.boundaryField() 
+   alphaf2.boundaryField()*phi2.boundaryField();
setSnGrad<fixedFluxPressureFvPatchScalarField>
(
    p_rgh_.boundaryFieldRef(),
    (
        phiHbyA.boundaryField() - phib
    )/(mesh_.magSf().boundaryField()*rAUf.boundaryField())
);
