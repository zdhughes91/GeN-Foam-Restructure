forAll(Beta_,precI)
{

    forAll(entries_, zoneI)
    {
        const word& name = entries_[zoneI].keyword();

        label zoneId = mesh_.cellZones().findZoneID(name);

        if(zoneId != -1)
        {

            const UList<label>& addr = mesh_.cellZones()[zoneId];
    
            IndirectList<scalar>(Beta_[precI].internalField(), addr)
            = scalar(BetaList_[zoneI][precI]);

            IndirectList<scalar>(lambda_[precI].internalField(), addr)
            = scalar(lambdaList_[zoneI][precI]);
        }

    }

    Beta_[precI].correctBoundaryConditions();
    lambda_[precI].correctBoundaryConditions();

}


