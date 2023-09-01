

forAll(fluxStar,energyI)
{
    bool adjust = true;
    forAll(groupsWoDF_, groupI)
    {
        if(energyI == groupsWoDF_[groupI])
        {
            adjust = false;
            Info << "Not calculating DF for group: " << groupsWoDF_[groupI] << endl;
        }
    }
    if(adjust)
    {
        volScalarField integrateYesNo
        (
            IOobject
            (
                "integrateYesNo",
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh_,
            dimensionedScalar("", dimensionSet(0,0,0,0,0,0,0),0.0),
            zeroGradientFvPatchScalarField::typeName
        );
        scalar totFluxSerpent = 0;
        scalar totFluxGeNFoam = 0;

        forAll(entries_, zoneI)
        {
            
            const word& name = entries_[zoneI].keyword();

            label zoneId = mesh_.cellZones().findZoneID(name);

            if(zoneId == -1)
            {

            }
            else
            {
                if(dfAdjustList_[zoneI]==true)
                {
        
                    forAll(mesh_.cellZones()[zoneId], cellIlocal)
                    {
        
                        label cellIglobal = mesh_.cellZones()[zoneId][cellIlocal];
            
                        integrateYesNo[cellIglobal] = 1.0;
            
                    }
                    totFluxSerpent += integralFluxList_[zoneI][energyI];
                
                }   
    
            }
        }

        totFluxGeNFoam = fvc::domainIntegrate(integrateYesNo*fluxStar[energyI]).value();

        integrateYesNo *= 0.0 ;
        forAll(entries_, zoneI)
        {

            const word& name = entries_[zoneI].keyword();

            Info << "Mesh zone: " << name  << endl;
    
            label zoneId = mesh_.cellZones().findZoneID(name);

            if(zoneId == -1)
            {
                Info << "Warning: mesh zone: " << name << " does not exists" << endl;
                Info << "Ignore this warning if it is a completely withdrawn control rod" << endl;
            }
            else
            {
                if(dfAdjustList_[zoneI]==true)
                {
    
                    forAll(mesh_.cellZones()[zoneId], cellIlocal)
                    {
        
                        label cellIglobal = mesh_.cellZones()[zoneId][cellIlocal];
    
                        integrateYesNo[cellIglobal] = 1.0;
    
                    }

                    discFactorList_[zoneI][energyI] =  (fvc::domainIntegrate(integrateYesNo*fluxStar[energyI]).value()/totFluxGeNFoam ) / (integralFluxList_[zoneI][energyI]/totFluxSerpent);
                }
                
                Info << "discFactorList: zone " <<  entries_[zoneI].keyword() << ", energy " << energyI << ": " << discFactorList_[zoneI][energyI]  << endl;
    
                integrateYesNo *= 0.0 ;     
    
            }
        }
    }

}


