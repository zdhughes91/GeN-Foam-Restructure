// This file sets constant values, that is, values that are not affected
// by feedbacks or at least are assumed to not change with time during
// the simulation.
 
forAll(nuSigmaEff_,energyI)
{
    forAll(entries_, zoneI)
    {
        label zone = zoneI;

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
            forAll(mesh_.cellZones()[zoneId], cellIlocal)
            {
                label cellIglobal = mesh_.cellZones()[zoneId][cellIlocal];
                zone = zoneI;

                //constant
                IV_[energyI][cellIglobal] = IVList_[zone][energyI];
                chiPrompt_[energyI][cellIglobal] = chiPromptList_[zone][energyI];
                chiDelayed_[energyI][cellIglobal] = chiDelayedList_[zone][energyI];
                BetaTot_[cellIglobal] = BetaTotList_[zone];

                fuelFraction_[cellIglobal] = fuelFractionList_[zone];
                secondaryPowerVolumeFraction_[cellIglobal] = secondaryPowerVolumeFractionList_[zone];
                fractionToSecondaryPower_[cellIglobal] = fractionToSecondaryPowerList_[zone];

            }
        }
    }
    IV_[energyI].correctBoundaryConditions();
    chiPrompt_[energyI].correctBoundaryConditions();
    chiDelayed_[energyI].correctBoundaryConditions();
    BetaTot_.correctBoundaryConditions();
    fuelFraction_.correctBoundaryConditions();
    secondaryPowerVolumeFraction_.correctBoundaryConditions();
    fractionToSecondaryPower_.correctBoundaryConditions();
}
 
//- Check that XSs have been specified for all regions
{
    forAll(mesh_.cellZones().names(), i)
    {
        bool zoneInEntries(false);
        word zoneName(mesh_.cellZones().names()[i]);
        forAll(entries_, j)
        {
            word entryName(entries_[j].keyword());
            if (zoneName == entryName)
            {
                zoneInEntries = true;
                break;
            }
        }
        if (!zoneInEntries)
        {
            FatalErrorInFunction 
                << "Zone: " << zoneName << " not found in neutronicProperties!"
                << exit(FatalError);
        }
    }
}
