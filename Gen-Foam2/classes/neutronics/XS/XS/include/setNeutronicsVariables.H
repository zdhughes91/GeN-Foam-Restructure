// Calculate Control Rod position
forAll(CRentries_,CRzoneI)
{
    //recalculate positions in case mesh has moved 
    label FOLLOWERzoneId = mesh_.cellZones().findZoneID(CRFollowerName_[CRzoneI]);

    if(FOLLOWERzoneId == -1)
    {
      Info << "Control rod with follower: " << CRFollowerName_(CRzoneI) << " does not exists" << endl;
    }

    scalar cellHeight = -100.0 ;
    scalar maxFollowerHeight = cellHeight ;

    volScalarField centersVSF(mesh_.C() & vector(0,0,1));

    Field<scalar> centersSF(centersVSF.internalField());

    forAll(mesh_.cellZones()[FOLLOWERzoneId], cellIlocal)
    {

        label cellIglobal = mesh_.cellZones()[FOLLOWERzoneId][cellIlocal];

        maxFollowerHeight = max(cellHeight,centersSF[cellIglobal]) ;

        cellHeight = maxFollowerHeight;

    }
    reduce(maxFollowerHeight, maxOp<scalar>());
    CRinitialPosition_[CRzoneI] = maxFollowerHeight-initialDistanceFromMeshCR_[CRzoneI];
    if(mesh_.time().value()<CRfinish_[CRzoneI])
    {
        CRposition_[CRzoneI] = (CRinitialPosition_[CRzoneI]-CRspeed_[CRzoneI]*(max(mesh_.time().value(),CRstart_[CRzoneI])-CRstart_[CRzoneI]));
    }
    else
    {
        CRposition_[CRzoneI] = (CRinitialPosition_[CRzoneI]-CRspeed_[CRzoneI]*(CRfinish_[CRzoneI]-CRstart_[CRzoneI]));
    }
    Info << "Control rod " << CRentries_[CRzoneI].keyword() << " in position: " << CRposition_[CRzoneI] << endl;
}

tmp<volVectorField> tDNeutroRad(new volVectorField(Disp));
volVectorField& DNeutroRad = tDNeutroRad.ref();

tmp<volScalarField> tDNeutroAx(new volScalarField(Disp & axialOrientation_));
volScalarField& DNeutroAx = tDNeutroAx.ref();

DNeutroRad.replace(0,DNeutroRad.component(0)*(1-axialOrientation_[0]));
DNeutroRad.replace(1,DNeutroRad.component(1)*(1-axialOrientation_[1]));
DNeutroRad.replace(2,DNeutroRad.component(2)*(1-axialOrientation_[2]));

DNeutroRad.correctBoundaryConditions();
DNeutroAx.correctBoundaryConditions();

volScalarField radExp(fvc::div(DNeutroRad)/2.0);
volScalarField axExp(fvc::grad(DNeutroAx)& axialOrientation_);

radExp.correctBoundaryConditions();
axExp.correctBoundaryConditions();

tDNeutroRad.clear();
DNeutroAx.clear();

volScalarField logT( log(max(0.01,(Tfuel/dimensionedScalar("",dimensionSet(0,0,0,1,0,0,0),TfuelRef_)))) );
volScalarField diffT( (pow(max(Tfuel,dimensionedScalar("",dimensionSet(0,0,0,1,0,0,0),SMALL)),0.5)-pow(dimensionedScalar("",dimensionSet(0,0,0,1,0,0,0),TfuelRef_),0.5)) );

forAll(nuSigmaEff_,energyI)
{

    bool parametrize = true;
    forAll(doNotParametrize_, groupI)
    {
        if(energyI == doNotParametrize_[groupI])
        {
            parametrize = false;
            Info << "Not Parametrizing XS for group: " << doNotParametrize_[groupI] << endl;
        }
    }

    forAll(entries_, zoneI)
    {
        label zone = zoneI;

        const word& name = entries_[zoneI].keyword();

        //Info << "Mesh zone: " << name  << endl;

        label zoneId = mesh_.cellZones().findZoneID(name);

        if(zoneId == -1)
        {
            Info << "WARNING! Mesh zone " << name << " does not exists" << endl;
            Info << "Ignore this warning if it is a completely withdrawn control rod" << endl;
        }
        else
        {
            //check if you are in a follower or CR
            bool isFollower = false;
            bool isCR = false;
            label CRzone = 0;

            forAll(CRentries_,CRzoneI)
            {
                const word& CRname = CRentries_[CRzoneI].keyword();
                
                if(name == CRname)
                {
                    isCR = true ;
                    CRzone = CRzoneI;
                }
                if(name == CRFollowerName_[CRzoneI])
                {
                    isFollower = true ;
                    CRzone = CRzoneI;
                }
            }

            forAll(mesh_.cellZones()[zoneId], cellIlocal)
            {
                label cellIglobal = mesh_.cellZones()[zoneId][cellIlocal];
                zone = zoneI;

                if(CRNumber_>0)
                {
                    if(isFollower && mesh_.C()[cellIglobal].z()>CRposition_[CRzone])
                    {
                        forAll(entries_, zoneII)
                        {
                            if(entries_[zoneII].keyword() == CRentries_[CRzone].keyword())
                            {
                                zone = zoneII;
                            }
                        }
                    }

                    if(isCR && mesh_.C()[cellIglobal].z()<CRposition_[CRzone])
                    {
                        forAll(entries_, zoneII)
                        {
                            if(entries_[zoneII].keyword() == CRFollowerName_[CRzone])
                            {
                                zone = zoneII;
                            }
                        }
                    }
                }

                // Adjustment of D and XS for current reactor conditions of:
                // rhoCool, TCool, fuelTemp, cladExp, axialExp and radialExp.

                D_[energyI][cellIglobal] = DList_[zone][energyI];
                nuSigmaEff_[energyI][cellIglobal] = nuSigmaEffList_[zone][energyI];
                sigmaPow_[energyI][cellIglobal] = sigmaPowList_[zone][energyI];
                sigmaDisapp_[energyI][cellIglobal] = sigmaDisappList_[zone][energyI];
                forAll(sigmaFromTo_,momentI)
                {
                  forAll(sigmaFromTo_[momentI],energyJ)
                  {
                      sigmaFromTo_[momentI][energyJ][energyI][cellIglobal] = sigmaFromToList_[zone][momentI][energyJ][energyI];
                  }
                }

                if (parametrize)
                {
                    D_[energyI][cellIglobal]
                 += rhoCoolDList_[zone][energyI]*(rhoCool[cellIglobal] - rhoCoolRef_)
                  + TCoolDList_[zone][energyI]*(TCool[cellIglobal] - TCoolRef_)
                  + (
                        fastNeutrons_
                      ? fuelTempDList_[zone][energyI]*logT[cellIglobal]
                      : fuelTempDList_[zone][energyI]*diffT[cellIglobal]
                    )
                  + cladExpDList_[zone][energyI]*(Tclad[cellIglobal] - TcladRef_)
                  + axialExpDList_[zone][energyI]*(axExp[cellIglobal])
                  + radialExpDList_[zone][energyI]*(radExp[cellIglobal]);

                    nuSigmaEff_[energyI][cellIglobal]
                 += rhoCoolNuSigmaEffList_[zone][energyI]*(rhoCool[cellIglobal] - rhoCoolRef_)
                  + TCoolNuSigmaEffList_[zone][energyI]*(TCool[cellIglobal] - TCoolRef_)
                  + (
                        fastNeutrons_
                      ? fuelTempNuSigmaEffList_[zone][energyI]*logT[cellIglobal]
                      : fuelTempNuSigmaEffList_[zone][energyI]*diffT[cellIglobal]
                    )
                  + cladExpNuSigmaEffList_[zone][energyI]*(Tclad[cellIglobal] - TcladRef_)
                  + axialExpNuSigmaEffList_[zone][energyI]*(axExp[cellIglobal])
                  + radialExpNuSigmaEffList_[zone][energyI]*(radExp[cellIglobal]);

                    sigmaPow_[energyI][cellIglobal]
                 += rhoCoolSigmaPowList_[zone][energyI]*(rhoCool[cellIglobal] - rhoCoolRef_)
                  + TCoolSigmaPowList_[zone][energyI]*(TCool[cellIglobal] - TCoolRef_)
                  + (
                        fastNeutrons_
                      ? fuelTempSigmaPowList_[zone][energyI]*logT[cellIglobal]
                      : fuelTempSigmaPowList_[zone][energyI]*diffT[cellIglobal]
                    )
                  + cladExpSigmaPowList_[zone][energyI]*(Tclad[cellIglobal] - TcladRef_)
                  + axialExpSigmaPowList_[zone][energyI]*(axExp[cellIglobal])
                  + radialExpSigmaPowList_[zone][energyI]*(radExp[cellIglobal]);

                    sigmaDisapp_[energyI][cellIglobal]
                 += rhoCoolSigmaDisappList_[zone][energyI]*(rhoCool[cellIglobal] - rhoCoolRef_)
                  + TCoolSigmaDisappList_[zone][energyI]*(TCool[cellIglobal] - TCoolRef_)
                  + (
                        fastNeutrons_
                      ? fuelTempSigmaDisappList_[zone][energyI]*logT[cellIglobal]
                      : fuelTempSigmaDisappList_[zone][energyI]*diffT[cellIglobal]
                    )
                  + cladExpSigmaDisappList_[zone][energyI]*(Tclad[cellIglobal] - TcladRef_)
                  + axialExpSigmaDisappList_[zone][energyI]*(axExp[cellIglobal])
                  + radialExpSigmaDisappList_[zone][energyI]*(radExp[cellIglobal]);



                  forAll(sigmaFromTo_,momentI)  
                  {
                    forAll(sigmaFromTo_[momentI],energyJ)
                    {
                        sigmaFromTo_[momentI][energyJ][energyI][cellIglobal]
                     += rhoCoolSigmaFromToList_[zone][momentI][energyJ][energyI]*(rhoCool[cellIglobal] - rhoCoolRef_)
                      + TCoolSigmaFromToList_[zone][momentI][energyJ][energyI]*(TCool[cellIglobal] - TCoolRef_)
                      + (
                            fastNeutrons_
                          ? fuelTempSigmaFromToList_[zone][momentI][energyJ][energyI]*logT[cellIglobal]
                          : fuelTempSigmaFromToList_[zone][momentI][energyJ][energyI]*diffT[cellIglobal]
                        )
                      + cladExpSigmaFromToList_[zone][momentI][energyJ][energyI]*(Tclad[cellIglobal] - TcladRef_)
                      + axialExpSigmaFromToList_[zone][momentI][energyJ][energyI]*(axExp[cellIglobal])
                      + radialExpSigmaFromToList_[zone][momentI][energyJ][energyI]*(radExp[cellIglobal]);
                    }
                  }
                }

                if(useGivenDiscFactors_ || ((mesh_.time().timeIndex()-mesh_.time().startTimeIndex())>0))
                {
                    discFactor_[energyI][cellIglobal] = discFactorList_[zone][energyI] ;
                }
                else
                {
                    discFactorList_[zone][energyI] = discFactor_[energyI][cellIglobal] ;
                }

            }//close forall cellIlocal in zoneId
        }//close else from zoneId check
    }//close forall zoneI in zone entries_


    D_[energyI].correctBoundaryConditions();
    nuSigmaEff_[energyI].correctBoundaryConditions();
    sigmaPow_[energyI].correctBoundaryConditions();
    sigmaDisapp_[energyI].correctBoundaryConditions();
    forAll(sigmaFromTo_, momentI)
    {
      forAll(sigmaFromTo_[momentI],energyJ)
      {
          sigmaFromTo_[momentI][energyJ][energyI].correctBoundaryConditions();
      }
    }
    discFactor_[energyI].correctBoundaryConditions();

}
