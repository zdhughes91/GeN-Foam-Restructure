Info<< "\nReading nuclear data\n" << endl;

PtrList<entry> entriesRadialExp(nuclearDataRadialExp_.lookup("zones"));
PtrList<entry> entriesAxialExp(nuclearDataAxialExp_.lookup("zones"));
PtrList<entry> entriesFuelTemp(nuclearDataFuelTemp_.lookup("zones"));
PtrList<entry> entriesRhoCool(nuclearDataRhoCool_.lookup("zones"));
PtrList<entry> entriesTCool(nuclearDataTCool_.lookup("zones"));
PtrList<entry> entriesCladExp(nuclearDataCladExp_.lookup("zones"));


//Nuclear data affected by feedbacks
//RadExp
PtrList<scalarField > DRadialExpList(zoneNumber_);
PtrList<scalarField > nuSigmaEffRadialExpList(zoneNumber_);
PtrList<scalarField > sigmaPowRadialExpList(zoneNumber_);
PtrList<scalarField > sigmaDisappRadialExpList(zoneNumber_);
List<PtrList<scalarSquareMatrix > > sigmaFromToRadialExpList(zoneNumber_);  //scattering matrix
//AxExp
PtrList<scalarField > DAxialExpList(zoneNumber_);
PtrList<scalarField > nuSigmaEffAxialExpList(zoneNumber_);
PtrList<scalarField > sigmaPowAxialExpList(zoneNumber_);
PtrList<scalarField > sigmaDisappAxialExpList(zoneNumber_);
List<PtrList<scalarSquareMatrix > > sigmaFromToAxialExpList(zoneNumber_);   //scattering matrix
//fuelTemp
PtrList<scalarField > DFuelTempList(zoneNumber_);
PtrList<scalarField > nuSigmaEffFuelTempList(zoneNumber_);
PtrList<scalarField > sigmaPowFuelTempList(zoneNumber_);
PtrList<scalarField > sigmaDisappFuelTempList(zoneNumber_);
List<PtrList<scalarSquareMatrix > > sigmaFromToFuelTempList(zoneNumber_);   //scattering matrix
//RhoCool
PtrList<scalarField > DRhoCoolList(zoneNumber_);
PtrList<scalarField > nuSigmaEffRhoCoolList(zoneNumber_);
PtrList<scalarField > sigmaPowRhoCoolList(zoneNumber_);
PtrList<scalarField > sigmaDisappRhoCoolList(zoneNumber_);
List<PtrList<scalarSquareMatrix > > sigmaFromToRhoCoolList(zoneNumber_);    //scattering matrix
//TCool
PtrList<scalarField > DTCoolList(zoneNumber_);
PtrList<scalarField > nuSigmaEffTCoolList(zoneNumber_);
PtrList<scalarField > sigmaPowTCoolList(zoneNumber_);
PtrList<scalarField > sigmaDisappTCoolList(zoneNumber_);
List<PtrList<scalarSquareMatrix > > sigmaFromToTCoolList(zoneNumber_);  //scattering matrix
//cladExp
PtrList<scalarField > DCladExpList(zoneNumber_);
PtrList<scalarField > nuSigmaEffCladExpList(zoneNumber_);
PtrList<scalarField > sigmaPowCladExpList(zoneNumber_);
PtrList<scalarField > sigmaDisappCladExpList(zoneNumber_);
List<PtrList<scalarSquareMatrix > > sigmaFromToCladExpList(zoneNumber_);    //scattering matrix

//Fill the data
forAll(entries_,zoneI)
{
    dictionary& dict = entries_[zoneI].dict();


    const word& name = entries_[zoneI].keyword();

    Info << "name: " << name << endl;

    fuelFractionList_.set(zoneI,new scalar(dict.lookupOrDefault("fuelFraction",1.0)));
    secondaryPowerVolumeFractionList_.set(zoneI,new scalar(dict.lookupOrDefault("secondaryPowerVolumeFraction",1.0)));
    fractionToSecondaryPowerList_.set(zoneI,new scalar(dict.lookupOrDefault("fractionToSecondaryPower",0.0)));

    dfAdjustList_.set(zoneI,new bool(dict.lookupOrDefault("dfAdjust",true)));

    //Nuclear data (assumed to be) not affected by feedbacks
    IVList_.set(zoneI,new scalarField("IV",dict,energyGroups_));
    chiPromptList_.set(zoneI,new scalarField("chiPrompt",dict,energyGroups_));
    chiDelayedList_.set(zoneI,new scalarField("chiDelayed",dict,energyGroups_));
    discFactorList_.set(zoneI,new scalarField("discFactor",dict,energyGroups_));
    integralFluxList_.set(zoneI,new scalarField("integralFlux",dict,energyGroups_));
    BetaList_.set(zoneI,new scalarField("Beta",dict,precGroups_));
    BetaTotList_.set(zoneI,new scalar(sum(BetaList_[zoneI])));
    lambdaList_.set(zoneI,new scalarField("lambda",dict,precGroups_));
    Info << "Reading nominal " << endl;
    //Nominal
    DList_.set(zoneI,new scalarField("D",dict,energyGroups_));
    nuSigmaEffList_.set(zoneI,new scalarField("nuSigmaEff",dict,energyGroups_));
    sigmaPowList_.set(zoneI,new scalarField("sigmaPow",dict,energyGroups_));
    sigmaDisappList_.set(zoneI,new scalarField("sigmaDisapp",dict,energyGroups_));
    sigmaFromToList_[zoneI].setSize(legendreMoments_);
    forAll(sigmaFromToList_[zoneI],momentI)
    {
      sigmaFromToList_[zoneI].set(momentI,new scalarSquareMatrix(dict.lookup("scatteringMatrix"+std::string("P")+Foam::name(momentI))));
    }
    Info << "Finishing reading nominal " << endl;
    //RadExp
    bool read = false;
    int perturbedEntry = 0;
    forAll(entriesRadialExp,entryI)
    {
      if(entriesRadialExp[entryI].keyword() == name)
      {
        read = true;
        perturbedEntry = entryI;
      }
    }
    if(read)
    {
      dictionary& dictRadialExp = entriesRadialExp[perturbedEntry].dict();
      DRadialExpList.set(zoneI,new scalarField("D",dictRadialExp,energyGroups_));
      nuSigmaEffRadialExpList.set(zoneI,new scalarField("nuSigmaEff",dictRadialExp,energyGroups_));
      sigmaPowRadialExpList.set(zoneI,new scalarField("sigmaPow",dictRadialExp,energyGroups_));
      sigmaDisappRadialExpList.set(zoneI,new scalarField("sigmaDisapp",dictRadialExp,energyGroups_));
      sigmaFromToRadialExpList[zoneI].setSize(legendreMoments_);
      forAll(sigmaFromToRadialExpList[zoneI],momentI)
      {
        sigmaFromToRadialExpList[zoneI].set(momentI,new scalarSquareMatrix(dictRadialExp.lookup("scatteringMatrix"+std::string("P")+Foam::name(momentI))));
      }
    }
    else
    {
      DRadialExpList.set(zoneI,new scalarField(DList_[zoneI]));
      nuSigmaEffRadialExpList.set(zoneI,new scalarField(nuSigmaEffList_[zoneI]));
      sigmaPowRadialExpList.set(zoneI,new scalarField(sigmaPowList_[zoneI]));
      sigmaDisappRadialExpList.set(zoneI,new scalarField(sigmaDisappList_[zoneI]));
      sigmaFromToRadialExpList[zoneI].setSize(legendreMoments_);
      forAll(sigmaFromToRadialExpList[zoneI],momentI)
      {
        sigmaFromToRadialExpList[zoneI].set(momentI,new scalarSquareMatrix(sigmaFromToList_[zoneI][momentI]));
      }
    }
    //AxExp
    read = false;
    perturbedEntry = 0;
    forAll(entriesAxialExp,entryI)
    {
      if(entriesAxialExp[entryI].keyword() == name)
      {
        read = true;
        perturbedEntry = entryI;
      }
    }
    if(read)
    {
      dictionary& dictAxialExp = entriesAxialExp[perturbedEntry].dict();
      DAxialExpList.set(zoneI,new scalarField("D",dictAxialExp,energyGroups_));
      nuSigmaEffAxialExpList.set(zoneI,new scalarField("nuSigmaEff",dictAxialExp,energyGroups_));
      sigmaPowAxialExpList.set(zoneI,new scalarField("sigmaPow",dictAxialExp,energyGroups_));
      sigmaDisappAxialExpList.set(zoneI,new scalarField("sigmaDisapp",dictAxialExp,energyGroups_));
      sigmaFromToAxialExpList[zoneI].setSize(legendreMoments_);
      forAll(sigmaFromToAxialExpList[zoneI],momentI)
      {
        sigmaFromToAxialExpList[zoneI].set(momentI,new scalarSquareMatrix(dictAxialExp.lookup("scatteringMatrix"+std::string("P")+Foam::name(momentI))));
      }
    }
    else
    {
      DAxialExpList.set(zoneI,new scalarField(DList_[zoneI]));
      nuSigmaEffAxialExpList.set(zoneI,new scalarField(nuSigmaEffList_[zoneI]));
      sigmaPowAxialExpList.set(zoneI,new scalarField(sigmaPowList_[zoneI]));
      sigmaDisappAxialExpList.set(zoneI,new scalarField(sigmaDisappList_[zoneI]));
      sigmaFromToAxialExpList[zoneI].setSize(legendreMoments_);
      forAll(sigmaFromToAxialExpList[zoneI],momentI)
      {
        sigmaFromToAxialExpList[zoneI].set(momentI,new scalarSquareMatrix(sigmaFromToList_[zoneI][momentI]));
      }
    }
    //fuelTemp
    read = false;
    perturbedEntry = 0;
    forAll(entriesFuelTemp,entryI)
    {
      if(entriesFuelTemp[entryI].keyword() == name)
      {
        read = true;
        perturbedEntry = entryI;
      }
    }
    if(read)
    {
      dictionary& dictFuelTemp = entriesFuelTemp[perturbedEntry].dict();
      DFuelTempList.set(zoneI,new scalarField("D",dictFuelTemp,energyGroups_));
      nuSigmaEffFuelTempList.set(zoneI,new scalarField("nuSigmaEff",dictFuelTemp,energyGroups_));
      sigmaPowFuelTempList.set(zoneI,new scalarField("sigmaPow",dictFuelTemp,energyGroups_));
      sigmaDisappFuelTempList.set(zoneI,new scalarField("sigmaDisapp",dictFuelTemp,energyGroups_));
      sigmaFromToFuelTempList[zoneI].setSize(legendreMoments_);
      forAll(sigmaFromToFuelTempList[zoneI],momentI)
      {
        sigmaFromToFuelTempList[zoneI].set(momentI,new scalarSquareMatrix(dictFuelTemp.lookup("scatteringMatrix"+std::string("P")+Foam::name(momentI))));
      }
    }
    else
    {
      DFuelTempList.set(zoneI,new scalarField(DList_[zoneI]));
      nuSigmaEffFuelTempList.set(zoneI,new scalarField(nuSigmaEffList_[zoneI]));
      sigmaPowFuelTempList.set(zoneI,new scalarField(sigmaPowList_[zoneI]));
      sigmaDisappFuelTempList.set(zoneI,new scalarField(sigmaDisappList_[zoneI]));
      sigmaFromToFuelTempList[zoneI].setSize(legendreMoments_);
      forAll(sigmaFromToFuelTempList[zoneI],momentI)
      {
        sigmaFromToFuelTempList[zoneI].set(momentI,new scalarSquareMatrix(sigmaFromToList_[zoneI][momentI]));
      }
    }
    //RhoCool
    read = false;
    perturbedEntry = 0;
    forAll(entriesRhoCool,entryI)
    {
      if(entriesRhoCool[entryI].keyword() == name)
      {
        read = true;
        perturbedEntry = entryI;
      }
    }
    if(read)
    {
      dictionary& dictRhoCool = entriesRhoCool[perturbedEntry].dict();
      DRhoCoolList.set(zoneI,new scalarField("D",dictRhoCool,energyGroups_));
      nuSigmaEffRhoCoolList.set(zoneI,new scalarField("nuSigmaEff",dictRhoCool,energyGroups_));
      sigmaPowRhoCoolList.set(zoneI,new scalarField("sigmaPow",dictRhoCool,energyGroups_));
      sigmaDisappRhoCoolList.set(zoneI,new scalarField("sigmaDisapp",dictRhoCool,energyGroups_));
      sigmaFromToRhoCoolList[zoneI].setSize(legendreMoments_);
      forAll(sigmaFromToRhoCoolList[zoneI],momentI)
      {
        sigmaFromToRhoCoolList[zoneI].set(momentI,new scalarSquareMatrix(dictRhoCool.lookup("scatteringMatrix"+std::string("P")+Foam::name(momentI))));
      }
    }
    else
    {
      DRhoCoolList.set(zoneI,new scalarField(DList_[zoneI]));
      nuSigmaEffRhoCoolList.set(zoneI,new scalarField(nuSigmaEffList_[zoneI]));
      sigmaPowRhoCoolList.set(zoneI,new scalarField(sigmaPowList_[zoneI]));
      sigmaDisappRhoCoolList.set(zoneI,new scalarField(sigmaDisappList_[zoneI]));
      sigmaFromToRhoCoolList[zoneI].setSize(legendreMoments_);
      forAll(sigmaFromToRhoCoolList[zoneI],momentI)
      {
        sigmaFromToRhoCoolList[zoneI].set(momentI,new scalarSquareMatrix(sigmaFromToList_[zoneI][momentI]));
      }
    }
    //TCool
    read = false;
    perturbedEntry = 0;
    forAll(entriesTCool,entryI)
    {
      if(entriesTCool[entryI].keyword() == name)
      {
        read = true;
        perturbedEntry = entryI;
      }
    }
    if(read)
    {
      dictionary& dictTCool = entriesTCool[perturbedEntry].dict();
      DTCoolList.set(zoneI,new scalarField("D",dictTCool,energyGroups_));
      nuSigmaEffTCoolList.set(zoneI,new scalarField("nuSigmaEff",dictTCool,energyGroups_));
      sigmaPowTCoolList.set(zoneI,new scalarField("sigmaPow",dictTCool,energyGroups_));
      sigmaDisappTCoolList.set(zoneI,new scalarField("sigmaDisapp",dictTCool,energyGroups_));
      sigmaFromToTCoolList[zoneI].setSize(legendreMoments_);
      forAll(sigmaFromToTCoolList[zoneI],momentI)
      {
        sigmaFromToTCoolList[zoneI].set(momentI,new scalarSquareMatrix(dictTCool.lookup("scatteringMatrix"+std::string("P")+Foam::name(momentI))));
      }
    }
    else
    {
      DTCoolList.set(zoneI,new scalarField(DList_[zoneI]));
      nuSigmaEffTCoolList.set(zoneI,new scalarField(nuSigmaEffList_[zoneI]));
      sigmaPowTCoolList.set(zoneI,new scalarField(sigmaPowList_[zoneI]));
      sigmaDisappTCoolList.set(zoneI,new scalarField(sigmaDisappList_[zoneI]));
      sigmaFromToTCoolList[zoneI].setSize(legendreMoments_);
      forAll(sigmaFromToTCoolList[zoneI],momentI)
      {
        sigmaFromToTCoolList[zoneI].set(momentI,new scalarSquareMatrix(sigmaFromToList_[zoneI][momentI]));
      }
    }
    //cladExp
    read = false;
    perturbedEntry = 0;
    forAll(entriesCladExp,entryI)
    {
      if(entriesCladExp[entryI].keyword() == name)
      {
        read = true;
        perturbedEntry = entryI;
      }
    }
    if(read)
    {
      dictionary& dictCladExp = entriesCladExp[perturbedEntry].dict();
      DCladExpList.set(zoneI,new scalarField("D",dictCladExp,energyGroups_));
      nuSigmaEffCladExpList.set(zoneI,new scalarField("nuSigmaEff",dictCladExp,energyGroups_));
      sigmaPowCladExpList.set(zoneI,new scalarField("sigmaPow",dictCladExp,energyGroups_));
      sigmaDisappCladExpList.set(zoneI,new scalarField("sigmaDisapp",dictCladExp,energyGroups_));
      sigmaFromToCladExpList[zoneI].setSize(legendreMoments_);
      forAll(sigmaFromToCladExpList[zoneI],momentI)
      {
        sigmaFromToCladExpList[zoneI].set(momentI,new scalarSquareMatrix(dictCladExp.lookup("scatteringMatrix"+std::string("P")+Foam::name(momentI))));
      }
    }
    else
    {
      DCladExpList.set(zoneI,new scalarField(DList_[zoneI]));
      nuSigmaEffCladExpList.set(zoneI,new scalarField(nuSigmaEffList_[zoneI]));
      sigmaPowCladExpList.set(zoneI,new scalarField(sigmaPowList_[zoneI]));
      sigmaDisappCladExpList.set(zoneI,new scalarField(sigmaDisappList_[zoneI]));
      sigmaFromToCladExpList[zoneI].setSize(legendreMoments_);
      forAll(sigmaFromToCladExpList[zoneI],momentI)
      {
        sigmaFromToCladExpList[zoneI].set(momentI,new scalarSquareMatrix(sigmaFromToList_[zoneI][momentI]));
      }
    }

    //Feedbacks Coefficients
    if(fastNeutrons_)
    {
        dimensionedScalar TrefRatio("",dimensionSet(0,0,0,0,0,0,0),(TfuelRef_/TfuelPerturbed_));
        fuelTempDList_.set(zoneI,new scalarField((DList_[zoneI]-DFuelTempList[zoneI])/log(TrefRatio).value()));
        fuelTempNuSigmaEffList_.set(zoneI,new scalarField((nuSigmaEffList_[zoneI]-nuSigmaEffFuelTempList[zoneI])/log(TrefRatio).value()));
        fuelTempSigmaPowList_.set(zoneI,new scalarField((sigmaPowList_[zoneI]-sigmaPowFuelTempList[zoneI])/log(TrefRatio).value()));
        fuelTempSigmaDisappList_.set(zoneI,new scalarField((sigmaDisappList_[zoneI]-sigmaDisappFuelTempList[zoneI])/log(TrefRatio).value()));
        fuelTempSigmaFromToList_[zoneI].setSize(legendreMoments_);
        forAll(fuelTempSigmaFromToList_[zoneI],momentI)
        {
          fuelTempSigmaFromToList_[zoneI].set(momentI,new scalarSquareMatrix(sigmaFromToList_[zoneI][momentI]));
          for(label energyJ=0; energyJ<energyGroups_ ; energyJ++)
          {
              for(label energyI=0; energyI<energyGroups_ ; energyI++)
              {
                  fuelTempSigmaFromToList_[zoneI][momentI][energyJ][energyI] = (sigmaFromToList_[zoneI][momentI][energyJ][energyI]-sigmaFromToFuelTempList[zoneI][momentI][energyJ][energyI])/log(TrefRatio).value();
              }
          }
        }
    }
    else
    {
        dimensionedScalar TrefDifference("",dimensionSet(0,0,0,0,0,0,0),(Foam::pow(TfuelRef_,0.5)-Foam::pow(TfuelPerturbed_,0.5)));
        fuelTempDList_.set(zoneI,new scalarField((DList_[zoneI]-DFuelTempList[zoneI])/TrefDifference.value()));
        fuelTempNuSigmaEffList_.set(zoneI,new scalarField((nuSigmaEffList_[zoneI]-nuSigmaEffFuelTempList[zoneI])/TrefDifference.value()));
        fuelTempSigmaPowList_.set(zoneI,new scalarField((sigmaPowList_[zoneI]-sigmaPowFuelTempList[zoneI])/TrefDifference.value()));
        fuelTempSigmaDisappList_.set(zoneI,new scalarField((sigmaDisappList_[zoneI]-sigmaDisappFuelTempList[zoneI])/TrefDifference.value()));
        fuelTempSigmaFromToList_[zoneI].setSize(legendreMoments_);
        forAll(fuelTempSigmaFromToList_[zoneI],momentI)
        {
          fuelTempSigmaFromToList_[zoneI].set(momentI,new scalarSquareMatrix(sigmaFromToList_[zoneI][momentI]));
          for(label energyJ=0; energyJ<energyGroups_ ; energyJ++)
          {
              for(label energyI=0; energyI<energyGroups_ ; energyI++)
              {
                  fuelTempSigmaFromToList_[zoneI][momentI][energyJ][energyI] = (sigmaFromToList_[zoneI][momentI][energyJ][energyI]-sigmaFromToFuelTempList[zoneI][momentI][energyJ][energyI])/TrefDifference.value();
              }
          }
        }
    }
    axialExpDList_.set(zoneI,new scalarField(-(DList_[zoneI]-DAxialExpList[zoneI])/AxExp_));
    axialExpNuSigmaEffList_.set(zoneI,new scalarField(-(nuSigmaEffList_[zoneI]-nuSigmaEffAxialExpList[zoneI])/AxExp_));
    axialExpSigmaPowList_.set(zoneI,new scalarField(-(sigmaPowList_[zoneI]-sigmaPowAxialExpList[zoneI])/AxExp_));
    axialExpSigmaDisappList_.set(zoneI,new scalarField(-(sigmaDisappList_[zoneI]-sigmaDisappAxialExpList[zoneI])/AxExp_));
    axialExpSigmaFromToList_[zoneI].setSize(legendreMoments_);
    forAll(axialExpSigmaFromToList_[zoneI],momentI)
    {
      axialExpSigmaFromToList_[zoneI].set(momentI,new scalarSquareMatrix(sigmaFromToList_[zoneI][momentI]));
      for(label energyJ=0; energyJ<energyGroups_ ; energyJ++)
      {
          for(label energyI=0; energyI<energyGroups_ ; energyI++)
          {
              axialExpSigmaFromToList_[zoneI][momentI][energyJ][energyI] = (-sigmaFromToList_[zoneI][momentI][energyJ][energyI]+sigmaFromToAxialExpList[zoneI][momentI][energyJ][energyI])/AxExp_;
          }
      }
    }    

    radialExpDList_.set(zoneI,new scalarField(-(DList_[zoneI]-DRadialExpList[zoneI])/RadExp_));
    radialExpNuSigmaEffList_.set(zoneI,new scalarField(-(nuSigmaEffList_[zoneI]-nuSigmaEffRadialExpList[zoneI])/RadExp_));
    radialExpSigmaPowList_.set(zoneI,new scalarField(-(sigmaPowList_[zoneI]-sigmaPowRadialExpList[zoneI])/RadExp_));
    radialExpSigmaDisappList_.set(zoneI,new scalarField(-(sigmaDisappList_[zoneI]-sigmaDisappRadialExpList[zoneI])/RadExp_));
    radialExpSigmaFromToList_[zoneI].setSize(legendreMoments_);
    forAll(radialExpSigmaFromToList_[zoneI],momentI)
    {
      radialExpSigmaFromToList_[zoneI].set(momentI,new scalarSquareMatrix(sigmaFromToList_[zoneI][momentI]));
      for(label energyJ=0; energyJ<energyGroups_ ; energyJ++)
      {
          for(label energyI=0; energyI<energyGroups_ ; energyI++)
          {
              radialExpSigmaFromToList_[zoneI][momentI][energyJ][energyI] = (-sigmaFromToList_[zoneI][momentI][energyJ][energyI]+sigmaFromToRadialExpList[zoneI][momentI][energyJ][energyI])/RadExp_;
          }
      }
    }  

    scalar deltaRho(rhoCoolRef_ - rhoCoolPerturbed_);
    rhoCoolDList_.set(zoneI,new scalarField((DList_[zoneI]-DRhoCoolList[zoneI])/deltaRho));
    rhoCoolNuSigmaEffList_.set(zoneI,new scalarField((nuSigmaEffList_[zoneI]-nuSigmaEffRhoCoolList[zoneI])/deltaRho));
    rhoCoolSigmaPowList_.set(zoneI,new scalarField((sigmaPowList_[zoneI]-sigmaPowRhoCoolList[zoneI])/deltaRho));
    rhoCoolSigmaDisappList_.set(zoneI,new scalarField((sigmaDisappList_[zoneI]-sigmaDisappRhoCoolList[zoneI])/deltaRho));
    rhoCoolSigmaFromToList_[zoneI].setSize(legendreMoments_);
    forAll(rhoCoolSigmaFromToList_[zoneI],momentI)
    {
      rhoCoolSigmaFromToList_[zoneI].set(momentI,new scalarSquareMatrix(sigmaFromToList_[zoneI][momentI]));
      for(label energyJ=0; energyJ<energyGroups_ ; energyJ++)
      {
          for(label energyI=0; energyI<energyGroups_ ; energyI++)
          {
              rhoCoolSigmaFromToList_[zoneI][momentI][energyJ][energyI] = (sigmaFromToList_[zoneI][momentI][energyJ][energyI]-sigmaFromToRhoCoolList[zoneI][momentI][energyJ][energyI])/deltaRho;
          }
      }
    }  

    scalar deltaTCool(TCoolRef_ - TCoolPerturbed_);
    TCoolDList_.set(zoneI,new scalarField((DList_[zoneI]-DTCoolList[zoneI])/deltaTCool));
    TCoolNuSigmaEffList_.set(zoneI,new scalarField((nuSigmaEffList_[zoneI]-nuSigmaEffTCoolList[zoneI])/deltaTCool));
    TCoolSigmaPowList_.set(zoneI,new scalarField((sigmaPowList_[zoneI]-sigmaPowTCoolList[zoneI])/deltaTCool));
    TCoolSigmaDisappList_.set(zoneI,new scalarField((sigmaDisappList_[zoneI]-sigmaDisappTCoolList[zoneI])/deltaTCool));
    TCoolSigmaFromToList_[zoneI].setSize(legendreMoments_);
    forAll(TCoolSigmaFromToList_[zoneI],momentI)
    {
      TCoolSigmaFromToList_[zoneI].set(momentI,new scalarSquareMatrix(sigmaFromToList_[zoneI][momentI]));
      for(label energyJ=0; energyJ<energyGroups_ ; energyJ++)
      {
          for(label energyI=0; energyI<energyGroups_ ; energyI++)
          {
              TCoolSigmaFromToList_[zoneI][momentI][energyJ][energyI] = (sigmaFromToList_[zoneI][momentI][energyJ][energyI]-sigmaFromToTCoolList[zoneI][momentI][energyJ][energyI])/deltaTCool;
          }
      }
    }  

    scalar deltaTClad(TcladRef_ - TcladPerturbed_);
    cladExpDList_.set(zoneI,new scalarField((DList_[zoneI]-DCladExpList[zoneI])/deltaTClad));
    cladExpNuSigmaEffList_.set(zoneI,new scalarField((nuSigmaEffList_[zoneI]-nuSigmaEffCladExpList[zoneI])/deltaTClad));
    cladExpSigmaPowList_.set(zoneI,new scalarField((sigmaPowList_[zoneI]-sigmaPowCladExpList[zoneI])/deltaTClad));
    cladExpSigmaDisappList_.set(zoneI,new scalarField((sigmaDisappList_[zoneI]-sigmaDisappCladExpList[zoneI])/deltaTClad));
    cladExpSigmaFromToList_[zoneI].setSize(legendreMoments_);
    forAll(cladExpSigmaFromToList_[zoneI],momentI)
    {
      cladExpSigmaFromToList_[zoneI].set(momentI,new scalarSquareMatrix(sigmaFromToList_[zoneI][momentI]));
      for(label energyJ=0; energyJ<energyGroups_ ; energyJ++)
      {
          for(label energyI=0; energyI<energyGroups_ ; energyI++)
          {
              cladExpSigmaFromToList_[zoneI][momentI][energyJ][energyI] = (sigmaFromToList_[zoneI][momentI][energyJ][energyI]-sigmaFromToCladExpList[zoneI][momentI][energyJ][energyI])/deltaTClad;
          }
      }
    }
}

forAll(CRentries_,CRzoneI)
{

    dictionary& CRdict = CRentries_[CRzoneI].dict();
    CRstart_.set(CRzoneI,new scalar(CRdict.lookupOrDefault("startTime",0.0)));
    CRfinish_.set(CRzoneI,new scalar(CRdict.lookupOrDefault("endTime",0.0)));
    CRspeed_.set(CRzoneI,new scalar(CRdict.lookupOrDefault("speed",0.0)));
    CRFollowerName_.set(CRzoneI,new word(CRdict.lookup("followerName")));
    initialDistanceFromMeshCR_.set(CRzoneI,new scalar(CRdict.lookupOrDefault("initialDistanceFromMeshCR",0.0)));

    label FOLLOWERzoneId = mesh.cellZones().findZoneID(CRFollowerName_[CRzoneI]);

    if(FOLLOWERzoneId == -1)
    {
      Info << "Control rod with follower: " << CRFollowerName_(CRzoneI) << " does not exists" << endl;
    }

    scalar cellHeight = -100.0 ;
    scalar maxFollowerHeight = cellHeight ;

    volScalarField centersVSF(mesh.C() & vector(0,0,1));

    Field<scalar> centersSF(centersVSF.internalField());

    forAll(mesh.cellZones()[FOLLOWERzoneId], cellIlocal)
    {

        label cellIglobal = mesh.cellZones()[FOLLOWERzoneId][cellIlocal];

        maxFollowerHeight = max(cellHeight,centersSF[cellIglobal]) ;

        cellHeight = maxFollowerHeight;

    }
    reduce(maxFollowerHeight, maxOp<scalar>());
    CRinitialPosition_.set(CRzoneI, new scalar(maxFollowerHeight-initialDistanceFromMeshCR_[CRzoneI]));
    if(mesh.time().value()<CRfinish_[CRzoneI])
    {
        CRinitialPosition_[CRzoneI] = (CRinitialPosition_[CRzoneI]-CRspeed_[CRzoneI]*(max(mesh.time().value(),CRstart_[CRzoneI])-CRstart_[CRzoneI]));
    }
    else
    {
        CRinitialPosition_[CRzoneI] = (CRinitialPosition_[CRzoneI]-CRspeed_[CRzoneI]*(CRfinish_[CRzoneI]-CRstart_[CRzoneI]));
    }
    CRposition_.set(CRzoneI, new scalar(CRinitialPosition_[CRzoneI]));
    Info << "Control rod " << CRentries_[CRzoneI].keyword() << " in position: " << CRposition_[CRzoneI] << endl;

}



forAll(nuSigmaEff_,energyI)
{
  forAll(nuSigmaEff_,energyJ)
  {
    label exchangeYesNo = 0;
    forAll(entries_,zoneI)
    {
      if(sigmaFromToList_[zoneI][0][energyJ][energyI]>5.0e-6) {exchangeYesNo = 1;}
    }
    sigmaFromToYesNo_[energyJ][energyI] = exchangeYesNo;
  }
}


Info << "Finished reading nuclear data " << endl;




