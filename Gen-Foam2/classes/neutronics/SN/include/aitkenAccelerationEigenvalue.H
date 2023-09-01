
bool aitkenYesNo = true;

if(aitkenIterCounter_==0) // only happening at first timeStep 
{

    keff0_=keff_;

}  
if(aitkenIterCounter_==1)
{

    keff1_=keff_;

}
if(aitkenIterCounter_==2)
{
    aitkenYesNo = true;

    keff2_=keff_;

    if(
    mag(keff2_-keff1_)/(max(mag(keff1_-keff0_),SMALL)) > 0.95 
    || (mag(keff2_-2*keff1_+keff0_)<SMALL)
    ) aitkenYesNo=false;


    if(aitkenYesNo==true)
    {
        Info << "Intra-step eigenvalue Aitken acceleration no: " << aitkenIterNo_ << endl;

        //keff2_=keff_;
        /*
        scalar aitkenScaling = 
        (
            ( 
                keff0_
                - (keff1_-keff0_)
                * (keff1_-keff0_)
                / (keff2_-2*keff1_+keff0_)
            )/keff2_
        );
        keff_ *= aitkenScaling ;
        */
        keff_ = keff0_ - pow((keff1_-keff0_),2) / (keff2_-2*keff1_+keff0_);

        keff0_=keff_;

        aitkenIterNo_++ ;
        aitkenIterCounter_ = 0;
        Info<< "Done "  << endl;
    }
    else
    {
        Info << "Aitken not performed because of diverging series" << endl;

        keff0_ = keff1_;
        keff1_ = keff2_;

        aitkenIterCounter_-- ;
    }

}
aitkenIterCounter_++ ;

