
adjointPrec_beta_ *= 0.0;
forAll(prec_,precI)
{
    adjointPrec_beta_ += Beta[precI]*prec_[precI];
}

