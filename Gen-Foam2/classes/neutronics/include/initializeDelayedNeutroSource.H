delayedNeutroSource_ *= 0.0;
forAll(prec_, precI)
{
    delayedNeutroSource_ += lambda[precI]*prec_[precI];
}

