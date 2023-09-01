//- Update and limit the static pressure
p_ = max(p_rgh_ + rho_*gh_, pMin_);

//- Limit p_rgh
p_rgh_ = p_ - rho_*gh_;

//- Update densities from change in p_rgh
rho1 += fluid1_.thermo().psi()*(p_rgh_ - p_rgh_0);
rho2 += fluid2_.thermo().psi()*(p_rgh_ - p_rgh_0);

//- Correct p_rgh for consistency with p and the updated densities
rho_ = (rho1*alpha1 + rho2*alpha2)/movingAlpha_;
p_rgh_ = p_ - rho_*gh_;
p_rgh_.correctBoundaryConditions();
