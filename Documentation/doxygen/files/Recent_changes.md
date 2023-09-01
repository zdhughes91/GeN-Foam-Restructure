# Recent changes in the case folder {#CHANGES}

* On September 18th 2022, the reactorState dictionary has been moved from constant to dictionary. We are aware that changes in the input deck should be minimized. However, we find that having reactorState under the uniform folder is much more logical and that it greatly simlifies the use of GeN-Foam. In particular, this allows to have a separate reactorState dictionary for every time step, which in turns allows for a natural restart of time-dependent simulations, notably for the point-kinetics solver.