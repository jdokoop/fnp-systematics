# fnp-systematics
Calculation of systematic uncertainties on the measured fraction of non-photonic electrons

The calculation of the fraction of non-photonic electrons, styled FNP, in p+p collisions depends quite sensitively on how we choose to model the production of etas, pizeros, and direct photons. We decided to use existing published measurements of the invariant yield of electrons from the above particles, and fit the measurements with a modified Hagedorn functional form. The particles for our simulations are sampled from these fits. 

In order to estimate the systematic uncertainties on FNP arising from our simulation of photonic electrons, we vary the fit to the published data to take into account the following sources of error
* The particle yield at low pT, where we rely on extrapolating the fits to published data
* Uncertainties on the shape of the published invariant yield, as quantified by their own systematics
* The pT-dependent ratio of etas to pions
