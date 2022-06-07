# riblets-linear-stability-model

Part of Vu's IIB Project - Capturing the onset of Kelvin-Helmholtz rollers over riblets.

Steps to analyse the stability of a given riblet geometry:

1. Open the "Stokes flow" folder.
2. Use stokes_shear.m to obtain porous boundary position and shear-driven permeability coefficients.
3. Use stokes_pressure.m to obtain pressure-driven coefficients - change ly to calculate these coefficients at different porous boundary positions.

4. Return to the main folder.
5. Run launchscript.m using these coefficients to obtain data file.
6. Plot results using plot_something.m scripts.