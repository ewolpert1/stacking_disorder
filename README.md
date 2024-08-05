This repository contains code written in fortran which simulates the stacking disorder in graphite. `graphite_stacking_ising.f90` simulates the stacking disorder of graphite using the Ising model:

$$    H=J\sum_{\langle ij\rangle}S_iS_j$$

where $`S_{i/j}=\pm1`$ and represents the slip directions between neighbouring layers such that when $`S_i=1(-1)`$ there is a slip in the "positive(negative)" direction *i.e.* AB(AC). `graphite_stacking_continuous.f90` models the stacking disorder using the following equation:

$$    H=J_1\sum_{\langle ij\rangle}\cos(2\pi \Delta x_{ij} )\cos(2\pi \Delta y_{ij})\cos(2\pi (\Delta x_{ij}-\Delta y_{ij}))\\
    -J_2\sum_{\langle\langle ij\rangle\rangle}\cos(2\pi \Delta x_{ij} )\cos(2\pi \Delta y_{ij})\cos(2\pi (\Delta x_{ij}-\Delta y_{ij})) $$


where ($`x_i,y_i`$) and ($`x_j,y_j`$) represent the fractional coordinates of neighbouring layers (or next-nearest layers), $`i`$ and $`j`$, of graphite and $`\Delta x_{ij}=x_i-x_j`$, $`\Delta y_{ij}=y_i-y_j`$. For positive values of $`J_1`$ and $`J_2`$, the Hamiltonian is minimised for stacking arrangements of graphene sheets in an ABAB fashion.

The parameters to used to generate the results from Fig 4a and 4b are give in the `parameters_Ising.txt` and `parameters_continuous.txt` respectively. These files are read in when running the two different fortran codes.

Running this code outputs cif files of graphite at different temperatures of the simulations. For the Ising/continuous model the cifs are labelled `phase_Ising_T.cif/phase_continuous_T.cif` where T is the temperature.

The size of the unit cell, interlayer spacing, interaction parameters, and number of layers in the model can all be changed by editing either the `parameters_Ising.txt` of `parameters_continuous.txt` files.

