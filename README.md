This repository contains code written in fortran which uses Monte Carlo simulations to investigate the stacking disorder in graphite with the continuous and Ising models. The parameters of the models were derived from PBE+MBD calculations along with three stacking-fault energies from the literature: 0.05, 0.16, and 0.21 $\mathrm{meV}\\mathrm{atom}^{-1}$ at the interface. `graphite_stacking_Ising.f90` simulates the stacking disorder of graphite using the Ising model:

$$    H=J\sum_{\langle ij\rangle}S_iS_j$$

where $`S_{i/j}=\pm1`$ and represents the slip directions between neighbouring layers such that when $`S_i=1(-1)`$ there is a slip in the "positive(negative)" direction *i.e.* AB(AC). `graphite_stacking_continuous.f90` models the stacking disorder using the following equation:

$$    H=J_1\sum_{\langle ij\rangle}\cos(\pi \Delta x_{ij} )\cos(\pi \Delta y_{ij})\cos(\pi (\Delta x_{ij}-\Delta y_{ij}))\\
    -J_2\sum_{\langle\langle ij\rangle\rangle}\cos(\pi \Delta x_{ij} )\cos(\pi \Delta y_{ij})\cos(\pi (\Delta x_{ij}-\Delta y_{ij})) $$


where ($`x_i,y_i`$) and ($`x_j,y_j`$) represent the fractional coordinates of neighbouring layers (or next-nearest layers), $`i`$ and $`j`$, of graphite and $`\Delta x_{ij}=x_i-x_j`$, $`\Delta y_{ij}=y_i-y_j`$. For positive values of $`J_1`$ and $`J_2`$, the Hamiltonian is minimised for stacking arrangements of graphene sheets in an ABAB fashion.

Parameter files are provided in `fault_0.05`, `fault_0.16` and `fault_0.21`. Each directory contains `parameters_continuous.txt` and `parameters_Ising.txt`. The coupling inputs are the per-layer-carbon coefficients divided by Boltzmann’s constant, \$J_\mathrm{code}=j/k_B$.

The calculations in the paper use 1000 layers, 10 million attempted moves per temperature and 987 temperature stages from 1000 K to approximately 0.05 K with three seeds: 1729, 2718, and 3141. 

The size of the unit cell, interlayer spacing, interaction parameters, and number of layers in the model can all be changed by editing the `parameters_Ising.txt` or `parameters_continuous.txt` files.

## Running the simulations

Compile the two models using:

```bash
gfortran -O3 -ffree-line-length-none -fcheck=bounds \
  -o continuous graphite_stacking_continuous.f90

gfortran -O3 -ffree-line-length-none -fcheck=bounds \
  -o ising graphite_stacking_Ising.f90
```

 Each simulation should run in a separate directory containing the appropriate parameter file, named `parameters_continuous.txt` or `parameters_Ising.txt`.

 The command-line arguments are:

```text
./continuous SEED NMOVES NSTAGES COOLING_RATIO [RESTART_FILE [START_T [FLOOR_T]]]
./ising      SEED NMOVES NSTAGES COOLING_RATIO [RESTART_FILE [START_T [FLOOR_T]]]
```

- `SEED`: integer used to initialise the random-number generator.
- `NMOVES`: number of attempted moves at each temperature.
- `NSTAGES`: number of temperature stages to perform.
- `COOLING_RATIO`: factor by which the temperature is multiplied after each stage.
- `RESTART_FILE`: optional binary file containing one saved configuration. An empty argument (`""`) starts from a random configuration.
- `START_T`: starting model temperature; defaults to 1000.
- `FLOOR_T`: minimum model temperature; defaults to zero.

The seed is required. If omitted, `NMOVES`, `NSTAGES` and `COOLING_RATIO` default to 1000000, 135 and 0.95, respectively. These defaults differ from the revised manuscript protocol.


## Outputs

At each temperature stage, the programs write:

- A CIF containing the current stacking configuration, named `phase_continuous_T.cif` or `phase_Ising_T.cif`, where `T` is the sampled model temperature.
- A row in `stats.csv` recording the stage, temperature, couplings, accumulated energy, independently recalculated energy, energy drift and acceptance fraction.
- The current configuration in `states.bin`, a raw binary trajectory containing one configuration per temperature stage.

The continuous-model trajectory stores the internal layer translation coordinates; the Ising trajectory stores the slip spins. The programs also write `initial_state.bin` and `initial_energy_K.txt` for the starting configuration.


