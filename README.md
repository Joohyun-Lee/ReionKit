# ReionKit: An Analysis Tool for Reionization Simulation

## Prerequisites
1. Numpy
2. Scipy
3. (Mpi4py): preferred
4. Healpy

## Modules

### Measuring mean free path (MFP) bubble sizes
1) Initialize a simulation object
2) Run analysis on a slice of the simulated ionization field (binary; np.bool dtype)

---
Example code
```
import numpy as np # type: ignore
import healpy as hp # type: ignore
import scipy # type: ignore
import scipy.ndimage # type: ignore
import time # type: ignore
from mpi4py import MPI # type: ignore

import simulation
from functions import binary_MFP_map

# change cosmological parameters!
h=0.677699966430664
omega_b=0.450000017881393E-01
omega_m=0.307114988565445E+00

sim = simulation.CosmoSim(boxsize=64, h=h, Omega_m=omega_m, Omega_lambda=1 - omega_m, Omega_b=omega_b)

path = 'whatever path of the ionization field (.npy with np.bool dtype)'
savepath = 'whatever path of a directory you want to save the outputs'
offset = i_z # layer number you want to calculate MFP bubble sizes

binary_MFP_map(CoDa, path, savepath,
               n_cell = 1024, num_sightlines = 192, max_reach = 2.5, resolution = 2,
               mpi_start_i = 0, mpi_rank = offset)
```

