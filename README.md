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
from functions import remove_similar_elements, healpix_direction_vectors, binary_MFP_map
```





