import numpy as np # type: ignore
import healpy as hp # type: ignore
import scipy # type: ignore
import scipy.ndimage # type: ignore
import time # type: ignore
from mpi4py import MPI # type: ignore

import simulation
from functions import remove_similar_elements, healpix_direction_vectors, binary_MFP_map


# CC=mpicc MPICC=mpicc pip install mpi4py --no-binary=mpi4py


h=0.677699966430664
omega_b=0.450000017881393E-01
omega_m=0.307114988565445E+00



# initialize simulation
CoDa = simulation.CosmoSim(boxsize=64, h=h, Omega_m=omega_m, Omega_lambda=1 - omega_m, Omega_b=omega_b)



mpi_comm = MPI.COMM_WORLD
mpi_size = mpi_comm.Get_size()
mpi_rank = mpi_comm.Get_rank()

print('mpi (size, rank): ', mpi_size, mpi_rank)



snapnum=90


path = './'
path = path + f'output_{snapnum :06d}_ion_reduced.npy'
savepath = f'./output_{snapnum :06d}'


# [0, 48]
# [96, 144]
# [192, 240] 
# [288, 336]
# [384, 432]
# [480, 528]
# [576, 624]
# [672, 720]
# [768, 816]
# [864, 912, 960]
# [1008-1024]
for offset in [96, 144]:
    binary_MFP_map(CoDa, path, savepath,
                   n_cell = 1024, num_sightlines = 192, max_reach = 2.5, resolution = 2,
                   mpi_start_i = 0, mpi_rank = offset + mpi_rank)

