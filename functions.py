import numpy as np  # type: ignore
import healpy as hp  # type: ignore
import time  # type: ignore
import scipy  # type: ignore
import simulation  # type: ignore


### remove elements with too little difference
### from an array of points that intersect with grid
def remove_similar_elements(arr, threshold):
    if len(arr) == 0:
        return arr

    # Initialize the result array with the first element
    result = [arr[0]]

    # Iterate through the sorted array and add elements
    # that are sufficiently different
    for i in range(1, len(arr)):
        if np.abs(arr[i] - result[-1]) >= threshold:
            result.append(arr[i])

    return np.array(result)


### HEALPix sightlines
def healpix_direction_vectors(num_sightlines=192):
    nside = hp.npix2nside(num_sightlines)

    pixel_indices = np.arange(hp.nside2npix(nside))

    theta, phi = hp.pix2ang(nside, pixel_indices)

    x = np.sin(theta) * np.cos(phi)
    y = np.sin(theta) * np.sin(phi)
    z = np.cos(theta)

    direction_vectors = np.vstack((x, y, z)).T

    return direction_vectors


### measure MFP bubble size
def binary_MFP_map(
    cosmosim: simulation.CosmoSim,
    path: str,
    savepath: str,
    n_cell: int = 1024,
    num_sightlines: int = 192,
    max_reach: float = 2.5,
    resolution: int = 2,
    mpi_start_i: int = 0,
    mpi_rank: int = 0,
):
    """
    Saves n_cell MFP bubble size arrays
    with n_cell^2 shape using mpi.
    Each MFP bubble size array is one layer of n_cell^3 array,
    using binary_map[:,:,mpi_start_i + mpi_rank].

    Parameters
    ----------
    - cosmosim: Parent cosmological simulation object
    - path: Path of binary_map numpy array
            with n_cell^3 shape (str)
    - savepath: Path to save n_cell MFP_map arrays
                with n_cell^2 shape (str)
    - n_cell: Number of catersian cells (int)
    - num_sightlines: Number of sightlines used for Healpix sightlines
      (int; 12, 48, 192, 768)
    - max_reach: Maximum reach of sightlines in n_cell unit (float)
    - resolution: How finely interpolate binary_map
                  on sightline; delta_x = cell_width/resolution (int)
    - mpi_start_i: Layer number to start mpi process
      (int; 0 - (n_cell-1))
    - mpi_rank: mpi_rank of this CPU
      (int; recommend mpi4py.MPI.COMM_WORLD.Get_rank())

    Outputs
    -------
    - Return: Nothing
    - Print: Progress with execution time
    - Save: MFP_{z :d}.npy and sigma_{z :d}.npy in savepath

    Example
    -------
        path = f'./output_000090_ion_reduced.npy'
        savepath = f'./output_000090'

        binary_MFP_map(cosmosim,
                       path,
                       savepath,
                       n_cell=1024,
                       num_sightlines=192,
                       max_reach=2.5,
                       resolution=4,
                       mpi_start_i=0,
                       mpi_rank = mpi4py.MPI.COMM_WORLD.Get_rank())
    """

    ### define sightline vectors
    direction_vectors = healpix_direction_vectors(num_sightlines)

    ### load data
    binary_map = np.ascontiguousarray(np.load(path))

    ### initialize MFP bubble size array 
    ### and standard deviation array
    MFP = np.ascontiguousarray(np.zeros((n_cell, n_cell), dtype=np.float32))
    sigma = np.ascontiguousarray(np.zeros((n_cell, n_cell), dtype=np.float32))

    start_time = time.time()

    ### measure MFP_size in each cell
    z = mpi_start_i + mpi_rank

    for x in range(n_cell):
        for y in range(n_cell):
            if binary_map[x, y, z]:
                MFP_t = np.zeros(num_sightlines, np.float64)

                for i in range(num_sightlines):
                    d1, d2, d3 = direction_vectors[i, :]

                    alpha = np.linspace(
                        0, max_reach * n_cell, int(max_reach * n_cell * resolution)
                    )  # skewer resolution = cell_width / resolution

                    xyz = np.c_[x + alpha * d1, y + alpha * d2, z + alpha * d3]

                    skewer = scipy.ndimage.map_coordinates(
                        binary_map,
                        xyz.T,
                        output=bool,
                        order=0,
                        mode="grid-wrap",
                        cval=0.0,
                        prefilter=False,
                    )

                    if np.sum(~skewer):
                        MFP_t[i] = alpha[~skewer][0]
                    else:
                        MFP_t[i] = n_cell * max_reach
                        print(x, y, z, i, 'hit max reach')

                MFP[x, y] = np.mean(MFP_t) * cosmosim.boxsize / n_cell / cosmosim.h
                sigma[x, y] = np.std(MFP_t) * cosmosim.boxsize / n_cell / cosmosim.h

    # save MFP bubble size in cMpc unit
    np.save(savepath + f"/MFP_{z :d}.npy", MFP)
    np.save(savepath + f"/sigma_{z :d}.npy", sigma)

    print(z, "done, time:", time.time() - start_time)
