import os
import h5py



def read_hdf5(fname, group, slices=None):
    """
    Read data from hdf5 file from a specific group.
    Parameters
    ----------
    fname : str
        String defining the path or file name.
    group : str
        Path to the group inside hdf5 file where data is located.
    slices : list of slice objects, optional
        (slice(start_1, end_1, step_1), ... , slice(start_N, end_N, step_N))
        defines slicing parameters for each axis of the data matrix.
    Returns
    -------
    ndarray
        Data.
    """
    fname = os.path.abspath(fname)
    f = h5py.File(fname, "r")
    arr = f[group]
    ndim = len(arr.shape)
    if slices is None:
        slices = [slice(None)] * ndim
    elif len(slices) < ndim:
        slices = list(slices) + [slice(None)] * (ndim - len(slices))
    arr = arr[slices]
    f.close()
    return arr


def read_aps_2bm(fname, slices, white_slices, dark_slices, path2white=None, path2dark=None):
    """
    Read APS 2-BM standard data format.
    Parameters
    ----------
    fname : str
        Path to hdf5 file.
    slices : list of slice objects, optional
        Specifies the slices to read of each dimension (projections, z, x) 
        of the data from a list of slice objects.
    white_slices : list of slice objects, optional
        Specifies the slices to read of each dimension (projections, z, x) 
        of the white from a list of slice objects.
    dark_slices : list of slice objects, optional
        Specifies the slices to read of each dimension (projections, z, x) 
        of the dark from a list of slice objects.
    path2white : str, optional
        Specifies the path to the file containing whites. 
        If None, whites from datafile will be used.
    path2dark : str, optional
        Specifies the path to the file containing darks. 
        If None, darks from datafile will be used.
    Returns
    -------
    ndarray
        3D tomographic data.
    ndarray
        3d flat field data.
    ndarray
        3D dark field data.
    """
    if path2white is None:
        path2white = fname
    if path2dark is None:
        path2dark = fname
    tomo_grp = os.path.join('exchange', 'data')
    flat_grp = os.path.join('exchange', 'data_white')
    dark_grp = os.path.join('exchange', 'data_dark')
    tomo = read_hdf5(fname, tomo_grp, slices=slices)
    flat = read_hdf5(path2white, flat_grp, slices=white_slices)
    dark = read_hdf5(path2white, dark_grp, slices=dark_slices)
    return tomo, flat, dark