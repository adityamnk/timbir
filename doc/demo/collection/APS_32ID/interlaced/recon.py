from __future__ import print_function
import os
import tomopy
import numpy as np
import sys
import reader
from TIMBIR_angles import gen_theta as _gen_theta_timbir
from TIMBIR_angles import calc_dropped_angles


# # Program settings
# MPI
_usempi = False
_nprocs = 1
_rank = 0
# multiprocessing
# only _nchunk is used in starting processes
# if _nchunk is None, _nchunk = (ndim - 1) // _ncore + 1
# where ndim = arr.shape[axis]
_ncore = None
_nchunk = None
# # Number of z slices in a single reconstruction run
# _z_recon_size = None


def set_mpi(usempi=True):
    global _nprocs
    global _rank
    if usempi:
        _usempi = True
        from mpi4py import MPI
        comm = MPI.COMM_WORLD
        _nprocs = comm.Get_size()
        _rank = comm.Get_rank()
    else:
        _usempi = False
        _nprocs = 1
        _rank = 0


def set_mp(ncore=None, nchunk=None):
    global _ncore, _nchunk
    _ncore = ncore
    _nchunk = nchunk


# def set_z_recon_size(z_recon_size=None):
#     global _z_recon_size
#     _z_recon_size = z_recon_size


def show_settings():
    print("rank =", _rank)
    print("nprocs =", _nprocs)
    print("ncore =", _ncore)
    print("nchunk =", _nchunk)
    # print("z_recon_size =", _z_recon_size)


def gen_theta(n, period=np.pi):
    return np.arange(0, period, period/n)


def gen_theta_timbir(K, N_theta, SlewSpeed=0, MinAcqTime=0, TotalNumCycles=1):
    thetas = _gen_theta_timbir(K, N_theta, TotalNumCycles)
    return thetas[~calc_dropped_angles(thetas, SlewSpeed*np.pi/180*MinAcqTime, verbose=True, TotalNumCycles=TotalNumCycles)]


def recon(io_paras, data_paras, rot_center=None, normalize=True, stripe_removal=10, phase_retrieval=False, 
            opt_center=False, diag_center=False, output="tiff"):
    # Input and output
    datafile = io_paras.get('datafile')
    path2white = io_paras.get('path2white', datafile)
    path2dark = io_paras.get('path2dark', path2white)
    out_dir = io_paras.get('out_dir')
    diag_cent_dir = io_paras.get('diag_cent_dir', out_dir+"/center_diagnose/")
    recon_dir = io_paras.get('recon_dir', out_dir+"/recon/")
    out_prefix = io_paras.get('out_prefix', "recon_")

    # Parameters of dataset
    NumCycles = data_paras.get('NumCycles', 1) # Number of cycles used for recon
    ProjPerCycle = data_paras.get('ProjPerCycle') # Number of projections per cycle, N_theta
    cycle_offset = data_paras.get('cycle_offset', 0) # Offset in output cycle number
    proj_start = data_paras.get('proj_start', 0) # Starting projection of reconstruction 
    proj_step = data_paras.get('proj_step')
    z_start = data_paras.get('z_start', 0)
    z_end = data_paras.get('z_end', z_start+1)
    z_step = data_paras.get('z_step')
    x_start = data_paras.get('x_start')
    x_end = data_paras.get('x_end', x_start+1)
    x_step = data_paras.get('x_step')
    white_start = data_paras.get('white_start')
    white_end = data_paras.get('white_end')
    dark_start = data_paras.get('dark_start')
    dark_end = data_paras.get('dark_end')

    rot_center_copy = rot_center

    for cycle in xrange(NumCycles):
        # Set start and end of each cycle
        projections_start = cycle * ProjPerCycle + proj_start
        projections_end = projections_start + ProjPerCycle
        slice1 = slice(projections_start, projections_end, proj_step)
        slice2 = slice(z_start, z_end, z_step)
        slice3 = slice(x_start, x_end, x_step)
        slices = (slice1, slice2, slice3)
        white_slices = (slice(white_start, white_end), slice2, slice3)
        dark_slices = (slice(dark_start, dark_end), slice2, slice3)
        print("Running cycle #%s (projs %s to %s)" 
            % (cycle, projections_start, projections_end))
        
        # Read HDF5 file.
        print("Reading datafile %s..." % datafile, end="")
        sys.stdout.flush()
        data, white, dark = reader.read_aps_2bm(datafile, slices, white_slices, dark_slices, 
                                        path2white=path2white, path2dark=path2dark)
        theta = gen_theta(data.shape[0])
        print("Done!")
        print("Data shape = %s;\nwhite shape = %s;\ndark shape = %s." 
            % (data.shape, white.shape, dark.shape))
        
        ## Normalize dataset using data_white and data_dark
        if normalize:
            print("Normalizing data ...")
            # white = white.mean(axis=0).reshape(-1, *data.shape[1:])
            # dark = dark.mean(axis=0).reshape(-1, *data.shape[1:])
            # data = (data - dark) / (white - dark)
            data = tomopy.normalize(data, white, dark, cutoff=None, ncore=_ncore, nchunk=None)[...]
    
        ## Remove stripes caused by dead pixels in the detector
        if stripe_removal:
            print("Removing stripes ...")
            data = tomopy.remove_stripe_fw(data, level=stripe_removal, wname='db5', sigma=2,
                                    pad=True, ncore=_ncore, nchunk=None)
            # data = tomopy.remove_stripe_ti(data, nblock=0, alpha=1.5, 
            #                                 ncore=None, nchunk=None)

#        # Show preprocessed projection
#        plt.figure("%s-prep" % projections_start)
#        plt.imshow(d.data[0,:,:], cmap=cm.Greys_r)
#        plt.savefig(out_dir+"/preprocess/%s-prep.jpg" 
#                    % projections_start)
#        # plt.show()
#        continue

        ## Phase retrieval
        if phase_retrieval:
            print("Retrieving phase ...")
            data = tomopy.retrieve_phase(data,
                        pixel_size=1e-4, dist=50, energy=20,
                        alpha=1e-3, pad=True, ncore=_ncore, nchunk=None)
        
        ## Determine and set the center of rotation 
        if opt_center or (rot_center == None):
            ### Using optimization method to automatically find the center
            # d.optimize_center()
            print("Optimizing center ...", end="")
            sys.stdout.flush()
            rot_center = tomopy.find_center(data, theta, ind=None, emission=True, init=None,
                                            tol=0.5, mask=True, ratio=1.)
            print("Done!")
            print("center = %s" % rot_center)
        if diag_center:
            ### Output the reconstruction results using a range of centers,
            ### and then manually find the optimal center.
            # d.diagnose_center()
            if not os.path.exists(diag_cent_dir):
                os.makedirs(diag_cent_dir)
            print("Testing centers ...", end="")
            sys.stdout.flush()
            tomopy.write_center(data, theta, dpath=diag_cent_dir, 
                                cen_range=[center_start, center_end, center_step], 
                                ind=None, emission=False, mask=False, ratio=1.)
            print("Done!")
        
        ## Flip odd frames
        if (cycle % 2):
            data[...] = data[...,::-1]
            rot_center = data.shape[-1] - rot_center_copy
        else:
            rot_center = rot_center_copy

        ## Reconstruction using FBP
        print("Running gridrec ...", end="")
        sys.stdout.flush()
        recon = tomopy.recon(data, theta, center=rot_center, emission=False, algorithm='gridrec',
                                # num_gridx=None, num_gridy=None, filter_name='shepp',
                                ncore=_ncore, nchunk=_nchunk)
        print("Done!")

        ## Collect background
        # if cycle == 0:
        #     bg = recon
        # elif cycle < 4:
        #     bg += recon
        # else:
        #     recon -= bg/4.

        # Write to stack of TIFFs.
        if not os.path.exists(recon_dir):
            os.makedirs(recon_dir)
        out_fname = recon_dir+"/"+out_prefix+"t_%d" % (cycle + cycle_offset)      
        if "hdf" in output: 
            hdf_fname = out_fname + ".hdf5"
            print("Writing reconstruction output file %s..." 
                 % hdf_fname, end="")
            sys.stdout.flush()
            tomopy.write_hdf5(recon, fname=hdf_fname, gname='exchange', overwrite=False)
            print("Done!")
        if "tif" in output:
            tiff_fname = out_fname + ".tiff"
            print("Writing reconstruction tiff files %s ..."
                    % tiff_fname, end="")
            sys.stdout.flush()
            tomopy.write_tiff_stack(recon, fname=tiff_fname, axis=0, digit=5, start=0, overwrite=False)
            print("Done!")
        if "bin" in output:
            bin_fname = out_fname + ".bin"
            print("Writing reconstruction to binary files %s..." 
                    % bin_fname, end="")
            sys.stdout.flush()
            recon.tofile(bin_fname)


def recon2(io_paras, data_paras, rot_center=None, normalize=True, stripe_removal=10, phase_retrieval=False, 
            opt_center=False, diag_center=False, output="tiff"):
    # Input and output
    datafile = io_paras.get('datafile')
    path2white = io_paras.get('path2white', datafile)
    path2dark = io_paras.get('path2dark', path2white)
    out_dir = io_paras.get('out_dir')
    diag_cent_dir = io_paras.get('diag_cent_dir', out_dir+"/center_diagnose/")
    recon_dir = io_paras.get('recon_dir', out_dir+"/recon/")
    out_prefix = io_paras.get('out_prefix', "recon_")
    
    # Parameters of dataset
    NumCycles = data_paras.get('NumCycles', 1) # Number of cycles used for recon
    ProjPerCycle = data_paras.get('ProjPerCycle') # Number of projections per cycle, N_theta
    cycle_offset = data_paras.get('cycle_offset', 0) # Offset in output cycle number
    proj_start = data_paras.get('proj_start', 0) # Starting projection of reconstruction 
    proj_step = data_paras.get('proj_step')
    z_start = data_paras.get('z_start', 0)
    z_end = data_paras.get('z_end', z_start+1)
    z_step = data_paras.get('z_step')
    x_start = data_paras.get('x_start')
    x_end = data_paras.get('x_end', x_start+1)
    x_step = data_paras.get('x_step')
    white_start = data_paras.get('white_start')
    white_end = data_paras.get('white_end')
    dark_start = data_paras.get('dark_start')
    dark_end = data_paras.get('dark_end')
    slice3 = slice(x_start, x_end, x_step)

    rot_center_copy = rot_center

    for cycle in xrange(NumCycles):
        # Set start and end of each cycle
        projections_start = cycle * ProjPerCycle + proj_start
        projections_end = projections_start + ProjPerCycle
        slice1 = slice(projections_start, projections_end, proj_step)
        # Distribute z slices to processes
        if z_step is None: # global z_step declaration is needed, because assignment is used below, therefore local variable is assumed by default.
            z_step = 1
        z_list = range(z_start, z_end, z_step)
        for i in range(_rank, len(z_list), _nprocs):
            z = z_list[i]
            slice2 = slice(z, z+1)
            slices = (slice1, slice2, slice3)
            white_slices = (slice(white_start, white_end), slice2, slice3)
            dark_slices = (slice(dark_start, dark_end), slice2, slice3)
            print("Running cycle #%s (projs %s to %s, z = %s) on process %s of %s" 
                % (cycle, projections_start, projections_end, z, _rank, _nprocs))
            
            # Read HDF5 file.
            print("Reading datafile %s..." % datafile, end="")
            sys.stdout.flush()
            data, white, dark = reader.read_aps_2bm(datafile, slices, white_slices, dark_slices, 
                                            path2white=path2white, path2dark=path2dark)
            theta = gen_theta(data.shape[0])
            print("Done!")
            print("Data shape = %s;\nwhite shape = %s;\ndark shape = %s." 
                % (data.shape, white.shape, dark.shape))
            
            
            # data = tomopy.focus_region(data, dia=1560, xcoord=1150, ycoord=1080, 
            #                 center=rot_center, pad=False, corr=True)
            # rot_center = None
            # print("Data shape = %s;\nwhite shape = %s;\ndark shape = %s." 
            #     % (data.shape, white.shape, dark.shape))

            ## Normalize dataset using data_white and data_dark
            if normalize:
                print("Normalizing data ...")
                # white = white.mean(axis=0).reshape(-1, *data.shape[1:])
                # dark = dark.mean(axis=0).reshape(-1, *data.shape[1:])
                # data = (data - dark) / (white - dark)
                data = tomopy.normalize(data, white, dark, cutoff=None, ncore=_ncore, nchunk=None)[...]
        
            ## Remove stripes caused by dead pixels in the detector
            if stripe_removal:
                print("Removing stripes ...")
                data = tomopy.remove_stripe_fw(data, level=stripe_removal, wname='db5', sigma=2,
                                        pad=True, ncore=_ncore, nchunk=None)
                # data = tomopy.remove_stripe_ti(data, nblock=0, alpha=1.5, 
                #                                 ncore=None, nchunk=None)

    #        # Show preprocessed projection
    #        plt.figure("%s-prep" % projections_start)
    #        plt.imshow(d.data[0,:,:], cmap=cm.Greys_r)
    #        plt.savefig(out_dir+"/preprocess/%s-prep.jpg" 
    #                    % projections_start)
    #        # plt.show()
    #        continue

            ## Phase retrieval
            if phase_retrieval:
                print("Retrieving phase ...")
                data = tomopy.retrieve_phase(data,
                                        pixel_size=6.5e-5, dist=33, energy=30,
                                        alpha=1e-3, pad=True, ncore=_ncore, nchunk=None)
            
            ## Determine and set the center of rotation 
            if opt_center: # or (rot_center == None):
                ### Using optimization method to automatically find the center
                # d.optimize_center()
                print("Optimizing center ...", end="")
                sys.stdout.flush()
                rot_center = tomopy.find_center(data, theta, ind=None, emission=True, init=None,
                                                tol=0.5, mask=True, ratio=1.)
                print("Done!")
                print("center = %s" % rot_center)
            if diag_center:
                ### Output the reconstruction results using a range of centers,
                ### and then manually find the optimal center.
                # d.diagnose_center()
                if not os.path.exists(diag_cent_dir):
                    os.makedirs(diag_cent_dir)
                print("Testing centers ...", end="")
                sys.stdout.flush()
                tomopy.write_center(data, theta, dpath=diag_cent_dir, 
                                    cen_range=[center_start, center_end, center_step], 
                                    ind=None, emission=False, mask=False, ratio=1.)
                print("Done!")

            ## Flip odd frames
            if (cycle % 2):
                data[...] = data[...,::-1]
                rot_center = data.shape[-1] - rot_center_copy
            else:
                rot_center = rot_center_copy

            ## Reconstruction using FBP
            print("Running gridrec ...", end="")
            sys.stdout.flush()
            recon = tomopy.recon(data, theta, center=rot_center, emission=False, algorithm='gridrec',
                                    # num_gridx=None, num_gridy=None, filter_name='shepp',
                                    ncore=_ncore, nchunk=_nchunk)
            print("Done!")

            ## Collect background
            # if cycle == 0:
            #     bg = recon
            # elif cycle < 4:
            #     bg += recon
            # else:
            #     recon -= bg/4.

            # Write to stack of TIFFs.
            if not os.path.exists(recon_dir):
                os.makedirs(recon_dir)
            out_fname = recon_dir+"/"+out_prefix+"t_%d_z_%d" % (cycle + cycle_offset, z)      
            if "hdf" in output: 
                hdf_fname = out_fname + ".hdf5"
                print("Writing reconstruction output file %s..." 
                     % hdf_fname, end="")
                sys.stdout.flush()
                tomopy.write_hdf5(recon, fname=hdf_fname, gname='exchange', overwrite=False)
                print("Done!")
            if "tif" in output:
                tiff_fname = out_fname + ".tiff"
                print("Writing reconstruction tiff files %s ..."
                        % tiff_fname, end="")
                sys.stdout.flush()
                tomopy.write_tiff(recon, fname=tiff_fname, overwrite=False)
                print("Done!")
            if "bin" in output:
                bin_fname = out_fname + ".bin"
                print("Writing reconstruction to binary files %s..." 
                        % bin_fname, end="")
                sys.stdout.flush()
                recon.tofile(bin_fname)
    if _usempi:
        comm.Barrier()
    if _rank == 0:
        print("All done!")


def center(io_paras, data_paras, center_start, center_end, center_step, diag_cycle=0, 
            mode='diag', normalize=True, stripe_removal=10, phase_retrieval=False):
    
    # Input and output
    datafile = io_paras.get('datafile')
    path2white = io_paras.get('path2white', datafile)
    path2dark = io_paras.get('path2dark', path2white)
    out_dir = io_paras.get('out_dir')
    diag_cent_dir = io_paras.get('diag_cent_dir', out_dir+"/center_diagnose/")
    recon_dir = io_paras.get('recon_dir', out_dir+"/recon/")
    out_prefix = io_paras.get('out_prefix', "recon_")

    # Parameters of dataset
    NumCycles = data_paras.get('NumCycles', 1) # Number of cycles used for recon
    ProjPerCycle = data_paras.get('ProjPerCycle') # Number of projections per cycle, N_theta
    cycle_offset = data_paras.get('cycle_offset', 0) # Offset in output cycle number
    proj_start = data_paras.get('proj_start', 0) # Starting projection of reconstruction 
    proj_step = data_paras.get('proj_step')
    z_start = data_paras.get('z_start', 0)
    z_end = data_paras.get('z_end', z_start+1)
    z_step = data_paras.get('z_step')
    x_start = data_paras.get('x_start')
    x_end = data_paras.get('x_end', x_start+1)
    x_step = data_paras.get('x_step')
    white_start = data_paras.get('white_start')
    white_end = data_paras.get('white_end')
    dark_start = data_paras.get('dark_start')
    dark_end = data_paras.get('dark_end')

    # TIMBIR parameters
    NumSubCycles = data_paras.get('NumSubCycles', 1) # Number of subcycles in one cycle, K
    SlewSpeed = data_paras.get('SlewSpeed', 0) # In deg/s
    MinAcqTime = data_paras.get('MinAcqTime', 0) # In s
    TotalNumCycles = data_paras.get('TotalNumCycles', 1) # Total number of cycles in the full scan data
    ProjPerRecon = data_paras.get('ProjPerRecon', ProjPerCycle) # Number of projections per reconstruction

    # Calculate thetas for interlaced scan
    theta = gen_theta_timbir(NumSubCycles, ProjPerCycle, SlewSpeed, MinAcqTime, TotalNumCycles)
    if ProjPerRecon is None:
        ProjPerCycle = theta.size//TotalNumCycles
    else:
        ProjPerCycle = ProjPerRecon

    print("Will use %s projections per reconstruction." % ProjPerCycle)

    # Set start and end of each subcycle
    projections_start = diag_cycle * ProjPerCycle + proj_start
    projections_end = projections_start + ProjPerCycle
    slice1 = slice(projections_start, projections_end, proj_step)
    slice2 = slice(z_start, z_end, z_step)
    slice3 = slice(x_start, x_end, x_step)
    slices = (slice1, slice2, slice3)
    white_slices = (slice(white_start, white_end), slice2, slice3)
    dark_slices = (slice(dark_start, dark_end), slice2, slice3)
    print("Running center diagnosis (projs %s to %s)" 
        % (projections_start, projections_end))
    
    # Read HDF5 file.
    print("Reading datafile %s..." % datafile, end="")
    sys.stdout.flush()
    data, white, dark = reader.read_aps_2bm(datafile, slices, white_slices, dark_slices, 
                                    path2white=path2white, path2dark=path2dark)
    data += 1
    # theta = gen_theta(data.shape[0])
    print("Done!")
    print("Data shape = %s;\nwhite shape = %s;\ndark shape = %s." 
        % (data.shape, white.shape, dark.shape))
    
    ## Normalize dataset using data_white and data_dark
    if normalize:
        data = tomopy.normalize(data, white, dark, cutoff=None, ncore=_ncore, nchunk=None)

    ## Remove stripes caused by dead pixels in the detector
    if stripe_removal:
        data = tomopy.remove_stripe_fw(data, level=stripe_removal, wname='db5', 
                                        sigma=2, pad=True, ncore=None, nchunk=None)
        # data = tomopy.remove_stripe_ti(data, nblock=0, alpha=1.5, 
        #                                 ncore=None, nchunk=None)
    
#        # Show preprocessed projection
#        plt.figure("%s-prep" % projections_start)
#        plt.imshow(d.data[0,:,:], cmap=cm.Greys_r)
#        plt.savefig(out_dir+"/preprocess/%s-prep.jpg" 
#                    % projections_start)
#        # plt.show()
#        continue

    ## Phase retrieval
    if phase_retrieval:
        data = tomopy.retrieve_phase(data,
                    pixel_size=1.1e-4, dist=5, energy=25.7,
                    alpha=1e-3, pad=True, ncore=_ncore, nchunk=None)
    
    ## Determine and set the center of rotation
    ### Using optimization method to automatically find the center
    # d.optimize_center()
    if 'opti' in mode:
        print("Optimizing center ...", end="")
        sys.stdout.flush()
        rot_center = tomopy.find_center(data, theta, ind=None, emission=True, init=None,
                                        tol=0.5, mask=True, ratio=1.)
        print("Done!")
        print("center = %s" % rot_center)
    ### Output the reconstruction results using a range of centers,
    ### and then manually find the optimal center.
    if 'diag' in mode:
        if not os.path.exists(diag_cent_dir):
            os.makedirs(diag_cent_dir)
        print("Testing centers ...", end="")
        sys.stdout.flush()
        tomopy.write_center(data, theta, dpath=diag_cent_dir, 
                            cen_range=[center_start, center_end, center_step], 
                            ind=None, emission=False, mask=False, ratio=1.)
        print("Done!")


def sweep_parameter(levels, sigmas, *args, **kwargs):
    for level in levels:
        for sigma in sigmas:
            print("Running level = %s, sigma = %s" % (level, sigma))
            kwargs['stripe_removal'] = level
            kwargs['stripe_sigma'] = sigma
            out_prefix = "recon_sr_l_%s_s_%s_" % (level, sigma)
            args[0]['out_prefix'] = out_prefix
            recon3(*args, **kwargs)


def recon3(io_paras, data_paras, rot_center=None, normalize=True, stripe_removal=10, stripe_sigma=2, phase_retrieval=False, 
            opt_center=False, diag_center=False, output="tiff", z_recon_size=None):
    # Input and output
    datafile = io_paras.get('datafile')
    path2white = io_paras.get('path2white', datafile)
    path2dark = io_paras.get('path2dark', path2white)
    out_dir = io_paras.get('out_dir')
    diag_cent_dir = io_paras.get('diag_cent_dir', out_dir+"/center_diagnose/")
    recon_dir = io_paras.get('recon_dir', out_dir+"/recon/")
    out_prefix = io_paras.get('out_prefix', "recon_")

    # Parameters of dataset
    NumCycles = data_paras.get('NumCycles', 1) # Number of cycles used for recon
    ProjPerCycle = data_paras.get('ProjPerCycle') # Number of projections per cycle, N_theta
    cycle_offset = data_paras.get('cycle_offset', 0) # Offset in output cycle number
    proj_start = data_paras.get('proj_start', 0) # Starting projection of reconstruction 
    proj_step = data_paras.get('proj_step')
    z_start = data_paras.get('z_start', 0)
    z_end = data_paras.get('z_end', z_start+1)
    z_step = data_paras.get('z_step')
    x_start = data_paras.get('x_start')
    x_end = data_paras.get('x_end', x_start+1)
    x_step = data_paras.get('x_step')
    white_start = data_paras.get('white_start')
    white_end = data_paras.get('white_end')
    dark_start = data_paras.get('dark_start')
    dark_end = data_paras.get('dark_end')

    # TIMBIR parameters
    NumSubCycles = data_paras.get('NumSubCycles', 1) # Number of subcycles in one cycle, K
    SlewSpeed = data_paras.get('SlewSpeed', 0) # In deg/s
    MinAcqTime = data_paras.get('MinAcqTime', 0) # In s
    TotalNumCycles = data_paras.get('TotalNumCycles', 1) # Total number of cycles in the full scan data
    ProjPerRecon = data_paras.get('ProjPerRecon', ProjPerCycle) # Number of projections per reconstruction

    # Calculate thetas for interlaced scan
    theta = gen_theta_timbir(NumSubCycles, ProjPerCycle, SlewSpeed, MinAcqTime, TotalNumCycles)
    if ProjPerRecon is None:
        ProjPerCycle = theta.size//TotalNumCycles
    else:
        ProjPerCycle = ProjPerRecon

    print("Will use %s projections per reconstruction." % ProjPerCycle)

    # Distribute z slices to processes
    if z_step is None: 
        z_step = 1
    
    z_pool = get_pool(z_start, z_end, z_step, z_chunk_size=z_recon_size, fmt='slice')

    slice3 = slice(x_start, x_end, x_step)

    rot_center_copy = rot_center

    for cycle in xrange(NumCycles):
        
        # Set start and end of each cycle
        projections_start = cycle * ProjPerCycle + proj_start
        projections_end = projections_start + ProjPerCycle
        slice1 = slice(projections_start, projections_end, proj_step)

        # Setup continuous output
        if "cont" in output:
            if not os.path.exists(recon_dir):
                os.makedirs(recon_dir)
            cont_fname = recon_dir+"/"+out_prefix+"t_%d_z_%d_%d.bin" \
                        % (cycle + cycle_offset, z_start, z_end)
            cont_file = file(cont_fname, 'wb')
        # Distribute z slices to processes
        for i in range(_rank, len(z_pool), _nprocs):
            slice2 = z_pool[i]
            slices = (slice1, slice2, slice3)
            white_slices = (slice(white_start, white_end), slice2, slice3)
            dark_slices = (slice(dark_start, dark_end), slice2, slice3)
            print("Running cycle #%s (projs %s to %s, z = %s - %s) on process %s of %s" 
                % (cycle, projections_start, projections_end, slice2.start, slice2.stop, _rank, _nprocs))
            
            # Read HDF5 file.
            print("Reading datafile %s..." % datafile, end="")
            sys.stdout.flush()
            data, white, dark = reader.read_aps_2bm(datafile, slices, white_slices, dark_slices, 
                                            path2white=path2white, path2dark=path2dark)
            # data += 1
            # theta = gen_theta(data.shape[0])
            print("Done!")
            print("Data shape = %s;\nwhite shape = %s;\ndark shape = %s." 
                % (data.shape, white.shape, dark.shape))
            
            
            # data = tomopy.focus_region(data, dia=1560, xcoord=1150, ycoord=1080, 
            #                 center=rot_center, pad=False, corr=True)
            # rot_center = None
            # print("Data shape = %s;\nwhite shape = %s;\ndark shape = %s." 
            #     % (data.shape, white.shape, dark.shape))

            ## Normalize dataset using data_white and data_dark
            if normalize:
                print("Normalizing data ...")
                # white = white.mean(axis=0).reshape(-1, *data.shape[1:])
                # dark = dark.mean(axis=0).reshape(-1, *data.shape[1:])
                # data = (data - dark) / (white - dark)
                data = tomopy.normalize(data, white, dark, cutoff=None, ncore=_ncore, nchunk=_nchunk)[...]
        
            ## Remove stripes caused by dead pixels in the detector
            if stripe_removal:
                print("Removing stripes ...")
                data = tomopy.remove_stripe_fw(data, level=stripe_removal, wname='db5', sigma=stripe_sigma,
                                        pad=True, ncore=_ncore, nchunk=_nchunk)
                # data = tomopy.remove_stripe_ti(data, nblock=0, alpha=1.5, 
                #                                 ncore=None, nchunk=None)

    #        # Show preprocessed projection
    #        plt.figure("%s-prep" % projections_start)
    #        plt.imshow(d.data[0,:,:], cmap=cm.Greys_r)
    #        plt.savefig(out_dir+"/preprocess/%s-prep.jpg" 
    #                    % projections_start)
    #        # plt.show()
    #        continue

            ## Phase retrieval
            if phase_retrieval:
                print("Retrieving phase ...")
                data = tomopy.retrieve_phase(data,
                                        pixel_size=1.1e-4, dist=6, energy=25.7,
                                        alpha=1e-2, pad=True, ncore=_ncore, nchunk=_nchunk)
            
            ## Determine and set the center of rotation 
            if opt_center: # or (rot_center == None):
                ### Using optimization method to automatically find the center
                # d.optimize_center()
                print("Optimizing center ...", end="")
                sys.stdout.flush()
                rot_center = tomopy.find_center(data, theta, ind=None, emission=True, init=None,
                                                tol=0.5, mask=True, ratio=1.)
                print("Done!")
                print("center = %s" % rot_center)
            if diag_center:
                ### Output the reconstruction results using a range of centers,
                ### and then manually find the optimal center.
                # d.diagnose_center()
                if not os.path.exists(diag_cent_dir):
                    os.makedirs(diag_cent_dir)
                print("Testing centers ...", end="")
                sys.stdout.flush()
                tomopy.write_center(data, theta, dpath=diag_cent_dir, 
                                    cen_range=[center_start, center_end, center_step], 
                                    ind=None, emission=False, mask=False, ratio=1.)
                print("Done!")

            ## Flip odd frames
#            if (cycle % 2):
#                data[...] = data[...,::-1]
#                rot_center = data.shape[-1] - rot_center_copy
#            else:
#                rot_center = rot_center_copy

            ## Reconstruction using FBP
            print("Running gridrec ...", end="")
            sys.stdout.flush()
            recon = tomopy.recon(data, theta[slice1], center=rot_center, emission=False, algorithm='gridrec',
                                    # num_gridx=None, num_gridy=None, filter_name='shepp',
                                    ncore=_ncore, nchunk=_nchunk)
            print("Done!")

            ## Collect background
            # if cycle == 0:
            #     bg = recon
            # elif cycle < 4:
            #     bg += recon
            # else:
            #     recon -= bg/4.

            # Write to stack of TIFFs.
            if not os.path.exists(recon_dir):
                os.makedirs(recon_dir)
            out_fname = recon_dir+"/"+out_prefix+"t_%d_z_" % (cycle + cycle_offset)      
            if "hdf" in output: 
                hdf_fname = out_fname + "%d_%d.hdf5" % (slice2.start, slice2.stop)
                print("Writing reconstruction output file %s..." 
                     % hdf_fname, end="")
                sys.stdout.flush()
                tomopy.write_hdf5(recon, fname=hdf_fname, gname='exchange', overwrite=False)
                print("Done!")
            if "tif" in output:
                if "stack" in output: # single stacked file for multiple z
                    tiff_fname = out_fname + "%d_%d.tiff" % (slice2.start, slice2.stop)
                    print("Writing reconstruction tiff files %s ..."
                            % tiff_fname, end="")
                    sys.stdout.flush()
                    tomopy.write_tiff(recon, fname=tiff_fname, overwrite=False)
                    print("Done!")
                    
                else: # separate files for different z
                    for iz, z in enumerate(range(slice2.start, slice2.stop, slice2.step)):
                        tiff_fname = out_fname + "%d.tiff" % z
                        print("Writing reconstruction tiff files %s ..."
                                % tiff_fname, end="")
                        sys.stdout.flush()
                        tomopy.write_tiff(recon[iz], fname=tiff_fname, overwrite=False)
                        print("Done!")
            if "bin" in output:
                bin_fname = out_fname + "%d_%d.bin" % (slice2.start, slice2.stop)
                print("Writing reconstruction to binary files %s..." 
                        % bin_fname, end="")
                sys.stdout.flush()
                recon.tofile(bin_fname)
            if "cont" in output:
                print("Writing reconstruction to binary files %s..." 
                        % cont_fname, end="")
                sys.stdout.flush()
                recon.tofile(cont_file)
                print("Done!")
        if "cont" in output:
            cont_file.close()
    
    if _usempi:
        comm.Barrier()
    if _rank == 0:
        print("All done!")


def get_pool(z_start, z_end, z_step=1, z_chunk_size=None, fmt='slice'):
    if fmt == 'list': 
        z_list = range(z_start, z_end, z_step)
        if z_chunk_size is not None:
            npool = len(z_list) // z_chunk_size + 1
            z_pool = [z_list[i:i+z_chunk_size] for i in range(npool)]
        else:
            z_pool = z_list
    else:
        nz = (z_end - z_start - 1) // z_step + 1
        if z_chunk_size is not None:
            npool = nz // z_chunk_size + bool(nz % z_chunk_size)
            z_pool = [slice(z_start + i * z_step * z_chunk_size, 
                            z_start + (i + 1) * z_step * z_chunk_size,
                            z_step) for i in range(npool - 1)]
            z_pool += [slice(z_start + (npool - 1) * z_step * z_chunk_size,
                            z_end, z_step)]
        else:
            z_pool = [slice(z_start, z_end, z_step)]
    return z_pool


def image_normalize(data, dtype=int, vmax=255, cutoffmin=None, cutoffmax=None):
    if cutoffmin is not None:
        data[data<cutoffmin] = cutoffmin
    if cutoffmax is not None:
        data[data>cutoffmax] = cutoffmax

    data = (data - data.min())/(data.max()-data.min())*vmax

    if dtype is int:
        data = np.array(data * (vmax+1), dtype=int)
        data[data==(vmax+1)] = vmax

    return data
