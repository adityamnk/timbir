# -*- coding: utf-8 -*-
import os


# Input and output
datafile = "/media/TOURO Desk Pro/S1A1_Heating_T627_N1_K8_nTheta3200_90DegPerSec_2msExpTime_edge_10x_50mm_Rolling_10umLuAG_1.8cmUnd_24.0mmGap_25.7keV_32ID_CHutch/proj_0336.hdf"
path2white = None #"/media/TOURO Desk Pro/S1A1_Heating_T627_N1_K8_nTheta3200_90DegPerSec_2msExpTime_edge_10x_50mm_Rolling_10umLuAG_1.8cmUnd_24.0mmGap_25.7keV_32ID_CHutch/proj_0337.hdf"
path2dark = None #"/media/TOURO Desk Pro/S1A1_Heating_T627_N1_K8_nTheta3200_90DegPerSec_2msExpTime_edge_10x_50mm_Rolling_10umLuAG_1.8cmUnd_24.0mmGap_25.7keV_32ID_CHutch/proj_0338.hdf"
out_dir = os.path.dirname(datafile)+"/"+os.path.basename(datafile).split('.')[0]+"/"
diag_cent_dir = out_dir+"/center_diagnose/"
recon_dir = out_dir+"/recon/"
out_prefix = "recon_norm_phase_"



# Parameters of dataset
NumCycles = 1 # Number of cycles used for recon
ProjPerCycle = 3200 # Number of projections per cycle, N_theta
cycle_offset = 0 # Offset in output cycle number
proj_start = 1 # Starting projection of reconstruction 
# proj_end = 19000 # Not implemented
proj_step = None
z_start = 0
z_end = 800
z_step = 50
x_start = 0
x_end = x_start+2560
x_step = None
white_start = 1
white_end = None
dark_start = 1
dark_end = None

# TIMBIR parameters
NumSubCycles = 8 # Number of subcycles in one cycle, K
SlewSpeed = 90
MinAcqTime = 1e-4
TotalNumCycles = 1 # Total number of cycles in the full scan data
ProjPerRecon = None # Number of projections per reconstruction

# Diagnose center
diag_cycle = 0
center_start = 550-x_start
center_end = 650-x_start
center_step = 5

# Set center of rotation
rot_center = 2560/2-x_start


# Program settings
usempi = 0
ncore = 24
nchunk = None
z_recon_size = None


from recon import *

set_mpi(usempi)
set_mp(ncore, nchunk)
# set_z_recon_size(z_recon_size)
show_settings()


from kwarg import kwarg

# Input and output
io_paras = kwarg(datafile = datafile,
                    path2white = path2white,
                    path2dark = path2dark,
                    out_dir = out_dir,
                    diag_cent_dir = diag_cent_dir,
                    recon_dir = recon_dir,
                    out_prefix = out_prefix)
# Parameters of dataset
data_paras = kwarg(NumCycles = NumCycles, # Number of cycles used for recon
                    ProjPerCycle = ProjPerCycle, # Number of projections per cycle, N_theta
                    cycle_offset = cycle_offset, # Offset in output cycle number
                    proj_start = proj_start, # Starting projection of reconstruction 
                    proj_step = proj_step,
                    z_start = z_start,
                    z_end = z_end,
                    z_step = z_step,
                    x_start = x_start,
                    x_end = x_end,
                    x_step = x_step,
                    white_start = white_start,
                    white_end = white_end,
                    dark_start = dark_start,
                    dark_end = dark_end,
                    # TIMBIR parameters
                    NumSubCycles = NumSubCycles, # Number of subcycles in one cycle, K
                    SlewSpeed = SlewSpeed,
                    MinAcqTime = MinAcqTime,
                    TotalNumCycles = TotalNumCycles, # Total number of cycles in the full scan data
                    ProjPerRecon = ProjPerRecon # Number of projections per reconstruction
                    )


if __name__ == "__main__":
    # center(io_paras, data_paras, center_start, center_end, center_step, diag_cycle=diag_cycle, normalize=False, stripe_removal=False, mode='diag')#+opti')
    # recon2(io_paras, data_paras, rot_center=rot_center, stripe_removal=20, normalize=True, phase_retrieval=False, output="tiff")
    recon3(io_paras, data_paras, rot_center=rot_center, stripe_removal=0, normalize=False, phase_retrieval=True, output="tiffstack", z_recon_size=z_recon_size)
    # sweep_parameter()
