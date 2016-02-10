# -*- coding: utf-8 -*-
import os


# Input and output
datafile = "/Volumes/HGST8/S6A1_Heating_T540_N1_K1_nTheta1024_720DegPerSec_200usExpTime_dimax_10x_50mm_Rolling_100umLuAG_1.8cmUnd_22.2mmGap_25.7keV_32ID_CHutch/proj_0022.hdf"
path2white = None # "/mnt/quest/projects/b1003/yue/data/APS14_AlCu_07.hdf"
path2dark = path2white
out_dir = os.path.dirname(datafile)+"/"+os.path.basename(datafile).split('.')[0]+"/"
diag_cent_dir = out_dir+"/center_diagnose/"
recon_dir = out_dir+"/recon/"
out_prefix = "recon_"


# Parameters of dataset
NumCycles = 1 # Number of cycles used for recon
ProjPerCycle = 1024 # Number of projections per cycle, N_theta
cycle_offset = 0 # Offset in output cycle number
proj_start = 0 # Starting projection of reconstruction 
# proj_end = 19000 # Not implemented
proj_step = None
z_start = 0
z_end = 1800
z_step = 50
x_start = 0
x_end = x_start+1200
x_step = None
white_start = 1
white_end = None
dark_start = 1
dark_end = None

# Diagnose center
diag_cycle = 1
center_start = 1325-x_start
center_end = 1345-x_start
center_step = 1

# Set center of rotation
rot_center = 610-x_start


# Program settings
usempi = 0
ncore = None
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
                    dark_end = dark_end)


if __name__ == "__main__":
    # center(io_paras, data_paras, center_start, center_end, center_step, diag_cycle=diag_cycle, normalize=False, stripe_removal=False, mode='diag+opti')
    # recon2(io_paras, data_paras, rot_center=rot_center, stripe_removal=20, normalize=True, phase_retrieval=False, output="tiff")
    recon3(io_paras, data_paras, rot_center=rot_center, stripe_removal=0, normalize=False, phase_retrieval=False, output="tiffstack", z_recon_size=z_recon_size)
    # sweep_parameter()
