import numpy as np
import glob


dir2data = "/Volumes/HGST5/APS1504_AlCu/recon/D_Coarsening_2/recon_2/"
prefix = "recon_"
dir2out = "/Volumes/HGST5/APS1504_AlCu/recon/D_Coarsening_2/combine_test/"


flist = glob.glob(dir2data+prefix+"*.bin")


def parse_filename(filename):
    t, z = filename.split('.')[0].split('_')[-3::2]
    return int(t), int(z)


fdict = {}
for f in flist:
    t, z = parse_filename(f)
    try:
        fdict[t].append(z)
    except KeyError:
        fdict[t] = [z]

for t, zlist in fdict.items():
    out_name = dir2out+prefix+"t_%d.bin" % t
    print("Opening output file %s" % out_name)
    fout = open(out_name, 'wb')
    for z in sorted(zlist):
        data_name = dir2data+prefix+"t_%d_z_%d.bin" % (t, z)
        print("Reading data file %s" % data_name)
        fin = open(data_name, 'rb')
        fout.write(fin.read())
        fin.close()
    fout.close()
