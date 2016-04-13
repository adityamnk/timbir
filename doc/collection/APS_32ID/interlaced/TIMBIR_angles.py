import numpy as np
import matplotlib.pyplot as plt
from sys import argv

# Bit-reversal function
# def _bit_reverse(n):
#     return int(bin(n)[:1:-1], 2)

def _bit_reverse(n, width=8):
    b = '{:0{width}b}'.format(n, width=width)
    return int(b[::-1], 2)

# bit_reverse = np.frompyfunc(_bit_reverse, 1, 1)
bit_reverse = np.vectorize(_bit_reverse, excluded='width')

# Function generating view angles for one cycle
def gen_theta_one_cycle(K, N_theta):
    n = np.arange(N_theta)
    width = len(bin(K))-3
#    return n*K//N_theta
#    return (n*K//N_theta) % K
#    return bit_reverse((n*K//N_theta) % K)
    return (n*K + bit_reverse((n*K//N_theta) % K, width=width)) * (np.pi/N_theta)

# Function generating view angles for mutiple cycle
def gen_theta(K, N_theta, TotalNumCycles=1):
    offset = 0
    thetas = np.array([],dtype=np.float)
    for i in range(TotalNumCycles):
        thetas = np.append(thetas, gen_theta_one_cycle(K, N_theta)+offset)
        offset += np.pi*K
    return thetas

def calc_dropped_angles(thetas, dtheta_min, verbose=True, TotalNumCycles=1):
    dthetas = thetas[1:] - thetas[:-1]
    dropped_angles = np.append([False], dthetas<dtheta_min)
    if verbose:
        print("Minimum angle between frames %s" % np.min(dthetas))
        print("Maximum angle between frames %s" % np.max(dthetas))
        if TotalNumCycles > 1:
            n_one_cycle = dthetas.size//TotalNumCycles
            print("Dropped %s frames in the first cycle, %s frames in subsequent cycles." % 
                (np.sum(dropped_angles[:n_one_cycle]), np.sum(dropped_angles[n_one_cycle:2*n_one_cycle])))
        print("Total number of dropped frames %s" % np.sum(dropped_angles))
    return dropped_angles
    
if __name__ == "__main__":
    if len(argv) >= 3:
        K = int(argv[1])
        N_theta = int(argv[2])
    else:
        K = 4
        N_theta = 16
    
    TotalNumCycles = 1    
    for i in range(3,len(argv)):
        if 'n_cycles=' in argv[i]:
            TotalNumCycles = int(argv[i].split('=')[-1])
    
    thetas = gen_theta(K, N_theta, TotalNumCycles)

    dtheta_min = 0
    for i in range(3,len(argv)):
        if 'dtheta_min=' in argv[i]:
            dtheta_min = float(argv[i].split('=')[-1])
            dropped_angles = calc_dropped_angles(thetas, dtheta_min, TotalNumCycles=TotalNumCycles)

    if not ('nowrap' in argv):
        thetas = thetas % np.pi
    if 'sort' in argv:
        indx = np.argsort(thetas)
        thetas = thetas[indx]
        if dtheta_min != 0:
        	dropped_angles = dropped_angles[indx]
    plt.plot(thetas, 'ro--', label='Viewing angles')
    if dtheta_min != 0:
        plt.plot(np.arange(len(thetas))[dropped_angles], thetas[dropped_angles], 'go', label='Dropped angles')

    plt.xlabel('Viewing angle index')
    plt.ylabel('Viewing angle (Rad)')
    plt.legend(loc='upper center', ncol=2, bbox_to_anchor=(0.5,1.12), fancybox=True, shadow=True)
    plt.show()
