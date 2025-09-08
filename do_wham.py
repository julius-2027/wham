import numpy as np
import matplotlib.pyplot as plt
import argparse
from wham import WHAM

kB = 0.001987191 # Boltzmann constant in kcal/(mol*K)
T_default = 310 # default temperature in Kelvin (approximate body temp in humans)

parser = argparse.ArgumentParser()

parser.add_argument('windows_file')
parser.add_argument('output_name')

parser.add_argument('-t_min', '--t_min', type=float, default=None)
parser.add_argument('-t_max', '--t_max', type=float, default=None)
parser.add_argument('-n_part', '--num_partitions', type=int, default=0)

parser.add_argument('-xi_min', '--xi_min', type=float, default=None)
parser.add_argument('-xi_max', '--xi_max', type=float, default=None)
parser.add_argument('-n_bins', '--num_bins', type=int, default=200)

parser.add_argument('-tol', '--tolerance', type=float, default=0.0001)
parser.add_argument('-kT', '--kT', type=float, default=kB*T_default)

parser.add_argument('-print_progress', action='store_true') #'store_true' means default is False
parser.add_argument('-plot', action='store_true')

args = parser.parse_args()

print('Input: ')
print(vars(args))

print('\n')
wham = WHAM(args.windows_file, kT=args.kT, print_progress=args.print_progress)

print('\n')
wham.load_data(t_min=args.t_min, t_max=args.t_max, print_progress=args.print_progress)

largest_t_start = max(wham.arr_t_start)
smallest_t_end = min(wham.arr_t_end)

n_part = int(args.num_partitions)
arr_t_part = np.linspace(largest_t_start, smallest_t_end, max(2,n_part+1))

# all data (excluding time values that are not present in all windows)
print('\nAll data: t='+str(arr_t_part[0])+' to t='+str(arr_t_part[-1]))
wham.select_data(t_min=arr_t_part[0], t_max=arr_t_part[-1])
wham.hist_data(xi_min=args.xi_min, xi_max=args.xi_max, n_bins=args.num_bins, print_progress=args.print_progress)
wham.combine_windows(tolerance=args.tolerance, print_progress=args.print_progress)
wham.pdf_to_pmf()

results = [wham.xi_binCenters, wham.pmf]

xi_min = wham.xi_binEdges[0]
xi_max = wham.xi_binEdges[-1]

# partitions
for i in range(n_part):
    print('\nPartition '+str(i+1)+' of '+str(n_part)+': t='+str(arr_t_part[i])+' to t='+str(arr_t_part[i+1]))
    wham.select_data(t_min=arr_t_part[i], t_max=arr_t_part[i+1])
    wham.hist_data(xi_min=xi_min, xi_max=xi_max, n_bins=args.num_bins, print_progress=args.print_progress)
    wham.combine_windows(tolerance=args.tolerance, print_progress=args.print_progress)
    wham.pdf_to_pmf()
    results.append(wham.pmf)

results = np.array(results)
np.savetxt(args.output_name, results.T)

if args.plot:
    plt.plot(results[0], results[1], label='all data', color='black', linestyle='--')

    for i in range(n_part):
        plt.plot(results[0], results[i+2], label='partition '+str(i+1)+' of '+str(n_part))

    plt.grid()
    plt.legend()
    plt.show()
