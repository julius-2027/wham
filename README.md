# wham
Python script implementing WHAM (Weighted Histogram Analysis Method). Based on Roux, B. (1995).[^fn1] Calculates PMF (Potential of Mean Force) along a single coordinate using data from umbrella sampling MD (Molecular Dynamics) simulations.

## Requirements
**NumPy** (https://numpy.org/). Used for loading and saving data files, histogramming the data, and implementing the WHAM algorithm.

**os** (https://docs.python.org/3/library/os.html). Built-in Python module. Used for efficiently reading the last line of simulation data files to determine the simulation end time.

**argparse** (https://docs.python.org/3/library/argparse.html). Built-in Python module. Used for parsing command line arguments in `do_wham.py`.

**Matplotlib** (https://matplotlib.org/). Used for plotting the output of the calculation if the optional `-plot` flag of `do_wham.py` is set.

`wham.py` must be included in the same directory as `do_wham.py` since `do_wham.py` imports the `WHAM` class from `wham.py`.

## Usage of do_wham.py
```
python do_wham.py <windows_file> <output_name>
```

`<windows_file>` specifies the path to a file describing the results of the umbrella sampling. The number of rows in `<windows_file>` should be equal to the number of windows. Each row should be formatted as follows:

```
<data_file> <window_center> 0 0 <k_spring> 0 0
```

where `<data_file>` contains a list of coordinate values over time in the MD simulation for the row's window:

```
<time in frame 1> <value of coordinate in frame 1>
...
<time in frame n> <value of coordinate in frame n>
```

The script assumes that harmonic restraints were applied to the coordinate of interest, $\xi$ (xi), in each window. These restraints have the form $w_i(\xi) = \frac{1}{2} k_i (\xi - \xi_i)^2$ where $\xi_i$ is specified by `<window center>`, and $k_i$ is specified by `<k_spring>`. The units of `<k_spring>` should be consistent with the units of energy chosen for `-kT` and `-tol`:

$\[k_\text{spring}\] = \frac{ \[\text{energy}\] }{ [\xi]^2 } $

`output_name` specifies a name (without file extension) for the file to which the output will be written.

Only the overlapping portions of MD simulations from different windows are considered for calculating output. That is, if one window was simulated longer than another window, then the additional simulation time in the longer window is thrown out. Thus, it is ideal to run simulations for the same time in each window. From here on, the simulation data obtained after exluding any non-overlapping portions will be referred to as *all data*

The output is written in the following format:

```
<center of bin 1> <PMF, all data> <PMF, partition 1 of n> ... <PMF, partition n of n>
...
<center of bin n> <PMF, all data> <PMF, partition 1 of n> ... <PMF, partition n of n>
```

The PMF over all data is calculated using the definition of all data above. If the optional argument `-n_part` is set to a value greater than its default of 0, then the simulation data from each window is partitioned into n_part equal-sized time intervals. The PMF is calculated for each partition using only the simulation data from that partition's corresponding time interval.

Under the assumption that the system is at equilibrium during the umbrella sampling simulations, we expect the PMF to be the same for each partition. Differences in the PMF for different partitions tell us about uncertainties in the PMF calculated over all simulation data. For example, the uncertainty may in `<PMF, all data>` at each bin may be calculated as the maximum absolute difference between any pair of values from `<PMF, partition 1 of n> ... <PMF, partition n of n>`. If the PMF shifts in one specific direction as the simulation progresses, this might tell us that the system is not in equilibrium after all.

### Optional arguments

`-t_min` Minimum simulation time. Any data with t < t_min is skipped during loading of data.

`-t_max` Maximum simulation time. Any data with t > t_max is skipped during loading of data.

`-n_part` Number of partitions along time. PMFs are calulated for each of the partitions in addition to the PMF over all data. Default is n_part = 0.

`-xi_min` Minimum value of the coordinate (xi) to be considered. This specifies the lowest bin edge. If no value is specified, then the lowest bin edge is the minimum value of xi present in the data.

`-xi_max` Maximum value of the coordinate (xi) to be considered. This specifies the highest bin edge. If no value is specified, then the highest bin edge is the maximum value of xi present in the data.

`-n_bins` Number of bins along the coordinate (xi). Default is n_bins = 200.

`-tol` Tolerance for the WHAM algorithm, units of kcal/mol. WHAM iterates two equations (Roux Eq. 8 and 9) to determine a set of free energy constants $F_i$ that is required to align the PMFs from neighboring windows. The iterations are continued until the difference in $F_i$ to the previous iteration is below the tolerance. Default is tolerance = 0.0001 kcal/mol

`-kT` Thermal energy, units of kcal/mol. Deafult is $k_B T = 0.001987191 \text{ kcal } \text{ mol}^{-1} \text{ K}^{-1} \cdot 310 \text{K} = 0.61602921 \text{ kcal } \text{ mol}^{-1}$, where $T = 310 \text{K}$ is the approximate body temperature in humans. 

`-print_progress` If this flag is given, then the progress of data loading and WHAM iterations will be printed out. Otherwise, only the time interval used for each partition is printed out.

`-plot` If this flag is given, then a plot of all calculated PMFs is produced using `matplotlib.pyplot` after the calculation is finished.

## Example data
`windows_reus1_T310_N17_alpha.dat` and the `Data` directory contain example data taken from REUS (Replica Exchange Umbrella Sampling) MD (Molecular Dynamics) simulations of the N17 peptide (N-terminal 17 amino acids) from the huntingtin protein.

N17 associates with lipid membranes, and tends to fold into an alpha helix in their presence. To probe this folding, the RMSD (root mean square deviation) of N17 from an alpha helical reference was chosen as the coordinate $\xi$ (xi) along which to calculate the PMF. Low values of RMSD indicate a folded state, and high values of RMSD indicate an unfolded state.

The simulation data included here is from an umbrella sampling with 16 evenly-spaced window centers ranging from 0.9 &angst; to 6.15 &angst;. In each of the simulations, N17 is associated with a 5-component model of a mammalian lipid membrane, and its RMSD is harmonically restrained around the window center with a spring constant of 8.0 $\text{kcal } \text{ mol}^{-1}$ &angst; $^{-1}$. An integration time step of 4 fs was used, and the RMSD of N17 was recorded every 100 steps (0.4 ps). The data in this repository includes only every 100th recorded RMSD value (spacing of 40 ps between RMSD values). Each window's simulation was run for over 1 $\mu \text{s}$ (250 million steps).

The simulations were performed using NAMD 2.14 (Nanoscale Molecular Dynamics) together with the CHARMM36m all-hydrogen parameters for proteins, the CHARMM36 all-hydrogen parameters for lipids, and the CHARMM36 TIP3P water model were used. The simulations were run under NPT (constant number of atoms, constant pressure, constant temperature) at a temperature of 310 K (body temperature) and a pressure of 1.01325 bar (atmospheric pressure).

## References
[^fn1]: B. Roux, “The calculation of the potential of mean force using computer simulations,” Computer Physics Communications, vol. 91, no. 1, pp. 275–282, Sep. 1995, doi: 10.1016/0010-4655(95)00053-I.
