# Lennard_Jones_fluid

This repo contains the necessary files to perform the simulation of a Lennard-Jones fluid in both the C and python programming languages.

## Physics background

The system is composed by $N$ identical particles in a box of volume $V$, whose value is decided by fixing the density of the system $\rho$.

Two generic particles $i$ and $j$ interact via the 6-12 Lennard-Jones potential, namely:

$$
V_{LJ}(r_{ij}) = 4 \left( \frac{1}{r_{ij}^{12}} - \frac{1}{r_{ij}^6} \right)
$$

where $r_{ij}$ is the modulus of the distance between the particles. All the distances are measured in units of $\sigma$ and all the energies in units of $\epsilon$, being $\sigma$ and $\epsilon$ the LJ parameters of the chosen chemical element.

Please note that, when $r_{ij}$ is greater than a certain cutoff, the particles are considered non-interacting (since the LJ potential goes to zero quite fast). A correction is introduced in the potential energy to take this into account.

To prevent finite-size effects in the system, the minimum image convention is employed.

The simulation is divided into the following steps:

1. *System creation*: the particles are initially disposed in an FCC crystal structure within the box, to avoid dangerous contacts, and with velocity equal to 0.

2. *Thermalization*: the system is thermalized by periodically rescaling the velocities in order to achieve the desired value of temperature (Berendsen thermostat). At the end of this phase, the potential must oscillate around a constant value.

3. *Production run*: the last step consists of an NVE production run

During the last phase, the 2-body radial distribution function $g(r)$ of the system is calculated.

## Technical details: C version

The scripts in the `lj_C_simulation` directory are written in C. All the necessary libraries should already be installed on your system.

To compile and run the simulation (only gcc is required), open the directory containing this file in the terminal and type the following commands:

```bash
cd lj_C_simulation
make
cd build
./exe
```

At the end of the execution, in the build directory you will find two additional files, namely `potential.txt` and `g.txt`: the first contains the potential in function of time during the equilibration; the second contains the 2-body radial distribution function in function of the distance.  
You can plot them with your favourite tool; if you have `gnuplot` installed, a simple way to visualize the results is to open the build directory in the terminal and type:

```bash
gnuplot
plot 'potential.txt' w l lt 1 lw 2
plot 'g.txt' w l lt 1 lw 2
```

Use `quit` to close `gnuplot`.

> **Warning**: these commands have been tested on Ubuntu 22.04, but they should work on any (not too bizarre) linux system. Nothing is guaranteed on Windows.

## Technical details: python version

The scripts in the `lj_py_simulation` directory are written in python. The only additional packages required to run the simulation are `numpy`, `numba` and `matplotlib`; if you have conda installed, the command:

```bash
conda create -n lj numba matplotlib
```

should create an environment with all the necessary packages.

After this, you can start the simulations with the commands:

```bash
conda activate lj
cd lj_py_simulation
python main.py
```

## Future improvements

* Calculate total energy during the production run
* Print the trajectory in a VMD-compatible format
* Add the neighbour list

## Contact

If you have any doubt or if you spot mistakes/bugs, feel free to contact [fabio.passi24@gmail.com](fabio.passi24@gmail.com)
