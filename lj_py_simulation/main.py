# Importing modules
import numpy as np
import matplotlib.pyplot as plt
from functions import *
import sys
import time
from multiprocessing import cpu_count

# Variables 
rho = 0.776                     # System density                
T = 0.85                        # System temperature
N = 256                         # Number of particles
x = np.zeros((N,3))             # Positions
v = np.zeros((N,3))             # Velocities
a = np.zeros((N,3))             # Accelerations (i.e. forces, since m = 1 in our units)
L = (N / rho)**(1/3)            # Box length
cutoff = 3                      # Cutoff for LJ interaction
dt = 0.003                      # Time step
t_count = 1                     # Counter of time steps
t_rescale = 20                  # Number of timesteps after which the thermostat acts
t_thermostat = 2000             # Number of timesteps for the thermostatting phases
t_run = 16000                   # Number of steps for production run
t_measure = 160                 # Number of timesteps after which I measure stuff
num_measure = 0                 # Number of measures
S = 100                         # Number of bins for the discretization of space in g(r) calculation
dr = L / (2*S)                  # Bin width for g(r)
g = np.zeros(S)                 # Radial distribution function
U_mean = 0                      # Average potential energy per particle
K_mean = 0                      # Average total kinetic energy
ncpu = cpu_count()              # Number of cpu available in the machine
print("\nThe number of employed threads is ", ncpu, "\n")

# Measure time
start = time.perf_counter()

# Position initialization
generate_positions(x, L)
print("FCC lattice : created")

# Calculation of initial forces
a = calculate_forces(x, cutoff, L)

# Status message
print( f"\nStarting thermostatting at temperature {T}:" )

# NVT equilibration with Berendsen thermostat
U = []                          # Potential energy list
t = []                          # List for x axis

while (t_count <= t_thermostat) :
    # Measure potential
    t.append(t_count * dt)
    U.append( V_LJ(x, cutoff, L, rho) )
    # Time integration
    x, v, a = velocity_verlet(x, v, a, dt, cutoff, L)
    # Thermostat
    if (t_count % t_rescale == 0) :
        v = rescale_velocity(v, T)
    # Increase time
    t_count += 1

print (f"\nThermostatting phase : finished")

# I check the kinetic energy before going on
trial = 2 * kinetic_energy(v) / (3 * v.shape[0])
print(f"Actual temperature : {trial:.3f}")
print(""); print (50 * "~"); print ("")
print("Starting the NVE simulation.\n")

# Production run without thermostat
t_count = 0                               # Reset counter
while (t_count <= t_run) :

    # Time integration
    x, v, a = velocity_verlet(x, v, a, dt, cutoff, L)
    
    # Measure
    if ( t_count % t_measure == 0 ) :
        U_mean += V_LJ(x, cutoff, L, rho)
        K_mean += kinetic_energy(v)
        # I measure the histogram for all particles and I sum it to my g vector
        hist = eval_g(x, L, dr, S)
        g = g + hist
        num_measure += 1
    
    # Increase time
    t_count += 1
    
    # Progress bar
    actual_prog = -1
    prog = np.floor(100*t_count/t_run).astype(int)
    idx = 1
    max_ast = np.rint(100/idx).astype(int)              # The maximum number of asterisks
    num_ast = np.floor(prog/idx).astype(int)            # Actual number of asterisks
    if prog != actual_prog :
        sys.stdout.write("\r" + "Calculating" + "\t" + "[" + "#" * num_ast + "." * (max_ast -num_ast) + "]" + "\t" + prog.astype(str) + "%" )
        actual_prog = prog

# Calculating mean energies, not total energies
U_mean /= num_measure
K_mean /= num_measure

# Results output
print("\n\n\nResults:")
print(f"<U> = {U_mean:.3f}")
print(f"<T> = {2*K_mean / (3*N):.3f}")
print("Number of measures = ", num_measure, "\n")
print(f"Initial parameters :\n\t-> T = {T:.3f}\n\t-> rho = {rho:.3f}\n\t-> N = {N}\n\t-> L = {L:.3f}")

# End time measurement
end = time.perf_counter()
print(f"\nExecution time = {end-start:.3f} s\n")

# Plot of the potential
fig, ax = plt.subplots()
ax.plot(t, U)
ax.set_xlabel('t')
ax.set_ylabel('U')
plt.show()

# Plot of the g(r)
g = g / dr                                              # Normalization to bin width
r = (np.arange(S) + 0.5) * dr                           # x axis values for distance
g = g / ( num_measure*N*4*3.1415*rho*np.square(r) )     # Normalization
fig, ax = plt.subplots()
ax.plot(r, g)
ax.set_xlabel("r")
ax.set_ylabel("g")
plt.show()
