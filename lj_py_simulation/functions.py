""" This module contains all the functions for the LJ fluid molecular dynamics script """

# Importing modules
import numpy as np
from numba import njit, prange


# Functions

def generate_positions (x, L) :
    """ This function generates the initial positions for the particles on a FCC
        NOTE: sorry for the awful code, the function was copied from the C version, as it can be seen from the code """
    # Basis vectors in the FCC unit cell
    r_FCC = np.array( [ [0.0, 0.0, 0.0], [0.0, 0.5, 0.5], [0.5, 0.0, 0.5], [0.5, 0.5, 0] ] )
    r_cell = np.zeros (3)           # Position of the unit cell in the cubic lattice 
    for c in range (1, 2*x.shape[0]):
        if (4 * c**3 >= x.shape[0]):
            break
    b = L / c               # Side length of the unit cell
    p = 0                   # Number of particles added
    # i,j,k : indeces of the cubic Bravais lattice
    # m,n : indeces of the unit cell
    for i in range (0, c):
        r_cell[0] = i
        for j in range (0, c):
            r_cell[1] = j
            for k in range (0, c):
                r_cell[2] = k
                for m in range (0,4):
                    if p < x.shape[0]:
                        x[p] = b * (r_cell + r_FCC[m])
                        p += 1



@njit(parallel = True)
def calculate_forces (x, cutoff, L) :
    """ This function calculates the forces acting on each particle and store them in the matrix a """

    a = np.zeros(x.shape)                       # I put a to zero, makes debugging easier
    dx = np.zeros(x.shape)                      # Matrix of the displacements
    r_2 = np.zeros(x.shape[0])                  # Matrix of the distances
    # I evaluate forces for each particle
    for i in prange(0, x.shape[0]) :
        dx = x - x[i]
        dx = dx - L * np.rint( dx/L )           # Minimum image convention
        r_2 = np.square(dx).sum(axis=1)         # Square distance between particle i and the other particles
        r_2_bool = (np.sqrt(r_2) < cutoff) 
        r_2_bool[i] = False                     # No self-interaction, otherwise the forces blow up
        r_2_eff = r_2[r_2_bool]                 # I take just particles within cutoff radius
        r_2_eff = r_2_eff.reshape(r_2_eff.size, 1)
        dx_eff = dx[r_2_bool]
        a[i] += (np.multiply(-24 * dx_eff, 2 * np.power(r_2_eff, -7) - np.power(r_2_eff, -4) ) ).sum(axis = 0)
    
    return a



@njit
def velocity_verlet (x, v, a, dt, cutoff, L) :
    """ This function upgrades positions and velocities of the particles using the velocity
        verlet algorithm. This function makes use of calculate_forces """
    v_dt_2 = np.zeros(v.shape)
    v_dt_2 = v + 0.5 * a * dt                       # Update velocities pt.1
    x = x + v_dt_2 * dt                             # Update position
    a = calculate_forces(x, cutoff, L)              # Update forces
    v = v_dt_2 + 0.5 * a * dt                       # Update velocities pt.2
    return x, v, a



@njit
def kinetic_energy (v) :
    """ This function evaluates the kinetic energy """
    return (np.square(v)).sum() * 0.5



@njit
def rescale_velocity (v, T) :
    """ This function is a thermostat: it rescales the velocity to guarantee that the 
        kinetic energy has a value compatible with the temperature T """
    K = kinetic_energy(v)                               # Kinetic energy
    T_i = 2 * K / (3*v.shape[0])                        # "Istantaneous" temperature
    alpha = np.sqrt( T/T_i )                            # Rescaling factor
    return v * alpha



@njit (parallel = True)
def V_LJ (x, cutoff, L, rho) :
    """ This function evaluates the value of the LJ potential with a given cutoff,
        applying a correction to take the excess energy into account """
    dx = np.zeros(x.shape)             # Matrix of the distances
    U = 0.0                            # Potential to be calculated

    for i in prange(0, x.shape[0]) :
        dx = x - x[i]
        dx = dx - L * np.rint(dx/L)         # Minimum image convention
        r_2 = np.square(dx).sum(axis=1)
        r_2_bool = (np.sqrt(r_2) < cutoff)
        r_2_bool[i] = False                 # No self-interaction, otherwise the forces blow up
        r_2_eff = r_2[r_2_bool]
        U += ( 4 * ( np.power(r_2_eff, -6) - np.power(r_2_eff, -3) ) ).sum()

    U = U/2                                                                 # Avoid double-counting
    U -= x.shape[0] * 8 * 3.1415 * rho / (3 * np.power (cutoff,3) )         # Cutoff excess energy correction
    return U/x.shape[0]                                                     # I return potential energy per particle



def eval_g (x, L, dr, S) :
    """ This function evaluates the radial distribution function """
    dx = np.zeros(x.shape)             					        # Matrix of the displacements
    r_2 = np.zeros(x.shape)            					        # Matrix of the distances
    hist = np.zeros(S)

    for i in prange(0, x.shape[0]) :
        dx = x - x[i]
        dx = dx - L * np.rint(dx/L)        				        # Minimum image convention
        r_2 = np.sqrt( np.square(dx).sum(axis=1) )
        r_2_bool = ( r_2 < L/2 ); r_2_bool[i] = False
        r_2_eff = r_2[r_2_bool]
        r_2_eff = r_2_eff / dr                      			# Distances in units of dr
        g_bin = np.floor(r_2_eff).astype(int)                   # Bin in which the distance is included
        tmp_1, _ = np.histogram(g_bin, np.arange(S+1))          # Histogram for particle i
        hist = hist + tmp_1   

    return hist
