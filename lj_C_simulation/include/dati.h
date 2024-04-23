#ifndef __DATI_H__
#define __DATI_H__


#define S 100           /* Bins for radial distribution function calculation */


double** allocate_double_matrix(int rows, int cols);

void free_matrix(double** matrix, int rows, int cols);

void generate_FCC(unsigned int N, double L, double** x);

void calculate_forces (unsigned int N, double** x, double** a, double cutoff, double L);

void velocity_verlet(unsigned int N, double** x, double** v, double** a, double dt, double cutoff, double L);

double kinetic_energy(unsigned int N, double** v);

void rescale_velocities(unsigned int N, double** v, double T);

double v_lenard_jones (unsigned int N, double** x, double L, double cutoff, double rho);

void eval_g (unsigned int N, double* g, double** x, double L);

void check_file_existence(FILE *pf);


#endif

