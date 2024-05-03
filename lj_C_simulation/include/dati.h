#ifndef __DATI_H__
#define __DATI_H__

/* Global variables */

extern const double rho;
extern const double T;
extern const unsigned int N;
extern const double cutoff;
extern const double dt;


/* Symbolic constants */

#define S 100


/* Headers */

double** allocate_double_matrix(int rows, int cols);

void free_matrix(double** matrix, int rows);

void generate_FCC(double L, double** x);

void calculate_forces (double** x, double** a, double L);

void velocity_verlet(double** x, double** v, double** a, double L);

double kinetic_energy(double** v);

void rescale_velocities(double** v);

double v_lenard_jones (double** x, double L);

void eval_g (double* g, double** x, double L);

void check_file_existence(FILE *pf);


#endif

