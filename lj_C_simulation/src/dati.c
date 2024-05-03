#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <omp.h>
#include "dati.h"

double** allocate_double_matrix(int rows, int cols) {
    /* Function to allocate a matrix of doubles of size (rows, cols) */

    double** matrix = (double**) malloc(rows * sizeof(double*));

    for (int i=0; i<rows; i++) {
        matrix[i] = (double*) malloc(cols * sizeof(double));
    }

    return matrix;
}



void free_matrix(double** matrix, int rows){
    /* Function to free a matrix */

    for (int i = 0; i < rows; i++) {
        if (matrix[i] != NULL ) {
            free(matrix[i]);
        } else {
            printf("Error: you are trying to free a location which is not initialized in memory!");
            exit(1);
        }
    }
    
	free(matrix);
}



void generate_FCC(double L, double** x)
{
  /*
    This function fills the matrix x[N][3] with the coordinates of
    particles in a FCC lattice within a box of side L whose center is the
    center of the coordinate system.
    The box is filled uniformly if N = 4 n^3 with n integer hence for
    N = 4     32    108    256    500    864   1372   2048   2916
    4000   5324 6912   8788  10976  13500  16384 ...
	Thanks to Dott. Giovanni Garberoglio for this function.
  */

  int i,j,k,m,n,p,c;
  
  /* position within the primary cell */
  const double rFCC[4][3] = {{0.0, 0.0, 0.0}, {0.0, 0.5, 0.5},
                             {0.5, 0.0, 0.5}, {0.5, 0.5, 0.0}};
  double b, rCell[3];
  
  for (c = 1; ; c++)
    if (4*c*c*c >= N)
      break;

  b = L / (double)c;            /* side length of the primary cell */
  p = 0;                        /* particles placed so far */
  for (i = 0; i < c; i++)
    {
      rCell[0] = i;
      for (j = 0; j < c; j++)
        {
          rCell[1] = j;
          for (k = 0; k < c; k++)
            {
              rCell[2] = k;
              for (m = 0; m < 4; m++) /* 4 particles in cell */
                if (p < N)
                  {
                    
                    /* add the com to each bead, and project to the real cell */
                    for(n=0;n<3;n++)
                      {
                        x[p][n] = b * (rCell[n] + rFCC[m][n]);
                        x[p][n] -= L * rint(x[p][n]/L);
                      }
                ++p;
                  }
            }
        }
    }  
}



void calculate_forces(double** x, double** a, double L){

	/* Initialize acceleration to zero */
	for (int i = 0; i < N; i++){
		for (int k = 0; k < 3; k++){
			a[i][k]=0;
		}
	}

	/* Calculate forces in parallel */
	#pragma omp parallel
	{
		#pragma omp for
		for (int i = 0; i < N; i++){
			double dx[3], r_2;
			for (int j = 0; j < N; j++){
				for(int k = 0; k < 3; k++){
					dx[k] = x[i][k]-x[j][k];				/* Distance between particles i and j */
					dx[k]= dx[k] - L*rint(dx[k]/L);			/* Minimum image convention */
				}
				/* Calculate forces only within the cutoff */
				r_2 = dx[0]*dx[0] + dx[1]*dx[1] + dx[2]*dx[2];
				if ( sqrt(r_2) <= cutoff && i!=j ) {
					for (int k = 0; k < 3; k++){
						a[i][k] += 24 * dx[k] * ( 2.0 * pow(r_2, -7) - pow(r_2, -4) );
					}
				}
			}
		}
	}

	/* Wait for every core to finish */
	#pragma omp barrier
}



void velocity_verlet(double** x, double** v, double** a, double L){
	/* Velocity-verlet algorithm */

	double v_dt_2[N][3];					/* Intermediate velocities */

	for (int i = 0; i < N; i++){
		for (int k = 0; k < 3; k++){
			/* First update of the velocities */
			v_dt_2[i][k]= v[i][k]+ a[i][k]*dt*0.5;
			/* Update positions */
			x[i][k]= x[i][k] +v_dt_2[i][k]*dt;
		}
	}
	/* Update forces */
	calculate_forces(x, a, L);
	/* Second update of the velocities */
	for (int i = 0; i < N; i++){
		for(int k = 0;k < 3; k++){
			v[i][k] = v_dt_2[i][k] + a[i][k]*dt*0.5;
		}
	}
}



double kinetic_energy(double** v){
	/* Function to calculate the average kinetic energy */

	double K = 0;
	for (int i = 0; i < N; i++){
		for (int k = 0; k < 3; k++){
			K += v[i][k] * v[i][k];
		}
	}
	K *= 0.5;
	K /= N;
	return K;
}



void rescale_velocities(double** v){
	/* Berendsen thermostat */

	int i,k;
	double alfa, K_i, T_i;
	/* Calculate rescaling factor */
	K_i = kinetic_energy(v);
	T_i = 2.0 * K_i / 3;
	alfa = sqrt(T/T_i);
	/* Rescale velocities */
	for (i = 0; i < N; i++){
		for (k = 0; k < 3; k++) {
			v[i][k] = alfa*v[i][k];
		}
	}
}



double v_lenard_jones(double** x, double L) {
	/* Function to calculate the LJ potential */
	double U_tot, dx[3], r_2;
	U_tot = 0;

	for (int i = 0; i < N; i++){
		for (int j = 0; j < i; j++){
			for (int k = 0; k < 3; k++){
				dx[k] = x[i][k]-x[j][k];				/* Distance between particles i and j */
				dx[k]= dx[k] - L*rint(dx[k]/L);			/* Minimum image convention */
			}
			r_2 = dx[0]*dx[0] + dx[1]*dx[1] + dx[2]*dx[2];
			if (sqrt(r_2) <= cutoff){
				U_tot += 4 *( pow(r_2, -6) - pow(r_2, -3) );
			}
		}
	}

	// I calculate the average potential energy per particle
	U_tot = U_tot/N;
	// Correction for the finite cutoff
	U_tot -= 8*(3.1415)*rho/(3*pow(cutoff,3));

	return U_tot;
}



void eval_g(double* g, double** x, double L){
	/* Evaluate histogram measuring the distances between all the couples of particles */

	int i, j, k, bin;
	double dx[3], dr_ij;
	double dr = L/(2*S);				/* dr is the width of a bin of the g histogram, with S bins between 0 and L/2 */

	for (i = 0; i < N; i++){
		for (j = 0; j < N; j++){
			if (i != j){
				for (k = 0 ;k < 3; k++){
					dx[k]=x[i][k]-x[j][k];
					dx[k]=dx[k]- L*rint(dx[k]/L);
				}
				dr_ij = sqrt(dx[0]*dx[0] + dx[1]*dx[1] + dx[2]*dx[2]);
				if (dr_ij < L/2) {
				    bin = floor(dr_ij / dr);		/* The bin is obtained from the distance measured in units of r */
				    g[bin] += 1.0 / dr;				/* I add +1, normalizing with the bin size dr */
				}
			}
		}
	}
}



void check_file_existence(FILE *pf) {
	if (pf == NULL) {
		fprintf(stderr, "\nError: a file was not correctly found.\nPlease check the paths you used.\n");
		exit(EXIT_FAILURE);
	}
}