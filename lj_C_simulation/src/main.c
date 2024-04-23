#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
#include <math.h>
#include "dati.h"

/* Units:
	-> mass = 1
	-> sigma = 1
	-> epsilon = 1
*/

int main(int argc, char *argv[]) {

	/* VARIABLES */
	int i,k;
	double rho= 0.776; 									/* Density */
    double T = 0.85; 									/* Temperature */
    unsigned int N = 256; 								/* Number of particles*/
	double** x = allocate_double_matrix(N,3); 	 		/* Positions */
	double** v = allocate_double_matrix(N,3);			/* Velocities */
	double** a = allocate_double_matrix(N,3);			/* Accelerations (= forces) */
	double cutoff = 3; 									/* Cutoff for the LJ interaction */
	double dt = 0.003; 									/* Time-step */
	int t = 1;	 										/* Number of time steps elapsed */
	int t_thermalize = 2000;							/* Number of steps for thermalization phase */
	int t_simulation = 16000;							/* Number of steps for the production run */
	int t_thermostat = 20;								/* Number of steps every which I apply the thermostat */
	int measure = 160; 									/* I perform a measure every "measure" steps */
	double potential; 									/* Variable used to store the potential value */
	double U_mean, K_mean; U_mean=0; K_mean=0; 			/* Average values of energues */
	int num_measure = 0;								/* Number of measures */
	double* g = malloc(S * sizeof(*g));					/* Radial distribution function */
	double L = pow ((float)N/rho , 1.0/3);				/* Box length */

	/* Progress bar stuff */
	int actual_prog = -1;
	int idx = 1;											
	int max_ast = rint(100*idx);			             /* The maximum number of asterisks */
	int prog, num_ast;


	/* SIMULATION */

	/* Start measuring time */
    struct timeval begin, end;
    gettimeofday(&begin, 0);

	/* I initialize to zero the velocities and the vector g */
	for (i=0; i<S; i++){
		g[i]=0;
	}

	for (i=0; i<N; i++){
		for (k=0;k<3;k++){
			v[i][k]=0;
		}
	}

	/* Generate initial positions */
	printf("\nGenerating initial position on FCC lattice.\n\n");
	generate_FCC (N, L, x);

	/* Evaluate initial forces */
	calculate_forces (N, x, a, cutoff, L);

	FILE *equil_data;
	equil_data = fopen ("potential.txt","w");
	check_file_existence(equil_data);
	fprintf (equil_data, "# Time\tPotential\n");

	/* Status message */
	printf("Starting thermostatting at temperature %.3lf.\n\n", T);

	/* Thermalization phase (NVT) */
	while (t <= t_thermalize) {

		/* Evaluate potential and print to file */
		potential = v_lenard_jones (N, x,L, cutoff,rho);
		fprintf (equil_data, "%lf	%lf\n", t * dt, potential);
		/* Time evolution */
		velocity_verlet (N, x, v, a, dt, cutoff, L);
		/* Thermostat */
		if (t % t_thermostat == 0) {
			rescale_velocities (N, v, T);
		}

		++t;
	}

	fclose (equil_data);

	/* Kinetic energy check after thermalization */
	double K = kinetic_energy(N, v);

	/* Status message */
	printf ("Thermostatting phase: finished.\nActual temperature: %.3f\n\n", 2*K/3);
	printf("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n\n");
    printf ("Starting the NVE simulation.\n\n");
	
    /* NVE production run and data acquisition */
    t=0;
	while (t <= t_simulation){

		/* Time evolution */
		velocity_verlet (N, x, v, a, dt, cutoff, L);

		/* Measure */
		if (t%measure == 0) {
			U_mean += v_lenard_jones (N, x, L, cutoff, rho);
			K_mean += kinetic_energy(N, v);
			eval_g (N, g, x, L);
			++num_measure;
		}

		// Progress bar
		prog = (int) (100.0* t / t_simulation);
		num_ast = prog/idx;					         		// Actual number of asterisks
		if (prog != actual_prog) {
			printf("\rCalculating\t[ ");
			for (i=0; i<num_ast; i++) {printf("#");}
			for (i=num_ast; i<max_ast; i++) {printf(".");}
			printf(" ]  %d%%", prog);
			actual_prog = prog;
			fflush(stdout);
		}

		t++;
	}
	printf("\n\n");

	/* I evaluate the radial distribution function from the histogram of the counts */
	double r[S], dr;														/* r contains the position vectors in which I evaluate g */
	dr = L/(2*S);
	for (i = 0; i < S; i++){
		r[i]= ( (2*i+1.0) / 2 ) * dr;										/* The positions are the "middle of the bins" */
		g[i]= g[i] / ( num_measure * N * 4 * 3.1415 * rho * r[i] * r[i] );	/* Normalization factor */
	}

	/* Print g(r) on file */
	FILE *g_file;
	g_file = fopen ("g.txt","w");
	check_file_existence(g_file);
	fprintf (g_file, "#r\tg(r)\n");

	for (i = 0; i < S; i++){
		fprintf (g_file, "%lf	%lf\n", r[i], g[i]);
	}

	fclose (g_file);

	/* Evaluation of average energies */
	U_mean /= num_measure;
	K_mean /= num_measure;

	/* Status message */
	printf("\nResults:\n");
	printf ("<U> = %.3lf \n", U_mean);
	printf ("<T> = %.3lf \n", 2*K_mean/3);
	printf ("Number of measures = %d \n", num_measure);
	printf ("\nInitial parameters:\n\t-> T= %.3lf\n\t-> rho= %.3lf\n\t-> N= %d\n\t-> L= %.3lf\n\n",T, rho, N, L);
	
	// Stop measuring time and calculate the elapsed time
    gettimeofday(&end, 0);
    long seconds = end.tv_sec - begin.tv_sec;
    long microseconds = end.tv_usec - begin.tv_usec;
    double elapsed = seconds + microseconds*1e-6;
	printf("Execution time = %.3f s\n\n", elapsed);

	// Free memory
	free_matrix(x, N, 3);
	free_matrix(v, N, 3);
	free_matrix(a, N, 3);
	free(g);

	return 0;
}
