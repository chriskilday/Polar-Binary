/* polar_binary.c */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <mpi.h>
#include "rebound.h"
#include "tree.h"
#include "output.h"

typedef struct reb_simulation simulation;
typedef struct reb_treecell treenode;

void heartbeat(simulation*);
void simulation_init(simulation*);
void initial_conditions(simulation*);
void write(simulation*);

const double TMAX = 1e-3;

char* out_tree = "tree";
const double OUTPUT_INTERVAL = 100;
const double BOXSIZE = 8.0;
const double DISC_MASS = 1e-2;

int main(int argc, char* argv[]) {
	simulation* const s = reb_create_simulation();

	simulation_init(s);
	reb_mpi_init(s);
//	initial_conditions(s);
//	if (s->mpi_id == 0) clean_files(out_tree);
//	MPI_Barrier(MPI_COMM_WORLD);

	reb_integrate(s, TMAX);

//	write_tree(s, out_tree);
	reb_mpi_finalize(s);
	reb_free_simulation(s);

	return 0;
}

void heartbeat(simulation* s) {
	reb_output_timing(s, 0);
	if (s->t >= TMAX)
		printf("\n");
}

void simulation_init(simulation* s) {
	s->integrator = REB_INTEGRATOR_LEAPFROG;
	s->gravity = REB_GRAVITY_TREE;
	s->boundary = REB_BOUNDARY_OPEN;
	s->opening_angle2 = 1.5;
	s->G = 1;
	s->softening = 0.02;
	s->dt = 3e-2;
//	s->heartbeat = heartbeat;

	reb_configure_box(s, 10.2, 2,2,1);
}

void initial_conditions(simulation* s) {
	struct reb_particle star = {0};
	star.m = 1.0;
	if (s->mpi_id == 0) reb_add(s, star);

	int n = 100;
	for (int i = 0; i < n; i++) {
		struct reb_particle pt = {0};
		double a = reb_random_powerlaw(s,BOXSIZE/10,BOXSIZE/2/1.2,-1.5);
		double phi = reb_random_uniform(s, 0.0, 2 * M_PI);
		double mu = star.m + DISC_MASS * (pow(a,-3./2.)-pow(BOXSIZE/10,-3./2.))/(pow(BOXSIZE/2./1.2,-3./2.)-pow(BOXSIZE/10.,-3./2.));
		double vkep = sqrt(s->G*mu/a);

		pt.x = a*cos(phi); pt.y = a*sin(phi);
		pt.z = a*reb_random_normal(s, 0.001);
		pt.vx = vkep*sin(phi); pt.vy = -vkep*cos(phi); pt.vz = 0;
		pt.m = DISC_MASS / (double)n;

		reb_add(s, pt);
	}
}

void write(simulation* s) {
	printf("N: %d\n", s->N);
}
