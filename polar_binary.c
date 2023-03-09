/* polar_binary.c */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "rebound.h"
#include "tree.h"
#include "output.h"
#include "communication_mpi.h"

typedef struct reb_simulation simulation;
typedef struct reb_treecell treenode;

void heartbeat(simulation* const);
void read_initial_conditions(int, char**);
void simulation_init(simulation*);
void initial_conditions(simulation*);


//const int NSTEPS = 1e1;

// Default initial conditions
const double TMAX = 8e2;
double MASS_RATIO = 0.0, A = 1.0, E = 0.0;
int N_INITIAL = 1e5;
double DISC_MASS = 2e-1;
//double INC_BASE = 0, OMEGA_BASE = 0.0;
//double DINC = M_PI * 1e-3, DOMEGA = M_PI * 1e-3;

char out_pt[512];
char* out_tree = "tree";
const double OUTPUT_INTERVAL = TMAX;
const double OUTPUT_TREE_INTERVAL = 50;
const double BOXSIZE = 10.0;
const double BOXN = 2;

int main(int argc, char* argv[]) {
	simulation* const s = reb_create_simulation();

	read_initial_conditions(argc, argv);
	simulation_init(s);

	if (s->mpi_id == 0) {
		initial_conditions(s);
		clean_files(out_pt, out_tree);
		write_header(out_pt, N_INITIAL, MASS_RATIO, A, E);
	}

	reb_communication_mpi_distribute_particles(s);
	move_to_com(s);
	MPI_Barrier(MPI_COMM_WORLD);

	reb_integrate(s, TMAX);

	MPI_Barrier(MPI_COMM_WORLD);
	reb_mpi_finalize(s);
	reb_free_simulation(s);
}

void heartbeat(simulation* const s) {
	move_to_com(s);

	if (reb_output_check(s, 10*s->dt))
		reb_output_timing(s, 0);

	if (reb_output_check(s, OUTPUT_INTERVAL))
		write_sim(s, out_pt, MASS_RATIO);

	if (s->t > 500 && reb_output_check(s, OUTPUT_TREE_INTERVAL))
		write_tree(s, out_tree, BOXN);

	if (s->t == TMAX && s->mpi_id == 0)
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
	s->heartbeat = heartbeat;

	reb_configure_box(s, BOXSIZE, BOXN,BOXN,1);
	reb_mpi_init(s);
}

void initialize_binary(simulation* s) {
	/* binary */
/*
	struct reb_particle star1 = {0}; star1.m = 1.0;
	struct reb_particle star2 = reb_tools_orbit_to_particle(s->G, star1, star1.m*MASS_RATIO,
									     A,E,0,0,M_PI/2,0);

	star1.hash = reb_hash("Star1");
	star2.hash = reb_hash("Star2");

	reb_add(s, star1);
	reb_add(s, star2);

	reb_move_to_com(s);
*/

	/* unary */
	struct reb_particle star = {0}; star.m = 1.0;
	star.hash = reb_hash("Star1");
	reb_add(s, star);
	reb_move_to_com(s);
}

void initialize_disc(simulation* s) {
//	double inner_a = 2 * (1+E)*A, outer_a = 5 * (1+E)*A;
	double inner_a = 1.02, outer_a = 10.2/2./1.2;
	double m = DISC_MASS/N_INITIAL;

//	struct reb_particle com = reb_get_com(s);
//	struct reb_particle com = {0}; com.m = 1+MASS_RATIO; com.x = 0; com.y = 0; com.z = 0;

	for (int i = 0; i < N_INITIAL; i++) {
/*
		double a = reb_random_uniform(s, inner_a, outer_a);
		double nu = reb_random_uniform(s, 0, 0.5*M_PI);
		double inc = reb_random_uniform(s, INC_BASE-DINC, INC_BASE+DINC);
		double Omega = reb_random_uniform(s, OMEGA_BASE-DOMEGA, OMEGA_BASE+DOMEGA);

		struct reb_particle pt = reb_tools_orbit_to_particle(s->G, com, m,a,0,inc,Omega,0,nu);
*/

		struct reb_particle pt = {0};
		double a = reb_random_powerlaw(s, inner_a, outer_a, -1.5);
		double phi = reb_random_uniform(s, 0, 2*M_PI);
		double mu = (1.0 + MASS_RATIO) + DISC_MASS * (pow(a, -3./2.) - pow(inner_a, -3./2.)) / (pow(outer_a, -3./2.) - pow(inner_a, -3./2.));
		double vkep = sqrt(s->G*mu/a);

		pt.x = a*cos(phi); pt.y = a*sin(phi);
		pt.z = a*reb_random_normal(s, 0.001);
		pt.vx = vkep*sin(phi); pt.vy = -vkep*cos(phi); pt.vz = 0;
		pt.m = m;

/*
		double a = inner_a + i * (outer_a-inner_a)/N_INITIAL;
		double nu = i * 2*M_PI/N_INITIAL;
		double inc = INC_BASE;
		double Omega = OMEGA_BASE;

		struct reb_particle pt = reb_tools_orbit_to_particle(s->G, com, m,a,0,inc,Omega,0,nu);
*/
/*
		struct reb_particle pt = {0};
		double a = inner_a + i * (outer_a-inner_a)/N_INITIAL;
		double phi = i * 2*M_PI/N_INITIAL;
		double mu = 1 + MASS_RATIO;
		double vkep = sqrt(s->G*mu/a);

		pt.x = a*cos(phi); pt.y = a*sin(phi);
		pt.z = 0;
		pt.vx = vkep*sin(phi); pt.vy = -vkep*cos(phi); pt.vz = 0;
		pt.m = m;
*/

		char hash_str[16]; sprintf(hash_str, "%d", i);
		pt.hash = reb_hash(hash_str);

		reb_add(s, pt);
	}
}

void initial_conditions(simulation* s) {
	initialize_binary(s);
	initialize_disc(s);
}

void read_initial_conditions(int argc, char* argv[]) {
	sprintf(out_pt, "%s", "particles.txt");

	for (int i = 1; i < argc; i++) {
		if (strcmp(argv[i], "-o") == 0)
			sprintf(out_pt, "%s", argv[i+1]);
		if (strcmp(argv[i], "-mr") == 0)
			MASS_RATIO = strtod(argv[i+1], NULL);
		if (strcmp(argv[i], "-e") == 0)
			E = strtod(argv[i+1], NULL);
	}
}
