/**
 * Self-gravitating disc with MPI
 *
 * A self-gravitating disc is integrated using
 * the leap frog integrator. Collisions are not resolved.
 * This program makes use of MPI. Note that you need 
 * to have MPI compilers (mpicc) installed. The code is using 
 * four root boxes to distribute to particles to one, two 
 * or four MPI nodes. How to efficiently run this code on
 * large clusters goes beyond this simple example and
 * almost certainly requires experimentation.
 */
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include <time.h>
#include <sys/time.h>
#include "rebound.h"
#include "tools.h"

#include "output.h"

int N_INITIAL = 1e4;
double DISC_MASS = 2e-1;
const int TMAX = 1e2;

int OUT_INTERVAL = 10;
int OUT_TREE_N = 0;

void heartbeat(struct reb_simulation* const r);

int main(int argc, char* argv[]){
    struct reb_simulation* const r = reb_create_simulation();
    // Setup constants
    r->integrator    = REB_INTEGRATOR_LEAPFROG;
    r->gravity    = REB_GRAVITY_TREE;
    r->boundary    = REB_BOUNDARY_OPEN;
    r->opening_angle2    = 1.5;        // This constant determines the accuracy of the tree code gravity estimate.
    r->G         = 1;        
    r->softening     = 0.02;        // Gravitational softening length
    r->dt         = 3e-2;        // Timestep
    const double boxsize = 10.2;
    // Setup root boxes for gravity tree.
    // Here, we use 2x2=4 root boxes (each with length 'boxsize')
    // This allows you to use up to 4 MPI nodes.
    reb_configure_box(r,boxsize,2,2,1);
	r->collisions = REB_COLLISION_TREE;

    // Initialize MPI
    // This can only be done after reb_configure_box.
    reb_mpi_init(r);

	if (r->mpi_id==0) {
		clean_files("particles.txt", "tree");
		write_header("particles.txt", N_INITIAL, 0, 1, 0);
	}


    // Setup particles only on master node
    // In the first timestep, the master node will 
    // distribute particles to other nodes. 
    // Note that this is not the most efficient method
    // for very large particle numbers.
    double disc_mass = DISC_MASS/r->mpi_num;    // Total disc mass
    int N = N_INITIAL/r->mpi_num;            // Number of particles
    // Initial conditions
    struct reb_particle star = {0};
    star.m         = 1;
	star.hash = reb_hash("Star1");
    if (r->mpi_id==0){
        reb_add(r, star);
    }
	int start_i = r->mpi_id*N;
    for (int i=start_i;i<start_i+N;i++){
        struct reb_particle pt = {0};
        double a    = reb_random_powerlaw(r, boxsize/10.,boxsize/2./1.2,-1.5);
        double phi     = reb_random_uniform(r, 0,2.*M_PI);
        pt.x         = a*cos(phi);
        pt.y         = a*sin(phi);
        pt.z         = a*reb_random_normal(r, 0.001);
        double mu     = star.m + disc_mass * (pow(a,-3./2.)-pow(boxsize/10.,-3./2.))/(pow(boxsize/2./1.2,-3./2.)-pow(boxsize/10.,-3./2.));
        double vkep     = sqrt(r->G*mu/a);
        pt.vx         =  vkep * sin(phi);
        pt.vy         = -vkep * cos(phi);
        pt.vz         = 0;
        pt.m         = disc_mass/(double)N;
	pt.r = 4e-5;

	char hash_str[16]; sprintf(hash_str, "%d", i);
	pt.hash = reb_hash(hash_str);
        reb_add(r, pt);
    }
    r->heartbeat = heartbeat;

#ifdef OPENGL
    // Hack to artificially increase particle array.
    // This cannot be done once OpenGL is activated. 
    r->allocatedN *=8;
    r->particles = realloc(r->particles,sizeof(struct reb_particle)*r->allocatedN);
#endif // OPENGL

    // Start the integration
    reb_integrate(r, TMAX);

    // Cleanup
    MPI_Barrier(MPI_COMM_WORLD);
    reb_mpi_finalize(r);
    reb_free_simulation(r); 
}

char out_tree[16];
void heartbeat(struct reb_simulation* const r){
	move_to_com(r);

	if (reb_output_check(r,10.0*r->dt))
		reb_output_timing(r,0);

	if (reb_output_check(r, OUT_INTERVAL))
		write_sim(r, "particles.txt", 0);

	if (r->t > 0 && reb_output_check(r, (double) TMAX / OUT_TREE_N)) {
		sprintf(out_tree, "tree%d", (int) (r->t / ((double) TMAX / OUT_TREE_N))-1);
		write_tree(r, out_tree);
	}

	if (r->t >= TMAX) {
		reb_output_timing(r,0);
		printf("\n");
	}
	
}
