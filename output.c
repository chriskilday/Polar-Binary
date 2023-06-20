/* output.c */

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include "output.h"
#include "tree.h"
#include "communication_mpi.h"

extern double A, E;
extern double MASS_RATIO;
extern double DISC_MASS;
extern double DISC_INNER, DISC_OUTER;
extern double MANUAL_INNER;
extern int N_INITIAL;

char* itostr(int i, char* out) {
	sprintf(out, "%d", i);
	return out;
}

// Reverse hash
char* hash_index(uint32_t hash, char* str) {
	if (hash == reb_hash("Star1")) return "Star1";
	if (hash == reb_hash("Star2")) return "Star2";

	for (int i = 0; i >= 0; i++) {
		char istr[16];
		if (hash == reb_hash(itostr(i, istr))) {
			sprintf(str, "%d", i);
			return str;
		}
	}

	return "-1";
}

// Output String for particle
void particle_str(struct reb_particle pt, char* pt_str, struct reb_simulation* const s, struct reb_particle com) {
	double mp = DISC_MASS / N_INITIAL;
	double inner_a = DISC_INNER * (1+E)*A, outer_a = DISC_OUTER * (1+E)*A;

	struct reb_particle primary = {0};
	primary.x = com.x; primary.y = com.y; primary.z = com.z;
	primary.vx = com.vx; primary.vy = com.vy; primary.vz = com.vz;

	double r = sqrt(pt.x*pt.x + pt.y*pt.y + pt.z*pt.z);
	primary.m = 1.0 + (DISC_MASS - 2*mp*(1/sqrt(inner_a) - 1/sqrt(MANUAL_INNER))) *
		(pow(r, -1.5) - pow(MANUAL_INNER, -1.5)) / (pow(outer_a, -1.5) - pow(MANUAL_INNER, -1.5));

	pt.sim = NULL;
	char hash_str[16];
	struct reb_orbit o = reb_tools_particle_to_orbit(s->G, pt, primary);
	sprintf(pt_str, "%s %f %f %f %f %f %f %f %f %f\n",
			hash_index(pt.hash, hash_str), pt.x, pt.y, pt.z, pt.vx, pt.vy, pt.vz, o.v, o.a, o.e);
}

// Output String for tree node
void tree_str(struct reb_treecell* node, char* node_str) {
	sprintf(node_str, "x: %f y: %f z: %f w: %f m: %f\n",
			node->x, node->y, node->z, node->w, node->m);
}

void write_particle(struct reb_particle pt, FILE* of, struct reb_simulation* const s, struct reb_particle com) {
	char pt_str[512]; particle_str(pt, pt_str, s, com);
	fprintf(of, pt_str);
}

void write_node(struct reb_treecell* node, FILE* out_tree) {
	if (node->pt >= 0) return;

	char node_str[512]; tree_str(node, node_str);
	fprintf(out_tree, node_str);

	for (int i = 0; i < 8; i++)
		if (node->oct[i] != NULL)
			write_node(node->oct[i], out_tree);
}

// Gather all particles on one processor and output
void write_particles(struct reb_simulation* const s, char* out_pt) {
	MPI_Barrier(MPI_COMM_WORLD);
	struct reb_particle com = get_com(s);

	for (int i = 0; i < s->mpi_num; i++)
		MPI_Scatter(s->particles_send_N, 1, MPI_INT, &(s->particles_recv_N[i]), 1, MPI_INT, i, MPI_COMM_WORLD);

	for (int i = 0; i < s->mpi_num; i++) {
		if (i == s->mpi_id) continue;

		while (s->particles_recv_Nmax[i] < s->particles_recv_N[i]) {
			s->particles_recv_Nmax[i] += 32;
			s->particles_recv[i] = realloc(s->particles_recv[i], sizeof(struct reb_particle)*s->particles_recv_Nmax[i]);
		}
	}

	MPI_Request request[s->mpi_num];
	for (int i = 0; i < s->mpi_num; i++) {
		if (i == s->mpi_id) continue;
		if (s->particles_recv_N[i] == 0) continue;

		MPI_Irecv(s->particles_recv[i], s->particles_recv_N[i], s->mpi_particle, i, i*s->mpi_num + s->mpi_id, MPI_COMM_WORLD, &(request[i]));
	}

	for (int i = 0; i < s->mpi_num; i++) {
		if (i == s->mpi_id) continue;
		if (s->particles_send_N[i] == 0) continue;

		MPI_Send(s->particles_send[i], s->particles_send_N[i], s->mpi_particle, i, s->mpi_id*s->mpi_num + i, MPI_COMM_WORLD);
	}

	for (int i = 0; i < s->mpi_num; i++) {
		if (i == s->mpi_id) continue;
		if (s->particles_recv_N[i] == 0) continue;

		MPI_Status status;
		MPI_Wait(&(request[i]), &status);
	}

	if (s->mpi_id == 0) {
		FILE* of = fopen(out_pt, "a");
		fprintf(of, "@@@ %f\n", s->t);

		for (int i = 0; i < s->N; i++)
			write_particle(s->particles[i], of, s, com);

		int n_cur = s->N;
		for (int i = 0; i < s->mpi_num; i++) {
			n_cur += s->particles_recv_N[i];
			for (int j = 0; j < s->particles_recv_N[i]; j++) {
				write_particle(s->particles_recv[i][j], of, s, com);
			}
		}
		fprintf(of, "### %d\n", n_cur);
		fclose(of);
	}

	MPI_Barrier(MPI_COMM_WORLD);
	for (int i = 0; i < s->mpi_num; i++) {
		s->particles_recv_N[i] = 0;
	}
}

// Write all particles in simulation
void write_sim(struct reb_simulation* const s, char* out_pt) {
	if (s->mpi_id != 0)
		for (int i = 0; i < s->N; i++)
			reb_communication_mpi_add_particle_to_send_queue(s, s->particles[i], 0);

	char out_pt_mpi[128]; sprintf(out_pt_mpi, "./%s/particles.txt", out_pt);
	write_particles(s, out_pt_mpi);
}

// Prepare files for output
void clean_files(char* ptfile) {
	char call[512];
	sprintf(call, "rm -rf %s pb.zip && mkdir %s", ptfile, ptfile);
	system(call);
}

// Write simulation parameters
void write_header(char* particle_file, int n, double mr, double a, double e) {
	char filename[128]; sprintf(filename, "%s/particles.txt", particle_file);
	FILE* pt_file = fopen(filename, "a");
	fprintf(pt_file, "*** Simulation\n*** N = %d\n*** MR = %f\tA = %f\tE = %f\n", n,mr,a,e);
	fclose(pt_file);
}

// get_com and move_to_com adapted for multiple processors
struct reb_particle get_com(struct reb_simulation* const s) {
	struct reb_particle com;
	if (s->mpi_num == 1)
		com = reb_get_com(s);
	else {
		if (s->mpi_id == 0) {
			struct reb_particle coms[s->mpi_num];
			coms[s->mpi_num-1] = reb_get_com(s);

			MPI_Request request[s->mpi_num-1];
			for (int i = 0; i < s->mpi_num-1; i++) {
				MPI_Irecv(&coms[i], 1, s->mpi_particle, i+1, (i+1)*s->mpi_num, MPI_COMM_WORLD, &(request[i]));
			}

			for (int i = 0; i < s->mpi_num-1; i++) {
				MPI_Status status;
				MPI_Wait(&(request[i]), &status);
			}

			double m_tot = 0;
			for (int i = 0; i < s->mpi_num; i++)
				m_tot += coms[i].m;

			struct reb_particle local_com = {0}; local_com.m = m_tot;
			for (int i = 0; i < s->mpi_num; i++) {
				double mr = coms[i].m/m_tot;
				local_com.x += coms[i].x * mr;
				local_com.y += coms[i].y * mr;
				local_com.z += coms[i].z * mr;
			}

			for (int i = 0; i < s->mpi_num-1; i++)
				MPI_Send(&local_com, 1, s->mpi_particle, i+1, s->mpi_num+(i+1), MPI_COMM_WORLD);

			com = local_com;
		} else {
			struct reb_particle local_com = reb_get_com(s);
			MPI_Send(&local_com, 1, s->mpi_particle, 0, s->mpi_id*s->mpi_num, MPI_COMM_WORLD);

			MPI_Request request;
			MPI_Status status;
			MPI_Irecv(&local_com, 1, s->mpi_particle, 0, s->mpi_num+s->mpi_id, MPI_COMM_WORLD, &request);
			MPI_Wait(&request, &status);

			com = local_com;		
		}
	}

	return com;
}

void move_to_com(struct reb_simulation* const s) {
	if (s->mpi_num == 1)
		reb_move_to_com(s);
	else {
		struct reb_particle com = get_com(s);
		for (int i = 0; i < s->N; i++) {
			s->particles[i].x -= com.x;
			s->particles[i].y -= com.y;
			s->particles[i].z -= com.z;
			s->particles[i].vx -= com.vx;
			s->particles[i].vy -= com.vy;
			s->particles[i].vz -= com.vz;
		}
	}
}
