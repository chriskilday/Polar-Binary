/* output.c */

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include "rebound.h"
#include "tree.h"
#include "communication_mpi.h"

#define ASCIIFILE "ascii"
#define ORBITFILE "orbit"

int N_PARTICLES_INITIAL = 10;
int N_PARTICLES;

char* hash_index(uint32_t, char*);

void clean_files(char* ptfile, char* treefile) {
	char call[512];
	sprintf(call, "rm -f %s %s*.txt", ptfile, treefile);
	system(call);
}

void write_header(char* particle_file, int n, double mr, double a, double e) {
	FILE* pt_file = fopen(particle_file, "a");
	fprintf(pt_file, "*** Simulation\n*** N = %d\n*** MR = %f\tA = %f\tE = %f\n", n,mr,a,e);
	fclose(pt_file);
}

void particle_str(struct reb_particle pt, struct reb_particle com, char* pt_str, struct reb_simulation* s) {
	char hash_str[16];
	//struct reb_orbit o = reb_tools_particle_to_orbit(s->G, pt, com);
	sprintf(pt_str, "%s %f %f %f %f %f %f %f\n", hash_index(pt.hash, hash_str), pt.x, pt.y, pt.z,0.,0.,0.,0.);
									      //o.a, o.e, o.inc, o.Omega);
}

void tree_str(struct reb_treecell* node, char* node_str) {
	sprintf(node_str, "x: %f y: %f z: %f w: %f m: %f\n",
			node->x, node->y, node->z, node->w, node->m);
}

void write_particle(struct reb_particle pt, struct reb_particle com, FILE* of, struct reb_simulation* s) {
	char pt_str[512]; particle_str(pt, com, pt_str, s);
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

int nout = 0;
void write_tree_helper(struct reb_simulation* s, char* treefile, int index) {
	char treefile_mpi[512];
	sprintf(treefile_mpi, "%s%d_%d.txt", treefile, nout++, index);

	FILE* out_tree = fopen(treefile_mpi, "a");
	write_node(s->tree_root[index], out_tree);
	fclose(out_tree);
}

void write_tree(struct reb_simulation* s, char* treefile, int boxsize) {
	if (s->mpi_num == 1)
		for (int i = 0; i < boxsize*boxsize; i++)
			write_tree_helper(s, treefile, i);
	else
		write_tree_helper(s, treefile, s->mpi_id);
}

char* itostr(int i, char* out) {
	sprintf(out, "%d", i);
	return out;
}

/*
void set_hash(struct reb_simulation* s) {
	for (int i = 0; i < N_PARTICLES; i++) {
		char istr[11];
		s->particles[i+2].hash = reb_hash(itostr(i, istr));
	}
}
*/

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

/*
void order_str(struct reb_simulation* s, char* out) {
	out[0] = '\0';
	char istr[16];
	for (int i = 0; i < N_PARTICLES; i++)
		sprintf(out+strlen(out), "%s ", hash_index(s->particles[i+2].hash, istr));
}

void write_hash_order(struct reb_simulation* s) {
	char asciifile[512], orbitfile[512];
	sprintf(asciifile, "%s_%d", ASCIIFILE, s->mpi_id);
	sprintf(orbitfile, "%s_%d", ORBITFILE, s->mpi_id);
	FILE *af = fopen(asciifile, "a"), *of = fopen(orbitfile, "a");

	char order[N_PARTICLES*12];
	order_str(s, order);
	fprintf(af, "^^^ %s\n", order); fprintf(of, "^^^ %s\n", order);

	fclose(af); fclose(of);
}

void write_periodic(struct reb_simulation* s) {
	write_hash_order(s);

	reb_output_ascii(s, ASCIIFILE);
	reb_output_orbits(s, ORBITFILE);
}
*/

void write_particles(struct reb_simulation* const s, char* out_pt, double mr) {
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
		struct reb_particle com = {0}; com.m = (1.0 + mr);
		FILE* of = fopen(out_pt, "a");
		fprintf(of, "@@@ %f\n", s->t);

		for (int i = 0; i < s->N; i++)
			write_particle(s->particles[i], com, of, s);

		int n_cur = s->N;
		for (int i = 0; i < s->mpi_num; i++) {
			n_cur += s->particles_recv_N[i];
			for (int j = 0; j < s->particles_recv_N[i]; j++) {
				write_particle(s->particles_recv[i][j], com, of, s);
			}
		}
		fprintf(of, "### %d\n", n_cur);
		fclose(of);
	}

	MPI_Barrier(MPI_COMM_WORLD);
	for (int i = 0; i < s->mpi_num; i++) {
		s->particles_send_N[i] = 0;
		s->particles_recv_N[i] = 0;
	}
}

void write_sim(struct reb_simulation* const s, char* out_pt, double mr) {
	if (s->mpi_id != 0)
		for (int i = 0; i < s->N; i++)
			reb_communication_mpi_add_particle_to_send_queue(s, s->particles[i], 0);

	write_particles(s, out_pt, mr);
}

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
