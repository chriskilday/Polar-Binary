/* output.h*/

void clean_files(char*, char*);
void write_header(char*, int, double, double, double);
void write_particles(struct reb_simulation* const, char* out_pt, double mr);
void write_particle(struct reb_particle, struct reb_particle, FILE*, struct reb_simulation*);
void write_tree(struct reb_simulation*, char*);
