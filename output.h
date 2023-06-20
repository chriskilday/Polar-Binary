/* output.h*/
#include "rebound.h"

void clean_files(char*);
void write_header(char*, int, double, double, double);
void write_sim(struct reb_simulation* const, char*);
struct reb_particle get_com(struct reb_simulation* const);
void move_to_com(struct reb_simulation* const);
