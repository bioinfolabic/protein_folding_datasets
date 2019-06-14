#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <sys/time.h>

#ifndef TRUE
#define TRUE 1
#define FALSE 0
#endif

#define FAIL 0
#define SUCCESS 1
#define SEQ_MAX_LENGTH 1000
#define N_DIM 3
#define RAND_MAX_MT 4294967295 // Mersenne Twister maximum 2^32-1
#define MAX_MOL 1026

#define sqr(x) ((x) * (x))
#define min(x, y) (((x) < (y))? (x) : (y))

//typedef double double; // Problem with atomic operations

/***   Structures   ***/
#if N_DIM == 2
	typedef struct {
		double x;
		double y;
	} VectorR;
#else
	typedef struct {
		double x;
		double y;
		double z;
	} VectorR;
#endif

/***#if N_DIM == 2
	typedef struct {
		double x;
		double y;
	} VectorR;
	typedef struct {
		float x;
		float y;
	} VectorF;
#else
	typedef struct {
		double x;
		double y;
		double z;
	} VectorR;
	typedef struct {
		float x;
		float y;
		float z;
	} VectorF;
#endif***/

typedef struct {
	double u1;
	double u2;
	double u3;
	double u4;
} Quaternion;


typedef struct {
	VectorR v_r;
	VectorR v_v;
	VectorR v_a;
} Particle;


/*
typedef struct {
	VectorR v_r;
	VectorR v_v;
	VectorR v_a;
	VectorR v_a1;
	VectorR v_a2;
	VectorR torque;
	Quaternion q_r;
	Quaternion q_v;
	Quaternion q_a;
	Quaternion q_a1;
	Quaternion q_a2;
} Particle;
*/

typedef struct {
	VectorR v;
	int ik;
	int jk;
} Constraint;

Particle *particles;
Particle *best_structure;
Constraint * constraint;
int **mMat;
double *lMat;

/***   From file   ***/
char sequence[SEQ_MAX_LENGTH];
double mass;
int n_mol;
int bond_len;
int prot_len;
int n_c;

double LV;
double dt;
double r_cut;
double c_T;
double temperature;
double temperature_steps;

int display_interval;
int step_limit;
int step2resc_vels;
char report_file;
int step2report;
char print_summary;
int print_summary_interval;
char print_summary2file;
char save_pathways;
int pathways_step;

double shake_cons_prec;
int shake_max_cycle;
int shake_step2shake;

/***   Simulation control   ***/
time_t t_i, delta_t;
struct timeval t_now, last_t, result_t;
int i_step;
double vel_mag;
int en_update;

/***   Energies   ***/
double uSum;
double uBond;
double uTorsion;
double uLJ;
double *uLJ1;
double kinetic_energy;
double total_energy;

/***   SHAKE   ***/
int nCycleR;
int nCycleV;

/***   Other   ***/
double current_temperature;
double density;
double bond_avg;
VectorR center_mass;
double rGH, rGP, rG;

/***   Best   ***/
double best_potencial_energy;
int best_step;

/***   Device variables   ***/
unsigned int blockSize;
Particle *d_particles;
double *d_LV;
double *d_uB;
double *d_uT;
double *d_uLJ;
double *d_uLJVector;
char *d_sequence;

Constraint * d_constraint;

/*
double *d_particles_r;
double *d_particles_v;
double *d_particles_a;
*/
