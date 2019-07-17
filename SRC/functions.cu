
extern "C" {
    #include "functions.h"
    #include "mt.h"
}

__constant__ int d_nmol;

/***   Timer   ***/
void initTimer() {
	t_i = time(NULL);
	gettimeofday(&last_t, NULL);
}

void finTimer() {
	delta_t = time(NULL) - t_i;
    //printf("%lf", (double)delta_t);
}

void timeval_subtract (struct timeval *result, struct timeval *x, struct timeval *y) {
	int nsec;

	/* Perform the carry for the later subtraction by updating y. */
	if (x->tv_usec < y->tv_usec) {
		nsec = (y->tv_usec - x->tv_usec) / 1000000 + 1;
		y->tv_usec -= 1000000 * nsec;
		y->tv_sec += nsec;
	}
	if (x->tv_usec - y->tv_usec > 1000000) {
		nsec = (x->tv_usec - y->tv_usec) / 1000000;
		y->tv_usec += 1000000 * nsec;
		y->tv_sec -= nsec;
	}

	/* Compute the time remaining to wait. tv_usec is certainly positive. */
	result->tv_sec = x->tv_sec - y->tv_sec;
	result->tv_usec = x->tv_usec - y->tv_usec;
}

/***   Files   ***/
int getParameter(const char *field_name, const char *token, const char *format, void *variable, FILE *fi) {
	char line[SEQ_MAX_LENGTH+12];
	char *p_token, *p_aux;
	int flag_done = 0;
	int len, i;

	do {
		if(fgets(line, SEQ_MAX_LENGTH+12, fi) == 0)
			return FAIL;

		//line[SEQ_MAX_LENGTH+10] = '\0';
		len = strlen(line);
		if(line[len-1] == '\n')
			line[len-1] = '\0';

		switch(flag_done) {
			case 0:
				if(strstr(line, "/*") != 0) {
					flag_done = 1;
					if(strstr(line, "*/") != 0)
						flag_done = 2;
				}
				break;
			case 1:
				if(strstr(line, "*/") != 0)
					flag_done = 2;
				break;
			case 2:
				if(strstr(line, "/*") != 0) {
					flag_done = 1;
					if(strstr(line, "*/") != 0)
						flag_done = 2;
				} else {
					flag_done = 0;
				}
				break;
		}
	} while(flag_done != 0);

	if(strstr(line, field_name) != 0) {
		p_token = strtok(line, token);

		while(p_token != NULL) {
			p_aux = p_token;
			p_token = strtok(NULL, token);
		}

		if(p_aux[0] == ' ') {

			for(i = 0; i < len - 1; i++) {
				p_aux[i] = p_aux[i + 1];
			}
		}
		//puts(p_aux);
		sscanf(p_aux, format, variable);
		return SUCCESS;

	} else {
		return FAIL;
	}
}

int loadFile(char **argv) {
	FILE *fi;
    	fi = fopen(argv[1], "r");

	if(fi != 0) {
		getParameter("sequence", "=", "%s", sequence, fi);
		getParameter("ProtLen", "=", "%d", &prot_len, fi);		
		getParameter("LV", "=", "%lf", &LV, fi);
		getParameter("stepLimit", "=", "%d", &step_limit, fi);
		getParameter("temperature", "=", "%lf", &temperature, fi);
		getParameter("savepathways", "=", "%c", &save_pathways, fi);
		getParameter("pathwaysstep", "=", "%d", &pathways_step, fi);

		fclose(fi);

		mass = 1.0; 
		n_mol = prot_len;
		bond_len = 1.0;
		n_c = (prot_len-1);
		dt = 0.0001;
		r_cut = 2.38;
		c_T = 0.01;
		display_interval = 100;
		step2resc_vels = 1;
		temperature_steps = 0.10;
		report_file = 'y';
		step2report = 'n';
		print_summary = 160;
		print_summary_interval = 'y';
		print_summary2file = 'y';
		shake_cons_prec = 1.0e-06;
		shake_max_cycle = 1; 
		shake_step2shake = 10;

				

        if(n_mol < 1025)
		{
			if(n_mol > 512) {
                blockSize = 1024;
            } else {
                if(n_mol > 256) {
                    blockSize = 512;
                } else {
                    if(n_mol > 128) {
                        blockSize = 256;
                    } else {
                        if(n_mol > 64) {
                            blockSize = 128;
                        } else {
                            if(n_mol > 32) {
                                blockSize = 64;
                            } else {
                                blockSize = 32;
                            }
                        }
                    }
                }
            }
        } 
		else {
            return FAIL;
        }

        if(print_summary != 'y' && print_summary != 'Y') print_summary_interval = step_limit - 1;
        if(report_file != 'y' && report_file != 'Y') step2report = step_limit - 1;
        if(save_pathways != 'y' && save_pathways != 'Y') pathways_step = step_limit - 1;
        en_update = min(print_summary_interval, step2report);
        en_update = min(en_update, pathways_step);

		return SUCCESS;
	} else {
		return FAIL;
	}
}

void putParameters() {
	printf("Sequence >>                         %s\n", sequence);
	printf("Mass >>                             %lf\n", mass);
	printf("Number of Particles >>              %d\n", n_mol);
	printf("Bond Length >>                      %d\n", bond_len);
	printf("Protein Length >>                   %d\n", prot_len);
	printf("Number of Constraints >>            %d\n", n_c);
	printf("Dimension of the Box >>             %lf\n", LV);
	printf("Time Step >>                        %lf\n", dt);
	printf("Cutoff Radius >>                    %lf\n", r_cut);
	printf("Skin >>                             %lf\n", skin);
	printf("Steps to Update Neighbour List >>   %d\n", step2neighbour);
	printf("cT >>                               %lf\n", c_T);
	printf("Display Interval >>                 %d\n", display_interval);
	printf("Step Limit >>                       %d\n", step_limit);
	printf("Steps to Rescale Velocities >>      %d\n", step2resc_vels);
	printf("Temperature >>                      %lf\n", temperature);
	printf("Temperature Steps >>                %lf\n", temperature_steps);
	printf("Report File >>                      %c\n", report_file);
	printf("Steps to Report >>                  %d\n", step2report);
	printf("Print Summary >>                    %c\n", print_summary);
	printf("Print Summary Interval >>           %d\n", print_summary_interval);
	printf("Print Summary to File >>            %c\n", print_summary2file);
	printf("Shake-consPrec >>                   %lf\n", shake_cons_prec);
	printf("Shake-maxCycle >>                   %d\n", shake_max_cycle);
	printf("Steps to Shake >>                   %d\n", shake_step2shake);
	printf("Save Pathways >>                    %c\n", save_pathways);
	printf("Pathways Step >>                    %d\n", pathways_step);
}

/***   Utilities   ***/
double randdouble(double max) {
	double ret;


	ret = fabs(((randomMT() + (RAND_MAX_MT/2)) / ((double)(RAND_MAX_MT)) * max));

	return ret;
}

int isUnique(Particle *p, int last) {
	int i;
	for(i = 0; i < last; i++) {
		if((p[i].v_r.x == p[last].v_r.x) && (p[i].v_r.y == p[last].v_r.y) && (p[i].v_r.z == p[last].v_r.z))
			return 0;
	}
	return 1;
}

void verifyBoundary1(VectorR *v) {
	if(v->x < 0) {
		v->x += LV;
	} else {
		if(v->x >= LV)
			v->x -= LV;
	}

	if(v->y < 0) {
		v->y += LV;
	} else {
		if(v->y >= LV)
			v->y -= LV;
	}

	if(v->z < 0) {
		v->z += LV;
	} else {
		if(v->z >= LV)
			v->z -= LV;
	}
}

void verifyBoundary2(VectorR *v) {
	if(v->x >= 0.5 * LV) {
		v->x -= LV;
	} else {
		if(v->x < -0.5 * LV)
			v->x += LV;
	}

	if(v->y >= 0.5 * LV) {
		v->y -= LV;
	} else {
		if(v->y < -0.5 * LV)
			v->y += LV;
	}

	if(v->z >= 0.5 * LV) {
		v->z -= LV;
	} else {
		if(v->z < -0.5 * LV)
			v->z += LV;
	}
}

__device__ void cudaVerifyBoundary1(VectorR *v, double LV) {
	if(v->x < 0) {
		v->x += LV;
	} else {
		if(v->x >= LV)
			v->x -= LV;
	}

	if(v->y < 0) {
		v->y += LV;
	} else {
		if(v->y >= LV)
			v->y -= LV;
	}

	if(v->z < 0) {
		v->z += LV;
	} else {
		if(v->z >= LV)
			v->z -= LV;
	}
}

__device__ void cudaVerifyBoundary2(VectorR *v, double LV) {
    if(v->x >= 0.5 * LV) {
		v->x -= LV;
	} else {
		if(v->x < -0.5 * LV)
			v->x += LV;
	}

	if(v->y >= 0.5 * LV) {
		v->y -= LV;
	} else {
		if(v->y < -0.5 * LV)
			v->y += LV;
	}

	if(v->z >= 0.5 * LV) {
		v->z -= LV;
	} else {
		if(v->z < -0.5 * LV)
			v->z += LV;
	}
}

/***   Initialize   ***/
void allocArrays() {
	int i;

	particles = (Particle *) malloc(n_mol * sizeof(Particle));

	best_structure = (Particle *) malloc(n_mol * sizeof(Particle));

	constraint = (Constraint *) malloc(n_c * sizeof(Constraint));

	mMat = (int **) malloc(n_mol * sizeof(int *));
	mMat[0] = (int *) malloc(n_mol * n_c * sizeof(int));
	for(i = 1; i < n_mol; i++)
		mMat[i] = mMat[0] + i * n_c;

	
}

void allocDevice() {
	cudaMemcpyToSymbol(d_nmol, &n_mol, sizeof(int));
    cudaMalloc((void **)&d_particles, blockSize * sizeof(Particle));
    cudaMalloc((void **)&d_uB, sizeof(double));
    cudaMalloc((void **)&d_uT, sizeof(double));
    cudaMalloc((void **)&d_uLJ, sizeof(double));
    cudaMalloc((void **)&d_uLJVector, blockSize * sizeof(double));
    cudaMalloc((void **)&d_sequence, n_mol * sizeof(char));
    cudaMalloc((void **)&d_neighbour_list, n_mol * n_mol * sizeof(int));
    cudaMalloc((void**)&d_constraint, n_c * sizeof(Constraint));
}

void setVelMag() {
	vel_mag = sqrt(N_DIM * (1. - 1./prot_len) * temperature);
}

void initCoords() {
	int i, j, conf_OK;
	double theta;
	double phi;
    double dist;

	particles[0].v_r.x = LV/2;
	particles[0].v_r.y = LV/2;
	particles[0].v_r.z = LV/2;
	//printf("(%lf, %lf, %lf)\n", particles[0].v_r.x, particles[0].v_r.y, particles[0].v_r.z);
	for(i = 1;i < n_mol; i++) {
		do {
			theta = randdouble(180.0) * M_PI / 180.0;
			phi = randdouble(360.0) * M_PI / 180.0;


			particles[i].v_r.x = particles[i-1].v_r.x + sin(theta) * cos(phi);
			particles[i].v_r.y = particles[i-1].v_r.y + sin(theta) * sin(phi);
			particles[i].v_r.z = particles[i-1].v_r.z + cos(theta);
            //printf("(%lf, %lf, %lf) theta %lf\n", particles[i].v_r.x, particles[i].v_r.y, particles[i].v_r.z, theta);

			verifyBoundary1(&particles[i].v_r);

            conf_OK = 1; //OK
			for (j = 0; j < i; j++)
			{
				dist = sqrt(sqr(particles[i].v_r.x - particles[j].v_r.x) + sqr(particles[i].v_r.y - particles[j].v_r.y) + sqr(particles[i].v_r.z - particles[j].v_r.z));
				if (dist < 1)
					conf_OK = 0; 
			}
            //printf("\t(%lf, %lf, %lf)\n", particles[i].v_r.x, particles[i].v_r.y, particles[i].v_r.z);
		} while((isUnique(particles, i) != 1) || (conf_OK != 1));
	}
}

void initVels() {
	int i;
	double theta;
	double phi;
	VectorR sum_v;

	sum_v.x = 0;
	sum_v.y = 0;
	sum_v.z = 0;

	/***   Generating Vectors   ***/
	setVelMag();
	for(i = 0; i < n_mol; i++) {
		theta = randdouble(180.0) * M_PI / 180.0;
		phi = randdouble(360.0) * M_PI / 360.0;

		particles[i].v_v.x = sin(theta) * cos(phi) * vel_mag;
		particles[i].v_v.y = sin(theta) * sin(phi) * vel_mag;
		particles[i].v_v.z = cos(theta) * vel_mag;
        //printf("(%lf, %lf, %lf)\n", particles[i].v_v.x/vel_mag, particles[i].v_v.y/vel_mag, particles[i].v_v.z/vel_mag);

		sum_v.x += particles[i].v_v.x;
		sum_v.y += particles[i].v_v.y;
		sum_v.z += particles[i].v_v.z;

	}

	/***   Center of Mass at Rest   ***/
	for(i = 0; i < n_mol; i++) {
		particles[i].v_v.x = particles[i].v_v.x - (sum_v.x / prot_len);
		particles[i].v_v.y = particles[i].v_v.y - (sum_v.y / prot_len);
		particles[i].v_v.z = particles[i].v_v.z - (sum_v.z / prot_len);

	}
}

void initAccs() {
	int i;

	for(i = 0; i < n_mol; i++) {
		particles[i].v_a.x = 0.;
		particles[i].v_a.y = 0.;
		particles[i].v_a.z = 0.;
    }
}

void buildConstMatrix() {
	int i, j;

	for(i = 0; i < n_mol; i++)
		for(j = 0; j < n_c; j++)
			mMat[i][j] = 0;

	for(i = 0; i < n_mol; i++) {
		if(i - 1 > 0)
			mMat[i - 1][i] = 2;
		if(i < n_c)
			mMat[i][i] = -2;
	}

	for(i = 0; i < n_c; i++) {
		/***   Distance sqr taken as 1   ***/
		constraint[i].ik = i;
		constraint[i].jk = i + 1;

        //printf("%d %d\n", constraint[i].ik, constraint[i].jk);
	}
}

void initMD() {
	allocArrays();
    allocDevice();
	initCoords();
	initVels();
	initAccs();
	buildConstMatrix();

    cudaMemcpy(d_sequence, sequence, n_mol * sizeof(char), cudaMemcpyHostToDevice);

}

/***   MD   ***/
__global__ void cudaUpdatePos(Particle *particles, double dt ,double LV) {
	//int threadId = threadIdx.x;

	for(int threadId = threadIdx.x; threadId < d_nmol; threadId += blockDim.x){
		particles[threadId].v_r.x += particles[threadId].v_v.x * dt + 0.5 * particles[threadId].v_a.x * sqr(dt);
		particles[threadId].v_r.y += particles[threadId].v_v.y * dt + 0.5 * particles[threadId].v_a.y * sqr(dt);
		particles[threadId].v_r.z += particles[threadId].v_v.z * dt + 0.5 * particles[threadId].v_a.z * sqr(dt);

        //cudaVerifyBoundary1(&particles[i].v_r);
        if(particles[threadId].v_r.x < 0) {
    		particles[threadId].v_r.x += LV;
    	} else {
    		if(particles[threadId].v_r.x >= LV)
    			particles[threadId].v_r.x -= LV;
    	}

    	if(particles[threadId].v_r.y < 0) {
    		particles[threadId].v_r.y += LV;
    	} else {
    		if(particles[threadId].v_r.y >= LV)
    			particles[threadId].v_r.y -= LV;
    	}

    	if(particles[threadId].v_r.z < 0) {
    		particles[threadId].v_r.z += LV;
    	} else {
    		if(particles[threadId].v_r.z >= LV)
    			particles[threadId].v_r.z -= LV;
    	}

        // First part of the velocity verlet
		particles[threadId].v_v.x += 0.5 * dt * particles[threadId].v_a.x;
		particles[threadId].v_v.y += 0.5 * dt * particles[threadId].v_a.y;
		particles[threadId].v_v.z += 0.5 * dt * particles[threadId].v_a.z;

        // Initialize the step acceleration
		particles[threadId].v_a.x = 0.;
		particles[threadId].v_a.y = 0.;
		particles[threadId].v_a.z = 0.;

	}
}


__global__ void cudaUpdateNeighbourList(Particle *particles, int *neighbour_list, double r_cut, double skin) {

	for(int i = blockIdx.x; i < d_nmol - 2; i += gridDim.x) {
		for(int j = threadIdx.x; j >= i + 2 && j < d_nmol; j += blockDim.x) {
			double r = sqrt(sqr(particles[j].v_r.x - particles[i].v_r.x) + sqr(particles[j].v_r.y - particles[i].v_r.y) + sqr(particles[j].v_r.z - particles[i].v_r.z));

			if(r <= (r_cut + skin))
				neighbour_list[i * d_nmol + j] = 1;
			else
            	neighbour_list[i * d_nmol + j] = 0;
		}
	}
}

__global__ void cudaBondEnergy(Particle *particles, double *d_uB, double LV) {
    // Restrict to one block. Limit: 1024 threads/amino acids (compute capability >= 2.x)
    extern __shared__ double d_uBond[];
    int blockSize = blockDim.x;
    int threadId = threadIdx.x;
    VectorR dr1, dr2, a1, a2;
    double c11, c12, c22, cd, f;
    //double c

    if(threadId < (d_nmol - 2)) {
        dr1.x = particles[threadId + 1].v_r.x - particles[threadId].v_r.x;
        dr1.y = particles[threadId + 1].v_r.y - particles[threadId].v_r.y;
        dr1.z = particles[threadId + 1].v_r.z - particles[threadId].v_r.z;
        dr2.x = particles[threadId + 2].v_r.x - particles[threadId + 1].v_r.x;
        dr2.y = particles[threadId + 2].v_r.y - particles[threadId + 1].v_r.y;
        dr2.z = particles[threadId + 2].v_r.z - particles[threadId + 1].v_r.z;

        cudaVerifyBoundary2(&dr1, LV);
        cudaVerifyBoundary2(&dr2, LV);

        c11 = dr1.x * dr1.x + dr1.y * dr1.y + dr1.z * dr1.z;
        c12 = dr1.x * dr2.x + dr1.y * dr2.y + dr1.z * dr2.z;
        c22 = dr2.x * dr2.x + dr2.y * dr2.y + dr2.z * dr2.z;

        cd = sqrt (c11 * c22);


        d_uBond[threadId] = c12;

        f = -1.0;
        a1.x = (f / cd) * ((c12 / c11) * dr1.x  - dr2.x);
        a1.y = (f / cd) * ((c12 / c11) * dr1.y  - dr2.y);
        a1.z = (f / cd) * ((c12 / c11) * dr1.z  - dr2.z);
        a2.x = (f / cd) * (dr1.x - (c12 / c22) * dr2.x);
        a2.y = (f / cd) * (dr1.y - (c12 / c22) * dr2.y);
        a2.z = (f / cd) * (dr1.z - (c12 / c22) * dr2.z);

        particles[threadId + 2].v_a.x = 0;
        particles[threadId + 2].v_a.y = 0;
        particles[threadId + 2].v_a.z = 0;
        __syncthreads();
        particles[threadId].v_a.x = a1.x;
        particles[threadId].v_a.y = a1.y;
        particles[threadId].v_a.z = a1.z;
        __syncthreads();
        particles[threadId + 1].v_a.x -= a1.x;
        particles[threadId + 1].v_a.y -= a1.y;
        particles[threadId + 1].v_a.z -= a1.z;
        particles[threadId + 1].v_a.x -= a2.x;
        particles[threadId + 1].v_a.y -= a2.y;
        particles[threadId + 1].v_a.z -= a2.z;
        __syncthreads();
        particles[threadId + 2].v_a.x += a2.x;
        particles[threadId + 2].v_a.y += a2.y;
        particles[threadId + 2].v_a.z += a2.z;

    } else {
        d_uBond[threadId] = 0;
    }

    if (blockSize >= 512) {
        if (threadId < 256)
            d_uBond[threadId] += d_uBond[threadId + 256];
        __syncthreads();
    }
    if (blockSize >= 256) {
        if (threadId < 128)
            d_uBond[threadId] += d_uBond[threadId + 128];
        __syncthreads();
    }
    if (blockSize >= 128) {
        if (threadId < 64)
            d_uBond[threadId] += d_uBond[threadId + 64];
        __syncthreads();
    }

    if (threadId < 32) {
        if (blockSize >= 64) d_uBond[threadId] += d_uBond[threadId + 32];
        if (blockSize >= 32) d_uBond[threadId] += d_uBond[threadId + 16];
        if (blockSize >= 16) d_uBond[threadId] += d_uBond[threadId + 8];
        if (blockSize >= 8) d_uBond[threadId] += d_uBond[threadId + 4];
        if (blockSize >= 4) d_uBond[threadId] += d_uBond[threadId + 2];
        if (blockSize >= 2) d_uBond[threadId] += d_uBond[threadId + 1];

        if(threadId == 0) {
            *d_uB = d_uBond[0];
        }
    }

}

__global__ void cudaTorsionEnergy(Particle *particles, double *d_uT, double LV) {
    extern __shared__ double d_uTorsion[];
    int blockSize = blockDim.x;
    int threadId = threadIdx.x;
    VectorR dr1, dr2, dr3, a1, a2;
    double c11, c12, c13, c22, c23, c33, pi, qia, qib, cr1, cr2, t1, t2, t3, t4, t5, t6, f;

    if(threadId < (d_nmol - 3)) {
        dr1.x = particles[threadId + 1].v_r.x - particles[threadId].v_r.x;
        dr1.y = particles[threadId + 1].v_r.y - particles[threadId].v_r.y;
        dr1.z = particles[threadId + 1].v_r.z - particles[threadId].v_r.z;
        dr2.x = particles[threadId + 2].v_r.x - particles[threadId + 1].v_r.x;
        dr2.y = particles[threadId + 2].v_r.y - particles[threadId + 1].v_r.y;
        dr2.z = particles[threadId + 2].v_r.z - particles[threadId + 1].v_r.z;
        dr3.x = particles[threadId + 3].v_r.x - particles[threadId + 2].v_r.x;
        dr3.y = particles[threadId + 3].v_r.y - particles[threadId + 2].v_r.y;
        dr3.z = particles[threadId + 3].v_r.z - particles[threadId + 2].v_r.z;

        cudaVerifyBoundary2(&dr1, LV);
        cudaVerifyBoundary2(&dr2, LV);
        cudaVerifyBoundary2(&dr3, LV);

        c11 = dr1.x * dr1.x + dr1.y * dr1.y + dr1.z * dr1.z;
        c12 = dr1.x * dr2.x + dr1.y * dr2.y + dr1.z * dr2.z;
        c13 = dr1.x * dr3.x + dr1.y * dr3.y + dr1.z * dr3.z;
        c22 = dr2.x * dr2.x + dr2.y * dr2.y + dr2.z * dr2.z;
        c23 = dr2.x * dr3.x + dr2.y * dr3.y + dr2.z * dr3.z;
        c33 = dr3.x * dr3.x + dr3.y * dr3.y + dr3.z * dr3.z;

        pi = c13 * c22 - c12 * c23;
        qia = c11 * c22 - c12 * c12;
        qib = c22 * c33 - c23 * c23;
        cr1 = c12 / c22;
        cr2 = c23 / c22;

        t1 = pi;
        t2 = c11 * c23 - c12 * c13;
        t3 = - qia;
        t4 = qib;
        t5 = c13 * c23 - c12 * c33;
        t6 = - t1;

        d_uTorsion[threadId] = (- 0.5) * c13;

        f = 0.5;

        a1.x = f * c22 * (t1 * dr1.x + t2 * dr2.x + t3 * dr3.x) / (sqrt(qia * qib) * qia);
        a1.y = f * c22 * (t1 * dr1.y + t2 * dr2.y + t3 * dr3.y) / (sqrt(qia * qib) * qia);
        a1.z = f * c22 * (t1 * dr1.z + t2 * dr2.z + t3 * dr3.z) / (sqrt(qia * qib) * qia);
        a2.x = f * c22 * (t4 * dr1.x + t5 * dr2.x + t6 * dr3.x) / (sqrt(qia * qib) * qib);
        a2.y = f * c22 * (t4 * dr1.y + t5 * dr2.y + t6 * dr3.y) / (sqrt(qia * qib) * qib);
        a2.z = f * c22 * (t4 * dr1.z + t5 * dr2.z + t6 * dr3.z) / (sqrt(qia * qib) * qib);

        particles[threadId].v_a.x += a1.x;
        particles[threadId].v_a.y += a1.y;
        particles[threadId].v_a.z += a1.z;
        __syncthreads();
        particles[threadId + 1].v_a.x += (- (1. + cr1) * a1.x) + cr2 * a2.x;
        particles[threadId + 1].v_a.y += (- (1. + cr1) * a1.y) + cr2 * a2.y;
        particles[threadId + 1].v_a.z += (- (1. + cr1) * a1.z) + cr2 * a2.z;
        __syncthreads();
        particles[threadId + 2].v_a.x += cr1 * a1.x + (- (1. + cr2)) * a2.x;
        particles[threadId + 2].v_a.y += cr1 * a1.y + (- (1. + cr2)) * a2.y;
        particles[threadId + 2].v_a.z += cr1 * a1.z + (- (1. + cr2)) * a2.z;
        __syncthreads();
        particles[threadId + 3].v_a.x += a2.x;
        particles[threadId + 3].v_a.y += a2.y;
        particles[threadId + 3].v_a.z += a2.z;

    } else {
        d_uTorsion[threadId] = 0.;
    }


    if (blockSize >= 1024) {
        if (threadId < 512)
            d_uTorsion[threadId] += d_uTorsion[threadId + 512];
        __syncthreads();
    }
    if (blockSize >= 512) {
        if (threadId < 256)
            d_uTorsion[threadId] += d_uTorsion[threadId + 256];
        __syncthreads();
    }
    if (blockSize >= 256) {
        if (threadId < 128)
            d_uTorsion[threadId] += d_uTorsion[threadId + 128];
        __syncthreads();
    }
    if (blockSize >= 128) {
        if (threadId < 64)
            d_uTorsion[threadId] += d_uTorsion[threadId + 64];
        __syncthreads();
    }

    if (threadId < 32) {
        if (blockSize >= 64) d_uTorsion[threadId] += d_uTorsion[threadId + 32];
        if (blockSize >= 32) d_uTorsion[threadId] += d_uTorsion[threadId + 16];
        if (blockSize >= 16) d_uTorsion[threadId] += d_uTorsion[threadId + 8];
        if (blockSize >= 8) d_uTorsion[threadId] += d_uTorsion[threadId + 4];
        if (blockSize >= 4) d_uTorsion[threadId] += d_uTorsion[threadId + 2];
        if (blockSize >= 2) d_uTorsion[threadId] += d_uTorsion[threadId + 1];

        if(threadId == 0) {
            *d_uT = d_uTorsion[0];

		}
    }

}

__global__ void cudaLJcomp(Particle *particles, int *d_neighbour_list, char *d_sequence, double LV, double *d_uLJVector, int i, int n_mol) {

	extern __shared__ double d_uLJComp[];

	int j = threadIdx.x;
	int blockSize = blockDim.x;
	double r2, u_LJ, f_LJ;
	VectorR dr1;

	if((j < n_mol) && (j > i + 1)) {

    		dr1.x = particles[i].v_r.x - particles[j].v_r.x;
    		dr1.y = particles[i].v_r.y - particles[j].v_r.y;
    		dr1.z = particles[i].v_r.z - particles[j].v_r.z;

    		cudaVerifyBoundary2(&dr1, LV);

    		r2 = (dr1.x * dr1.x) + (dr1.y * dr1.y) + (dr1.z * dr1.z);
    		u_LJ = 4. * (pow(r2, -6) - pow(r2, -3));
    		f_LJ = 24. * (2. * pow(r2, -7) - pow(r2, -4));


            
            if( (d_sequence[i] != 'A') || (d_sequence[j] != 'A') ) {
    			u_LJ = 0.5 * u_LJ;
    			f_LJ = 0.5 * f_LJ;
    		}

    		atomicAdd(&particles[i].v_a.x, (dr1.x * f_LJ));
    		atomicAdd(&particles[i].v_a.y, (dr1.y * f_LJ));
    		atomicAdd(&particles[i].v_a.z, (dr1.z * f_LJ));
    		atomicAdd(&particles[j].v_a.x, -(dr1.x * f_LJ));
    		atomicAdd(&particles[j].v_a.y, -(dr1.y * f_LJ));
    		atomicAdd(&particles[j].v_a.z, -(dr1.z * f_LJ));

    		d_uLJComp[j] = u_LJ;
			
		
    } else {

		d_uLJComp[j] = 0.;
	}
	
    __syncthreads(); 


    if (blockSize >= 1024) {
        if (j < 512){
            d_uLJComp[j] += d_uLJComp[j + 512];}
        __syncthreads();
    }
    if (blockSize >= 512) {
        if (j < 256){
            d_uLJComp[j] += d_uLJComp[j + 256];}
        __syncthreads();
    }
    if (blockSize >= 256) {
        if (j < 128){
            d_uLJComp[j] += d_uLJComp[j + 128];}
        __syncthreads();
    }
    if (blockSize >= 128) {
        if (j < 64){
            d_uLJComp[j] += d_uLJComp[j + 64];}
        __syncthreads();
    }
	
    
        
	if (blockSize >= 64){
		if(j < 32){
			d_uLJComp[j] += d_uLJComp[j + 32];}
	}
	__syncthreads();
	if (blockSize >= 32){
		if(j < 16){
			d_uLJComp[j] += d_uLJComp[j + 16];}
	}
	__syncthreads();
	if (blockSize >= 16){
		if(j < 8){
			d_uLJComp[j] += d_uLJComp[j + 8];}
	}
	__syncthreads();
	if (blockSize >= 8){
		if(j < 4){
			d_uLJComp[j] += d_uLJComp[j + 4];} 
	}
	__syncthreads();
	if (blockSize >= 4){
		if(j < 2){
			d_uLJComp[j] += d_uLJComp[j + 2];}
	}
	__syncthreads();
	if (blockSize >= 2){
		if(j < 1){
			d_uLJComp[j] += d_uLJComp[j + 1];}
	}
	__syncthreads();
	if (blockSize >= 1){
		d_uLJVector[i] = d_uLJComp[0];
	}

	

}

__global__ void cudaLJ(Particle *particles, int *d_neighbour_list, char *d_sequence, double LV, double *d_uT, double *d_uB, double *d_uLJ, double *d_uLJVector, int step, int n_mol){

	extern __shared__ double d_u[];

	// Other Variables
	unsigned int threadId = threadIdx.x;
	unsigned int blockId = blockIdx.x;; 
	unsigned int blockSize = blockDim.x;

	
	if(blockId == 2) {
		if(threadId < n_mol-2)
		{
			__syncthreads();
			if(threadId == 0){
					cudaLJcomp <<< 1, blockSize, (blockSize + 1) * sizeof(double) >>> (particles, d_neighbour_list, d_sequence, LV, d_uLJVector, threadId, n_mol);
			}
			__syncthreads();
			if(threadId < 20){
				if(threadId > 0){
					cudaLJcomp <<< 1, blockSize, (blockSize + 1) * sizeof(double) >>> (particles, d_neighbour_list, d_sequence, LV, d_uLJVector, threadId, n_mol);}
			}
			__syncthreads();
			if(threadId < 40){
				if(threadId > 19){
					cudaLJcomp <<< 1, blockSize, (blockSize + 1) * sizeof(double) >>> (particles, d_neighbour_list, d_sequence, LV, d_uLJVector, threadId, n_mol);}
			}
			__syncthreads();
			if(threadId < 60){
				if(threadId > 39){
					cudaLJcomp <<< 1, blockSize, (blockSize + 1) * sizeof(double) >>> (particles, d_neighbour_list, d_sequence, LV, d_uLJVector, threadId, n_mol);}
			}
			__syncthreads();
			if(threadId < 80){
				if(threadId > 59){
					cudaLJcomp <<< 1, blockSize, (blockSize + 1) * sizeof(double) >>> (particles, d_neighbour_list, d_sequence, LV, d_uLJVector, threadId, n_mol);}
			}
			__syncthreads();
			if(threadId < 100){
				if(threadId > 79){
					cudaLJcomp <<< 1, blockSize, (blockSize + 1) * sizeof(double) >>> (particles, d_neighbour_list, d_sequence, LV, d_uLJVector, threadId, n_mol);}
			}
			__syncthreads();
			if(threadId < 120){
				if(threadId > 99){
					cudaLJcomp <<< 1, blockSize, (blockSize + 1) * sizeof(double) >>> (particles, d_neighbour_list, d_sequence, LV, d_uLJVector, threadId, n_mol);}
			}
			__syncthreads();
			if(threadId < 140){
				if(threadId > 119){
					cudaLJcomp <<< 1, blockSize, (blockSize + 1) * sizeof(double) >>> (particles, d_neighbour_list, d_sequence, LV, d_uLJVector, threadId, n_mol);}
			}
			__syncthreads();
			if(threadId < 160){
				if(threadId > 139){
					cudaLJcomp <<< 1, blockSize, (blockSize + 1) * sizeof(double) >>> (particles, d_neighbour_list, d_sequence, LV, d_uLJVector, threadId, n_mol);}
			}
			__syncthreads();
			if(threadId < 180){
				if(threadId > 159){
					cudaLJcomp <<< 1, blockSize, (blockSize + 1) * sizeof(double) >>> (particles, d_neighbour_list, d_sequence, LV, d_uLJVector, threadId, n_mol);}
			}
			__syncthreads();
			if(threadId < 200){
				if(threadId > 179){
					cudaLJcomp <<< 1, blockSize, (blockSize + 1) * sizeof(double) >>> (particles, d_neighbour_list, d_sequence, LV, d_uLJVector, threadId, n_mol);}
			}
			__syncthreads();
			if(threadId < 220){
				if(threadId > 199){
					cudaLJcomp <<< 1, blockSize, (blockSize + 1) * sizeof(double) >>> (particles, d_neighbour_list, d_sequence, LV, d_uLJVector, threadId, n_mol);}
			}
			__syncthreads();
			if(threadId > 119){
					cudaLJcomp <<< 1, blockSize, (blockSize + 1) * sizeof(double) >>> (particles, d_neighbour_list, d_sequence, LV, d_uLJVector, threadId, n_mol);
			}
			__syncthreads();
			d_u[threadId] = d_uLJVector[threadId];
			//printf("d_u[%d] = %.8lf\n",threadId, d_u[threadId] );

		}
		else {
			//printf("else - threadId: %d\n", threadId);
			d_u[threadId] = 0.;
		}
	}

}

__global__ void cudaEnergy(Particle *particles, int *d_neighbour_list, char *d_sequence, double LV, double *d_uT, double *d_uB, double *d_uLJ, double *d_uLJVector, int step) {

	extern __shared__ double d_u[];

	// Other Variables
	unsigned int threadId = threadIdx.x;
	unsigned int blockId = blockIdx.x;
	unsigned int blockSize = blockDim.x;
	VectorR dr1, dr2, dr3, a1, a2;
	double cd, c11, c12, c13, c22, c23, c33, pi, qia, qib, cr1, cr2, t1, t2, t3, t4, t5, t6, f;
    //double c

	// ============================================================================================================================ //
	// Bond Energy
	if(blockId == 0) {
		if(threadId < (d_nmol - 2)) {
			dr1.x = particles[threadId + 1].v_r.x - particles[threadId].v_r.x;
			dr1.y = particles[threadId + 1].v_r.y - particles[threadId].v_r.y;
			dr1.z = particles[threadId + 1].v_r.z - particles[threadId].v_r.z;
			dr2.x = particles[threadId + 2].v_r.x - particles[threadId + 1].v_r.x;
			dr2.y = particles[threadId + 2].v_r.y - particles[threadId + 1].v_r.y;
			dr2.z = particles[threadId + 2].v_r.z - particles[threadId + 1].v_r.z;

			cudaVerifyBoundary2(&dr1, LV);  
			cudaVerifyBoundary2(&dr2, LV);

			c11 = dr1.x * dr1.x + dr1.y * dr1.y + dr1.z * dr1.z;
			c12 = dr1.x * dr2.x + dr1.y * dr2.y + dr1.z * dr2.z;
			c22 = dr2.x * dr2.x + dr2.y * dr2.y + dr2.z * dr2.z;

			cd = sqrt (c11 * c22);


			d_u[threadId] = c12;


			f = -1.0;
			a1.x = (f / cd) * ((c12 / c11) * dr1.x  - dr2.x);
			a1.y = (f / cd) * ((c12 / c11) * dr1.y  - dr2.y);
			a1.z = (f / cd) * ((c12 / c11) * dr1.z  - dr2.z);
			a2.x = (f / cd) * (dr1.x - (c12 / c22) * dr2.x);
			a2.y = (f / cd) * (dr1.y - (c12 / c22) * dr2.y);
			a2.z = (f / cd) * (dr1.z - (c12 / c22) * dr2.z);

			atomicAdd(&particles[threadId].v_a.x, a1.x);
			atomicAdd(&particles[threadId].v_a.y, a1.y);
			atomicAdd(&particles[threadId].v_a.z, a1.z);
			atomicAdd(&particles[threadId + 1].v_a.x, -(a1.x + a2.x));
			atomicAdd(&particles[threadId + 1].v_a.y, -(a1.y + a2.y));
			atomicAdd(&particles[threadId + 1].v_a.z, -(a1.z + a2.z));
			atomicAdd(&particles[threadId + 2].v_a.x, a2.x);
			atomicAdd(&particles[threadId + 2].v_a.y, a2.y);
			atomicAdd(&particles[threadId + 2].v_a.z, a2.z);

		} else {
			d_u[threadId] = 0.;
		}
	}

	// Torsion Energy
	if(blockId == 1) {
		if(threadId < (d_nmol - 3)) {
			dr1.x = particles[threadId + 1].v_r.x - particles[threadId].v_r.x;
			dr1.y = particles[threadId + 1].v_r.y - particles[threadId].v_r.y;
			dr1.z = particles[threadId + 1].v_r.z - particles[threadId].v_r.z;
			dr2.x = particles[threadId + 2].v_r.x - particles[threadId + 1].v_r.x;
			dr2.y = particles[threadId + 2].v_r.y - particles[threadId + 1].v_r.y;
			dr2.z = particles[threadId + 2].v_r.z - particles[threadId + 1].v_r.z;
			dr3.x = particles[threadId + 3].v_r.x - particles[threadId + 2].v_r.x;
			dr3.y = particles[threadId + 3].v_r.y - particles[threadId + 2].v_r.y;
			dr3.z = particles[threadId + 3].v_r.z - particles[threadId + 2].v_r.z;



			cudaVerifyBoundary2(&dr1, LV);
			cudaVerifyBoundary2(&dr2, LV);
			cudaVerifyBoundary2(&dr3, LV);

			c11 = dr1.x * dr1.x + dr1.y * dr1.y + dr1.z * dr1.z;
			c12 = dr1.x * dr2.x + dr1.y * dr2.y + dr1.z * dr2.z;
			c13 = dr1.x * dr3.x + dr1.y * dr3.y + dr1.z * dr3.z;
			c22 = dr2.x * dr2.x + dr2.y * dr2.y + dr2.z * dr2.z;
			c23 = dr2.x * dr3.x + dr2.y * dr3.y + dr2.z * dr3.z;
			c33 = dr3.x * dr3.x + dr3.y * dr3.y + dr3.z * dr3.z;

			pi = c13 * c22 - c12 * c23;
			qia = c11 * c22 - c12 * c12;
			qib = c22 * c33 - c23 * c23;
			cr1 = c12 / c22;
			cr2 = c23 / c22;

			t1 = pi;
			t2 = c11 * c23 - c12 * c13;
			t3 = - qia;
			t4 = qib;
			t5 = c13 * c23 - c12 * c33;
			t6 = - t1;

			d_u[threadId] = (- 0.5) * c13;


			f = 0.5;

			a1.x = f * c22 * (t1 * dr1.x + t2 * dr2.x + t3 * dr3.x) / (sqrt(qia * qib) * qia);
			a1.y = f * c22 * (t1 * dr1.y + t2 * dr2.y + t3 * dr3.y) / (sqrt(qia * qib) * qia);
			a1.z = f * c22 * (t1 * dr1.z + t2 * dr2.z + t3 * dr3.z) / (sqrt(qia * qib) * qia);
			a2.x = f * c22 * (t4 * dr1.x + t5 * dr2.x + t6 * dr3.x) / (sqrt(qia * qib) * qib);
			a2.y = f * c22 * (t4 * dr1.y + t5 * dr2.y + t6 * dr3.y) / (sqrt(qia * qib) * qib);
			a2.z = f * c22 * (t4 * dr1.z + t5 * dr2.z + t6 * dr3.z) / (sqrt(qia * qib) * qib);


			atomicAdd(&particles[threadId].v_a.x, a1.x);
			atomicAdd(&particles[threadId].v_a.y, a1.y);
			atomicAdd(&particles[threadId].v_a.z, a1.z);
			atomicAdd(&particles[threadId + 1].v_a.x, (- (1. + cr1) * a1.x) + cr2 * a2.x);
			atomicAdd(&particles[threadId + 1].v_a.y, (- (1. + cr1) * a1.y) + cr2 * a2.y);
			atomicAdd(&particles[threadId + 1].v_a.z, (- (1. + cr1) * a1.z) + cr2 * a2.z);
			atomicAdd(&particles[threadId + 2].v_a.x, cr1 * a1.x + (- (1. + cr2)) * a2.x);
			atomicAdd(&particles[threadId + 2].v_a.y, cr1 * a1.y + (- (1. + cr2)) * a2.y);
			atomicAdd(&particles[threadId + 2].v_a.z, cr1 * a1.z + (- (1. + cr2)) * a2.z);
			atomicAdd(&particles[threadId + 3].v_a.x, a2.x);
			atomicAdd(&particles[threadId + 3].v_a.y, a2.y);
			atomicAdd(&particles[threadId + 3].v_a.z, a2.z);

		} else {
			d_u[threadId] = 0.;
		}
    }


	if(blockId == 2) {
		if(threadId > (d_nmol - 1)) {

		}

	}


    // ============================================================================================================================ //

	__syncthreads();
	cudaDeviceSynchronize();

    if (blockSize >= 1024) {
        if (threadId < 512)
            d_u[threadId] += d_u[threadId + 512];
        __syncthreads();
    }
    if (blockSize >= 512) {
        if (threadId < 256)
            d_u[threadId] += d_u[threadId + 256];
        __syncthreads();
    }
    if (blockSize >= 256) {
        if (threadId < 128)
            d_u[threadId] += d_u[threadId + 128];
        __syncthreads();
    }
    if (blockSize >= 128) {
        if (threadId < 64)
            d_u[threadId] += d_u[threadId + 64];
        __syncthreads();
    }

    if (threadId < 32) {
        if (blockSize >= 64) d_u[threadId] += d_u[threadId + 32];
		__syncthreads();
        if (blockSize >= 32) d_u[threadId] += d_u[threadId + 16];
		__syncthreads();
        if (blockSize >= 16) d_u[threadId] += d_u[threadId + 8];
		__syncthreads();
        if (blockSize >= 8) d_u[threadId] += d_u[threadId + 4];
		__syncthreads();
        if (blockSize >= 4) d_u[threadId] += d_u[threadId + 2];
		__syncthreads();
        if (blockSize >= 2) d_u[threadId] += d_u[threadId + 1];

        if(threadId == 0) {
			switch (blockId) {
				case 0:
					*d_uB = d_u[0];
					break;
				case 1:
					*d_uT = d_u[0];
					break;
				case 2:
					*d_uLJ = d_u[0];

					break;
			}
		}
    }
}

__global__ void cudaUpdateVelocities(Particle *particles, double dt) {

	for(int i = threadIdx.x; i < d_nmol; i += blockDim.x) {
		particles[i].v_v.x += 0.5 * dt * particles[i].v_a.x;
		particles[i].v_v.y += 0.5 * dt * particles[i].v_a.y;
		particles[i].v_v.z += 0.5 * dt * particles[i].v_a.z;
	}
}

void step() {
    //int i, j;

    cudaMemcpy(d_particles, particles, n_mol * sizeof(Particle), cudaMemcpyHostToDevice);

	/***   Update positions of all atoms   ***/
	cudaUpdatePos <<< 1, blockSize >>> (d_particles, dt, LV);
    cudaDeviceSynchronize();

	cudaUpdateVelocities <<< 1, blockSize >>> (d_particles, dt);
    cudaDeviceSynchronize();

	cudaLJ <<< 3, blockSize, (blockSize + 1) * sizeof(double) >>> (d_particles, d_neighbour_list, d_sequence, LV, d_uT, d_uB, d_uLJ, d_uLJVector, i_step, n_mol);
	cudaDeviceSynchronize();
    cudaEnergy <<< 3, blockSize, (blockSize + 1) * sizeof(double) >>> (d_particles, d_neighbour_list, d_sequence, LV, d_uT, d_uB, d_uLJ, d_uLJVector, i_step);
    cudaDeviceSynchronize();

    

    cudaMemcpy(particles, d_particles, n_mol * sizeof(Particle), cudaMemcpyDeviceToHost);


	cudaMemcpy(&uBond, d_uB, sizeof(double), cudaMemcpyDeviceToHost);
	cudaMemcpy(&uTorsion, d_uT, sizeof(double), cudaMemcpyDeviceToHost);
	cudaMemcpy(&uLJ, d_uLJ, sizeof(double), cudaMemcpyDeviceToHost);
    cudaDeviceSynchronize();
}

__global__ void cudaShakeRelaxation(Particle *particles, Constraint *constraint, int n_c, double shake_cons_prec, double LV) {
	VectorR dr, dv;
	double cDev, cDevR, cDevV, g, ga;
	int changed, m1, m2, maxCycle;

	maxCycle = 200;
	cDevR = 0;
	cDevV = 0;



	for(int i = blockIdx.x; i < d_nmol; i += gridDim.x) {
		int nCycleR = 0;
		changed = 1;

		while(nCycleR < maxCycle && changed) {
			nCycleR++;
			changed = 0;
			cDev = 0;

			for(int m = threadIdx.x; m < n_c; m += blockDim.x) {
				m1 = constraint[m].ik;
				m2 = constraint[m].jk;

				dr.x = particles[m1].v_r.x - particles[m2].v_r.x;
				dr.y = particles[m1].v_r.y - particles[m2].v_r.y;
				dr.z = particles[m1].v_r.z - particles[m2].v_r.z;

				cudaVerifyBoundary2(&dr, LV);

				g = (sqr(dr.x) + sqr(dr.y) + sqr(dr.z) - 1.) / 4.;
				ga = fabs(g);
				cDev = (cDev > ga) ? cDev: ga;

				if(ga > shake_cons_prec) {
					changed = 1;
					atomicAdd(&particles[m1].v_r.x, -(g * dr.x));
					atomicAdd(&particles[m1].v_r.y, -(g * dr.y));
					atomicAdd(&particles[m1].v_r.z, -(g * dr.z));

					atomicAdd(&particles[m2].v_r.x, (g * dr.x));
					atomicAdd(&particles[m2].v_r.y, (g * dr.y));
					atomicAdd(&particles[m2].v_r.z, (g * dr.z));

				}
			}
		}
		cDevR = (cDev > cDevR) ? cDev: cDevR;

		int nCycleV = 0;
		changed = 1;

		while(nCycleV < maxCycle && changed) {
			nCycleV++;
			changed = 0;
			cDev = 0;
			for(int m = threadIdx.x; m < n_c; m += blockDim.x) {
				m1 = constraint[m].ik;
				m2 = constraint[m].jk;

				dr.x = particles[m1].v_r.x - particles[m2].v_r.x;
				dr.y = particles[m1].v_r.y - particles[m2].v_r.y;
				dr.z = particles[m1].v_r.z - particles[m2].v_r.z;

				cudaVerifyBoundary2(&dr, LV);

				dv.x = particles[m1].v_v.x - particles[m2].v_v.x;
				dv.y = particles[m1].v_v.y - particles[m2].v_v.y;
				dv.z = particles[m1].v_v.z - particles[m2].v_v.z;

				g = ((dr.x * dv.x) + (dr.y * dv.y) + (dr.z * dv.z)) / 2.;
				ga = fabs(g);
				cDev = (cDev > ga) ? cDev: ga;

				if(ga > shake_cons_prec) {
					changed = 1;
                    atomicAdd(&particles[m1].v_v.x, -(g * dr.x));
					atomicAdd(&particles[m1].v_v.y, -(g * dr.y));
					atomicAdd(&particles[m1].v_v.z, -(g * dr.z));
					atomicAdd(&particles[m2].v_v.x, (g * dr.x));
					atomicAdd(&particles[m2].v_v.y, (g * dr.y));
					atomicAdd(&particles[m2].v_v.z, (g * dr.z));
				}
			}
		}
		cDevV = (cDev > cDevV) ? cDev : cDevV;
	}
}

void shakeRelaxation() {
	VectorR dr, dv;
	double cDev, cDevR, cDevV, g, ga;
	int changed, m, m1, m2, maxCycle, i;

	maxCycle = 200;
	cDevR = 0;
	cDevV = 0;



	for(i = 0; i < n_mol; i++) {
		nCycleR = 0;
		changed = 1;

		while(nCycleR < maxCycle && changed) {
			nCycleR++;
			changed = 0;
			cDev = 0;

			for(m = 0; m < n_c; m++) {
				m1 = constraint[m].ik;
				m2 = constraint[m].jk;

				dr.x = particles[m1].v_r.x - particles[m2].v_r.x;
				dr.y = particles[m1].v_r.y - particles[m2].v_r.y;
				dr.z = particles[m1].v_r.z - particles[m2].v_r.z;

				verifyBoundary2(&dr);

				g = (sqr(dr.x) + sqr(dr.y) + sqr(dr.z) - 1.) / 4.;
				ga = fabs(g);
				cDev = (cDev > ga) ? cDev: ga;

				if(ga > shake_cons_prec) {
					changed = 1;
					particles[m1].v_r.x -= (g * dr.x);
					particles[m1].v_r.y -= (g * dr.y);
					particles[m1].v_r.z -= (g * dr.z);

					particles[m2].v_r.x += (g * dr.x);
					particles[m2].v_r.y += (g * dr.y);
					particles[m2].v_r.z += (g * dr.z);

				}
			}
		}

		cDevR = (cDev > cDevR) ? cDev: cDevR;

		nCycleV = 0;
		changed = 1;

		while(nCycleV < maxCycle && changed) {
			nCycleV++;
			changed = 0;
			cDev = 0;
			for(m = 0; m < n_c; m++) {
				m1 = constraint[m].ik;
				m2 = constraint[m].jk;

				dr.x = particles[m1].v_r.x - particles[m2].v_r.x;
				dr.y = particles[m1].v_r.y - particles[m2].v_r.y;
				dr.z = particles[m1].v_r.z - particles[m2].v_r.z;

				verifyBoundary2(&dr);

				dv.x = particles[m1].v_v.x - particles[m2].v_v.x;
				dv.y = particles[m1].v_v.y - particles[m2].v_v.y;
				dv.z = particles[m1].v_v.z - particles[m2].v_v.z;

				g = ((dr.x * dv.x) + (dr.y * dv.y) + (dr.z * dv.z)) / 2.;
				ga = fabs(g);
				cDev = (cDev > ga) ? cDev: ga;

				if(ga > shake_cons_prec) {
					changed = 1;
					particles[m1].v_v.x -= (g * dr.x);
					particles[m1].v_v.y -= (g * dr.y);
					particles[m1].v_v.z -= (g * dr.z);
					particles[m2].v_v.x += (g * dr.x);
					particles[m2].v_v.y += (g * dr.y);
					particles[m2].v_v.z += (g * dr.z);
				}
			}
		}

		cDevV = (cDev > cDevV) ? cDev : cDevV;
	}
}

void berendsenThermostat() {
	double sum, lambda, temp;
	int i;

	sum = 0;
	for(i = 0; i < n_mol; i++) {
		sum += sqr(particles[i].v_v.x) + sqr(particles[i].v_v.y) + sqr(particles[i].v_v.z);
	}

	temp = sum / (3. * (n_mol - 1));
	lambda = sqrt(1. + (dt / c_T) * (temperature / temp - 1.));

	for(i = 0; i < n_mol; i++) {
		particles[i].v_v.x *= lambda;
		particles[i].v_v.y *= lambda;
		particles[i].v_v.z *= lambda;
	}
}

void calcStatus() {
	int i;
	double sum;

	sum = 0.;
	for(i = 0; i < n_mol; i++) {
		sum += sqr(particles[i].v_v.x) + sqr(particles[i].v_v.y) + sqr(particles[i].v_v.z);
	}

	kinetic_energy = 0.5 * sum;
	current_temperature = sum / (3. * (n_mol - 1));
	density = n_mol / (LV * sqr(LV));

	uSum = uTorsion + uBond + uLJ;
	total_energy = uSum + kinetic_energy;

	sum = 0.;
	for(i = 0; i < n_mol - 1; i++) {
		sum += sqrt(sqr(particles[i].v_r.x - particles[i + 1].v_r.x) + sqr(particles[i].v_r.y - particles[i + 1].v_r.y) + sqr(particles[i].v_r.z - particles[i + 1].v_r.z));
 	}
 	bond_avg = sum / (n_mol - 1);
}

void calcRG() {
	int i;
	int h = 0, p = 0;
	double x_avg = 0., y_avg = 0., z_avg = 0.;
	double x_avgH = 0., y_avgH = 0., z_avgH = 0.;
	double x_avgP = 0., y_avgP = 0., z_avgP = 0.;

	rGH = 0.;
	rGP = 0.;
	rG = 0.;
	for(i = 0; i < n_mol; i++) {
		x_avg += particles[i].v_r.x;
		y_avg += particles[i].v_r.y;
		z_avg += particles[i].v_r.z;

		if(sequence[i] == 'A') {
			x_avgH += particles[i].v_r.x;
			y_avgH += particles[i].v_r.y;
			z_avgH += particles[i].v_r.z;

			h++;
		} else {
			x_avgP += particles[i].v_r.x;
			y_avgP += particles[i].v_r.y;
			z_avgP += particles[i].v_r.z;

			p++;
		}
	}

	x_avg = x_avg / n_mol;
	y_avg = y_avg / n_mol;
	z_avg = z_avg / n_mol;

	x_avgH = x_avgH / h;
	y_avgH = y_avgH / h;
	z_avgH = z_avgH / h;

	x_avgP = x_avgP / p;
	y_avgP = y_avgP / p;
	z_avgP = z_avgP / p;

	for(i = 0; i < n_mol; i++) {
		rG += sqr(particles[i].v_r.x - x_avg) + sqr(particles[i].v_r.y - y_avg) + sqr(particles[i].v_r.z - z_avg);

		if(sequence[i] == 'A')
			rGH += sqr(particles[i].v_r.x - x_avgH) + sqr(particles[i].v_r.y - y_avgH) + sqr(particles[i].v_r.z - z_avgH);
		else
			rGP += sqr(particles[i].v_r.x - x_avgP) + sqr(particles[i].v_r.y - y_avgP) + sqr(particles[i].v_r.z - z_avgP);
	}

	rG = sqrt(rG / n_mol);
	rGH = sqrt(rGH / h);
	rGP = sqrt(rGP / p);
}

void evaluate() {
	int i;

	calcStatus();

	if(bond_avg <= 1.) {
		if(i_step == 0) {
			best_potencial_energy = uSum;
			best_step = i_step;
		} else {
			if(uSum < best_potencial_energy) {
				best_potencial_energy = uSum;
				best_step = i_step;

				for(i = 0; i < n_mol; i++) {
					best_structure[i].v_r.x = particles[i].v_r.x;
					best_structure[i].v_r.y = particles[i].v_r.y;
					best_structure[i].v_r.z = particles[i].v_r.z;
				}
			}
		}
	}

	calcRG();
}

void calcCenterMass() {
	int i;
	VectorR r;

	r.x = 0.;
	r.y = 0.;
	r.z = 0.;
	for(i = 0; i < n_mol; i++) {
		r.x += particles[i].v_r.x;
		r.y += particles[i].v_r.y;
		r.z += particles[i].v_r.z;
	}

	center_mass.x = r.x / n_mol;
	center_mass.y = r.y / n_mol;
	center_mass.z = r.z / n_mol;
}

void printSummary(char **argv) {
	FILE *fo;
    char file_name[200];

    //cudaMemcpy(particles, d_particles, n_mol * sizeof(Particle), cudaMemcpyDeviceToHost);

	calcCenterMass();
	system("clear");
	printf("Step: %d\n\nTemp = %lf\nTotal Lennard-Jones Potential = %lf\nTotal Torsion Potential = %lf\nTotal Chain Angle Potential = %lf\nTotal Potential Energy = %lf\n", i_step, current_temperature, uLJ, uTorsion, uBond, uSum);
	printf("Center of Mass = (%.2lf, %.2lf, %.2lf)\n", center_mass.x, center_mass.y, center_mass.z);
	printf("Radius of Gyration - Hydrophobic = %lf\nRadius of Gyration - Polar = %lf\nRadius of Gyration - All = %lf\n", rGH, rGP, rG);
	printf("Bond Length Average = %lf\n", bond_avg);

	if(print_summary2file == 'y' || print_summary2file == 'Y') {
		sprintf(file_name, "%s_summary_%s.txt", argv[2], argv[3]);
        fo = fopen(file_name, "a+");
		if(i_step == 0)
			fprintf(fo, "Step\tTemperature\tU_LJ\tU_Torsion\tU_ChainAngles\tU_Total\trGH\trGP\trG\n");

		fprintf(fo, "%d\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\n", i_step, temperature, uLJ, uTorsion, uBond, uSum, rGH, rGP, rG);
		fclose(fo);
	}
}


void savePathways(char **argv) {
	FILE *fo;
	int i;
    char file_name[200];


    sprintf(file_name, "%s_%s.txt", argv[2], argv[3]);
	fo = fopen(file_name, "a+");

   	fprintf(fo, "N\tx\ty\tz\n");
	for (i = 0; i < n_mol; i++) {
		fprintf(fo, "%d\t%lf\t%lf\t%lf\n", i, particles[i].v_r.x, particles[i].v_r.y, particles[i].v_r.z);
	}

	fprintf(fo, "\n\nPotential Energy = %lf\nStep = %d\n", uSum, i_step);
	fprintf(fo, "uLJ = %lf\nTorsion = %lf\nBond = %lf\n",  uLJ, uTorsion, uBond); 
	fprintf(fo, "rGAll = %lf\nrGH = %lf\nrGP = %lf\n\n\n", rG, rGH, rGP);

	fclose(fo);

}

void freeArrays() {
	free(particles);
	free(best_structure);
	free(constraint);
	free((void *)mMat[0]);
	free((void *)mMat);
}

void freeDevice() {
    cudaFree(d_particles);
    cudaFree(d_sequence);
    cudaFree(d_neighbour_list);
	cudaFree(d_uB);
	cudaFree(d_uT);
	cudaFree(d_uLJ);
    cudaFree(d_uLJVector);
}



void finishSim(char **argv) {
	finTimer();
	savePathways(argv);
	freeArrays();
    freeDevice();
}
