/* 1 - Iniciar temporizador
 * 2 - Ler parametros
 * 3 - Iniciar MD
 * 4 - Finalizar MD e temporizador
 */
#include "functions.h"


int main(int argc, char **argv) {

	initTimer();
	int i_seed = (int)*argv[4] - 48;
	printf("i_seed = %d\n", i_seed);
	srand(time(NULL));
	seedMT();

	int device = (int)*argv[5] - 48;
	cudaSetDevice(device);
	printf("GPU: %d\n", device);

    if(loadFile(argv) == SUCCESS) {
        putParameters();

   		/***   Initiallize atoms' status and LJ potential tables   ***/
    	/***   Set parameters controlling the simulation   ***/
    	initMD();

		/***   Iterations   ***/
		for(i_step = 0; i_step < step_limit; i_step++) {
		//for(i_step = 0; i_step < 15; i_step++) {
			step();

			if(i_step % shake_step2shake == 0)
				shakeRelaxation();

			berendsenThermostat();

			evaluate();

			/* APENAS SALVAR OS PATHWAYS PARA PEGAR O TEMPO DOS PATHWAYS
			if((print_summary == 'y' || print_summary == 'Y') && (i_step % print_summary_interval == 0))
				printSummary(argv);

			if((report_file == 'y' || report_file == 'Y') && (i_step % step2report == 0))
				report(argv);
			*/

			if((save_pathways == 'y' || save_pathways == 'Y') && (i_step % pathways_step == 0) && (i_step != 0) || (i_step == 1) )
				savePathways(argv);

		}

		finishSim(argv);

    } else {
        printf("Error: could not load file.");
    }
}
