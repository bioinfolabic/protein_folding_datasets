#include "func_MD.h"


int main (int argc, char **argv)
{
	sscanf(argv[2],"%d",&ID);
   
	int finish, k = 0;
	srand(time(0));

	ti = time(NULL);				//initial time (used to calculate the processing time)


	initial_message();
	GetParameters(argv);			// lê e coleta os dado do arquivo inicial
	SetParams();					// faz algumas contas de velocidade e temperatura
	PrintParameters();				// apenas da print nos parâmetros coletados do arquivo e calculados, não é pra atrasar o programa
	Initialize_MD();				// inicializa os dados dos aminoácidos da proteína

	finish = 0;
	while (finish == 0)
	{
		finish = Step();		// faz todas as contas para calcular as posiçẽos de cada dobramento
		k++;
	}	//main function MD
	// time_ex();
	freeArrays();					// apenas libera tudo o que foi alocado no programa
	
	return 1;
}

