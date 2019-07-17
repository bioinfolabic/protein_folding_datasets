# Copyright (c) 2019, Universidade Tecnológica Federal do Paraná - Bioinformatics and Computacional Intelligence Laboratory
#All rights reserved.
#
#The 3-clause BSD license is applied to this software.
ab:
	cd SRC && python AB_sequence.py 2gb1 

gpu:
	cd SRC && make -f Makefile_PATHMOLD-AB

cpu:
	cd SRC_CPU && make -f Makefile_CPU 

run_gpu:
	cd SRC && ./a.out ../INPUT/2GB1_56.in ../OUTPUT/pathways56 0 0

run_cpu:
	cd SRC_CPU && ./a.out ../INPUT/2GB1_56.in 0 && mv pathways56_0.txt ../OUTPUT/pathways56_0.txt

visualize:
	cd SRC && python pathway_print_multi-subplot.py && mv folding.mp4 ../OUTPUT/folding.mp4 && mv img ../OUTPUT/


