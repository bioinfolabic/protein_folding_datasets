all:
	nvcc --gpu-architecture=compute_61 --device-c main.c functions.cu
	nvcc --gpu-architecture=compute_61 main.o functions.o

run_13FIBO:
	./a.out Proteinas/13FIBO/13FIBO_13.in Proteinas/13FIBO/13FIBO_13_pathway 00 > 13FIBO_13.txt 99 0; \

run_2GB1:
	./a.out Proteinas/2GB1/2GB1_56.in Proteinas/2GB1/2GB1_56_pathway 00 > 2GB1_56.txt 0 0; \

run_1PLC:
	./a.out Proteinas/1PLC/1PLC_99.in Proteinas/1PLC/1PLC_99_pathway 00 > 1PLC_99.txt 0 0; \

run_5NAZ:
	./a.out Proteinas/5NAZ/5NAZ_229.in Proteinas/5NAZ/5NAZ_229_pathway 00 > 5NAZ_229.txt 7 0; \

clean:
	rm a.out
	rm main.o
	rm functions.o
