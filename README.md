# protein_folding_datasets

 This work reports datasets for the study of protein folding dynamics, corresponding to the spatiotemporal data of folding trajectories obtained by an in silico method. Four datasets are described: the first one is based on the Fibonacci sequence (13FIBO) and the other three on real biological sequences (**2GB1*,**1PLC** and **5NAZ**). Each dataset is composed of **1000 different trajectories data**, each containing structural information of the protein during the folding process including the spatial coordinates of each amino acid at each time step, and the free energy and radius of gyration values for each structure.

An animation of a 1PCL protein folding path is shown below.
![example_dataset](https://github.com/bioinfolabic/protein_folding_datasets/blob/master/Images/gif_exemplo.gif)


The programs used to generate the datasets and charts for the analyzes are described below.

## Path DM_cesar
### ▪ Sub-path create_datasets: 

This folder contains the program used to generate the data of folding trajectories using Molecular Dynamics (*define.h*,*function_MD.c* and *main.c*), which consist of files in the C ++ language. It is worth mentioning that there is no need to modify these codes because the variables are all parameterized in input files separately from the DM code.

### ▪ executa0.sh: 

The script *executa0.sh* allows to generate 1000 different data of folding trajectories, as shown below:
```
#!/bin/bash

gcc -o executavel main.c func_MD.c -lm
for i in {1..1001}
do
   ./executavel md_test_2gb1_albertsclassification.in $i
done
```
where the first line compiles the DM code and each loop generates one folding trajectory.
the input file "md_test_2gb1_albertsclassification.in" is the protein data, and "$ i" represents the ID of that structure, being able to assume values from 1 to 1000. The program saves the data in a folder whose name is the number of amino acids of the protein. To run the program:

program                         |language |where to run                |running the program                  |virtual machine        |
--------------------------------|---------|----------------------------|-------------------------------------|-----------------------|
create_datasets                 |C++      |linux terminal              |./executa0.sh                        |-| 


#### ▪▪ md_test_2gb1_alberts_classification.in:
Input file for the previous program. For each protein, the following variables were modified:

 Protein  |  Sequence | nMol | ProtLen |  LV  |  nC  |
--------- | ----------|------|---------|------|------|
13FIB     | see paper |  13  |    13   |  26  |  12  |
2GB1      | see paper |  56  |    56   |  112 |  55  |
1PCL      | see paper |  99  |    99   |  198 |  98  |

where **Sequence** is the AB sequence of the protein, **nMol** and ProtLen are the number of amino acids (the size of the protein), **LV** is the size of the box edge where the protein is contained and **nC** is the number of amino acids minus one.


The following figure represents one of the program output files, which contains the beginning of one pathway.

![example_dataset](https://github.com/bioinfolabic/protein_folding_datasets/blob/master/Images/format_dataset.png)


### ▪ old_versions:
Old versions of the "create_datasets" program.
### ▪ pathways_test:
Examples of output from the "create_datasets" program using a fibonacci sequence. 





## Path Images_folding
### ▪ pathway_print_multi-subplot.py:
A python program that creates images of the structure folding path. The program saves the images in .png format, where an image is ploted for each state of pathway, so if the folding path has 1000 states, 1000 images will be created. Before you run the program,  replace the following variables with the desired input files:
```
path_pathways = '/home/bruna/teste/13_FIBO/' # folder in which the datasets are
filename = 'pathways13_999.txt'              # name of the dataset from which the images will be created
filesequencia = 'seq_13.txt'                 # file containing the AB sequence of the protein
path_save = 'img_13'                         # folder where the images will be saved
```
Running the program:

program                         |language |where to run                |running the program                  |virtual machine        |
--------------------------------|---------|----------------------------|-------------------------------------|-----------------------|
pathway_print_multi-subplot.py  |python 2.7|visual stidio code terminal |python pathway_print_multi-subplot.py|virtualenv: imageio, matplotlib|

An example of the image produced by the program is shown below.

![example_dataset](https://github.com/bioinfolabic/protein_folding_datasets/blob/master/Images/exemplo_img_56_1000.png)

To create a protein folding video, the name of the image files must follow a numerical sequence, as shown in the example:

![example_file_image](https://github.com/bioinfolabic/protein_folding_datasets/blob/master/Images/exemplo_arquivo_imagem.png)

with the images files, open the terminal in the folder, and enter  the follow command:

```
ffmpeg -r 3 -f image2 -s 720x480 -start_number 0 -i %d.png -vframes 1000 -vcodec libx264 -crf 25 -pix_fmt yuv420p folding.mp4
```

where 3 the number of frames per second, 720x480 the resolution, "%d.png" the name of the images, where "%d" assumes values between 0 and 1000, and "folding.mp4" is the resulting video.


### ▪ Sub-Path images:
Folder containing the examples of images obtained through the previous program.
### ▪ remove_lines.py:
Python program that removes all the contents of a file from a specific line. This program is very useful if you run the "create_datasets" program twice and the pathways override.





## Path Programa Heamap Kabsch RMSD
An analysis was accomplished to observe the spatial differences between the initial and final structures of the folding trajectories using heatmaps. For each dataset, all the 1,000 protein structures were compared to each other, at the beginning and the end of the trajectories. The comparison of pairs of structures was accomplished with the Kabsch algorithm, normalized between zero and one. Then, computed values were the plot in a heatmap.

### ▪ dataset_heatmap_kabsch.py:

This program creates a .npy matrix of spatial differences between the initial and final structures of the folding trajectories. Each point of the horizontal and vertical axes of the matrix represents a protein structure at a given point of the trajectory (either initial of final point). The higher number, the closer to 1 value it is in the Kabsch scale,  meaning that the structures tend to be more different. The opposite holds, meaning similar spatial structures.
Before run the program, change the following data by the desired input files:
```
number_pathway = 1000                                                 # number of datasets / pathways
path_pathways = '/home/bruna/teste/13FIBO/'                           # file folder
filename = 'ABBABBABABBAB_pathwaystep3000_'+ str(i) +'_pathways.txt'  # name of datasets / pathways
```
Running the program:

  program                         |language |where to run              |running the program                  |virtual machine        |
--------------------------------|---------|----------------------------|-------------------------------------|-----------------------|
dataset_heatmap_kabsch.py       |python 2.7|visual stidio code terminal |python dataset_heatmap_kabsch.py     |virtualenv: numpy|

### ▪ generate_dataset:
Examples of outputs from the previous program.
### ▪ other programs
The other programs are required for dataset_heatmap_kabsch.py to work.



## Path Programs
### ▪ rg.py:
For each pathway this program calculates the values of the RgAll, RgH, and RgP, and saves it to a .txt file. Before run the program, change the following data by the desired input files:
```
n_arquivos = 1000                                                               # number of pathways
s = "ABBABBABABBAB"                                                             # AB protein sequence
with open("/home/bruna/heatmap/13_FIBO_rG/rg"+str(i)+".txt", "w") as saida:     # output file, file containing the new rGs
with open("/home/bruna/heatmap/13_FIBO/pathways13_"+str(i)+".txt", "r") as f:   # pathways from which the rG will be calculated
```
Running the program:

program                         |language |where to run                |running the program                  |virtual machine        |
--------------------------------|---------|----------------------------|-------------------------------------|-----------------------|
rg.py                           |python 2.7|visual stidio code terminal |python rg.py                         |-|

The following image represents the output of the program:

![exemplo_rg](https://github.com/bioinfolabic/protein_folding_datasets/blob/master/Images/exemplo_rg.png)


### ▪ graphics_rg.py:
This program calculates the average Rg (RgH, RgP and RgALL) of all pathways, then plots a graph of the rG values per iteration. Before running the program, change the input file:
```
with open("/home/bruna/Downloads/13_fibonacci/pathways13_"+str(i)+".txt", "r") as input: # file generated by the program "create datasets"
```
Running the program:

program             |language |where to run                |running the program                  |virtual machine        |
--------------------|---------|----------------------------|-------------------------------------|-----------------------|
graphics_rg.py|python 2.7|visual stidio code terminal |python graphics_rg.py   |virtualenv: numpy, matplotlib, seaborn|

The following image represents the output of the program:
![example_rG](https://github.com/bioinfolabic/protein_folding_datasets/blob/master/Images/average_Rg_1PCL.png)

### ▪ graphics_ep.py:
This program calculates the average potential energy value of all pathways of three proteins (13FIBO, 2GB1 and 1PLC) and plots those values in a graph of the potential energy of each protein by step. Before running the program, change the input files:
```
with open("/home/bruna/heatmap/13FIBO/13_fibonacci_"+str(i+1)+".txt", "r") as input: # pathways of the first protein (13FIB)
with open("/home/bruna/heatmap/2GB1/pathways56_"+str(i)+".txt", "r") as input:       # pathways of the second protein (2GB1)
with open("/home/bruna/heatmap/1PCL_99/pathways99_"+str(i)+".txt", "r") as input:    # pathways of the third protein (1PLC)
with open("/home/bruna/heatmap/5NAZ_229/pathways229_"+str(i)+".txt", "r") as input:    # pathways of the fourth protein (5NAZ)
```
Running the program:

program             |language |where to run                |running the program                  |virtual machine        |
--------------------|---------|----------------------------|-------------------------------------|-----------------------|
graphics_ep.py|python 2.7|visual stidio code terminal |python graphics_ep.py   |virtualenv: numpy, matplotlib, seaborn|

The following image represents the output of the program:
![example_Ep](https://github.com/bioinfolabic/protein_folding_datasets/blob/master/Images/average_Ep%20(1).png)

### ▪ graphics_heatmap.py:
This program uses the .npy array generated by the program "dataset_heatmap_kabsch.py" to create two heatmaps, one comparing the initial structures of the pathways and another comparing the final structures.
To create the initial heatmap, comment on the part of the program that says "HEATMAP END" and change the input file:
```
a = np.load('/home/bruna/heatmap/generate_dataset/heatmap_13_fibonacci_begin.npy') # file generated by the program "dataset_heatmap_kabsch.py"
```
To create the final heatmap, comment on the part of the program that says "HEATMAP BEGIN" and change the input file:
```
b = np.load('/home/bruna/heatmap/generate_dataset/heatmap_13_fibonacci_end.npy') # file generated by the program 
```
Running the program:

program             |language |where to run                |running the program                  |virtual machine        |
--------------------|---------|----------------------------|-------------------------------------|-----------------------|
graphics_heatmap.py|python 2.7|visual stidio code terminal|python graphics_heatmap.py|virtualenv: numpy, matplotlib, seaborn|

The following image represents the output of the program:

![example_heatmap](https://github.com/bioinfolabic/protein_folding_datasets/blob/master/Images/heatmap_1PCY_begin.png)


## Path Results
Contains all datasets and graphs generated during this work.


## All data that needs to be changed to place the input files are indicated in the program with the following comment:
```
#######################################################################################################
########################################### CHANGE HERE ###############################################
#######################################################################################################
```
or
```
########################################################################################### CHANGE HERE
```
