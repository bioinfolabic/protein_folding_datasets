# protein_folding_datasets

![example_dataset](https://github.com/bioinfolabic/protein_folding_datasets/blob/master/Images/gif_exemplo.gif)

  This work presents datasets for the Protein Folding Problem (PFP), they consist in data of folding trajectories obtained by an in silico method.
  Three datasets are introduced in this paper, where one is based on the Fibonacci sequence (13FIBO) and the other two on biological sequences (2GB1 and 1PCY). Each dataset is composed of 1000 different trajectories data, which contains information on the structures of the protein during the folding process. In addition, it is also present the energy and radius rotation values for each structure.
All datasets proposed in this work are available in https://mega.nz/#F!O4wBHQiB!mpM81jK9ycbSKPvdP57avw




## DM cesar
### ▪ create_datasets: 
Within this part there is a C ++ program which, through the TAL input file, manufactures 1000 different initial structures for each protein and calculates its folding paths (the datasets), as exemplified in the image below, which contains the beginning of a dataset of one of the 1000 initial structures of a protein.

![example_dataset](https://github.com/bioinfolabic/protein_folding_datasets/blob/master/Images/format_dataset.png)

### ▪ md_test_2gb1_alberts classification:
Input file for the previous program. For each protein, the following variables were modified:

proteína  | sequência |nMol|ProtLen|  LV |
--------- | ----------|----|-------|-----|
13FIB     | see paper | 13 |   13  |  26 |
2GB1      | see paper | 56 |   56  | 112 |
1PCY      | see paper | 99 |   99  | 198 |

### ▪ old_versions:
Old versions of the "create_datasets" program.
### ▪ executa0.sh:
Old versions of the "create_datasets" program.
### ▪ pathways_test:
Examples of output from the "create_datasets" program using a fibonacci sequence. 





## Images_folding
### ▪ pathway_print_multi-subplot.py:
A python program that creates images of the structure's folding path from a dataset. The program variables are the protein size, the "AB" sequence and the folding path of the structure.  The program saves the images in .png format in a folder, being an image for each configuration of the structure, so if the folding path has 1000 configurations, 1000 images will be made. Below is an example of the image produced by the program.

![example_dataset](https://github.com/bioinfolabic/protein_folding_datasets/blob/master/Images/exemplo_img_56_1000.png)

### ▪ images:
Folder containing the examples of images obtained through the previous program.
### ▪ remove_lines.py:
Python program that removes all the contents of a file from a specific line. This program is very useful if you run the "create_datasets" program twice and the pathways override.





## Program Heamap Kabsch RMSD
### ▪ dataset_heatmap_kabsch.py:
Python program that, starting from the 1000 datasets of a protein, creates a .npy array from which the initial and final heatmap will be created.
### ▪ generate_dataset:
Examples of outputs from the previous program.





## Programs
### ▪ rg.py:
For each dataset it calculates all the spinning spokes and saves it to a .txt file. This program is useful because it produces more accurate results.
### ▪ switch_rg.py:
Cria novos arquivos .txt usando os datasets das proteínas mas com os raios de giração calculados pelo programa anterior.
### ▪ graphics.py:
This program is responsible for creating the initial and final heatmaps of each protein, as well as creating the graph of the potential energy per step and the graph of the turn-by-step rays. Below are examples of the heatmap, graph of potential energy and graph of the spinning rays respectively.

![example_heatmap](https://github.com/bioinfolabic/protein_folding_datasets/blob/master/Images/heatmap_1PCY_begin.png)

![example_Ep](https://github.com/bioinfolabic/protein_folding_datasets/blob/master/Images/average_Ep.png)

![example_rG](https://github.com/bioinfolabic/protein_folding_datasets/blob/master/Images/average_Rg_1PCY.png)





## Results
Contains all datasets and graphs generated during this work.
