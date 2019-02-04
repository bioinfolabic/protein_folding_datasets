# protein_folding_datasets

![example_dataset](https://github.com/bioinfolabic/protein_folding_datasets/blob/master/Images/gif_exemplo.gif)

  This work presents datasets for the Protein Folding Problem (PFP), they consist in data of folding trajectories obtained by an in silico method.
  Three datasets are introduced in this paper, where one is based on the Fibonacci sequence (13FIBO) and the other two on biological sequences (2GB1 and 1PCY). Each dataset is composed of 1000 different trajectories data, which contains information on the structures of the protein during the folding process. In addition, it is also present the energy and radius rotation values for each structure.
All datasets proposed in this work are available in https://mega.nz/#F!O4wBHQiB!mpM81jK9ycbSKPvdP57avw




## DM cesar
### ▪ create_datasets: 
Within this part there is a C ++ program which, through an input file, creates 1000 different initial structures for each protein and calculates its folding paths (the datasets).
To create a thousand different initial structures for each protein was used a random function initialized from a "seed" inserted in the program generated initial coodenates creating new structures. This seed is the identification number (ID) of each new structure. To run the program with a thousand different IDs an .sh file (executa0.sh) was created, as follows:

```
#!/bin/bash

for i in {1..1001}
do
   gcc -o executavel main.c func_MD.c -lm
   ./executavel md_test_2gb1_albertsclassification.in $i
done
```

Being "md_test_2gb1_albertsclassification.in" the input file with the protein data, and "$ i" represents the ID of that structure, being able to assume values from 1 to 1000.
The following figure represents one of the program output files, which contains the beginning of a dataset of one of the 1000 initial structures of a protein.

![example_dataset](https://github.com/bioinfolabic/protein_folding_datasets/blob/master/Images/format_dataset.png)

### ▪ md_test_2gb1_alberts classification:
Input file for the previous program. For each protein, the following variables were modified:

 Protein  |  Sequence | nMol | ProtLen |  LV  |  nC  |
--------- | ----------|------|---------|------|------|
13FIB     | see paper |  13  |    13   |  26  |  12  |
2GB1      | see paper |  56  |    56   |  112 |  55  |
1PCY      | see paper |  99  |    99   |  198 |  98  |

Being "Sequence" the AB sequence of the protein, nMol and ProtLen the number of amino acids (the size of the protein), LV twice the size of it and nC is the number of amino acids minus one.
### ▪ old_versions:
Old versions of the "create_datasets" program.
### ▪ pathways_test:
Examples of output from the "create_datasets" program using a fibonacci sequence. 





## Images_folding
### ▪ pathway_print_multi-subplot.py:
A python program that creates images of the structure's folding path from a dataset. The program saves the images in .png format, being an image for each configuration of the structure, so if the folding path has 1000 configurations, 1000 images will be made. Below is an example of the image produced by the program.

![example_dataset](https://github.com/bioinfolabic/protein_folding_datasets/blob/master/Images/exemplo_img_56_1000.png)

To run the program simply replace the variables indicated in it, they are:

```
path_pathways = '/home/bruna/DM_leandro/pathways/1PCY/'       # folder in which the datasets are
filename = 'pathways99_1000.txt'                              # name of the dataset from which the images will be created
filesequencia = 'seq_99.txt'                                  # file containing the AB sequence of the protein
path_save = 'imagens/img99/'                                  # folder where the images will be saved
```

To create a protein folding video, the name of the image files must follow a numerical sequence, as shown in the example:

![example_file_image](https://github.com/bioinfolabic/protein_folding_datasets/blob/master/Images/exemplo_arquivo_imagem.png)

Having the images, using linux just open the terminal in the folder containing the images and enter the command:

```
ffmpeg -r 3 -f image2 -s 720x480 -start_number 0 -i% d.png -vframes 1000 -vcodec libx264 -crf 25 -pix_fmt yuv420p folding.mp4
```

Being 3 the number of frames per second, 720x480 the resolution, "%d.png" the name of the images, where "%d" assumes values between 0 and 1000, and "folding.mp4" is the resulting video.


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
