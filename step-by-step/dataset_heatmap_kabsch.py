# -*- coding: utf-8 -*-


import csv
import pickle
import StringIO
import re
import numpy as np
from importlib import import_module
import kabsch_algorithm as kabsch
#kabsch = import_module("metrics.kabsch_algorithm")
import sys


def open_file(filename):
	f = open(filename, 'r')
	return f

def close_file(f):
	f.close()

def read_file(f):
	count = 0
	#dataset = []
	data		 = [] 
	vector_ep	= []
	vector_rgh   = []
	vector_rgp   = []
	vector_rgall = []
	f.readline()
	str_ep	= "Potential Energy = "
	str_rgh   = "rGH = "
	str_rgp   = "rGP = "
	str_rgall = "rGAll = "
	for lines_read, line in enumerate(f): # lines_read - the number of the line
										  # line	   - the content in the current lines_read variable
		test  = re.findall(r'[a-zA-Z]+', line) #Check if there are letters  
		test2 = re.findall(r'[0-9]+', line)	#Check if there are numbers   
		# Begin Test
		#print("Lineread: ",lines_read)
		#print("Test letters content: ", test)
		#print("Test letters size ", len(test))
		#print("Test number contant: ", test2) 
		#print("Test number size ", len(test2)) 
		# End Test  		
		#print type(line)

		ep	  = line.find(str_ep)
		rgh   = line.find(str_rgh)
		rgp   = line.find(str_rgp)
		rgall = line.find(str_rgall)

		if(len(test) == 0 and len(test2) > 0): # copy the line that has the coord information   
			aux = (np.genfromtxt(StringIO.StringIO(line), delimiter="\t"))
			data = np.append(aux,data)
		elif(ep != -1):
			aux = float(line.replace(str_ep,''))
			vector_ep.append(aux)
		elif(rgh != -1):
			aux = float(line.replace(str_rgh,''))
			vector_rgh.append(aux)
		elif(rgp != -1):
			aux = float(line.replace(str_rgp,''))
			vector_rgp.append(aux)
		elif(rgall != -1):
			aux = float(line.replace(str_rgall,''))
			vector_rgall.append(aux)
		else:
			aux = 0

	return np.array(data), np.asarray(vector_ep), np.asarray(vector_rgh), np.asarray(vector_rgp), np.asarray(vector_rgall)


def format_dataset(dataset, shape):			 # dataset = raw dataset, shape = (numere of sample, sequence, coord x-y-z)
	#print np.array(dataset).shape

	dataset = dataset.reshape((shape))	   # make the reshape 
	dataset = dataset[:,::-1,:]			  # invert to first aminio acid to the last
	dataset = dataset[::-1,:,:]			  # invert to the first fold to the last fold
	dataset = dataset[:,:,1:4]			   # remove the amino acid position 
	return dataset

def generate_dataset_label(dataset):
	amount_data = dataset.shape[0]
	label = dataset[1:amount_data,:,:]
	dataset = dataset[0:(amount_data-1),:,:]
	return dataset, label

def save_pickl(data, namefile):
	pickle.dump( data, open( namefile + ".p", "wb" ))

def save_dataset(dataset, label, namefile, path_save):
	save_pickl(dataset, path_save + namefile)
	save_pickl(label  , path_save + "label_"+namefile)

def get_shape(dataset):
	sequence_size = int(dataset[0]) + 1
	dimension	 = 4 # amino acid sequence; x; y; z	 
	total_element = dataset.shape[0]
	folding	   = total_element / (sequence_size * dimension) 
	return [folding, sequence_size, dimension]

def save_heatmap(dataset_name, kabsch_matrix):
	np.save(dataset_name, kabsch_matrix)

def create_kabsch_matrix(all_structure, number_pathway):
	#kabsch = import_module("metrics.kabsch_algorithm")      # it is need to create a __init__.py and compile it
	kabsch_matrix = np.ndarray(shape=(number_pathway,number_pathway), dtype=float)

	# Calcular o kabasch entre as estruturas
	for i in xrange(all_structure.shape[0]):
		for j in xrange(all_structure.shape[0]):
			print all_structure[i].shape
			kabsch_matrix[i,j] = kabsch.get_kabsch(all_structure[i],all_structure[j])
	return kabsch_matrix

def	mean_pathway_kabsch(total_structure_begin):
	pathway_amount  = total_structure_begin.shape[0]
	fold_amount     = total_structure_begin.shape[1]
	print "reading all pathway..."
	mean_path = []
	for i in xrange(fold_amount):
		all_structure = []
		for j in xrange(pathway_amount):
			all_structure.append(total_structure_begin[j,i])
		all_structure = np.asarray(all_structure)
		kabsch_matrix = create_kabsch_matrix(all_structure)
		mean_path.append(kabsch_matrix)
	mean_path  = np.asarray(mean_path)
	print mean_path.shape
	print np.mean(mean_path,axis=0).shape
	return np.mean(mean_path,axis=0)

def main():
	all_structure = []
	all_structure_begin = []
	total_structure_begin = []
	all_ep        = []
	all_rgh 	  = []
	all_rgp 	  = []
	all_rgall 	  = []

	####################################################################################################################################################################################
	#################################################################################### CHANGE HERE ###################################################################################
	####################################################################################################################################################################################
	number_pathway = 500
	for i in xrange(1,number_pathway+1): 
		path_pathways = '/home/bruna/DM/Proteinas/5NAZ/5NAZ/' # file folder
		filename      = '5NAZ_229_pathway_'+ str(i) +'.txt' # name of each pathway
		################################################################################################################################################################################
		print i
		path_save     = 'datasets/' # destination folder to save

		f = open_file(path_pathways + filename)
		dataset, vector_ep, vector_rgh, vector_rgp, vector_rgall = read_file(f)
		close_file(f)

		shape 		= get_shape(dataset)
		dataset		= format_dataset(dataset, shape)

		all_structure_begin.append(dataset[0])
		all_structure.append(dataset[dataset.shape[0]-1])
		print vector_ep[dataset.shape[0]-1]
		all_ep.append(vector_ep[dataset.shape[0]-1])
		all_rgh.append(vector_rgh[dataset.shape[0]-1])
		all_rgp.append(vector_rgp[dataset.shape[0]-1])
		all_rgall.append(vector_rgall[dataset.shape[0]-1])
		
		total_structure_begin.append(dataset)

		print filename
		#print "shape: ", shape
		#print "new shape: ", dataset.shape
		#print "energy p: ", vector_ep.shape
		#dataset, label = generate_dataset_label(dataset)
		#dataset_name   = filename[0:filename.find(".")]
		#save_dataset(dataset, label, dataset_name, path_save)
		#print ("\ndataset saved! ", dataset_name)
		#print "\n"
		#print ("\ndataset shape: ",dataset.shape)
		#print ("dataset sequence shape: ", dataset[0].shape)
		#print ("first sample of the dataset: ",dataset[0,:,:])


	all_structure = np.asarray(all_structure)
	all_structure_begin = np.asarray(all_structure_begin)
	new_all_structure = []
	
	# If you want to remove some structures
	'''
	count = 0
	new_all_structure.append(all_structure[0])
	for i in xrange(all_structure.shape[0]):
		if(count == 9):
			new_all_structure.append(all_structure[i])
			count = 0
		count = count + 1 
	all_structure = np.asarray(new_all_structure)
	'''
	
	print all_structure.shape
	
	# Calculate kabasch between structures
	kabsch_matrix       = create_kabsch_matrix(all_structure, number_pathway)
	kabsch_matrix_begin = create_kabsch_matrix(all_structure_begin, number_pathway)
	print "kabsch_matrix: ", kabsch_matrix.shape
	print "kabsch_matrix_begin: ", kabsch_matrix.shape

	# Save similarity matrix between final structures
	dataset_name   = filename[0:filename.find("_")]
	dataset_name   = "generate_dataset/heatmap_" + dataset_name + "_end.npy"
	save_heatmap(dataset_name, kabsch_matrix)

	# Save similarity matrix between initial structures
	dataset_name   = filename[0:filename.find("_")]
	dataset_name   = "generate_dataset/heatmap_" + dataset_name + "_begin.npy"
	save_heatmap(dataset_name, kabsch_matrix_begin)

	all_ep = np.asarray(all_ep)
	#kabsch_matrix_mean = np.mean(kabsch_matrix, axis=1)
	''' aux = []
	for i in xrange(kabsch_matrix.shape[0]):
		mean_value = 0.0
		for j in xrange(kabsch_matrix.shape[1]):
			if(i != j):
				mean_value = mean_value + kabsch_matrix[i,j]
		mean_value = mean_value/(kabsch_matrix.shape[1]-1.)
		aux.append(mean_value)
	aux = np.asarray(aux)
	for i in xrange(aux.shape[0]):
		print("%f\t%f" %(aux[i],all_ep[i]))
	#print kabsch_matrix_mean
	#print kabsch_matrix_mean.shape
	#print all_ep

	print 'end'

	# Salvar matriz de similaridade entre as estruturas todas as estruturas de todos os pathways
	dataset_name   = filename[0:filename.find("_")]
	dataset_name   = "generate_dataset/heatmap_" + dataset_name + "_all.dat"
	total_structure_begin = np.asarray(total_structure_begin)
	kabsch_matrix_all = mean_pathway_kabsch(total_structure_begin)
	save_heatmap(dataset_name, kabsch_matrix_all) '''
	

if __name__ == "__main__":
	main()
