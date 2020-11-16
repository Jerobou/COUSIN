#!/usr/bin/env python3

################################################################################
################################################################################
################################ Imports #######################################
################################################################################
################################################################################

import sys,re,os,cgi,csv, getopt, argparse, time, math, shutil
from copy import deepcopy
from random import randrange, uniform
from numpy import *
import scipy as sp
import numpy as np
import pandas as pd

from sklearn.cluster import KMeans
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler

import matplotlib
matplotlib.use('Agg')

import matplotlib.pyplot as plt
from scipy.stats import gaussian_kde
from scipy.stats import mannwhitneyu

import seaborn as sns

import cousin_functions

################################################################################
################################################################################
################################################################################
################## TAKING DATA AND PARAMETERS AS INPUT #########################
################################################################################
################################################################################
################################################################################


################################################################################
############################## PARAMETERS ######################################
################################################################################
parameters_value = '''


################################################################################
################################################################################
################################################################################
############################## COUSIN V 0.4 ####################################
################################################################################
################################################################################
################################################################################

################################################################################
PARAMETERS TAKEN AS INPUT :
'''

################################################################################
######################### USER CHOICES FUNCTION ################################
################################################################################
## Here, parameters as defined above are given a value thanks to user choice.

start_time = time.time()

if __name__ == '__main__' :
	parser = argparse.ArgumentParser(prog="Mandatory arguments", formatter_class=argparse.RawTextHelpFormatter)

	parser.add_argument("-s", nargs= "?", type = str,  metavar = "<path_to_input_seq_file> The path to your input file containing DNA or AA sequences")
	parser.add_argument("-c", nargs = "?", type = str,  metavar = "<path_to_input_cu_table_file> The path to your input file containing the Codon Usage Table")
	parser.add_argument("-g", nargs= "?", default = 1, type = int,  metavar = "<number> The genetic code used (1= standard)")
	parser.add_argument("-X", default = "no", nargs= "?", type = str,  metavar = "<Yes|No> Perform a specific COUSIN calculation with deletion of amino acids with no CUB in the reference")
	parser.add_argument("-S", default = "no", nargs= "?", type = str,  metavar = "Perform routine simulations")
	parser.add_argument("-J", default = "no", nargs= "?", type = str,  metavar = "Create vectors of occurences, frequencies and RSCU value for each query")
	parser.add_argument("-b", default = "nothing", type = str,  metavar = "<path_to_optimal_codons_file> The path to your input file containing a list of optimal codons.")
	parser.add_argument("-o", nargs= "?", default = "output", type = str, metavar = "<path_to_output_directory> The path of you output directory which will contain the results of your analysis.")

	sub_parser = parser.add_subparsers(title ="7 differents steps", description="Perform a basic calculation on queries ; Create a Codon Usage Table from the query dataset ; Optimize sequences thru a random, and random-guided or a \"one amino acid - one codon\" method ; Do a simulation on the query dataset ; Do an analysis by pattern on the query dataset ; Do a clustering step on each variable ; Compare two Codon Usage Tables altogether.")

	calculation_step = sub_parser.add_parser("calculation", help = "Perform basic calculation.")
	calculation_step.set_defaults(which = 'calculation')
	calculation_args = calculation_step.add_argument_group('Mandatory parameters for the basic calculation step')
	calculation_args.add_argument("-s", required = True, type = str,  metavar = "<path_to_input_seq_file> The path to your input file containing DNA or AA sequences")
	calculation_args.add_argument("-c", required = True, type = str,  metavar = "<path_to_input_cu_table_file> The path to your input file containing the Codon Usage Table")
	calculation_args.add_argument("-g", default = 1, type = int,  metavar = "<number> The genetic code used (1= standard)")
	calculation_args.add_argument("-X", default = "no", nargs= "?", type = str,  metavar = "<Yes|No> Perform a specific COUSIN calculation with deletion of amino acids with no CUB in the reference")
	calculation_args.add_argument("-S", default = "no", nargs= "?", type = str,  metavar = "Perform routine simulations")
	calculation_args.add_argument("-J", default = "no", nargs= "?", type = str,  metavar = "Create vectors of occurences, frequencies and RSCU value for each query")
	calculation_args.add_argument("-b", default = "nothing", type = str,  metavar = "<path_to_optimal_codons_file> The path to your input file containing a list of optimal codons.")
	calculation_args.add_argument("-o", default = "output", type = str, metavar = "<path_to_output_directory> The path of you output directory which will contain the results of your analysis.")

	create_table_step = sub_parser.add_parser("create_table", help = "Create a codon Usage Table from the input queries.")
	create_table_step.set_defaults(which = 'create_table')
	create_table_args = create_table_step.add_argument_group('Mandatory parameters for the creation of a Codon Usage Table')
	create_table_args.add_argument("-s", required = True, type = str,  metavar = "<path_to_input_seq_file> The path to your input file containing DNA or AA sequences")
	create_table_args.add_argument("-g", default = 1, type = int,  metavar = "<number> The genetic code used (1= standard)")
	create_table_args.add_argument("-X", default = "no", nargs= "?", type = str,  metavar = "<Yes|No> Perform a specific COUSIN calculation with deletion of amino acids with no CUB in the reference")
	create_table_args.add_argument("-S", default = "no", nargs= "?", type = str,  metavar = "Perform routine simulations")
	create_table_args.add_argument("-J", default = "no", nargs= "?", type = str,  metavar = "Create vectors of occurences, frequencies and RSCU value for each query")
	create_table_args.add_argument("-b", default = "nothing", type = str,  metavar = "<path_to_optimal_codons_file> The path to your input file containing a list of optimal codons.")
	create_table_args.add_argument("-o", default = "output", type = str, metavar = "<path_to_output_directory> The path of you output directory which will contain the results of your analysis.")

	pattern_analysis_step = sub_parser.add_parser("pattern_analysis", help = "Perform an analysis by pattern.")
	pattern_analysis_step.set_defaults(which='pattern_analysis')
	pattern_analysis_args = pattern_analysis_step.add_argument_group('Mandatory parameters for pattern analysis step')
	pattern_analysis_args.add_argument("-s", required = True, type = str,  metavar = "<path_to_input_seq_file> The path to your input file containing DNA or AA sequences")
	pattern_analysis_args.add_argument("-c", required = True, type = str,  metavar = "<path_to_input_cu_table_file> The path to your input file containing the Codon Usage Table")
	pattern_analysis_args.add_argument("-z", default = "nothing", type = str, metavar = "<path_to_pattern_file> The path to your input file containing patterns four your analysis. Must be in a .csv format with a separator (optional)")
	pattern_analysis_args.add_argument("-g", default = 1, type = int,  metavar = "<number> The genetic code used (1= standard)")
	pattern_analysis_args.add_argument("-X", default = "no", nargs= "?", type = str,  metavar = "<Yes|No> Perform a specific COUSIN calculation with deletion of amino acids with no CUB in the reference")
	pattern_analysis_args.add_argument("-S", default = "no", nargs= "?", type = str,  metavar = "Perform routine simulations")
	pattern_analysis_args.add_argument("-J", default = "no", nargs= "?", type = str,  metavar = "Create vectors of occurences, frequencies and RSCU value for each query")
	pattern_analysis_args.add_argument("-b", default = "nothing", type = str,  metavar = "<path_to_optimal_codons_file> The path to your input file containing a list of optimal codons.")
	pattern_analysis_args.add_argument("-o", default = "output", type = str, metavar = "<path_to_output_directory> The path of you output directory which will contain the results of your analysis.")
	pattern_analysis_args.add_argument("-x", default = "No", type = str, metavar = "If you want to display graphical output")

	optimization_step = sub_parser.add_parser("optimization", help = "Perform an optimization analysis")
	optimization_step.set_defaults(which='optimization')
	optimization_args = optimization_step.add_argument_group('Mandatory parameters for optimization step')
	optimization_args.add_argument("-s", required = True, type = str,  metavar = "<path_to_input_seq_file> The path to your input file containing DNA or AA sequences")
	optimization_args.add_argument("-k", default = "DNA", type = str,  metavar = "<DNA/AA> The format of your sequences")
	optimization_args.add_argument("-c", required = True, type = str,  metavar = "<path_to_input_cu_table_file> The path to your input file containing the Codon Usage Table")
	optimization_args.add_argument("-g", default = 1, type = int,  metavar = "<number> The genetic code used (1= standard)")
	optimization_args.add_argument("-X", default = "no", nargs= "?", type = str,  metavar = "<Yes|No> Perform a specific COUSIN calculation with deletion of amino acids with no CUB in the reference")
	optimization_args.add_argument("-S", default = "no", nargs= "?", type = str,  metavar = "Perform routine simulations")
	optimization_args.add_argument("-J", default = "no", nargs= "?", type = str,  metavar = "Create vectors of occurences, frequencies and RSCU value for each query")
	optimization_args.add_argument("-b", default = "nothing", type = str,  metavar = "<path_to_optimal_codons_file> The path to your input file containing a list of optimal codons.")
	optimization_args.add_argument("-o", default = "output", type = str, metavar = "<path_to_output_directory> The path of you output directory which will contain the results of your analysis.")
	optimization_args.add_argument("-x", default = 1, type = int, metavar = "The kind of optimization wanted") # between 1 and 8
	optimization_args.add_argument("-n", type = int, default = 1, metavar = "The number of optimization repetitions ") #between 0 and 1000 maybe ?
	optimization_args.add_argument("-r", type = int, default = 0, metavar = "Reverse the frequencies found in the codon usage table (optional)") #Yes or no

	simulation_analysis_step = sub_parser.add_parser("simulation_analysis", help = "Perform a simulation analysis")
	simulation_analysis_step.set_defaults(which='simulation_analysis')
	simulation_analysis_args = simulation_analysis_step.add_argument_group('Mandatory parameters for simulation step')
	simulation_analysis_args.add_argument("-s", required = True, type = str,  metavar = "<path_to_input_seq_file> The path to your input file containing DNA or AA sequences")
	simulation_analysis_args.add_argument("-c", required = True, type = str,  metavar = "<path_to_input_cu_table_file> The path to your input file containing the Codon Usage Table")
	simulation_analysis_args.add_argument("-g", default = 1, type = int,  metavar = "<number> The genetic code used (1= standard)")
	simulation_analysis_args.add_argument("-X", default = "no", nargs= "?", type = str,  metavar = "<Yes|No> Perform a specific COUSIN calculation with deletion of amino acids with no CUB in the reference")
	simulation_analysis_args.add_argument("-S", default = "no", nargs= "?", type = str,  metavar = "Perform routine simulations")
	simulation_analysis_args.add_argument("-J", default = "no", nargs= "?", type = str,  metavar = "Create vectors of occurences, frequencies and RSCU value for each query")
	simulation_analysis_args.add_argument("-b", default = "nothing", type = str,  metavar = "<path_to_optimal_codons_file> The path to your input file containing a list of optimal codons.")
	simulation_analysis_args.add_argument("-n", type = int, default = 500, metavar = "The number of sequences created during simulation")
	simulation_analysis_args.add_argument("-o", default = "output", type = str, metavar = "<path_to_output_directory> The path of you output directory which will contain the results of your analysis.")

	clustering_analysis_step = sub_parser.add_parser("clustering_analysis", help = "Perform a PCA and a clustering analysis on COUSIN tool variables")
	clustering_analysis_step.set_defaults(which='clustering_analysis')
	clustering_analysis_args = clustering_analysis_step.add_argument_group('Mandatory parameters for clustering_analysis step')
	clustering_analysis_args.add_argument("-s", required = True, type = str,  metavar = "<path_to_input_seq_file> The path to your input file containing DNA or AA sequences")
	clustering_analysis_args.add_argument("-c", required= True, type = str,  metavar = "<path_to_input_cu_table_file> The path to your input file containing the Codon Usage Table")
	clustering_analysis_args.add_argument("-g", default = 1, type = int,  metavar = "<number> The genetic code used (1= standard)")
	clustering_analysis_args.add_argument("-X", default = "no", nargs= "?", type = str,  metavar = "<Yes|No> Perform a specific COUSIN calculation with deletion of amino acids with no CUB in the reference")
	clustering_analysis_args.add_argument("-S", default = "no", nargs= "?", type = str,  metavar = "Perform routine simulations")
	clustering_analysis_args.add_argument("-J", default = "no", nargs= "?", type = str,  metavar = "Create vectors of occurences, frequencies and RSCU value for each query")
	clustering_analysis_args.add_argument("-b", default = "nothing", type = str,  metavar = "<path_to_optimal_codons_file> The path to your input file containing a list of optimal codons.")
	clustering_analysis_args.add_argument("-w", default = "variables", type = str,  metavar = "The data on where clustering should be done (codons or variables)")
	clustering_analysis_args.add_argument("-n", default = 2, type = str,  metavar = "The number of clusters you want to create (not mandatory ==> automatically detected by x-means algorithm) / If you indicate the path to a pattern file, it could be use to create predefined clusters.")
	clustering_analysis_args.add_argument("-o", default = "output", type = str, metavar = "<path_to_output_directory> The path of you output directory which will contain the results of your analysis.")

	compare_data_step = sub_parser.add_parser("compare_data", help = "Compares two datasets thru calculation of euclidian distances")
	compare_data_step.set_defaults(which='compare_data')
	compare_data_args = compare_data_step.add_argument_group('Mandatory parameters for the codon usage tables comparison step')
	compare_data_args.add_argument("-C", required= True, type = str,  metavar = "<path_to_input_cu_table_file> The path to your input file containing the first Codon Usage Table")
	compare_data_args.add_argument("-c", required= True, type = str,  metavar = "<path_to_input_cu_table_file> The path to your input file containing the first Codon Usage Table")
	compare_data_args.add_argument("-d", required= True, type = str,  metavar = "<path_to_input_cu_table_file> The path to your input file containing the second Codon Usage Table")
	compare_data_args.add_argument("-g", default = 1, type = int,  metavar = "<number> The genetic code used (1= standard)")
	compare_data_args.add_argument("-o", default = "output", type = str, metavar = "<path_to_output_directory> The path of you output directory which will contain the results of your analysis.")

	args = vars(parser.parse_args())

if "s" in args :
	if args["s"] != None :
		with open(args["s"], 'r') as my_file_seq :
			input_file_seq = my_file_seq.read()

if "k" in args :
	input_file_seq_type = args["k"]
else :
	input_file_seq_type = "DNA"

genet_code_type = args["g"]

perform_COUSIN_with_deletion = args["X"]

perform_vector_analysis = args["J"]


if args["b"] != "nothing" :
	with open(args["b"], 'r') as optimal_codons_file_path :
		optimal_codons_file = optimal_codons_file_path.read()
else :
	optimal_codons_file = "nothing"

output_dir = args["o"]

routine_simulation = args["S"]

cousin_functions.create_output_dir(output_dir, args["which"], routine_simulation, perform_vector_analysis)

genetic_code_ref = {}
optimal_codons_list = {}

def compare_data(genet_code_type, output_dir) :

	compare_type = args["C"]

	print ('### You have chosen the "compare_data" step. ###')

	with open(args["c"], 'r') as my_file_dataset :
		input_file_dataset_1 = my_file_dataset.read()

	with open(args["d"], 'r') as my_file_dataset :
		input_file_dataset_2 = my_file_dataset.read()

	if compare_type == "cvc" :

		print ('## Comparing two codon usage tables. ##')

		cousin_functions.sequential_order__data_comparison = []
		cousin_functions.results_for_file_data_comparison = {}

		genetic_code_1 = cousin_functions.get_genetic_code(genet_code_type, output_dir)
		genetic_code_2 = cousin_functions.get_genetic_code(genet_code_type, output_dir)

		cousin_functions.get_CU_table_content(input_file_dataset_1, genetic_code_1)
		cousin_functions.get_CU_table_content(input_file_dataset_2, genetic_code_2)
		cousin_functions.compare_data(genetic_code_1, genetic_code_2)

		results_data_comparison = cousin_functions.send_results_data_comparison()
		results_data_comparison_list = cousin_functions.send_results_data_comparison_list()
		if not os.path.isdir(output_dir + "/data_comparison") :
			os.system('mkdir ./' + output_dir + "/data_comparison")
		cousin_functions.write_file(output_dir + "/data_comparison/comparison_results.txt", results_data_comparison, results_data_comparison_list)

	elif compare_type == "svs" :

		print ('## Comparing two datasets of sequences. ##')


		genetic_code_ref = cousin_functions.get_genetic_code(genet_code_type, output_dir)

		liste_sequences_1 = {}
		liste_sequences_2 = {}
		cousin_functions.content_seq = {}
		cousin_functions.parse_input_file_seq(input_file_dataset_1, input_file_seq_type, output_dir)
		for header_seq in cousin_functions.content_seq :
			liste_sequences_1[header_seq] = cousin_functions.seq_to_dict(cousin_functions.content_seq[header_seq], header_seq, genetic_code_ref, output_dir)

		cousin_functions.content_seq = {}
		cousin_functions.parse_input_file_seq(input_file_dataset_2, input_file_seq_type, output_dir)

		for header_seq in cousin_functions.content_seq :
			liste_sequences_2[header_seq] = cousin_functions.seq_to_dict(cousin_functions.content_seq[header_seq], header_seq, genetic_code_ref, output_dir)

		for header_seq in liste_sequences_1 :
			for header_seq_2 in liste_sequences_2 :

				cousin_functions.sequential_order__data_comparison = []
				cousin_functions.results_for_file_data_comparison = {}

				cousin_functions.compare_data(liste_sequences_1[header_seq], liste_sequences_2[header_seq_2])

				results_data_comparison = cousin_functions.send_results_data_comparison()
				results_data_comparison_list = cousin_functions.send_results_data_comparison_list()
				if not os.path.isdir(output_dir + "/data_comparison") :
					os.system('mkdir ./' + output_dir + "/data_comparison")
				cousin_functions.write_file(output_dir + "/data_comparison/comparison_results_" + header_seq[:50] + "_vs_" + header_seq_2[:50] + ".txt", results_data_comparison, results_data_comparison_list)



	elif compare_type == "cvs" :

		print ('## Comparing a codon usage table to a dataset of sequences. ##')


		genetic_code_ref = cousin_functions.get_genetic_code(genet_code_type, output_dir)
		genetic_code_1 = cousin_functions.get_genetic_code(genet_code_type, output_dir)
		liste_sequences_2 = {}
		cousin_functions.get_CU_table_content(input_file_dataset_1, genetic_code_1)

		cousin_functions.parse_input_file_seq(input_file_dataset_2, input_file_seq_type, output_dir)
		for header_seq in cousin_functions.content_seq :
			liste_sequences_2[header_seq] = cousin_functions.seq_to_dict(cousin_functions.content_seq[header_seq], header_seq, genetic_code_ref, output_dir)

			cousin_functions.sequential_order__data_comparison = []
			cousin_functions.results_for_file_data_comparison = {}

			cousin_functions.compare_data(genetic_code_1, liste_sequences_2[header_seq])

			results_data_comparison = cousin_functions.send_results_data_comparison()
			results_data_comparison_list = cousin_functions.send_results_data_comparison_list()
			if not os.path.isdir(output_dir + "/data_comparison") :
				os.system('mkdir ./' + output_dir + "/data_comparison")
			cousin_functions.write_file(output_dir + "/data_comparison/comparison_results_cu_table_vs_" + header_seq[:50] + ".txt", results_data_comparison, results_data_comparison_list)


	elif compare_type == "svc" :

		print ('## Comparing a dataset of sequences to a codon usage table. ##')


		genetic_code_ref = cousin_functions.get_genetic_code(genet_code_type, output_dir)
		liste_sequences_1 = {}
		genetic_code_2 = cousin_functions.get_genetic_code(genet_code_type, output_dir)
		cousin_functions.get_CU_table_content(input_file_dataset_2, genetic_code_2)

		cousin_functions.parse_input_file_seq(input_file_dataset_1, input_file_seq_type, output_dir)
		for header_seq in cousin_functions.content_seq :
			liste_sequences_1[header_seq] = cousin_functions.seq_to_dict(cousin_functions.content_seq[header_seq], header_seq, genetic_code_ref, output_dir)

			cousin_functions.sequential_order__data_comparison = []
			cousin_functions.results_for_file_data_comparison = {}

			cousin_functions.compare_data(liste_sequences_1[header_seq], genetic_code_2)

			results_data_comparison = cousin_functions.send_results_data_comparison()
			results_data_comparison_list = cousin_functions.send_results_data_comparison_list()
			if not os.path.isdir(output_dir + "/data_comparison") :
				os.system('mkdir ./' + output_dir + "/data_comparison")
			cousin_functions.write_file(output_dir + "/data_comparison/comparison_results_cu_table_vs_" + header_seq[:50] + ".txt", results_data_comparison, results_data_comparison_list)

	print ('### The "Compare data" step is done ! ###')
	if os.stat("./" + output_dir + "/log/log_cousin.txt").st_size != 0 :
		print ("### IMPORTANT ### Some errors / warning arised during you analysis. Please check the log file to see what went wrong.") 

def prepare_data(input_file_seq, genet_code_type, optimal_codons_file, output_dir) :

	global genetic_code_ref
	global optimal_codons_list

	if args["which"] == "create_table" :

		print ('#### You have chosen the "create_table" step. ####')

		print ("## creating a Codon Usage Table. ##")

		cousin_functions.parse_input_file_seq(input_file_seq, input_file_seq_type, output_dir)
		genetic_code_ref = cousin_functions.get_genetic_code(genet_code_type, output_dir)
		cousin_functions.create_cu_table(genetic_code_ref, output_dir)
		optimal_codons_list = cousin_functions.calc_optimal_codons(genetic_code_ref, optimal_codons_file)

		print ('## Your Codon Usage Table has been created in a kazusa-style format. ##')
		if os.stat("./" + output_dir + "/log/log_cousin.txt").st_size != 0 :
			print ("### IMPORTANT ### Some errors / warnings arised during you analysis. Please check the log file to see what went wrong.") 

	else :

		with open(args["c"], 'r') as my_file_CU :
			input_file_CU_table = my_file_CU.read()
		cousin_functions.parse_input_file_seq(input_file_seq, input_file_seq_type, output_dir)
		genetic_code_ref = cousin_functions.get_genetic_code(genet_code_type, output_dir)
		cousin_functions.get_CU_table_content(input_file_CU_table, genetic_code_ref)
		optimal_codons_list = cousin_functions.calc_optimal_codons(genetic_code_ref, optimal_codons_file)

	if input_file_seq_type == "DNA" :

		print ('### Doing the mandatory "calculation step" of COUSIN. ###')
		cousin_functions.unvoid_empty_frequencies(genetic_code_ref)
		cousin_functions.get_GC_and_AT_max_codons(genetic_code_ref)
		cousin_functions.calculate_input_values(input_file_seq_type, genetic_code_ref, optimal_codons_list, perform_COUSIN_with_deletion, perform_vector_analysis, output_dir)
		if routine_simulation == "yes" or routine_simulation == "Yes" or routine_simulation == "YES" or routine_simulation == "Y" or routine_simulation == "y" :
			print ('## Doing the basic simulation analysis. ##')
			cousin_functions.seq_simulations(genetic_code_ref,500,100,200,150, optimal_codons_list, perform_COUSIN_with_deletion, output_dir)
			cousin_functions.treat_simulation_data(output_dir)
		results = cousin_functions.send_results()
		results_list = cousin_functions.send_results_list()
		print ('### Calculation step done !###')

		if not os.path.isdir(output_dir + "/results") :
			os.system('mkdir ./' + output_dir + "/results")
		cousin_functions.write_file(output_dir + "/results/results.txt", results, results_list)
	elif input_file_seq_type == "AA" :
		print ('## Preparing AA data for optimization step (if you chosed AA as a -k argument for another step, an error will occur). ##')
		cousin_functions.unvoid_empty_frequencies(genetic_code_ref)
		cousin_functions.get_GC_and_AT_max_codons(genetic_code_ref)
		print ('## Preparation of AA data done ! ##')

	else :
		print("### ERROR ### Please enter a correct sequence format (option -k ; AA for amino acids sequence(s) and DNA for nucleotides sequence(s)).")
		print ("No COUSIN analysis has been done. Exiting...")

	if args["which"] == "calculation" : 
		if os.stat("./" + output_dir + "/log/log_cousin.txt").st_size != 0 :
			print ("### IMPORTANT ### Some errors / warning arised during you analysis. Please check the log file to see what went wrong.") 

#cousin_functions.do_wilcox_test(set_of_seq["cousin_value"], set_of_seq["cousin_value_muta_bias"])

def do_optimization(output_dir) :

	print ('### You have chosen the "optimization" additional step. ####')

	input_file_seq_type = args["k"]
	optimization_type = args["x"]
	optimization_number = args["n"]
	reverse_frequencies = args["r"]

	print ('## Doing optimization... ##')

	cousin_functions.optimize(genetic_code_ref, optimal_codons_list, perform_COUSIN_with_deletion, optimization_type, optimization_number, input_file_seq_type, reverse_frequencies, output_dir)
	results_optimization = cousin_functions.send_results_optimization()
	results_optimization_list = cousin_functions.send_results_optimization_list()
	if not os.path.isdir(output_dir + "/optimization") :
		os.system('mkdir ./' + output_dir + "/optimization")
	cousin_functions.write_file_optimization(output_dir + "/optimization/optimization_results.txt", results_optimization, results_optimization_list)
	
	print ("## Optimization done ! ##")
	if os.stat("./" + output_dir + "/log/log_cousin.txt").st_size != 0 :
		print ("### IMPORTANT ### Some errors / warning arised during you analysis. Please check the log file to see what went wrong.") 

def do_pattern_analysis(output_dir) :

	print ('### You have chosen the "pattern_analysis" additional step. ####')

	print ('## Doing pattern analysis using the pattern file given with the "-z" parameter... ##')


	file_pattern = args["z"]

	cousin_functions.pattern_analysis_step(file_pattern)
	cousin_functions.treat_pattern_analysis_step_data(output_dir)

	print ("## Pattern analysis done ! ##")
	if os.stat("./" + output_dir + "/log/log_cousin.txt").st_size != 0 :
		print ("### IMPORTANT ### Some errors / warning arised during you analysis. Please check the log file to see what went wrong.") 


def do_simulation_analysis(output_dir) :

	print ('### You have chosen the "simulation_analysis" additional step. ####')

	print ('## Doing simulation analysis... ##')

	simulation_number = args["n"]
	if not os.path.isdir(output_dir + "/simulation_results") :
		os.system('mkdir ./' + output_dir + "/simulation_results")
	if not os.path.isdir(output_dir + "/simulation_graphs") :
		os.system('mkdir ./' + output_dir + "/simulation_graphs")    
	cousin_functions.input_seq_simulation_analysis(genetic_code_ref, simulation_number, optimal_codons_list, perform_COUSIN_with_deletion, output_dir)
	cousin_functions.treat_input_simulation_analysis_data(output_dir)
	results_simulation = cousin_functions.send_results_simulation()
	results_simulation_list = cousin_functions.send_results_simulation_list()
	cousin_functions.write_file_simulation(output_dir + "/simulation_results/simulation_results.txt", results_simulation, results_simulation_list)

	print ('## Simulation analysis done ! ##')
	if os.stat("./" + output_dir + "/log/log_cousin.txt").st_size != 0 :
		print ("### IMPORTANT ### Some errors / warning arised during you analysis. Please check the log file to see what went wrong.") 


def do_clustering_analysis(output_dir) :

	print ('### You have chosen the "clustering_analysis" additional step. ####')


	print ('## Doing clustering analysis... ##')

	try:

		number_of_clusters = int(args["n"])

	except ValueError:

		number_of_clusters = args["n"]

	type_of_analysis = args["w"]

	results = cousin_functions.send_results()

	cousin_functions.perform_clustering_analysis(results, genetic_code_ref, number_of_clusters, type_of_analysis, output_dir)

	results_clustering = cousin_functions.send_results_clustering()
	results_clustering_list = cousin_functions.send_results_clustering_list()

	if not os.path.isdir(output_dir + "/clustering") :
		os.system('mkdir ./' + output_dir + "/clustering")
	cousin_functions.write_file_clustering(output_dir + "/clustering/clustering_results.txt", results_clustering, results_clustering_list)

	print ('## Clustering analysis done ! ##')
	if os.stat("./" + output_dir + "/log/log_cousin.txt").st_size != 0 :
		print ("### IMPORTANT ### Some errors / warning arised during you analysis. Please check the log file to see what went wrong.") 




print ("""

#################################################################################
#               Welcome on COUSIN (COdon Usage Similarity INdex)!               #
#                                     V 1.0                                     #
#         Created by Jerome Bourret, Samuel Alizon and Ignacio G. Bravo         #
#               Any question ? Contact me : jerome.bourret@ird.fr               #
#                                    Enjoy !                                    #
#################################################################################

	""")

if args["which"] == "compare_data" :
	compare_data(genet_code_type, output_dir)
else :
	prepare_data(input_file_seq, genet_code_type, optimal_codons_file, output_dir)

	if args["which"] == "optimization" :
		do_optimization(output_dir)
	elif args["which"] == "simulation_analysis" :
		do_simulation_analysis(output_dir)
	elif args["which"] == "pattern_analysis" :
		do_pattern_analysis(output_dir)
	elif args["which"] == "clustering_analysis" :
		do_clustering_analysis(output_dir)
	elif args["which"] != "calculation" and args["which"] != "create_table" : 
		print ("### ERROR ### You have chosen a wrong option (typo ?). Only the basic calculation step has been done.")
		print ("Exiting COUSIN without doing any option...")

end_time = time.time()

elapsed_time = round((end_time - start_time), 3)
elapsed_time_hours = round(int(elapsed_time // 3600), 3)
elapsed_time = round((elapsed_time % 3600), 3)
elapsed_time_minutes = round(int(elapsed_time // 60), 3)
elapsed_time_seconds = round((elapsed_time % 60), 3)

print ("This COUSIN analysis took " + str(elapsed_time_hours) + " HH " + str(elapsed_time_minutes) + " MM " + str(elapsed_time_seconds) + " SS.")

