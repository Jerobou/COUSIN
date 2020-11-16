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
from scipy.cluster import hierarchy
from scipy.cluster.hierarchy import dendrogram, linkage

import seaborn as sns

import pyclustering

from pyclustering.samples.definitions import SIMPLE_SAMPLES, FCPS_SAMPLES;

from pyclustering.cluster.center_initializer import kmeans_plusplus_initializer;

from pyclustering.cluster import cluster_visualizer;
from pyclustering.cluster.xmeans import xmeans, splitting_type;
from pyclustering.cluster.kmeans import kmeans;

from pyclustering.utils import read_sample;

################################################################################
################################################################################
################################################################################
############# READING CU TABLE AND SEQUENCES (user input) ######################
################################################################################
################################################################################
################################################################################

################################################################################
########################## FASTA FILES PARSER ##################################
################################################################################

def replace_forbidden_char(string, forbidden_char, allowed_char):
	for element in forbidden_char :
		if element in string :
			string = string.replace(element, allowed_char)
	return  string

def parse_input_file_seq (input_file_seq, input_file_seq_type, output_dir) :
	lines = input_file_seq.split("\n")
	line_seq = ''
	header_seq = 'sequence'
	cpt = 0
	forbidden_chars = set('\/\\\?\%\*\:\|\"\<\>\.')
	alphabet_DNA = set('AaTtUuGgCc')
	alphabet_AA = set('AaRrDdNnCcEeQqGgHhIiLlKkMmFfPpSsTtWwYyVvXx\*')
	

	if not os.path.isdir("./" + output_dir + "/log"):
		os.makedirs("./" + output_dir + "/log")

	log_file = open("./" + output_dir + "/log/log_cousin.txt", "a")

	for line in lines :
		
		if line != '' :

			if line[0] == '>' and line_seq == '' :
				header_seq = line[1:].replace('\t', ' ') #POUR EVITER D'AVOIR DES "\t" qui gènent le fichier de résultats final. Autre solution ? Car peut être gourmand en énergie
				if any((c in forbidden_chars) for c in header_seq): 
					log_file.write("## WARNING : the header \"" + header_seq + "\" contains forbidden characters (such as \"\/\", \"\\\", \"\?\", \"\%\", \"\*\", \"\:\", \"\|\", \"\"\", \"\<\", \"\>\", \"\.\"). These characters have been replaced by an underscore. ## \n\n")
					header_seq = replace_forbidden_char(header_seq, ['/', '\\', '?', '%', '*', ':', '|','"', '<', '>', '.'], "_")
			elif line[0] == '>' and line_seq != '' :
				if header_seq not in content_seq :


					if (len(line_seq) ) == 0 :
						log_file.write("## SEQUENCE ERROR : the sequence related to the header " + header_seq + "is empty. (check your file) ##\n")
						header_seq = line[1:].replace('\t', ' ') #POUR EVITER D'AVOIR DES "\t" qui gènent le fichier de résultats final. Autre solution ? Car peut être gourmand en énergie
						
						if any((c in forbidden_chars) for c in header_seq): 
							log_file.write("## WARNING : the header \"" + header_seq + "\" contains forbidden characters (such as \"\/\", \"\\\", \"\?\", \"\%\", \"\*\", \"\:\", \"\|\", \"\"\", \"\<\", \"\>\", \"\.\"). These characters have been replaced by an underscore. ## \n\n")
							header_seq = replace_forbidden_char(header_seq, ['/', '\\', '?', '%', '*', ':', '|','"', '<', '>', '.'], "_")
						line_seq = ''
					

					if input_file_seq_type == "DNA" and (len(line_seq) % 3) != 0 :
							
						log_file.write("## SEQUENCE ERROR : the sequence related to the header " + header_seq + " doesn't have a number that can be divided by 3. (check your file) ##\n")
						header_seq = line[1:].replace('\t', ' ') #POUR EVITER D'AVOIR DES "\t" qui gènent le fichier de résultats final. Autre solution ? Car peut être gourmand en énergie
						
						if any((c in forbidden_chars) for c in header_seq):  
							log_file.write("## WARNING : the header \"" + header_seq + "\" contains forbidden characters (such as \"\/\", \"\\\", \"\?\", \"\%\", \"\*\", \"\:\", \"\|\", \"\"\", \"\<\", \"\>\", \"\.\"). These characters have been replaced by an underscore. ## \n\n")
							header_seq = replace_forbidden_char(header_seq, ['/', '\\', '?', '%', '*', ':', '|','"', '<', '>', '.'], "_")
						line_seq = ''
					

					else :
						cpt = 0
						cpt_check_char = 0
						if input_file_seq_type == "DNA" : 
							for c in line_seq:
								if c not in alphabet_DNA :
									cpt_check_char =+ 1
							if cpt_check_char > 0 : 
								log_file.write("## SEQUENCE ERROR : the sequence with the header " + header_seq  + " has characters that don't belong in the DNA alphabet ##\n")
							else : 
								content_seq[header_seq] = line_seq.upper()
						elif input_file_seq_type == "AA" : 

							for c in line_seq:
								if c not in alphabet_AA :
									cpt_check_char =+ 1
							if cpt_check_char > 0 : 
								log_file.write("## SEQUENCE ERROR : the sequence with the header " + header_seq  + " has characters that don't belong in the AA alphabet ##\n")
							else : 
								line_seq = line_seq.replace("X", "*") #replace X character, which is the stop codon by *
								content_seq[header_seq] = line_seq.upper()			
								

						header_seq = line[1:].replace('\t', ' ') #POUR EVITER D'AVOIR DES "\t" qui gènent le fichier de résultats final. Autre solution ? Car peut être gourmand en énergie
						if any((c in forbidden_chars) for c in header_seq): 
							log_file.write("## WARNING : the header \"" + header_seq + "\" contains forbidden characters (such as \"\/\", \"\\\", \"\?\", \"\%\", \"\*\", \"\:\", \"\|\", \"\"\", \"\<\", \"\>\", \"\.\"). These characters have been replaced by an underscore. ## \n\n")
							header_seq = replace_forbidden_char(header_seq, ['/', '\\', '?', '%', '*', ':', '|','"', '<', '>', '.'], "_")
						line_seq = ''
				else :
					cpt += 1
					log_file.write("## SEQUENCE ERROR : the header " + header_seq + " is already taken. Sequence ignored (check your file) ##\n")
					header_seq = line[1:].replace('\t', ' ') #POUR EVITER D'AVOIR DES "\t" qui gènent le fichier de résultats final. Autre solution ? Car peut être gourmand en énergie
					if any((c in forbidden_chars) for c in header_seq): 
						log_file.write("## WARNING : the header \"" + header_seq + "\" contains forbidden characters (such as \"\/\", \"\\\", \"\?\", \"\%\", \"\*\", \"\:\", \"\|\", \"\"\", \"\<\", \"\>\", \"\.\"). These characters have been replaced by an underscore. ## \n\n")
						header_seq = replace_forbidden_char(header_seq, ['\/', '\\', '\?', '\%', '\*', '\:', '\|','\"', '\<', '\>', '\.'], "_")
					line_seq = ''
			


			else :
				line_seq = line_seq + line
	

	if header_seq not in content_seq :

		if (len(line_seq) ) == 0 :
			log_file.write("## SEQUENCE ERROR : the sequence related to the header " + header_seq + "is empty. (check your file) ##\n")
		
		if input_file_seq_type == "DNA" and (len(line_seq) % 3) != 0 : 
			log_file.write("## SEQUENCE ERROR : the sequence related to the header " + header_seq + "doesn't have a number that can be divided by 3. (check your file) ##\n")
		
		else :
			cpt_check_char = 0
			if input_file_seq_type == "DNA" : 
				for c in line_seq:
					if c not in alphabet_DNA :
						cpt_check_char =+ 1
				if cpt_check_char > 0 : 
					log_file.write("## SEQUENCE ERROR : the sequence with the header " + header_seq + " has characters that don't belong in the DNA alphabet ##\n")
				else : 
					content_seq[header_seq] = line_seq.upper()
			elif input_file_seq_type == "AA" : 

				for c in line_seq:
					if c not in alphabet_AA :
						cpt_check_char =+ 1
				if cpt_check_char > 0 : 
					log_file.write("## SEQUENCE ERROR : the sequence with the header " + header_seq + " has characters that don't belong in the AA alphabet ##\n")
				else : 
					line_seq = line_seq.replace("X", "*") #replace X character, which is the stop codon by *
					content_seq[header_seq] = line_seq.upper()						
	else :
		log_file.write("## SEQUENCE ERROR : the header " + header_seq + " is already taken. Sequence ignored (check your file) ##\n")

	log_file.close()

################################################################################
################# DEFINITION OF GENETIC CODE USED FOR ANALYSIS #################
################################################################################

def get_genetic_code(number, output_dir) :

	if not os.path.isdir("./" + output_dir + "/log"):
		os.makedirs("./" + output_dir + "/log")

	log_file = open("./" + output_dir + "/log/log_cousin.txt", "a")

	if number == 1 : #The Standard Code (transl_table=1)
		genetic_code = {'F': {'UUU': [], 'UUC': []}, 'S': {'UCU': [], 'UCC': [], 'UCA': [], 'UCG': [], 'AGU': [], 'AGC': []}, 'Y': {'UAU': [], 'UAC': []}, 'C': {'UGU': [], 'UGC': []}, 'L': {'UUA': [], 'UUG': [], 'CUU': [], 'CUC': [], 'CUA': [], 'CUG': []}, '*': {'UAA': [], 'UGA': [], 'UAG': []}, 'W': {'UGG': []}, 'P': {'CCU': [], 'CCC': [], 'CCA': [], 'CCG': []}, 'H': {'CAU': [], 'CAC': []}, 'R': {'CGU': [], 'CGC': [], 'CGA': [], 'CGG': [], 'AGA': [], 'AGG': []}, 'Q': {'CAA': [], 'CAG': []}, 'I': {'AUU': [], 'AUC': [], 'AUA': []}, 'T': {'ACU': [], 'ACC': [], 'ACA': [], 'ACG': []}, 'N': {'AAU': [], 'AAC': []}, 'K': {'AAA': [], 'AAG': []}, 'M': {'AUG': []}, 'V': {'GUU': [], 'GUC': [], 'GUA': [], 'GUG': []}, 'A': {'GCU': [], 'GCC': [], 'GCA': [], 'GCG': []}, 'D': {'GAU': [], 'GAC': []}, 'G': {'GGU': [], 'GGC': [], 'GGA': [], 'GGG': []}, 'E': {'GAA': [], 'GAG': []}}

	elif number == 2 : # The Vertebrate Mitochondrial Code (transl_table=2)
		genetic_code = {'F': {'UUU': [], 'UUC': []}, 'S': {'UCU': [], 'UCC': [], 'UCA': [], 'UCG': [], 'AGU': [], 'AGC': []}, 'Y': {'UAU': [], 'UAC': []}, 'C': {'UGU': [], 'UGC': []}, 'L': {'UUA': [], 'UUG': [], 'CUU': [], 'CUC': [], 'CUA': [], 'CUG': []}, '*': {'UAA': [], 'UAG': [], 'AGA': [], 'AGG': []}, 'W': {'UGA': [], 'UGG': []}, 'P': {'CCU': [], 'CCC': [], 'CCA': [], 'CCG': []}, 'H': {'CAU': [], 'CAC': []}, 'R': {'CGU': [], 'CGC': [], 'CGA': [], 'CGG': []}, 'Q': {'CAA': [], 'CAG': []}, 'I': {'AUU': [], 'AUC': []}, 'T': {'ACU': [], 'ACC': [], 'ACA': [], 'ACG': []}, 'N': {'AAU': [], 'AAC': []}, 'M': {'AUA': [], 'AUG': []}, 'K': {'AAA': [], 'AAG': []}, 'V': {'GUU': [], 'GUC': [], 'GUA': [], 'GUG': []}, 'A': {'GCU': [], 'GCC': [], 'GCA': [], 'GCG': []}, 'D': {'GAU': [], 'GAC': []}, 'G': {'GGU': [], 'GGC': [], 'GGA': [], 'GGG': []}, 'E': {'GAA': [], 'GAG': []}}

	elif number == 3 : # 3. The Yeast Mitochondrial Code (transl_table=3)
		genetic_code = {'F': {'UUU': [], 'UUC': []}, 'S': {'UCU': [], 'UCC': [], 'UCA': [], 'UCG': [], 'AGU': [], 'AGC': []}, 'Y': {'UAU': [], 'UAC': []}, 'C': {'UGU': [], 'UGC': []}, 'L': {'UUA': [], 'UUG': []}, '*': {'UAA': [], 'UAG': []}, 'W': {'UGA': [], 'UGG': []}, 'T': {'CUU': [], 'CUC': [], 'CUA': [], 'CUG': [], 'ACU': [], 'ACC': [], 'ACA': [], 'ACG': []}, 'P': {'CCU': [], 'CCC': [], 'CCA': [], 'CCG': []}, 'H': {'CAU': [], 'CAC': []}, 'R': {'CGU': [], 'CGC': [], 'CGA': [], 'CGG': [], 'AGA': [], 'AGG': []}, 'Q': {'CAA': [], 'CAG': []}, 'I': {'AUU': [], 'AUC': []}, 'N': {'AAU': [], 'AAC': []}, 'M': {'AUA': [], 'AUG': []}, 'K': {'AAA': [], 'AAG': []}, 'V': {'GUU': [], 'GUC': [], 'GUA': [], 'GUG': []}, 'A': {'GCU': [], 'GCC': [], 'GCA': [], 'GCG': []}, 'D': {'GAU': [], 'GAC': []}, 'G': {'GGU': [], 'GGC': [], 'GGA': [], 'GGG': []}, 'E': {'GAA': [], 'GAG': []}}

	elif number == 4 : # The Mold, Protozoan, and Coelenterate Mitochondrial Code and the Mycoplasma/Spiroplasma Code (transl_table=4)
		genetic_code = {'F': {'UUU': [], 'UUC': []}, 'S': {'UCU': [], 'UCC': [], 'UCA': [], 'UCG': [], 'AGU': [], 'AGC': []}, 'Y': {'UAU': [], 'UAC': []}, 'C': {'UGU': [], 'UGC': []}, 'L': {'UUA': [], 'UUG': [], 'CUU': [], 'CUC': [], 'CUA': [], 'CUG': []}, '*': {'UAA': [], 'UAG': []}, 'W': {'UGA': [], 'UGG': []}, 'P': {'CCU': [], 'CCC': [], 'CCA': [], 'CCG': []}, 'H': {'CAU': [], 'CAC': []}, 'R': {'CGU': [], 'CGC': [], 'CGA': [], 'CGG': [], 'AGA': [], 'AGG': []}, 'Q': {'CAA': [], 'CAG': []}, 'I': {'AUU': [], 'AUC': [], 'AUA': []}, 'T': {'ACU': [], 'ACC': [], 'ACA': [], 'ACG': []}, 'N': {'AAU': [], 'AAC': []}, 'K': {'AAA': [], 'AAG': []}, 'M': {'AUG': []}, 'V': {'GUU': [], 'GUC': [], 'GUA': [], 'GUG': []}, 'A': {'GCU': [], 'GCC': [], 'GCA': [], 'GCG': []}, 'D': {'GAU': [], 'GAC': []}, 'G': {'GGU': [], 'GGC': [], 'GGA': [], 'GGG': []}, 'E': {'GAA': [], 'GAG': []}}

	elif number == 5 : # The Invertebrate Mitochondrial Code (transl_table=5)
		genetic_code = {'F': {'UUU': [], 'UUC': []}, 'S': {'UCU': [], 'UCC': [], 'UCA': [], 'UCG': [], 'AGU': [], 'AGC': [], 'AGA': [], 'AGG': []}, 'Y': {'UAU': [], 'UAC': []}, 'C': {'UGU': [], 'UGC': []}, 'L': {'UUA': [], 'UUG': [], 'CUU': [], 'CUC': [], 'CUA': [], 'CUG': []}, '*': {'UAA': [], 'UAG': []}, 'W': {'UGA': [], 'UGG': []}, 'P': {'CCU': [], 'CCC': [], 'CCA': [], 'CCG': []}, 'H': {'CAU': [], 'CAC': []}, 'R': {'CGU': [], 'CGC': [], 'CGA': [], 'CGG': []}, 'Q': {'CAA': [], 'CAG': []}, 'I': {'AUU': [], 'AUC': []}, 'T': {'ACU': [], 'ACC': [], 'ACA': [], 'ACG': []}, 'N': {'AAU': [], 'AAC': []}, 'M': {'AUA': [], 'AUG': []}, 'K': {'AAA': [], 'AAG': []}, 'V': {'GUU': [], 'GUC': [], 'GUA': [], 'GUG': []}, 'A': {'GCU': [], 'GCC': [], 'GCA': [], 'GCG': []}, 'D': {'GAU': [], 'GAC': []}, 'G': {'GGU': [], 'GGC': [], 'GGA': [], 'GGG': []}, 'E': {'GAA': [], 'GAG': []}}

	elif number == 6 : # The Ciliate, Dasycladacean and Hexamita Nuclear Code (transl_table=6)
		genetic_code = {'F': {'UUU': [], 'UUC': []}, 'S': {'UCU': [], 'UCC': [], 'UCA': [], 'UCG': [], 'AGU': [], 'AGC': []}, 'Y': {'UAU': [], 'UAC': []}, 'C': {'UGU': [], 'UGC': []}, 'L': {'UUA': [], 'UUG': [], 'CUU': [], 'CUC': [], 'CUA': [], 'CUG': []}, 'Q': {'UAA': [], 'UAG': [], 'CAA': [], 'CAG': []}, '*': {'UGA': []}, 'W': {'UGG': []}, 'P': {'CCU': [], 'CCC': [], 'CCA': [], 'CCG': []}, 'H': {'CAU': [], 'CAC': []}, 'R': {'CGU': [], 'CGC': [], 'CGA': [], 'CGG': [], 'AGA': [], 'AGG': []}, 'I': {'AUU': [], 'AUC': [], 'AUA': []}, 'T': {'ACU': [], 'ACC': [], 'ACA': [], 'ACG': []}, 'N': {'AAU': [], 'AAC': []}, 'K': {'AAA': [], 'AAG': []}, 'M': {'AUG': []}, 'V': {'GUU': [], 'GUC': [], 'GUA': [], 'GUG': []}, 'A': {'GCU': [], 'GCC': [], 'GCA': [], 'GCG': []}, 'D': {'GAU': [], 'GAC': []}, 'G': {'GGU': [], 'GGC': [], 'GGA': [], 'GGG': []}, 'E': {'GAA': [], 'GAG': []}}

	elif number == 9 : # The Echinoderm Mitochondrial Code (transl_table=9)
		genetic_code = {'F': {'UUU': [], 'UUC': []}, 'S': {'UCU': [], 'UCC': [], 'UCA': [], 'UCG': [], 'AGU': [], 'AGC': [], 'AGA': [], 'AGG': []}, 'Y': {'UAU': [], 'UAC': []}, 'C': {'UGU': [], 'UGC': []}, 'L': {'UUA': [], 'UUG': [], 'CUU': [], 'CUC': [], 'CUA': [], 'CUG': []}, '*': {'UAA': [], 'UAG': []}, 'W': {'UGA': [], 'UGG': []}, 'P': {'CCU': [], 'CCC': [], 'CCA': [], 'CCG': []}, 'H': {'CAU': [], 'CAC': []}, 'R': {'CGU': [], 'CGC': [], 'CGA': [], 'CGG': []}, 'Q': {'CAA': [], 'CAG': []}, 'I': {'AUU': [], 'AUC': [], 'AUA': []}, 'T': {'ACU': [], 'ACC': [], 'ACA': [], 'ACG': []}, 'N': {'AAU': [], 'AAC': [], 'AAA': []}, 'M': {'AUG': []}, 'K': {'AAG': []}, 'V': {'GUU': [], 'GUC': [], 'GUA': [], 'GUG': []}, 'A': {'GCU': [], 'GCC': [], 'GCA': [], 'GCG': []}, 'D': {'GAU': [], 'GAC': []}, 'G': {'GGU': [], 'GGC': [], 'GGA': [], 'GGG': []}, 'E': {'GAA': [], 'GAG': []}}

	elif number == 10 : # The Euplotid Nuclear Code (transl_table=10)
		genetic_code = {'F': {'UUU': [], 'UUC': []}, 'S': {'UCU': [], 'UCC': [], 'UCA': [], 'UCG': [], 'AGU': [], 'AGC': [], 'AGA': [], 'AGG': []}, 'Y': {'UAU': [], 'UAC': []}, 'C': {'UGU': [], 'UGC': []}, 'L': {'UUA': [], 'UUG': [], 'CUU': [], 'CUC': [], 'CUA': [], 'CUG': []}, '*': {'UAA': [], 'UAG': []}, 'W': {'UGA': [], 'UGG': []}, 'P': {'CCU': [], 'CCC': [], 'CCA': [], 'CCG': []}, 'H': {'CAU': [], 'CAC': []}, 'R': {'CGU': [], 'CGC': [], 'CGA': [], 'CGG': []}, 'Q': {'CAA': [], 'CAG': []}, 'I': {'AUU': [], 'AUC': [], 'AUA': []}, 'T': {'ACU': [], 'ACC': [], 'ACA': [], 'ACG': []}, 'N': {'AAU': [], 'AAC': [], 'AAA': []}, 'M': {'AUG': []}, 'K': {'AAG': []}, 'V': {'GUU': [], 'GUC': [], 'GUA': [], 'GUG': []}, 'A': {'GCU': [], 'GCC': [], 'GCA': [], 'GCG': []}, 'D': {'GAU': [], 'GAC': []}, 'G': {'GGU': [], 'GGC': [], 'GGA': [], 'GGG': []}, 'E': {'GAA': [], 'GAG': []}}

	elif number == 11 : # The Bacterial "Code" (transl_table=11)
		genetic_code = {'F': {'UUU': [], 'UUC': []}, 'S': {'UCU': [], 'UCC': [], 'UCA': [], 'UCG': [], 'AGU': [], 'AGC': []}, 'Y': {'UAU': [], 'UAC': []}, 'C': {'UGU': [], 'UGC': []}, 'L': {'UUA': [], 'UUG': [], 'CUU': [], 'CUC': [], 'CUA': [], 'CUG': []}, '*': {'UAA': [], 'UGA': [], 'UAG': []}, 'W': {'UGG': []}, 'P': {'CCU': [], 'CCC': [], 'CCA': [], 'CCG': []}, 'H': {'CAU': [], 'CAC': []}, 'R': {'CGU': [], 'CGC': [], 'CGA': [], 'CGG': [], 'AGA': [], 'AGG': []}, 'Q': {'CAA': [], 'CAG': []}, 'I': {'AUU': [], 'AUC': [], 'AUA': []}, 'T': {'ACU': [], 'ACC': [], 'ACA': [], 'ACG': []}, 'N': {'AAU': [], 'AAC': []}, 'K': {'AAA': [], 'AAG': []}, 'M': {'AUG': []}, 'V': {'GUU': [], 'GUC': [], 'GUA': [], 'GUG': []}, 'A': {'GCU': [], 'GCC': [], 'GCA': [], 'GCG': []}, 'D': {'GAU': [], 'GAC': []}, 'G': {'GGU': [], 'GGC': [], 'GGA': [], 'GGG': []}, 'E': {'GAA': [], 'GAG': []}}

	elif number == 12 : # The Alternative Yeast Nuclear Code (transl_table=12)
		genetic_code = {'F': {'UUU': [], 'UUC': []}, 'S': {'UCU': [], 'UCC': [], 'UCA': [], 'UCG': [], 'CUG': [], 'AGU': [], 'AGC': []}, 'Y': {'UAU': [], 'UAC': []}, 'C': {'UGU': [], 'UGC': []}, 'L': {'UUA': [], 'UUG': [], 'CUU': [], 'CUC': [], 'CUA': []}, '*': {'UAA': [], 'UGA': [], 'UAG': []}, 'W': {'UGG': []}, 'P': {'CCU': [], 'CCC': [], 'CCA': [], 'CCG': []}, 'H': {'CAU': [], 'CAC': []}, 'R': {'CGU': [], 'CGC': [], 'CGA': [], 'CGG': [], 'AGA': [], 'AGG': []}, 'Q': {'CAA': [], 'CAG': []}, 'I': {'AUU': [], 'AUC': [], 'AUA': []}, 'T': {'ACU': [], 'ACC': [], 'ACA': [], 'ACG': []}, 'N': {'AAU': [], 'AAC': []}, 'K': {'AAA': [], 'AAG': []}, 'M': {'AUG': []}, 'V': {'GUU': [], 'GUC': [], 'GUA': [], 'GUG': []}, 'A': {'GCU': [], 'GCC': [], 'GCA': [], 'GCG': []}, 'D': {'GAU': [], 'GAC': []}, 'G': {'GGU': [], 'GGC': [], 'GGA': [], 'GGG': []}, 'E': {'GAA': [], 'GAG': []}}

	elif number == 13 : # The Ascidian Mitochondrial Code (transl_table=13)
		genetic_code = {'F': {'UUU': [], 'UUC': []}, 'S': {'UCU': [], 'UCC': [], 'UCA': [], 'UCG': [], 'AGU': [], 'AGC': []}, 'Y': {'UAU': [], 'UAC': []}, 'C': {'UGU': [], 'UGC': []}, 'L': {'UUA': [], 'UUG': [], 'CUU': [], 'CUC': [], 'CUA': [], 'CUG': []}, '*': {'UAA': [], 'UAG': []}, 'W': {'UGA': [], 'UGG': []}, 'P': {'CCU': [], 'CCC': [], 'CCA': [], 'CCG': []}, 'H': {'CAU': [], 'CAC': []}, 'R': {'CGU': [], 'CGC': [], 'CGA': [], 'CGG': []}, 'Q': {'CAA': [], 'CAG': []}, 'I': {'AUU': [], 'AUC': []}, 'T': {'ACU': [], 'ACC': [], 'ACA': [], 'ACG': []}, 'N': {'AAU': [], 'AAC': []}, 'M': {'AUA': [], 'AUG': []}, 'K': {'AAA': [], 'AAG': []}, 'G': {'AGA': [], 'AGG': [], 'GGU': [], 'GGC': [], 'GGA': [], 'GGG': []}, 'V': {'GUU': [], 'GUC': [], 'GUA': [], 'GUG': []}, 'A': {'GCU': [], 'GCC': [], 'GCA': [], 'GCG': []}, 'D': {'GAU': [], 'GAC': []}, 'E': {'GAA': [], 'GAG': []}}

	elif number == 14 : # 14. The Flatworm Mitochondrial Code (transl_table=14)
		genetic_code = {'F': {'UUU': [], 'UUC': []}, 'S': {'UCU': [], 'UCC': [], 'UCA': [], 'UCG': [], 'AGU': [], 'AGC': [], 'AGA': [], 'AGG': []}, 'Y': {'UAU': [], 'UAC': [], 'UAA': []}, 'C': {'UGU': [], 'UGC': []}, 'L': {'UUA': [], 'UUG': [], 'CUU': [], 'CUC': [], 'CUA': [], 'CUG': []}, 'W': {'UGA': [], 'UGG': []}, '*': {'UAG': []}, 'P': {'CCU': [], 'CCC': [], 'CCA': [], 'CCG': []}, 'H': {'CAU': [], 'CAC': []}, 'R': {'CGU': [], 'CGC': [], 'CGA': [], 'CGG': []}, 'Q': {'CAA': [], 'CAG': []}, 'I': {'AUU': [], 'AUC': [], 'AUA': []}, 'T': {'ACU': [], 'ACC': [], 'ACA': [], 'ACG': []}, 'N': {'AAU': [], 'AAC': [], 'AAA': []}, 'M': {'AUG': []}, 'K': {'AAG': []}, 'V': {'GUU': [], 'GUC': [], 'GUA': [], 'GUG': []}, 'A': {'GCU': [], 'GCC': [], 'GCA': [], 'GCG': []}, 'D': {'GAU': [], 'GAC': []}, 'G': {'GGU': [], 'GGC': [], 'GGA': [], 'GGG': []}, 'E': {'GAA': [], 'GAG': []}}

	elif number == 15 : # Blepharisma Nuclear Code (transl_table=15)
		genetic_code = {'F': {'UUU': [], 'UUC': []}, 'S': {'UCU': [], 'UCC': [], 'UCA': [], 'UCG': [], 'AGU': [], 'AGC': []}, 'Y': {'UAU': [], 'UAC': []}, 'C': {'UGU': [], 'UGC': []}, 'L': {'UUA': [], 'UUG': [], 'CUU': [], 'CUC': [], 'CUA': [], 'CUG': []}, '*': {'UAA': [], 'UGA': []}, 'Q': {'UAG': [], 'CAA': [], 'CAG': []}, 'W': {'UGG': []}, 'P': {'CCU': [], 'CCC': [], 'CCA': [], 'CCG': []}, 'H': {'CAU': [], 'CAC': []}, 'R': {'CGU': [], 'CGC': [], 'CGA': [], 'CGG': [], 'AGA': [], 'AGG': []}, 'I': {'AUU': [], 'AUC': [], 'AUA': []}, 'T': {'ACU': [], 'ACC': [], 'ACA': [], 'ACG': []}, 'N': {'AAU': [], 'AAC': []}, 'K': {'AAA': [], 'AAG': []}, 'M': {'AUG': []}, 'V': {'GUU': [], 'GUC': [], 'GUA': [], 'GUG': []}, 'A': {'GCU': [], 'GCC': [], 'GCA': [], 'GCG': []}, 'D': {'GAU': [], 'GAC': []}, 'G': {'GGU': [], 'GGC': [], 'GGA': [], 'GGG': []}, 'E': {'GAA': [], 'GAG': []}}
	
	elif number == 16 : # Chlorophycean Mitochondrial Code (transl_table=16)
		genetic_code = {'F': {'UUU': [], 'UUC': []}, 'S': {'UCU': [], 'UCC': [], 'UCA': [], 'UCG': [], 'AGU': [], 'AGC': []}, 'Y': {'UAU': [], 'UAC': []}, 'C': {'UGU': [], 'UGC': []}, 'L': {'UUA': [], 'UUG': [], 'CUU': [], 'CUC': [], 'CUA': [], 'CUG': [], 'UAG': []}, '*': {'UAA': [], 'UGA': []}, 'W': {'UGG': []}, 'P': {'CCU': [], 'CCC': [], 'CCA': [], 'CCG': []}, 'H': {'CAU': [], 'CAC': []}, 'R': {'CGU': [], 'CGC': [], 'CGA': [], 'CGG': [], 'AGA': [], 'AGG': []}, 'Q': {'CAA': [], 'CAG': []}, 'I': {'AUU': [], 'AUC': [], 'AUA': []}, 'T': {'ACU': [], 'ACC': [], 'ACA': [], 'ACG': []}, 'N': {'AAU': [], 'AAC': []}, 'K': {'AAA': [], 'AAG': []}, 'M': {'AUG': []}, 'V': {'GUU': [], 'GUC': [], 'GUA': [], 'GUG': []}, 'A': {'GCU': [], 'GCC': [], 'GCA': [], 'GCG': []}, 'D': {'GAU': [], 'GAC': []}, 'G': {'GGU': [], 'GGC': [], 'GGA': [], 'GGG': []}, 'E': {'GAA': [], 'GAG': []}}

	elif number == 21 : # Trematode Mitochondrial Code (transl_table=21)
		genetic_code = {'F': {'UUU': [], 'UUC': []}, 'S': {'UCU': [], 'UCC': [], 'UCA': [], 'UCG': [], 'AGU': [], 'AGC': [], 'AGA': [], 'AGG': []}, 'Y': {'UAU': [], 'UAC': []}, 'C': {'UGU': [], 'UGC': []}, 'L': {'UUA': [], 'UUG': [], 'CUU': [], 'CUC': [], 'CUA': [], 'CUG': []}, '*': {'UAA': [], 'UAG': []}, 'W': {'UGG': [], 'UGA': []}, 'P': {'CCU': [], 'CCC': [], 'CCA': [], 'CCG': []}, 'H': {'CAU': [], 'CAC': []}, 'R': {'CGU': [], 'CGC': [], 'CGA': [], 'CGG': []}, 'Q': {'CAA': [], 'CAG': []}, 'I': {'AUU': [], 'AUC': []}, 'T': {'ACU': [], 'ACC': [], 'ACA': [], 'ACG': []}, 'N': {'AAU': [], 'AAC': [],'AAA': []}, 'K': {'AAG': []}, 'M': {'AUG': [], 'AUA': []}, 'V': {'GUU': [], 'GUC': [], 'GUA': [], 'GUG': []}, 'A': {'GCU': [], 'GCC': [], 'GCA': [], 'GCG': []}, 'D': {'GAU': [], 'GAC': []}, 'G': {'GGU': [], 'GGC': [], 'GGA': [], 'GGG': []}, 'E': {'GAA': [], 'GAG': []}}

	elif number == 22 : # Scenedesmus obliquus Mitochondrial Code (transl_table=22)
		genetic_code = {'F': {'UUU': [], 'UUC': []}, 'S': {'UCU': [], 'UCC': [], 'UCG': [], 'AGU': [], 'AGC': []}, 'Y': {'UAU': [], 'UAC': []}, 'C': {'UGU': [], 'UGC': []}, 'L': {'UUA': [], 'UAG': [], 'UUG': [], 'CUU': [], 'CUC': [], 'CUA': [], 'CUG': []}, '*': {'UAA': [], 'UGA': [], 'UCA': []}, 'W': {'UGG': []}, 'P': {'CCU': [], 'CCC': [], 'CCA': [], 'CCG': []}, 'H': {'CAU': [], 'CAC': []}, 'R': {'CGU': [], 'CGC': [], 'CGA': [], 'CGG': [], 'AGA': [], 'AGG': []}, 'Q': {'CAA': [], 'CAG': []}, 'I': {'AUU': [], 'AUC': [], 'AUA': []}, 'T': {'ACU': [], 'ACC': [], 'ACA': [], 'ACG': []}, 'N': {'AAU': [], 'AAC': []}, 'K': {'AAA': [], 'AAG': []}, 'M': {'AUG': []}, 'V': {'GUU': [], 'GUC': [], 'GUA': [], 'GUG': []}, 'A': {'GCU': [], 'GCC': [], 'GCA': [], 'GCG': []}, 'D': {'GAU': [], 'GAC': []}, 'G': {'GGU': [], 'GGC': [], 'GGA': [], 'GGG': []}, 'E': {'GAA': [], 'GAG': []}}

	elif number == 23 : # Thraustochytrium Mitochondrial Code (transl_table=23)
		genetic_code = {'F': {'UUU': [], 'UUC': []}, 'S': {'UCU': [], 'UCC': [], 'UCA': [], 'UCG': [], 'AGU': [], 'AGC': []}, 'Y': {'UAU': [], 'UAC': []}, 'C': {'UGU': [], 'UGC': []}, 'L': {'UUG': [], 'CUU': [], 'CUC': [], 'CUA': [], 'CUG': []}, '*': {'UUA': [],'UAA': [], 'UGA': [], 'UAG': []}, 'W': {'UGG': []}, 'P': {'CCU': [], 'CCC': [], 'CCA': [], 'CCG': []}, 'H': {'CAU': [], 'CAC': []}, 'R': {'CGU': [], 'CGC': [], 'CGA': [], 'CGG': [], 'AGA': [], 'AGG': []}, 'Q': {'CAA': [], 'CAG': []}, 'I': {'AUU': [], 'AUC': [], 'AUA': []}, 'T': {'ACU': [], 'ACC': [], 'ACA': [], 'ACG': []}, 'N': {'AAU': [], 'AAC': []}, 'K': {'AAA': [], 'AAG': []}, 'M': {'AUG': []}, 'V': {'GUU': [], 'GUC': [], 'GUA': [], 'GUG': []}, 'A': {'GCU': [], 'GCC': [], 'GCA': [], 'GCG': []}, 'D': {'GAU': [], 'GAC': []}, 'G': {'GGU': [], 'GGC': [], 'GGA': [], 'GGG': []}, 'E': {'GAA': [], 'GAG': []}}

	elif number == 24 : # Pterobranchia Mitochondrial Code (transl_table=24)
		genetic_code = {'F': {'UUU': [], 'UUC': []}, 'S': {'AGA': [], 'UCU': [], 'UCC': [], 'UCA': [], 'UCG': [], 'AGU': [], 'AGC': []}, 'Y': {'UAU': [], 'UAC': []}, 'C': {'UGU': [], 'UGC': []}, 'L': {'UUA': [], 'UUG': [], 'CUU': [], 'CUC': [], 'CUA': [], 'CUG': []}, '*': {'UAA': [], 'UAG': []}, 'W': {'UGA': [], 'UGG': []}, 'P': {'CCU': [], 'CCC': [], 'CCA': [], 'CCG': []}, 'H': {'CAU': [], 'CAC': []}, 'R': {'CGU': [], 'CGC': [], 'CGA': [], 'CGG': []}, 'Q': {'CAA': [], 'CAG': []}, 'I': {'AUU': [], 'AUC': [], 'AUA': []}, 'T': {'ACU': [], 'ACC': [], 'ACA': [], 'ACG': []}, 'N': {'AAU': [], 'AAC': []}, 'K': {'AGG': [],'AAA': [], 'AAG': []}, 'M': {'AUG': []}, 'V': {'GUU': [], 'GUC': [], 'GUA': [], 'GUG': []}, 'A': {'GCU': [], 'GCC': [], 'GCA': [], 'GCG': []}, 'D': {'GAU': [], 'GAC': []}, 'G': {'GGU': [], 'GGC': [], 'GGA': [], 'GGG': []}, 'E': {'GAA': [], 'GAG': []}}

	elif number == 25 : # Candidate Division SR1 and Gracilibacteria Code (transl_table=25)
		genetic_code = {'F': {'UUU': [], 'UUC': []}, 'S': {'UCU': [], 'UCC': [], 'UCA': [], 'UCG': [], 'AGU': [], 'AGC': []}, 'Y': {'UAU': [], 'UAC': []}, 'C': {'UGU': [], 'UGC': []}, 'L': {'UUA': [], 'UUG': [], 'CUU': [], 'CUC': [], 'CUA': [], 'CUG': []}, '*': {'UAA': [], 'UAG': []}, 'W': {'UGG': []}, 'P': {'CCU': [], 'CCC': [], 'CCA': [], 'CCG': []}, 'H': {'CAU': [], 'CAC': []}, 'R': {'CGU': [], 'CGC': [], 'CGA': [], 'CGG': [], 'AGA': [], 'AGG': []}, 'Q': {'CAA': [], 'CAG': []}, 'I': {'AUU': [], 'AUC': [], 'AUA': []}, 'T': {'ACU': [], 'ACC': [], 'ACA': [], 'ACG': []}, 'N': {'AAU': [], 'AAC': []}, 'K': {'AAA': [], 'AAG': []}, 'M': {'AUG': []}, 'V': {'GUU': [], 'GUC': [], 'GUA': [], 'GUG': []}, 'A': {'GCU': [], 'GCC': [], 'GCA': [], 'GCG': []}, 'D': {'GAU': [], 'GAC': []}, 'G': {'UGA': [], 'GGU': [], 'GGC': [], 'GGA': [], 'GGG': []}, 'E': {'GAA': [], 'GAG': []}}
	
	elif number == 26 : # Pachysolen tannophilus Nuclear Code (transl_table=26)
		genetic_code = {'F': {'UUU': [], 'UUC': []}, 'S': {'UCU': [], 'UCC': [], 'UCA': [], 'UCG': [], 'AGU': [], 'AGC': []}, 'Y': {'UAU': [], 'UAC': []}, 'C': {'UGU': [], 'UGC': []}, 'L': {'UUA': [], 'UUG': [], 'CUU': [], 'CUC': [], 'CUA': []}, '*': {'UAA': [], 'UGA': [], 'UAG': []}, 'W': {'UGG': []}, 'P': {'CCU': [], 'CCC': [], 'CCA': [], 'CCG': []}, 'H': {'CAU': [], 'CAC': []}, 'R': {'CGU': [], 'CGC': [], 'CGA': [], 'CGG': [], 'AGA': [], 'AGG': []}, 'Q': {'CAA': [], 'CAG': []}, 'I': {'AUU': [], 'AUC': [], 'AUA': []}, 'T': {'ACU': [], 'ACC': [], 'ACA': [], 'ACG': []}, 'N': {'AAU': [], 'AAC': []}, 'K': {'AAA': [], 'AAG': []}, 'M': {'AUG': []}, 'V': {'GUU': [], 'GUC': [], 'GUA': [], 'GUG': []}, 'A': {'CUG': [], 'GCU': [], 'GCC': [], 'GCA': [], 'GCG': []}, 'D': {'GAU': [], 'GAC': []}, 'G': {'GGU': [], 'GGC': [], 'GGA': [], 'GGG': []}, 'E': {'GAA': [], 'GAG': []}}	
	
	elif number == 27 : # Karyorelict Nuclear Code (transl_table=27) REFAIRE A PARTIR DE LA
		genetic_code = {'F': {'UUU': [], 'UUC': []}, 'S': {'UCU': [], 'UCC': [], 'UCA': [], 'UCG': [], 'AGU': [], 'AGC': []}, 'Y': {'UAU': [], 'UAC': []}, 'C': {'UGU': [], 'UGC': []}, 'L': {'UUA': [], 'UUG': [], 'CUU': [], 'CUC': [], 'CUA': [], 'CUG': []}, '*': {}, 'W': {'UGA': [], 'UGG': []}, 'P': {'CCU': [], 'CCC': [], 'CCA': [], 'CCG': []}, 'H': {'CAU': [], 'CAC': []}, 'R': {'CGU': [], 'CGC': [], 'CGA': [], 'CGG': [], 'AGA': [], 'AGG': []}, 'Q': {'UAA': [], 'UAG': [], 'CAA': [], 'CAG': []}, 'I': {'AUU': [], 'AUC': [], 'AUA': []}, 'T': {'ACU': [], 'ACC': [], 'ACA': [], 'ACG': []}, 'N': {'AAU': [], 'AAC': []}, 'K': {'AAA': [], 'AAG': []}, 'M': {'AUG': []}, 'V': {'GUU': [], 'GUC': [], 'GUA': [], 'GUG': []}, 'A': {'GCU': [], 'GCC': [], 'GCA': [], 'GCG': []}, 'D': {'GAU': [], 'GAC': []}, 'G': {'GGU': [], 'GGC': [], 'GGA': [], 'GGG': []}, 'E': {'GAA': [], 'GAG': []}}	

	elif number == 28 : # Condylostoma Nuclear Code (transl_table=28)
		genetic_code = {'F': {'UUU': [], 'UUC': []}, 'S': {'UCU': [], 'UCC': [], 'UCA': [], 'UCG': [], 'AGU': [], 'AGC': []}, 'Y': {'UAU': [], 'UAC': []}, 'C': {'UGU': [], 'UGC': []}, 'L': {'UUA': [], 'UUG': [], 'CUU': [], 'CUC': [], 'CUA': [], 'CUG': []}, '*': {}, 'W': {'UGA': [], 'UGG': []}, 'P': {'CCU': [], 'CCC': [], 'CCA': [], 'CCG': []}, 'H': {'CAU': [], 'CAC': []}, 'R': {'CGU': [], 'CGC': [], 'CGA': [], 'CGG': [], 'AGA': [], 'AGG': []}, 'Q': {'UAA': [], 'UAG': [], 'CAA': [], 'CAG': []}, 'I': {'AUU': [], 'AUC': [], 'AUA': []}, 'T': {'ACU': [], 'ACC': [], 'ACA': [], 'ACG': []}, 'N': {'AAU': [], 'AAC': []}, 'K': {'AAA': [], 'AAG': []}, 'M': {'AUG': []}, 'V': {'GUU': [], 'GUC': [], 'GUA': [], 'GUG': []}, 'A': {'GCU': [], 'GCC': [], 'GCA': [], 'GCG': []}, 'D': {'GAU': [], 'GAC': []}, 'G': {'GGU': [], 'GGC': [], 'GGA': [], 'GGG': []}, 'E': {'GAA': [], 'GAG': []}}	

	elif number == 29 : # Mesodinium Nuclear Code (transl_table=29)
		genetic_code = {'F': {'UUU': [], 'UUC': []}, 'S': {'UCU': [], 'UCC': [], 'UCA': [], 'UCG': [], 'AGU': [], 'AGC': []}, 'Y': {'UAA': [], 'UAG': [],'UAU': [], 'UAC': []}, 'C': {'UGU': [], 'UGC': []}, 'L': {'UUA': [], 'UUG': [], 'CUU': [], 'CUC': [], 'CUA': [], 'CUG': []}, '*': {'UGA': []}, 'W': {'UGG': []}, 'P': {'CCU': [], 'CCC': [], 'CCA': [], 'CCG': []}, 'H': {'CAU': [], 'CAC': []}, 'R': {'CGU': [], 'CGC': [], 'CGA': [], 'CGG': [], 'AGA': [], 'AGG': []}, 'Q': {'CAA': [], 'CAG': []}, 'I': {'AUU': [], 'AUC': [], 'AUA': []}, 'T': {'ACU': [], 'ACC': [], 'ACA': [], 'ACG': []}, 'N': {'AAU': [], 'AAC': []}, 'K': {'AAA': [], 'AAG': []}, 'M': {'AUG': []}, 'V': {'GUU': [], 'GUC': [], 'GUA': [], 'GUG': []}, 'A': {'GCU': [], 'GCC': [], 'GCA': [], 'GCG': []}, 'D': {'GAU': [], 'GAC': []}, 'G': {'GGU': [], 'GGC': [], 'GGA': [], 'GGG': []}, 'E': {'GAA': [], 'GAG': []}}

	elif number == 30 : # Peritrich Nuclear Code (transl_table=30)
		genetic_code = {'F': {'UUU': [], 'UUC': []}, 'S': {'UCU': [], 'UCC': [], 'UCA': [], 'UCG': [], 'AGU': [], 'AGC': []}, 'Y': {'UAU': [], 'UAC': []}, 'C': {'UGU': [], 'UGC': []}, 'L': {'UUA': [], 'UUG': [], 'CUU': [], 'CUC': [], 'CUA': [], 'CUG': []}, '*': {'UGA': []}, 'W': {'UGG': []}, 'P': {'CCU': [], 'CCC': [], 'CCA': [], 'CCG': []}, 'H': {'CAU': [], 'CAC': []}, 'R': {'CGU': [], 'CGC': [], 'CGA': [], 'CGG': [], 'AGA': [], 'AGG': []}, 'Q': {'UAA': [], 'UAG': [], 'CAA': [], 'CAG': []}, 'I': {'AUU': [], 'AUC': [], 'AUA': []}, 'T': {'ACU': [], 'ACC': [], 'ACA': [], 'ACG': []}, 'N': {'AAU': [], 'AAC': []}, 'K': {'AAA': [], 'AAG': []}, 'M': {'AUG': []}, 'V': {'GUU': [], 'GUC': [], 'GUA': [], 'GUG': []}, 'A': {'GCU': [], 'GCC': [], 'GCA': [], 'GCG': []}, 'D': {'GAU': [], 'GAC': []}, 'G': {'GGU': [], 'GGC': [], 'GGA': [], 'GGG': []}, 'E': {'GAA': [], 'GAG': []}}	

	elif number == 31 : # Blastocrithidia Nuclear Code (transl_table=31)
		genetic_code = {'F': {'UUU': [], 'UUC': []}, 'S': {'UCU': [], 'UCC': [], 'UCA': [], 'UCG': [], 'AGU': [], 'AGC': []}, 'Y': {'UAU': [], 'UAC': []}, 'C': {'UGU': [], 'UGC': []}, 'L': {'UUA': [], 'UUG': [], 'CUU': [], 'CUC': [], 'CUA': [], 'CUG': []}, '*': {}, 'W': {'UGG': [],'UGG': []}, 'P': {'CCU': [], 'CCC': [], 'CCA': [], 'CCG': []}, 'H': {'CAU': [], 'CAC': []}, 'R': {'CGU': [], 'CGC': [], 'CGA': [], 'CGG': [], 'AGA': [], 'AGG': []}, 'Q': {'CAA': [], 'CAG': []}, 'I': {'AUU': [], 'AUC': [], 'AUA': []}, 'T': {'ACU': [], 'ACC': [], 'ACA': [], 'ACG': []}, 'N': {'AAU': [], 'AAC': []}, 'K': {'AAA': [], 'AAG': []}, 'M': {'AUG': []}, 'V': {'GUU': [], 'GUC': [], 'GUA': [], 'GUG': []}, 'A': {'GCU': [], 'GCC': [], 'GCA': [], 'GCG': []}, 'D': {'GAU': [], 'GAC': []}, 'G': {'GGU': [], 'GGC': [], 'GGA': [], 'GGG': []}, 'E': {'UAA': [], 'UAG': [], 'GAA': [], 'GAG': []}}	

	elif number == 33 : # Cephalodiscidae Mitochondrial UAA-Tyr Code (transl_table=33)
		genetic_code = {'F': {'UUU': [], 'UUC': []}, 'S': {'AGA': [], 'UCU': [], 'UCC': [], 'UCA': [], 'UCG': [], 'AGU': [], 'AGC': []}, 'Y': {'UAA': [], 'UAU': [], 'UAC': []}, 'C': {'UGU': [], 'UGC': []}, 'L': {'UUA': [], 'UUG': [], 'CUU': [], 'CUC': [], 'CUA': [], 'CUG': []}, '*': {'UAG': []}, 'W': {'UGA': [], 'UGG': []}, 'P': {'CCU': [], 'CCC': [], 'CCA': [], 'CCG': []}, 'H': {'CAU': [], 'CAC': []}, 'R': {'CGU': [], 'CGC': [], 'CGA': [], 'CGG': []}, 'Q': {'CAA': [], 'CAG': []}, 'I': {'AUU': [], 'AUC': [], 'AUA': []}, 'T': {'ACU': [], 'ACC': [], 'ACA': [], 'ACG': []}, 'N': {'AAU': [], 'AAC': []}, 'K': {'AGG': [], 'AAA': [], 'AAG': []}, 'M': {'AUG': []}, 'V': {'GUU': [], 'GUC': [], 'GUA': [], 'GUG': []}, 'A': {'GCU': [], 'GCC': [], 'GCA': [], 'GCG': []}, 'D': {'GAU': [], 'GAC': []}, 'G': {'GGU': [], 'GGC': [], 'GGA': [], 'GGG': []}, 'E': {'GAA': [], 'GAG': []}}			
	
	else :
		log_file.write("## GENETIC CODE ERROR, you must enter a valid number for your genetic code (1~6 / 9~16 / 21~31 / 33 ##")
		print ("## GENETIC CODE ERROR, you must enter a valid number for your genetic code (1~6 / 9~16 / 21~31 / 33, exiting... ##") 
		sys.exit()
	return genetic_code

#To check
################################################################################
########################## CU TABLE PARSER #####################################
################################################################################

def get_CU_table_content(CU_table, genetic_code) :

	global occ_tot
	global occ_tot_aa
	global occ_tot_aa_4_sim
	CU_codon = []
	CU_table_line = []

	split_line = re.compile('\n').split(CU_table)
	if split_line != None :
		i = 0
		while i < len(split_line) :
			if split_line[i] != '' :
				CU_table_line.append(split_line[i].replace("T","U"))
			i += 1

	i = 0
	while i < len(CU_table_line) :
		split_codon = re.findall('([A-z]{3}[^\(]+\({1}\s*[0-9]*\){1})', CU_table_line[i])
		if split_codon :
			j = 0
			while j <len(split_codon) :
				CU_codon.append(split_codon[j])
				j += 1
		i += 1

	i = 0
	while i < len(CU_codon) :
		codon = re.search('(^[A-z]{3})', CU_codon[i]).group(1)
		occ_codon = re.search('\(\s*([0-9]*)\)', CU_codon[i]).group(1)
		for amino_acid in genetic_code :
			if codon in genetic_code[amino_acid] :
				genetic_code[amino_acid][codon].append(int(occ_codon))
		i += 1

	for amino_acid in genetic_code :
		occ_tot_aa[amino_acid] = 0
		for codon in genetic_code[amino_acid] :
			if amino_acid not in occ_tot_aa :
				occ_tot_aa[amino_acid] = genetic_code[amino_acid][codon][0]
			else : 
				occ_tot_aa[amino_acid] += genetic_code[amino_acid][codon][0]
		if amino_acid not in occ_tot_ref :
			occ_tot_ref[amino_acid] = occ_tot_aa[amino_acid] # pour la fonction unvoid_empty_frequencies
		occ_tot += occ_tot_aa[amino_acid]
		if amino_acid != '*' :
			occ_tot_aa_4_sim += occ_tot_aa[amino_acid]
		for codon in genetic_code[amino_acid] :
			genetic_code[amino_acid][codon].append(round((genetic_code[amino_acid][codon][0] / occ_tot_aa[amino_acid]), 3))


#checked

################################################################################
################# FUNCTION USED TO AVOID 0 VALUES IN CU TABLE ##################
################################################################################

def unvoid_empty_frequencies(genetic_code) :
	for amino_acid in genetic_code :
		cpt = 0
		for codon in genetic_code[amino_acid] :
			if genetic_code[amino_acid][codon][0] == 0 and genetic_code[amino_acid][codon][1] == 0 :
				cpt += 1
				genetic_code[amino_acid][codon][0] = (1 / (61 * occ_tot))
				genetic_code[amino_acid][codon][1] = (1 / (61 * occ_tot))
		for codon in genetic_code[amino_acid] :
			if genetic_code[amino_acid][codon][1] > (1 / (61 * occ_tot)) and genetic_code[amino_acid][codon][0] > (1 / (61 * occ_tot)) :
				genetic_code[amino_acid][codon][1] = round(genetic_code[amino_acid][codon][1]  - (( cpt * (1 / (61 * occ_tot)))/(len(genetic_code[amino_acid]) - cpt )), 3)
				
#Seems okay

################################################################################
############ FUNCTION USED TO GET REVERSE FREQUENCIES ON CU TABLE ##############
################################################################################

def get_reverse_frequencies(genetic_code) :
	genetic_code_rev = {}
	genetic_code_rev = deepcopy(genetic_code)
	for amino_acid in genetic_code_rev :
		denominator = 0
		for codon in genetic_code_rev[amino_acid] :
			denominator += (1/genetic_code[amino_acid][codon][1])
		for codon in genetic_code_rev[amino_acid] :
			genetic_code_rev[amino_acid][codon][1] = round((((1/genetic_code_rev[amino_acid][codon][1]) / denominator)), 3)
		for codon in genetic_code_rev[amino_acid] :
			genetic_code_rev[amino_acid][codon][0] = round(genetic_code_rev[amino_acid][codon][1] * occ_tot_ref[amino_acid])

	return genetic_code_rev

#To check
# NOT GOOD, ROUND FUNCTION DOES NOT WORK AT 2.5

################################################################################
## FUNCTION USED TO GET CODONS WITH GC-AT MAX CONTENT FREQUENCIES ON CU TABLE ##
################################################################################

def get_GC_and_AT_max_codons(genetic_code) :

	for amino_acid in genetic_code :
		ngc = 0
		nat = 0
		total_number_gc = 0
		total_number_at = 0
		for codon in genetic_code[amino_acid] :

			nb_gc = codon.count('G') + codon.count('C')
			nb_at = codon.count('A') + codon.count('U')

			if nb_gc > ngc :
				codon_gc_max[amino_acid] = {}
				codon_gc_max[amino_acid][codon] = []
				ngc = nb_gc
			elif nb_gc == ngc :
				if amino_acid not in codon_gc_max :
					codon_gc_max[amino_acid] = {}
				codon_gc_max[amino_acid][codon] = []

			if nb_at > nat :
				codon_at_max[amino_acid] = {}
				codon_at_max[amino_acid][codon] = []
				nat = nb_at
			elif nb_at == nat :
				if amino_acid not in codon_at_max :
					codon_at_max[amino_acid] = {}
				codon_at_max[amino_acid][codon] = []

		for codon in genetic_code[amino_acid] :
			for codon_gc in codon_gc_max[amino_acid] :
				if codon == codon_gc :
					codon_gc_max[amino_acid][codon_gc].append(genetic_code[amino_acid][codon][0])
					total_number_gc += codon_gc_max[amino_acid][codon_gc][0]
		for codon_gc in codon_gc_max[amino_acid] :
			if total_number_gc <= 0 :
				codon_gc_max[amino_acid][codon_gc].append(0)
			else :
				codon_gc_max[amino_acid][codon_gc].append(codon_gc_max[amino_acid][codon_gc][0] / total_number_gc)
				codon_gc_max[amino_acid][codon_gc].append((1 - (codon_gc_max[amino_acid][codon_gc][0] / total_number_gc)))

		for codon in genetic_code[amino_acid] :
			for codon_at in codon_at_max[amino_acid] :
				if codon == codon_at :
					codon_at_max[amino_acid][codon_at].append(genetic_code[amino_acid][codon][0])
					total_number_at += codon_at_max[amino_acid][codon_at][0]
		for codon_at in codon_at_max[amino_acid] :
			if total_number_at <= 0 :
				codon_at_max[amino_acid][codon_at].append(0)
			else :
				codon_at_max[amino_acid][codon_at].append(codon_at_max[amino_acid][codon_at][0] / total_number_at)
				codon_at_max[amino_acid][codon_at].append((1 - (codon_at_max[amino_acid][codon_at][0] / total_number_at)))

#TO CHECK BUT PROBABLY NEED A SPECIFIC DICT FOR GENETIC CODE.

################################################################################
######## FUNCTION USED TO TRANSFORM A STRING SEQUENCE IN A DICTIONNARY ONE #####
################################################################################

def seq_to_dict(sequence, header, genetic_code_ref, output_dir) :

	total_by_amino_acid = {}

	if not os.path.isdir("./" + output_dir + "/log"):
		os.makedirs("./" + output_dir + "/log")

	log_file = open("./" + output_dir + "/log/log_cousin.txt", "a")

	genetic_code_seq = dict.fromkeys(genetic_code_ref)
	for amino_acid in genetic_code_seq :
		genetic_code_seq[amino_acid] = dict.fromkeys(genetic_code_ref[amino_acid])
		for codon in genetic_code_seq[amino_acid] :
			genetic_code_seq[amino_acid][codon] = []
			genetic_code_seq[amino_acid][codon].append(0)


	i = 0
	length_seq = len(sequence) / 3
	cpt = 0
	while i < length_seq :
		codon_seq = sequence[:3]
		sequence  = sequence[3:]
		for amino_acid in genetic_code_seq :
			if codon_seq in genetic_code_seq[amino_acid] :
				if amino_acid == "*" and len(sequence) > 0 :
					log_file.write("## SEQUENCE ERROR : the sequence " + header  + " has a STOP codon in another place than at the end of it (check your file) ##\n")
					try :
						del (content_seq[header])
						return 'deleted'
					except :
						log_file.write("## Oops, this data seems already deleted / isn't part of your original dataset. ##\n")
						return 'deleted'
				genetic_code_seq[amino_acid][codon_seq][0] += 1
				if amino_acid not in total_by_amino_acid :
					total_by_amino_acid[amino_acid] = 1
				else :
					total_by_amino_acid[amino_acid] += 1
				cpt += 1
		i += 1

	for amino_acid in genetic_code_seq :
		for codon in genetic_code_seq[amino_acid] :
			if amino_acid in total_by_amino_acid :
				genetic_code_seq[amino_acid][codon].append(genetic_code_seq[amino_acid][codon][0] / total_by_amino_acid[amino_acid])
			else :
				genetic_code_seq[amino_acid][codon].append(0)

	log_file.close()
	return genetic_code_seq

#TO CHECK

################################################################################
################################################################################
################################################################################
############## SIMULATION OF SEQUENCES (thanks to CU table) ####################
################################################################################
################################################################################
################################################################################

################################################################################
########## FUNCTION USED TO CREATE SIMULATED SEQUENCES (random guided) #########
################################################################################

def create_sequences_from_CU_table(size, genetic_code_simulated, genetic_code_ref) :


	i = 0
	total_by_amino_acid = {}

	global occ_tot
	while i < size :
		separateur_aa = 0
		random_choice_aa = uniform(0, 1)
		for amino_acid in genetic_code_ref :
			if amino_acid != '*' :
				separateur_aa += occ_tot_ref[amino_acid] / occ_tot_aa_4_sim # MODIFIER POUR QU'ON AIT UN BON RATIO DANS LE CADRE D'UNE AUTRE TABLE. APPAREMMENT CA VA CAR LE NOMBRE d'ACIDES AMINES NE DOIT PAS CHANGER
				if random_choice_aa <= separateur_aa :
					separateur_codon = 0
					random_choice_codon = uniform(0, 1)
					for codon in genetic_code_ref[amino_acid] :
						separateur_codon += genetic_code_ref[amino_acid][codon][1]
						if random_choice_codon <= separateur_codon :
							genetic_code_simulated[amino_acid][codon][0] += 1
							break
					break
		i += 1

	for amino_acid in genetic_code_simulated :
		for codon in genetic_code_simulated[amino_acid] :
			if amino_acid not in total_by_amino_acid :
				total_by_amino_acid[amino_acid] = genetic_code_simulated[amino_acid][codon][0]
			else :
				total_by_amino_acid[amino_acid] += genetic_code_simulated[amino_acid][codon][0]

	for amino_acid in genetic_code_simulated :
		for codon in genetic_code_simulated[amino_acid] :
			if total_by_amino_acid[amino_acid] > 0 :
				if len(genetic_code_simulated[amino_acid][codon]) > 1 :
					genetic_code_simulated[amino_acid][codon][1] = (genetic_code_simulated[amino_acid][codon][0] / total_by_amino_acid[amino_acid])
				else :
					genetic_code_simulated[amino_acid][codon].append(genetic_code_simulated[amino_acid][codon][0] / total_by_amino_acid[amino_acid])
			else :
				if len(genetic_code_simulated[amino_acid][codon]) > 1 :
					genetic_code_simulated[amino_acid][codon][1] = 0
				else :
					genetic_code_simulated[amino_acid][codon].append(0)

	return genetic_code_simulated

################################################################################
################################################################################
################################################################################
############################# CALCULATION METHODS  #############################
################################################################################
################################################################################
################################################################################

################################################################################
############## FUNCTION USED TO GET GC CONTENT OF A SEQUENCE  ##################
################################################################################

def get_freq_GC_seq(sequence) :
	i = 0
	GC_percent = []
	nb_1st = 0
	nb_2nd = 0
	nb_3th = 0
	nb_gc = sequence.count('G') + sequence.count('C')
	GC_percent.append(100*(nb_gc / len(sequence)))

	seq_length = (len(sequence)/3)
	while i < seq_length :
		codon_seq = sequence[:3]
		sequence = sequence[3:]
		nb_1st += codon_seq[0].count('G') + codon_seq[0].count('C')
		nb_2nd += codon_seq[1].count('G') + codon_seq[1].count('C')
		nb_3th += codon_seq[2].count('G') + codon_seq[2].count('C')
		i += 1
	GC_percent.append(100*(nb_1st / (seq_length)))
	GC_percent.append(100*(nb_2nd / (seq_length)))
	GC_percent.append(100*(nb_3th / (seq_length)))

	return GC_percent

################################################################################
### FUNCTION USED TO GET ATGC CONTENT OF A SEQUENCE (only 3rd base returned) ###
################################################################################

def get_freq_ATGC_seq(sequence) :
	i = 0
	ATGC_percent = []
	nb_A = 0
	nb_T = 0
	nb_G = 0
	nb_C = 0

	seq_length = (len(sequence)/3)
	while i < seq_length :
		codon_seq = sequence[:3]
		sequence = sequence[3:]
		nb_A += codon_seq[2].count('A')
		nb_T += codon_seq[2].count('T')
		nb_T += codon_seq[2].count('U')
		nb_G += codon_seq[2].count('G')
		nb_C += codon_seq[2].count('C')
		i += 1
	ATGC_percent.append(100*(nb_A / (seq_length)))
	ATGC_percent.append(100*(nb_T / (seq_length)))
	ATGC_percent.append(100*(nb_G / (seq_length)))
	ATGC_percent.append(100*(nb_C / (seq_length)))

	return ATGC_percent

################################################################################
####### FUNCTION USED TO GET THE FREQUENCES OF THE SIMULATED SEQUENCE  #########
################################################################################

# def calc_freq_codon_req(genetic_code) : #VERIFIER SI FREQUENCES BIEN CONSTRUITES
#
#     freq_codon_req = {}
#
#     for amino_acid in genetic_code :
#         total_occ = 0
#         for codon in genetic_code[amino_acid] :
#                 total_occ += genetic_code[amino_acid][codon][0]
#         if total_occ == 0 :
#             for codon in genetic_code[amino_acid] :
#                 freq_codon_req[codon] = 0
#         else :
#             for codon in genetic_code[amino_acid] :
#                 if total_occ > 0 :
#                     freq_codon_req[codon] = genetic_code[amino_acid][codon][0] / total_occ
#                 else :
#                     freq_codon_req[codon] = 0 # A VERIFIER ABSOLUMENT
#
#     return freq_codon_req

def check_if_aa_in_req(genetic_code) : #VERIFIER SI FREQUENCES BIEN CONSTRUITES

	is_aa_in_req = {}

	for amino_acid in genetic_code :
		total_occ = 0
		for codon in genetic_code[amino_acid] :
				total_occ += genetic_code[amino_acid][codon][0]
		if total_occ == 0 :
			is_aa_in_req[amino_acid] = False
		else :
			is_aa_in_req[amino_acid] = True

	return is_aa_in_req

################################################################################
##FUNCTION USED TO GET OCCURENCES OF EACH AMINO ACID IN THE SIMULATED SEQUENCE##
################################################################################

def get_occ_aa (genetic_code, dict_occ_aa) :
	for amino_acid in genetic_code :
		dict_occ_aa[amino_acid] = 0

	for amino_acid in genetic_code :
		for codon in genetic_code[amino_acid] :
			if amino_acid not in dict_occ_aa :
				dict_occ_aa[amino_acid] = genetic_code[amino_acid][codon][0]
			else :
				dict_occ_aa[amino_acid] += genetic_code[amino_acid][codon][0]

def get_freq_aa (genetic_code, dict_occ_aa) : #Voir si fusionner avec la fonction d'au dessus, car on fait presque la même chose sauf qu'on rallonge le temps de calcul (donc voir si on garde la fonction de calcul des distances euclidiennes pour fusionner)
	total_aa = 0
	for amino_acid in genetic_code :
		if amino_acid != '*' :
			dict_occ_aa[amino_acid] = 0

	for amino_acid in genetic_code :
		if amino_acid != '*' :
			for codon in genetic_code[amino_acid] :
				if amino_acid not in dict_occ_aa :
					dict_occ_aa[amino_acid] = genetic_code[amino_acid][codon][0]
					total_aa += genetic_code[amino_acid][codon][0]
				else :
					dict_occ_aa[amino_acid] += genetic_code[amino_acid][codon][0]
					total_aa += genetic_code[amino_acid][codon][0]

	for amino_acid in dict_occ_aa :
		dict_occ_aa[amino_acid] = dict_occ_aa[amino_acid] / total_aa

################################################################################
################################################################################
################################################################################
############################### OPTIMIZER ######################################
################################################################################
################################################################################
################################################################################

################################################################################
############ FUNCTION USED TO TRANSFORM DNA SEQUENCES IN RNA ONES  #############
################################################################################

def DNA_to_RNA(sequence) :
	sequence = sequence.replace("T", "U")
	return sequence

################################################################################
############ FUNCTION USED TO TRANSFORM RNA SEQUENCES IN DNA ONES  #############
################################################################################

def RNA_to_DNA(sequence) :
	sequence = sequence.replace("U", "T")
	return sequence

################################################################################
############ FUNCTION USED TO TRANSFORM RNA SEQUENCES IN AA ONES  ##############
################################################################################

def RNA_to_AA(sequence, genetic_code) : # Besoin du code génétique de la référence pour construire la séquence en acides aminés
	seq_aa = ''
	seq_length = (len(sequence))
	i = 0
	while i < seq_length :
		codon_seq = sequence[:3]
		sequence = sequence[3:]
		for amino_acid in genetic_code :
			if codon_seq in genetic_code[amino_acid] :
				seq_aa += amino_acid
		i += 1

	return seq_aa

################################################################################
############ FUNCTION USED TO OPTIMIZE SEQUENCE IN A RANDOM ONE  ###############
################################################################################

def random_simulation(sequence, genetic_code) :
	seq_simulated = ''
	for c in sequence :
		for amino_acid in genetic_code :
			if c == amino_acid :
				random_choice = uniform(0, 1)
				separateur = 1/len(genetic_code[amino_acid])
				n = 1
				while n <= len(genetic_code[amino_acid]) :
					if random_choice <= (n*separateur) and random_choice > ((n-1)*separateur) :
						random_select = n-1
					n += 1

				list_codon = list(genetic_code[amino_acid])
				choice = list_codon[random_select]

				seq_simulated += choice
	return seq_simulated

################################################################################
######## FUNCTION USED TO OPTIMIZE SEQUENCE IN A RANDOM-GUIDED ONE  ############
################################################################################

def random_guided_simulation(sequence, genetic_code) :

	seq_simulated = ''
	for c in sequence :
		for amino_acid in genetic_code :
			separateur = 0
			if c == amino_acid :
				random_choice = uniform(0, 1)
				for codon in genetic_code[amino_acid] :
					separateur += genetic_code[amino_acid][codon][1]
					if random_choice <= separateur :
						seq_simulated += codon
						break
	return seq_simulated

################################################################################
##### FUNCTION USED TO OPTIMIZE SEQUENCE IN A RANDOM ONE (MAX GC CONTENT)  #####
################################################################################

def gc_max_random_simulation(sequence, genetic_code) :
	seq_simulated = ''
	for c in sequence :
		for amino_acid in codon_gc_max :
			if c == amino_acid :
				random_choice = uniform(0, 1)
				separateur = 1/len(codon_gc_max[amino_acid])
				n = 1
				while n <= len(codon_gc_max[amino_acid]) :
					if random_choice <= (n*separateur) and random_choice > ((n-1)*separateur) :
						random_select = n-1
					n += 1
				list_codon = list(codon_gc_max[amino_acid])
				choice = list_codon[random_select]
				seq_simulated += choice
	return seq_simulated

################################################################################
## FUNCTION USED TO OPTIMIZE SEQUENCE IN A RANDOM-GUIDED ONE (MAX GC CONTENT) ##
################################################################################

def gc_max_random_guided_simulation(sequence, genetic_code) :

	seq_simulated = ''
	for c in sequence :
		for amino_acid in codon_gc_max :
			separateur = 0
			if c == amino_acid :
				random_choice = uniform(0, 1)
				for codon in codon_gc_max[amino_acid] :
					separateur += codon_gc_max[amino_acid][codon][1]
					if random_choice <= separateur :
						 seq_simulated += codon
						 break
	return seq_simulated

################################################################################
##### FUNCTION USED TO OPTIMIZE SEQUENCE IN A RANDOM ONE (MAX AT CONTENT)  #####
################################################################################

def at_max_random_simulation(sequence, genetic_code) :

	seq_simulated = ''
	for c in sequence :
		for amino_acid in codon_at_max :
			if c == amino_acid :
				random_choice = uniform(0, 1)
				separateur = 1/len(codon_at_max[amino_acid])
				n = 1
				while n <= len(codon_at_max[amino_acid]) :
					if random_choice <= (n*separateur) and random_choice > ((n-1)*separateur) :
						random_select = n-1
					n += 1
				list_codon = list(codon_at_max[amino_acid])
				choice = list_codon[random_select]

				seq_simulated += choice
	return seq_simulated

################################################################################
## FUNCTION USED TO OPTIMIZE SEQUENCE IN A RANDOM-GUIDED ONE (MAX AT CONTENT) ##
################################################################################

def at_max_random_guided_simulation(sequence, genetic_code) :
	seq_simulated = ''
	for c in sequence :
		for amino_acid in codon_at_max :
			separateur = 0
			if c == amino_acid :
				random_choice = uniform(0, 1)
				for codon in codon_at_max[amino_acid] :
					separateur += codon_at_max[amino_acid][codon][1]
					if random_choice <= separateur :
						 seq_simulated += codon
						 break
	return seq_simulated

################################################################################
##### FUNCTION USED TO OPTIMIZE SEQUENCE IN A ONE AA - ONE CODON ONE (MAX) #####
################################################################################

def max_simulation(sequence, genetic_code) :
	best_choice = {}
	for amino_acid in genetic_code :
		best_number = 0
		for codon in genetic_code[amino_acid] :
			if genetic_code[amino_acid][codon][0] > best_number :
				best_number = genetic_code[amino_acid][codon][0]
				best_choice[amino_acid] = codon
			elif genetic_code[amino_acid][codon][0] == best_number and len(genetic_code[amino_acid]) > 1 :
				print ("## IMPORTANT : we have found ties in your CU table for " + amino_acid + " please note that only one of them is taken in account. ##" )
			elif len(genetic_code[amino_acid]) == 1 :
				best_choice[amino_acid] = codon

	seq_simulated = ''
	for c in sequence :
		for amino_acid in best_choice :
			if c == amino_acid :
				seq_simulated += best_choice[amino_acid]
	return seq_simulated

################################################################################
##### FUNCTION USED TO OPTIMIZE SEQUENCE IN A ONE AA - ONE CODON ONE (MIN) #####
################################################################################

def min_simulation(sequence, genetic_code) :
	worst_choice = {}
	for amino_acid in genetic_code :
		worst_number = 1
		for codon in genetic_code[amino_acid] :
			if genetic_code[amino_acid][codon][1] < worst_number :
				worst_number = genetic_code[amino_acid][codon][1]
				worst_choice[amino_acid] = codon
			elif genetic_code[amino_acid][codon][1] == worst_number and len(genetic_code[amino_acid]) > 1 :
				worst_number = genetic_code[amino_acid][codon][1]
				worst_choice[amino_acid] = codon
				print ("## IMPORTANT : we have found ties in your CU table for " + amino_acid + " please note that only one of them is taken in account. ##" )
			elif len(genetic_code[amino_acid]) == 1 :
				worst_choice[amino_acid] = codon

	seq_simulated = ''
	for c in sequence :
		for amino_acid in worst_choice :
			if c == amino_acid :
				seq_simulated += worst_choice[amino_acid]
	return seq_simulated

################################################################################
############## FUNCTION USED TO GENERATE OPTIMIZED SEQUENCES ###################
################################################################################

def optimize(genetic_code, optimal_codons_list, perform_COUSIN_with_boundaries, optimization_type, optimization_number, input_file_seq_type, reverse_frequencies, output_dir) :

	sequential_order_results_optimization.append("firstline")
	results_for_file_optimization["firstline"] = "Header" + "\t" + "Length" + "\t" + "%A(3)" + "\t" + "%T(3)" + "\t" + "%G(3)" + "\t" + "%C(3)" + "\t" + "%GC(all)" + "\t" + "%GC(1)" + "\t" + "%GC(2)" + "\t" + "%GC(3)" + "\t" + "COUSIN_18" + "\t" + "COUSIN_59" + "\t" + "CAI_59" + "\t" + "CAI_18" + "\t" + "ENC" + "\t" + "SCUO" + "\t" + "FOP" + "\t" + "CBI" + "\t" + "ICDI" + "\t" + "Scaled_Chi" + "\t" "GRAVY" + "\t" + "AROMA" + "\t" + "eucl_dist_vs_ref (a.a)" + "\t" + "eucl_dist_vs_H0 (a.a)"
	if not os.path.isdir("./" + output_dir + "/optimization"):
		os.makedirs("./" + output_dir + "/optimization")

	freq_aa_ref = {}

	if reverse_frequencies == 1 : # VERIFIER CETTE HISTOIRE DE REVERSE FREQUENCIES. FAUT IL AVOIR COMME PARAM GENETIC_CODE POUR, PAR EXEMPLE, calc_freq_codon_ref ?
		genetic_code_rev = get_reverse_frequencies(genetic_code)
		get_GC_and_AT_max_codons(genetic_code)
		get_freq_aa(genetic_code, freq_aa_ref)
		calc_freq_equ_muta_bias(genetic_code_rev) # A VERIFIER
		genetic_code_nucl_comp = generate_table_nucl_comp (genetic_code_rev, output_dir)
		deviation_score = calc_deviation_scores(genetic_code)
		deviation_score_nucl_comp = calc_deviation_scores(genetic_code_nucl_comp)
		w_scores = calc_adaptative_values(genetic_code_rev)

	else :
		get_GC_and_AT_max_codons(genetic_code)
		get_freq_aa(genetic_code, freq_aa_ref)
		calc_freq_equ_muta_bias(genetic_code)
		genetic_code_nucl_comp = generate_table_nucl_comp (genetic_code, output_dir)
		deviation_score = calc_deviation_scores(genetic_code)
		deviation_score_nucl_comp = calc_deviation_scores(genetic_code_nucl_comp)
		w_scores = calc_adaptative_values(genetic_code)

	if input_file_seq_type == 'AA' :
		for header_seq in content_seq :
			freq_aa_req = {}
			if optimization_number != '' :
				i = 0
				seq_length = len(content_seq[header_seq]) * 3
				while  i < optimization_number :
					if reverse_frequencies == 0 :

						if optimization_type == 0 :
							seq_simulated = random_simulation(content_seq[header_seq], genetic_code)

						elif optimization_type == 1 :
							seq_simulated = random_guided_simulation(content_seq[header_seq], genetic_code)

						elif optimization_type == 2 :
							seq_simulated = gc_max_random_simulation(content_seq[header_seq], genetic_code)

						elif optimization_type == 3 :
							seq_simulated = gc_max_random_guided_simulation(content_seq[header_seq], genetic_code)

						elif optimization_type == 4 :
							seq_simulated = at_max_random_simulation(content_seq[header_seq], genetic_code)

						elif optimization_type == 5 :
							seq_simulated = at_max_random_guided_simulation(content_seq[header_seq], genetic_code)

						elif optimization_type == 6 :
							seq_simulated = max_simulation(content_seq[header_seq], genetic_code)

						elif optimization_type == 7 :
							seq_simulated = min_simulation(content_seq[header_seq], genetic_code)

						ATGC_percent = get_freq_ATGC_seq(seq_simulated)
						GC_percent = get_freq_GC_seq(seq_simulated)
						dict_seq_simulated = seq_to_dict(seq_simulated, header_seq + "_" + str(i), genetic_code, output_dir) # ATTENTION, LA FORMULE SEQ_TO_DICT POURRAIT ESSAYER DE SUPPRIMER LE DICO DE LA SEQUENCE SIMULEE (CA SERAIT BIZARRE MAIS ATTENTION)

						if dict_seq_simulated != 'deleted' :

							aa_in_req = check_if_aa_in_req(dict_seq_simulated)
							deviation_score_aa = calc_deviation_scores_aa(genetic_code, aa_in_req)

							get_freq_aa(dict_seq_simulated, freq_aa_req)
							RSCU_data = calc_RSCU(dict_seq_simulated)

							CAI_cod = calc_CAI_59(dict_seq_simulated, genetic_code, w_scores)
							CAI_aa = calc_CAI_18(dict_seq_simulated, genetic_code, w_scores)
							COUSIN = calc_COUSIN_18(dict_seq_simulated, genetic_code, deviation_score, aa_in_req, perform_COUSIN_with_boundaries)
							COUSIN_pond = calc_COUSIN_59(dict_seq_simulated, genetic_code, deviation_score, deviation_score_aa, aa_in_req, perform_COUSIN_with_boundaries)
	
							ENC = calc_ENC(dict_seq_simulated)
							SCUO = calc_SCUO(dict_seq_simulated, genetic_code)
							FOP = calc_FOP(dict_seq_simulated, optimal_codons_list)
							CBI = calc_CBI(dict_seq_simulated, optimal_codons_list)
							ICDI = calc_ICDI(dict_seq_simulated, RSCU_data, aa_in_req)
							scaled_chi = calc_scaled_chi(dict_seq_simulated)

							GRAVY = calc_GRAVY(dict_seq_simulated)
							AROMA = calc_AROMA(dict_seq_simulated)
							eucl_dist_aa_que_vs_ref = calc_eucl_dist_que_vs_ref(freq_aa_req, freq_aa_ref, aa_in_req)
							eucl_dist_aa_que_vs_H0 = calc_eucl_dist_que_vs_H0(freq_aa_req, aa_in_req)

							if (header_seq + "_" + str(i)) not in results_for_file_optimization :
								sequential_order_results_optimization.append(header_seq + "_" + str(i))
								results_for_file_optimization[header_seq + "_" + str(i)] = header_seq + "_" + str(i) + "\t" + str(seq_length) + "\t" + str(round(ATGC_percent[0],3)) + "\t" + str(round(ATGC_percent[1],3)) + "\t" + str(round(ATGC_percent[2],3)) + "\t" + str(round(ATGC_percent[3],3)) + "\t" + str(round(GC_percent[0],3)) + "\t" + str(round(GC_percent[1],3)) + "\t" + str(round(GC_percent[2],3)) + "\t" + str(round(GC_percent[3],3)) + "\t" + str(round(COUSIN,3))  + "\t" + str(round(COUSIN_pond,3)) + "\t" + str(round(CAI_cod,3))  + "\t" + str(round(CAI_aa,3)) + "\t" + str(round(ENC,3)) + "\t" + str(round(SCUO,3)) +  "\t" + str(round(FOP,3)) +  "\t" + str(round(CBI,3)) +  "\t" + str(round(ICDI,3)) +   "\t" + str(round(scaled_chi,3)) + "\t" + str(round(GRAVY,3)) + "\t" + str(round(AROMA,3)) + "\t" + str(round(eucl_dist_aa_que_vs_ref,3)) + "\t" + str(round(eucl_dist_aa_que_vs_H0,3))
							# if (header_seq + "_" + str(i)) not in results_for_file_optimization :
							#     sequential_order_results_optimization.append(header_seq + "_" + str(i))
							#     results_for_file_optimization[header_seq + "_" + str(i)] = header_seq + "_" + str(i) + "\t" + str(seq_length) + "\t" + str(round(ATGC_percent[0],3)) + "\t" + str(round(ATGC_percent[1],3)) + "\t" + str(round(ATGC_percent[2],3)) + "\t" + str(round(ATGC_percent[3],3)) + "\t" + str(round(GC_percent[0],3)) + "\t" + str(round(GC_percent[1],3)) + "\t" + str(round(GC_percent[2],3)) + "\t" + str(round(GC_percent[3],3)) + "\t" + str(round(COUSIN,3))  + "\t" + str(round(COUSIN_pond,3)) + "\t" + str(round(COUSIN_muta_bias,3)) + "\t" + str(round(COUSIN_muta_bias_pond,3))  + "\t" + str(round(CAI_cod,3))  + "\t" + str(round(CAI_aa,3)) + "\t" + str(round(ENC,3)) + "\t" + str(round(SCUO,3)) +  "\t" + str(round(FOP,3)) +  "\t" + str(round(CBI,3)) +  "\t" + str(round(ICDI,3)) +   "\t" + str(round(scaled_chi,3)) + "\t" + str(round(GRAVY,3)) + "\t" + str(round(AROMA,3)) + "\t" + str(round(eucl_dist_aa_que_vs_ref,3))


							optimized_sequences = open("./" + output_dir + "/optimization/optimized_sequences.fasta", "a")
							optimized_sequences.write(">" + header_seq + "_" + str(i) + "\n" + RNA_to_DNA(seq_simulated) + "\n\n")
							optimized_sequences.close()

					elif reverse_frequencies == 1 :

						if optimization_type == 0 :
							seq_simulated = random_simulation(content_seq[header_seq], genetic_code_rev)

						elif optimization_type == 1 :
							seq_simulated = random_guided_simulation(content_seq[header_seq], genetic_code_rev)

						elif optimization_type == 2 :
							seq_simulated = gc_max_random_simulation(content_seq[header_seq], genetic_code_rev)

						elif optimization_type == 3 :
							seq_simulated = gc_max_random_guided_simulation(content_seq[header_seq], genetic_code_rev)

						elif optimization_type == 4 :
							seq_simulated = at_max_random_simulation(content_seq[header_seq], genetic_code_rev)

						elif optimization_type == 5 :
							seq_simulated = at_max_random_guided_simulation(content_seq[header_seq], genetic_code_rev)

						elif optimization_type == 6 :
							seq_simulated = max_simulation(content_seq[header_seq], genetic_code_rev)

						elif optimization_type == 7 :
							seq_simulated = min_simulation(content_seq[header_seq], genetic_code_rev)

						ATGC_percent = get_freq_ATGC_seq(seq_simulated)
						GC_percent = get_freq_GC_seq(seq_simulated)
						dict_seq_simulated = seq_to_dict(seq_simulated, header_seq + "_" + str(i), genetic_code, output_dir)

						if dict_seq_simulated != 'deleted' :

							aa_in_req = check_if_aa_in_req(dict_seq_simulated)
							deviation_score_aa = calc_deviation_scores_aa(genetic_code, aa_in_req)

							get_freq_aa(dict_seq_simulated, freq_aa_req)
							RSCU_data = calc_RSCU(dict_seq_simulated)

							CAI_cod = calc_CAI_59(dict_seq_simulated, genetic_code, w_scores)
							CAI_aa = calc_CAI_18(dict_seq_simulated, genetic_code, w_scores)
							COUSIN = calc_COUSIN_18(dict_seq_simulated, genetic_code, deviation_score, aa_in_req, perform_COUSIN_with_boundaries)
							COUSIN_pond = calc_COUSIN_59(dict_seq_simulated, genetic_code, deviation_score, deviation_score_aa, aa_in_req, perform_COUSIN_with_boundaries)
							
							ENC = calc_ENC(dict_seq_simulated)
							SCUO = calc_SCUO(dict_seq_simulated, genetic_code)
							FOP = calc_FOP(dict_seq_simulated, optimal_codons_list)
							CBI = calc_CBI(dict_seq_simulated, optimal_codons_list)
							ICDI = calc_ICDI(dict_seq_simulated, RSCU_data, aa_in_req)
							scaled_chi = calc_scaled_chi(dict_seq_simulated)

							GRAVY = calc_GRAVY(dict_seq_simulated)
							AROMA = calc_AROMA(dict_seq_simulated)
							eucl_dist_aa_que_vs_ref = calc_eucl_dist_que_vs_ref(freq_aa_req, freq_aa_ref, aa_in_req)
							eucl_dist_aa_que_vs_H0 = calc_eucl_dist_que_vs_H0(freq_aa_req, aa_in_req)


							if (header_seq + "_" + str(i)) not in results_for_file_optimization :
								sequential_order_results_optimization.append(header_seq + "_" + str(i))
								results_for_file_optimization[header_seq + "_" + str(i)] = header_seq + "_" + str(i) + "\t" + str(seq_length) + "\t" + str(round(ATGC_percent[0],3)) + "\t" + str(round(ATGC_percent[1],3)) + "\t" + str(round(ATGC_percent[2],3)) + "\t" + str(round(ATGC_percent[3],3)) + "\t" + str(round(GC_percent[0],3)) + "\t" + str(round(GC_percent[1],3)) + "\t" + str(round(GC_percent[2],3)) + "\t" + str(round(GC_percent[3],3)) + "\t" + str(round(COUSIN,3))  + "\t" + str(round(COUSIN_pond,3)) + "\t" + str(round(CAI_cod,3))  + "\t" + str(round(CAI_aa,3)) + "\t" + str(round(ENC,3)) + "\t" + str(round(SCUO,3)) +  "\t" + str(round(FOP,3)) +  "\t" + str(round(CBI,3)) +  "\t" + str(round(ICDI,3)) + "\t" + str(round(scaled_chi,3)) + "\t" + str(round(GRAVY,3)) + "\t" + str(round(AROMA,3)) + "\t" + str(round(eucl_dist_aa_que_vs_ref,3)) + "\t" + str(round(eucl_dist_aa_que_vs_H0,3))
					
							optimized_sequences = open("./" + output_dir + "/optimization/optimized_sequences.fasta", "a")
							optimized_sequences.write(">" + header_seq + "_" + str(i) + "\n" + RNA_to_DNA(seq_simulated) + "\n\n")
							optimized_sequences.close()
					i += 1

	elif input_file_seq_type == "DNA" :
		for header_seq in content_seq :
			freq_aa_req = {}
			seq_length = len(content_seq[header_seq])
			seq_DNA_to_RNA = DNA_to_RNA(content_seq[header_seq])
			seq_RNA_to_AA = RNA_to_AA(seq_DNA_to_RNA, genetic_code)

			ATGC_percent = get_freq_ATGC_seq(content_seq[header_seq])
			GC_percent = get_freq_GC_seq(content_seq[header_seq])

			if header_seq not in results_for_file_optimization :
				sequential_order_results_optimization.append(header_seq)
				results_for_file_optimization[header_seq] = header_seq + "\t" + str(len(content_seq[header_seq])) + "\t" + str(round(ATGC_percent[0],3)) + "\t" + str(round(ATGC_percent[1],3)) + "\t" + str(round(ATGC_percent[2],3)) + "\t" + str(round(ATGC_percent[3],3)) + "\t" + str(round(GC_percent[0],3)) + "\t" + str(round(GC_percent[1],3)) + "\t" + str(round(GC_percent[2],3)) + "\t" + str(round(GC_percent[3],3)) + "\t" + str(round(input_seq_value_COUSIN[header_seq],3))  + "\t" + str(round(input_seq_value_COUSIN_pond[header_seq],3))  + "\t" + str(round(input_seq_value_CAI_cod[header_seq],3))  + "\t" + str(round(input_seq_value_CAI_aa[header_seq],3)) + "\t" + str(round(input_seq_value_ENC[header_seq],3)) + "\t" + str(round(input_seq_value_SCUO[header_seq],3)) + "\t" + str(round(input_seq_value_FOP[header_seq],3)) + "\t" + str(round(input_seq_value_CBI[header_seq],3)) + "\t" + str(round(input_seq_value_ICDI[header_seq],3)) + "\t" + str(round(input_seq_value_scaled_chi[header_seq],3)) + "\t" + str(round(input_seq_value_GRAVY[header_seq],3)) + "\t" + str(round(input_seq_value_AROMA[header_seq],3)) + "\t" + str(round(input_seq_value_eucl_dist_aa_que_vs_ref[header_seq],3)) + "\t" + str(round(input_seq_value_eucl_dist_aa_que_vs_H0[header_seq],3))

			optimized_sequences = open("./" + output_dir + "/optimization/optimized_sequences.fasta", "a")
			optimized_sequences.write(">" + header_seq + "\n" + content_seq[header_seq] + "\n\n")
			optimized_sequences.close()

			if optimization_number != '' :
				i = 0
				seq_length = len(content_seq[header_seq])
				while  i < optimization_number :
					if reverse_frequencies == 0 :

						if optimization_type == 0 :
							seq_simulated = random_simulation(seq_RNA_to_AA, genetic_code)

						elif optimization_type == 1 :
							seq_simulated = random_guided_simulation(seq_RNA_to_AA, genetic_code)

						elif optimization_type == 2 :
							seq_simulated = gc_max_random_simulation(seq_RNA_to_AA, genetic_code)

						elif optimization_type == 3 :
							seq_simulated = gc_max_random_guided_simulation(seq_RNA_to_AA, genetic_code)

						elif optimization_type == 4 :
							seq_simulated = at_max_random_simulation(seq_RNA_to_AA, genetic_code)

						elif optimization_type == 5 :
							seq_simulated = at_max_random_guided_simulation(seq_RNA_to_AA, genetic_code)

						elif optimization_type == 6 :
							seq_simulated = max_simulation(seq_RNA_to_AA, genetic_code)

						elif optimization_type == 7 :
							seq_simulated = min_simulation(seq_RNA_to_AA, genetic_code)

						ATGC_percent = get_freq_ATGC_seq(seq_simulated)
						GC_percent = get_freq_GC_seq(seq_simulated)
						dict_seq_simulated = seq_to_dict(seq_simulated, header_seq + "_" + str(i), genetic_code, output_dir)
						if dict_seq_simulated != 'deleted' :
							aa_in_req = check_if_aa_in_req(dict_seq_simulated)

							deviation_score_aa = calc_deviation_scores_aa(genetic_code, aa_in_req)
							get_freq_aa(dict_seq_simulated, freq_aa_req)
							RSCU_data = calc_RSCU(dict_seq_simulated)

							CAI_cod = calc_CAI_59(dict_seq_simulated, genetic_code, w_scores)
							CAI_aa = calc_CAI_18(dict_seq_simulated, genetic_code, w_scores)
							COUSIN = calc_COUSIN_18(dict_seq_simulated, genetic_code, deviation_score, aa_in_req, perform_COUSIN_with_boundaries)
							COUSIN_pond = calc_COUSIN_59(dict_seq_simulated, genetic_code, deviation_score, deviation_score_aa, aa_in_req, perform_COUSIN_with_boundaries)
							#COUSIN_muta_bias = calc_COUSIN_18(dict_seq_simulated, genetic_code_nucl_comp, deviation_score_nucl_comp, aa_in_req)
							#COUSIN_muta_bias_pond = calc_COUSIN_59(dict_seq_simulated, genetic_code_nucl_comp, deviation_score_nucl_comp, aa_in_req)

							ENC = calc_ENC(dict_seq_simulated)
							SCUO = calc_SCUO(dict_seq_simulated, genetic_code)
							FOP = calc_FOP(dict_seq_simulated, optimal_codons_list)
							CBI = calc_CBI(dict_seq_simulated, optimal_codons_list)
							ICDI = calc_ICDI(dict_seq_simulated, RSCU_data, aa_in_req)
							scaled_chi = calc_scaled_chi(dict_seq_simulated)

							GRAVY = calc_GRAVY(dict_seq_simulated)
							AROMA = calc_AROMA(dict_seq_simulated)
							eucl_dist_aa_que_vs_ref = calc_eucl_dist_que_vs_ref(freq_aa_req, freq_aa_ref, aa_in_req)
							eucl_dist_aa_que_vs_H0 = calc_eucl_dist_que_vs_H0(freq_aa_req, aa_in_req)

							if (header_seq + "_" + str(i)) not in results_for_file_optimization :
								sequential_order_results_optimization.append(header_seq + "_" + str(i))
								results_for_file_optimization[header_seq + "_" + str(i)] = header_seq + "_" + str(i) + "\t" + str(seq_length) + "\t" + str(round(ATGC_percent[0],3)) + "\t" + str(round(ATGC_percent[1],3)) + "\t" + str(round(ATGC_percent[2],3)) + "\t" + str(round(ATGC_percent[3],3)) + "\t" + str(round(GC_percent[0],3)) + "\t" + str(round(GC_percent[1],3)) + "\t" + str(round(GC_percent[2],3)) + "\t" + str(round(GC_percent[3],3)) + "\t" + str(round(COUSIN,3))  + "\t" + str(round(COUSIN_pond,3)) + "\t" + str(round(CAI_cod,3))  + "\t" + str(round(CAI_aa,3)) + "\t" + str(round(ENC,3)) + "\t" + str(round(SCUO,3)) +  "\t" + str(round(FOP,3)) +  "\t" + str(round(CBI,3)) +  "\t" + str(round(ICDI,3)) + "\t" + str(round(scaled_chi,3)) + "\t" + str(round(GRAVY,3)) + "\t" + str(round(AROMA,3)) + "\t" + str(round(eucl_dist_aa_que_vs_ref,3)) + "\t" + str(round(eucl_dist_aa_que_vs_H0,3))
							# if (header_seq + "_" + str(i)) not in results_for_file_optimization :
							#     sequential_order_results_optimization.append(header_seq + "_" + str(i))
							#     results_for_file_optimization[header_seq + "_" + str(i)] = header_seq + "_" + str(i) + "\t" + str(seq_length) + "\t" + str(round(ATGC_percent[0],3)) + "\t" + str(round(ATGC_percent[1],3)) + "\t" + str(round(ATGC_percent[2],3)) + "\t" + str(round(ATGC_percent[3],3)) + "\t" + str(round(GC_percent[0],3)) + "\t" + str(round(GC_percent[1],3)) + "\t" + str(round(GC_percent[2],3)) + "\t" + str(round(GC_percent[3],3)) + "\t" + str(round(COUSIN,3))  + "\t" + str(round(COUSIN_pond,3)) + "\t" + str(round(COUSIN_muta_bias,3)) + "\t" + str(round(COUSIN_muta_bias_pond,3))  + "\t" + str(round(CAI_cod,3))  + "\t" + str(round(CAI_aa,3)) + "\t" + str(round(ENC,3)) + "\t" + str(round(SCUO,3)) +  "\t" + str(round(FOP,3)) +  "\t" + str(round(CBI,3)) +  "\t" + str(round(ICDI,3)) + "\t" + str(round(scaled_chi,3)) + "\t" + str(round(GRAVY,3)) + "\t" + str(round(AROMA,3)) + "\t" + str(round(eucl_dist_aa_que_vs_ref,3)) + "\t" + str(round(eucl_dist_aa_que_vs_H0,3))

							optimized_sequences = open("./" + output_dir + "/optimization/optimized_sequences.fasta", "a")
							optimized_sequences.write(">" + header_seq + "_" + str(i) + "\n" + RNA_to_DNA(seq_simulated) + "\n\n")
							optimized_sequences.close()

					elif reverse_frequencies == 1 :

						if optimization_type == 0 :
							seq_simulated = random_simulation(seq_RNA_to_AA, genetic_code_rev)

						elif optimization_type == 1 :
							seq_simulated = random_guided_simulation(seq_RNA_to_AA, genetic_code_rev)

						elif optimization_type == 2 :
							seq_simulated = gc_max_random_simulation(seq_RNA_to_AA, genetic_code_rev)

						elif optimization_type == 3 :
							seq_simulated = gc_max_random_guided_simulation(seq_RNA_to_AA, genetic_code_rev)

						elif optimization_type == 4 :
							seq_simulated = at_max_random_simulation(seq_RNA_to_AA, genetic_code_rev)

						elif optimization_type == 5 :
							seq_simulated = at_max_random_guided_simulation(seq_RNA_to_AA, genetic_code_rev)

						elif optimization_type == 6 :
							seq_simulated = max_simulation(seq_RNA_to_AA, genetic_code_rev)

						elif optimization_type == 7 :
							seq_simulated = min_simulation(seq_RNA_to_AA, genetic_code_rev)

						ATGC_percent = get_freq_ATGC_seq(seq_simulated)
						GC_percent = get_freq_GC_seq(seq_simulated)
						dict_seq_simulated = seq_to_dict(seq_simulated, header_seq + "_" + str(i), genetic_code, output_dir)

						if dict_seq_simulated != 'deleted ' :

							aa_in_req = check_if_aa_in_req(dict_seq_simulated)
							deviation_score_aa = calc_deviation_scores_aa(genetic_code, aa_in_req)
							get_freq_aa(dict_seq_simulated, freq_aa_req)
							RSCU_data = calc_RSCU(dict_seq_simulated)

							CAI_cod = calc_CAI_59(dict_seq_simulated, genetic_code, w_scores)
							CAI_aa = calc_CAI_18(dict_seq_simulated, genetic_code, w_scores)
							COUSIN = calc_COUSIN_18(dict_seq_simulated, genetic_code, deviation_score, aa_in_req, perform_COUSIN_with_boundaries)
							COUSIN_pond = calc_COUSIN_59(dict_seq_simulated, genetic_code, deviation_score, deviation_score_aa, aa_in_req, perform_COUSIN_with_boundaries)
							#COUSIN_muta_bias = calc_COUSIN_18(dict_seq_simulated, genetic_code_rev_nucl_comp, deviation_score_nucl_comp, aa_in_req)
							#COUSIN_muta_bias_pond = calc_COUSIN_59(dict_seq_simulated, genetic_code_rev_nucl_comp, deviation_score_nucl_comp, aa_in_req)
							ENC = calc_ENC(dict_seq_simulated)

							SCUO = calc_SCUO(dict_seq_simulated, genetic_code)
							FOP = calc_FOP(dict_seq_simulated, optimal_codons_list)
							CBI = calc_CBI(dict_seq_simulated, optimal_codons_list)
							ICDI = calc_ICDI(dict_seq_simulated, RSCU_data, aa_in_req)
							scaled_chi = calc_scaled_chi(dict_seq_simulated)

							GRAVY = calc_GRAVY(dict_seq_simulated)
							AROMA = calc_AROMA(dict_seq_simulated)
							eucl_dist_aa_que_vs_ref = calc_eucl_dist_que_vs_ref(freq_aa_req, freq_aa_ref, aa_in_req)
							eucl_dist_aa_que_vs_H0 = calc_eucl_dist_que_vs_H0(freq_aa_req, aa_in_req)

							if (header_seq + "_" + str(i)) not in results_for_file_optimization :
								sequential_order_results_optimization.append(header_seq + "_" + str(i))
								results_for_file_optimization[header_seq + "_" + str(i)] = header_seq + "_" + str(i) + "\t" + str(seq_length) + "\t" + str(round(ATGC_percent[0],3)) + "\t" + str(round(ATGC_percent[1],3)) + "\t" + str(round(ATGC_percent[2],3)) + "\t" + str(round(ATGC_percent[3],3)) + "\t" + str(round(GC_percent[0],3)) + "\t" + str(round(GC_percent[1],3)) + "\t" + str(round(GC_percent[2],3)) + "\t" + str(round(GC_percent[3],3)) + "\t" + str(round(COUSIN,3))  + "\t" + str(round(COUSIN_pond,3)) + "\t" + str(round(CAI_cod,3))  + "\t" + str(round(CAI_aa,3)) + "\t" + str(round(ENC,3)) + "\t" + str(round(SCUO,3)) +  "\t" + str(round(FOP,3)) +  "\t" + str(round(CBI,3)) +  "\t" + str(round(ICDI,3)) + "\t" + str(round(scaled_chi,3)) + "\t" + str(round(GRAVY,3)) + "\t" + str(round(AROMA,3)) + "\t" + str(round(eucl_dist_aa_que_vs_ref,3)) + "\t" + str(round(eucl_dist_aa_que_vs_H0,3))
							# if (header_seq + "_" + str(i)) not in results_for_file_optimization :
							#     sequential_order_results_optimization.append(header_seq + "_" + str(i))
							#     results_for_file_optimization[header_seq + "_" + str(i)] = header_seq + "_" + str(i) + "\t" + str(seq_length) + "\t" + str(round(ATGC_percent[0],3)) + "\t" + str(round(ATGC_percent[1],3)) + "\t" + str(round(ATGC_percent[2],3)) + "\t" + str(round(ATGC_percent[3],3)) + "\t" + str(round(GC_percent[0],3)) + "\t" + str(round(GC_percent[1],3)) + "\t" + str(round(GC_percent[2],3)) + "\t" + str(round(GC_percent[3],3)) + "\t" + str(round(COUSIN,3))  + "\t" + str(round(COUSIN_pond,3)) + "\t" + str(round(COUSIN_muta_bias,3)) + "\t" + str(round(COUSIN_muta_bias_pond,3))  + "\t" + str(round(CAI_cod,3))  + "\t" + str(round(CAI_aa,3)) + "\t" + str(round(ENC,3)) + "\t" + str(round(SCUO,3)) +  "\t" + str(round(FOP,3)) +  "\t" + str(round(CBI,3)) +  "\t" + str(round(ICDI,3)) + "\t" + str(round(scaled_chi,3)) + "\t" + str(round(GRAVY,3)) + "\t" + str(round(AROMA,3)) + "\t" + str(round(eucl_dist_aa_que_vs_ref,3)) + "\t" + str(round(eucl_dist_aa_que_vs_H0,3))
							optimized_sequences = open("./" + output_dir + "/optimization/optimized_sequences.fasta", "a")
							optimized_sequences.write(">" + header_seq + "_" + str(i) + "\n" + RNA_to_DNA(seq_simulated) + "\n\n")
							optimized_sequences.close()
					i += 1

################################################################################
################################################################################
############################# CAI CALCULATION ##################################
################################################################################
################################################################################

################################################################################
############################# CALCULATION WITH CODONS ##########################
################################################################################

def calc_adaptative_values(genetic_code) :

	w_scores = {}
	best_freq_codon_ref = {}

	for amino_acid in genetic_code :
		best_freq_value = 0
		for codon in genetic_code[amino_acid] :
			if genetic_code[amino_acid][codon][1] > best_freq_value :
				best_freq_codon_ref[amino_acid] = genetic_code[amino_acid][codon][1]
				best_freq_value = genetic_code[amino_acid][codon][1]

	for amino_acid in genetic_code :
		if len(genetic_code[amino_acid]) > 1 and amino_acid != '*' :
			for codon in genetic_code[amino_acid] :
				w_scores[codon] = (genetic_code[amino_acid][codon][1] / best_freq_codon_ref[amino_acid])

	return w_scores

def calc_CAI_59(genetic_code, genetic_code_ref, w_scores) :
	CAI = 0
	length_seq_CAI = 0

	for amino_acid in genetic_code :
		if len(genetic_code[amino_acid]) > 1 and amino_acid != '*' :
			for codon in genetic_code[amino_acid] :
				length_seq_CAI += genetic_code[amino_acid][codon][0]
	for amino_acid in genetic_code :
		if len(genetic_code[amino_acid]) > 1 and amino_acid != '*' :
			for codon in genetic_code[amino_acid] :
				if genetic_code[amino_acid][codon][0] > 0 :
					if codon in w_scores :
						CAI += (math.log(w_scores[codon]) * genetic_code[amino_acid][codon][0])
					else :
						CAI += (math.log(0.5) * genetic_code[amino_acid][codon][0])
	if length_seq_CAI > 0 :
		CAI = math.exp(CAI * (1/(length_seq_CAI)))
	else :
		print ("## IMPORTANT (CAI index): error, your sequence seems empty or isn't compounded of amino acids associated with synonymous codons. CAI score = 0.\n ##")
		CAI = 0

	return CAI

################################################################################
################### CALCULATION WITH AMINO ACIDS FAMILIES ######################
################################################################################


def calc_CAI_18(genetic_code, genetic_code_ref, w_scores) :
	CAI = 0
	length_seq_amino_acid = {}
	sum_cod_aa_family = {}
	CAI_by_amino_acid = {}

	for amino_acid in genetic_code :
		if len(genetic_code[amino_acid]) > 1 and amino_acid != '*' :
			for codon in genetic_code[amino_acid] :
				if genetic_code[amino_acid][codon][0] > 0 :
					if amino_acid not in length_seq_amino_acid :
						length_seq_amino_acid[amino_acid] = genetic_code[amino_acid][codon][0]
					else :
						length_seq_amino_acid[amino_acid] += genetic_code[amino_acid][codon][0]

			for codon in genetic_code[amino_acid] :
				if genetic_code[amino_acid][codon][0] > 0 :
					if codon in w_scores :
						if amino_acid not in CAI_by_amino_acid :
							CAI_by_amino_acid[amino_acid] = (math.log(w_scores[codon]) * genetic_code[amino_acid][codon][0])
						else :
							CAI_by_amino_acid[amino_acid] += (math.log(w_scores[codon]) * genetic_code[amino_acid][codon][0])
					else :
						if amino_acid not in CAI_by_amino_acid :
							CAI_by_amino_acid[amino_acid] = (math.log(0.5) * genetic_code[amino_acid][codon][0])
						else :
							CAI_by_amino_acid[amino_acid] += (math.log(0.5) * genetic_code[amino_acid][codon][0])

	for amino_acid in genetic_code :
		if amino_acid in CAI_by_amino_acid and amino_acid in length_seq_amino_acid :
			CAI += math.log(math.exp(CAI_by_amino_acid[amino_acid] * (1/(length_seq_amino_acid[amino_acid]))))

	if len(CAI_by_amino_acid) > 0 :
		CAI = math.exp(CAI * (1/len(CAI_by_amino_acid)))
	else :
		print ("## IMPORTANT (CAI index): error, your sequence seems empty or isn't compounded of amino acids associated with synonymous codons. CAI score = 0.\n ##")
		CAI = 0

	return CAI

################################################################################
################################################################################
########################## COUSIN CALCULATION ##################################
################################################################################
################################################################################

################################################################################
############################# CALCULATION WITH CODONS ##########################
################################################################################

def calc_COUSIN_59(genetic_code, genetic_code_ref, deviation_score, deviation_score_aa, is_aa_in_req, perform_COUSIN_with_boundaries) : # A VERIFIER A 100%

	COUSIN = 0
	cpt_nb_amino_acid_req = 0

	weight_req = {}
	weight_ref = {}
	sum_COUSIN_aa_by_aa = 0
	occ_tot_in_req = 0
	occ_tot_in_ref = 0
	occ_aa_in_req = {}
	occ_aa_in_ref = {}
	freq_aa_ref = {}
	freq_aa_req = {}
	weight_aa_ref = {}
	weight_aa_req = {}

	COUSIN = 0

	if perform_COUSIN_with_boundaries == "Yes" or perform_COUSIN_with_boundaries == "YES" or perform_COUSIN_with_boundaries == "yes" or perform_COUSIN_with_boundaries == "Y" or perform_COUSIN_with_boundaries == "y"  :

		for amino_acid in genetic_code_ref :
			if len(genetic_code_ref[amino_acid]) > 1 and amino_acid != '*' :
					if (is_aa_in_req[amino_acid] == True) :

						for codon in genetic_code_ref[amino_acid] :

							occ_tot_in_req += genetic_code[amino_acid][codon][0]
							occ_tot_in_ref += genetic_code_ref[amino_acid][codon][0]

							if amino_acid not in weight_req and amino_acid not in occ_aa_in_req :
								occ_aa_in_req[amino_acid] = genetic_code[amino_acid][codon][0]
								weight_req[amino_acid] = (genetic_code[amino_acid][codon][1] * deviation_score[codon])
							else :
								occ_aa_in_req[amino_acid] += genetic_code[amino_acid][codon][0]
								weight_req[amino_acid] += (genetic_code[amino_acid][codon][1] * deviation_score[codon])

							if amino_acid not in weight_ref and amino_acid not in occ_aa_in_ref :
								occ_aa_in_ref[amino_acid] = genetic_code_ref[amino_acid][codon][0]
								weight_ref[amino_acid] = (genetic_code_ref[amino_acid][codon][1] * deviation_score[codon])
							else :
								occ_aa_in_ref[amino_acid] += genetic_code_ref[amino_acid][codon][0]
								weight_ref[amino_acid] += (genetic_code_ref[amino_acid][codon][1] * deviation_score[codon])

		total_freq = 0
		for amino_acid in genetic_code_ref :
			if len(genetic_code_ref[amino_acid]) > 1 and amino_acid != '*' and amino_acid in weight_req and amino_acid in weight_ref :
					if (is_aa_in_req[amino_acid] == True) :
						weight_by_amino_acid_req = (weight_req[amino_acid]) 
						weight_by_amino_acid_ref = (weight_ref[amino_acid])
						if weight_ref[amino_acid] != 0 :

							if (weight_by_amino_acid_req / weight_by_amino_acid_ref) > 4 :
								sum_COUSIN_aa_by_aa += (occ_aa_in_req[amino_acid] / occ_tot_in_req) *  (4)
							elif (weight_by_amino_acid_req / weight_by_amino_acid_ref) < - 3 :
								sum_COUSIN_aa_by_aa += ((occ_aa_in_req[amino_acid] / occ_tot_in_req)) *  (- 3)
							else :
								sum_COUSIN_aa_by_aa += round(( (occ_aa_in_req[amino_acid] / occ_tot_in_req) * (weight_by_amino_acid_req / weight_by_amino_acid_ref)),3)
						else :
							sum_COUSIN_aa_by_aa += 0
							weight_ref.pop(amino_acid)
							weight_req.pop(amino_acid)

		if weight_req != 0 and weight_ref != 0 and len(weight_req) == len(weight_ref) :
			COUSIN = (sum_COUSIN_aa_by_aa)
		else :
			COUSIN = 0
		return COUSIN

	else :

		for amino_acid in genetic_code_ref :
			if len(genetic_code_ref[amino_acid]) > 1 and amino_acid != '*' :
					if (is_aa_in_req[amino_acid] == True) :

						for codon in genetic_code_ref[amino_acid] :

							occ_tot_in_req += genetic_code[amino_acid][codon][0]
							occ_tot_in_ref += genetic_code_ref[amino_acid][codon][0]

							if amino_acid not in weight_req and amino_acid not in occ_aa_in_req :
								occ_aa_in_req[amino_acid] = genetic_code[amino_acid][codon][0]
								weight_req[amino_acid] = (genetic_code[amino_acid][codon][1] * deviation_score[codon])
							else :
								occ_aa_in_req[amino_acid] += genetic_code[amino_acid][codon][0]
								weight_req[amino_acid] += (genetic_code[amino_acid][codon][1] * deviation_score[codon])

							if amino_acid not in weight_ref and amino_acid not in occ_aa_in_ref :
								occ_aa_in_ref[amino_acid] = genetic_code_ref[amino_acid][codon][0]
								weight_ref[amino_acid] = (genetic_code_ref[amino_acid][codon][1] * deviation_score[codon])
							else :
								occ_aa_in_ref[amino_acid] += genetic_code_ref[amino_acid][codon][0]
								weight_ref[amino_acid] += (genetic_code_ref[amino_acid][codon][1] * deviation_score[codon])

		total_freq = 0
		for amino_acid in genetic_code_ref :
			if len(genetic_code_ref[amino_acid]) > 1 and amino_acid != '*' and amino_acid in weight_req and amino_acid in weight_ref :
					if (is_aa_in_req[amino_acid] == True) :
						weight_by_amino_acid_req = (weight_req[amino_acid]) 
						weight_by_amino_acid_ref = (weight_ref[amino_acid])
						if weight_ref[amino_acid] != 0 :
							sum_COUSIN_aa_by_aa += round(( (occ_aa_in_req[amino_acid] / occ_tot_in_req) * (weight_by_amino_acid_req / weight_by_amino_acid_ref)),3)
						else :
							sum_COUSIN_aa_by_aa += 0
							weight_ref.pop(amino_acid)
							weight_req.pop(amino_acid)

		if weight_req != 0 and weight_ref != 0 and len(weight_req) == len(weight_ref) :
			COUSIN = (sum_COUSIN_aa_by_aa)
		else :
			COUSIN = 0
		return COUSIN

################################################################################
#################### CALCULATION WITH AMINO ACIDS FAMILIES #####################
################################################################################


def calc_deviation_scores (genetic_code_ref) :

	deviation_score = {}

	for amino_acid in genetic_code_ref :
		if amino_acid not in freq_equ :
			freq_equ[amino_acid] = (1 / (len(genetic_code_ref[amino_acid])))

	for amino_acid in genetic_code_ref :
		if len(genetic_code_ref[amino_acid]) > 1 and amino_acid != '*' :
			for codon in genetic_code_ref[amino_acid] :
				if codon not in deviation_score :
					deviation_score[codon] = round((genetic_code_ref[amino_acid][codon][1] - freq_equ[amino_acid]),3)

	return deviation_score

def calc_deviation_scores_aa (genetic_code_ref, is_aa_in_req) :

	deviation_score_aa = {}
	occ_aa = {}
	occ_tot_aa = 0
	number_of_aa = 0

	for amino_acid in genetic_code_ref : 
		if amino_acid in is_aa_in_req and len(genetic_code_ref[amino_acid]) > 1 and amino_acid != '*' : 
			number_of_aa += 1
			for codon in genetic_code_ref[amino_acid] : 
				if amino_acid not in occ_aa : 
					occ_aa[amino_acid] = genetic_code_ref[amino_acid][codon][0]
				else : 
					occ_aa[amino_acid] += genetic_code_ref[amino_acid][codon][0]
				occ_tot_aa += genetic_code_ref[amino_acid][codon][0]
	
	freq_equ_aa = (1/number_of_aa)

	for amino_acid in genetic_code_ref : 
		if amino_acid in is_aa_in_req and len(genetic_code_ref[amino_acid]) > 1 and amino_acid != '*' : 

			freq_aa_ref = (occ_aa[amino_acid] / occ_tot_aa)
			deviation_score_aa[amino_acid] = freq_aa_ref - freq_equ_aa

	return deviation_score_aa


def calc_COUSIN_18(genetic_code, genetic_code_ref, deviation_score, is_aa_in_req, perform_COUSIN_with_boundaries) : #Changer genetic_code par freq_codon_seq (appel de la requête) et remplacer genet code par genetic code (dans le cas où reverse)

	weight_req = {}
	weight_ref = {}
	sum_COUSIN_aa_by_aa = 0

	COUSIN = 0

	if perform_COUSIN_with_boundaries == "Yes" or perform_COUSIN_with_boundaries == "YES" or perform_COUSIN_with_boundaries == "yes" or perform_COUSIN_with_boundaries == "Y" or perform_COUSIN_with_boundaries == "y" :

		for amino_acid in genetic_code_ref :
			if len(genetic_code_ref[amino_acid]) > 1 and amino_acid != '*' :
					if (is_aa_in_req[amino_acid] == True) :
						for codon in genetic_code_ref[amino_acid] :
							if amino_acid not in weight_req :
								weight_req[amino_acid] = (genetic_code[amino_acid][codon][1] * deviation_score[codon])
							else :
								weight_req[amino_acid] += (genetic_code[amino_acid][codon][1] * deviation_score[codon])
							if amino_acid not in weight_ref :
								weight_ref[amino_acid] = (genetic_code_ref[amino_acid][codon][1] * deviation_score[codon])
							else :
								weight_ref[amino_acid] += (genetic_code_ref[amino_acid][codon][1] * deviation_score[codon])

						if weight_ref[amino_acid] != 0 :
							if (weight_req[amino_acid]/weight_ref[amino_acid]) > 4 :
								sum_COUSIN_aa_by_aa += 4
							elif (weight_req[amino_acid]/weight_ref[amino_acid]) < - 3 :
								sum_COUSIN_aa_by_aa += - 3
							else :
								sum_COUSIN_aa_by_aa += round((weight_req[amino_acid]/weight_ref[amino_acid]),3)
						else :
							weight_ref.pop(amino_acid)
							weight_req.pop(amino_acid)

		if weight_req != 0 and weight_ref != 0 and len(weight_req) == len(weight_ref) :
			COUSIN = (sum_COUSIN_aa_by_aa / len(weight_req))
		else :
			COUSIN = 0
		return COUSIN

	else :

		for amino_acid in genetic_code_ref :
			if len(genetic_code_ref[amino_acid]) > 1 and amino_acid != '*' :
				if (is_aa_in_req[amino_acid] == True) :
					for codon in genetic_code_ref[amino_acid] :
						if amino_acid not in weight_req :
							weight_req[amino_acid] = (genetic_code[amino_acid][codon][1] * deviation_score[codon])
						else :
							weight_req[amino_acid] += (genetic_code[amino_acid][codon][1] * deviation_score[codon])
						if amino_acid not in weight_ref :
							weight_ref[amino_acid] = (genetic_code_ref[amino_acid][codon][1] * deviation_score[codon])
						else :
							weight_ref[amino_acid] += (genetic_code_ref[amino_acid][codon][1] * deviation_score[codon])

					if weight_ref[amino_acid] != 0 :
						sum_COUSIN_aa_by_aa += round((weight_req[amino_acid]/weight_ref[amino_acid]),3)
					else :
						weight_ref.pop(amino_acid)
						weight_req.pop(amino_acid)

		if weight_req != 0 and weight_ref != 0 and len(weight_req) == len(weight_ref) :
			COUSIN = (sum_COUSIN_aa_by_aa / len(weight_req))
		else :
			COUSIN = 0
		return COUSIN


def calc_freq_equ_muta_bias (genetic_code) :

	total_freq = 0
	for amino_acid in genetic_code :

		if amino_acid != '*' :

			for codon in genetic_code[amino_acid] :

				if codon[2] not in freq_nuc :
					freq_nuc[codon[2]] = genetic_code[amino_acid][codon][0]
					total_freq += genetic_code[amino_acid][codon][0]

				else :
					freq_nuc[codon[2]] += genetic_code[amino_acid][codon][0]
					total_freq += genetic_code[amino_acid][codon][0]


	for nucleotide in freq_nuc :
		freq_nuc[nucleotide] = freq_nuc[nucleotide] / total_freq


def generate_table_nucl_comp (genetic_code, output_dir) : #A CHECKER : CE QU'IL SE PASSE EN FONCTION D'UN CODON = 0. ICI, CALCUL EN POSSITION DES FREQUENCES DE CHAQUE CODON AU SEIN DE LA TABLE, CAR PRISE EN COMTPE DES CODONS = 1 (dont l'occurence doit être elle aussi modifier)

	sequential_order_cu_table = ["UUU", "UCU", "UAU", "UGU", "UUC", "UCC", "UAC", "UGC", "UUA", "UCA", "UAA", "UGA", "UUG", "UCG", "UAG", "UGG", "CUU", "CCU", "CAU", "CGU", "CUC", "CCC", "CAC", "CGC", "CUA", "CCA", "CAA", "CGA", "CUG", "CCG", "CAG", "CGG", "AUU", "ACU", "AAU", "AGU", "AUC", "ACC", "AAC", "AGC", "AUA", "ACA", "AAA", "AGA", "AUG", "ACG", "AAG", "AGG", "GUU", "GCU", "GAU", "GGU", "GUC", "GCC", "GAC", "GGC", "GUA", "GCA", "GAA", "GGA", "GUG", "GCG", "GAG", "GGG"]

	genetic_code_nucl_comp = deepcopy(genetic_code)

	total_occ_aa = {}
	total_occ = 0
	sum_frequencies_aa = {}
	sum_frequencies = 0
	total = 0
	cpt = 0

	for amino_acid in genetic_code :
		for codon in genetic_code[amino_acid] :
			if amino_acid not in total_occ_aa :
				total_occ_aa[amino_acid] = genetic_code[amino_acid][codon][0]
			else :
				total_occ_aa[amino_acid] += genetic_code[amino_acid][codon][0]
			if amino_acid not in sum_frequencies_aa :
				sum_frequencies_aa[amino_acid] = freq_nuc[codon[0]] * freq_nuc[codon[1]] * freq_nuc[codon[2]]
			else :
				sum_frequencies_aa[amino_acid] += freq_nuc[codon[0]] * freq_nuc[codon[1]] * freq_nuc[codon[2]]

		total_occ += total_occ_aa[amino_acid]
		sum_frequencies += sum_frequencies_aa[amino_acid]

	for amino_acid in genetic_code :
		temp_value = (occ_tot_aa[amino_acid] * (sum_frequencies_aa[amino_acid] / sum_frequencies))
		for codon in genetic_code[amino_acid] :
			genetic_code_nucl_comp[amino_acid][codon][1] = round((freq_nuc[codon[0]] * freq_nuc[codon[1]] * freq_nuc[codon[2]] / sum_frequencies_aa[amino_acid]),3)
			genetic_code_nucl_comp[amino_acid][codon][0] = round(genetic_code_nucl_comp[amino_acid][codon][1] * total_occ_aa[amino_acid])
			total += round(genetic_code_nucl_comp[amino_acid][codon][1] * genetic_code_nucl_comp[amino_acid][codon][0])

	return genetic_code_nucl_comp

def calc_ENC (genetic_code) : # A REVOIR CAR BUGGUEE SUR L'ADAPTATION A UN NOUVEAU CODE GENETIQUE

	ENC = 0
	aa_alone = 0

	category = {}
	sum_pow_frequencies = {}
	homozygosity = {}
	sum_homozygosity = {}
	N_value = {}
	mean = {}
	amount = {}
	tot_occ_aa = {}

	for amino_acid in genetic_code :
		if amino_acid != "*" :
			if len(genetic_code[amino_acid]) > 1 :
				for codon in genetic_code[amino_acid] :
					if amino_acid not in tot_occ_aa :
						tot_occ_aa[amino_acid] = genetic_code[amino_acid][codon][0]
					else :
						tot_occ_aa[amino_acid] += genetic_code[amino_acid][codon][0]
				if len(genetic_code[amino_acid]) not in category :
					category[len(genetic_code[amino_acid])] = {}
					category[len(genetic_code[amino_acid])][amino_acid] = {}
				else :
					category[len(genetic_code[amino_acid])][amino_acid] = {}
				if len(genetic_code[amino_acid]) not in amount :
					amount[len(genetic_code[amino_acid])] = 1
				else :
					amount[len(genetic_code[amino_acid])] += 1
			else :
				aa_alone += 1

	keyList = sorted(category.keys())

	for number in category :
		for amino_acid in category[number] :
			if tot_occ_aa[amino_acid] > 0 :
				for codon in genetic_code[amino_acid] :
					if genetic_code[amino_acid][codon][0] > 0 :
						if amino_acid not in sum_pow_frequencies :
							sum_pow_frequencies[amino_acid] = pow(genetic_code[amino_acid][codon][1], 2) # carré des fréquences au sein de l'aa
						else :
							sum_pow_frequencies[amino_acid] += pow(genetic_code[amino_acid][codon][1], 2) # carré des fréquences au sein de l'aa

				if (tot_occ_aa[amino_acid] - 1) > 0 : # CONDITION A VERIFIER : EMPECHE DIVISION PAR 0 SI UN SEUL CODON POUR UN ACIDE AMINE (evenement rare). OR SI UN SEUL ACIDE AMINE / UN SEUL CODON ==> homozygosity DE 1
					homozygosity[amino_acid] = (((tot_occ_aa[amino_acid] * sum_pow_frequencies[amino_acid]) - 1) / (tot_occ_aa[amino_acid] - 1))
					if ( homozygosity[amino_acid] == 0 ) :
						amount[number] -= 1
				else :
					homozygosity[amino_acid] = 1
			else :
				amount[number] -= 1

	for number in category :
		sum_homozygosity[number] = 0
		N_value[number] = None
		for amino_acid in category[number] :
			if amino_acid in homozygosity :
				sum_homozygosity[number] += homozygosity[amino_acid]
		if amount[number] != 0 :
			mean[number] = (sum_homozygosity[number] / amount[number])
		else :
			mean[number] = None
		if mean[number] != None and mean[number] != 0  :
			N_value[number] = (amount[number] / mean[number])
		else :
			N_value[number] = None

	for i in keyList :
		if N_value[i] == None :
			current_number = i
			h = 0
			j = 0
			previous_number = "NaN"
			next_number = "NaN"
			while h < current_number :
				if h in keyList :
					previous_number = h
				h += 1
			while j <= max(keyList) :
				if j > current_number and j in keyList :
					next_number = j
					break
				j += 1
			if previous_number == "NaN" or next_number == "NaN" :
				ENC = 0
				return ENC
			else :
				if N_value[previous_number] != None and N_value[next_number] != None :
					if ((N_value[previous_number] + N_value[next_number]) / 2) > i :
						N_value[i] = i
					else :
						N_value[i] = (N_value[previous_number] + N_value[next_number]) / 2
				else :
					ENC = 0
					return ENC

	for number in category :
		ENC += N_value[number]
	ENC += aa_alone
	return ENC


def calc_SCUO(genetic_code, genetic_code_ref) : #A VERIFIER (résultats legerement different de ceux de l'outil et je ne vois pas pourquoi.)

	entropy = {}
	entropy_max = {}
	SCUO = {}
	composition_ratio = {}
	total_occ_aa = {}
	total_occ = 0
	SCUO_score = 0
	cpt_null_codon = {}


	for amino_acid in genetic_code :

		if (len(genetic_code_ref[amino_acid]) > 1 and amino_acid != '*') :
			total_occ_aa[amino_acid] = 0
			for codon in genetic_code[amino_acid] :
				total_occ_aa[amino_acid] += genetic_code[amino_acid][codon][0]
				total_occ += genetic_code[amino_acid][codon][0]
				if genetic_code[amino_acid][codon][1] != 0 : # A vérifier dans le cas freq = 0 car ce n'est pas precise sur l'article de SCUO
					if amino_acid not in entropy :
						if amino_acid not in cpt_null_codon :
							cpt_null_codon[amino_acid] = 0
						entropy[amino_acid] = ( - round(genetic_code[amino_acid][codon][1],3) ) * (math.log(round(genetic_code[amino_acid][codon][1],3)))
					else :
						entropy[amino_acid] += ( - round(genetic_code[amino_acid][codon][1],3) ) * (math.log(round(genetic_code[amino_acid][codon][1],3)))
				else :
					if amino_acid not in cpt_null_codon :
						cpt_null_codon[amino_acid] = 1
					else :
						cpt_null_codon[amino_acid] += 1

			if amino_acid in entropy :
				cpt_null_codon[amino_acid] = len(genetic_code[amino_acid]) - cpt_null_codon[amino_acid]
				entropy_max[amino_acid] = ( - math.log((1/cpt_null_codon[amino_acid])))
				if entropy_max[amino_acid] != 0 :
					SCUO[amino_acid] = ((round(entropy_max[amino_acid],3) - round(entropy[amino_acid],3)) / round(entropy_max[amino_acid],3))


	if total_occ > 0 :
		for amino_acid in genetic_code :

			if (len(genetic_code_ref[amino_acid]) > 1 and amino_acid != '*') :
				composition_ratio[amino_acid] = total_occ_aa[amino_acid] / total_occ
				if amino_acid in SCUO :
					SCUO_score += ( composition_ratio[amino_acid] * SCUO[amino_acid])
	else :
		print ("###IMPORTANT (SCUO): your sequence seems empty or isn't compounded of amino acids associated with synonymous codons. SCUO score = 0.\n ###")
		SCUO_score = 0

	return SCUO_score

def calc_optimal_codons(genetic_code, optimal_codons_file) :
	best_codons = ""
	optimal_codons_list = {}

	if optimal_codons_file == "nothing" :
		print("## IMPORTANT : You didn't give a list indicating optimal codons. Determining optimal codons with the reference (ones with greatest frequency) ##")
		for amino_acid in genetic_code :
			if amino_acid != "*" and len(genetic_code[amino_acid]) > 1 :
				best_number = 0
				for codon in genetic_code[amino_acid] :
					if genetic_code[amino_acid][codon][0] > best_number :
						best_number = genetic_code[amino_acid][codon][0]
						best_codons = codon
					elif genetic_code[amino_acid][codon][0] == best_number and best_number > 0 :
						best_codons += "," + codon
				if "," in best_codons :
					multiple_codons = best_codons.split(",")
					print ("## Note : We have found ties in your codon usage table for " + amino_acid + " please note that all tied codons will be taken in the calculation of the FOP ##" )
					for codon in multiple_codons :
						optimal_codons_list[codon] = ""
				else :
					optimal_codons_list[best_codons] = ""

	else :
		if "," in optimal_codons_file :
			optimal_codons_file = optimal_codons_file.replace("T","U")
			multiple_codons = optimal_codons_file.split(",")
			for codon in multiple_codons :
				result = re.match("^[A|U|C|G|a|u|c|g]{3}", codon)
				if result == None :
					print("## error on your optimal codon dataset ; it must be codons in a DNA/RNA alphabet. A selection of optimal codons is now done with reference dataset. ##")
					optimal_codons_file = "nothing"
					optimal_codons_list = calc_optimal_codons(genetic_code, optimal_codons_file)

			optimal_codons_list = {}
			for amino_acid in genetic_code :
				cpt = 0
				if amino_acid != "*" and len(genetic_code[amino_acid]) > 1 :
					for codon in genetic_code[amino_acid] :
						if codon in multiple_codons :
							optimal_codons_list[codon] = ""
							cpt += 1
					if cpt > 1 :
						print ("## Warning, the amino acid " + amino_acid + " contains more than one optimal codon. Please be aware of this specificity. All of the related codons are kept as optimal for the analysis. ##")
				else :
					for codon in genetic_code[amino_acid] :
						if codon in multiple_codons :
							print(" ## error on your optimal codon dataset, the codon " + codon + " is either a STOP codon or a codon without synonymous. Codon ejected of the optimal codons list. ##")
		else :
			print("## error on your optimal codon dataset, each codon must be separated by \",\" characters (and no spaces between codons). A selection of optimal codons is now done with reference dataset. ##")
			optimal_codons_file = "nothing"
			optimal_codons_list = calc_optimal_codons(genetic_code, optimal_codons_file)

	return optimal_codons_list

def calc_FOP(genetic_code, optimal_codons) :

	occ_optimal_codons = 0
	occ_tot_FOP = 0

	for amino_acid in genetic_code :
		if amino_acid != "*" and len(genetic_code[amino_acid]) > 1 :
			for codon in genetic_code[amino_acid] :
				if codon in optimal_codons :
					occ_optimal_codons += genetic_code[amino_acid][codon][0]
					occ_tot_FOP += genetic_code[amino_acid][codon][0]
				else :
					occ_tot_FOP += genetic_code[amino_acid][codon][0]

	if occ_tot_FOP != 0 :
		FOP = (occ_optimal_codons / occ_tot_FOP)
		return FOP
	else :
		print ("## IMPORTANT (FOP): your sequence seems empty, FOP score = 0 ##")
		FOP = 0
		return FOP

def calc_CBI(genetic_code, optimal_codons) :
	expected_number = 0
	tot_occ_optimal_codons = 0
	occ_tot_CBI = 0
	occ_codon_optimal = {}
	occ_aa_CBI = {}
	number_of_opt_codons = {}

	CBI = 0


	for amino_acid in genetic_code :
		if amino_acid != "*" and len(genetic_code[amino_acid]) > 1 :
			for codon in genetic_code[amino_acid] :
				if codon in optimal_codons :
					if amino_acid not in occ_codon_optimal :
						occ_codon_optimal[amino_acid] = genetic_code[amino_acid][codon][0]
					else :
						occ_codon_optimal[amino_acid] += genetic_code[amino_acid][codon][0]
					if amino_acid not in number_of_opt_codons :
						number_of_opt_codons[amino_acid] = 1
					else :
						number_of_opt_codons[amino_acid] += 1
					if amino_acid not in occ_aa_CBI :
						occ_aa_CBI[amino_acid] = genetic_code[amino_acid][codon][0]
					else :
						occ_aa_CBI[amino_acid] += genetic_code[amino_acid][codon][0]
					tot_occ_optimal_codons += genetic_code[amino_acid][codon][0]
					occ_tot_CBI += genetic_code[amino_acid][codon][0]
				else :
					if amino_acid not in occ_aa_CBI :
						occ_aa_CBI[amino_acid] = genetic_code[amino_acid][codon][0]
					else :
						occ_aa_CBI[amino_acid] += genetic_code[amino_acid][codon][0]
					occ_tot_CBI += genetic_code[amino_acid][codon][0]

	for amino_acid in occ_aa_CBI :
		if amino_acid in occ_codon_optimal :
			#expected_number += (occ_aa_CBI[amino_acid]) * ((occ_codon_optimal[amino_acid])/(len(genetic_code[amino_acid]))) #Ce qui semble indiqué, de manière floue, dans l'artcile de Roth et al. En réalité, le numérateur semble être le nombre qualitatif de codons optimaux (1 ou plus donc)
			expected_number += (occ_aa_CBI[amino_acid]) * ((number_of_opt_codons[amino_acid])/(len(genetic_code[amino_acid]))) #Avec cette version, cela semble plus logique. Le "expected number" indique le nombre d'acides aminés encodés par les codons optimaux en absence de biais d'usage des codons.


	if occ_tot_CBI - expected_number != 0 :
		CBI = ((tot_occ_optimal_codons - expected_number) / (occ_tot_CBI - expected_number))
	return CBI


def calc_RSCU(genetic_code)  :

	RSCU = {}
	occ_aa_RSCU = {}

	for amino_acid in genetic_code :
		if amino_acid != "*" and len(genetic_code[amino_acid]) > 1 :
			for codon in genetic_code[amino_acid] :
				if amino_acid not in occ_aa_RSCU :
					occ_aa_RSCU[amino_acid] = genetic_code[amino_acid][codon][0]
				else :
					occ_aa_RSCU[amino_acid] += genetic_code[amino_acid][codon][0]

	for amino_acid in genetic_code :
		if amino_acid != "*" and len(genetic_code[amino_acid]) > 1 :
			if occ_aa_RSCU[amino_acid] > 0 :
				for codon in genetic_code[amino_acid] :
					if amino_acid not in RSCU :
						RSCU[amino_acid] = {}
					RSCU[amino_acid][codon] = ((genetic_code[amino_acid][codon][0]) / (((1)/(len(genetic_code[amino_acid])))*(occ_aa_RSCU[amino_acid])))
			else :
				for codon in genetic_code[amino_acid] :
					if amino_acid not in RSCU :
						RSCU[amino_acid] = {}
					RSCU[amino_acid][codon] = 0

	return RSCU


def calc_ICDI(genetic_code, RSCU_data, is_aa_in_req) :

	S_value = {}
	S_value_tot = 0
	cpt_aa = 0

	for amino_acid in genetic_code :
		if amino_acid != "*" and len(genetic_code[amino_acid]) > 1 and is_aa_in_req[amino_acid] == True :
			cpt_aa += 1
			sum_RSCU = 0
			for codon in genetic_code[amino_acid] :
					sum_RSCU += pow((RSCU_data[amino_acid][codon] - 1),2)
			S_value[amino_acid] = ((1)/((len(genetic_code[amino_acid])*(len(genetic_code[amino_acid]) - 1)))) * sum_RSCU
			S_value_tot += S_value[amino_acid]

	ICDI = (((1) / (cpt_aa)) * (S_value_tot))
	return ICDI

def calc_scaled_chi(genetic_code) :

	sum_codons = None
	total_occ_scaled_chi = 0
	scaled_chi = 0

	for amino_acid in genetic_code :
		if amino_acid != "*" and len(genetic_code[amino_acid]) > 1 :
			for codon in genetic_code[amino_acid] :
				if (genetic_code[amino_acid][codon][0] > 0) :
					if sum_codons == None :
						sum_codons = ((genetic_code[amino_acid][codon][0] - (len(genetic_code[amino_acid]) - 1)) / (len(genetic_code[amino_acid]) - 1))
						total_occ_scaled_chi += genetic_code[amino_acid][codon][0]
					else :
						sum_codons += ((genetic_code[amino_acid][codon][0] - (len(genetic_code[amino_acid]) - 1)) / (len(genetic_code[amino_acid]) - 1))
						total_occ_scaled_chi += genetic_code[amino_acid][codon][0]

	scaled_chi = (sum_codons / total_occ_scaled_chi)

	return scaled_chi


def calc_GRAVY (genetic_code) :

	GRAVY_score = 0
	total_occ_aa_gravy = 0

	for amino_acid in genetic_code :
		if amino_acid != '*' :
			for codon in genetic_code[amino_acid] :
				GRAVY_score += (genetic_code[amino_acid][codon][0] * hydropathy_value[amino_acid])
				total_occ_aa_gravy += genetic_code[amino_acid][codon][0]

	GRAVY_score = (GRAVY_score / total_occ_aa_gravy)

	return GRAVY_score

def calc_AROMA (genetic_code) :

	AROMA_score = 0
	occ_aa_aromatic = {}
	total_occ_aa_AROMA = 0

	for amino_acid in genetic_code :
		if amino_acid != '*' :
			if amino_acid == 'F' or amino_acid == 'Y' or amino_acid == 'W' :
				for codon in genetic_code[amino_acid] :
					if amino_acid not in occ_aa_aromatic :
						occ_aa_aromatic[amino_acid] = genetic_code[amino_acid][codon][0]
						total_occ_aa_AROMA += genetic_code[amino_acid][codon][0]
					else :
						occ_aa_aromatic[amino_acid] += genetic_code[amino_acid][codon][0]
						total_occ_aa_AROMA += genetic_code[amino_acid][codon][0]
			else :
				for codon in genetic_code[amino_acid] :
					total_occ_aa_AROMA += genetic_code[amino_acid][codon][0]

	for amino_acid in occ_aa_aromatic :
		AROMA_score += (occ_aa_aromatic[amino_acid] / total_occ_aa_AROMA)

	return AROMA_score

################################################################################
################################################################################
################################################################################
############################ R CALCULATIONS ####################################
################################################################################
################################################################################
################################################################################

matrice = {}
mean = {}
sd = {}
rnorm = {}

quantile_99_left_x = {}
quantile_99_right_x = {}
quantile_95_left_x = {}
quantile_95_right_x = {}

quantile_99_left_y = {}
quantile_99_right_y = {}
quantile_95_left_y = {}
quantile_95_right_y = {}

density_result = {}

def get_matrice_and_quantiles(input_value, dict_values, header, header_file, results_file, output_dir) :

	global firstline

	global matrice
	global quantile_99_left_x
	global quantile_99_right_x
	global quantile_95_left_x
	global quantile_95_right_x
	global quantile_99_left_y
	global quantile_99_right_y
	global quantile_95_left_y
	global quantile_95_right_y

	matrice[header] = np.array(dict_values)

	quantile_99_left_x[header] = np.percentile(matrice[header], 0.5)
	quantile_99_right_x[header] = np.percentile(matrice[header], 99.5)
	quantile_95_left_x[header] = np.percentile(matrice[header], 2.5)
	quantile_95_right_x[header] = np.percentile(matrice[header], 97.5)

	if firstline :
		csv_fill = "\t" + csv_name
		results_file["firstline"] += csv_fill
		firstline = False

	if input_value < quantile_99_left_x[header] and input_value < quantile_95_left_x[header] and input_value < quantile_99_right_x[header] and input_value < quantile_95_right_x[header] :
		results_file[header_file] += "\t-2"
	elif input_value >= quantile_99_left_x[header] and input_value < quantile_95_left_x[header] and input_value < quantile_99_right_x[header] and input_value < quantile_95_right_x[header] :
		results_file[header_file] += "\t-1"
	elif input_value > quantile_99_left_x[header] and input_value >= quantile_95_left_x[header] and input_value < quantile_99_right_x[header] and input_value <= quantile_95_right_x[header] :
		results_file[header_file] += "\t0"
	elif input_value > quantile_99_left_x[header] and input_value > quantile_95_left_x[header] and input_value <= quantile_99_right_x[header] and input_value > quantile_95_right_x[header] :
		results_file[header_file] += "\t1"
	elif input_value > quantile_99_left_x[header] and input_value > quantile_95_left_x[header] and input_value > quantile_99_right_x[header] and input_value > quantile_95_right_x[header] :
		results_file[header_file] += "\t2"
	else : 
		if not os.path.isdir("./" + output_dir + "/log"):
			os.makedirs("./" + output_dir + "/log")

		log_file = open("./" + output_dir + "/log/log_cousin.txt", "a")
		log_file.write("## SIMULATION STEP ERROR (step : "+ csv_name +") : 95% and 99% intervals limits seem to be merged (at least on one side). Please consider to check your related sequence and intervals limits manually (" + str(header) + ")\n")
		results_file[header_file] += "\tNA"
	

################################################################################
################################################################################
################################################################################
############################### R GRAPHS #######################################
################################################################################
################################################################################
################################################################################

################################################################################
######################## CREATE DENSITY CURVES #################################
################################################################################

def create_graphs(output_dir, graph_dir, input_value, legends, lengths, type_of_data, *args) : #AJOUTER INPUT VALUE DANS LES ARGS ET TRAITER CORRECTEMENT. VALEUR UNIQUE DE TYPE FLOAT # pourrait être une liste incrémentée de manière ephémère au sein de la fonction qui appelle get_matrice_and_quantiles

	density_result = []
	xlim_min_graph = None
	xlim_max_graph = None
	ylim_min_graph = None
	ylim_max_graph = None
	x_min = []
	x_max = []
	y_min = []
	y_max = []

	create_first_plot = 0

	i = 0
	if input_value == None :
		sns.set(rc={'figure.figsize':(14,10)})
		sns.set(color_codes=True)
		sns.set_style("whitegrid")

		set_of_colors_colorblind = sns.color_palette("colorblind")
		set_of_colors_colorblind.extend(["#b15928","#9b59b6","#fb9a99","#ffff99","#6a3d9a","#1f78b4","#95a5a6","#34495e", "#2ecc71"])
		sns.set_palette(sns.color_palette(set_of_colors_colorblind, n_colors=15))

		while i < len(args) :

			sns_plot = sns.distplot(args[i], hist=False, kde_kws={"label":legends[i] + " : " + str(lengths[i]) + " " + type_of_data +  "\n" + "95% : [ " + str(round(np.percentile(args[i], 2.5),3)) + " ; " + str(round(np.percentile(args[i], 97.5),3)) + " ]" + "\n" + "99% : [ " + str(round(np.percentile(args[i], 0.5),3)) + " ; " + str(round(np.percentile(args[i], 99.5),3)) + " ]" +  "\n"  })
			plt.axvline(x=np.percentile(args[i], 2.5), ls = '--', color = sns_plot.get_lines()[-1].get_color())
			plt.axvline(x=np.percentile(args[i], 97.5), ls = '--', color = sns_plot.get_lines()[-1].get_color())
			box = sns_plot.get_position() # get position of figure
			plt.subplots_adjust(hspace=0.4, wspace=0.4)
			plt.legend(fontsize = 12)
			plt.xticks(fontsize=14)
			plt.yticks(fontsize=14)
			sns_plot.set_position([box.x0, box.y0, box.width * 0.90, box.height])
			plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
			plt.setp(sns_plot.get_legend().get_texts(), fontsize='16')
			i += 1

	elif input_value != None :
		sns.set(rc={'figure.figsize':(14,10)})
		sns.set(color_codes=True)
		sns.set_style("whitegrid")

		set_of_colors_colorblind = sns.color_palette("colorblind")
		set_of_colors_colorblind.extend(["#b15928","#9b59b6","#fb9a99","#ffff99","#6a3d9a","#1f78b4","#95a5a6","#34495e", "#2ecc71"])
		sns.set_palette(sns.color_palette(set_of_colors_colorblind, n_colors=15))

		while i < len(args) :

			sns_plot = sns.distplot(args[i], hist=False, kde_kws={"label":legends[i] + " : " + str(lengths[i]) + " " + type_of_data +  "\n" + "95% : [ " + str(round(np.percentile(args[i], 2.5),3)) + " ; " + str(round(np.percentile(args[i], 97.5),3)) + " ]" + "\n" + "99% : [ " + str(round(np.percentile(args[i], 0.5),3)) + " ; " + str(round(np.percentile(args[i], 99.5),3)) + " ]" +  "\n"  });
			plt.axvline(x=np.percentile(args[i], 2.5), ls = '--', color = sns_plot.get_lines()[-1].get_color())
			plt.axvline(x=np.percentile(args[i], 97.5), ls = '--', color = sns_plot.get_lines()[-1].get_color())
			box = sns_plot.get_position() 
			plt.subplots_adjust(hspace=0.4, wspace=0.4)
			plt.xticks(fontsize=14)
			plt.yticks(fontsize=14)
			sns_plot.set_position([box.x0, box.y0, box.width * 0.84, box.height])
			plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
			plt.setp(sns_plot.get_legend().get_texts(), fontsize='16')
			i += 1

		xlim_min_graph, xlim_max_graph = sns_plot.get_xlim()
		ylim_min_graph, ylim_max_graph = sns_plot.get_ylim()

		sns_plot.set_xlim([xlim_min_graph,xlim_max_graph])
		sns_plot.set_ylim([ylim_min_graph,ylim_max_graph])

		if input_value < xlim_min_graph or input_value > xlim_max_graph :
			if input_value < xlim_min_graph :

				sns_plot.plot((xlim_min_graph + abs(((xlim_max_graph - xlim_min_graph)/50))), (ylim_min_graph + abs(((ylim_max_graph - ylim_min_graph)/50))), color= "red", marker="$\u2190$ ", markersize=25)
				sns_plot.annotate(round(input_value,3), xy=((xlim_min_graph + abs(((xlim_max_graph - xlim_min_graph)/50))), (ylim_min_graph + abs(((ylim_max_graph - ylim_min_graph)/50)))), xycoords='data', xytext=(xlim_min_graph,(ylim_min_graph + abs(((ylim_max_graph - ylim_min_graph)/10)))), horizontalalignment="center", arrowprops=dict(arrowstyle="-", connectionstyle="arc3"), bbox=dict(boxstyle="round", fc="w"),)

			else :

				sns_plot.plot((xlim_max_graph - abs(((xlim_max_graph - xlim_min_graph)/50))), (ylim_min_graph + abs(((ylim_max_graph - ylim_min_graph)/50))), color= "red", marker="$\u2192$ ", markersize=25)
				sns_plot.annotate(round(input_value,3), xy=((xlim_max_graph - abs(((xlim_max_graph - xlim_min_graph)/50))), (ylim_min_graph + abs(((ylim_max_graph - ylim_min_graph)/50)))), xycoords='data', xytext=(xlim_max_graph,(ylim_min_graph + abs(((ylim_max_graph - ylim_min_graph)/10)))), horizontalalignment="center", arrowprops=dict(arrowstyle="-", connectionstyle="arc3"), bbox=dict(boxstyle="round", fc="w"),)

		else :
			sns_plot.plot(input_value, (ylim_min_graph + abs(((ylim_max_graph - ylim_min_graph)/50))), color= "black", marker="$\u2193$ ", markersize=25)
			sns_plot.annotate(round(input_value,3), xy=((input_value, (ylim_min_graph + abs(((ylim_max_graph - ylim_min_graph)/50))))), xycoords='data', xytext=(input_value,(ylim_min_graph + abs(((ylim_max_graph - ylim_min_graph)/10)))), horizontalalignment="center", arrowprops=dict(arrowstyle="-", connectionstyle="arc3"), bbox=dict(boxstyle="round", fc="w"),)

	plt.suptitle(step_title, fontsize=22)
	plt.title(title_name, fontsize = 18)
	sns_plot.set_xlabel(absc_name, fontsize = 18)
	sns_plot.set_ylabel(ord_name, fontsize = 18)
	fig = sns_plot.get_figure()
	fig.savefig("./" + output_dir + "/" + graph_dir + "/" + file_name + '.pdf')
	plt.close(fig)


	i = 0

colors = ["#812321","#50d885","#5e2f8f","#73c85e","#4d68d7","#acba3b","#5d84f4","#e2ac3d","#5e55ba","#79a033","#b370d6","#3e7f27","#db6ac6","#4bc98d","#983082","#37d8b0","#d54a4a","#36dee6","#b84229","#4493d9","#b7891a","#6c81de","#8ac166","#57286e","#5cb66b","#d14b75","#47b795","#892254","#3f9b66","#d2669e","#2f6b2a","#d179bb","#7fb26d","#35418b","#bcaf4f","#6464b0","#c16420","#9792e2","#686611","#bb82cd","#aeb867","#d35760","#827e34","#ba5b6e","#d7b16d","#894e19","#cd8f40","#ce7058","#c2894e","#df7d52"]

################################################################################
################################################################################
################################################################################
########################## LAUNCHING FUNCTIONS #################################
################################################################################
################################################################################
################################################################################

################################################################################
################################ VARIABLES #####################################
################################################################################

genetic_code_ref = {}
genetic_code_ref_rev = {}
genetic_code_sim = {}
genetic_code_sim_nucl_comp = {}
content_seq = {}
codon_gc_max = {}
codon_at_max = {}
occ_codon = {}
occ_codon_seq = {}
w_scores = {}
freq_equ = {}
diff_codon_seq = {}
occ_tot_ref = {}
occ_tot = 0
occ_tot_aa_4_sim = 0
occ_tot_aa = {}
occ_aa_seq = {}
dict_freq_equ_muta_bias = {}

input_seq_value_CAI_cod = {}
input_seq_value_CAI_aa = {}
input_seq_value_COUSIN = {}
input_seq_value_COUSIN_pond = {}
input_seq_value_ENC = {}
input_seq_value_SCUO = {}
input_seq_value_FOP = {}
input_seq_value_CBI = {}
input_seq_value_ICDI = {}
input_seq_value_scaled_chi = {}
input_seq_value_GRAVY = {}
input_seq_value_AROMA = {}
input_seq_value_eucl_dist_aa_que_vs_ref = {}
input_seq_value_eucl_dist_aa_que_vs_H0 = {}


freq_nuc = {}

hydropathy_value =  {'F': 2.800 , 'S': -0.800 , 'Y': -1.300 , 'C': 2.500 , 'L': 3.800 , 'W': -0.900 , 'H': -3.200 , 'R': -4.500 , 'Q': -3.500 , 'I': 4.500 , 'T': -0.700 , 'N': -3.500 , 'K': -3.900 , 'M': 1.900 , 'V': 4.200 , 'A': 1.800 , 'D': -3.500 , 'G': -0.400 , 'E': -3.500 , 'P' : -1.600 }

#########################lancement fonctions####################################

CAI_codon_calc = {}
CAI_aa_calc = {}
COUSIN_codon_calc = {}
COUSIN_calc = {}
COUSIN_calc_pond = {}
ENC_calc = {}
SCUO_calc = {}
FOP_calc = {}
CBI_calc = {}
ICDI_calc = {}
scaled_chi_calc = {}

input_sim_CAI_codon_calc = {}
input_sim_CAI_aa_calc = {}
input_sim_COUSIN_calc = {}
input_sim_COUSIN_pond_calc = {}

input_sim_ENC = {}
input_sim_SCUO = {}
input_sim_FOP = {}
input_sim_CBI = {}
input_sim_ICDI = {}
input_sim_scaled_chi = {}

input_sim_CAI_codon_calc_rand_aa = {}
input_sim_CAI_aa_calc_rand_aa = {}
input_sim_COUSIN_calc_rand_aa = {}
input_sim_COUSIN_pond_calc_rand_aa = {}

input_sim_ENC_rand_aa = {}
input_sim_SCUO_rand_aa = {}
input_sim_FOP_rand_aa = {}
input_sim_CBI_rand_aa = {}
input_sim_ICDI_rand_aa = {}
input_sim_scaled_chi_rand_aa = {}

################################################################################
############### FUNCTION : CALCULATE VALUES FROM INPUT #########################
################################################################################

set_of_seq = {"cousin_value" : [], "cousin_value_pond" : [], "CAI_codon" : [], "CAI_cod_famil" : [], "ENC_value" : [], "SCUO_value" : [], "FOP_value" : [], "CBI_value" : [], "ICDI_value" : [], "ICDI_value_nucl_comp" : [], "scaled_chi_value" : [], "GRAVY_value" : [], "AROMA_value" : [], "eucl_dist_aa_value_que_vs_ref" : [], "eucl_dist_aa_value_que_vs_H0" : []}

results_for_file = {}
results_for_file_optimization = {}
results_for_file_simulation = {}
results_for_file_clustering = {}
results_for_file_data_comparison = {}

sequential_order_results = []
sequential_order_results_optimization = []
sequential_order_results_simulation = []
sequential_order_results_clustering = []
sequential_order__data_comparison = []

def send_results () :
	return results_for_file

def send_results_optimization () :
	return results_for_file_optimization

def send_results_simulation () :
	return results_for_file_simulation

def send_results_clustering () :
	return results_for_file_clustering

def send_results_data_comparison () :
	return results_for_file_data_comparison

def send_results_list () :
	return sequential_order_results

def send_results_optimization_list () :
	return sequential_order_results_optimization

def send_results_simulation_list () :
	return sequential_order_results_simulation

def send_results_clustering_list () :
	return sequential_order_results_clustering

def send_results_data_comparison_list () :
	return sequential_order__data_comparison

def write_file (path_file, data, data_list) :

	results_file = open(path_file, 'a')

	for header_seq in data_list :
		if "\n" not in data[header_seq] :
			results_file.write(data[header_seq] + "\n")
		else :
			results_file.write(data[header_seq])


	results_file.close()

def write_file_optimization (path_file, data, data_list) :

	results_file = open(path_file, 'a')

	for header_seq in data_list :
		if "\n" not in data[header_seq] :
			results_file.write(data[header_seq] + "\n")
		else :
			results_file.write(data[header_seq])

	results_file.close()

def write_file_simulation (path_file, data, data_list) :

	results_file = open(path_file, 'a')

	for header_seq in data_list :
		if "\n" not in data[header_seq] :
			results_file.write(data[header_seq] + "\n")
		else :
			results_file.write(data[header_seq])

	results_file.close()

def write_file_clustering (path_file, data, data_list) :

	results_file = open(path_file, 'a')

	for header_seq in data_list :
		if "\n" not in data[header_seq] :
			results_file.write(data[header_seq] + "\n")
		else :
			results_file.write(data[header_seq])

	results_file.close()

def create_data_vector(cu_dict, di_cu_dict, RSCU_data, header, output_dir) :

	sequential_order_amino_acid = ["M","W","N","D","C","Q","E","H","K","F","Y","I","A","G","P","T","V","R","L","S","*"]
	
	sequential_order_codon = ["AUG", "UGG", "AAU", "AAC", "GAU", "GAC", "UGU", "UGC", "CAA", "CAG", "GAA", "GAG", "CAU", "CAC", "AAA", "AAG", "UUU", "UUC", "UAU", "UAC", "AUU", "AUC", "AUA", "GCU", "GCC", "GCA", "GCG", "GGU", "GGC", "GGA", "GGG", "CCU", "CCC", "CCA", "CCG", "ACU", "ACC", "ACA", "ACG", "GUU", "GUC", "GUA", "GUG", "CGU", "CGC", "CGA", "CGG", "AGA", "AGG", "UUA", "UUG", "CUU", "CUC", "CUA", "CUG", "UCU", "UCC", "UCA", "UCG", "AGU", "AGC", "UAG", "UAA", "UGA"]
	total_list_dicodon = [str(first_codon) + "_" + str(second_codon) for first_codon in sequential_order_codon for second_codon in sequential_order_codon]

	if os.path.isfile(output_dir + "/vector_data/vector_occurences.txt") :

		vector_file = open(output_dir + "/vector_data/vector_occurences.txt", "a")
		vector_file.write(str(header) + "\t")

		for amino_acid in sequential_order_amino_acid :
			if amino_acid in cu_dict :
				for codon in sequential_order_codon :
					if codon in cu_dict[amino_acid] :
						vector_file.write(str(cu_dict[amino_acid][codon][0]) + "\t")

		vector_file.write("\n")
		vector_file.close()

	else :

		vector_file = open(output_dir + "/vector_data/vector_occurences.txt", "w")
		vector_file.write("Header" + "\t" + "AUG" + "\t" + "UGG" + "\t" + "AAU" + "\t" + "AAC" + "\t" + "GAU" + "\t" + "GAC" + "\t" + "UGU" + "\t" + "UGC" + "\t" + "CAA" + "\t" + "CAG" + "\t" + "GAA" + "\t" + "GAG" + "\t" + "CAU" + "\t" + "CAC" + "\t" + "AAA" + "\t" + "AAG" + "\t" + "UUU" + "\t" + "UUC" + "\t" + "UAU" + "\t" + "UAC" + "\t" + "AUU" + "\t" + "AUC" + "\t" + "AUA" + "\t" + "GCU" + "\t" + "GCC" + "\t" + "GCA" + "\t" + "GCG" + "\t" + "GGU" + "\t" + "GGC" + "\t" + "GGA" + "\t" + "GGG" + "\t" + "CCU" + "\t" + "CCC" + "\t" + "CCA" + "\t" + "CCG" + "\t" + "ACU" + "\t" + "ACC" + "\t" + "ACA" + "\t" + "ACG" + "\t" + "GUU" + "\t" + "GUC" + "\t" + "GUA" + "\t" + "GUG" + "\t" + "CGU" + "\t" + "CGC" + "\t" + "CGA" + "\t" + "CGG" + "\t" + "AGA" + "\t" + "AGG" + "\t" + "UUA" + "\t" + "UUG" + "\t" + "CUU" + "\t" + "CUC" + "\t" + "CUA" + "\t" + "CUG" + "\t" + "UCU" + "\t" + "UCC" + "\t" + "UCA" + "\t" + "UCG" + "\t" + "AGU" + "\t" + "AGC" + "\t" + "UAG" + "\t" + "UAA" + "\t" + "UGA" + "\n")
		vector_file.write(str(header) + "\t")

		for amino_acid in sequential_order_amino_acid :
			if amino_acid in cu_dict :
				for codon in sequential_order_codon :
					if codon in cu_dict[amino_acid] :
						vector_file.write(str(cu_dict[amino_acid][codon][0]) + "\t")

		vector_file.write("\n")
		vector_file.close()


	if os.path.isfile(output_dir + "/vector_data/vector_frequencies.txt") :

		vector_file = open(output_dir + "/vector_data/vector_frequencies.txt", "a")
		vector_file.write(str(header) + "\t")

		for amino_acid in sequential_order_amino_acid :
			if amino_acid in cu_dict :
				for codon in sequential_order_codon :
					if codon in cu_dict[amino_acid] :
						vector_file.write(str(round(cu_dict[amino_acid][codon][1],3)) + "\t")

		vector_file.write("\n")
		vector_file.close()

	else :
		
		vector_file = open(output_dir + "/vector_data/vector_frequencies.txt", "w")
		vector_file.write("Header" + "\t" + "AUG" + "\t" + "UGG" + "\t" + "AAU" + "\t" + "AAC" + "\t" + "GAU" + "\t" + "GAC" + "\t" + "UGU" + "\t" + "UGC" + "\t" + "CAA" + "\t" + "CAG" + "\t" + "GAA" + "\t" + "GAG" + "\t" + "CAU" + "\t" + "CAC" + "\t" + "AAA" + "\t" + "AAG" + "\t" + "UUU" + "\t" + "UUC" + "\t" + "UAU" + "\t" + "UAC" + "\t" + "AUU" + "\t" + "AUC" + "\t" + "AUA" + "\t" + "GCU" + "\t" + "GCC" + "\t" + "GCA" + "\t" + "GCG" + "\t" + "GGU" + "\t" + "GGC" + "\t" + "GGA" + "\t" + "GGG" + "\t" + "CCU" + "\t" + "CCC" + "\t" + "CCA" + "\t" + "CCG" + "\t" + "ACU" + "\t" + "ACC" + "\t" + "ACA" + "\t" + "ACG" + "\t" + "GUU" + "\t" + "GUC" + "\t" + "GUA" + "\t" + "GUG" + "\t" + "CGU" + "\t" + "CGC" + "\t" + "CGA" + "\t" + "CGG" + "\t" + "AGA" + "\t" + "AGG" + "\t" + "UUA" + "\t" + "UUG" + "\t" + "CUU" + "\t" + "CUC" + "\t" + "CUA" + "\t" + "CUG" + "\t" + "UCU" + "\t" + "UCC" + "\t" + "UCA" + "\t" + "UCG" + "\t" + "AGU" + "\t" + "AGC" + "\t" + "UAG" + "\t" + "UAA" + "\t" + "UGA" + "\n")
		vector_file.write(str(header) + "\t")

		for amino_acid in sequential_order_amino_acid :
			if amino_acid in cu_dict :
				for codon in sequential_order_codon :
					if codon in cu_dict[amino_acid] :
						vector_file.write(str(round(cu_dict[amino_acid][codon][1],3)) + "\t")

		vector_file.write("\n")
		vector_file.close()

	if os.path.isfile(output_dir + "/vector_data/di_vector_occurences.txt") :

		vector_file = open(output_dir + "/vector_data/di_vector_occurences.txt", "a")
		vector_file.write(str(header) + "\t")

		for first_codon in sequential_order_codon :
			for second_codon in sequential_order_codon :
				if first_codon in di_cu_dict : 
					if second_codon in di_cu_dict[first_codon] :
						vector_file.write(str(di_cu_dict[first_codon][second_codon][0]) + "\t")
				else : 
					vector_file.write("NA" + "\t")

		vector_file.write("\n")
		vector_file.close()

	else :

		vector_file = open(output_dir + "/vector_data/di_vector_occurences.txt", "w")
		vector_file.write("Header" + "\t" + "\t".join(total_list_dicodon) + "\n")
		vector_file.write(str(header) + "\t")

		for first_codon in sequential_order_codon :
			for second_codon in sequential_order_codon :
				if first_codon in di_cu_dict : 
					if second_codon in di_cu_dict[first_codon] :
						vector_file.write(str(di_cu_dict[first_codon][second_codon][0]) + "\t")
				else : 
					vector_file.write("NA" + "\t")

		vector_file.write("\n")
		vector_file.close()


	if os.path.isfile(output_dir + "/vector_data/di_vector_frequencies.txt") :

		vector_file = open(output_dir + "/vector_data/di_vector_frequencies.txt", "a")
		vector_file.write(str(header) + "\t")

		for first_codon in sequential_order_codon :
			for second_codon in sequential_order_codon :
				if first_codon in di_cu_dict : 
					if second_codon in di_cu_dict[first_codon] :
						vector_file.write(str(round(di_cu_dict[first_codon][second_codon][1],4)) + "\t")
				else : 
					vector_file.write("NA" + "\t")
		vector_file.write("\n")
		vector_file.close()

	else :
		
		vector_file = open(output_dir + "/vector_data/di_vector_frequencies.txt", "w")
		vector_file.write("Header" + "\t" + "\t".join(total_list_dicodon) + "\n")
		vector_file.write(str(header) + "\t")

		for first_codon in sequential_order_codon :
			for second_codon in sequential_order_codon :
				if first_codon in di_cu_dict : 
					if second_codon in di_cu_dict[first_codon] :
						vector_file.write(str(round(di_cu_dict[first_codon][second_codon][1],4)) + "\t")
				else : 
					vector_file.write("NA" + "\t")

		vector_file.write("\n")
		vector_file.close()

	if os.path.isfile(output_dir + "/vector_data/vector_RSCU.txt") :

		vector_file = open(output_dir + "/vector_data/vector_RSCU.txt", "a")
		vector_file.write(str(header) + "\t")

		for amino_acid in sequential_order_amino_acid :
			if amino_acid in RSCU_data :
				for codon in sequential_order_codon :
					if codon in RSCU_data[amino_acid] :
						vector_file.write(str(round(RSCU_data[amino_acid][codon],3)) + "\t")
		vector_file.write("\n")
		vector_file.close()

	else :
		
		vector_file = open(output_dir + "/vector_data/vector_RSCU.txt", "w")
		vector_file.write("Header" + "\t" + "AUG" + "\t" + "UGG" + "\t" + "AAU" + "\t" + "AAC" + "\t" + "GAU" + "\t" + "GAC" + "\t" + "UGU" + "\t" + "UGC" + "\t" + "CAA" + "\t" + "CAG" + "\t" + "GAA" + "\t" + "GAG" + "\t" + "CAU" + "\t" + "CAC" + "\t" + "AAA" + "\t" + "AAG" + "\t" + "UUU" + "\t" + "UUC" + "\t" + "UAU" + "\t" + "UAC" + "\t" + "AUU" + "\t" + "AUC" + "\t" + "AUA" + "\t" + "GCU" + "\t" + "GCC" + "\t" + "GCA" + "\t" + "GCG" + "\t" + "GGU" + "\t" + "GGC" + "\t" + "GGA" + "\t" + "GGG" + "\t" + "CCU" + "\t" + "CCC" + "\t" + "CCA" + "\t" + "CCG" + "\t" + "ACU" + "\t" + "ACC" + "\t" + "ACA" + "\t" + "ACG" + "\t" + "GUU" + "\t" + "GUC" + "\t" + "GUA" + "\t" + "GUG" + "\t" + "CGU" + "\t" + "CGC" + "\t" + "CGA" + "\t" + "CGG" + "\t" + "AGA" + "\t" + "AGG" + "\t" + "UUA" + "\t" + "UUG" + "\t" + "CUU" + "\t" + "CUC" + "\t" + "CUA" + "\t" + "CUG" + "\t" + "UCU" + "\t" + "UCC" + "\t" + "UCA" + "\t" + "UCG" + "\t" + "AGU" + "\t" + "AGC" + "\t" + "UAG" + "\t" + "UAA" + "\t" + "UGA" + "\n")
		vector_file.write(str(header) + "\t")

		for amino_acid in sequential_order_amino_acid :
			if amino_acid in RSCU_data :
				for codon in sequential_order_codon :
					if codon in RSCU_data[amino_acid] :
						vector_file.write(str(round(RSCU_data[amino_acid][codon],3)) + "\t")

		vector_file.write("\n")
		vector_file.close()



def calculate_input_values(input_file_seq_type, genetic_code_ref, optimal_codons_list, perform_COUSIN_with_boundaries, perform_vector_analysis, output_dir) :

	global input_seq_value_CAI_cod
	global input_seq_value_CAI_aa
	global input_seq_value_COUSIN_cod
	global input_seq_value_COUSIN
	global input_seq_value_ENC
	global results_for_file
	global sequential_order_results
	global results_for_file_simulation
	global sequential_order_results_simulation

	freq_aa_ref = {}

	get_freq_aa(genetic_code_ref, freq_aa_ref)
	calc_freq_equ_muta_bias(genetic_code_ref)

	sequential_order_results.append("firstline")
	results_for_file["firstline"] = "Header" + "\t" + "length" + "\t" + "%A(3)" + "\t" + "%T(3)" + "\t" + "%G(3)" + "\t" + "%C(3)" + "\t" + "%GC(all)" + "\t" + "%GC(1)" + "\t" + "%GC(2)" + "\t" + "%GC(3)" + "\t" + "COUSIN_18" + "\t" + "COUSIN_59" + "\t" + "CAI_59" + "\t" + "CAI_18" + "\t" + "ENC" + "\t" + "SCUO" + "\t" + "FOP" + "\t" + "CBI" + "\t" + "ICDI" + "\t" + "Scaled_Chi" + "\t" + "GRAVY" + "\t" + "AROMA" + "\t" + "eucl_dist_vs_ref (a.a)" + "\t" + "eucl_dist_vs_H0 (a.a)"
	if input_file_seq_type == "DNA" :
		genetic_code_ref_nucl_comp = generate_table_nucl_comp (genetic_code_ref, output_dir)
		deviation_score = calc_deviation_scores(genetic_code_ref)
		deviation_score_nucl_comp = calc_deviation_scores(genetic_code_ref_nucl_comp)
		w_scores = calc_adaptative_values (genetic_code_ref)


		for header_seq in content_seq.copy() :

			occ_aa = {}
			freq_aa_req = {}
			ATGC_percent = get_freq_ATGC_seq(content_seq[header_seq])
			GC_percent = get_freq_GC_seq(content_seq[header_seq])
			seq_DNA_to_RNA = DNA_to_RNA(content_seq[header_seq])
			dicodon_dict = create_dicodon_dict(genetic_code_ref)
			seq_dict = seq_to_dict(seq_DNA_to_RNA, header_seq, genetic_code_ref, output_dir)
			dicodon_dict_seq = seq_to_dicodon_dict(seq_DNA_to_RNA, header_seq, dicodon_dict, output_dir)

			if seq_dict != 'deleted' :

				aa_in_req = check_if_aa_in_req(seq_dict)
				deviation_score_aa = calc_deviation_scores_aa(genetic_code_ref, aa_in_req)
				get_freq_aa(seq_dict, freq_aa_req)
				RSCU_data = calc_RSCU(seq_dict)
				if perform_vector_analysis == "Yes" or perform_vector_analysis == "yes"  or perform_vector_analysis == "YES" or perform_vector_analysis == "Y" or perform_vector_analysis == "y": 
					create_data_vector(seq_dict, dicodon_dict_seq, RSCU_data, header_seq[:50], output_dir)

				input_seq_value_CAI_cod[header_seq] = calc_CAI_59(seq_dict, genetic_code_ref, w_scores)
				input_seq_value_CAI_aa[header_seq] = calc_CAI_18(seq_dict, genetic_code_ref, w_scores)
				input_seq_value_COUSIN[header_seq] = calc_COUSIN_18(seq_dict, genetic_code_ref, deviation_score, aa_in_req, perform_COUSIN_with_boundaries)
				input_seq_value_COUSIN_pond[header_seq] = calc_COUSIN_59(seq_dict, genetic_code_ref, deviation_score, deviation_score_aa, aa_in_req, perform_COUSIN_with_boundaries)
				input_seq_value_ENC[header_seq] = calc_ENC(seq_dict)
				input_seq_value_SCUO[header_seq] = calc_SCUO(seq_dict, genetic_code_ref)
				input_seq_value_FOP[header_seq] = calc_FOP(seq_dict, optimal_codons_list)
				input_seq_value_CBI[header_seq] = calc_CBI(seq_dict, optimal_codons_list)
				input_seq_value_ICDI[header_seq] = calc_ICDI(seq_dict, RSCU_data, aa_in_req)
				input_seq_value_scaled_chi[header_seq] = calc_scaled_chi(seq_dict)
				input_seq_value_GRAVY[header_seq] = calc_GRAVY(seq_dict)
				input_seq_value_AROMA[header_seq] = calc_AROMA(seq_dict)
				input_seq_value_eucl_dist_aa_que_vs_ref[header_seq] = calc_eucl_dist_que_vs_ref(freq_aa_req, freq_aa_ref, aa_in_req)
				input_seq_value_eucl_dist_aa_que_vs_H0[header_seq] = calc_eucl_dist_que_vs_H0(freq_aa_req, aa_in_req)



				set_of_seq["CAI_codon"].append(input_seq_value_CAI_cod[header_seq])
				set_of_seq["CAI_cod_famil"].append(input_seq_value_CAI_aa[header_seq])
				set_of_seq["cousin_value"].append(input_seq_value_COUSIN[header_seq])
				set_of_seq["cousin_value_pond"].append(input_seq_value_COUSIN_pond[header_seq])
				set_of_seq["ENC_value"].append(input_seq_value_ENC[header_seq])
				set_of_seq["SCUO_value"].append(input_seq_value_SCUO[header_seq])
				set_of_seq["FOP_value"].append(input_seq_value_FOP[header_seq])
				set_of_seq["CBI_value"].append(input_seq_value_CBI[header_seq])
				set_of_seq["ICDI_value"].append(input_seq_value_ICDI[header_seq])
				set_of_seq["scaled_chi_value"].append(input_seq_value_scaled_chi[header_seq])
				set_of_seq["GRAVY_value"].append(input_seq_value_GRAVY[header_seq])
				set_of_seq["AROMA_value"].append(input_seq_value_AROMA[header_seq])
				set_of_seq["eucl_dist_aa_value_que_vs_ref"].append(input_seq_value_eucl_dist_aa_que_vs_ref[header_seq])
				set_of_seq["eucl_dist_aa_value_que_vs_H0"].append(input_seq_value_eucl_dist_aa_que_vs_H0[header_seq])


				sequential_order_results.append(header_seq)
				if header_seq not in results_for_file :
					results_for_file[header_seq] = header_seq + "\t" + str(len(content_seq[header_seq])) + "\t" + str(round(ATGC_percent[0],3)) + "\t" + str(round(ATGC_percent[1],3)) + "\t" + str(round(ATGC_percent[2],3)) + "\t" + str(round(ATGC_percent[3],3)) + "\t" + str(round(GC_percent[0],3)) + "\t" + str(round(GC_percent[1],3)) + "\t" + str(round(GC_percent[2],3)) + "\t" + str(round(GC_percent[3],3)) + "\t" + str(round(input_seq_value_COUSIN[header_seq],3))  + "\t" + str(round(input_seq_value_COUSIN_pond[header_seq],3)) + "\t" + str(round(input_seq_value_CAI_cod[header_seq],3))  + "\t" + str(round(input_seq_value_CAI_aa[header_seq],3)) + "\t" + str(round(input_seq_value_ENC[header_seq],3)) + "\t" + str(round(input_seq_value_SCUO[header_seq],3)) + "\t" + str(round(input_seq_value_FOP[header_seq],3)) + "\t"+ str(round(input_seq_value_CBI[header_seq],3)) + "\t"+ str(round(input_seq_value_ICDI[header_seq],3)) +"\t"+ str(round(input_seq_value_scaled_chi[header_seq],3)) + "\t" + str(round(input_seq_value_GRAVY[header_seq],3)) + "\t" + str(round(input_seq_value_AROMA[header_seq],3)) + "\t" + str(round(input_seq_value_eucl_dist_aa_que_vs_ref[header_seq],3)) + "\t" + str(round(input_seq_value_eucl_dist_aa_que_vs_H0[header_seq],3))
				else :
					results_for_file[header_seq] += header_seq + "\t" + str(len(content_seq[header_seq])) + "\t" + str(round(ATGC_percent[0],3)) + "\t" + str(round(ATGC_percent[1],3)) + "\t" + str(round(ATGC_percent[2],3)) + "\t" + str(round(ATGC_percent[3],3)) + "\t" + str(round(GC_percent[0],3)) + "\t" + str(round(GC_percent[1],3)) + "\t" + str(round(GC_percent[2],3)) + "\t" + str(round(GC_percent[3],3)) + "\t" + str(round(input_seq_value_COUSIN[header_seq],3))  + "\t" + str(round(input_seq_value_COUSIN_pond[header_seq],3))  + "\t" + str(round(input_seq_value_CAI_cod[header_seq],3))  + "\t" + str(round(input_seq_value_CAI_aa[header_seq],3)) + "\t" + str(round(input_seq_value_ENC[header_seq],3)) + "\t" + str(round(input_seq_value_SCUO[header_seq],3)) + "\t" + str(round(input_seq_value_FOP[header_seq],3)) + "\t"+ str(round(input_seq_value_CBI[header_seq],3)) + "\t"+ str(round(input_seq_value_ICDI[header_seq],3)) + "\t"+ str(round(input_seq_value_scaled_chi[header_seq],3)) + "\t" + str(round(input_seq_value_GRAVY[header_seq],3)) + "\t" + str(round(input_seq_value_AROMA[header_seq],3)) + "\t" + str(round(input_seq_value_eucl_dist_aa_que_vs_ref[header_seq],3)) + "\t" + str(round(input_seq_value_eucl_dist_aa_que_vs_H0[header_seq],3))

	elif input_file_seq_type == "AA" :
		print ("## Note : Since your sequences seem to be in an amino acid format, COUSIN can't perform the Codon Usage analysis ##")
	
	results_for_file_simulation = results_for_file.copy()
	sequential_order_results_simulation = sequential_order_results



file_name = ""
absc_name = ""
ord_name = ""
title_name = ""
seq_name = ""
csv_name = ""
header_name = ""
step_title = ""
firstline = True


def treat_simulation_data(output_dir) :

	### CREATION OF DENSITY CURVES
	global file_name
	global absc_name
	global ord_name
	global title_name
	global csv_name
	global firstline
	global step_title

	step_title = "Graph on simulated sequences (basic simulation)"
	graph_dir = "basic_simulation_graphs"

	list_size = [100, 200, 150]

	i = 0
	x = 1
	firstline = True
	legends_names = []
	legends_length = []
	size_tmp = 0
	for size in list_size :
		file_name = "simulation_dens_COUSIN_18"
		absc_name = "$\mathrm{COUSIN}_{18}\/\/\mathrm{score}$"
		ord_name = "density score "
		title_name = "$\mathrm{COUSIN}_{18}\/\/\mathrm{SCORE}$"
		csv_name = "COUSIN 18"
		type_of_data = "a.a"

		matrice[size] = {}
		quantile_99_left_x[size] = {}
		quantile_99_right_x[size] = {}
		quantile_95_left_x[size] = {}
		quantile_95_right_x[size] = {}
		quantile_99_left_y[size] = {}
		quantile_99_right_y[size] = {}
		quantile_95_left_y[size] = {}
		quantile_95_right_y[size] = {}

		size_tmp = size_tmp + size
		csv_name = "COUSIN_18 (" + str(size_tmp) + " a.a)"
		graph_name = "step " + str(x)
		x += 1
		legends_names.append(graph_name)
		legends_length.append(size_tmp)

		for header_seq in content_seq :
			get_matrice_and_quantiles(input_seq_value_COUSIN[header_seq], COUSIN_calc[size],size, header_seq, results_for_file, output_dir)
		firstline = True

	create_graphs(output_dir, graph_dir, None, legends_names, legends_length, type_of_data, matrice[100], matrice[200], matrice[150] )

	i = 0
	x = 1
	firstline = True
	legends_names = []
	legends_length = []
	size_tmp = 0
	for size in list_size :
		file_name = "simulation_dens_COUSIN_59"
		absc_name = "$\mathrm{COUSIN}_{59}\/\/\mathrm{score}$"
		ord_name = "density score "
		title_name = "$\mathrm{COUSIN}_{59}\/\/\mathrm{SCORE}$"
		csv_name = "COUSIN_59"
		type_of_data = "a.a"

		matrice[size] = {}
		quantile_99_left_x[size] = {}
		quantile_99_right_x[size] = {}
		quantile_95_left_x[size] = {}
		quantile_95_right_x[size] = {}
		quantile_99_left_y[size] = {}
		quantile_99_right_y[size] = {}
		quantile_95_left_y[size] = {}
		quantile_95_right_y[size] = {}

		size_tmp = size_tmp + size
		csv_name = "COUSIN_59 (" + str(size_tmp) + " a.a)"
		graph_name = "step " + str(x)
		x += 1
		legends_names.append(graph_name)
		legends_length.append(size_tmp)
		#
		for header_seq in content_seq :
			get_matrice_and_quantiles(input_seq_value_COUSIN_pond[header_seq], COUSIN_calc_pond[size],size, header_seq, results_for_file, output_dir)
		firstline = True

	create_graphs(output_dir, graph_dir, None, legends_names, legends_length, type_of_data, matrice[100], matrice[200], matrice[150] )

	i = 0
	x = 1
	firstline = True
	legends_names = []
	legends_length = []
	size_tmp = 0
	for size in list_size :
		file_name = "simulation_dens_CAI_59"
		absc_name = "$\mathrm{CAI}_{59}\/\/\mathrm{score}$"
		ord_name = "density score  "
		title_name = "$\mathrm{CAI}_{59}\/\/\mathrm{SCORE}$"
		csv_name = "CAI_59"
		type_of_data = "a.a"

		matrice[size] = {}
		quantile_99_left_x[size] = {}
		quantile_99_right_x[size] = {}
		quantile_95_left_x[size] = {}
		quantile_95_right_x[size] = {}
		quantile_99_left_y[size] = {}
		quantile_99_right_y[size] = {}
		quantile_95_left_y[size] = {}
		quantile_95_right_y[size] = {}

		size_tmp = size_tmp + size
		csv_name = "CAI_59 (" + str(size_tmp) + " a.a)"
		graph_name = "step " + str(x)
		x += 1
		legends_names.append(graph_name)
		legends_length.append(size_tmp)

		for header_seq in content_seq :
			get_matrice_and_quantiles(input_seq_value_CAI_cod[header_seq], CAI_codon_calc[size],size, header_seq, results_for_file, output_dir)
		firstline = True

	create_graphs(output_dir, graph_dir, None, legends_names, legends_length, type_of_data, matrice[100], matrice[200], matrice[150])

	i = 0
	x = 1
	firstline = True
	legends_names = []
	legends_length = []
	size_tmp = 0
	for size in list_size :
		file_name = "simulation_dens_CAI_18"
		absc_name = "$\mathrm{CAI}_{18}\/\/\mathrm{score}$"
		ord_name = "density score  "
		title_name = "$\mathrm{CAI}_{18}\/\/\mathrm{SCORE}$"
		csv_name = "CAI_18"
		type_of_data = "a.a"

		matrice[size] = {}
		quantile_99_left_x[size] = {}
		quantile_99_right_x[size] = {}
		quantile_95_left_x[size] = {}
		quantile_95_right_x[size] = {}
		quantile_99_left_y[size] = {}
		quantile_99_right_y[size] = {}
		quantile_95_left_y[size] = {}
		quantile_95_right_y[size] = {}

		size_tmp = size_tmp + size
		csv_name = "CAI_18 (" + str(size_tmp) + " a.a)"
		graph_name = "step " + str(x)
		x += 1
		legends_names.append(graph_name)
		legends_length.append(size_tmp)

		for header_seq in content_seq :
			get_matrice_and_quantiles(input_seq_value_CAI_aa[header_seq], CAI_aa_calc[size],size, header_seq, results_for_file, output_dir)
		firstline = True

	create_graphs(output_dir, graph_dir, None, legends_names, legends_length, type_of_data, matrice[100], matrice[200], matrice[150])

	i = 0
	x = 1
	firstline = True
	legends_names = []
	legends_length = []
	size_tmp = 0
	for size in list_size :
		file_name = "simulation_dens_ENC"
		absc_name = "ENC value "
		ord_name = "density score  "
		title_name = "Density curves for ENC measure"
		csv_name = "ENC"
		type_of_data = "a.a"

		matrice[size] = {}
		quantile_99_left_x[size] = {}
		quantile_99_right_x[size] = {}
		quantile_95_left_x[size] = {}
		quantile_95_right_x[size] = {}
		quantile_99_left_y[size] = {}
		quantile_99_right_y[size] = {}
		quantile_95_left_y[size] = {}
		quantile_95_right_y[size] = {}

		size_tmp = size_tmp + size
		csv_name = "ENC (" + str(size_tmp) + " a.a)"
		graph_name = "step " + str(x)
		x += 1
		legends_names.append(csv_name)
		legends_length.append(size_tmp)

		for header_seq in content_seq :
			get_matrice_and_quantiles(input_seq_value_ENC[header_seq], ENC_calc[size],size, header_seq, results_for_file, output_dir)
		firstline = True

	create_graphs(output_dir, graph_dir, None, legends_names, legends_length, type_of_data, matrice[100], matrice[200], matrice[150] )

	i = 0
	x = 1
	firstline = True
	legends_names = []
	legends_length = []
	size_tmp = 0
	for size in list_size :
		file_name = "simulation_dens_SCUO"
		absc_name = "SCUO value "
		ord_name = "density score  "
		title_name = "Density curves for SCUO measure"
		csv_name = "SCUO"
		type_of_data = "a.a"

		matrice[size] = {}
		quantile_99_left_x[size] = {}
		quantile_99_right_x[size] = {}
		quantile_95_left_x[size] = {}
		quantile_95_right_x[size] = {}
		quantile_99_left_y[size] = {}
		quantile_99_right_y[size] = {}
		quantile_95_left_y[size] = {}
		quantile_95_right_y[size] = {}

		size_tmp = size_tmp + size
		csv_name = "SCUO (" + str(size_tmp) + " a.a)"
		graph_name = "step " + str(x)
		x += 1
		legends_names.append(csv_name)
		legends_length.append(size_tmp)

		for header_seq in content_seq :
			get_matrice_and_quantiles(input_seq_value_SCUO[header_seq], SCUO_calc[size],size, header_seq, results_for_file, output_dir)
		firstline = True

	create_graphs(output_dir, graph_dir, None, legends_names, legends_length, type_of_data, matrice[100], matrice[200], matrice[150] )

	i = 0
	x = 1
	firstline = True
	legends_names = []
	legends_length = []
	size_tmp = 0
	for size in list_size :
		file_name = "simulation_dens_FOP"
		absc_name = "FOP value "
		ord_name = "density score  "
		title_name = "Density curves for FOP measure"
		csv_name = "FOP"
		type_of_data = "a.a"

		matrice[size] = {}
		quantile_99_left_x[size] = {}
		quantile_99_right_x[size] = {}
		quantile_95_left_x[size] = {}
		quantile_95_right_x[size] = {}
		quantile_99_left_y[size] = {}
		quantile_99_right_y[size] = {}
		quantile_95_left_y[size] = {}
		quantile_95_right_y[size] = {}

		size_tmp = size_tmp + size
		csv_name = "FOP (" + str(size_tmp) + " a.a)"
		graph_name = "step " + str(x)
		x += 1
		legends_names.append(csv_name)
		legends_length.append(size_tmp)

		for header_seq in content_seq :
			get_matrice_and_quantiles(input_seq_value_FOP[header_seq], FOP_calc[size],size, header_seq, results_for_file, output_dir)
		firstline = True

	create_graphs(output_dir, graph_dir, None, legends_names, legends_length, type_of_data, matrice[100], matrice[200], matrice[150] )

	i = 0
	x = 1
	firstline = True
	legends_names = []
	legends_length = []
	size_tmp = 0
	for size in list_size :
		file_name = "simulation_dens_CBI"
		absc_name = "CBI value "
		ord_name = "density score  "
		title_name = "Density curves for CBI measure"
		csv_name = "CBI"
		type_of_data = "a.a"

		matrice[size] = {}
		quantile_99_left_x[size] = {}
		quantile_99_right_x[size] = {}
		quantile_95_left_x[size] = {}
		quantile_95_right_x[size] = {}
		quantile_99_left_y[size] = {}
		quantile_99_right_y[size] = {}
		quantile_95_left_y[size] = {}
		quantile_95_right_y[size] = {}

		size_tmp = size_tmp + size
		csv_name = "CBI (" + str(size_tmp) + " a.a)"
		graph_name = "step " + str(x)
		x += 1
		legends_names.append(csv_name)
		legends_length.append(size_tmp)

		for header_seq in content_seq :
			get_matrice_and_quantiles(input_seq_value_CBI[header_seq], CBI_calc[size],size, header_seq, results_for_file, output_dir)
		firstline = True

	create_graphs(output_dir, graph_dir, None, legends_names, legends_length, type_of_data, matrice[100], matrice[200], matrice[150] )

	i = 0
	x = 1
	firstline = True
	legends_names = []
	legends_length = []
	size_tmp = 0
	for size in list_size :
		file_name = "simulation_dens_ICDI"
		absc_name = "ICDI value "
		ord_name = "density score  "
		title_name = "Density curves for ICDI measure"
		csv_name = "ICDI"
		type_of_data = "a.a"

		matrice[size] = {}
		quantile_99_left_x[size] = {}
		quantile_99_right_x[size] = {}
		quantile_95_left_x[size] = {}
		quantile_95_right_x[size] = {}
		quantile_99_left_y[size] = {}
		quantile_99_right_y[size] = {}
		quantile_95_left_y[size] = {}
		quantile_95_right_y[size] = {}

		size_tmp = size_tmp + size
		csv_name = "ICDI (" + str(size_tmp) + " a.a)"
		graph_name = "step " + str(x)
		x += 1
		legends_names.append(csv_name)
		legends_length.append(size_tmp)

		for header_seq in content_seq :
			get_matrice_and_quantiles(input_seq_value_ICDI[header_seq], ICDI_calc[size],size, header_seq, results_for_file, output_dir)
		firstline = True

	create_graphs(output_dir, graph_dir, None, legends_names, legends_length, type_of_data, matrice[100], matrice[200], matrice[150] )

	i = 0
	x = 1
	firstline = True
	legends_names = []
	legends_length = []
	size_tmp = 0
	for size in list_size :
		file_name = "simulation_dens_scaled_chi"
		absc_name = "scaled chi value "
		ord_name = "density score  "
		title_name = "Density curves for scaled chi measure"
		csv_name = "scaled_chi"
		type_of_data = "a.a"

		matrice[size] = {}
		quantile_99_left_x[size] = {}
		quantile_99_right_x[size] = {}
		quantile_95_left_x[size] = {}
		quantile_95_right_x[size] = {}
		quantile_99_left_y[size] = {}
		quantile_99_right_y[size] = {}
		quantile_95_left_y[size] = {}
		quantile_95_right_y[size] = {}

		size_tmp = size_tmp + size
		csv_name = "scaled_chi (" + str(size_tmp) + " a.a)"
		graph_name = "step " + str(x)
		x += 1
		legends_names.append(csv_name)
		legends_length.append(size_tmp)

		for header_seq in content_seq :
			get_matrice_and_quantiles(input_seq_value_scaled_chi[header_seq], scaled_chi_calc[size],size, header_seq, results_for_file, output_dir)
		firstline = True

	create_graphs(output_dir, graph_dir, None, legends_names, legends_length, type_of_data, matrice[100], matrice[200], matrice[150] )


def treat_input_simulation_analysis_data(output_dir) :

	### CREATION OF DENSITY CURVES
	global file_name
	global absc_name
	global ord_name
	global title_name
	global csv_name
	global firstline
	global header_name
	global step_title

	graph_dir = "simulation_graphs"

	i = 0
	firstline = True
	legends_names = []
	legends_length = {}
	for header_seq in content_seq :

		file_name = str(header_seq)[:210] + "_dens_COUSIN_18"
		absc_name = "$\mathrm{COUSIN}_{18}\/\/\mathrm{score}$"
		ord_name = "density score "
		title_name = "$\mathrm{COUSIN}_{18}\/\/\mathrm{score}$"
		header_name = header_seq
		step_title = "Graph on simulated data from the sequence " + str(header_seq)[:20]
		type_of_data = "a.a"

		header_seq_list = []
		j = 0
		while j < 2 :

			if header_seq not in legends_length :
				legends_length[header_seq] = []
			legends_length[header_seq].append(len(content_seq[header_seq]))
			tmp_header_seq_list = header_seq + "_" + str(j)
			header_seq_list.append(tmp_header_seq_list)

			matrice[header_seq_list[j]] = {}
			quantile_99_left_x[header_seq_list[j]] = {}
			quantile_99_right_x[header_seq_list[j]] = {}
			quantile_95_left_x[header_seq_list[j]] = {}
			quantile_95_right_x[header_seq_list[j]] = {}
			quantile_99_left_y[header_seq_list[j]] = {}
			quantile_99_right_y[header_seq_list[j]] = {}
			quantile_95_left_y[header_seq_list[j]] = {}
			quantile_95_right_y[header_seq_list[j]] = {}
			density_result[header_seq_list[j]] = {}

			if i == 0 :
				firstline = True

			if j == 0 :
				csv_name = "COUSIN_18 (simu_analysis)"
				graph_name = "same a.a content"
				legends_names.append(graph_name)
				get_matrice_and_quantiles(input_seq_value_COUSIN[header_seq], input_sim_COUSIN_calc[header_seq],header_seq_list[j], header_seq, results_for_file_simulation, output_dir)
			else :
				csv_name = "COUSIN_18 (a.a random guided) (simu_analysis)"
				graph_name = "a.a random-guided"
				legends_names.append(graph_name)
				get_matrice_and_quantiles(input_seq_value_COUSIN[header_seq], input_sim_COUSIN_calc_rand_aa[header_seq],header_seq_list[j], header_seq, results_for_file_simulation, output_dir)

			j += 1
		if i == 0 :
			firstline = True
		csv_name = "COUSIN_18 U_test"
		output_dir_wilcox = "simulation_results/Wilcoxon_Mann_Whitney_results"
		perform_wilcoxon_mann_whitney_u_test_simulation_step(input_sim_COUSIN_calc[header_seq], input_sim_COUSIN_calc_rand_aa[header_seq], csv_name, header_seq, results_for_file_simulation, output_dir, output_dir_wilcox)

		if i < 10 :
			create_graphs(output_dir, graph_dir, input_seq_value_COUSIN[header_seq], legends_names, legends_length[header_seq], type_of_data, matrice[header_seq_list[0]], matrice[header_seq_list[1]])
			i+= 1

	i = 0
	firstline = True
	legends_names = []
	legends_length = {}
	for header_seq in content_seq :

		file_name = str(header_seq)[:210] + "_dens_COUSIN_59"
		absc_name = "$\mathrm{COUSIN}_{59}\/\/\mathrm{score}$"
		ord_name = "density score "
		title_name = "$\mathrm{COUSIN}_{59}\/\/\mathrm{score}$"
		header_name = header_seq
		step_title = "Graph on simulated data from the sequence " + str(header_seq)[:20]
		type_of_data = "a.a"


		header_seq_list = []
		j = 0
		while j < 2 :

			if header_seq not in legends_length :
				legends_length[header_seq] = []
			legends_length[header_seq].append(len(content_seq[header_seq]))
			tmp_header_seq_list = header_seq + "_" + str(j)
			header_seq_list.append(tmp_header_seq_list)

			matrice[header_seq_list[j]] = {}
			quantile_99_left_x[header_seq_list[j]] = {}
			quantile_99_right_x[header_seq_list[j]] = {}
			quantile_95_left_x[header_seq_list[j]] = {}
			quantile_95_right_x[header_seq_list[j]] = {}
			quantile_99_left_y[header_seq_list[j]] = {}
			quantile_99_right_y[header_seq_list[j]] = {}
			quantile_95_left_y[header_seq_list[j]] = {}
			quantile_95_right_y[header_seq_list[j]] = {}
			density_result[header_seq_list[j]] = {}

			if i == 0 :
				firstline = True

			if j == 0 :
				csv_name = "COUSIN_59 (simu_analysis)"
				graph_name = "same a.a content"
				legends_names.append(graph_name)
				get_matrice_and_quantiles(input_seq_value_COUSIN_pond[header_seq], input_sim_COUSIN_pond_calc[header_seq],header_seq_list[j], header_seq, results_for_file_simulation, output_dir)
			else :
				csv_name = "COUSIN_59 (a.a random guided) (simu_analysis)"
				graph_name = "a.a random-guided"
				legends_names.append(graph_name)
				get_matrice_and_quantiles(input_seq_value_COUSIN_pond[header_seq], input_sim_COUSIN_pond_calc_rand_aa[header_seq],header_seq_list[j], header_seq, results_for_file_simulation, output_dir)

			j += 1
		if i == 0 :
			firstline = True
		csv_name = "COUSIN_59 U_test"
		output_dir_wilcox = "simulation_results/Wilcoxon_Mann_Whitney_results"
		perform_wilcoxon_mann_whitney_u_test_simulation_step(input_sim_COUSIN_pond_calc[header_seq], input_sim_COUSIN_pond_calc_rand_aa[header_seq], csv_name, header_seq, results_for_file_simulation, output_dir, output_dir_wilcox)
	
		if i < 10 :
			create_graphs(output_dir, graph_dir, input_seq_value_COUSIN_pond[header_seq], legends_names, legends_length[header_seq], type_of_data, matrice[header_seq_list[0]], matrice[header_seq_list[1]])
			i+=1

	i = 0
	firstline = True
	legends_names = []
	legends_length = {}
	for header_seq in content_seq :

		file_name = str(header_seq)[:210] + "_dens_CAI_59"
		absc_name = "$\mathrm{CAI}_{59}\/\/\mathrm{score}$"
		ord_name = "density score "
		title_name = "$\mathrm{CAI}_{59}\/\/\mathrm{score}$"
		header_name = header_seq
		step_title = "Graph on simulated data from the sequence " + str(header_seq)[:20]
		type_of_data = "a.a"

		header_seq_list = []
		j = 0
		while j < 2 :
			if header_seq not in legends_length :
				legends_length[header_seq] = []
			legends_length[header_seq].append(len(content_seq[header_seq]))
			tmp_header_seq_list = header_seq + "_" + str(j)
			header_seq_list.append(tmp_header_seq_list)

			matrice[header_seq_list[j]] = {}
			quantile_99_left_x[header_seq_list[j]] = {}
			quantile_99_right_x[header_seq_list[j]] = {}
			quantile_95_left_x[header_seq_list[j]] = {}
			quantile_95_right_x[header_seq_list[j]] = {}
			quantile_99_left_y[header_seq_list[j]] = {}
			quantile_99_right_y[header_seq_list[j]] = {}
			quantile_95_left_y[header_seq_list[j]] = {}
			quantile_95_right_y[header_seq_list[j]] = {}
			density_result[header_seq_list[j]] = {}

			if i == 0 :
				firstline = True

			if j == 0 :
				csv_name = "CAI_59 (simu_analysis)"
				graph_name = "same a.a content"
				legends_names.append(graph_name)
				get_matrice_and_quantiles(input_seq_value_CAI_cod[header_seq], input_sim_CAI_codon_calc[header_seq],header_seq_list[j], header_seq, results_for_file_simulation, output_dir)
			else :
				csv_name = "CAI_59 (a.a random guided) (simu_analysis)"
				graph_name = "a.a random-guided"
				legends_names.append(graph_name)
				get_matrice_and_quantiles(input_seq_value_CAI_cod[header_seq], input_sim_CAI_codon_calc_rand_aa[header_seq],header_seq_list[j], header_seq, results_for_file_simulation, output_dir)

			j += 1
		if i == 0 :
			firstline = True
		csv_name = "CAI_59 U_test"
		output_dir_wilcox = "simulation_results/Wilcoxon_Mann_Whitney_results"
		perform_wilcoxon_mann_whitney_u_test_simulation_step(input_sim_CAI_codon_calc[header_seq], input_sim_CAI_codon_calc_rand_aa[header_seq], csv_name, header_seq, results_for_file_simulation, output_dir, output_dir_wilcox)

		if i < 10 :
			create_graphs(output_dir, graph_dir, input_seq_value_CAI_cod[header_seq], legends_names, legends_length[header_seq], type_of_data, matrice[header_seq_list[0]], matrice[header_seq_list[1]])
			i+= 1

	i = 0
	firstline = True
	legends_names = []
	legends_length = {}
	for header_seq in content_seq :

		file_name = str(header_seq)[:210] + "_dens_CAI_18"
		absc_name = "$\mathrm{CAI}_{18}\/\/\mathrm{score}$"
		ord_name = "density score "
		title_name = "$\mathrm{CAI}_{18}\/\/\mathrm{score}$"
		header_name = header_seq
		step_title = "Graph on simulated data from the sequence " + str(header_seq)[:20]
		type_of_data = "a.a"

		header_seq_list = []
		j = 0
		while j < 2 :
			if header_seq not in legends_length :
				legends_length[header_seq] = []
			legends_length[header_seq].append(len(content_seq[header_seq]))
			tmp_header_seq_list = header_seq + "_" + str(j)
			header_seq_list.append(tmp_header_seq_list)

			matrice[header_seq_list[j]] = {}
			quantile_99_left_x[header_seq_list[j]] = {}
			quantile_99_right_x[header_seq_list[j]] = {}
			quantile_95_left_x[header_seq_list[j]] = {}
			quantile_95_right_x[header_seq_list[j]] = {}
			quantile_99_left_y[header_seq_list[j]] = {}
			quantile_99_right_y[header_seq_list[j]] = {}
			quantile_95_left_y[header_seq_list[j]] = {}
			quantile_95_right_y[header_seq_list[j]] = {}
			density_result[header_seq_list[j]] = {}

			if i == 0 :
				firstline = True

			if j == 0 :
				csv_name = "CAI_18 (simu_analysis)"
				graph_name = "same a.a content"
				legends_names.append(graph_name)
				get_matrice_and_quantiles(input_seq_value_CAI_aa[header_seq], input_sim_CAI_aa_calc[header_seq],header_seq_list[j], header_seq, results_for_file_simulation, output_dir)
			else :
				csv_name = "CAI_18 (a.a random guided) (simu_analysis)"
				graph_name = "a.a random-guided"
				legends_names.append(graph_name)
				get_matrice_and_quantiles(input_seq_value_CAI_aa[header_seq], input_sim_CAI_aa_calc_rand_aa[header_seq],header_seq_list[j], header_seq, results_for_file_simulation, output_dir)

			j += 1
		if i == 0 :
			firstline = True
		csv_name = "CAI_18 U_test"
		output_dir_wilcox = "simulation_results/Wilcoxon_Mann_Whitney_results"
		perform_wilcoxon_mann_whitney_u_test_simulation_step(input_sim_CAI_aa_calc[header_seq], input_sim_CAI_aa_calc_rand_aa[header_seq], csv_name, header_seq, results_for_file_simulation, output_dir, output_dir_wilcox)

		if i < 10 :
			create_graphs(output_dir, graph_dir, input_seq_value_CAI_aa[header_seq], legends_names, legends_length[header_seq], type_of_data, matrice[header_seq_list[0]], matrice[header_seq_list[1]])
			i+= 1

	i = 0
	firstline = True
	legends_names = []
	legends_length = {}
	for header_seq in content_seq :

		file_name = str(header_seq)[:210] + "_dens_ENC"
		absc_name = "ENC score "
		ord_name = "density score "
		title_name = "ENC SCORE"
		header_name = header_seq
		step_title = "Graph on simulated data from the sequence " + str(header_seq)[:20]
		type_of_data = "a.a"

		header_seq_list = []
		j = 0
		while j < 2 :
			if header_seq not in legends_length :
				legends_length[header_seq] = []
			legends_length[header_seq].append(len(content_seq[header_seq]))
			tmp_header_seq_list = header_seq + "_" + str(j)
			header_seq_list.append(tmp_header_seq_list)

			matrice[header_seq_list[j]] = {}
			quantile_99_left_x[header_seq_list[j]] = {}
			quantile_99_right_x[header_seq_list[j]] = {}
			quantile_95_left_x[header_seq_list[j]] = {}
			quantile_95_right_x[header_seq_list[j]] = {}
			quantile_99_left_y[header_seq_list[j]] = {}
			quantile_99_right_y[header_seq_list[j]] = {}
			quantile_95_left_y[header_seq_list[j]] = {}
			quantile_95_right_y[header_seq_list[j]] = {}
			density_result[header_seq_list[j]] = {}

			if i == 0 :
				firstline = True

			if j == 0 :
				csv_name = "ENC (simu_analysis)"
				graph_name = "same a.a content"
				legends_names.append(graph_name)
				get_matrice_and_quantiles(input_seq_value_ENC[header_seq], input_sim_ENC[header_seq],header_seq_list[j], header_seq, results_for_file_simulation, output_dir)
			else :
				csv_name = "ENC (a.a random guided) (simu_analysis)"
				graph_name = "a.a random-guided"
				legends_names.append(graph_name)
				get_matrice_and_quantiles(input_seq_value_ENC[header_seq], input_sim_ENC_rand_aa[header_seq],header_seq_list[j], header_seq, results_for_file_simulation, output_dir)

			j += 1
		if i == 0 :
			firstline = True
		csv_name = "ENC U_test"
		output_dir_wilcox = "simulation_results/Wilcoxon_Mann_Whitney_results"
		perform_wilcoxon_mann_whitney_u_test_simulation_step(input_sim_ENC[header_seq], input_sim_ENC_rand_aa[header_seq], csv_name, header_seq, results_for_file_simulation, output_dir, output_dir_wilcox)

		if i < 10 :
			create_graphs(output_dir, graph_dir, input_seq_value_ENC[header_seq], legends_names, legends_length[header_seq], type_of_data, matrice[header_seq_list[0]], matrice[header_seq_list[1]])
			i+= 1
	   
	i = 0
	firstline = True
	legends_names = []
	legends_length = {}
	for header_seq in content_seq :

		file_name = str(header_seq)[:210] + "_dens_SCUO"
		absc_name = "SCUO score "
		ord_name = "density score "
		title_name = "SCUO SCORE"
		header_name = header_seq
		step_title = "Graph on simulated data from the sequence " + str(header_seq)[:20]
		type_of_data = "a.a"

		header_seq_list = []
		j = 0
		while j < 2 :
			if header_seq not in legends_length :
				legends_length[header_seq] = []
			legends_length[header_seq].append(len(content_seq[header_seq]))
			tmp_header_seq_list = header_seq + "_" + str(j)
			header_seq_list.append(tmp_header_seq_list)

			matrice[header_seq_list[j]] = {}
			quantile_99_left_x[header_seq_list[j]] = {}
			quantile_99_right_x[header_seq_list[j]] = {}
			quantile_95_left_x[header_seq_list[j]] = {}
			quantile_95_right_x[header_seq_list[j]] = {}
			quantile_99_left_y[header_seq_list[j]] = {}
			quantile_99_right_y[header_seq_list[j]] = {}
			quantile_95_left_y[header_seq_list[j]] = {}
			quantile_95_right_y[header_seq_list[j]] = {}
			density_result[header_seq_list[j]] = {}

			if i == 0 :
				firstline = True

			if j == 0 :
				csv_name = "SCUO (simu_analysis)"
				graph_name = "same a.a content"
				legends_names.append(graph_name)
				get_matrice_and_quantiles(input_seq_value_SCUO[header_seq], input_sim_SCUO[header_seq],header_seq_list[j], header_seq, results_for_file_simulation, output_dir)
			else :
				csv_name = "SCUO (a.a random guided) (simu_analysis)"
				graph_name = "a.a random-guided"
				legends_names.append(graph_name)
				get_matrice_and_quantiles(input_seq_value_SCUO[header_seq], input_sim_SCUO_rand_aa[header_seq],header_seq_list[j], header_seq, results_for_file_simulation, output_dir)

			j += 1
		if i == 0 :
			firstline = True
		csv_name = "SCUO U_test"
		output_dir_wilcox = "simulation_results/Wilcoxon_Mann_Whitney_results"
		perform_wilcoxon_mann_whitney_u_test_simulation_step(input_sim_SCUO[header_seq], input_sim_SCUO_rand_aa[header_seq], csv_name, header_seq, results_for_file_simulation, output_dir, output_dir_wilcox)

		if i < 10 :
			create_graphs(output_dir, graph_dir, input_seq_value_SCUO[header_seq], legends_names, legends_length[header_seq], type_of_data, matrice[header_seq_list[0]], matrice[header_seq_list[1]])
			i+= 1
		
	i = 0
	firstline = True
	legends_names = []
	legends_length = {}
	for header_seq in content_seq :

		file_name = str(header_seq)[:210] + "_dens_FOP"
		absc_name = "FOP score "
		ord_name = "density score "
		title_name = "FOP SCORE"
		header_name = header_seq
		step_title = "Graph on simulated data from the sequence " + str(header_seq)[:20]
		type_of_data = "a.a"

		header_seq_list = []
		j = 0
		while j < 2 :
			if header_seq not in legends_length :
				legends_length[header_seq] = []
			legends_length[header_seq].append(len(content_seq[header_seq]))
			tmp_header_seq_list = header_seq + "_" + str(j)
			header_seq_list.append(tmp_header_seq_list)

			matrice[header_seq_list[j]] = {}
			quantile_99_left_x[header_seq_list[j]] = {}
			quantile_99_right_x[header_seq_list[j]] = {}
			quantile_95_left_x[header_seq_list[j]] = {}
			quantile_95_right_x[header_seq_list[j]] = {}
			quantile_99_left_y[header_seq_list[j]] = {}
			quantile_99_right_y[header_seq_list[j]] = {}
			quantile_95_left_y[header_seq_list[j]] = {}
			quantile_95_right_y[header_seq_list[j]] = {}
			density_result[header_seq_list[j]] = {}

			if i == 0 :
				firstline = True

			if j == 0 :
				csv_name = "FOP (simu_analysis)"
				graph_name = "same a.a content"
				legends_names.append(graph_name)
				get_matrice_and_quantiles(input_seq_value_FOP[header_seq], input_sim_FOP[header_seq],header_seq_list[j], header_seq, results_for_file_simulation, output_dir)
			else :
				csv_name = "FOP (a.a random guided) (simu_analysis)"
				graph_name = "a.a random-guided"
				legends_names.append(graph_name)
				get_matrice_and_quantiles(input_seq_value_FOP[header_seq], input_sim_FOP_rand_aa[header_seq],header_seq_list[j], header_seq, results_for_file_simulation, output_dir)

			j += 1
		if i == 0 :
			firstline = True
		csv_name = "FOP U_test"
		output_dir_wilcox = "simulation_results/Wilcoxon_Mann_Whitney_results"
		perform_wilcoxon_mann_whitney_u_test_simulation_step(input_sim_FOP[header_seq], input_sim_FOP_rand_aa[header_seq], csv_name, header_seq, results_for_file_simulation, output_dir, output_dir_wilcox)

		if i < 10 :
			create_graphs(output_dir, graph_dir, input_seq_value_FOP[header_seq], legends_names, legends_length[header_seq], type_of_data, matrice[header_seq_list[0]], matrice[header_seq_list[1]])
			i+= 1
		
	i = 0
	firstline = True
	legends_names = []
	legends_length = {}
	for header_seq in content_seq :

		file_name = str(header_seq)[:210] + "_dens_CBI"
		absc_name = "CBI score "
		ord_name = "density score "
		title_name = "CBI SCORE"
		header_name = header_seq
		step_title = "Graph on simulated data from the sequence " + str(header_seq)[:20]
		type_of_data = "a.a"

		header_seq_list = []
		j = 0
		while j < 2 :
			if header_seq not in legends_length :
				legends_length[header_seq] = []
			legends_length[header_seq].append(len(content_seq[header_seq]))
			tmp_header_seq_list = header_seq + "_" + str(j)
			header_seq_list.append(tmp_header_seq_list)

			matrice[header_seq_list[j]] = {}
			quantile_99_left_x[header_seq_list[j]] = {}
			quantile_99_right_x[header_seq_list[j]] = {}
			quantile_95_left_x[header_seq_list[j]] = {}
			quantile_95_right_x[header_seq_list[j]] = {}
			quantile_99_left_y[header_seq_list[j]] = {}
			quantile_99_right_y[header_seq_list[j]] = {}
			quantile_95_left_y[header_seq_list[j]] = {}
			quantile_95_right_y[header_seq_list[j]] = {}
			density_result[header_seq_list[j]] = {}

			if i == 0 :
				firstline = True

			if j == 0 :
				csv_name = "CBI (simu_analysis)"
				graph_name = "same a.a content"
				legends_names.append(graph_name)
				get_matrice_and_quantiles(input_seq_value_CBI[header_seq], input_sim_CBI[header_seq],header_seq_list[j], header_seq, results_for_file_simulation, output_dir)
			else :
				csv_name = "CBI (a.a random_guided) (simu_analysis)"
				graph_name = "a.a random-guided"
				legends_names.append(graph_name)
				get_matrice_and_quantiles(input_seq_value_CBI[header_seq], input_sim_CBI_rand_aa[header_seq],header_seq_list[j], header_seq, results_for_file_simulation, output_dir)

			j += 1
		if i == 0 :
			firstline = True
		csv_name = "CBI U_test"
		output_dir_wilcox = "simulation_results/Wilcoxon_Mann_Whitney_results"
		perform_wilcoxon_mann_whitney_u_test_simulation_step(input_sim_CBI[header_seq], input_sim_CBI_rand_aa[header_seq], csv_name, header_seq, results_for_file_simulation, output_dir, output_dir_wilcox)

		if i < 10 :
			create_graphs(output_dir, graph_dir, input_seq_value_CBI[header_seq], legends_names, legends_length[header_seq], type_of_data, matrice[header_seq_list[0]], matrice[header_seq_list[1]])
			i+= 1
		
	i = 0
	firstline = True
	legends_names = []
	legends_length = {}
	for header_seq in content_seq :

		file_name = str(header_seq)[:210] + "_dens_ICDI"
		absc_name = "ICDI score "
		ord_name = "density score "
		title_name = "ICDI SCORE"
		header_name = header_seq
		step_title = "Graph on simulated data from the sequence " + str(header_seq)[:20]
		type_of_data = "a.a"

		header_seq_list = []
		j = 0
		while j < 2 :
			if header_seq not in legends_length :
				legends_length[header_seq] = []
			legends_length[header_seq].append(len(content_seq[header_seq]))
			tmp_header_seq_list = header_seq + "_" + str(j)
			header_seq_list.append(tmp_header_seq_list)

			matrice[header_seq_list[j]] = {}
			quantile_99_left_x[header_seq_list[j]] = {}
			quantile_99_right_x[header_seq_list[j]] = {}
			quantile_95_left_x[header_seq_list[j]] = {}
			quantile_95_right_x[header_seq_list[j]] = {}
			quantile_99_left_y[header_seq_list[j]] = {}
			quantile_99_right_y[header_seq_list[j]] = {}
			quantile_95_left_y[header_seq_list[j]] = {}
			quantile_95_right_y[header_seq_list[j]] = {}
			density_result[header_seq_list[j]] = {}

			if i == 0 :
				firstline = True

			if j == 0 :
				csv_name = "ICDI (simu_analysis)"
				graph_name = "same a.a content"
				legends_names.append(graph_name)
				get_matrice_and_quantiles(input_seq_value_ICDI[header_seq], input_sim_ICDI[header_seq],header_seq_list[j], header_seq, results_for_file_simulation, output_dir)
			else :
				csv_name = "ICDI (a.a random guided) (simu_analysis)"
				graph_name = "a.a random-guided"
				legends_names.append(graph_name)
				get_matrice_and_quantiles(input_seq_value_ICDI[header_seq], input_sim_ICDI_rand_aa[header_seq],header_seq_list[j], header_seq, results_for_file_simulation, output_dir)

			j += 1
		if i == 0 :
			firstline = True
		csv_name = "ICDI U_test"
		output_dir_wilcox = "simulation_results/Wilcoxon_Mann_Whitney_results"
		perform_wilcoxon_mann_whitney_u_test_simulation_step(input_sim_ICDI[header_seq], input_sim_ICDI_rand_aa[header_seq], csv_name, header_seq, results_for_file_simulation, output_dir, output_dir_wilcox)

		if i < 10 :
			create_graphs(output_dir, graph_dir, input_seq_value_ICDI[header_seq], legends_names, legends_length[header_seq], type_of_data, matrice[header_seq_list[0]], matrice[header_seq_list[1]])
			i+= 1
		
	i = 0
	firstline = True
	legends_names = []
	legends_length = {}
	for header_seq in content_seq :

		file_name = str(header_seq)[:210] + "_dens_scaled_chi"
		absc_name = "scaled_chi score "
		ord_name = "density score "
		title_name = "scaled_chi SCORE"
		header_name = header_seq
		step_title = "Graph on simulated data from the sequence " + str(header_seq)[:20]
		type_of_data = "a.a"

		header_seq_list = []
		j = 0
		while j < 2 :
			if header_seq not in legends_length :
				legends_length[header_seq] = []
			legends_length[header_seq].append(len(content_seq[header_seq]))
			tmp_header_seq_list = header_seq + "_" + str(j)
			header_seq_list.append(tmp_header_seq_list)

			matrice[header_seq_list[j]] = {}
			quantile_99_left_x[header_seq_list[j]] = {}
			quantile_99_right_x[header_seq_list[j]] = {}
			quantile_95_left_x[header_seq_list[j]] = {}
			quantile_95_right_x[header_seq_list[j]] = {}
			quantile_99_left_y[header_seq_list[j]] = {}
			quantile_99_right_y[header_seq_list[j]] = {}
			quantile_95_left_y[header_seq_list[j]] = {}
			quantile_95_right_y[header_seq_list[j]] = {}
			density_result[header_seq_list[j]] = {}

			if i == 0 :
				firstline = True

			if j == 0 :
				csv_name = "scaled_chi (simu_analysis)"
				graph_name = "same a.a content"
				legends_names.append(graph_name)
				get_matrice_and_quantiles(input_seq_value_scaled_chi[header_seq], input_sim_scaled_chi[header_seq],header_seq_list[j], header_seq, results_for_file_simulation, output_dir)
			else :
				csv_name = "scaled_chi (a.a random guided) (simu_analysis)"
				graph_name = "a.a random-guided"
				legends_names.append(graph_name)
				get_matrice_and_quantiles(input_seq_value_scaled_chi[header_seq], input_sim_scaled_chi_rand_aa[header_seq],header_seq_list[j], header_seq, results_for_file_simulation, output_dir)
			j += 1
		if i == 0 :
			firstline = True
		csv_name = "Scaled_Chi U_test"
		output_dir_wilcox = "simulation_results/Wilcoxon_Mann_Whitney_results"
		perform_wilcoxon_mann_whitney_u_test_simulation_step(input_sim_scaled_chi[header_seq], input_sim_scaled_chi_rand_aa[header_seq], csv_name, header_seq, results_for_file_simulation, output_dir, output_dir_wilcox)

		if i < 10 :
			create_graphs(output_dir, graph_dir, input_seq_value_scaled_chi[header_seq], legends_names, legends_length[header_seq], type_of_data, matrice[header_seq_list[0]], matrice[header_seq_list[1]])
			i+= 1
		

def treat_pattern_analysis_step_data(output_dir) :

	global file_name
	global absc_name
	global ord_name
	global title_name
	global csv_name
	global firstline
	global header_name
	global step_title

	graph_dir = "pattern_analysis_graphs"


	for pattern_number in patterns :

		for pattern in patterns[pattern_number] :
			if len(patterns[pattern_number][pattern]) > 0 :
				step_title = "Graph on sequences containing the pattern " + str(pattern)[:20]

				i = 0
				if len(patterns[pattern_number][pattern]["cousin_value"]) > 1 :
					for header_seq in content_seq :
						firstline = True
						legends_names = []
						legends_length = []

						file_name = str(pattern_number) + "_" + str(pattern) + "_dens_COUSIN_18"
						absc_name = "$\mathrm{COUSIN}_{18}\/\/\mathrm{score}$"
						ord_name = "density score "
						title_name = "$\mathrm{COUSIN}_{18}\/\/\mathrm{score}$"
						header_name = header_seq
						legends_names.append(pattern)
						legends_length.append(len(patterns[pattern_number][pattern]["cousin_value"]))
						type_of_data = "sequences"

						matrice[header_seq] = {}
						quantile_99_left_x[header_seq] = {}
						quantile_99_right_x[header_seq] = {}
						quantile_95_left_x[header_seq] = {}
						quantile_95_right_x[header_seq] = {}
						quantile_99_left_y[header_seq] = {}
						quantile_99_right_y[header_seq] = {}
						quantile_95_left_y[header_seq] = {}
						quantile_95_right_y[header_seq] = {}
						density_result[header_seq] = {}

						csv_name = str(pattern) + " (pattern analysis) "
						graph_name = str(pattern)
						legends_names.append(graph_name)
						get_matrice_and_quantiles(input_seq_value_COUSIN[header_seq], patterns[pattern_number][pattern]["cousin_value"],header_seq, header_seq, results_for_file, output_dir)
						if i == 0 :
							create_graphs(output_dir, graph_dir, None, legends_names, legends_length, type_of_data, matrice[header_seq])

						firstline = True
						legends_names = []
						legends_length = []

						file_name = str(pattern_number) + "_" + str(pattern) + "_dens_COUSIN_59"
						absc_name = "$\mathrm{COUSIN}_{59}\/\/\mathrm{score}$"
						ord_name = "density score "
						title_name = "$\mathrm{COUSIN}_{59}\/\/\mathrm{score}$"
						header_name = header_seq
						legends_names.append(pattern)
						legends_length.append(len(patterns[pattern_number][pattern]["cousin_value_pond"]))
						type_of_data = "sequences"

						matrice[header_seq] = {}
						quantile_99_left_x[header_seq] = {}
						quantile_99_right_x[header_seq] = {}
						quantile_95_left_x[header_seq] = {}
						quantile_95_right_x[header_seq] = {}
						quantile_99_left_y[header_seq] = {}
						quantile_99_right_y[header_seq] = {}
						quantile_95_left_y[header_seq] = {}
						quantile_95_right_y[header_seq] = {}
						density_result[header_seq] = {}

						csv_name = str(pattern) + " (pattern analysis) "
						graph_name = str(pattern)
						legends_names.append(graph_name)
						get_matrice_and_quantiles(input_seq_value_COUSIN_pond[header_seq], patterns[pattern_number][pattern]["cousin_value_pond"],header_seq, header_seq, results_for_file, output_dir)

						if i == 0 :
							create_graphs(output_dir, graph_dir, None, legends_names, legends_length, type_of_data, matrice[header_seq])

						firstline = True
						legends_names = []
						legends_length = []

						file_name = str(pattern_number) + "_" + str(pattern) + "_dens_CAI_59"
						absc_name = "$\mathrm{CAI}_{59}\/\/\mathrm{score}$"
						ord_name = "density score "
						title_name = "$\mathrm{CAI}_{59}\/\/\mathrm{score}$"
						header_name = header_seq
						legends_names.append(pattern)
						legends_length.append(len(patterns[pattern_number][pattern]["CAI_codon"]))
						type_of_data = "sequences"

						matrice[header_seq] = {}
						quantile_99_left_x[header_seq] = {}
						quantile_99_right_x[header_seq] = {}
						quantile_95_left_x[header_seq] = {}
						quantile_95_right_x[header_seq] = {}
						quantile_99_left_y[header_seq] = {}
						quantile_99_right_y[header_seq] = {}
						quantile_95_left_y[header_seq] = {}
						quantile_95_right_y[header_seq] = {}
						density_result[header_seq] = {}

						csv_name = str(pattern) + " (pattern analysis) "
						graph_name = str(pattern)
						legends_names.append(graph_name)
						get_matrice_and_quantiles(input_seq_value_CAI_cod[header_seq], patterns[pattern_number][pattern]["CAI_codon"],header_seq, header_seq,  results_for_file, output_dir)

						if i == 0 :
							create_graphs(output_dir, graph_dir, None, legends_names, legends_length, type_of_data, matrice[header_seq])

						firstline = True
						legends_names = []
						legends_length = []

						file_name = str(pattern_number) + "_" + str(pattern) + "_dens_CAI_18"
						absc_name = "$\mathrm{CAI}_{18}\/\/\mathrm{score}$"
						ord_name = "density score "
						title_name = "$\mathrm{CAI}_{18}\/\/\mathrm{score}$"
						header_name = header_seq
						legends_names.append(pattern)
						legends_length.append(len(patterns[pattern_number][pattern]["CAI_cod_famil"]))
						type_of_data = "sequences"

						matrice[header_seq] = {}
						quantile_99_left_x[header_seq] = {}
						quantile_99_right_x[header_seq] = {}
						quantile_95_left_x[header_seq] = {}
						quantile_95_right_x[header_seq] = {}
						quantile_99_left_y[header_seq] = {}
						quantile_99_right_y[header_seq] = {}
						quantile_95_left_y[header_seq] = {}
						quantile_95_right_y[header_seq] = {}
						density_result[header_seq] = {}

						csv_name = str(pattern) + " (pattern analysis) "
						graph_name = str(pattern)
						legends_names.append(graph_name)
						get_matrice_and_quantiles(input_seq_value_CAI_aa[header_seq], patterns[pattern_number][pattern]["CAI_cod_famil"],header_seq, header_seq, results_for_file, output_dir)

						if i == 0 :
							create_graphs(output_dir, graph_dir, None, legends_names, legends_length, type_of_data, matrice[header_seq])

						firstline = True
						legends_names = []
						legends_length = []

						file_name = str(pattern_number) + "_" + str(pattern) + "_dens_ENC"
						absc_name = "ENC value "
						ord_name = "density score  "
						title_name = "Density curves for ENC measure"
						header_name = header_seq
						legends_names.append(pattern)
						legends_length.append(len(patterns[pattern_number][pattern]["ENC_value"]))
						type_of_data = "sequences"

						matrice[header_seq] = {}
						quantile_99_left_x[header_seq] = {}
						quantile_99_right_x[header_seq] = {}
						quantile_95_left_x[header_seq] = {}
						quantile_95_right_x[header_seq] = {}
						quantile_99_left_y[header_seq] = {}
						quantile_99_right_y[header_seq] = {}
						quantile_95_left_y[header_seq] = {}
						quantile_95_right_y[header_seq] = {}
						density_result[header_seq] = {}

						csv_name = pattern
						legends_names.append(csv_name)
						get_matrice_and_quantiles(input_seq_value_ENC[header_seq], patterns[pattern_number][pattern]["ENC_value"],header_seq, header_seq, results_for_file, output_dir)

						if i == 0 :
							create_graphs(output_dir, graph_dir, None, legends_names, legends_length, type_of_data, matrice[header_seq])

						firstline = True
						legends_names = []
						legends_length = []

						file_name = str(pattern_number) + "_" + str(pattern) + "_dens_SCUO"
						absc_name = "SCUO value "
						ord_name = "density score  "
						title_name = "Density curves for SCUO measure"
						header_name = header_seq
						legends_names.append(pattern)
						legends_length.append(len(patterns[pattern_number][pattern]["SCUO_value"]))
						type_of_data = "sequences"

						matrice[header_seq] = {}
						quantile_99_left_x[header_seq] = {}
						quantile_99_right_x[header_seq] = {}
						quantile_95_left_x[header_seq] = {}
						quantile_95_right_x[header_seq] = {}
						quantile_99_left_y[header_seq] = {}
						quantile_99_right_y[header_seq] = {}
						quantile_95_left_y[header_seq] = {}
						quantile_95_right_y[header_seq] = {}
						density_result[header_seq] = {}

						csv_name = pattern
						legends_names.append(csv_name)
						get_matrice_and_quantiles(input_seq_value_SCUO[header_seq], patterns[pattern_number][pattern]["SCUO_value"],header_seq, header_seq, results_for_file, output_dir)

						if i == 0 :
							create_graphs(output_dir, graph_dir, None, legends_names, legends_length, type_of_data, matrice[header_seq])

						firstline = True
						legends_names = []
						legends_length = []

						file_name = str(pattern_number) + "_" + str(pattern) + "_dens_FOP"
						absc_name = "FOP value "
						ord_name = "density score  "
						title_name = "Density curves for FOP measure"
						header_name = header_seq
						legends_names.append(pattern)
						legends_length.append(len(patterns[pattern_number][pattern]["FOP_value"]))
						type_of_data = "sequences"

						matrice[header_seq] = {}
						quantile_99_left_x[header_seq] = {}
						quantile_99_right_x[header_seq] = {}
						quantile_95_left_x[header_seq] = {}
						quantile_95_right_x[header_seq] = {}
						quantile_99_left_y[header_seq] = {}
						quantile_99_right_y[header_seq] = {}
						quantile_95_left_y[header_seq] = {}
						quantile_95_right_y[header_seq] = {}
						density_result[header_seq] = {}

						csv_name = pattern
						legends_names.append(csv_name)
						get_matrice_and_quantiles(input_seq_value_FOP[header_seq], patterns[pattern_number][pattern]["FOP_value"],header_seq, header_seq, results_for_file, output_dir)

						if i == 0 :
							create_graphs(output_dir, graph_dir, None, legends_names, legends_length, type_of_data, matrice[header_seq])

						firstline = True
						legends_names = []
						legends_length = []

						file_name = str(pattern_number) + "_" + str(pattern) + "_dens_CBI"
						absc_name = "CBI value "
						ord_name = "density score  "
						title_name = "Density curves for CBI measure"
						header_name = header_seq
						legends_names.append(pattern)
						legends_length.append(len(patterns[pattern_number][pattern]["CBI_value"]))
						type_of_data = "sequences"

						matrice[header_seq] = {}
						quantile_99_left_x[header_seq] = {}
						quantile_99_right_x[header_seq] = {}
						quantile_95_left_x[header_seq] = {}
						quantile_95_right_x[header_seq] = {}
						quantile_99_left_y[header_seq] = {}
						quantile_99_right_y[header_seq] = {}
						quantile_95_left_y[header_seq] = {}
						quantile_95_right_y[header_seq] = {}
						density_result[header_seq] = {}

						csv_name = pattern
						legends_names.append(csv_name)
						get_matrice_and_quantiles(input_seq_value_CBI[header_seq], patterns[pattern_number][pattern]["CBI_value"],header_seq, header_seq, results_for_file, output_dir)

						if i == 0 :
							create_graphs(output_dir, graph_dir, None, legends_names, legends_length, type_of_data, matrice[header_seq])

						firstline = True
						legends_names = []
						legends_length = []

						file_name = str(pattern_number) + "_" + str(pattern) + "_dens_ICDI"
						absc_name = "ICDI value "
						ord_name = "density score  "
						title_name = "Density curves for ICDI measure"
						header_name = header_seq
						legends_names.append(pattern)
						legends_length.append(len(patterns[pattern_number][pattern]["ICDI_value"]))
						type_of_data = "sequences"

						matrice[header_seq] = {}
						quantile_99_left_x[header_seq] = {}
						quantile_99_right_x[header_seq] = {}
						quantile_95_left_x[header_seq] = {}
						quantile_95_right_x[header_seq] = {}
						quantile_99_left_y[header_seq] = {}
						quantile_99_right_y[header_seq] = {}
						quantile_95_left_y[header_seq] = {}
						quantile_95_right_y[header_seq] = {}
						density_result[header_seq] = {}

						csv_name = pattern
						legends_names.append(csv_name)
						get_matrice_and_quantiles(input_seq_value_ICDI[header_seq], patterns[pattern_number][pattern]["ICDI_value"],header_seq, header_seq, results_for_file, output_dir)

						if i == 0 :
							create_graphs(output_dir, graph_dir, None, legends_names, legends_length, type_of_data, matrice[header_seq])

						firstline = True
						legends_names = []
						legends_length = []

						file_name = str(pattern_number) + "_" + str(pattern) + "_dens_scaled_chi"
						absc_name = "scaled_chi value "
						ord_name = "density score  "
						title_name = "Density curves for scaled_chi measure"
						header_name = header_seq
						legends_names.append(pattern)
						legends_length.append(len(patterns[pattern_number][pattern]["scaled_chi_value"]))
						type_of_data = "sequences"

						matrice[header_seq] = {}
						quantile_99_left_x[header_seq] = {}
						quantile_99_right_x[header_seq] = {}
						quantile_95_left_x[header_seq] = {}
						quantile_95_right_x[header_seq] = {}
						quantile_99_left_y[header_seq] = {}
						quantile_99_right_y[header_seq] = {}
						quantile_95_left_y[header_seq] = {}
						quantile_95_right_y[header_seq] = {}
						density_result[header_seq] = {}

						csv_name = pattern
						legends_names.append(csv_name)
						get_matrice_and_quantiles(input_seq_value_scaled_chi[header_seq], patterns[pattern_number][pattern]["scaled_chi_value"],header_seq, header_seq, results_for_file, output_dir)

						if i == 0 :
							create_graphs(output_dir, graph_dir, None, legends_names, legends_length, type_of_data, matrice[header_seq])
							i += 1
				else :
					if not os.path.isdir("./" + output_dir + "/log"):
						os.makedirs("./" + output_dir + "/log")

					log_file = open("./" + output_dir + "/log/log_cousin.txt", "a")


					log_file.write("## PATTERN ANALYSIS ERROR : the pattern \"" + str(pattern) + "\" associated to the group \"" + str(pattern_number) + "\" contains only one query. Please change your pattern to make it broader." + "\n")



################################################################################
##################### FUNCTION : CREATE CLUSTERING PLOTS #######################
################################################################################

def do_wilcox_test(matrix_1, matrix_2) :

	matrix_calc_1 = ro.r.c(ro.FloatVector(matrix_1))
	matrix_calc_2 = ro.r.c(ro.FloatVector(matrix_2))

	test = ro.r("wilcox.test")(matrix_calc_1, matrix_calc_2, paired=ro.r("TRUE"), mu=1)

################################################################################
##################### FUNCTION : DO THE CALCULATION STEP #######################
################################################################################

calculation_CAI_cod = {}
calculation_CAI_aa = {}
patterns_by_groups = {}
patterns = {}

def pattern_analysis_step(file_pattern) :

	if file_pattern == "nothing" : #DO WITHOUT PATTERN
		print ("## Note : since you didn't give a pattern file, no pattern analysis will be done on this dataset. That said, a density graph of all CDSs scores will be given. ##")
		pattern_number = ""
		patterns[pattern_number] = {}
		individual_pattern = "no_pattern"
		patterns[pattern_number][individual_pattern] = {}

		for header_seq in content_seq :

			if "CAI_codon" not in patterns[pattern_number][individual_pattern] :
				patterns[pattern_number][individual_pattern]["CAI_codon"] = []
			patterns[pattern_number][individual_pattern]["CAI_codon"].append(input_seq_value_CAI_cod[header_seq])

			if "CAI_cod_famil" not in patterns[pattern_number][individual_pattern] :
				patterns[pattern_number][individual_pattern]["CAI_cod_famil"] = []
			patterns[pattern_number][individual_pattern]["CAI_cod_famil"].append(input_seq_value_CAI_aa[header_seq])

			if "cousin_value" not in patterns[pattern_number][individual_pattern] :
				patterns[pattern_number][individual_pattern]["cousin_value"] = []
			patterns[pattern_number][individual_pattern]["cousin_value"].append(input_seq_value_COUSIN[header_seq])

			if "cousin_value_pond" not in patterns[pattern_number][individual_pattern] :
				patterns[pattern_number][individual_pattern]["cousin_value_pond"] = []
			patterns[pattern_number][individual_pattern]["cousin_value_pond"].append(input_seq_value_COUSIN_pond[header_seq])

			if "ENC_value" not in patterns[pattern_number][individual_pattern] :
				patterns[pattern_number][individual_pattern]["ENC_value"] = []
			patterns[pattern_number][individual_pattern]["ENC_value"].append(input_seq_value_ENC[header_seq])

			if "SCUO_value" not in patterns[pattern_number][individual_pattern] :
				patterns[pattern_number][individual_pattern]["SCUO_value"] = []
			patterns[pattern_number][individual_pattern]["SCUO_value"].append(input_seq_value_SCUO[header_seq])

			if "FOP_value" not in patterns[pattern_number][individual_pattern] :
				patterns[pattern_number][individual_pattern]["FOP_value"] = []
			patterns[pattern_number][individual_pattern]["FOP_value"].append(input_seq_value_FOP[header_seq])

			if "CBI_value" not in patterns[pattern_number][individual_pattern] :
				patterns[pattern_number][individual_pattern]["CBI_value"] = []
			patterns[pattern_number][individual_pattern]["CBI_value"].append(input_seq_value_CBI[header_seq])

			if "ICDI_value" not in patterns[pattern_number][individual_pattern] :
				patterns[pattern_number][individual_pattern]["ICDI_value"] = []
			patterns[pattern_number][individual_pattern]["ICDI_value"].append(input_seq_value_ICDI[header_seq])

			if "scaled_chi_value" not in patterns[pattern_number][individual_pattern] :
				patterns[pattern_number][individual_pattern]["scaled_chi_value"] = []
			patterns[pattern_number][individual_pattern]["scaled_chi_value"].append(input_seq_value_scaled_chi[header_seq])

	else :

		with open(file_pattern, 'r') as file_pattern :
			list_pattern = file_pattern.read()
			file_pattern.close()

		split_line = re.compile('\n').split(list_pattern)
		i = 0
		while i < len(split_line) :
			if split_line[i] != '' :
				patterns[i] = {}
				patterns_by_groups[i] = {}
				split_pattern = re.compile("\t|;|,|\s").split(split_line[i])
				j = 0
				while j < len(split_pattern) :
					patterns_by_groups[i][split_pattern[j]] = {}
					patterns[i][split_pattern[j]] = {}
					j += 1
			i += 1
	for header_seq in content_seq :
		for pattern_number in patterns :
			for individual_pattern in patterns[pattern_number] :
				search_pattern = re.match('.*' + individual_pattern + ".*", header_seq)
				if search_pattern != None :

					if "CAI_codon" not in patterns[pattern_number][individual_pattern] :
						patterns[pattern_number][individual_pattern]["CAI_codon"] = []
					patterns[pattern_number][individual_pattern]["CAI_codon"].append(input_seq_value_CAI_cod[header_seq])

					if "CAI_cod_famil" not in patterns[pattern_number][individual_pattern] :
						patterns[pattern_number][individual_pattern]["CAI_cod_famil"] = []
					patterns[pattern_number][individual_pattern]["CAI_cod_famil"].append(input_seq_value_CAI_aa[header_seq])

					if "cousin_value" not in patterns[pattern_number][individual_pattern] :
						patterns[pattern_number][individual_pattern]["cousin_value"] = []
					patterns[pattern_number][individual_pattern]["cousin_value"].append(input_seq_value_COUSIN[header_seq])

					if "cousin_value_pond" not in patterns[pattern_number][individual_pattern] :
						patterns[pattern_number][individual_pattern]["cousin_value_pond"] = []
					patterns[pattern_number][individual_pattern]["cousin_value_pond"].append(input_seq_value_COUSIN_pond[header_seq])

					if "ENC_value" not in patterns[pattern_number][individual_pattern] :
						patterns[pattern_number][individual_pattern]["ENC_value"] = []
					patterns[pattern_number][individual_pattern]["ENC_value"].append(input_seq_value_ENC[header_seq])

					if "SCUO_value" not in patterns[pattern_number][individual_pattern] :
						patterns[pattern_number][individual_pattern]["SCUO_value"] = []
					patterns[pattern_number][individual_pattern]["SCUO_value"].append(input_seq_value_SCUO[header_seq])

					if "FOP_value" not in patterns[pattern_number][individual_pattern] :
						patterns[pattern_number][individual_pattern]["FOP_value"] = []
					patterns[pattern_number][individual_pattern]["FOP_value"].append(input_seq_value_FOP[header_seq])

					if "CBI_value" not in patterns[pattern_number][individual_pattern] :
						patterns[pattern_number][individual_pattern]["CBI_value"] = []
					patterns[pattern_number][individual_pattern]["CBI_value"].append(input_seq_value_CBI[header_seq])

					if "ICDI_value" not in patterns[pattern_number][individual_pattern] :
						patterns[pattern_number][individual_pattern]["ICDI_value"] = []
					patterns[pattern_number][individual_pattern]["ICDI_value"].append(input_seq_value_ICDI[header_seq])

					if "scaled_chi_value" not in patterns[pattern_number][individual_pattern] :
						patterns[pattern_number][individual_pattern]["scaled_chi_value"] = []
					patterns[pattern_number][individual_pattern]["scaled_chi_value"].append(input_seq_value_scaled_chi[header_seq])


				for pattern_group in patterns_by_groups :
					if pattern_group != pattern_number :
						for individual_pattern_by_groups in patterns_by_groups[pattern_group] :
							if individual_pattern not in patterns_by_groups[pattern_group][individual_pattern_by_groups] :
								patterns_by_groups[pattern_group][individual_pattern_by_groups][individual_pattern] = {}

		for pattern_group in patterns_by_groups :
			for individual_pattern_by_groups in patterns_by_groups[pattern_group] :
				search_pattern_by_groups = re.match('.*' + individual_pattern_by_groups + ".*", header_seq)
				if search_pattern_by_groups != None :
					for individual_pattern in patterns_by_groups[pattern_group][individual_pattern_by_groups] :
						search_pattern = re.match('.*' + individual_pattern + ".*", header_seq)
						if search_pattern != None :

							if "CAI_codon" not in patterns_by_groups[pattern_group][individual_pattern_by_groups][individual_pattern] :
								patterns_by_groups[pattern_group][individual_pattern_by_groups][individual_pattern]["CAI_codon"] = []
							patterns_by_groups[pattern_group][individual_pattern_by_groups][individual_pattern]["CAI_codon"].append(input_seq_value_CAI_cod[header_seq])

							if "CAI_cod_famil" not in patterns_by_groups[pattern_group][individual_pattern_by_groups][individual_pattern] :
								patterns_by_groups[pattern_group][individual_pattern_by_groups][individual_pattern]["CAI_cod_famil"] = []
							patterns_by_groups[pattern_group][individual_pattern_by_groups][individual_pattern]["CAI_cod_famil"].append(input_seq_value_CAI_aa[header_seq])

							if "cousin_value" not in patterns_by_groups[pattern_group][individual_pattern_by_groups][individual_pattern] :
								patterns_by_groups[pattern_group][individual_pattern_by_groups][individual_pattern]["cousin_value"] = []
							patterns_by_groups[pattern_group][individual_pattern_by_groups][individual_pattern]["cousin_value"].append(input_seq_value_COUSIN[header_seq])

							if "cousin_value_pond" not in patterns_by_groups[pattern_group][individual_pattern_by_groups][individual_pattern] :
								patterns_by_groups[pattern_group][individual_pattern_by_groups][individual_pattern]["cousin_value_pond"] = []
							patterns_by_groups[pattern_group][individual_pattern_by_groups][individual_pattern]["cousin_value_pond"].append(input_seq_value_COUSIN_pond[header_seq])

							if "ENC_value" not in patterns_by_groups[pattern_group][individual_pattern_by_groups][individual_pattern] :
								patterns_by_groups[pattern_group][individual_pattern_by_groups][individual_pattern]["ENC_value"] = []
							patterns_by_groups[pattern_group][individual_pattern_by_groups][individual_pattern]["ENC_value"].append(input_seq_value_ENC[header_seq])

							if "SCUO_value" not in patterns_by_groups[pattern_group][individual_pattern_by_groups][individual_pattern] :
								patterns_by_groups[pattern_group][individual_pattern_by_groups][individual_pattern]["SCUO_value"] = []
							patterns_by_groups[pattern_group][individual_pattern_by_groups][individual_pattern]["SCUO_value"].append(input_seq_value_SCUO[header_seq])

							if "FOP_value" not in patterns_by_groups[pattern_group][individual_pattern_by_groups][individual_pattern] :
								patterns_by_groups[pattern_group][individual_pattern_by_groups][individual_pattern]["FOP_value"] = []
							patterns_by_groups[pattern_group][individual_pattern_by_groups][individual_pattern]["FOP_value"].append(input_seq_value_FOP[header_seq])

							if "CBI_value" not in patterns_by_groups[pattern_group][individual_pattern_by_groups][individual_pattern] :
								patterns_by_groups[pattern_group][individual_pattern_by_groups][individual_pattern]["CBI_value"] = []
							patterns_by_groups[pattern_group][individual_pattern_by_groups][individual_pattern]["CBI_value"].append(input_seq_value_CBI[header_seq])

							if "ICDI_value" not in patterns_by_groups[pattern_group][individual_pattern_by_groups][individual_pattern] :
								patterns_by_groups[pattern_group][individual_pattern_by_groups][individual_pattern]["ICDI_value"] = []
							patterns_by_groups[pattern_group][individual_pattern_by_groups][individual_pattern]["ICDI_value"].append(input_seq_value_ICDI[header_seq])

							if "scaled_chi_value" not in patterns_by_groups[pattern_group][individual_pattern_by_groups][individual_pattern] :
								patterns_by_groups[pattern_group][individual_pattern_by_groups][individual_pattern]["scaled_chi_value"] = []
							patterns_by_groups[pattern_group][individual_pattern_by_groups][individual_pattern]["scaled_chi_value"].append(input_seq_value_scaled_chi[header_seq])

################################################################################
############### FUNCTION : CREATE A CU_TABLE FROM INPUT DATA ###################
################################################################################

def create_cu_table(genetic_code_ref, output_dir) : # VERIFIER CODON STOP + /3

	sequential_order_cu_table = ["UUU", "UCU", "UAU", "UGU", "UUC", "UCC", "UAC", "UGC", "UUA", "UCA", "UAA", "UGA", "UUG", "UCG", "UAG", "UGG", "CUU", "CCU", "CAU", "CGU", "CUC", "CCC", "CAC", "CGC", "CUA", "CCA", "CAA", "CGA", "CUG", "CCG", "CAG", "CGG", "AUU", "ACU", "AAU", "AGU", "AUC", "ACC", "AAC", "AGC", "AUA", "ACA", "AAA", "AGA", "AUG", "ACG", "AAG", "AGG", "GUU", "GCU", "GAU", "GGU", "GUC", "GCC", "GAC", "GGC", "GUA", "GCA", "GAA", "GGA", "GUG", "GCG", "GAG", "GGG"]
	tot_aa = 0
	allowed_chars = set('AaTtCcGgUu')


	global occ_tot
	global occ_tot_aa_4_sim
	global occ_tot_aa

	dicodon_tot_codon = {}
	dicodon_tot_codon_by_amino_acid = {}

	dict_codons_frequencies = {}
	best_freq_value = {}


	if not os.path.isdir("./" + output_dir + "/log"):
		os.makedirs("./" + output_dir + "/log")

	log_file = open("./" + output_dir + "/log/log_cousin.txt", "a")

	dicodon_dict = {}

	for amino_acid in genetic_code_ref :
		if amino_acid != "*" : 
			for first_codon in genetic_code_ref[amino_acid] : 

				if first_codon not in dicodon_dict : 
					dicodon_dict[first_codon] = {}

				for amino_acid in genetic_code_ref : 
					for second_codon in genetic_code_ref[amino_acid]:
						if second_codon not in dicodon_dict[first_codon] : 
							dicodon_dict[first_codon][second_codon] = None

	for header_seq in content_seq.copy() :
			i = 0
			sequence = DNA_to_RNA(content_seq[header_seq])
			length_seq = len(sequence) / 3

			for c in sequence:
				if c not in allowed_chars :
					log_file.write("## SEQUENCE ERROR : the sequence " + header_seq  + " contains non-nucleotide characters (check your file) (create table) ## \n")
					try :
						del (content_seq[header_seq])
						break
					except :
						log_file.write("## Oops, this data" + header_seq + " seems already deleted / isn't part of your original dataset. This error should not appear. ##\n")
						break

			while i < length_seq :

				codon_seq = sequence[:3]
				sequence  = sequence[3:] # We suppress the first codons, since it won't be analyzed anymore

				for amino_acid in genetic_code_ref :
					if codon_seq in genetic_code_ref[amino_acid] :
						if amino_acid == "*" and len(sequence) > 0 :
							log_file.write("## SEQUENCE ERROR : the sequence " + header_seq  + " has a STOP codon in another place than at the end of it (check your file) (create table) ## \n")
							try :
								del (content_seq[header_seq])
								i = length_seq
								break
							except :
								log_file.write("## Oops, this data" + header_seq + " seems already deleted / isn't part of your original dataset. ## \n")
								i = length_seq
								break
				i += 1

					
	for header_seq in content_seq :
			i = 0
			sequence = DNA_to_RNA(content_seq[header_seq])
			length_seq = len(sequence) / 3
			while i < length_seq :

				dicodon_seq = sequence[:6] # We take the two codons
				first_codon = dicodon_seq[:int(len(dicodon_seq)/2)]
				second_codon = dicodon_seq[int(len(dicodon_seq)/2):]
				codon_seq = sequence[:3]
				sequence  = sequence[3:]

				if i < length_seq - 1 :

					if dicodon_dict[first_codon][second_codon] == None : 
						dicodon_dict[first_codon][second_codon] = 1
					else :
						dicodon_dict[first_codon][second_codon] += 1

				for amino_acid in genetic_code_ref :
					
					if amino_acid not in occ_tot_aa :
						occ_tot_aa[amino_acid] = 0
					if codon_seq in genetic_code_ref[amino_acid] :
						occ_tot_aa[amino_acid] += 1
						if not genetic_code_ref[amino_acid][codon_seq] :
							genetic_code_ref[amino_acid][codon_seq].append(1)
						else :
							genetic_code_ref[amino_acid][codon_seq][0] += 1
						tot_aa += 1

				i += 1

	for amino_acid in genetic_code_ref :
		if amino_acid not in occ_tot_ref :
			occ_tot_ref[amino_acid] = occ_tot_aa[amino_acid] # pour la fonction unvoid_empty_frequencies
		occ_tot += occ_tot_aa[amino_acid]
		if amino_acid != '*' :
			occ_tot_aa_4_sim += occ_tot_aa[amino_acid]
		for codon in genetic_code_ref[amino_acid] :
			if not genetic_code_ref[amino_acid][codon] :
				genetic_code_ref[amino_acid][codon].append(0)
			if occ_tot_aa[amino_acid] != 0 :
				genetic_code_ref[amino_acid][codon].append(round((genetic_code_ref[amino_acid][codon][0] / occ_tot_aa[amino_acid]), 3))
				if amino_acid not in best_freq_value : 
					best_freq_value[amino_acid] = genetic_code_ref[amino_acid][codon][1]
				else : 
					if genetic_code_ref[amino_acid][codon][1] > best_freq_value[amino_acid] : 
						best_freq_value[amino_acid] = genetic_code_ref[amino_acid][codon][1]

			else :
				genetic_code_ref[amino_acid][codon].append(0)
				dict_codons_frequencies[codon] = 0

				if amino_acid not in best_freq_value : 
					best_freq_value[amino_acid] = genetic_code_ref[amino_acid][codon][1]
				else : 
					if genetic_code_ref[amino_acid][codon][1] > best_freq_value[amino_acid] : 
						best_freq_value[amino_acid] = genetic_code_ref[amino_acid][codon][1]



	for amino_acid in genetic_code_ref :
		for codon in genetic_code_ref[amino_acid] :
			if best_freq_value[amino_acid] != 0 : 
				dict_codons_frequencies[codon] = round((genetic_code_ref[amino_acid][codon][1] / best_freq_value[amino_acid]), 3)
			else :
				dict_codons_frequencies[codon] = 0




	cpt = 0


	for first_codon in dicodon_dict :
		for second_codon in dicodon_dict[first_codon] : 

			if dicodon_dict[first_codon][second_codon] == None : 
				dicodon_dict[first_codon][second_codon] = []
				dicodon_dict[first_codon][second_codon].append(0)
				dicodon_dict[first_codon][second_codon].append(0) # get the frequency of this dicodon
				if first_codon not in dicodon_tot_codon : 
					dicodon_tot_codon[first_codon] = 1 # = 1 only to avoid division by 0. Need to check if it'll be important later.

			else :				
				value = dicodon_dict[first_codon][second_codon]
				dicodon_dict[first_codon][second_codon] = []
				dicodon_dict[first_codon][second_codon].append(value)
				dicodon_dict[first_codon][second_codon].append(value / (occ_tot - len(content_seq))) # get the frequency of this dicodon. len(content_seq) allow to suppress the last aa of each sequence
				if first_codon not in dicodon_tot_codon : 
					dicodon_tot_codon[first_codon] = value
				else : 
					dicodon_tot_codon[first_codon] += value

				# for amino_acid in genetic_code_ref : 
				# 	if second_codon in genetic_code_ref[amino_acid] : 
				# 		if first_codon not in dicodon_tot_codon_by_amino_acid : 
				# 			dicodon_tot_codon_by_amino_acid[first_codon] = {}
				# 			dicodon_tot_codon_by_amino_acid[first_codon][amino_acid] = value
				# 		else : 
				# 			if amino_acid not in dicodon_tot_codon_by_amino_acid[first_codon] : 
				# 				dicodon_tot_codon_by_amino_acid[first_codon][amino_acid] = value
				# 			else : 
				# 				dicodon_tot_codon_by_amino_acid[first_codon][amino_acid] += value

	if not os.path.isdir("./" + output_dir + "/codon_usage_tables"):
		os.makedirs("./" + output_dir + "/codon_usage_tables")

	cu_table = open("./" + output_dir + "/codon_usage_tables/cu_table.txt", "a")

	for codon in sequential_order_cu_table :
		for amino_acid in genetic_code_ref :
			if codon in genetic_code_ref[amino_acid] :
				number_of_spaces_occ = " " * (6 - len(str(genetic_code_ref[amino_acid][codon][0])))
				number_of_spaces_freq = " " * (5 - len(str(round(((genetic_code_ref[amino_acid][codon][0]/tot_aa) * 1000), 1))))
				if cpt < 4 :
					cu_table.write(codon + number_of_spaces_freq + str(round(((genetic_code_ref[amino_acid][codon][0]/tot_aa) * 1000), 1)) + "(" + number_of_spaces_occ + str(genetic_code_ref[amino_acid][codon][0]) + ")  ")
					cpt += 1
				elif cpt == 4 :
					cu_table.write("\n")
					cu_table.write(codon + number_of_spaces_freq + str(round(((genetic_code_ref[amino_acid][codon][0]/tot_aa) * 1000), 1)) + "(" + number_of_spaces_occ + str(genetic_code_ref[amino_acid][codon][0]) + ")  ")
					cpt = 1

	di_cu_table = open("./" + output_dir + "/codon_usage_tables/di_cu_table.txt", "a")

	for first_codon in sequential_order_cu_table :
		for second_codon in sequential_order_cu_table :  
			if first_codon in dicodon_dict : 
				if second_codon in dicodon_dict[first_codon] :
					number_of_spaces_occ = " " * (6 - len(str(dicodon_dict[first_codon][second_codon][0])))
					number_of_spaces_freq = " " * (5 - len(str(round(((dicodon_dict[first_codon][second_codon][0] / dicodon_tot_codon[first_codon]) * 1000), 1))))
					if cpt < 4 :
						di_cu_table.write(first_codon + "_" + second_codon + number_of_spaces_freq + str(round(((dicodon_dict[first_codon][second_codon][0] / dicodon_tot_codon[first_codon]) * 1000), 1)) + "(" + number_of_spaces_occ + str(dicodon_dict[first_codon][second_codon][0]) + ")  ")
						cpt += 1
					elif cpt == 4 :
						di_cu_table.write("\n")
						di_cu_table.write(first_codon + "_" + second_codon + number_of_spaces_freq + str(round(((dicodon_dict[first_codon][second_codon][0] / dicodon_tot_codon[first_codon]) * 1000), 1)) + "(" + number_of_spaces_occ + str(dicodon_dict[first_codon][second_codon][0]) + ")  ")
						cpt = 1
		di_cu_table.write("\n")

	cu_table.close()
	di_cu_table.close()
	log_file.close()

	create_cut_graph(genetic_code_ref, dict_codons_frequencies, "frequencies_of_synonymous_codons", output_dir)

	return genetic_code_ref

################################################################################
################################################################################
############################# MISCELLANOUS #####################################
################################################################################
################################################################################

################################################################################
#################### Creating output directory #################################
################################################################################

def create_output_dir(directory, step_used, routine_simulation, perform_vector_analysis) :

	path = "./"
	if not os.path.isdir(directory):
		os.makedirs(directory)

	else :
		sys.exit("Error, your directory already exists. Please change the name of your directory.")

	if not os.path.isdir(path + directory + "/codon_usage_tables"):
		os.makedirs(path + directory + "/codon_usage_tables")

	if routine_simulation == "yes" or routine_simulation == "Yes" or routine_simulation == "YES" or routine_simulation == "Y" or routine_simulation == "y" :
		if not os.path.isdir(path + directory + "/basic_simulation_graphs"):
			os.makedirs(path + directory + "/basic_simulation_graphs")

	if not os.path.isdir(path + directory + "/log"):
		os.makedirs(path + directory + "/log")
	log_file = open(path + directory + "/log/log_cousin.txt", "w")
	log_file.write("")
	log_file.close()

	if perform_vector_analysis == "Yes" or perform_vector_analysis == "yes" or perform_vector_analysis == "YES" or perform_vector_analysis == "Y" or perform_vector_analysis == "y" : 
		if not os.path.isdir(path + directory + "/vector_data"):
			os.makedirs(path + directory + "/vector_data")

	if step_used == "create_table" :
		if not os.path.isdir(path + directory + "/codon_usage_tables"):
			os.makedirs(path + directory + "/codon_usage_tables")
		cu_table = open(path + directory + "/codon_usage_tables/cu_table.txt", "w")
		cu_table.write("")
		cu_table.close()

	if step_used == "optimization" :
		if not os.path.isdir(path + directory + "/optimization"):
			os.makedirs(path + directory + "/optimization")
		optimized_sequences = open(path + directory + "/optimization/optimized_sequences.fasta", "w") #NEED TO CHANGE VAR NAME TO A BETTER ONE AND CHECK wITH OPTIMIZE FUNCTION\
		optimized_sequences.write("")
		optimized_sequences.close()

	if step_used == "clustering_analysis" :
		if not os.path.isdir(path + directory + "/clustering"):
			os.makedirs(path + directory + "/clustering")

	if step_used == "pattern_analysis" :
		if not os.path.isdir(path + directory + "/pattern_analysis_graphs"):
			os.makedirs(path + directory + "/pattern_analysis_graphs")

	if step_used == "simulation_analysis" :
		if not os.path.isdir(path + directory + "/simulation_graphs"):
			os.makedirs(path + directory + "/simulation_graphs")
		if not os.path.isdir(path + directory + "/simulation_results"):
			os.makedirs(path + directory + "/simulation_results")
		if not os.path.isdir(path + directory + "/simulation_results/Wilcoxon_Mann_Whitney_results"):
			os.makedirs(path + directory + "/simulation_results/Wilcoxon_Mann_Whitney_results")

################################################################################
#################### FUNCTION USED TO GENERATE THE SEQUENCES ###################
################################################################################

### CAUTION : THIS FUNCTION LAUNCHES THE CAI AND COUSIN CALCULATION (DESCRIBED BELOW)

def seq_simulations(genetic_code_ref, number_of_repetitions, size_1, size_2, size_3, optimal_codons_list, perform_COUSIN_with_boundaries, output_dir) :

	i = 0

	CAI_codon_calc[size_1] = []
	CAI_codon_calc[size_2] = []
	CAI_codon_calc[size_3] = []

	CAI_aa_calc[size_1] = []
	CAI_aa_calc[size_2] = []
	CAI_aa_calc[size_3] = []

	COUSIN_calc[size_1] = []
	COUSIN_calc[size_2] = []
	COUSIN_calc[size_3] = []

	COUSIN_calc_pond[size_1] = []
	COUSIN_calc_pond[size_2] = []
	COUSIN_calc_pond[size_3] = []

	ENC_calc[size_1] = []
	ENC_calc[size_2] = []
	ENC_calc[size_3] = []

	SCUO_calc[size_1] = []
	SCUO_calc[size_2] = []
	SCUO_calc[size_3] = []

	FOP_calc[size_1] = []
	FOP_calc[size_2] = []
	FOP_calc[size_3] = []

	CBI_calc[size_1] = []
	CBI_calc[size_2] = []
	CBI_calc[size_3] = []

	ICDI_calc[size_1] = []
	ICDI_calc[size_2] = []
	ICDI_calc[size_3] = []

	scaled_chi_calc[size_1] = []
	scaled_chi_calc[size_2] = []
	scaled_chi_calc[size_3] = []

	occ_aa_sim = {}

	calc_freq_equ_muta_bias(genetic_code_ref)
	genetic_code_ref_nucl_comp = generate_table_nucl_comp(genetic_code_ref, output_dir)
	deviation_score = calc_deviation_scores(genetic_code_ref)
	deviation_score_nucl_comp = calc_deviation_scores(genetic_code_ref_nucl_comp)
	w_scores = calc_adaptative_values(genetic_code_ref)

	while i < number_of_repetitions :
		occ_aa_sim[i] = {}
		genetic_code_sim[i] = dict.fromkeys(genetic_code_ref)
		for amino_acid in genetic_code_sim[i] :
			genetic_code_sim[i][amino_acid] = dict.fromkeys(genetic_code_ref[amino_acid])
			for codon in genetic_code_sim[i][amino_acid] :
				genetic_code_sim[i][amino_acid][codon] = []
				genetic_code_sim[i][amino_acid][codon].append(0)
		genetic_code_sim_nucl_comp[i] = dict.fromkeys(genetic_code_ref_nucl_comp)
		for amino_acid in genetic_code_sim_nucl_comp[i] :
			genetic_code_sim_nucl_comp[i][amino_acid] = dict.fromkeys(genetic_code_ref_nucl_comp[amino_acid])
			for codon in genetic_code_sim_nucl_comp[i][amino_acid] :
				genetic_code_sim_nucl_comp[i][amino_acid][codon] = []
				genetic_code_sim_nucl_comp[i][amino_acid][codon].append(0)

		genetic_code_sim[i] = create_sequences_from_CU_table(size_1, genetic_code_sim[i], genetic_code_ref)
		aa_in_req = check_if_aa_in_req(genetic_code_sim[i])
		deviation_score_aa = calc_deviation_scores_aa(genetic_code_ref, aa_in_req)
		RSCU_data = calc_RSCU(genetic_code_sim[i])

		CAI_codon_calc[size_1].append(calc_CAI_59(genetic_code_sim[i], genetic_code_ref, w_scores))
		CAI_aa_calc[size_1].append(calc_CAI_18(genetic_code_sim[i], genetic_code_ref, w_scores))
		COUSIN_calc[size_1].append(calc_COUSIN_18(genetic_code_sim[i], genetic_code_ref, deviation_score, aa_in_req, perform_COUSIN_with_boundaries))
		COUSIN_calc_pond[size_1].append(calc_COUSIN_59(genetic_code_sim[i], genetic_code_ref, deviation_score, deviation_score_aa, aa_in_req, perform_COUSIN_with_boundaries))
	
		ENC_calc[size_1].append(calc_ENC(genetic_code_sim[i]))
		SCUO_calc[size_1].append(calc_SCUO(genetic_code_sim[i], genetic_code_ref))
		FOP_calc[size_1].append(calc_FOP(genetic_code_sim[i], optimal_codons_list))
		CBI_calc[size_1].append(calc_CBI(genetic_code_sim[i], optimal_codons_list))
		ICDI_calc[size_1].append(calc_ICDI(genetic_code_sim[i], RSCU_data, aa_in_req))
		scaled_chi_calc[size_1].append(calc_scaled_chi(genetic_code_sim[i]))

		genetic_code_sim[i] = create_sequences_from_CU_table(size_2, genetic_code_sim[i], genetic_code_ref)
		aa_in_req = check_if_aa_in_req(genetic_code_sim[i])
		deviation_score_aa = calc_deviation_scores_aa(genetic_code_ref, aa_in_req)

		RSCU_data = calc_RSCU(genetic_code_sim[i])

		CAI_codon_calc[size_2].append(calc_CAI_59(genetic_code_sim[i], genetic_code_ref, w_scores))
		CAI_aa_calc[size_2].append(calc_CAI_18(genetic_code_sim[i], genetic_code_ref, w_scores))
		COUSIN_calc[size_2].append(calc_COUSIN_18(genetic_code_sim[i], genetic_code_ref, deviation_score, aa_in_req, perform_COUSIN_with_boundaries))
		COUSIN_calc_pond[size_2].append(calc_COUSIN_59(genetic_code_sim[i], genetic_code_ref, deviation_score, deviation_score_aa, aa_in_req, perform_COUSIN_with_boundaries))
		
		ENC_calc[size_2].append(calc_ENC(genetic_code_sim[i]))
		SCUO_calc[size_2].append(calc_SCUO(genetic_code_sim[i], genetic_code_ref))
		FOP_calc[size_2].append(calc_FOP(genetic_code_sim[i], optimal_codons_list))
		CBI_calc[size_2].append(calc_CBI(genetic_code_sim[i], optimal_codons_list))
		ICDI_calc[size_2].append(calc_ICDI(genetic_code_sim[i], RSCU_data, aa_in_req))
		scaled_chi_calc[size_2].append(calc_scaled_chi(genetic_code_sim[i]))

		genetic_code_sim[i] = create_sequences_from_CU_table(size_3, genetic_code_sim[i], genetic_code_ref)
		aa_in_req = check_if_aa_in_req(genetic_code_sim[i])
		deviation_score_aa = calc_deviation_scores_aa(genetic_code_ref, aa_in_req)

		RSCU_data = calc_RSCU(genetic_code_sim[i])

		CAI_codon_calc[size_3].append(calc_CAI_59(genetic_code_sim[i], genetic_code_ref, w_scores))
		CAI_aa_calc[size_3].append(calc_CAI_18(genetic_code_sim[i], genetic_code_ref, w_scores))
		COUSIN_calc[size_3].append(calc_COUSIN_18(genetic_code_sim[i], genetic_code_ref, deviation_score, aa_in_req, perform_COUSIN_with_boundaries))
		COUSIN_calc_pond[size_3].append(calc_COUSIN_59(genetic_code_sim[i], genetic_code_ref, deviation_score, deviation_score_aa, aa_in_req, perform_COUSIN_with_boundaries))
		
		ENC_calc[size_3].append(calc_ENC(genetic_code_sim[i]))
		SCUO_calc[size_3].append(calc_SCUO(genetic_code_sim[i], genetic_code_ref))
		FOP_calc[size_3].append(calc_FOP(genetic_code_sim[i], optimal_codons_list))
		CBI_calc[size_3].append(calc_CBI(genetic_code_sim[i], optimal_codons_list))
		ICDI_calc[size_3].append(calc_ICDI(genetic_code_sim[i], RSCU_data, aa_in_req))
		scaled_chi_calc[size_3].append(calc_scaled_chi(genetic_code_sim[i]))


		i += 1

def input_seq_simulation_analysis(genetic_code_ref, number_of_repetitions, optimal_codons_list, perform_COUSIN_with_boundaries, output_dir) :

	calc_freq_equ_muta_bias(genetic_code_ref)
	genetic_code_ref_nucl_comp = generate_table_nucl_comp(genetic_code_ref, output_dir)
	deviation_score = calc_deviation_scores(genetic_code_ref)
	deviation_score_nucl_comp = calc_deviation_scores(genetic_code_ref_nucl_comp)
	w_scores = calc_adaptative_values(genetic_code_ref)

	for header_seq in content_seq :

		input_sim_CAI_codon_calc[header_seq] = []
		input_sim_CAI_aa_calc[header_seq] = []
		input_sim_COUSIN_calc[header_seq] = []
		input_sim_COUSIN_pond_calc[header_seq] = []
		# input_sim_COUSIN_calc_muta_bias[header_seq] = []
		# input_sim_COUSIN_calc_muta_bias_pond[header_seq] = []

		input_sim_ENC[header_seq] = []
		input_sim_SCUO[header_seq] = []
		input_sim_FOP[header_seq] = []
		input_sim_CBI[header_seq] = []
		input_sim_ICDI[header_seq] = []
		input_sim_scaled_chi[header_seq] = []

		input_sim_CAI_codon_calc_rand_aa[header_seq] = []
		input_sim_CAI_aa_calc_rand_aa[header_seq] = []
		input_sim_COUSIN_calc_rand_aa[header_seq] = []
		input_sim_COUSIN_pond_calc_rand_aa[header_seq] = []
		# input_sim_COUSIN_calc_muta_bias_rand_aa[header_seq] = []
		# input_sim_COUSIN_calc_muta_bias_pond_rand_aa[header_seq] = []

		input_sim_ENC_rand_aa[header_seq] = []
		input_sim_SCUO_rand_aa[header_seq] = []
		input_sim_FOP_rand_aa[header_seq] = []
		input_sim_CBI_rand_aa[header_seq] = []
		input_sim_ICDI_rand_aa[header_seq] = []
		input_sim_scaled_chi_rand_aa[header_seq] = []


		occ_aa_sim = {}

		calc_freq_equ_muta_bias(genetic_code_ref)
		i = 0
		seq_dict = {}

		rna_seq = DNA_to_RNA(content_seq[header_seq])
		aa_seq = RNA_to_AA(rna_seq, genetic_code_ref)

		while (i < number_of_repetitions) :

			seq_simulated = random_guided_simulation(aa_seq, genetic_code_ref)
			seq_dict = seq_to_dict(seq_simulated, "simulation_data", genetic_code_ref, output_dir)

			if seq_dict != 'deleted' :

				aa_in_req = check_if_aa_in_req(seq_dict)
				deviation_score_aa = calc_deviation_scores_aa(genetic_code_ref, aa_in_req)

				RSCU_data = calc_RSCU(seq_dict)

				input_sim_CAI_codon_calc[header_seq].append(calc_CAI_59(seq_dict, genetic_code_ref, w_scores))
				input_sim_CAI_aa_calc[header_seq].append(calc_CAI_18(seq_dict, genetic_code_ref, w_scores))
				input_sim_COUSIN_calc[header_seq].append(calc_COUSIN_18(seq_dict, genetic_code_ref, deviation_score, aa_in_req, perform_COUSIN_with_boundaries))
				input_sim_COUSIN_pond_calc[header_seq].append(calc_COUSIN_59(seq_dict, genetic_code_ref, deviation_score, deviation_score_aa, aa_in_req, perform_COUSIN_with_boundaries))
				
				input_sim_ENC[header_seq].append(calc_ENC(seq_dict))
				input_sim_SCUO[header_seq].append(calc_SCUO(seq_dict, genetic_code_ref))
				input_sim_FOP[header_seq].append(calc_FOP(seq_dict, optimal_codons_list))
				input_sim_CBI[header_seq].append(calc_CBI(seq_dict, optimal_codons_list))
				input_sim_ICDI[header_seq].append(calc_ICDI(seq_dict, RSCU_data, aa_in_req))
				input_sim_scaled_chi[header_seq].append(calc_scaled_chi(seq_dict))

				genetic_code_sim_rand_aa = dict.fromkeys(genetic_code_ref)

				for amino_acid in genetic_code_sim_rand_aa :
					genetic_code_sim_rand_aa[amino_acid] = dict.fromkeys(genetic_code_ref[amino_acid])
					for codon in genetic_code_sim_rand_aa[amino_acid] :
						genetic_code_sim_rand_aa[amino_acid][codon] = []
						genetic_code_sim_rand_aa[amino_acid][codon].append(0)

				create_sequences_from_CU_table(len(aa_seq), genetic_code_sim_rand_aa, genetic_code_ref)
				aa_in_req = check_if_aa_in_req(genetic_code_sim_rand_aa)
				deviation_score_aa = calc_deviation_scores_aa(genetic_code_ref, aa_in_req)



				input_sim_CAI_codon_calc_rand_aa[header_seq].append(calc_CAI_59(genetic_code_sim_rand_aa, genetic_code_ref, w_scores))
				input_sim_CAI_aa_calc_rand_aa[header_seq].append(calc_CAI_18(genetic_code_sim_rand_aa, genetic_code_ref, w_scores))
				input_sim_COUSIN_calc_rand_aa[header_seq].append(calc_COUSIN_18(genetic_code_sim_rand_aa, genetic_code_ref, deviation_score, aa_in_req, perform_COUSIN_with_boundaries))
				input_sim_COUSIN_pond_calc_rand_aa[header_seq].append(calc_COUSIN_59(genetic_code_sim_rand_aa, genetic_code_ref, deviation_score, deviation_score_aa, aa_in_req, perform_COUSIN_with_boundaries))
				# input_sim_COUSIN_calc_muta_bias_rand_aa[header_seq].append(calc_COUSIN_18(genetic_code_sim_rand_aa, genetic_code_ref_nucl_comp, deviation_score_nucl_comp, aa_in_req))
				# input_sim_COUSIN_calc_muta_bias_pond_rand_aa[header_seq].append(calc_COUSIN_59(genetic_code_sim_rand_aa, genetic_code_ref_nucl_comp, deviation_score_nucl_comp, aa_in_req))
				input_sim_ENC_rand_aa[header_seq].append(calc_ENC(seq_dict))
				input_sim_SCUO_rand_aa[header_seq].append(calc_SCUO(seq_dict, genetic_code_ref))
				input_sim_FOP_rand_aa[header_seq].append(calc_FOP(seq_dict, optimal_codons_list))
				input_sim_CBI_rand_aa[header_seq].append(calc_CBI(seq_dict, optimal_codons_list))
				input_sim_ICDI_rand_aa[header_seq].append(calc_ICDI(seq_dict, RSCU_data, aa_in_req))
				input_sim_scaled_chi_rand_aa[header_seq].append(calc_scaled_chi(seq_dict))

				i += 1

def calc_eucl_dist_que_vs_ref(genetic_code_req, genetic_code_ref, is_aa_in_req) :
	eucl_dist = 0
	for amino_acid in genetic_code_ref :
		if is_aa_in_req[amino_acid] == True : 
			eucl_dist += math.pow((genetic_code_ref[amino_acid] - genetic_code_req[amino_acid]), 2)

	eucl_dist = math.sqrt(eucl_dist)
	return eucl_dist

def calc_eucl_dist_que_vs_H0(genetic_code_req, is_aa_in_req) :
	eucl_dist = 0
	cpt = 0

	for amino_acid in genetic_code_req :
		if is_aa_in_req[amino_acid] == True : 
			cpt += 1 
	freq_value_H0 = (1/cpt)

	for amino_acid in genetic_code_req : 
		eucl_dist += math.pow((genetic_code_req[amino_acid]) - freq_value_H0, 2)

	eucl_dist = math.sqrt(eucl_dist)
	return eucl_dist


def perform_clustering_analysis(results, genetic_code_ref, number_of_clusters, type_of_analysis, output_dir) :

	results_line = {}
	results_element = {}
	variable_PCA = None
	size_line = 0
	results_element_for_PCA = []

	sns.set(rc={'figure.figsize':(14,10)})
	sns.set(color_codes=True)
	sns.set_style("whitegrid")

	set_of_colors_colorblind = sns.color_palette("colorblind")
	set_of_colors_colorblind.extend(["#b15928","#9b59b6","#fb9a99","#ffff99","#6a3d9a","#1f78b4","#95a5a6","#34495e", "#2ecc71"])
	sns.set_palette(sns.color_palette(set_of_colors_colorblind, n_colors=15))

	if type_of_analysis == "variables" : 

		if len(results) > 1 :
			for line in results :
				if size_line == 0 :
					size_line = len(results[line].split("\t"))
				results_element[line] = results[line].split("\t")
				if line == "firstline" :
					headers_for_PCA = results_element[line]
				else :
					results_element_for_PCA.append(results_element[line])


			results_for_PCA = pd.DataFrame(results_element_for_PCA, columns = headers_for_PCA)
			results_for_PCA = results_for_PCA.set_index("Header")
			results_for_PCA =  results_for_PCA.drop(results_for_PCA.columns[22:], axis=1) #all

			headers_for_PCA = headers_for_PCA[1:23]


	elif type_of_analysis == "GC_content" : 

		if len(results) > 1 :
			for line in results :
				if size_line == 0 :
					size_line = len(results[line].split("\t"))
				results_element[line] = results[line].split("\t")
				if line == "firstline" :
					headers_for_PCA = results_element[line]
				else :
					results_element_for_PCA.append(results_element[line])


			results_for_PCA = pd.DataFrame(results_element_for_PCA, columns = headers_for_PCA)
			results_for_PCA = results_for_PCA.set_index("Header")
			results_for_PCA =  results_for_PCA.drop(results_for_PCA.columns[0:5], axis=1) #all
			results_for_PCA =  results_for_PCA.drop(results_for_PCA.columns[4:], axis=1) 

			headers_for_PCA = list(results_for_PCA)

	elif type_of_analysis == "ATGC_content" : 

		if len(results) > 1 :
			for line in results :
				if size_line == 0 :
					size_line = len(results[line].split("\t"))
				results_element[line] = results[line].split("\t")
				if line == "firstline" :
					headers_for_PCA = results_element[line]
				else :
					results_element_for_PCA.append(results_element[line])


			results_for_PCA = pd.DataFrame(results_element_for_PCA, columns = headers_for_PCA)
			results_for_PCA = results_for_PCA.set_index("Header")
			results_for_PCA =  results_for_PCA.drop(results_for_PCA.columns[0:1], axis=1) #all
			results_for_PCA =  results_for_PCA.drop(results_for_PCA.columns[4:], axis=1) 

			headers_for_PCA = list(results_for_PCA)



	elif type_of_analysis == "all_indexes" : 

		if len(results) > 1 :
			for line in results :
				if size_line == 0 :
					size_line = len(results[line].split("\t"))
				results_element[line] = results[line].split("\t")
				if line == "firstline" :
					headers_for_PCA = results_element[line]
				else :
					results_element_for_PCA.append(results_element[line])


			results_for_PCA = pd.DataFrame(results_element_for_PCA, columns = headers_for_PCA)
			results_for_PCA = results_for_PCA.set_index("Header")
			results_for_PCA =  results_for_PCA.drop(results_for_PCA.columns[0:9], axis=1) #all

			headers_for_PCA = list(results_for_PCA)


	elif type_of_analysis == "basic_analysis" : 

		if len(results) > 1 :
			for line in results :
				if size_line == 0 :
					size_line = len(results[line].split("\t"))
				results_element[line] = results[line].split("\t")
				if line == "firstline" :
					headers_for_PCA = results_element[line]
				else :
					results_element_for_PCA.append(results_element[line])


			results_for_PCA = pd.DataFrame(results_element_for_PCA, columns = headers_for_PCA)
			results_for_PCA = results_for_PCA.set_index("Header")
			results_for_PCA =  results_for_PCA.drop(results_for_PCA.columns[0:5], axis=1) #all
			results_for_PCA =  results_for_PCA.drop(results_for_PCA.columns[8:], axis=1) 

			headers_for_PCA = list(results_for_PCA)


	elif type_of_analysis == "codons" : 

		results = {}
		results["firstline"] = "Header"
		list_codons = []

		for amino_acid in genetic_code_ref : 
			if len(genetic_code_ref[amino_acid]) > 1 and amino_acid != '*' : 
				for codon in genetic_code_ref[amino_acid] : 
					list_codons.append(codon)
					results["firstline"] += "\t" + codon


		
		for header_seq in content_seq.copy() :

			occ_aa = {}
			freq_aa_req = {}
			ATGC_percent = get_freq_ATGC_seq(content_seq[header_seq])
			GC_percent = get_freq_GC_seq(content_seq[header_seq])
			seq_DNA_to_RNA = DNA_to_RNA(content_seq[header_seq])
			seq_dict = seq_to_dict(seq_DNA_to_RNA, header_seq, genetic_code_ref, output_dir)

			if seq_dict != 'deleted' :
				for amino_acid in seq_dict : 
					if len(seq_dict[amino_acid]) > 1 and amino_acid != '*' : 
						for codon in seq_dict[amino_acid] : 
							for element in list_codons : 
								if codon == element :  
									if header_seq not in results : 
										results[header_seq] = str(header_seq)
										results[header_seq] += "\t" + str(seq_dict[amino_acid][codon][1])
									else : 
										results[header_seq] += "\t" + str(seq_dict[amino_acid][codon][1])


		if len(results) > 1 :
			for line in results :
				if size_line == 0 :
					size_line = len(results[line].split("\t"))
				results_element[line] = results[line].split("\t")
				if line == "firstline" :

					headers_for_PCA = results_element[line]
				else :
					results_element_for_PCA.append(results_element[line])

			results_for_PCA = pd.DataFrame(results_element_for_PCA, columns = headers_for_PCA)
			results_for_PCA = results_for_PCA.set_index("Header")


		headers_for_PCA = headers_for_PCA[1:]

	if len(results) > 1 :

		n = min(len(results_for_PCA.index),len(results_for_PCA.columns))
		print (n)
		pca = PCA(n_components = n)
		X = StandardScaler().fit_transform(results_for_PCA)
		pca.fit(X)

		pca_vector_1 = pca.components_[0]
		pca_vector_2 = pca.components_[1]
		pca_vector_3 = pca.components_[2]
		pca_vector_4 = pca.components_[3]

		PCA_data = pca.transform(X)

		list_of_PC = []
		list_headers_PC = []
		i = 1
		for PC in pca.components_ :

			list_headers_PC.append("PC " + str(i))
			list_of_PC.append(PC.tolist())
			i += 1

		PCA_results_by_var = pd.DataFrame(pca.components_, index = list_headers_PC, columns = headers_for_PCA)
		PCA_results_by_var.to_csv("./" + output_dir + "/clustering/PCA_results_by_var.txt", sep='\t', encoding="utf-8")

		sns_plot = sns.heatmap(PCA_results_by_var)

		fig = sns_plot.get_figure()

		var_plot = fig.add_subplot(111)

		plt.suptitle('Impact of each variable on PCs', fontsize = 22)
		var_plot.set_xlabel('PC component', fontsize = 18)
		var_plot.set_ylabel('variable', fontsize = 18)
		plt.xticks(fontsize=14, rotation=90)
		plt.yticks(fontsize=14)

		fig.savefig("./" + output_dir + "/clustering/PCA_impact_variable_graphs.pdf")

		plt.close(fig)


		list_of_PC_scores = []
		list_headers_PC_scores = []
		i = 1
		for PC_score in pca.explained_variance_ratio_ :
			list_headers_PC_scores.append("PC " + str(i))
			list_of_PC_scores.append(PC_score.tolist())
			i += 1

		variance_explained_PCA = pd.DataFrame(list_of_PC_scores, columns = ["%"], index = list_headers_PC_scores)
		variance_explained_PCA.to_csv("./" + output_dir + "/clustering/variance_explained_PCA.txt", sep='\t', encoding="utf-8")

		sns_plot = sns.barplot(x = list_headers_PC_scores , y = list_of_PC_scores)

		total_number_of_PC = float(len(list_of_PC_scores))

		for bar in sns_plot.patches:
			height = bar.get_height()
			sns_plot.text(bar.get_x()+bar.get_width()/2., height + 3, '{:1.2f}'.format(height/total_number_of_PC), ha="center")

		fig = sns_plot.get_figure()

		bar_plot = fig.add_subplot(111)

		plt.suptitle('Variance explained by PC (%)', fontsize = 22)
		bar_plot.set_xlabel('PC component', fontsize = 18)
		bar_plot.set_ylabel('%', fontsize = 18)
		plt.legend(fontsize = 12)
		plt.xticks(fontsize=14, rotation=90)
		plt.yticks(fontsize=14)

		fig.savefig("./" + output_dir + "/clustering/PCA_variance_graphs.pdf")

		plt.close(fig)



		PCA_transform_1 = pca.transform(X)[:,0] # see 'prcomp(my_data)$x' in R
		PCA_transform_2 = pca.transform(X)[:,1]
		PCA_transform_3 = pca.transform(X)[:,2] # see 'prcomp(my_data)$x' in R
		PCA_transform_4 = pca.transform(X)[:,3]

		if isinstance(number_of_clusters, int) : 

			is_numeric = True

			if number_of_clusters > 0 : 

				initial_centers = kmeans_plusplus_initializer(PCA_data, number_of_clusters).initialize();

				kmeans_instance = kmeans(PCA_data, initial_centers)
				kmeans_instance.process()

				clusters = np.array(kmeans_instance.get_clusters())
				centers = np.array(kmeans_instance.get_centers())

			else : 

				xmeans_instance = xmeans(PCA_data);

				xmeans_instance.process()

				clusters = np.array(xmeans_instance.get_clusters())
				centers = np.array(xmeans_instance.get_centers())

				number_of_clusters = len(np.unique(clusters))


		elif isinstance(number_of_clusters, str) : 

			is_numeric = False

			with open(number_of_clusters, 'r') as file_pattern :
				list_pattern = file_pattern.read()
				file_pattern.close()

			split_line = re.compile('\n').split(list_pattern)
			i = 0
			while i < len(split_line) :

				if split_line[i] != '' :
					split_pattern = re.compile(";").split(split_line[i])

				i += 1	

			clusters = [[] for _ in split_pattern]

			pos = 0
			for header_seq in results_for_PCA.index.tolist() :
				pattern_found = False
				i = 0
				for pattern in split_pattern :
					search_pattern = re.match('.*' + split_pattern[i] + ".*", header_seq)
					if search_pattern != None :
						clusters[i].append(pos)
						pattern_found = True

						break

					i += 1

				if pattern_found == False : 
					if len(clusters) == len(split_pattern) and split_pattern[-1] != "No group" : 
						clusters.append([])
						split_pattern.append("No group")
					clusters[-1].append(pos)

				pos += 1

			number_of_clusters = len(clusters)

			clusters = np.array(clusters)

			cluster_list_pattern = [None] * len(PCA_data)

			i = 0

			while i < len(clusters) :
				j = 0
				while j < len(clusters[i]) : 
					cluster_list_pattern[clusters[i][j]] = split_pattern[i]
					j += 1
				i += 1

		cluster_list = [None] * len(PCA_data)

		i = 0

		while i < len(clusters) :
			j = 0
			while j < len(clusters[i]) : 
				cluster_list[clusters[i][j]] = i
				j += 1
			i += 1

		i = 0
		liste_headers_for_cluster = []
		liste_clusters = []

		for header in results_for_PCA.index :
			liste_headers_for_cluster.append(str(header))
			liste_clusters.append(str(cluster_list[i]))
			i+= 1

		i = 0
		liste_plot = []
		liste_legends = []
		while i < number_of_clusters :
			liste_plot.append("Cluster " + str(i))
			i += 1

		if is_numeric == True : 
			clustering_results_for_csv = pd.DataFrame(cluster_list, columns=["Cluster group"], index = results_for_PCA.index)
			clustering_results_for_csv.to_csv("./" + output_dir + "/clustering/clustering_results.txt", sep='\t', encoding="utf-8")

		else : 
			clustering_results_for_csv = pd.DataFrame(cluster_list_pattern, columns=["Cluster group"], index = results_for_PCA.index)
			clustering_results_for_csv.to_csv("./" + output_dir + "/clustering/clustering_results.txt", sep='\t', encoding="utf-8")
		clustering_results = pd.DataFrame( np.column_stack([results_for_PCA.index, PCA_data[:, 0], PCA_data[:, 1], PCA_data[:, 2], PCA_data[:, 3], cluster_list]), columns=['Header', 'PCA1', 'PCA2', 'PCA3', 'PCA4', "Cluster group"])

		i = 0

		# print (clustering_results)
		# while i < number_of_clusters :
		# 	current_cluster = clustering_results.loc[clustering_results['Cluster group'] == i]
		# 	if current_cluster.empty == True : 
		# 		print ("We found an empty cluster during clustering cleaning step. This shouldn't happen. Check if all your data has been analyzed!")
		# 		clustering_results.drop(clustering_results['Cluster group'] == i)
		# 	i += 1 
		
		cm = plt.get_cmap('gist_rainbow')

		NUM_COLORS = number_of_clusters

		x = np.arange(10)
		number_of_colors = [i+x+(i*x)**2 for i in range(number_of_clusters)]

		colors = matplotlib.cm.rainbow(np.linspace(0, 1, len(number_of_colors)))

		lim_min_PCA1 = (min(min(clustering_results['PCA1']), min([value * max(PCA_transform_1) * 1.1 for value in pca_vector_1]))) - (abs(0.1*(min(min(clustering_results['PCA1']), min([value * max(PCA_transform_1) * 1.1 for value in pca_vector_1])))))
		lim_min_PCA2 = (min(min(clustering_results['PCA2']), min([value * max(PCA_transform_2) * 1.1 for value in pca_vector_2]))) - (abs(0.1*(min(min(clustering_results['PCA2']), min([value * max(PCA_transform_2) * 1.1 for value in pca_vector_2])))))
		lim_min_PCA3 = (min(min(clustering_results['PCA3']), min([value * max(PCA_transform_3) * 1.1 for value in pca_vector_3]))) - (abs(0.1*(min(min(clustering_results['PCA3']), min([value * max(PCA_transform_3) * 1.1 for value in pca_vector_3])))))
		lim_min_PCA4 = (min(min(clustering_results['PCA4']), min([value * max(PCA_transform_4) * 1.1 for value in pca_vector_4]))) - (abs(0.1*(min(min(clustering_results['PCA4']), min([value * max(PCA_transform_4) * 1.1 for value in pca_vector_4])))))

		lim_max_PCA1 = (max(max(clustering_results['PCA1']), max([value * max(PCA_transform_1) * 1.1 for value in pca_vector_1]))) + (abs(0.1*(max(max(clustering_results['PCA1']), max([value * max(PCA_transform_1) * 1.1 for value in pca_vector_1])))))
		lim_max_PCA2 = (max(max(clustering_results['PCA2']), max([value * max(PCA_transform_2) * 1.1 for value in pca_vector_2]))) + (abs(0.1*(max(max(clustering_results['PCA2']), max([value * max(PCA_transform_2) * 1.1 for value in pca_vector_2])))))
		lim_max_PCA3 = (max(max(clustering_results['PCA3']), max([value * max(PCA_transform_3) * 1.1 for value in pca_vector_3]))) + (abs(0.1*(max(max(clustering_results['PCA3']), max([value * max(PCA_transform_3) * 1.1 for value in pca_vector_3])))))
		lim_max_PCA4 = (max(max(clustering_results['PCA4']), max([value * max(PCA_transform_4) * 1.1 for value in pca_vector_4]))) + (abs(0.1*(max(max(clustering_results['PCA4']), max([value * max(PCA_transform_4) * 1.1 for value in pca_vector_4])))))

		fig = plt.figure()
		cluster_plot = fig.add_subplot(111)

		i = 0

		print (number_of_clusters)

		while i < number_of_clusters :
			current_cluster = clustering_results.loc[clustering_results['Cluster group'] == i]
			if is_numeric == True : 
				plt.scatter(current_cluster["PCA1"].tolist(), current_cluster["PCA2"].tolist(), color = colors[i], s=30, label="Cluster " + str(i))
				plt.scatter(centers[:, 0], centers[:, 1], c='black', s=50, alpha=0.5, marker="+")
			else : 
				plt.scatter(current_cluster["PCA1"].tolist(), current_cluster["PCA2"].tolist(), color = colors[i], s=30, label=str(np.unique(cluster_list_pattern)[i]))

			i += 1

		for j in range(len(pca_vector_1)):
			plt.text(pca_vector_1[j]*max(PCA_transform_1)*1.1, pca_vector_2[j]*max(PCA_transform_2)*1.1, list(results_for_PCA.columns.values)[j], color='black', fontsize=14, fontweight = "bold")



		length_lines_PCA_1 = max(abs(min([value * max(PCA_transform_1) * 1.1 for value in pca_vector_1])) - abs(0.05*min([value * max(PCA_transform_1) * 1.1 for value in pca_vector_1])), abs(max([value * max(PCA_transform_1) * 1.1 for value in pca_vector_1])) + abs(0.05*max([value * max(PCA_transform_1) * 1.1 for value in pca_vector_1])))
		length_lines_PCA_2 = max(abs(min([value * max(PCA_transform_2) * 1.1 for value in pca_vector_2])) - abs(0.05*min([value * max(PCA_transform_2) * 1.1 for value in pca_vector_2])), abs(max([value * max(PCA_transform_2) * 1.1 for value in pca_vector_2])) + abs(0.05*max([value * max(PCA_transform_2) * 1.1 for value in pca_vector_2])))
		length_lines_PCA_3 = max(abs(min([value * max(PCA_transform_3) * 1.1 for value in pca_vector_3])) - abs(0.05*min([value * max(PCA_transform_3) * 1.1 for value in pca_vector_3])), abs(max([value * max(PCA_transform_3) * 1.1 for value in pca_vector_3])) + abs(0.05*max([value * max(PCA_transform_3) * 1.1 for value in pca_vector_3])))
		length_lines_PCA_4 = max(abs(min([value * max(PCA_transform_4) * 1.1 for value in pca_vector_4])) - abs(0.05*min([value * max(PCA_transform_4) * 1.1 for value in pca_vector_4])), abs(max([value * max(PCA_transform_4) * 1.1 for value in pca_vector_4])) + abs(0.05*max([value * max(PCA_transform_4) * 1.1 for value in pca_vector_4])))

		plt.plot([ - (length_lines_PCA_1 / 2) , (length_lines_PCA_1 / 2)], [0,0], linestyle = '--', color = "black")
		plt.plot([0,0], [- (length_lines_PCA_2 / 2) , (length_lines_PCA_2 / 2)], linestyle = '--', color = "black")

		plt.suptitle('K-means analysis on PCA results', fontsize = 22)
		cluster_plot.set_xlabel('PCA1', fontsize = 18)
		cluster_plot.set_ylabel('PCA2', fontsize = 18)
		cluster_plot.set_xlim(lim_min_PCA1,lim_max_PCA1)
		cluster_plot.set_ylim(lim_min_PCA2,lim_max_PCA2)
		plt.legend(fontsize = 12)
		plt.xticks(fontsize=14)
		plt.yticks(fontsize=14)

		box = cluster_plot.get_position() # get position of figure
		cluster_plot.set_position([box.x0, box.y0, box.width * 0.90, box.height])
		plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
		plt.setp(cluster_plot.get_legend().get_texts(), fontsize='16')
		fig.savefig("./" + output_dir + "/clustering/PCA_1_2_clustering_graph.pdf")
		plt.close(fig)


		fig = plt.figure()
		cluster_plot = fig.add_subplot(111)


		i = 0

		while i < number_of_clusters :
			current_cluster = clustering_results.loc[clustering_results['Cluster group'] == i]
			if is_numeric == True : 
				plt.scatter(current_cluster["PCA1"].tolist(), current_cluster["PCA3"].tolist(), color = colors[i], s=30, label="Cluster " + str(i))
				plt.scatter(centers[:, 0], centers[:, 1], c='black', s=50, alpha=0.5, marker="+")
			else : 
				plt.scatter(current_cluster["PCA1"].tolist(), current_cluster["PCA3"].tolist(), color = colors[i], s=30, label=str(np.unique(cluster_list_pattern)[i]))

			i += 1

		for j in range(len(pca_vector_1)):
			plt.text(pca_vector_1[j]*max(PCA_transform_1)*1.1, pca_vector_3[j]*max(PCA_transform_3)*1.1, list(results_for_PCA.columns.values)[j], color='black', fontsize=14, fontweight = "bold")

		plt.plot([ - (length_lines_PCA_1 / 2) , (length_lines_PCA_1 / 2)], [0,0], linestyle = '--', color = "black")
		plt.plot([0,0], [- (length_lines_PCA_3 / 2) , (length_lines_PCA_3 / 2)], linestyle = '--', color = "black")

		plt.suptitle('K-means analysis on PCA results', fontsize = 22)
		cluster_plot.set_xlabel('PCA1', fontsize = 18)
		cluster_plot.set_ylabel('PCA3', fontsize = 18)
		cluster_plot.set_xlim(lim_min_PCA1,lim_max_PCA1)
		cluster_plot.set_ylim(lim_min_PCA3,lim_max_PCA3)
		plt.legend(fontsize = 12)
		plt.xticks(fontsize=14)
		plt.yticks(fontsize=14)
		box = cluster_plot.get_position() # get position of figure
		cluster_plot.set_position([box.x0, box.y0, box.width * 0.90, box.height])
		plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
		plt.setp(cluster_plot.get_legend().get_texts(), fontsize='16')
		fig.savefig("./" + output_dir + "/clustering/PCA_1_3_clustering_graph.pdf")
		plt.close(fig)


		fig = plt.figure()
		cluster_plot = fig.add_subplot(111)



		i = 0

		while i < number_of_clusters :
			current_cluster = clustering_results.loc[clustering_results['Cluster group'] == i]
			if is_numeric == True : 
				plt.scatter(current_cluster["PCA1"].tolist(), current_cluster["PCA4"].tolist(), color = colors[i], s=30, label="Cluster " + str(i))
				plt.scatter(centers[:, 0], centers[:, 1], c='black', s=50, alpha=0.5, marker="+")
			else : 
				plt.scatter(current_cluster["PCA1"].tolist(), current_cluster["PCA4"].tolist(), color = colors[i], s=30, label=str(np.unique(cluster_list_pattern)[i]))

			i += 1

		for j in range(len(pca_vector_1)):
			plt.text(pca_vector_1[j]*max(PCA_transform_1)*1.1, pca_vector_4[j]*max(PCA_transform_4)*1.1, list(results_for_PCA.columns.values)[j], color='black', fontsize=14, fontweight = "bold")

		plt.plot([ - (length_lines_PCA_1 / 2) , (length_lines_PCA_1 / 2)], [0,0], linestyle = '--', color = "black")
		plt.plot([0,0], [- (length_lines_PCA_4 / 2) , (length_lines_PCA_4 / 2)], linestyle = '--', color = "black")

		plt.suptitle('K-means analysis on PCA results', fontsize = 22)
		cluster_plot.set_xlabel('PCA1', fontsize = 18)
		cluster_plot.set_ylabel('PCA4', fontsize = 18)
		cluster_plot.set_xlim(lim_min_PCA1,lim_max_PCA1)
		cluster_plot.set_ylim(lim_min_PCA4,lim_max_PCA4)
		plt.legend(fontsize = 12)
		plt.xticks(fontsize=14)
		plt.yticks(fontsize=14)
		box = cluster_plot.get_position() # get position of figure
		cluster_plot.set_position([box.x0, box.y0, box.width * 0.90, box.height])
		plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
		plt.setp(cluster_plot.get_legend().get_texts(), fontsize='16')
		fig.savefig("./" + output_dir + "/clustering/PCA_1_4_clustering_graph.pdf")
		plt.close(fig)


		fig = plt.figure()
		cluster_plot = fig.add_subplot(111)


		i = 0

		while i < number_of_clusters :
			current_cluster = clustering_results.loc[clustering_results['Cluster group'] == i]
			if is_numeric == True : 
				plt.scatter(current_cluster["PCA2"].tolist(), current_cluster["PCA3"].tolist(), color = colors[i], s=30, label="Cluster " + str(i))
				plt.scatter(centers[:, 0], centers[:, 1], c='black', s=50, alpha=0.5, marker="+")
			else : 
				plt.scatter(current_cluster["PCA2"].tolist(), current_cluster["PCA3"].tolist(), color = colors[i], s=30, label=str(np.unique(cluster_list_pattern)[i]))

			i += 1

		for j in range(len(pca_vector_2)):
			plt.text(pca_vector_2[j]*max(PCA_transform_2)*1.1, pca_vector_3[j]*max(PCA_transform_3)*1.1, list(results_for_PCA.columns.values)[j], color='black', fontsize=14, fontweight = "bold")

		plt.plot([ - (length_lines_PCA_2 / 2) , (length_lines_PCA_2 / 2)], [0,0], linestyle = '--', color = "black")
		plt.plot([0,0], [- (length_lines_PCA_3 / 2) , (length_lines_PCA_3 / 2)], linestyle = '--', color = "black")

		plt.suptitle('K-means analysis on PCA results', fontsize = 22)
		cluster_plot.set_xlabel('PCA2', fontsize = 18)
		cluster_plot.set_ylabel('PCA3', fontsize = 18)
		cluster_plot.set_xlim(lim_min_PCA2,lim_max_PCA2)
		cluster_plot.set_ylim(lim_min_PCA3,lim_max_PCA3)
		plt.legend(fontsize = 12)
		plt.xticks(fontsize=14)
		plt.yticks(fontsize=14)
		box = cluster_plot.get_position() # get position of figure
		cluster_plot.set_position([box.x0, box.y0, box.width * 0.90, box.height])
		plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
		plt.setp(cluster_plot.get_legend().get_texts(), fontsize='16')
		fig.savefig("./" + output_dir + "/clustering/PCA_2_3_clustering_graph.pdf")
		plt.close(fig)

		fig = plt.figure()
		cluster_plot = fig.add_subplot(111)


		i = 0

		while i < number_of_clusters :
			current_cluster = clustering_results.loc[clustering_results['Cluster group'] == i]
			if is_numeric == True : 
				plt.scatter(current_cluster["PCA2"].tolist(), current_cluster["PCA4"].tolist(), color = colors[i], s=30, label="Cluster " + str(i))
				plt.scatter(centers[:, 0], centers[:, 1], c='black', s=50, alpha=0.5, marker="+")
			else : 
				plt.scatter(current_cluster["PCA2"].tolist(), current_cluster["PCA4"].tolist(), color = colors[i], s=30, label=str(np.unique(cluster_list_pattern)[i]))

			i += 1

		for j in range(len(pca_vector_2)):
			plt.text(pca_vector_2[j]*max(PCA_transform_2)*1.1, pca_vector_4[j]*max(PCA_transform_4)*1.1, list(results_for_PCA.columns.values)[j], color='black', fontsize=14, fontweight = "bold")

		plt.plot([ - (length_lines_PCA_2 / 2) , (length_lines_PCA_2 / 2)], [0,0], linestyle = '--', color = "black")
		plt.plot([0,0], [- (length_lines_PCA_4 / 2) , (length_lines_PCA_4 / 2)], linestyle = '--', color = "black")

		plt.suptitle('K-means analysis on PCA results', fontsize = 22)
		cluster_plot.set_xlabel('PCA2', fontsize = 18)
		cluster_plot.set_ylabel('PCA4', fontsize = 18)
		cluster_plot.set_xlim(lim_min_PCA2,lim_max_PCA2)
		cluster_plot.set_ylim(lim_min_PCA3,lim_max_PCA3)
		plt.legend(fontsize = 12)
		plt.xticks(fontsize=14)
		plt.yticks(fontsize=14)
		box = cluster_plot.get_position() # get position of figure
		cluster_plot.set_position([box.x0, box.y0, box.width * 0.90, box.height])
		plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
		plt.setp(cluster_plot.get_legend().get_texts(), fontsize='16')
		fig.savefig("./" + output_dir + "/clustering/PCA_2_4_clustering_graph.pdf")
		plt.close(fig)


		fig = plt.figure()
		cluster_plot = fig.add_subplot(111)


		i = 0

		while i < number_of_clusters :
			current_cluster = clustering_results.loc[clustering_results['Cluster group'] == i]
			if is_numeric == True : 
				plt.scatter(current_cluster["PCA3"].tolist(), current_cluster["PCA4"].tolist(), color = colors[i], s=30, label="Cluster " + str(i))
				plt.scatter(centers[:, 0], centers[:, 1], c='black', s=50, alpha=0.5, marker="+")
			else : 
				plt.scatter(current_cluster["PCA3"].tolist(), current_cluster["PCA4"].tolist(), color = colors[i], s=30, label=str(np.unique(cluster_list_pattern)[i]))

			i += 1

		for j in range(len(pca_vector_3)):
			plt.text(pca_vector_3[j]*max(PCA_transform_3)*1.1, pca_vector_4[j]*max(PCA_transform_4)*1.1, list(results_for_PCA.columns.values)[j], color='black', fontsize=14, fontweight = "bold")
		
		plt.plot([ - (length_lines_PCA_3 / 2) , (length_lines_PCA_3 / 2)], [0,0], linestyle = '--', color = "black")
		plt.plot([0,0], [- (length_lines_PCA_4 / 2) , (length_lines_PCA_4 / 2)], linestyle = '--', color = "black")

		plt.suptitle('K-means analysis on PCA results', fontsize = 22)
		cluster_plot.set_xlabel('PCA3', fontsize = 18)
		cluster_plot.set_ylabel('PCA4', fontsize = 18)
		cluster_plot.set_xlim(lim_min_PCA3,lim_max_PCA3)
		cluster_plot.set_ylim(lim_min_PCA4,lim_max_PCA4)
		plt.legend(fontsize = 12)
		plt.xticks(fontsize=14)
		plt.yticks(fontsize=14)
		box = cluster_plot.get_position() # get position of figure
		cluster_plot.set_position([box.x0, box.y0, box.width * 0.90, box.height])
		plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
		plt.setp(cluster_plot.get_legend().get_texts(), fontsize='16')
		fig.savefig("./" + output_dir + "/clustering/PCA_3_4_clustering_graph.pdf")
		plt.close(fig)

		PCA_1_results_by_var = PCA_results_by_var.loc[["PC 1"]].to_dict('r')[0]
		PCA_2_results_by_var = PCA_results_by_var.loc[["PC 2"]].to_dict('r')[0]
		PCA_3_results_by_var = PCA_results_by_var.loc[["PC 3"]].to_dict('r')[0]
		PCA_4_results_by_var = PCA_results_by_var.loc[["PC 4"]].to_dict('r')[0]

		if type_of_analysis == "codons" : 
			create_cut_graph(genetic_code_ref, PCA_1_results_by_var, "PCA_1",  output_dir)
			create_cut_graph(genetic_code_ref, PCA_2_results_by_var, "PCA_2",  output_dir)
			create_cut_graph(genetic_code_ref, PCA_3_results_by_var, "PCA_3",  output_dir)
			create_cut_graph(genetic_code_ref, PCA_4_results_by_var, "PCA_4",  output_dir)

		cluster_group = clustering_results.pop("Cluster group")
		lut = dict(zip(sorted(cluster_group.unique()), colors))
		row_colors = cluster_group.map(lut)
		fig = plt.figure()

		results_for_PCA_normalized = results_for_PCA[results_for_PCA.columns].astype(float)
		results_for_PCA_normalized = (results_for_PCA_normalized-results_for_PCA_normalized.mean())/results_for_PCA_normalized.std()

		results_for_PCA_normalized =  results_for_PCA_normalized.drop(results_for_PCA_normalized.columns[results_for_PCA_normalized.nunique() == 0], axis = 1)

		hierarchical_clustering_graph = sns.clustermap(results_for_PCA_normalized.reset_index(drop=True), metric="euclidean", method="average", cmap="mako", robust=True, xticklabels=True, yticklabels=False, linewidths=0.0001, row_colors=row_colors)
		plt.setp(hierarchical_clustering_graph.ax_heatmap.xaxis.get_majorticklabels(), fontsize=7)

		hierarchical_clustering_graph.savefig("./" + output_dir + "/clustering/hierarchical_clustering_graph.pdf")
		leaf_names = results_for_PCA_normalized.index.tolist()
		linkage_tree_matrix = hierarchy.linkage(results_for_PCA_normalized, method="average", metric="euclidean")

		scipy_tree = hierarchy.to_tree(linkage_tree_matrix, False)
		newick_tree = getNewick(scipy_tree, "", scipy_tree.dist, leaf_names)
		file_newick_tree = open("./" + output_dir + "/clustering/hierarchical_clustering_tree.nwk", "w")
		file_newick_tree.write(newick_tree)
		file_newick_tree.close()

def getNewick(node, newick, parentdist, leaf_names):
    if node.is_leaf():
        return "%s:%.2f%s" % (leaf_names[node.id], parentdist - node.dist, newick)
    else:
        if len(newick) > 0:
            newick = "):%.2f%s" % (parentdist - node.dist, newick)
        else:
            newick = ");"
        newick = getNewick(node.get_left(), newick, node.dist, leaf_names)
        newick = getNewick(node.get_right(), ",%s" % (newick), node.dist, leaf_names)
        newick = "(%s" % (newick)
        return newick

def compare_data(genetic_code_1, genetic_code_2) :

	sequential_order_amino_acid = ["F","L","I","M","V","S","P","T","A","Y","*","H","Q","N","K","D","E","C","W","R","G"]
	sequential_order_codon = ["UUU", "UCU", "UAU", "UGU", "UUC", "UCC", "UAC", "UGC", "UUA", "UCA", "UAA", "UGA", "UUG", "UCG", "UAG", "UGG", "CUU", "CCU", "CAU", "CGU", "CUC", "CCC", "CAC", "CGC", "CUA", "CCA", "CAA", "CGA", "CUG", "CCG", "CAG", "CGG", "AUU", "ACU", "AAU", "AGU", "AUC", "ACC", "AAC", "AGC", "AUA", "ACA", "AAA", "AGA", "AUG", "ACG", "AAG", "AGG", "GUU", "GCU", "GAU", "GGU", "GUC", "GCC", "GAC", "GGC", "GUA", "GCA", "GAA", "GGA", "GUG", "GCG", "GAG", "GGG"]

	sequential_order__data_comparison.append("firstline")
	results_for_file_data_comparison["firstline"] = "### CODONS : EUCLIDIAN DISTANCES ###\n\n"
	occ_aa_genetic_code_1 = {}
	occ_aa_genetic_code_2 = {}
	occ_tot_genetic_code_1 = 0
	occ_tot_genetic_code_2 = 0
	euclidian_distance_codon = {}
	euclidian_distance_amino_acid = {}

	for amino_acid in sequential_order_amino_acid :

		if amino_acid not in occ_aa_genetic_code_1 :
			occ_aa_genetic_code_1[amino_acid] = 0

		if amino_acid not in occ_aa_genetic_code_2 :
			occ_aa_genetic_code_2[amino_acid] = 0

		for codon in sequential_order_codon :

			if codon in genetic_code_1[amino_acid] and codon in genetic_code_2[amino_acid] :

				occ_aa_genetic_code_1[amino_acid] += genetic_code_1[amino_acid][codon][0]
				occ_aa_genetic_code_2[amino_acid] += genetic_code_2[amino_acid][codon][0]
				occ_tot_genetic_code_1 += genetic_code_1[amino_acid][codon][0]
				occ_tot_genetic_code_2 += genetic_code_2[amino_acid][codon][0]

				euclidian_distance_codon[codon] = math.pow((genetic_code_1[amino_acid][codon][1] - genetic_code_2[amino_acid][codon][1]), 2)
				euclidian_distance_codon[codon] = math.sqrt(euclidian_distance_codon[codon])
				sequential_order__data_comparison.append(codon)
				results_for_file_data_comparison[codon] = codon + "\t" + str(round(euclidian_distance_codon[codon],3)) + "\n"

	sequential_order__data_comparison.append("secondline")
	results_for_file_data_comparison["secondline"] = "\n### AMINO ACIDS : EUCLIDIAN DISTANCES ###\n\n"
	for amino_acid in sequential_order_amino_acid :

		if amino_acid in occ_aa_genetic_code_1 and amino_acid in occ_aa_genetic_code_2 :

			euclidian_distance_amino_acid[amino_acid] = math.pow(((occ_aa_genetic_code_1[amino_acid]/occ_tot_genetic_code_1) - (occ_aa_genetic_code_2[amino_acid]/occ_tot_genetic_code_2)), 2)
			euclidian_distance_amino_acid[amino_acid] = math.sqrt(euclidian_distance_amino_acid[amino_acid])
			sequential_order__data_comparison.append(amino_acid)
			results_for_file_data_comparison[amino_acid] = amino_acid + "\t" + str(round(euclidian_distance_amino_acid[amino_acid],3)) + "\n"

	return results_for_file_data_comparison


def perform_wilcoxon_mann_whitney_u_test_simulation_step(dataset_1, dataset_2, csv_name, header, results_file, output_dir, output_dir_wilcox) :
	global firstline

	if firstline :
		csv_fill = "\t" + csv_name
		results_file["firstline"] += csv_fill
		firstline = False    
	wilcox_file = open("./" + output_dir + "/" + output_dir_wilcox + "/" + header[:210] + ".txt", "a")

	alpha = 0.05


	stat, p = mannwhitneyu(dataset_1, dataset_2)
	
	if p > alpha : 
		wilcox_file.write("\n\n" + "##############################" + "\n"+ csv_name + "\n\n" + "test statistics : " + str(stat) + " \n" + "p_value : " + str(p) + "\n")
		results_file[header] += "\t" + "H0 accepted (p-value : " + str(p) + ")"
	
	else :
		wilcox_file.write("\n\n" + "##############################" + "\n"+ csv_name + "\n\n" + "test statistics : " + str(stat) + " \n" + "p_value : " + str(p) + "\n")
		results_file[header] += "\t" + "H0 rejected (p-value : " + str(p) + ")"

def abline(slope, intercept):
	axes = plt.gca()
	x_vals = np.array(axes.get_xlim())
	y_vals = intercept + slope * x_vals
	plt.plot(x_vals, y_vals, '--')


def create_dicodon_dict(genetic_code_ref) : 

	dicodon_dict = {}

	for amino_acid in genetic_code_ref :
		if amino_acid != "*" : 
			for first_codon in genetic_code_ref[amino_acid] : 

				if first_codon not in dicodon_dict : 
					dicodon_dict[first_codon] = {}

				for amino_acid in genetic_code_ref : 
					for second_codon in genetic_code_ref[amino_acid]:
						if second_codon not in dicodon_dict[first_codon] : 
							dicodon_dict[first_codon][second_codon] = None

	return dicodon_dict		

def create_diamino_acid_dict(genetic_code_ref) : 

	diamino_acid_dict = {}

	for first_amino_acid in genetic_code_ref :
		if first_amino_acid != "*" : 
			if first_amino_acid not in diamino_acid_dict : 
				diamino_acid_dict[first_amino_acid] = {}
				for second_amino_acid in genetic_code_ref : 
					dicodon_dict[first_amino_acid][second_amino_acid] = None

	return dicodon_dict							


def create_cut_graph(cut_dict, gradient_variable, name_of_cut_table, output_dir) :


	sequential_order_cu_table = ["UUU", "UCU", "UAU", "UGU", "UUC", "UCC", "UAC", "UGC", "UUA", "UCA", "UAA", "UGA", "UUG", "UCG", "UAG", "UGG", "CUU", "CCU", "CAU", "CGU", "CUC", "CCC", "CAC", "CGC", "CUA", "CCA", "CAA", "CGA", "CUG", "CCG", "CAG", "CGG", "AUU", "ACU", "AAU", "AGU", "AUC", "ACC", "AAC", "AGC", "AUA", "ACA", "AAA", "AGA", "AUG", "ACG", "AAG", "AGG", "GUU", "GCU", "GAU", "GGU", "GUC", "GCC", "GAC", "GGC", "GUA", "GCA", "GAA", "GGA", "GUG", "GCG", "GAG", "GGG"]

	cm = sns.color_palette("coolwarm", 64)

	df = pd.DataFrame()

	df["1_col"] = ["UUU" , "UUC" , "UUA" , "UUG" , "CUU" , "CUC" , "CUA" , "CUG" , "AUU" , "AUC" , "AUA" , "AUG" , "GUU" , "GUC" , "GUA" , "GUG"]
	df["1_col_aa"] = [None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None]
	df["1_col_var"] = [None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None]

	df["2_col"] = ["UCU" , "UCC" , "UCA" , "UCG" , "CCU" , "CCC" , "CCA" , "CCG" , "ACU" , "ACC" , "ACA" , "ACG" , "GCU" , "GCC" , "GCA" , "GCG"]
	df["2_col_aa"] = [None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None]
	df["2_col_var"] = [None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None]

	df["3_col"] = ["UAU" , "UAC" , "UAA" , "UAG" , "CAU" , "CAC" , "CAA" , "CAG" , "AAU" , "AAC" , "AAA" , "AAG" , "GAU" , "GAC" , "GAA" , "GAG"]
	df["3_col_aa"] = [None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None]
	df["3_col_var"] = [None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None]

	df["4_col"] = ["UGU" , "UGC" , "UGA" , "UGG" , "CGU" , "CGC" , "CGA" , "CGG" , "AGU" , "AGC" , "AGA" , "AGG" , "GGU" , "GGC" , "GGA" , "GGG"]
	df["4_col_aa"] = [None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None]
	df["4_col_var"] = [None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None]

	cpt = 0
	cpt_fill_1 = 0
	cpt_fill_2 = 0
	cpt_fill_3 = 0
	cpt_fill_4 = 0

	for codon in sequential_order_cu_table : 

		cpt += 1

		for amino_acid in cut_dict : 

			if codon in cut_dict[amino_acid] : 
				
				if cpt == 1 : 

					df["1_col_aa"][cpt_fill_1] = amino_acid
					if codon in gradient_variable : 
						df["1_col_var"][cpt_fill_1] = round(gradient_variable[codon],3)
					else : 
						df["1_col_var"][cpt_fill_1] = 0
					cpt_fill_1 += 1

				elif cpt == 2 : 
					

					df["2_col_aa"][cpt_fill_2] = amino_acid
					if codon in gradient_variable : 
						df["2_col_var"][cpt_fill_2] = round(gradient_variable[codon],3)
					else : 
						df["2_col_var"][cpt_fill_2] = 0
					cpt_fill_2 += 1

				elif cpt == 3 : 

					df["3_col_aa"][cpt_fill_3] = amino_acid
					if codon in gradient_variable : 
						df["3_col_var"][cpt_fill_3] = round(gradient_variable[codon],3)
					else : 
						df["3_col_var"][cpt_fill_3] = 0
					cpt_fill_3 += 1

				elif cpt == 4 : 
					
					df["4_col_aa"][cpt_fill_4] = amino_acid
					if codon in gradient_variable : 
						df["4_col_var"][cpt_fill_4] = round(gradient_variable[codon],3)
					else : 
						df["4_col_var"][cpt_fill_4] = 0
					cpt = 0
					cpt_fill_4 += 1


	dataframe_values = pd.DataFrame()
	dataframe_values["col1"] = df["1_col_var"]
	dataframe_values["col2"] = df["2_col_var"]
	dataframe_values["col3"] = df["3_col_var"]
	dataframe_values["col4"] = df["4_col_var"]

	dataframe_values = dataframe_values.astype(float)


	df["1_col"] = [str(element) + " " + "(" for element in df["1_col"]]
	df["1_col_aa"] = [str(element) + ")" + "\t" for element in df["1_col_aa"]]
	df["1_col_var"] = [str(element) for element in df["1_col_var"]]

	df["2_col"] = [str(element) +  " " + "(" for element in df["2_col"]]
	df["2_col_aa"] = [str(element) +  ")" + "\t" for element in df["2_col_aa"]]
	df["2_col_var"] = [str(element) for element in df["2_col_var"]]

	df["3_col"] = [str(element) +  " " + "(" for element in df["3_col"]]
	df["3_col_aa"] = [str(element) +  ")" + "\t" for element in df["3_col_aa"]]
	df["3_col_var"] = [str(element) for element in df["3_col_var"]]

	df["4_col"] = [str(element) +  " " + "(" for element in df["4_col"]]
	df["4_col_aa"] = [str(element) +  ")" + "\t" for element in df["4_col_aa"]]
	df["4_col_var"] = [str(element) for element in df["4_col_var"]]

	dataframe_table = pd.DataFrame()

	dataframe_table["col1"] = df['1_col'].str.cat(df['1_col_aa'], sep ="")
	dataframe_table["col1"] = dataframe_table['col1'].str.cat(df['1_col_var'], sep ="")

	dataframe_table["col2"] = df['2_col'].str.cat(df['2_col_aa'], sep ="")
	dataframe_table["col2"] = dataframe_table['col2'].str.cat(df['2_col_var'], sep ="")

	dataframe_table["col3"] = df['3_col'].str.cat(df['3_col_aa'], sep ="")
	dataframe_table["col3"] = dataframe_table['col3'].str.cat(df['3_col_var'], sep ="")

	dataframe_table["col4"] = df['4_col'].str.cat(df['4_col_aa'], sep ="")
	dataframe_table["col4"] = dataframe_table['col4'].str.cat(df['4_col_var'], sep ="")

	dataframe_table.columns = [''] * len(dataframe_table.columns)
	dataframe_table.reset_index(drop=True, inplace=True)

	dataframe_values.columns = [''] * len(dataframe_values.columns)
	dataframe_values.reset_index(drop=True, inplace=True)


	sns_plot = sns.heatmap(dataframe_values, annot=dataframe_table, fmt="", cmap='RdYlGn', cbar=False, annot_kws={"size": 12})
	sns_plot.tick_params(left=False, bottom=False)
	sns_plot.set_xticks([])
	sns_plot.set_yticks([])
	fig = sns_plot.get_figure()
	fig.savefig(output_dir + "/codon_usage_tables/" + str(name_of_cut_table) + ".pdf")
	plt.close(fig)


def seq_to_dicodon_dict(sequence, header, dicodon_dict, output_dir) :

	dicodon_dict_seq = deepcopy(dicodon_dict)

	i = 0
	length_seq = len(sequence) / 3
	cpt = 0
	while i < length_seq - 1 : # -1 because there is no need to search for a dicodon on the last codon.
		dicodon_seq = sequence[:6] # We take the two codons
		
		first_codon = dicodon_seq[:int(len(dicodon_seq)/2)]
		second_codon = dicodon_seq[int(len(dicodon_seq)/2):]

		if dicodon_dict_seq[first_codon][second_codon] == None : 
			dicodon_dict_seq[first_codon][second_codon] = 1
		else :
			dicodon_dict_seq[first_codon][second_codon] += 1

		sequence  = sequence[3:] # We suppress the first codons, since it won't be analyzed anymore
		i += 1

	for first_codon in dicodon_dict :
		for second_codon in dicodon_dict[first_codon] : 

			if dicodon_dict_seq[first_codon][second_codon] == None : 
				dicodon_dict_seq[first_codon][second_codon] = []
				dicodon_dict_seq[first_codon][second_codon].append(0)
				dicodon_dict_seq[first_codon][second_codon].append(0) # get the frequency of this dicodon
			else : 
				value = dicodon_dict_seq[first_codon][second_codon]
				dicodon_dict_seq[first_codon][second_codon] = []
				dicodon_dict_seq[first_codon][second_codon].append(value)
				dicodon_dict_seq[first_codon][second_codon].append(value/(length_seq - 1)) # get the frequency of this dicodon

	return dicodon_dict_seq
