#!/usr/bin/env python3
#-*-coding: utf-8-*-

#MODULES###########################################################################

import shutil
import sys
import os

from Bio import SeqIO

#FUNCTIONS########################################################################

def fasta_conversion(input_file, output_file):

	"""
	Function whose arguments are an input file and an output file. 
	It turns the genbank file (input_file) into a fasta one (output_file).
	"""

	with open(input_file, "r") as input_handle:
		for record in SeqIO.parse(input_handle, "genbank"):
			for feature in record.features:

				#We take the locus tag, the accession and the protein sequence if
				# this exists. 

				if feature.type == 'CDS':
					try:						
						if feature.qualifiers['translation'][0] != " ":
							sys.stdout=open(output_file,'a')
							print (">"+feature.qualifiers['locus_tag'][0]+"@"+ 
								record.name)
							print(feature.qualifiers['translation'][0])
						
					except:
						pass

		sys.stdout.close()
		sys.stdout = open("/dev/stdout", "w")		
							

def fasta(input_file, output_file):

	"""
	Function whose arguments are an input file and an output file. 
	If the input file has fasta format, it will be copied into the output file.
	"""
		
	f = open(input_file, 'r')
	content = f.read()
	if content.startswith(">"):
		shutil.copy(input_file, output_file)
	f.close()


	