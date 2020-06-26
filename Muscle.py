#!/usr/bin/env python3
#-*-coding: utf-8-*-

#MODULES###########################################################################

import os
import shutil
import sys

from Bio import SeqIO
from Bio import Phylo

#FUNCTIONS##########################################################################


def Muscle(input_file, output_file):

	"""
	This function generates an alignment between the sequences present in the BlastP
	results filtered.This alignment is introduced in the output file. 
	Arguments: BlastP results filtered file and output file name.
	"""

	muscle="muscle -in "+input_file+" -out "+output_file+" 2>/dev/null "
	os.system(muscle)


def tree(input_file, output_file):

	"""
	This function generates a tree with the sequences present in the alignment.
	This tree created is introduced in the output file. 
	Arguments: alignment and output file name. 
	"""

	tree="muscle -maketree -in "+input_file+" -out "+output_file+ \
		" -cluster neighborjoining "+ " 2>/dev/null "
	os.system(tree)


def Phylo_maketree(input_file, output_file):

	"""
	This function creates a tree that is more visual than the previous one. This new
	tree is stored in the output file. 
	Arguments: tree in newick format, name of the output file.
	"""

	sys.stdout = open(output_file, 'w')
	tree=Phylo.read(input_file, 'newick')
	Phylo.draw_ascii(tree)
	sys.stdout.close()
	sys.stdout = open("/dev/stdout", "w")	
		



	
		
	




