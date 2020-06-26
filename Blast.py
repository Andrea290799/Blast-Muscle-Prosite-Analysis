#!/usr/bin/env python3
# -*-coding: utf-8-*-

#MODULES###########################################################################

from subprocess import call
import os
import shutil
import sys

from Bio import SeqIO


#FUNCTIONS########################################################################


def default(arg_num, value):

    """
    This funtion returns a value introduced as an argument. If it has not been 
    introduced, it returns a default value. 
    Arguments: script argument number, default value.
    """

    if sys.argv[arg_num].isdigit():
        variable = str(sys.argv[arg_num])
    elif (not sys.argv[arg_num].isdigit() and sys.argv[arg_num] != "-h" 
		and sys.argv[arg_num] != "--help"):
        variable = str(value)
    return variable


def Blast(query, subject, evalue, ID, output_file):

    """
    This funtion makes a BlastP analysis of the sequence(s) introduced onto a
	database also introduced. It returns the result inside an output_file. 
    Arguments: file with the sequences, database (both files in fasta format),
	e-value, the name of the query, output file. 
    """

    comands_arguments_blast = "blastp -query " + query + " -subject " + subject + \
        ' -evalue '+evalue + ' -outfmt "6 sseqid pident qcovs evalue sseq " -out '+ \
		output_file
    os.system(comands_arguments_blast)


def blast_results_to_fasta(input, output, identity, coverage):

	"""
	This function filters the BlastP results according its identity and its
	coverage. The results are stored in the output file, with fasta format. 
	Arguments: BlastP results' file, the output file, identity and coverage. 
	"""

	file = open(input, 'r')

	for line in file.readlines():
		fields = line.rstrip().split("\t")
		name = fields[0]
		id = fields[1]
		cov = fields[2]
		seq = fields[4]

		sys.stdout = open(output, 'a')

		if id >= identity and cov >= coverage:
			print(">"+name+"\n"+seq)
			sys.stdout.close()
		sys.stdout = open("/dev/stdout", "w")


def subject_selection(total_subjects, subjects_blast, output_file, ID):

    """
    This function takes the resulting subjects of BlastP (they are incomplete) and 
	searchs for the complete ones in the initial subject file that was used for 
	doing BlastP.
    Arguments: file which has all the sequences of the genbanks, BlastP results 
	filtered and transformed to fasta format, output file and query ID.
    """

    file_to_prosite = open(output_file, 'a')

    with open(total_subjects, 'r') as multifasta:
        for recordM in SeqIO.parse(multifasta, "fasta"):
            with open(subjects_blast, 'r') as blast_def:
                for recordB in SeqIO.parse(blast_def, "fasta"):

                    # We check if the sequence IDs are the same, and we take the 
					# sequence and its ID from the file where there are all the 
					# sequences extracted from the genbanks.

                    if str(recordB.id) == str(recordM.id):
                        file_to_prosite.write(
                            ">"+str(recordM.id)+"\n"+str(recordM.seq)+"\n")
    file_to_prosite.close()
