#!/usr/bin/env python3
# -*-coding: utf-8-*-

#MODULES###########################################################################

import os
import re
import sys

from Bio import SeqIO
from Bio.ExPASy import Prosite, Prodoc


def translator(domain):

    """
    This function changes the regular expression format from Prosite to the Python
    one of the string given. 
    """

    dom = (str(domain).upper().replace(".", "").replace("-", "").replace("X", "."
		).replace("{", "[^").replace("}", "]").replace("<", "^").replace(">", "$"
		).replace("(", "{").replace(")", "}"))
    return dom


def results(output_file, id, name, accession, description, pattern, sequence):

    """
    This function writes in the output file the results if there is any match 
    between the domain and the sequence. 
    Arguments: output file; domain ID, name, accession, description, pattern;
    and the sequence to analyze in each case. 
    """

    sys.stdout = open(output_file, 'a')

    print("\nRESULTS "+id + ":")
    print("------------------------------------------------------------------")
    print("OBTAINED DOMAIN: ")
    print("Name: "+name)
    print("Accesion: "+accession)
    print("Description: "+description)
    print("Found pattern: "+pattern)
    print("Protein's domain: " +
          str(re.findall(pattern, sequence)[0]).replace('"', ""))
    match = re.search(pattern, sequence)
    print("Domain position: " + str(match.start())+"-"+str(match.end()))
    print("------------------------------------------------------------------\n")


def no_results(output_file, id):

    """
    This function prints the results if there isn't any match between the domain
	and the sequence. 
    Arguments: output file name, domain ID
    """

    sys.stdout = open(output_file, 'a')

    print("\nRESULTS "+id + ":")
    print("---------------------------------------------------------------------")
    print("No matching domain found.")
    print("---------------------------------------------------------------------")


def domains_finder(input_file, database, output_file):
	
    """
    This function chacks sequence by sequence if there is any matched domain
	present in the prosite database (.dat). 
    Arguments: file which contais the sequences, database, output_file name. 
    """

    check = 0

    with open(input_file) as handle2:
        for recordBlast in SeqIO.parse(handle2, "fasta"):
            check = 0

            # Variable which contains the sequence
            subject_seqs = str(recordBlast.seq)

            with open(database, 'r') as handle:
                for recordD in Prosite.parse(handle):

                    # Variable which contains the "translated" domain
                    a = translator(recordD.pattern)

                    if re.search(a, subject_seqs) and recordD.pattern != "":

                        results(output_file, recordBlast.id, recordD.name,
                                recordD.accession, recordD.description, a, 
								subject_seqs)
                        check = 1

                if check == 0:
                    no_results(output_file, recordBlast.id)

    sys.stdout.close()
    sys.stdout = open("/dev/stdout", "w")
