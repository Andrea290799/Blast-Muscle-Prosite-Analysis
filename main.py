#!/usr/bin/env python3
# -*-coding: utf-8-*-
# Andrea Escolar Peña 2020


#MODULES###########################################################################

from datetime import datetime
from datetime import date
import os
import shutil
import sys

from Bio import SeqIO

import To_fasta_conversion as Tfc
import Blast
import Muscle
import PROSITE
import Blast_Graph as BG
import Prosite_Graph as PG

#FUNCTIONS########################################################################


def help():

    """
    This function prints usage when called.
    """

    print("USAGE")
    print("----------------------------------------------------------------------")
    print(sys.argv[0]+" [query_file] [genbanks_directory/] [identity] [coverage]" +
          " [evalue] [prosite_file]\n")
    print("\033[1;01mprosite_file\033[00m -> file that contains prosite database to analize." +
          " It has to be located in the same directory as this script.")
    print("\033[1;01mquery_file\033[00m -> file that contains query sequences to analize." +
          " It has to be located in the same directory as this script.")
    print("\033[1;01mgenbanks_directory\033[00m -> directory that contains genbanks files to " +
          " analize. It is necessary to include a '/' after de directory's name")
    print("It has to be located in the same directory as this script.\n")
    print("-----For using the default value, introduce a non numeric character.-----")
    print("\033[1;01midentity\033[00m -> minimum identity percentage required. 0 <= number <= 100.")
    print("\033[1;01mcoverage\033[00m -> minimum coverage percentage required. 0 <= number <= 100.")
    print("\033[1;01mevalue\033[00m -> minimum evalue percentage required. 0 <= number <= 100.\n")
    print("----------------------------------------------------------------------")


def detailed_help():

    """
    This function prints more detailed information of usage when called.
    """

    print("This program will be used to analize some query proteins using " +
          "BlastP; to obtain a MSA and a phylogenetic trees (N-J) using "+
          "MUSCLE and to find proteins' domains from Prosite database "+
          "that are present in the BlasP results. Results from second and"+
          " third analysis are also graphicated.\n ")
    print("A results folder will be created. The results of the different"+
    " queries are organized in their corresponding folder,")
    print("regarding their names. The data used is storaged in the DATA "+
    "folder, inside the results folder generated.")
    print("Results folder: Blast_Muscle_Prosite_results_[date_time].\n")
    help()
    print("Write "+sys.argv[0]+" [-h] or " +
          sys.argv[0]+" [--help] to print this message.\n")


def arguments():

    """
    This function detects whether the number of arguments introduced is 
    correct or not.
    """

    if len(sys.argv) != 7:
        print("Wrong number of arguments.")
        print("Arguments introduced: "+str(len(sys.argv)-1))
        print("Arguments to introduce: 6.\n")

        help()

        sys.exit()


def ask_for_help():

    """
    This funtion prints detailed help if called.
    """

    for e in sys.argv:
        if e == "-h" or e == "--help":
            
            detailed_help()
            
            sys.exit()


def file_exists(file):

    """
    This function checks whether the input file exist.
    Argument: input file.
    """

    if not os.path.isfile(file):
        print("Input file "+file+" does not exist.\n")
       
        help()
        
        sys.exit()


def directory_exists(directory):

    """
    This function checks whether the input directory exist.
    Argument: input directory.
    """

    if not os.path.isdir(directory):
        print("Input directory "+directory+" does not exist.\n")
       
        help()
        
        sys.exit()


def file_format(file, element):

    """
    This function checks if the input file has an appropiate format.
    Arguments: input file, element that defines the format.
    """

    f = open(file, 'r').read()

    if not f.startswith(element):
        print("Input file "+file+" does not have an appropiate format.\n")
       
        help()
       
        sys.exit()


def number_range(number):

    """
    This function checks if the introduced number is positive.
    """

    tmp_number = number.replace("-", "").replace(".", "", 1)

    if str(tmp_number).isdigit():
        number = float(number)
        if number < 0:
            number = str(number)
            print("Introduced number "+number+"  is not positive.\n")
            help()
            sys.exit()
        elif number > 100:
            number = str(number)
            print("Introduced number "+number+" > 100.\n")
            help()
            sys.exit()


def how_many_sequences(input_file):

    """
    This function returns the number of sequences in the input_file.
    """

    i = 0
    with open(input_file, "r") as input_handle:
        for record in SeqIO.parse(input_handle, "fasta"):
            i += 1
    return i


def only_gbks(filenames_list, dirpath):

    """
    This function creates a list of genbakns present in a directory given a 
    list of files in this directory and its path. It also informs about if
    some file does not have genbank format. 
    """

    GBK_list = []

    for element in filenames_list:

        GBK = str(dirpath)+element

        gbk = open(str(dirpath)+element, 'r').read()

        if gbk.startswith("LOCUS"):
            GBK_list.append(element)
        else:
            print("\t\033[0;31mThe file "+GBK +
                  " does not have genbank format, so it won't be used for the analysis.\033[00m\n")

    return GBK_list

def end(identity, coverage, evalue, date_time):

    """
    This function prints this message when the programme is finished. 
    Arguments: identity, coverage, date and time variable
    """

    print("\nThe total analysis has succesfully finished!\n")
    print("Parametres used:\n")
    print("\t-Identity: "+identity)
    print("\t-Coverage: "+coverage)
    print("\t-Evalue: "+evalue+"\n")

    print("All the results and data are stored in Blast_Muscle_Prosite_results_" +
          date_time+" in the current directory.\n")


def main():

    """
    This program will be used to analize some query proteins using BlastP; to 
    obtain a MSA and a phylogenetic trees (N-J) using MUSCLE and Phylo and to 
    find proteins' domains from Prosite database that are present in the BlasP
    results. Results from second and third analysis are also graphicated.
    A results folder will be created. The results of the different queries are 
    organized in their corresponding folder, regarding their names. The data used
    is storaged in the DATA folder, inside the results folder generated.
    """

    print("")

    # We check if all is coorect in order to continue
    ask_for_help()
    arguments()
    file_exists(sys.argv[1])
    directory_exists(sys.argv[2])
    number_range(sys.argv[3])
    number_range(sys.argv[4])
    number_range(sys.argv[5])
    file_exists(sys.argv[6])
    file_format(sys.argv[1], ">")
    file_format(sys.argv[6], "CC")

    # We create the directories where data and results will be storaged
    now = datetime.now()
    time = (format(now.hour) + ":" +
            format(now.minute) + ":" + format(now.second))

    date_time = str(now.date())+"@"+str(time)
    data = "Blast_Muscle_Prosite_results_"+date_time+"/DATA"
    os.makedirs(data)

    querys_input = sys.argv[1]
    GBKs_directory = sys.argv[2]

    # We copy the fasta input file into another file, for making it easier
    Tfc.fasta(querys_input, "query.fa")

    # We will use it later for the counter shown in the screen while waiting
    i = how_many_sequences("query.fa")
    total = str(i)
    count = 1

    # These lists will be used for changing the accession code for the 
    # species name in the phylogenetic tree
    code_list = []
    species_list = []

    for (dirpath, dirnames, filenames) in os.walk(GBKs_directory):

        # This list contains only the genbanks in the directory
        GBK_list = only_gbks(filenames, dirpath)

        # We refill the lists of codes and species with the data of each genbank.
        for GBK in GBK_list:

            GBK = str(dirpath)+GBK

            with open(GBK, "r") as handle:
                for recordG in SeqIO.parse(handle, "genbank"):
                    code = recordG.name
                    species = recordG.annotations["organism"].replace(" ", "_")
                    break

            code_list.append(code)
            species_list.append(species)

            # We convert each genbank into a fasta format. They will be all
            # together in this file.
            Tfc.fasta_conversion(GBK, 'subject.fa')

    # We define the different values of the variables below
    identity = Blast.default(3, 30)
    coverage = Blast.default(4, 50)
    evalue = Blast.default(5, 0.00001)

    # This loop is necessary for distinguishing what query we are working with.
    # We will obtain 1 result/query of each analysis
    with open("query.fa", "r") as input_handle:
        for record in SeqIO.parse(input_handle, "fasta"):

            # We create the results directories, depending on each query
            results_query = ("Blast_Muscle_Prosite_results_"+date_time 
            + "/" +record.id)

            try:
                os.makedirs(results_query+"/Blast")
                os.mkdir(results_query+"/Muscle")
                os.mkdir(results_query+"/Trees")
                os.mkdir(results_query+"/Prosite")
                os.mkdir(results_query+"/Graphics")

            except:
                pass

            print("Sequence "+record.id +
                  " analysis started ("+str(count)+"/"+total+").")
            print(
                "-----------------------------------------------------------------")
            count += 1

            # We separate one query/iteration and we introduce it in a temporal file
            # to do BlastP.
            sequence = str(">"+record.id+"\n"+record.seq)
            tmp_file = open("file.fa", "w")
            tmp_file.write(sequence)
            tmp_file.close()
            # BlastP
            Blast.Blast("file.fa", "subject.fa", evalue, record.id,
                        "blast_results_"+record.id+".tsv")
            # We filter the results by coverage and identity and we create a fasta 
            # format file
            Blast.blast_results_to_fasta(
                'blast_results_'+record.id+".tsv", "blast_def_"+record.id+".fa", 
                identity, coverage)

            os.remove("file.fa")
            # We select from the subject.fa (initial multifasta) the hits of blast, 
            # for analize them with the prosite module
            Blast.subject_selection("subject.fa", "blast_def_"+record.id +
                                    ".fa", "blast_complete_sequences_"+
                                    record.id+".fa", record.id)

            # We add the query sequence in both files
            to_prosite = open("blast_complete_sequences_"+record.id+".fa"
            , 'a')
            to_prosite.write(str(sequence))
            to_prosite.close()

            blast_output = open("blast_def_"+record.id+".fa", 'a')
            blast_output.write(str(sequence))
            blast_output.close()

            print("\tBlast analysis finished.")

            # Alignment with Muscle
            Muscle.Muscle("blast_def_"+record.id+".fa", "alignment_"+record.id)
            try:

                # Tree with Muscle, newick format for output file
                Muscle.tree("alignment_"+record.id, "tree_"+record.id+".nw")

                tree = open("tree_"+record.id+".nw", 'r').read()

                if tree != '':
                    tree_species = open(
                        "tree_species_"+record.id+".nw", 'w')

                # We replace the accession code by the species name
                for i in range(0, len(code_list)):
                    tree = tree.replace(code_list[i], species_list[i])

                tree_species.write(tree)
                tree_species.close()

                # We move the files to their directories
                Muscle.Phylo_maketree(
                    "tree_species_"+record.id+".nw", "phylo_tree_species_"+
                    record.id)

                shutil.move("tree_"+record.id+".nw",
                            results_query+"/Trees/"+"tree_"+record.id+".nw")
                shutil.move("tree_species_"+record.id+".nw",
                            results_query+"/Trees/"+"tree_species_"+record.id+
                            ".nw")
                shutil.move("phylo_tree_species_"+record.id,
                            results_query+"/Trees/"+"phylo_tree_species_"+
                            record.id)

            except:
                print(
                    "\t\033[0;31mIt is impossible to make a tree with the "+
                    " file alignment_"+record.id+"\033[00m")

            print("\tMuscle analysis finished.")

            # Domains analysis
            PROSITE.domains_finder(
                "blast_complete_sequences_"+record.id+".fa", sys.argv[6],
                "prosite_results_"+record.id)

            print("\tPROSITE analysis finished.")
            print(
                "------------------------------------------------------------")

            try:
                # Graphing BlastP results
                BG.graphic_blast('blast_results_'+record.id+".tsv",
                                 date_time, record.id, identity, coverage)

                print("\tBlast's results graphic finished.")

            except:
                print(
                    "\t\033[0;31mIt is impossible to graph BlastP results because "+
                    "there aren't enough hits.\033[00m")

            # Graphing domain analysis results
            PG.graphic_prosite("blast_complete_sequences_" +
                               record.id+".fa", "prosite.dat", date_time, 
                               record.id)
                               
            print("\tProsite's results graphic finished.\n\n")

            # We move the files to their directories
            shutil.move('blast_results_'+record.id+".tsv",
                        results_query+"/Blast/"+'blast_results_'+record.id+
                        ".tsv")
            shutil.move("blast_def_"+record.id+".fa",
                        results_query+"/Blast/"+"blast_def_"+
                        record.id+".fa")
            shutil.move("alignment_"+record.id,
                        results_query+"/Muscle/"+"alignment_"+
                        record.id)
            shutil.move("prosite_results_"+record.id,
                        results_query+"/Prosite/"+"prosite_results_"+
                        record.id)
            shutil.move("blast_complete_sequences_"+record.id+".fa",
                        results_query+"/Blast/"+"blast_complete_sequences_"+
                        record.id+".fa")

    shutil.move("query.fa", data+"/query.fa")
    shutil.move("subject.fa", data+"/subject.fa")

    shutil.copy(sys.argv[1], data+"/"+sys.argv[1])
    shutil.copytree(sys.argv[2], data+"/"+sys.argv[2])
    shutil.copy(sys.argv[6], data+"/"+sys.argv[6])

    end(identity, coverage, evalue, date_time)

if __name__ == '__main__':
    main()
