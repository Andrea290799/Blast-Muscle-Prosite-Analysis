# Blast-Muscle-Prosite-Analysis


**DESCRIPTION**

This program will be used to analize some query proteins using BlastP; to obtain a MSA and a phylogenetic trees (N-J) using MUSCLE and to find proteins' domains from Prosite database that are present in the BlasP results. Results from second and third analysis are also graphicated.

A results folder will be created. The results of the different queries are organized in their corresponding folder regarding their names. The data used is storaged in the DATA  folder, inside the results folder generated.

Results folder: Blast_Muscle_Prosite_results_[date_time].


**INSTALATION**

Muscle, Matplolib and BioPython.


**USAGE**

    ./main.py [query_file] [genbanks_directory/] [identity] [coverage] [evalue] [prosite_file]


prosite -> file that contains prosite database to analize. It has to be located in the same directory as this script.

query_file -> file that contains query sequences to analize. It has to be located in the same directory as this script.

genbanks_directory -> directory that contains genbanks files to analize. It is necessary to include a '\' after de directory's name. It has to be located in the same directory as this script.

-----------------For using the default value, introduce a non numeric character.---------------------

identity -> minimum identity percentage required. 0 <= number <= 100
coverage -> minimum coverage percentage required. 0 <= number <= 100
evalue -> minimum evalue percentage required. 0 <= number <= 100.

