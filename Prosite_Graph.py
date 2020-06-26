#!/usr/bin/env python3
#-*-coding: utf-8-*-

#MODULES###########################################################################

import os
import re
import sys

from Bio import SeqIO
from Bio.ExPASy import Prosite, Prodoc
from matplotlib import colors as mcolors
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import numpy as np

from PROSITE import translator

#FUNCTIONS#########################################################################


def graphic_prosite(input, database, date_time, ID):

	"""
	This function graphs the complete sequences of BlastP results, taking into
	account that its size is the same for all of them: 1000 aminoacids. Domains 
	present in the sequence are also represented by rectangles of differents
	colors. The names of the hits are in the right.  
	Arguments: input file that contains the complete sequences of the BlastP 
	results, date_time for saving the image and the query ID.  
	"""

	# We define the axes
	fig = plt.figure()
	axes = fig.subplots(nrows=1, ncols=1)

	# We obtain the number of sequences
	with open(input) as handle2:
		l=0
		for recordBlast in SeqIO.parse(handle2, "fasta"):
			l+=1

	with open(input) as handle2:
		dom_list=[]
		for recordBlast in SeqIO.parse(handle2, "fasta"):
			
			# This variable contains the different sequences
			subject_seq=str(recordBlast.seq) 
						
			with open(database, 'r') as handle:
				for recordD in Prosite.parse(handle):

					# This variable contains the different domains, translated
					a=translator(recordD.pattern)
								
					if re.search(a, subject_seq) and recordD.pattern != "":
						dom_list.append(a)

				# Drawing the proteins					
				x=np.linspace(0,1000,1000)
				y=[l]*1000
				l-=1
				axes.plot(x, y, 'red')

	# Colors list
	colors = list(sorted(mcolors.BASE_COLORS.keys()))
	colors2 = list(sorted(mcolors.TABLEAU_COLORS.keys()))
	colors3=list(sorted(mcolors.CSS4_COLORS.keys()))
	colors.extend(colors2)
	colors.extend(colors3)
	colors.remove('w')
	colors.extend(colors)

	patches=[]
	blast_names=[]

	l=0

	# We obtain the sequences number in the input file and we create a hits 
	# names list
	with open(input) as handle3:
		for recordBlast in SeqIO.parse(handle3, "fasta"):
			l+=1
			blast_names.append(recordBlast.id)
					
		lmax=l
					
	with open(input) as handle4:
		for recordBlast in SeqIO.parse(handle4, "fasta"):

			short=sorted(list(set(dom_list))) # We remove repeated domains
					
			subject_seq=str(recordBlast.seq)
						
			for i in range(0, len(short)):
						
				e= re.search(str(short[i]), subject_seq)
						
				if e:

					start_end=[] # It will contain the start and end of each domain
					start_end.append(round(e.start()*(1000/len(recordBlast.seq)),
						1))
					start_end.append(round(e.end()*(1000/len(recordBlast.seq)),
						1)) 

					x=np.linspace(float(start_end[0]), float(start_end[1]),
						int(start_end[1]-start_end[0]))
						
					y=[l]*int(start_end[1]-start_end[0])
					axes.plot(x, y, color= colors[i], linewidth=3)										
			l-=1
				
	# Printing the hits names
	for i in range(1, lmax+1):
		# If I don't distinguish between these two cases the name of the query is
		# cut
		if i == 1:
			plt.text(1110,i, blast_names[-i], fontsize=4, ha='center')
		else:
			plt.text(1110,i, blast_names[-i].replace(blast_names[-i][blast_names[
				-i].rfind("@"):], ""), fontsize=4, ha='center')

	# Graph fit		
	fig.subplots_adjust(bottom=0.3, wspace=0.33, hspace= 1)

	only_colors=[]

	# We generate a list of colors used
	j=0
	with open(input) as handle4:
		for recordBlast in SeqIO.parse(handle4, "fasta"):
			subject_seq=str(recordBlast.seq)
			for i in range(0, len(short)):
				e= re.search(short[i], subject_seq)
				if e:
					only_colors.append(colors[j])
					j+=1

	# Legend
	for i in range(len(short)):
		patches.append(mpatches.Patch(color=only_colors[i], label= short[i]))

	axes.legend(handles=patches,loc='lower center', fontsize=5, 
	     bbox_to_anchor=(0.5, -0.5),fancybox=False, shadow=False, ncol=3)

	# Axes and title setting
	axes.set_title(ID)
	axes.set_xticks(range(0, 1050, 250))
	axes.set_yticks([])
	axes.set_xlabel("Protein length", fontsize= 12)
	plt.box(False)
		
	plt.savefig("Blast_Muscle_Prosite_results_"+date_time+"/"+ID+ \
		"/Graphics/Prosite_results_graph_"+ID+".pdf")
	


