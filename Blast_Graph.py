#!/usr/bin/env python3
# -*-coding: utf-8-*-

#MODULES###########################################################################

import os

from Bio import SeqIO
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import numpy as np

#FUNCTIONS#########################################################################
	
def graphic_blast(input_file, date_time, ID, ident, cove):

	"""
	This function graphs the BlastP results. Legth represents coverage, and colors,
	identity:
		Red: identity < 50
		Yellow: 50 < identity < 90
		Green: identity > 90
	
	Hits names are in the left. 

	Arguments: input file that contains BlastP hits, date_time varibale for saving
	the image, protein ID, identity, coverage. 
	"""

	# We define the axes
	fig = plt.figure()
	axes = fig.subplots(nrows=1, ncols=1)

	# We iniciate the necessary inicial (_i) lists
	blast_names_i=[]
	identity_i= []
	coverage_i = []

	# We iniciate the definitive lists for drawing
	identity=[] 
	blast_names=[]

	# We filter the input file by intensity and coverage and we refill the inicial
	# lists
	file = open(input_file, 'r')

	for line in file.readlines():
		fields = line.rstrip().split("\t")
		name = fields[0]
		id = fields[1]
		cov  = fields[2]

		if id >= ident and cov >= cove:

			identity_i.append(id)
			coverage_i.append(cov)
			blast_names_i.append(name)

			# We order the list so it's sorted by coverage	
			coverage=sorted(coverage_i)

	# We order the other lists regarding the coverage list order
	for i in range(len(coverage)):
		for j in range(len(coverage)):
			if coverage_i[j]==coverage[i]:
				identity.append(identity_i[j])
				blast_names.append(blast_names_i[j])
				coverage_i[j]=""

	l = len(coverage)

	# We draw. Space sponsored by Art Academy. 	
	for i in range(0, l):

		y = [i]*int(coverage[i])

		# As we said before, the x axis represent the coverage		
		x =np.linspace(0,float(coverage[i]), int(coverage[i]))

		# Hits names
		plt.text(-10,i, blast_names[i].replace(blast_names[i][blast_names[i
		].rfind("@"):], ""), fontsize=3, ha='center')
				
		if float(identity[i]) >= 90:
			axes.plot(x, y, 'green') 
				
		elif float(identity[i]) < 50:
			axes.plot(x, y, 'red')
		else:
			axes.plot(x, y, 'orange') 
				
	# Graph fit			
	fig.subplots_adjust(bottom=0.3, wspace=0.33)

	# We set the legend
	red_patch = mpatches.Patch(color='red', label='Identity < 50')
	green_patch = mpatches.Patch(color='green', label='Identity >= 90')
	orange_patch = mpatches.Patch(color='orange', label='50 < Identity < 90')

	axes.legend(handles=[green_patch,orange_patch,red_patch],loc='upper center', 
				bbox_to_anchor=(0.5, -0.2),fancybox=False, shadow=False, ncol=3)

	# Axes and title setting
	axes.set_title(ID)
	axes.set_xticks(range(0, 150, 50))
	axes.set_yticks([])
	axes.set_xlabel("Coverage", fontsize= 12)
	plt.box(False)


	plt.savefig("Blast_Muscle_Prosite_results_"+date_time+"/"+ID+ \
		"/Graphics/Blast_results_graph_"+ID+".pdf")
	



		
