##########################################################
# Python Script to compute Protein Metabolite relation +
# Uses HMDB database (locally stored as XML files)
# Requires HMDB id of metabolites
# Uses parsing of XML files
# Author: Paurush Praveen
# Contact: praveen@cosbi.eu
##########################################################

#!/usr/bin/python

import glob
import csv
import sys
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import xml.etree.ElementTree as ET

#from itertools import groupby
#from xml.dom.minidom import parse


#yxml = []
#yxml.append(glob.glob("*.xml"))
# Read Input files and IDs
colnames = ['Compound', 'ID']
data = pd.read_csv('/home/praveen/Paurush/COSBI/Alzherimers/Metabolic_data/AD_metabolites.csv', delimiter="\t", names=colnames)
IDs=list(data.ID)

# Function to extract different information
def extractInfo(info,IDs, outfile):
	if info=='pathway':	
		x=tablePathway(IDs)
		writeFile(x, outfile)
	elif info=='disease':
		x=tableDisease(IDs)
		writeFile(x, outfile)
	elif info=='proteins':
		x=tableProteins(IDs)
		writeFile(x, outfile)
	else:
		print 'The information in not available'
		print 'The valid available information is'
		print '<pathway>'
		print '<disease>'
		print '<proteins>'
		print 'Please enter one of these values as <info> parameter (1)'
		exit
		
###############################
# Pathway annotation
def tablePathway(IDs):
	pathlist = []
	idlist = []
	for id in IDs:
		if id != 'XXX':
			print(id)
			xmlfile=id+'.xml'
			root=ET.parse(xmlfile).getroot()
			for node in root.findall('pathways'):
				pathlist.append(node.text)
				idlist.append(id)
	res=[idlist,pathlist]
	return(res)	


###############################
# Protein annotation
def tableProteins(IDs):
	genelist = []
	idlist = []
	uniprotlist = []
	for id in IDs:
		if id != 'XXX':
			print(id)
			xmlfile=id+'.xml'
			root=ET.parse(xmlfile).getroot()
			for node in root.findall('protein_associations'):
				for proteins in node.findall('protein'):
					for genes in proteins.findall('gene_name'):
						genelist.append(genes.text)
						idlist.append(id)
					for uprots in proteins.findall('uniprot_id'):
						uniprotlist.append(uprots.text)
	genelist=['None' if v is None else v for v in genelist]
	res=[idlist,genelist,uniprotlist]
	return(res)	
###########################
# Disease annotation
def tableDisease(IDs):
	dislist = []
	idlist = []
	for id in IDs:
		if id != 'XXX':
			print(id)
			xmlfile=id+'.xml'
			root=ET.parse(xmlfile).getroot()
			for node in root.findall('diseases'):
				dislist.append(node.text)
				idlist.append(id)
	res=[idlist,dislist]
	return(res)	

#######################
def writeFile(d, outfile):
	a=zip(*d)
	csv_file = open(outfile, "wb")
	writer = csv.writer(csv_file)
	writer.writerows(a)

######################
#x=tableProteins(IDs)
#writeFile(x, 'outfile.csv')
#y=tableDisease(IDs)
#writeFile(y, 'outfiledis.csv')
#z=tablePathway(IDs)
#writeFile(z, 'outfilepath.csv')




