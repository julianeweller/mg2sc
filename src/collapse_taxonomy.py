import anndata
import pandas as pd
import numpy as np
import scipy.sparse

def taxonomy_compare(A, B):
	'''
	A helper function to compare the level of two entries in the kraken hierarchy. Returns
	1 if A is higher in the hierarchy, 0 if A and B are equal, -1 if B is higher.
	
	Input
	-----
	A, B : ``str``
		Two strings taken from column four of hierarchy.txt, capturing the taxonomy 
		levels of the records.
	'''
	#the order of the taxonomy levels
	order = 'RDKPCOFGS'
	#get the levels of A and B in this hierarchy
	levA = order.find(A[0])
	levB = order.find(B[0])
	if levA < levB:
		#A is higher in the hierarchy
		return 1
	elif levA > levB:
		#B is higher in the hierarchy
		return -1
	else:
		#A and B are at the same level, so need more detailed investigation
		#need to turn to positions [1:], which capture the more detailed sub-hierarchy
		#(and are regular integers, if present)
		if A[1:] < B[1:]:
			#A is higher in the hierarchy
			return 1
		elif A[1:] > B[1:]:
			#B is higher in the hierarchy
			return -1
		else:
			#they're the same
			return 0

def identify_level(hierarchy, target):
	'''
	A helper function that identifies each record's corresponding entry of the appropriate
	level, for the entire kraken hierarchy file. Input as in collapse_taxonomy(). Outputs
	a dictionary with the record ID as the key and the appropriate hierarchy level as 
	the value, with entries higher in the hierarchy than the desired target level having 
	"" as the value.
	'''
	#load the file
	hierarchy = pd.read_table(hierarchy, comment='#', header=None)
	#we start off with no value
	tax_value = ''
	translation = {0:''}
	for level, tax_id, tax_name in zip(hierarchy.iloc[:,3], hierarchy.iloc[:,4], hierarchy.iloc[:,5]):
		#need to strip off all the left side space from the name
		tax_name = tax_name.lstrip()
		#check how our encountered level compares to the target level
		comp = taxonomy_compare(level, target)
		if comp == 1:
			#our record is higher in the hierarchy than the target level
			#so it can't have a value of the desired level assigned
			tax_value = ''
		elif comp == 0:
			#our record is of the exact correct hierarchy level
			#so it's the new value for anything lower in the hierarchy
			tax_value = tax_name
		#if comp == -1, then it's lower and the current tax_value should be usedd
		#one way or another, can store the value now
		translation[tax_id] = tax_value
	#and once we go across the whole thing, return the dictionary
	return translation

def collapse_taxonomy(adata, hierarchy, target):
	'''
	Collapse the single cell kraken pipeline counts based on taxonomy hierarchy, with all
	information more specific than the desired depth summed up. Returns an AnnData with the
	variable space being the summed counts at the specified point in the taxonomy hierarchy.
	
	Input
	-----
	adata : ``AnnData``
		Single cell kraken pipeline output, with the kraken taxonomy IDs in 
		``adata.var['gene_ids']``.
	hierarchy : filename
		Path to the ``hierarchy.txt`` file provided as part of the output of the kraken 
		pipeline.
	target : ``str``
		The desired hierarchy level to sum up to. Accepted values are ``["domain", "kingdom", 
		"phylum", "class", "order", "family", "genus", "species"]`` or any value of the 
		fourth column of ``hierarchy``.
	'''
	#if the user provides the target as a full name, extract and capitalise the first letter
	if len(target)>2:
		target = target[0].upper()
	#obtain the translation of the taxonomy IDs to the appropriate target level
	translation = identify_level(hierarchy, target)
	#use this translation to get the target level values for the IDs in the object
	levels = np.array([translation[i] for i in adata.var['gene_ids']])
	#prepare the output var names and collapse the count matrix accordingly
	#(remove '' as it's of no interest)
	var_names = np.unique(levels)
	var_names = var_names[var_names != '']
	#this goes better on CSC
	holder = adata.X.tocsc()
	c = scipy.sparse.csr_matrix(np.hstack([np.sum(holder[:,levels==n],1) for n in var_names]))
	#create output object
	bdata = anndata.AnnData(c, obs=adata.obs)
	bdata.var_names = var_names
	return bdata