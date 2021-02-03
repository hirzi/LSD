#!/usr/bin/python

# BAM2AlleleFreqs version
version = "BAM2AlleleFreqs version: 19_03_2019"

############################################## This script calculates allele frequencies from a list of mpileup files ##############################################
## The input is a list of mpileup files produced from BAM files (>samtools mpileup example.bam) and the output is a dict of allele frequencies for input into the sumstatscalc.py script ##
## BAM file from which mpileup file is produced must be SORTED ##

# Import modules
import pandas
import re
import numpy
#import os

# Set working directory
#os.chdir("/Users/luqman/Documents/Hirzi/ETHPHD/Analysis/PoolSimulator")

# Example usage:
# Calculate allele counts for single mpileup file (whether individual or poolseq).
#allele_freqs_single = calc_allele_freqs_mpileup(filename = "LOW_GTGAAA_realn_recal_sorted_scaffold4_size532381.txt" , start_position = 2500, end_position = 12500, quality_filter = 20, region = "single")
# Calculate allele counts for a single population (takes as input multiple individuals).
#file_list_e = "HIGH_ACAGTG_realn_recal_sorted_scaffold4_size532381.txt,HIGH_CAGATC_realn_recal_sorted_scaffold4_size532381.txt,HIGH_CTTGTA_realn_recal_sorted_scaffold4_size532381.txt,HIGH_ACAGTG_realn_recal_sorted_scaffold4_size532381.txt,HIGH_CAGATC_realn_recal_sorted_scaffold4_size532381.txt,HIGH_CTTGTA_realn_recal_sorted_scaffold4_size532381.txt" # first line of "test_individuals_pop.filelist"
#allele_freqs_df = calc_pop_allele_freqs_mpileup(file_list_e, no_inds = 6, start_position = 2500, end_position = 12500, quality_filter = 20, region = "single")	
# Calculate allele frequencies for/across multiple (individual data) populations.
#allele_freqs_df = calc_pooled_allele_freqs_mpileup("test_6individuals_pop.filelist", no_pops = 2, pop_size = [6,6], pooled = False, start_position = 2500, end_position = 12500, quality_filter = 20, min_allele_count = 1, min_cov = 10, max_cov = 200, error_method = "4", region = "single")	
# Calculate allele frequencies for/across multiple (poolseq) populations.
#allele_freqs_df = calc_pooled_allele_freqs_mpileup("scaffold4_size532381_filelist.txt", no_pops = 6, pooled = True, start_position = 2500, end_position = 12500, quality_filter = 20, min_allele_count = 1, min_cov = 10, max_cov = 200, error_method = "4", region = "single")	
# Final wrapper
#raw_allele_freqs = mpileup2lsd(allele_freqs_df = allele_freqs_df, no_pops = no_pops, region = region, window_size = window_size, start_position = start_position, end_position = end_position)
#raw_allele_freqs = mpileup2lsd(allele_freqs_df = allele_freqs_df, no_pops = 2, region = "single", window_size = 0, start_position = 2500, end_position = 12500)

#%%

############################################## NOTES: SAMPLE REQUIREMENTS ##############################################

# Note 0: This script returns the quality-filtered minor allele counts
# Note 1: If error_method set to "0" or "1", minimum allele count filter does not apply, even if set (since these two error methods don't consider minimum allele count in the lsd_high script). Minimum allele count filter only applies when using error methods "2" and "3".
# Note 2: The quality filter removes reads with Phred quality scores of less than or equal to the defined tag. E.g. -q 20 removes all reads with quality scores less than and equal to 20.

#%%

############################################## NOTES: SAMPLE REQUIREMENTS ##############################################

#1) Tajima's D requires there to be >= 4 individuals per population.
#2) Variant class frequencies (from the SFS, e.g. singletons, doubletons, tripletons etc) require n individuals (polarised data) or 2n individuals (polarised data) per population for calculation, where n is the variant class frequency desired. By default, only doubletons and tripletons are outputted (singletons are by default excluded because they are too sensitive to error). Adjust accordingly depending on needs. 
#3) Negative FST and DXY values are collapsed to 0.

#%%

### First we define the functions that modify and filter the input mpileup file.
## For details on the mpileup format, refer to: https://en.wikipedia.org/wiki/Pileup_format
# To account for indels in the reads column
def remove_indels(string):
	string_list = list(string)
	# We find the positions of the indels (these are denoted by +/-xACTGn, where + = insertion, - = deletion, x = size of indel, ACTGn = the bases of the indel))
	indel_positions_iter = re.finditer(r"[\-|\+]{1}\d+", string)
	indel_positions = [m.start() for m in indel_positions_iter]
	# We find the size of the the indels
	indel_size = re.findall(r"[\-|\+]{1}\d+",string)
	indel_size_int = [int(re.search(r'\d+', x).group()) for x in indel_size]
	# We make a list of the positions in the reads string to remove (i.e. of the indel)
	indel_indexes_raw = []
	for i in range(len(indel_size)):	
		indexes = [x for x in range(indel_positions[i],indel_positions[i]+indel_size_int[i]+len(indel_size[i]))]
		indel_indexes_raw.append(indexes)
	indel_indexes = [item for sublist in indel_indexes_raw for item in sublist]
	# We delete these characters from the string, making sure to delete from highest index to lowest
	for index in sorted(indel_indexes, reverse=True):
	    del string_list[index]
	string_no_indels = ''.join(string_list)
	return string_no_indels
# To account for caret and the succeeding single ASCII quality characters
def remove_carets(string):
	string_list = list(string)
	caret_positions_iter = re.finditer(r"\^", string)
	caret_positions = [m.start() for m in caret_positions_iter]
	# To account for cases where the caret indicates a preceding caret's quality score and NOT a start read position, we must include the following condition (i.e. in cases where ...^^...):
	carets_as_qualityscore = []		
	for count, i in enumerate(caret_positions):
		for j in caret_positions[count+1:count+2]:
	#		print(str(i)+"_"+str(j))
			if j-i == 1:
				carets_as_qualityscore.append(j)
	caret_positions = [x for x in caret_positions if x not in carets_as_qualityscore]
	caret_positions.extend([x+1 for x in caret_positions])
	start_read_chars = list(set(caret_positions))
	for index in sorted(start_read_chars, reverse=True):
		del string_list[index]
	string_no_carets = ''.join(string_list)
	return string_no_carets
# To account for n and * characters. This function removes the corresponding quality scores and the calling of this function must precede the calling of the remove_n_asterisks function below
def remove_n_asterisks_quality(string,quality):
	quality_list = list(quality)
	asterisk_positions_iter = re.finditer(r"[*]", string)
	n_positions_iter = re.finditer(r"N", string)
	asterisk_positions = [m.start() for m in asterisk_positions_iter]
	n_positions = [m.start() for m in n_positions_iter]
	undefined_positions = asterisk_positions
	undefined_positions.extend(n_positions)
	for index in sorted(undefined_positions, reverse=True):
		del quality_list[index]	
	quality_no_undefined = ''.join(quality_list)
	return quality_no_undefined
# To account for n and * characters. This function removes these characters from the reads column
def remove_n_asterisks(string):
	string_list = list(string)
	asterisk_positions_iter = re.finditer(r"[*]", string)
	n_positions_iter = re.finditer(r"N", string)
	asterisk_positions = [m.start() for m in asterisk_positions_iter]
	n_positions = [m.start() for m in n_positions_iter]
	undefined_positions = asterisk_positions
	undefined_positions.extend(n_positions)
	for index in sorted(undefined_positions, reverse=True):
		del string_list[index]	
	string_no_undefined = ''.join(string_list)
	return string_no_undefined
# To check that the length of quality scores matches that of the read bases
def len_test(x):
	if len(x["quality_mod"]) != len(list(x["reads_mod"])):
		output = "mismatch"
	else:
		output = "match"
	return output
# To select the reads that satisfy the defined quality condition.
def reads_quality_filter(reads_string, quality_string, quality):
#	filtered_reads = list(numpy.array(reads_string)[(quality < numpy.array(quality_string))])
	filtered_reads = [read for read, quality_score in zip(reads_string, quality_string) if quality_score > quality]
	return filtered_reads
# To select minor allele count, where minor_allele argument defines the (n+1)th most frequent allele, i.e. count_minor_allele(ACTG_sorted_count,minor_allele = 1) will return the count for the 2nd most frequent allele.
def count_minor_allele(ACTG_sorted_count,minor_allele):
	minor_count = ACTG_sorted_count[minor_allele]
	return minor_count
# To find the minor allele index (in the meta-population), where -> 0:A, 1:C, 2:T, 3:G. If reference allele is given, we take the most frequent of the 3 alternate alleles. If no reference allele is given (i.e. ref = N), we select the meta-population minor allele.
def index_minor_allele(ref_allele, ACTG_count):
	ref_str = str(ref_allele)
	if ref_str == "A" or ref_str == "C" or ref_str == "T" or ref_str == "G":
		if ref_str == "A":
			ref_allele_index = 0
		elif ref_str == "C":
			ref_allele_index = 1
		elif ref_str == "T":
			ref_allele_index = 2
		elif ref_str== "G":
			ref_allele_index = 3
		ACTG_count_mod = [x for x in ACTG_count]
		del ACTG_count_mod[ref_allele_index]
		sorted_ACTG_count = sorted(list(ACTG_count_mod), reverse = True)
		minor_allele_count = sorted_ACTG_count[0]
		minor_allele_index = ACTG_count.index(minor_allele_count)
	else:
		sorted_ACTG_count = sorted(list(ACTG_count), reverse = True)
		minor_allele_count = sorted_ACTG_count[1]
		minor_allele_index = ACTG_count.index(minor_allele_count)
	return minor_allele_index
# To calculate and output (population-specific) allele frequencies. Here the minor index is (0:A, 1:C, 2:T, 3:G) respectively
def allele_freqs(read_counts, minor_index, minor_allele):
	minor_count = read_counts[minor_index]
	filtered_depth = sum(read_counts)
	if minor_count > 0: # i.e. if the site is a polymorphic site (when minor_allele refers to the first minor allele)
		allele_freq = float(minor_count / filtered_depth)
	else:
		allele_freq = 0
	return allele_freq
# To apply a minimum allele count filter (population-specific). Note: a minimum quality filter only applies to error methods 2,3 and 4.
def min_allele_count_filter(read_counts, min_allele_count, error_method):
	sorted_read_counts = sorted(list(read_counts), reverse = True)
	if error_method == "2":
		if (sorted_read_counts[2] > min_allele_count) or (sorted_read_counts[3] > min_allele_count) or (0 < sorted_read_counts[1] <= min_allele_count):		
			output = "FAIL"
		else:
			output = "PASS"
	elif error_method == "3" or error_method == "4":
		if (sorted_read_counts[1] == sorted_read_counts[2] and sorted_read_counts[1] > 0) or (0 < sorted_read_counts[1] <= min_allele_count):
			output = "FAIL"
		else:
			output = "PASS"
	return output	

### This function, which depends on the pre-defined functions above, takes a SINGLE mpileup file and returns the allele counts in a pandas dataframe format. The start and end positions are only mandatory for the 'single' option.
def calc_allele_freqs_mpileup(filename, quality_filter, region, start_position = 0, end_position = 0):
	raw_full_df = pandas.read_csv(filename, sep='\t', header=None, keep_default_na=False)	# keep_default_na=False is important as read characters could be of the form "NA"
	raw_full_df.columns = ['scaffold', 'position', 'reference', 'depth', 'reads', 'quality']
	if region == "single":
		# Specify region in scaffold or chromosome, according to start and end position
		raw_df_temp = raw_full_df.loc[(raw_full_df['position'] >= start_position) & (raw_full_df['position'] <= end_position)]
		raw_df = raw_df_temp.reset_index(drop=True)
		del(raw_df_temp)
	elif region == "multiple":
		raw_df = raw_full_df
	del(raw_full_df)
	# Define a unique position number, allowing for multiple scaffolds
	# 'position' now refers to a unique position identifier, by concatenating the scaffold name with the scaffold position number. 'position original' refers to the original (non-unique) scaffold position number.
	raw_df.rename(columns = {'position':'position_original'}, inplace = True)
	raw_df["position"] = raw_df["scaffold"].astype(str) + "_pos" + raw_df["position_original"].astype(str)
	cols = ['scaffold', 'position_original', 'position', 'reference', 'depth', 'reads', 'quality']	
	raw_df = raw_df[cols]
	# Make sure datatype in pandas dataframe is in the appropriate datatype (specifically we want the reads and quality columns to be in string format)
	raw_df["reads"] = raw_df["reads"].astype(str)
	raw_df["quality"] = raw_df["quality"].astype(str)	
	# Make everything upper case
	raw_df["reads_mod"] = raw_df["reads"].apply(lambda x: x.upper())	
	# Account for indels
	raw_df["reads_mod"] = raw_df.apply(lambda x: remove_indels(x["reads_mod"]),axis=1)	
	# Account for caret and following ASCII quality characters
	raw_df["reads_mod"] = raw_df.apply(lambda x: remove_carets(x["reads_mod"]),axis=1)	
	# Remove special characters
	raw_df["reads_mod"] = raw_df["reads_mod"].apply(lambda x: re.sub(r'[^ACTGN*]', '', x))	
	# Account for n and * characters. This function removes the corresponding quality scores and the calling of this function must precede the calling of the remove_n_asterisks function below
	raw_df["quality"] = raw_df.apply(lambda x: remove_n_asterisks_quality(x["reads_mod"],x["quality"]),axis=1)
	# Account for n and * characters. This function removes these characters from the reads column
	raw_df["reads_mod"] = raw_df.apply(lambda x: remove_n_asterisks(x["reads_mod"]),axis=1)
	# Convert quality scores
	raw_df["quality_mod"] = raw_df["quality"].apply(lambda x: [(ord(y)-33) for y in list(x)])
	# Convert reads_mod column from string to list
	raw_df["reads_mod"] = raw_df["reads_mod"].apply(lambda x: list(x))
	# Check that length of quality scores matches that of read bases
	raw_df["check"] = raw_df.apply(len_test,axis=1)
	errors = raw_df.loc[raw_df['check'] == "mismatch"]
	if len(errors) != 0:
		print("Mismatch error!")
		error_rows = raw_df.loc[raw_df['check'] == "mismatch"]
		print(error_rows)
	# We select the reads that satisfy the defined quality condition.
	# To avoid a potential bug when inputting list data into a new pandas column. See: https://github.com/pandas-dev/pandas/pull/18577. Also consider the tags 'raw' and 'reduce'
	try:
		raw_df["filtered_reads"] = raw_df.apply(lambda x: reads_quality_filter(x["reads_mod"],x["quality_mod"],quality_filter),axis=1, reduce = True)	
	except ValueError:
		raw_df["filtered_reads"] = ''
		raw_df["filtered_reads"] = raw_df.apply(lambda x: reads_quality_filter(x["reads_mod"],x["quality_mod"],quality_filter),axis=1)	
	# Count the (filtered) reads
	countA = list(raw_df["filtered_reads"].apply(lambda x: x.count("A")))
	countC = list(raw_df["filtered_reads"].apply(lambda x: x.count("C")))
	countT = list(raw_df["filtered_reads"].apply(lambda x: x.count("T")))
	countG = list(raw_df["filtered_reads"].apply(lambda x: x.count("G")))
	countACTG = (list(zip(countA,countC,countT,countG)))
	raw_df["read_counts"] = pandas.Series(countACTG, index=raw_df.index)	
	# Remove columns for streamlined output file
	final_df = raw_df.drop("reads_mod", axis=1)
	final_df.drop("quality_mod", axis=1, inplace=True)
	final_df.drop("check", axis=1, inplace=True)
	final_df.drop("filtered_reads", axis=1, inplace=True)
	return final_df		

### For individual data. This function takes a LIST of individual mpileup files (that belong to the same population), and returns the population allele counts in a pandas dataframe format. It expands the above (single file) function to a many files (populations) function.
# Note, even though the input here are individuals (not pooled individuals), the corresponding entries in the pandas dataframe are referred to as "pops" rather than "inds", for compatibility and legacy reasons.
def calc_pop_allele_freqs_mpileup(file_list_e, no_inds, quality_filter, region, start_position = 0, end_position = 0):
	# Print program version number
	print(version)
	# Input data. For individual data, the file list is a row of populations, with each row (population) comprising of n comma-separated individuals 
	filelist = file_list_e.split(",")
	# Define empty list and dicts	
	positions_master=[]
	master={}
	master_trimmed={}
	# Read as pandas dataframe the mpileup files
	for count,file in enumerate(filelist):
		print("Reading file: " + str(file))
		pop_df = calc_allele_freqs_mpileup(file, quality_filter, region, start_position, end_position)
		master["Pop"+str(count+1)] = pop_df
		positions = set(list(pop_df["position"]))
		positions_master.append(positions)		
	# Calculate the common (intersection of) sites between the 6 populations
	across_pop_sites = list(set.intersection(*positions_master))	# This finds the intersection of a list of sets, where *a_list is list expansion
	# Trim the mpileup files retaining only these common sites
	for pop in range(no_inds):
		pop_df = master["Pop"+str(pop+1)]
		pop_df_trimmed_raw = pop_df[pop_df["position"].isin(across_pop_sites)]
		pop_df_trimmed = pop_df_trimmed_raw.reset_index(drop=True)
		master_trimmed["Pop"+str(pop+1)] = pop_df_trimmed	
	# Concatenate the 6 populations' mpileup files into one file
	popall_df = master_trimmed["Pop1"] 	
	popall_df.rename(columns={'read_counts': 'read_counts1'}, inplace=True)	
	for i in range(2,no_inds+1):
		popall_df.loc[:,"read_counts"+str(i)] = master_trimmed["Pop"+str(i)]["read_counts"]
	# Add a column of unique position numbers (integers)
	# 'position_reordered' gives a unique integer position number to each site. We need this for later region selection. We can't use 'position' because these are not unique, and cannot use 'position_original' because these are strings, not integers (see line 175).	
	popall_df["position_reordered"] = numpy.asarray(list(range(1,len(popall_df)+1)))
	# We split the read counts to A,C,T,G counts (columns), to facilitate downstream calculations
	for i in range(1,no_inds+1):
		popall_df[["Pop"+str(i)+"A","Pop"+str(i)+"C","Pop"+str(i)+"T","Pop"+str(i)+"G"]] = pandas.DataFrame(popall_df["read_counts"+str(i)].values.tolist(), index= popall_df.index)
	# We calculate the meta-population ACTG counts, to find the minor allele (in case reference is not given)
	popall_df["AllPops_A"] = popall_df["Pop1A"]
	popall_df["AllPops_C"] = popall_df["Pop1C"]
	popall_df["AllPops_T"] = popall_df["Pop1T"]
	popall_df["AllPops_G"] = popall_df["Pop1G"]
	for base in ["A","C","T","G"]:
		for i in range(2,no_inds+1):
			popall_df["AllPops_" + str(base)] += popall_df["Pop"+str(i)+str(base)]	
	popall_df["AllPops_ACTG"] = popall_df.apply(lambda x: list([x["AllPops_A"],x["AllPops_C"],x["AllPops_T"],x["AllPops_G"]]),axis=1) 
	# Remove the separate ACTG count and read count columns.
	for base in ["A","C","T","G"]:
		popall_df.drop("AllPops_" + str(base), axis=1, inplace=True)
		for i in range(1,no_inds+1):                                 			   
			popall_df.drop("Pop"+str(i)+str(base), axis=1, inplace=True)
	for pop in range(no_inds):	
		popall_df.drop("read_counts"+str(pop+1), axis=1, inplace=True)
	# Remove and rename columns for streamlined and compatible final output file
	final_df = popall_df.drop("position_reordered", axis=1)
	final_df.rename(columns={'AllPops_ACTG': 'read_counts'}, inplace=True)			
	return final_df

### This function takes a LIST (pooled data) or LIST OF LIST (individual population data) of mpileup files and returns the allele frequencies in a pandas dataframe format. It expands the above (single/population file) function to many files (populations) function, and calculates the corresponding allele frequencies.
def calc_pooled_allele_freqs_mpileup(file_list, no_pops, pooled, quality_filter, min_allele_count, error_method, region, pop_size = None, min_cov = 0, max_cov = 999, start_position = 0, end_position = 0):
	# Print program version number
	print(version)
	# Input data
	with open(file_list) as f:
		filelist = f.read().splitlines()
	# Define empty list and dicts	
	positions_master=[]
	master={}
	master_trimmed={}
	# Read as pandas dataframe the mpileup files
	for count,file in enumerate(filelist):
		print("Reading file: " + str(file))
		if pooled == True:
			pop_df = calc_allele_freqs_mpileup(file, quality_filter, region, start_position, end_position)
		elif pooled == False:
			no_inds = pop_size[count]
			pop_df = calc_pop_allele_freqs_mpileup(file, no_inds, quality_filter, region, start_position, end_position)
		master["Pop"+str(count+1)] = pop_df
		positions = set(list(pop_df["position"]))
		positions_master.append(positions)
	# Calculate the common (intersection of) sites between the 6 populations
	across_pop_sites = list(set.intersection(*positions_master))	# This finds the intersection of a list of sets, where *a_list is list expansion
	# Trim the mpileup files retaining only these common sites
	for pop in range(no_pops):
		pop_df = master["Pop"+str(pop+1)]
		pop_df_trimmed_raw = pop_df[pop_df["position"].isin(across_pop_sites)]
		pop_df_trimmed = pop_df_trimmed_raw.reset_index(drop=True)
		master_trimmed["Pop"+str(pop+1)] = pop_df_trimmed	
	# Concatenate the 6 populations' mpileup files into one file
	popall_df = master_trimmed["Pop1"]
	popall_df.rename(columns={'depth': 'depth1'}, inplace=True)
	popall_df.rename(columns={'read_counts': 'read_counts1'}, inplace=True)
	for i in range(2,no_pops+1):
		popall_df.loc[:,"depth"+str(i)] = master_trimmed["Pop"+str(i)]["depth"]
		popall_df.loc[:,"read_counts"+str(i)] = master_trimmed["Pop"+str(i)]["read_counts"]
	# Add a column of unique position numbers (integers)
	# 'position_reordered' gives a unique integer position number to each site. We need this for later region selection. We can't use 'position' because these are not unique, and cannot use 'position_original' because these are strings, not integers (see line 175).	
	popall_df["position_reordered"] = numpy.asarray(list(range(1,len(popall_df)+1)))
	# We split the read counts to A,C,T,G counts (columns), to facilitate downstream calculations
	for i in range(1,no_pops+1):
		popall_df[["Pop"+str(i)+"A","Pop"+str(i)+"C","Pop"+str(i)+"T","Pop"+str(i)+"G"]] = pandas.DataFrame(popall_df["read_counts"+str(i)].values.tolist(), index= popall_df.index)
	# We calculate the meta-population ACTG counts, to find the minor allele (in case reference is not given)
	popall_df["AllPops_A"] = popall_df["Pop1A"]
	popall_df["AllPops_C"] = popall_df["Pop1C"]
	popall_df["AllPops_T"] = popall_df["Pop1T"]
	popall_df["AllPops_G"] = popall_df["Pop1G"]
	for base in ["A","C","T","G"]:
		for i in range(2,no_pops+1):
			popall_df["AllPops_" + str(base)] += popall_df["Pop"+str(i)+str(base)]	
	popall_df["AllPops_ACTG"] = popall_df.apply(lambda x: list([x["AllPops_A"],x["AllPops_C"],x["AllPops_T"],x["AllPops_G"]]),axis=1) 
	# Remove the separate ACTG count columns
	for base in ["A","C","T","G"]:
		popall_df.drop("AllPops_" + str(base), axis=1, inplace=True)
		for i in range(1,no_pops+1):                                 			   
			popall_df.drop("Pop"+str(i)+str(base), axis=1, inplace=True)      
	# Count minor allele across all populations (to remove invariant sites and for the minor_index later)         
	popall_df["AllPops_ACTG_sort"] = popall_df["AllPops_ACTG"].apply(lambda x: sorted(x,reverse=True))
	popall_df["minor_count"] = popall_df["AllPops_ACTG_sort"].map(lambda x: count_minor_allele(x,1))
	# Remove invariant sites
	popall_segsites_df_temp = popall_df.loc[(popall_df["minor_count"] > 0)]
	popall_segsites_df = popall_segsites_df_temp.reset_index(drop=True)
	# To find the (meta-population) minor allele index, where -> 0:A, 1:C, 2:T, 3:G. This minor allele is used calculation of allele frequencies. If reference allele is given, we take the most frequent of the 3 alternate alleles. If no reference allele is given (i.e. ref = N), we select the meta-population minor allele.
	popall_segsites_df["minor_index"] = popall_segsites_df.apply(lambda x: index_minor_allele(x["reference"], x["AllPops_ACTG"]),axis=1)
	# Apply minimum allele count filter (this is the population-specific minimum allele count; NOT the meta-population minor allele count). Only applied if/when using error methods 2, 3 or 4.
	if error_method == "0" or error_method == "1":
		popall_segsites_filtered_df = popall_segsites_df.copy()
	elif error_method == "2" or error_method == "3" or error_method == "4":
		for pop in range(no_pops):		
			popall_segsites_df["min_allele_filter"+str(pop+1)] = popall_segsites_df.apply(lambda x: min_allele_count_filter(x["read_counts"+str(pop+1)],min_allele_count, error_method), axis=1)
		popall_segsites_filtered_df_temp = popall_segsites_df.copy()
		for pop in range(no_pops):
			popall_segsites_filtered_df_temp = popall_segsites_filtered_df_temp.loc[popall_segsites_filtered_df_temp["min_allele_filter"+str(pop+1)] == "PASS"]
		popall_segsites_filtered_df = popall_segsites_filtered_df_temp.reset_index(drop=True)	
	# Apply min and max depth filter (applies per sample depth)
	popall_segsites_and_depth_filtered_df_temp = popall_segsites_filtered_df	
	for i in range(1,no_pops+1):
		popall_segsites_and_depth_filtered_df_temp = popall_segsites_and_depth_filtered_df_temp.loc[(popall_segsites_and_depth_filtered_df_temp["depth"+str(i)] >= min_cov) & (popall_segsites_and_depth_filtered_df_temp["depth"+str(i)] <= max_cov)]
	popall_segsites_and_depth_filtered_df = popall_segsites_and_depth_filtered_df_temp.reset_index(drop=True)	
	# Calculating and outputting allele frequencies (of meta-population minor allele; NOT population specific minor allele). USE THIS DF FOR TROUBLESHOOTING PURPOSES
	for pop in range(no_pops):
		popall_segsites_and_depth_filtered_df["allele_freq"+str(pop+1)] = popall_segsites_and_depth_filtered_df.apply(lambda x: allele_freqs(x["read_counts"+str(pop+1)], x["minor_index"], 1),axis=1)
		popall_segsites_and_depth_filtered_df["allele_freq"+str(pop+1)] = popall_segsites_and_depth_filtered_df["allele_freq"+str(pop+1)].astype(float)			
	# Remove columns for streamlined final output file
	final_df = popall_segsites_and_depth_filtered_df.drop("quality", axis=1)
	for pop in range(no_pops):	
		final_df.drop("read_counts"+str(pop+1), axis=1, inplace=True)
		final_df.drop("depth"+str(pop+1), axis=1, inplace=True)
		if error_method == "2" or error_method == "3" or error_method == "4":
			final_df.drop("min_allele_filter"+str(pop+1), axis=1, inplace=True)
	final_df.drop("minor_count", axis=1, inplace=True)
	final_df.drop("reads", axis=1, inplace=True)
	final_df.drop("AllPops_ACTG_sort", axis=1, inplace=True)
	# Reorder columns
	# cols = final_df.columns.tolist()
	cols = ['scaffold', 'position_original', 'position', 'position_reordered', 'reference', 'AllPops_ACTG', 'minor_index']
	for pop in range(no_pops):
		cols.append("allele_freq"+str(pop+1))
	final_df = final_df[cols]	
	return final_df
	### CONSIDER: An alternative approach to filtering directly in the df here is to output a list of sites to be discarded, which can be handled and removed downstream by sumstatscalc.py, as you have done so for lsd_high.py.  

### This is the wrapper (main) function that performs everything; taking the output of the above function as the input and returning allele frequencies as the output in the appropriate format for the sumstatscalc.py script
def mpileup2lsd(allele_freqs_df, no_pops, region, window_size, start_position = 0, end_position = 0):
	# We define open dictionaries, into which we'll write our iteration (N/A here) and allele_freq/index_discarded_sites results. The order/organisation follows from that of the lsd_high.py simulator
	results = {}
	allelefreqs_master_array = {}
	index_discarded_sites_master_array = {}
	allele_freqs = {}
	# Specify region in scaffold or chromosome, according to start and end position
	if window_size > 0:
		if region == "single":
			region_allele_freqs_df = allele_freqs_df.loc[(allele_freqs_df['position_original'] >= start_position) & (allele_freqs_df['position_original'] <= end_position)]	
		elif region == "multiple":
			region_allele_freqs_df = allele_freqs_df.loc[(allele_freqs_df['position_reordered'] >= start_position) & (allele_freqs_df['position_reordered'] <= end_position)]
	elif window_size == 0:	
		region_allele_freqs_df = allele_freqs_df
	# Convert allele frequencies dataframes into list format
	for pop in range(no_pops):
		allele_freqs["Population"+str(pop+1)] = list(region_allele_freqs_df["allele_freq"+str(pop+1)])
	# Append data to dicts
	allelefreqs_master_array["Iteration_1"] = allele_freqs
	index_discarded_sites_master_array["Iteration_1"] = []
	results["Allele_freqs"] = allelefreqs_master_array
	results["index_discarded_sites"] = index_discarded_sites_master_array

	# Return final results dict
	return results
	
#%%
