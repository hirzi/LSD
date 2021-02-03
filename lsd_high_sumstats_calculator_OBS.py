#!/usr/bin/python

# lsd_high_sumstats_calculator_OBS version
version = "lsd_high_sumstats_calculator_OBS version: 30_03_2020"

############################################## ALLELE FREQUENCIES CONVERTER (OBSERVED) AND SUMSTATS CALCULATOR ##############################################

### NOTES AND EXAMPLES ###

# This is the main (master) function that runs the BAM2AlleleFreqs (observed data) converter and summary statistics calculator functions. Use this to run analyses for the observed data #
# It takes as input a list of mpileup files (produced from BAM files via: samtools mpileup example.bam) and produces an output of summary statistics #
# BAM file from which mpileup file is produced must be sorted (samtools sort -T /tmp/aln.sorted -o aln.sorted.bam aln.bam) #
# EXTENSION: This version has been extended to produce sumstat files for multiple windows in a scaffold.

## EXAMPLE USAGE
# cd /Users/luqman/Documents/Hirzi/ETHPHD/Analysis/PoolSimulator
# Example usage (GENOME-WIDE, NEUTRAL REGIONS, 5000kbMakerFlanked-1000kbOtherFlanked, anchored scaffolds > 50kb): python lsd_high_sumstats_calculator_OBS.py genomeWide_5000_1000_anchored50_filelist.txt -d 40 -d 40 -d 40 -d 40 -d 40 -d 40 -q 0 -m 2 -o GW_Neutral_5000_1000_anchored_50kb -f ABC- -r multiple --mindepth 10 --maxdepth 200 --windowSize 5000 --scaffold --pooled
# Example usage (ALL (FULLY) NEUTRAL SCAFFOLDS, CHR 1): python lsd_high_sumstats_calculator_OBS.py Chr1_filelist.txt -d 40 -d 40 -d 40 -d 40 -d 40 -d 40 -q 0 -m 2 -o Chr1 -f ABC- -r multiple --windowSize 5000 --mindepth 10 --maxdepth 100 --scaffold --pooled
# Example usage (ALL (FULLY) NEUTRAL SCAFFOLDS, CHR 13): python lsd_high_sumstats_calculator_OBS.py Chr13_filelist.txt -d 40 -d 40 -d 40 -d 40 -d 40 -d 40 -q 0 -m 2 -o Chr13 -f ABC- -r multiple --windowSize 5000 --mindepth 10 --maxdepth 100 --scaffold --pooled
# Example usage using command line input: python lsd_high_sumstats_calculator_OBS.py scaffold4_size532381_filelist.txt -d 40 -d 40 -d 40 -d 40 -d 40 -d 40 -q 20 -m 1 -o scaffold4_size532381_minDP10maxDP200 -f ABC --startPos 1 --endPos 532381 --mindepth 10 --maxdepth 200 --windowSize 10000 -r single --pooled
# Example usage using command line input: python lsd_high_sumstats_calculator_OBS.py scaffold1_size1318325_1_614794_filelist.txt -d 40 -d 40 -d 40 -d 40 -d 40 -d 40 -q 20 -m 1 -o scaffold1_size1318325_1_614794_minDP10maxDP200 -f ABC -r single --mindepth 10 --maxdepth 200 --pooled
# Example usage, individual (non-pooled) data, with 6 individuals per population: python lsd_high_sumstats_calculator_OBS.py test_6individuals_pop.filelist -d 6 -d 6 -q 20 -m 1 -o scaffold4_size532381_test_individuals_pop -f ABC -r single
# Example usage, individual (non-pooled) data, with 6 individuals per population (10kb windows): python lsd_high_sumstats_calculator_OBS.py test_6individuals_pop.filelist -d 6 -d 6 -q 20 -m 1 -o scaffold4_size532381_test_individuals_pop_windows -f ABC --startPos 1 --endPos 532381 --windowSize 10000 -r single
# Script won't work with < 6 individuals per pop (unpolarised data) as is because it by default outputs variant frequencies up till the tripleton class. Either change outputted frequency classes, or ignore if you anyway won't work with such small sample sizes.
# Script won't work with < 4 individuals per pop Tajima's D require at least this number of individuals per population.

## NOTES
# Note 0: A window size of 0 (default) will calculate the allele frequencies and summary statistics for the region defined between the start_position and the end_position. A window size > 0 will calculate allele frequencies and summary statistics on a per window basis.
# Note 1: Here, the minimum allele count is defined as follows: If allele frequency is greater than minimum allele count, then we keep. If equal to or below minimum allele count, we discard. E.g. a minimum allele count = 2 means any alleles with counts of 2 or below are discarded (need a minimum of 3 to keep).
# Note 2: The quality filter removes reads with Phred quality scores of less than or equal to the defined tag. E.g. -q 20 removes all reads with quality scores less than and equal to 20.
# Note 3: The error rate in this script defines the rate at which a individual base read is read incorrectly during the random sampling process of the pooling process (random sampling n coverage times with replacement), per site. It is NOT exactly analogous/comparable to genotyping or SNP call error rate (it is much lower).
# Note 4: This error simulation operates in addition to the sampling bias that will result from stochasticity from the random sampling.
# Note 5: No genotyping or SNP calling is performed in this script. Everything is kept in term of allele frequencies and probabilies.
# Note 6: msms, ms, cosi2, etc (i.e. most coalescent simulators) work with haploid samples. To use with diploids, simply take 2N individuals.
# Note 7: At present, we are simulating without recombination. We make the assumption that within simulated loci, there is no recombination; hoever between loci, there is free recombination. This is a simplified approach and should work well assuming you use relatively short loci.

### Note X: If you recieve the error: "ValueError: could not broadcast input array from shape (i) into shape (j)", where i and j are numbers, this points to lines 163-167 in the BAM2AlleleFreqs.py script. There appears to be a bug in pandas, when outputting data of a certain shape into a new column. See: https://github.com/pandas-dev/pandas/pull/18577

## NOTES: SAMPLE REQUIREMENTS
# Note 1: Tajima's D requires there to be >= 4 individuals per population.
# Note 2: Variant class frequencies (from the SFS, e.g. singletons, doubletons, tripletons etc) require n individuals (polarised data) or 2n individuals (polarised data) per population for calculation, where n is the variant class frequency desired. By default, only doubletons and tripletons are outputted (singletons are by default excluded because they are too sensitive to error). Adjust accordingly depending on needs. 
# Note 3: Negative FST and DXY values are collapsed to 0.
# Note 4: In the case of two populations, private S will always be the same. This is because a site that is private in pop1 must necessarily be private in pop2.

##  We import the necessary modules
import argparse
import os
import pandas
import BAM2AlleleFreqs
import sumstatscalc
#import csv
#import glob
#import numpy
from datetime import datetime

# os.chdir("/Users/luqman/Documents/Hirzi/ETHPHD/Analysis/PoolSimulator")

#%%

### PARSE ARGUMENTS ###

## We parse arguments from the command line using argparse 
parser = argparse.ArgumentParser(description="LSD formatter and summary statistics calculator for observed sequence data (version 0.1, Luqman 2020)", add_help=True)

# We allow for two input options:
"""
# 1. We parse arguments from a supplied parameter file:
if glob.glob("/Users/luqman/Documents/Hirzi/ETHPHD/Analysis/PoolSimulator/*.params"):     # Glob will return a list of files matching the wildcard-accessible path. If there are no files, it will return an empty list. This is really just os.listdir and fnmatch.filter together.
	parser.add_argument("paramfile", action="store", help="Parameters may either be entered in the command line as arguments or provided in a separate parameter file. If a .params file is provided in the working directory, this will be used. If no parameter file provided in the working directory, retype -h for full argument options (full command-line argument options are only shown when no .params file exists in the working path)")
	parser.add_argument("-o", action="store", default="output", dest="output", help="prefix for output files")
	args = parser.parse_args()	
	paramfile = args.paramfile	
	with open(paramfile) as paramfile:     # We parse arguments from from the input parameter file 
		reader = list(csv.reader(paramfile, delimiter="\t"))
		arg_dict = {}
		for row in reader:
			arg_dict[row[0]] = row[1]
	
	##  We associate the parsed arguments with our global variables (defined in the functions below)
	filename = arg_dict["filename"]
#	no_iterations = int(arg_dict["no_iterations"])
	pop_size = [int(x) for x in arg_dict["pop_size"].replace(" ","").split(",")]
	span =  arg_dict["across"]	
	plotSFS = bool(arg_dict["plotSFS"] == "false")
	plotPi = bool(arg_dict["plotPi"] == "false")
	output_prefix = arg_dict["output"]
	output_format = arg_dict["output_format"]
	quality_filter = int(arg_dict["quality_filter"])
	min_allele_count = int(arg_dict["min_allele_count"])
	start_position = int(arg_dict["start_position"])
	end_position = int(arg_dict["end_position"])
"""

# 2. Or we parse arguments directly from the command line:	
#else:
parser.add_argument("filename", action="store", help="The input file, which is a list of mpileup files.")
#parser.add_argument("-n", action="store", dest="nrep", default=1, type=int, help="The number of independent iterations (alternatively, can be interpreted as the number of loci simulated). N/A for observed data")
parser.add_argument("-d", action="append", required=True, dest="demesizes", help="The sample sizes of the different populations.")
parser.add_argument("-a", action="store", dest="across", default="segsites", help="Define whether to calculate Fst summary statistics across only segregating sites (segsites) or across all sites (allsites).")
parser.add_argument("--polarised", action="store_true", default=False, dest="polarised", help="Define whether to consider the simulated sequences as polarised (ancestral:0 - derived:1) or unpolarised (major:most frequent allele - minor: 2nd most frequent allele).")
parser.add_argument("--pooled", action="store_true", default=False, dest="pooled", help="Define whether data (BAM, mpileup files) is individual data or pooled data.")
parser.add_argument("--error_method", action="store", dest="error_method", default="4", help="Error calculation method (five options: 0, 1, 2, 3 and 4). Error method 0 doesn't account for error, only for stochastic sampling bias from the pooling process. Error methods 1, 2, 3 and 4 account for stochastic sampling bias from the pooling process in addition to randomly including error at each read per site at (probability) error_rate. Error method 1: If the error leads to a change to alleles 0 or 1, the read with the error is retained and appended to the total counts. If the error leads to inobservable alleles 2 and 3, the reads are discarded (no sites discarded in this method). Error method 2: All reads are retained, regardless of whether the error leads to alleles 0, 1, 2 or 3. If allele 2 and/or 3 is above the minimum allele count, we consider it a triallelic site and discard the site. If they are equal to or below the minimum allele count, and allele 1 is above the minimum allele count, we consider the site biallelic (for alleles 0 and 1) , and retain the original coverage (no reads are discarded in this method). Error method 3: All reads are retained, regardless of whether the error leads to alleles 0, 1, 2 or 3. This method counts the two most frequent alleles (the ancestral + the most frequent minor allele between allele 1, allele 2 and allele 3), given that the 1) most frequent minor allele is > the minumum allele count, and 2) the the most frequent minor allele is NOT equal to the second most frequent minor allele (otherwise site is considered triallelic and discarded). If these two points are satisfied, the two least frequent alleles are ignored and we consider the site biallelic, whilst retaining the original coverage (no reads are discarded here). Error method 4: An approximation of error method 3, which results in approximately 5x shorter runtimes. The approximation assumes that all retained erroneous sites have (min_allele_count + 1) number of errors. Sites with less number of errors are filtered by the minimum allele criterion, sites with more errors are exponentially much less likely, and to the closest approximation, results in negligable number of entries.")
parser.add_argument("--plotSFS", action="store_true", default=False, dest="plotSFS", help="To output plots of site frequency spectrum (one for each population).")
parser.add_argument("--plotPi", action="store_true", default=False, dest="plotPi", help="To output plots of nucleotide diversity, pi, across site positions (one for each population).")
parser.add_argument("-o", action="store", default="summary_stats_temp.txt", dest="output", help="Prefix for output files.")
parser.add_argument("-f", action="store", default="ABC", dest="output_format", help="Format of output file. 'Full' outputs a convenient, human-readable summary output while 'ABC' outputs the file in the format required for ABCToolbox.")
parser.add_argument("-q", action="store", dest="quality_filter", default=0, type=int, help="Remove reads with PHRED quality-score less than this value.")
parser.add_argument("-m", action="store", dest="min_allele_count", default=0, type=int, help="Apply a minimum allele count filter.")
parser.add_argument("-r", action="store", dest="region", help="Select 'single' when calculating summary statistics for a single region or scaffold, or 'multiple' when calculating summary statistics for multiple scaffolds/regions. The 'single' option allows defining custom sub-regions within a scaffold (e.g. between sites 100 to 90000); the 'multiple' option analyses all sites within all given regions/scaffolds.")
parser.add_argument("--mindepth", action="store", dest="min_cov", default=0, type=int, help="Minimum depth filter")
parser.add_argument("--maxdepth", action="store", dest="max_cov", default=999, type=int, help="Maximum depth filter")
parser.add_argument("--startPos", action="store", dest="start_position", default=0, type=int, help="The start position (bp) of the region of interest. Only called (is mandatory) when program is run for a single region; ignored if running for multiple regions, as in this case, all sites of all regions are analysed.")
parser.add_argument("--endPos", action="store", dest="end_position", default=0, type=int, help="The end position (bp) of the region of interest. Only called (is mandatory) when program is run for a single region; ignored if running for multiple regions, as in this case, all sites of all regions are analysed.")
parser.add_argument("--windowSize", action="store", dest="window_size", default=0, type=int, help="Splits the scaffold into windows. Summary statistics will be calculated per window.")
parser.add_argument("--windowStep", action="store", dest="window_step", default=None, type=int, help="Defines the length (number of sites) till the subsequent window. If window_step < window_size, windows will be overlapping; if window_step = window_size, windows will be non-overlapping and adjacent; if window_step > window_size, windows will be non-overlapping with gaps equal to the difference.")
#parser.add_argument("--windowMissing", action="store", dest="missing_data_filter", default=0.5, type=float, help="(Scaffold) windows with a certain fraction of missing sites will be excluded from calculation of summary statistics.")
parser.add_argument("--baseJump", action="store", dest="base_jump", default=100, type=int, help="When sites are missing in the scaffold (for defining windows), this parameter defines the jump size to the next available site.")
parser.add_argument("--scaffold", action="store_true", default=False, dest="scaffold_specific", help="If set as true, calculates allele frequencies and summary statistics per scaffold (scaffolds not merged).")

args = parser.parse_args()

##  We associate the parsed arguments with our global variables (defined in the functions below)
filename = args.filename
#no_iterations = args.nrep
pop_size = [int(x) for x in args.demesizes]
span = args.across
polarised = args.polarised
pooled = args.pooled
error_method = args.error_method
plotSFS = args.plotSFS
plotPi = args.plotPi
output_prefix = args.output
output_format = args.output_format
quality_filter = args.quality_filter
min_allele_count = args.min_allele_count
region = args.region
min_cov = args.min_cov
max_cov = args.max_cov
start_position = args.start_position
end_position = args.end_position
window_size = args.window_size
window_step = args.window_step
#missing_data_filter = args.missing_data_filter
base_jump = args.base_jump
scaffold_specific = args.scaffold_specific

# We define the indirect arguments:
scaffold_seq_length = end_position - start_position
no_pops = len(pop_size)

##  To test if this works:
print("Filename is: " + str(filename))
#print("Number of iterations is: " +str(no_iterations))
print("List of respective population sample sizes are: " + str(pop_size))
print("Summary statistics will be calculated for " + str(region) + " region(s)/scaffold(s)")
if scaffold_specific == True:
	print("Summary statistics will be calculated per scaffold (scaffolds will not be merged)")
if region == "single" and start_position != 0 and end_position != 0:
	print("Summary statistics will be calculated for the region between " + str(start_position) + "bp and " + str(end_position) + "bp")
	print("The length of the scaffold sequence is: " + str(scaffold_seq_length) + "bp")
if region == "single" and start_position == 0 and end_position == 0:	
	print("Summary statistics will be calculated for the whole scaffold" )
if window_size > 0:
	print("Window size is " + str(window_size) + "bp")
	print("Summary statistics are calculated per " + str(window_size) + "bp windows")
	print("Windows are generated in " + str(window_step) + "bp steps")
print("Remove reads of quality < " + str(quality_filter))
print("Error method: " + str(error_method))
print("The observed data is polarised: " + str(polarised))
print("The observed data is pooled: " + str(pooled))
print("Applying a filter on minimum allele count. Alleles with counts </= " + str(min_allele_count) + " will be removed")
print("Minimum depth filter: " + str(min_cov))
print("Maximum depth filter: " + str(max_cov))
#print("Fst and Dxy statistics are calculated across " + str(span))
print("Plots of SFS are produced: " + str(plotSFS))
print("Plots of pi are produced: " + str(plotPi))
print("Prefix for output files is: " + str(output_prefix))
print("Output format is: " + str(output_format))

#%%

### DEFINE THE FUNCTIONS ###

# This function returns a set (list) of the common sites, among all populations
def shared_sites(file_list, region, no_pops, pooled, scaffold_specific, list_of_scaffolds=[]):
	with open(file_list) as f:
		file_list_e = f.read().splitlines()
	if pooled == True:
		filelist = file_list_e
	elif pooled == False:
		filelist_temp = [x.split(",") for x in file_list_e]
		filelist = [item for sublist in filelist_temp for item in sublist]
	# Define empty list and dicts	
	positions_master=[]
	master={}	
	# Read as pandas dataframe the mpileup files
	for count,file in enumerate(filelist):
#		print("Reading file: " + str(file))
		test_df = pandas.read_csv(file, sep='\t', header=None)
		test_df.columns = ['scaffold', 'position_original', 'reference', 'depth', 'reads', 'quality']
		#pop_df = calc_allele_freqs_mpileup(file, start_position, end_position, quality_filter)
		if region == "multiple":
			test_df["position"] = test_df["scaffold"].astype(str) + "_pos" + test_df["position_original"].astype(str)
			cols = ['scaffold', 'position_original', 'position', 'reference', 'depth', 'reads', 'quality']	
			test_df = test_df[cols]
		elif region == "single":
			test_df.rename(columns = {'position_original':'position'}, inplace = True)
		master["Pop"+str(count+1)] = test_df
		positions = set(list(test_df["position"]))
		positions_master.append(positions)		
	# Calculate the common (intersection of) sites between the 6 populations
	across_pop_sites = list(set.intersection(*positions_master))	# This finds the intersection of a list of sets, where *a_list is list expansion
	across_pop_sites.sort()
	if region == "single":
		print("Number of common base pairs (across all populations) in scaffold is: " + str(len(across_pop_sites)))
		return across_pop_sites
	if region == "multiple" and scaffold_specific == False:
		across_pop_sites = list(range(1,len(across_pop_sites)+1))  # This line follows that we have added and defined the unique position index: "position_reordered" = numpy.asarray(list(range(1,len(popall_df)+1))) in the dataframe, following the same intersection step. Hence we can match according to this index list.
		print("Number of common base pairs (across all populations) in scaffold is: " + str(len(across_pop_sites)))
		return across_pop_sites
	elif region == "multiple" and scaffold_specific == True:
		# Produce a dataframe of scaffold, position_original and position by trimming the mpileup files to retain only the common (intersected) sites
		master_trimmed={}
		for pop in range(no_pops):
			pop_df = master["Pop"+str(pop+1)]
			pop_df_trimmed_raw = pop_df[pop_df["position"].isin(across_pop_sites)]
			pop_df_trimmed = pop_df_trimmed_raw.reset_index(drop=True)
			master_trimmed["Pop"+str(pop+1)] = pop_df_trimmed	
		# Remember, here we just need the intersecting positions (which is identical and included in all population dataframes), so we can just choose 1 population
		common_sites_df = master_trimmed["Pop1"] 	
		# Drop unnecessary columns
		common_sites_df.drop("reads", axis=1, inplace=True)
		common_sites_df.drop("quality", axis=1, inplace=True)
		common_sites_df.drop("depth", axis=1, inplace=True)
		common_sites_df.drop("reference", axis=1, inplace=True)
	#	common_sites_df_copy = common_sites_df.copy()
		# Index by scaffold, to allow easy partioning of dataframe by scaffold later
		common_sites_df.set_index(keys=['scaffold'], drop=False,inplace=True)
#		names = common_sites_df['scaffold'].unique().tolist() # to avoid bugs (e.g. in the case that there is a mismatch in the scaffold list of across_pop_sites_per_scaffold and allele_freqs_df_scaffold, which may arise in the case that some scaffolds have no variant sites), we loop and extract only over scaffolds which have variant sites (taken from the allele_freqs_df_scaffold scaffold list)
		# Partition dataframe by scaffold
		common_sites_df_list_scaffold_partitions = []
		for scaffold in list_of_scaffolds: # this names_scaffols is taken from names_scaffolds = allele_freqs_df['scaffold'].unique().tolist()
	#		print(scaffold)
			common_sites_df_list_scaffold_partitions.append(common_sites_df.loc[common_sites_df.scaffold==str(scaffold)])
		# Produce a LIST of across pop sites per scaffold	
		across_pop_sites_per_scaffold_list = []
		for i in common_sites_df_list_scaffold_partitions:
			across_pop_sites_per_scaffold = list(i["position_original"])
			across_pop_sites_per_scaffold_list.append(across_pop_sites_per_scaffold)
		return across_pop_sites_per_scaffold_list

# NOTE, THIS FUNCTION HAS BEEN DEPRECIATED AND REPLACED BY THE FUNCTION BELOW. This function returns the start and end positions of the desired windows, defining FIXED regions (e.g. if 10k windows: 1-10000, 10001-20000, 20001 - 30000, ...etc. ), given a defined start position, end position, and window size. Tolerance for amount of missing data can be set *NOTE if this tolerance is <1, different summary statistics will be affected differently (e.g. S, watteson's theta correlate directly with length of sequence, while e.g. pi, Fst, etc. are divided over sequence length)
# Note 0: This function hasn't been updated to work with multiple regions/scaffolds.
# Note 1: positions are just that positions, NOT python indexes, and thus they start from 1.
# Note 2: the start_position and end_position arguments are the start and end positions/sites of the e.g. scaffold.
#def make_windows_fixedRegions(start_position, end_position, window_size, missing_data_filter, base_jump):
#	positions={}
#	# Make list of start and end positions, based on window size and region of interest
#	list_start_position = []
#	list_end_position = []
#	while start_position + window_size <= end_position: # or <
#		list_start_position.append(start_position)
#		start_position = start_position + window_size
#		list_end_position.append(start_position - 1)
#	# Account for the cases where start position is missing (by searching for the next available starting position)
#	list_start_position_data_raw =[]
#	for start_position in list_start_position:
#		if start_position in across_pop_sites:
#			list_start_position_data_raw.append(start_position)
#		else:
#			while start_position not in across_pop_sites:
#				start_position = start_position + base_jump
#			list_start_position_data_raw.append(start_position)
#	# Account for the cases where end position is missing (by searching for the previous available ending position)
#	list_end_position_data_raw =[]
#	for end_position in list_end_position:
#		if end_position in across_pop_sites:
#			list_end_position_data_raw.append(end_position)
#		else:
#			while end_position not in across_pop_sites:
#				end_position = end_position - base_jump
#			list_end_position_data_raw.append(end_position)
#	# Remove windows with no or little (defined by missing_data_filter tag) data.
#	size_end_start_list = []
#	for i,j in zip(list_start_position_data_raw,list_end_position_data_raw):
#		size_end_start_list.append(len(across_pop_sites[i:j]))	
#	size_end_start = numpy.array(size_end_start_list)
#	gaps_to_remove = numpy.where(size_end_start < (window_size * missing_data_filter))[0]
#	# Remove start and end positions which define the windows (gaps) in scaffold with insufficient data
#	list_start_position_data = []
#	for index, start_pos in enumerate(list_start_position_data_raw):
#		if index not in gaps_to_remove:
#			list_start_position_data.append(start_pos)
#	list_end_position_data = []
#	for index, end_pos in enumerate(list_end_position_data_raw):
#		if index not in gaps_to_remove:
#			list_end_position_data.append(end_pos)
#	# Check output and return
#	if len(list_start_position_data) != len(list_end_position_data):
#		print("Error: Different number of start and end positions")
#	positions["start_positions"] = list_start_position_data
#	positions["end_positions"] = list_end_position_data
#	return positions

# This function returns the start and end positions of the desired windows, defining windows of FIXED length (e.g. if 10k windows: the 1st 10000 sites, the 2nd 10000 sites, the 3rd 10000 sites, etc), given a defined start position, end position, and window size.
# Note 1: positions are just that positions, NOT python indexes, and thus they start from 1.
# Note 2: the start_index and end_index arguments are the start and end indexes of the total number of shared sites.
# Note 3: window_step is defined as the length (number of sites) till the subsequent window. If window_step < window_size, windows will be overlapping, if window_step = window_size, windows will be non-overlapping and adjacent, if window_step > window_size, windows will be non-overlapping with gaps equa to the difference.		
def make_windows_fixedLength(across_pop_sites, region, start_position, end_position, start_index, end_index, window_size, window_step = None):
	positions = {}
	list_start_indexes = []
	list_end_indexes = []
	while start_index + window_size <= end_index:
		list_start_indexes.append(start_index)
		if window_step != None:
			final_index = start_index + window_size
			list_end_indexes.append(final_index)
			start_index = start_index + window_step
		elif window_step == None:
			start_index = start_index + window_size
			list_end_indexes.append(start_index)		
	list_start_positions_all = []
	for index in list_start_indexes:
		start_pos = across_pop_sites[index]
		list_start_positions_all.append(start_pos)
	list_end_positions_all = []
	for index in list_end_indexes:
		end_pos = across_pop_sites[index]
		list_end_positions_all.append(end_pos)
	if region == "single":
		list_start_positions_trim_start = [x for x in list_start_positions_all if x >= start_position]	
		num_trimmed_start = len(list_start_positions_all) - len(list_start_positions_trim_start)
		list_end_positions_trim_start = list_end_positions_all[num_trimmed_start:]
		list_end_positions_trim_start_end = [x for x in list_end_positions_trim_start if x <= end_position]
		num_trimmed_end = len(list_end_positions_trim_start) - len(list_end_positions_trim_start_end)
		if num_trimmed_end == 0:
			list_start_positions_trim_start_end = list_start_positions_trim_start
		elif num_trimmed_end > 0:
			list_start_positions_trim_start_end = list_start_positions_trim_start[:(-num_trimmed_end)]	
		positions["start_positions"] = list_start_positions_trim_start_end
		positions["end_positions"] = list_end_positions_trim_start_end
	elif region == "multiple":
		positions["start_positions"] = list_start_positions_all
		positions["end_positions"] = list_end_positions_all	
	return positions

#%%

### HERE IS THE MAIN FUNCTION ###

##  Finally, we run the BAM2AlleleFreqs and sumstats_calc functions:
if __name__ == "__main__":

	startTime = datetime.now()
	
	print(version)

	if window_size == 0:
		
		# Obtain the set (list) of the common sites, among all populations
		across_pop_sites_all = shared_sites(file_list = filename, region = region, no_pops = no_pops, pooled = pooled, scaffold_specific = False)
		
		# Calculates sequence length from which to calculate summary statistics
		if start_position != 0 and end_position != 0:
			start_pos_index = across_pop_sites_all.index(start_position)
			end_pos_index = across_pop_sites_all.index(end_position)
			seq_length = len(across_pop_sites_all[start_pos_index:end_pos_index])
		elif start_position == 0 and end_position == 0:
			start_position = across_pop_sites_all[0]
			end_position = across_pop_sites_all[-1]
			seq_length = len(across_pop_sites_all)
			
		# Generate allele frequency dataframe for selected sites, from a list of mpileup files
		allele_freqs_df = BAM2AlleleFreqs.calc_pooled_allele_freqs_mpileup(file_list = filename, no_pops = no_pops, pooled = pooled, quality_filter = quality_filter, min_allele_count = min_allele_count, error_method = error_method, region = region, pop_size = pop_size, min_cov = min_cov, max_cov = max_cov, start_position = start_position, end_position = end_position)
		
		if scaffold_specific == False:			
			if region == "single":
				print("Calculating allele frequencies and summary statistics for scaffold region: " + str(start_position) + "bp - " + str(end_position) + "bp")
			elif region == "multiple":
				print("Calculating allele frequencies and summary statistics for multiple (merged) scaffolds, of sum total length of: " + str(seq_length) + "bp")
			
			# Convert allele frequencies dataframe into input format for LSD-High (lsd_high)
			raw_allele_freqs = BAM2AlleleFreqs.mpileup2lsd(allele_freqs_df = allele_freqs_df, no_pops = no_pops, region = region, window_size = window_size, start_position = start_position, end_position = end_position)		# note here, as we want, start and end position arguments are not called, since window == 0 always in this clause. Region selection already performed in above step.
		
			# Calculate the summary statistics
			sumstatscalc.sumstats_calc (raw_allele_freqs, pop_size = pop_size, seq_length = seq_length, span = span, polarised = polarised, output = output_format, plotSFS = plotSFS, plotPi = plotPi)
		
			# Output summary statistic files
			if region == "single":
				os.rename("summary_stats_temp.txt", output_prefix + "_" + str(start_position) + "_" + str(end_position) + "_M" + str(min_allele_count) + "_Q" + str(quality_filter) + ".obs")
			elif region == "multiple":
				os.rename("summary_stats_temp.txt", output_prefix + "_allRegions" + "_M" + str(min_allele_count) + "_Q" + str(quality_filter) + ".obs")
				
		elif scaffold_specific == True:		
			# First we split dataframe into multiple scaffold-specific dataframes
			allele_freqs_df.set_index(keys=['scaffold'], drop=False,inplace=True)
			names_scaffolds = allele_freqs_df['scaffold'].unique().tolist()
			# We can then loop over the scaffold-specific dataframes
			allele_freqs_df_per_scaffold = []
			for scaffold in names_scaffolds:
			#	print(scaffold)
				allele_freqs_df_per_scaffold.append(allele_freqs_df.loc[allele_freqs_df.scaffold==str(scaffold)])
	
			# Obtain the set (list) of the common sites, among all populations, per scaffold
			across_pop_sites_per_scaffold = shared_sites(file_list = filename, region = region, no_pops = no_pops, pooled = pooled, scaffold_specific = True, list_of_scaffolds=names_scaffolds)	
	
			# Then we loop (the sumstats calculation, in windows) over all scaffolds
			output_names = []
			for i in range(len(across_pop_sites_per_scaffold)):
				# Define start and end indexes and positions per scaffold, as well as the sequence length from which to calculate summary statistics.
				start_index = 0
				end_index = len(across_pop_sites_per_scaffold[i]) - 1
				start_position = across_pop_sites_per_scaffold[i][0]
				end_position = across_pop_sites_per_scaffold[i][-1]
				seq_length = len(across_pop_sites_per_scaffold[i])

				allele_freqs_df_scaffold = allele_freqs_df_per_scaffold[i]
	
				print("Calculating allele frequencies and summary statistics for scaffold: " + str(names_scaffolds[i]) + "...")
						
				# Select window region from allele frequency dataframe and convert this into input format for LSD-High (lsd_high). Define region = "single" since we're calculating on a PER scaffold basis.
				raw_allele_freqs = BAM2AlleleFreqs.mpileup2lsd(allele_freqs_df = allele_freqs_df_scaffold, no_pops = no_pops, region = "single", window_size = window_size, start_position = start_position, end_position = end_position)
			
				# Calculate the summary statistics
				sumstatscalc.sumstats_calc (raw_allele_freqs, pop_size = pop_size, seq_length = seq_length, span = span, polarised = polarised, output = output_format, plotSFS = plotSFS, plotPi = plotPi)

				# Output summary statistic files
				output_name = output_prefix + "_" + names_scaffolds[i] + ".obs"
				output_names.append(output_name)
				os.rename("summary_stats_temp.txt", output_name)	
			
		# Should you wish to output the sumstats of the multiple scaffolds to the same file (instead of one output file per region), you can include the code below. Default setting as is.
		if region == "multiple":
			with open (output_prefix + "_start" + str(across_pop_sites_all[0]) + "_end" + str(across_pop_sites_all[-1]) + "_M" + str(min_allele_count) + "_Q" + str(quality_filter) + ".obs", "w") as master_file:
				for file in output_names[:1]:    
					f = open(file, "r")
					master_file.write(f.read())
					master_file.writelines("\n")
					f.close
					os.remove(file)
				for file in output_names[1:]:    
					f = open(file, "r")
					sumstat_line = f.read().splitlines()[1] + "\n"
					master_file.writelines(sumstat_line)
					f.close
					os.remove(file)
			if scaffold_specific == True:
				with open(output_prefix + "_start" + str(across_pop_sites_all[0]) + "_end" + str(across_pop_sites_all[-1]) + "_M" + str(min_allele_count) + "_Q" + str(quality_filter) + ".scaffoldlist", "w") as names_file:
					for scaffold_name in output_names:
						names_file.write("%s\n" % scaffold_name)
			
	elif window_size > 0:
		
		# Obtain the set (list) of the common sites, among all populations
		across_pop_sites_all = shared_sites(file_list = filename, region = region, no_pops = no_pops, pooled = pooled, scaffold_specific = False)
	
		# Define overall start and end positions
		if start_position == 0 and end_position == 0:
			start_position = across_pop_sites_all[0]
			end_position = across_pop_sites_all[-1]
	
		# Generate total allele frequency dataframe, from a list of mpileup files
		allele_freqs_df = BAM2AlleleFreqs.calc_pooled_allele_freqs_mpileup(file_list = filename, no_pops = no_pops, pooled = pooled, quality_filter = quality_filter, min_allele_count = min_allele_count, error_method = error_method, region = region, pop_size = pop_size, min_cov = min_cov, max_cov = max_cov, start_position = start_position, end_position = end_position)
	
		if scaffold_specific == False:
			start_index = 0
			end_index = len(across_pop_sites_all) - 1
	
			# Obtain the start and end positions of the desired windows, given a defined start position, end position, and window size. Accounts for missing bases in the scaffold (data).
			positions = make_windows_fixedLength(across_pop_sites = across_pop_sites_all, region = region, start_position = start_position, end_position = end_position, start_index = start_index, end_index = end_index, window_size = window_size, window_step = window_step)	
			start_positions = positions["start_positions"]
			end_positions = positions["end_positions"]
	
			# Calculate summary statistics across windows, window-by-window
			output_names = []
			for start_pos, end_pos in zip(start_positions,end_positions):
	
				print("Calculating allele frequencies and summary statistics for all scaffolds (merged); window " + str(start_pos) + "-" + str(end_pos) + "...")
						
				# Select window region from allele frequency dataframe and convert this into input format for LSD-High (lsd_high)
				raw_allele_freqs = BAM2AlleleFreqs.mpileup2lsd(allele_freqs_df = allele_freqs_df, no_pops = no_pops, region = region, window_size = window_size, start_position = start_pos, end_position = end_pos)
	
				# Calculates sequence length from which to calculate summary statistics. This is the window size minus the missing bases
				start_pos_index = across_pop_sites_all.index(start_pos)
				end_pos_index = across_pop_sites_all.index(end_pos)
				seq_length = len(across_pop_sites_all[start_pos_index:end_pos_index])
				
				# Calculate the summary statistics
				sumstatscalc.sumstats_calc (raw_allele_freqs, pop_size = pop_size, seq_length = seq_length, span = span, polarised = polarised, output = output_format, plotSFS = plotSFS, plotPi = plotPi)
			
				# Output summary statistic files
				output_name = output_prefix + "_window_" + str(start_pos) + "_" + str(end_pos) + ".obs"
				output_names.append(output_name)
				os.rename("summary_stats_temp.txt", output_name)
	
		elif scaffold_specific == True:
			# First we split dataframe into multiple scaffold-specific dataframes
			allele_freqs_df.set_index(keys=['scaffold'], drop=False,inplace=True)
			names_scaffolds = allele_freqs_df['scaffold'].unique().tolist()
			# We can then loop over the scaffold-specific dataframes
			allele_freqs_df_per_scaffold = []
			for scaffold in names_scaffolds:
			#	print(scaffold)
				allele_freqs_df_per_scaffold.append(allele_freqs_df.loc[allele_freqs_df.scaffold==str(scaffold)])
	
			# Obtain the set (list) of the common sites, among all populations, per scaffold
			across_pop_sites_per_scaffold = shared_sites(file_list = filename, region = region, no_pops = no_pops, pooled = pooled, scaffold_specific = True, list_of_scaffolds=names_scaffolds)	
	
			# Then we loop (the sumstats calculation, in windows) over all scaffolds
			output_names = []
			for i in range(len(across_pop_sites_per_scaffold)):
				# Define start and end indexes and positions per scaffold
				start_index = 0
				end_index = len(across_pop_sites_per_scaffold[i]) - 1
				start_position = across_pop_sites_per_scaffold[i][0]
				end_position = across_pop_sites_per_scaffold[i][-1]
				# Obtain the start and end positions of the desired windows, given a defined start position, end position, and window size, PER SCAFFOLD. Accounts for missing bases in the scaffold (data).
				positions = make_windows_fixedLength(across_pop_sites = across_pop_sites_per_scaffold[i], region = region, start_position = start_position, end_position = end_position, start_index = start_index, end_index = end_index, window_size = window_size, window_step = window_step)	
				start_positions = positions["start_positions"]
				end_positions = positions["end_positions"]
				allele_freqs_df_scaffold = allele_freqs_df_per_scaffold[i]
	
				# Calculate summary statistics across windows, window-by-window			
				for start_pos, end_pos in zip(start_positions,end_positions):  # need to have scaffold specific start and end positions!
	
					print("Calculating allele frequencies and summary statistics for scaffold: " + str(names_scaffolds[i]) + "; window " + str(start_pos) + "-" + str(end_pos) + "...")
							
					# Select window region from allele frequency dataframe and convert this into input format for LSD-High (lsd_high). Define region = "single" since we're calculating on a PER scaffold basis.
					raw_allele_freqs = BAM2AlleleFreqs.mpileup2lsd(allele_freqs_df = allele_freqs_df_scaffold, no_pops = no_pops, region = "single", window_size = window_size, start_position = start_pos, end_position = end_pos)
	
					# Calculates sequence length from which to calculate summary statistics. This is the window size minus the missing bases
					start_pos_index = across_pop_sites_per_scaffold[i].index(start_pos)
					end_pos_index = across_pop_sites_per_scaffold[i].index(end_pos)
					seq_length = len(across_pop_sites_per_scaffold[i][start_pos_index:end_pos_index])
					
					# Calculate the summary statistics
					sumstatscalc.sumstats_calc (raw_allele_freqs, pop_size = pop_size, seq_length = seq_length, span = span, polarised = polarised, output = output_format, plotSFS = plotSFS, plotPi = plotPi)
	
					# Output summary statistic files
					output_name = output_prefix + "_" + names_scaffolds[i] + "_window_" + str(start_pos) + "_" + str(end_pos) + ".obs"
					output_names.append(output_name)
					os.rename("summary_stats_temp.txt", output_name)	
			
		# Should you wish to output the sumstats of the multiple windows and scaffolds to the same file (instead of one output file per region), you can include the code below. Default setting as is.
		with open (output_prefix + "_start" + str(across_pop_sites_all[0]) + "_end" + str(across_pop_sites_all[-1]) + "_windowsize" + str(window_size) + "_M" + str(min_allele_count) + "_Q" + str(quality_filter) + ".obs", "w") as master_file:
			for file in output_names[:1]:    
				f = open(file, "r")
				master_file.write(f.read())
				master_file.writelines("\n")
				f.close
				os.remove(file)
			for file in output_names[1:]:    
				f = open(file, "r")
				sumstat_line = f.read().splitlines()[1] + "\n"
				master_file.writelines(sumstat_line)
				f.close
				os.remove(file)
		with open(output_prefix + "_start" + str(across_pop_sites_all[0]) + "_end" + str(across_pop_sites_all[-1]) + "_windowsize" + str(window_size) + "_M" + str(min_allele_count) + "_Q" + str(quality_filter) + ".scaffoldlist", "w") as names_file:
			for scaffold_name in output_names:
				names_file.write("%s\n" % scaffold_name)

	else:
		print("Window size must be defined as 0 or a positive integer")	

	endTime = datetime.now()
	total_runtime = "Total runtime:" + str(endTime - startTime)
	
	print(total_runtime)
	
#%%
