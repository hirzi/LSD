#!/usr/bin/python

# sumstatscalc version
version = "sumstatscalc version: 08_03_2019"

############################################## BELOW WE DEFINE FUNCTIONS THAT CALCULATE RELEVANT SUMMARY STATISTICS ##############################################

# These functions take the output of the lsd_high function (allele frequencies) as the input #

# Load dependencies and set working directory/path

# import os
import numpy
import scipy
import scipy.stats
import matplotlib.pyplot as plt          # for plotting SFS (required for plot_SFS function)
from numpy import euler_gamma            # for calculating the harmonic mean using the digamma method (optionally required for the mean_Watterson_theta and mean_Tajima_D functions)
from scipy.special import digamma        # for calculating the harmonic mean using the digamma method (optionally required for the mean_Watterson_theta and mean_Tajima_D functions)
from fractions import Fraction           # for calculating the harmonic mean using the classical method (optionally required for the mean_Watterson_theta and mean_Tajima_D functions)
#from datetime import datetime            # for getting time stamps

# os.chdir("/Users/luqman/Documents/Hirzi/ETHPHD/Analysis/PoolSimulator")

#%%

############################################## NOTES: SAMPLE REQUIREMENTS ##############################################

#1) Tajima's D requires there to be >= 4 individuals per population.
#2) Variant class frequencies (from the SFS, e.g. singletons, doubletons, tripletons etc) require n individuals (polarised data) or 2n individuals (polarised data) per population for calculation, where n is the variant class frequency desired. By default, only doubletons and tripletons are outputted (singletons are by default excluded because they are too sensitive to error). Adjust accordingly depending on needs. 
#3) Negative FST and DXY values are collapsed to 0.

#%%

############################################## FUNCTIONS (GENERAL) ##############################################
# THIS FUNCTION GIVES THE NUMBER OF POLYMORPHIC LOCI (ITERATIONS). I.E. IT REMOVES FIXED (UNINFORMATIVE) LOCI FROM THE SIMULATION AND REDEFINES THE NUMBER OF ITERATIONS

def redefine_no_iterations(data):
	no_iterations = len(data["Allele_freqs"])
	return no_iterations


#%%	

############################################## FUNCTIONS (GENERAL STATISTICS) ##############################################
# THESE FUNTIONS HAVE BEEN CODED TO WORK WITH SINGLE POPULATION (pop = "no"), MULTIPLE POPULATION (pop = "yes") AND POPULATION PAIR (e.g. for Fst; pop = "pairwise") POPULATION STATISTICS.

## Calculate Mean ##
def mean(data, pop = "yes", axis = 0):
	if pop == "no":	
		mean = numpy.mean(data, axis)
		return mean
	elif pop == "yes":
		mean_array = {}
		for pop in range(1,len(data)+1,1):
			pop_data = data["Population" + str(pop)]
			mean = ["NA" if isinstance(pop_data, str) == True else numpy.mean(pop_data, axis)]
			mean_array["Population" + str(pop)] = mean[0]
		return mean_array
	elif pop == "pairwise":
		mean_array = {}
		no_pops = int((1+(1+8*len(data))**0.5) / 2)   # this formula (positive quadratic root) gives the number of populations from the dimensions (length) of the data. We do this so that we don't have to add another argument in the function
		for i in range(1,no_pops+1,1):
			for j in range(i+1,no_pops+1,1):
				pop_data = data["Population" + str(i) + "_" + str(j)]
				mean = numpy.mean(pop_data, axis)
				mean_array["Population" + str(i) + "_" + str(j)] = mean
		return mean_array		
		
## Calculate Variance ##		
def variance(data, pop = "yes"):	
	if pop == "no":	
		mean_ = mean(data, pop = "no")
		var = (sum([(i - mean_)**2 for i in data])) / (len(data) - 1)	
		return var
	elif pop == "yes":
		var_array = {}
		for pop in range(1,len(data)+1,1):
			pop_data = data["Population" + str(pop)]
			mean_array = mean(data, pop = "yes")
			mean_ = mean_array["Population" + str(pop)]
			if isinstance(pop_data, str) == True or isinstance(mean, str) == True:
				var = "NA"
			else:			
				var = (sum([(i - j)**2 for i,j in zip(pop_data, mean_)])) / (len(pop_data) - 1)
			var_array["Population" + str(pop)] = var
		return var_array
	elif pop == "pairwise":
		var_array = {}
		no_pops = int((1+(1+8*len(data))**0.5) / 2)   # this formula (positive quadratic solution) gives the number of populations from the dimensions (length) of the data. We do this so that we don't have to add another argument in the function
		for i in range(1,no_pops+1,1):
			for j in range(i+1,no_pops+1,1):
				pop_data = data["Population" + str(i) + "_" + str(j)]
				mean_array = mean(data, pop = "pairwise")
				mean_ = mean_array["Population" + str(i) + "_" + str(j)]
				var = (sum([(i - mean_)**2 for i in pop_data])) / (len(pop_data) - 1)
				var_array["Population" + str(i) + "_" + str(j)] = var
		return var_array

## Calculate Standard Deviation ##
def standard_deviation(data, pop = "yes"):
	if pop == "no":
		sd = (variance(data, pop = "no"))**0.5
		return sd
	elif pop == "yes":
		sd_array = {}
		for pop in range(1,len(data)+1,1):		
			pop_data = data["Population" + str(pop)]
			if isinstance(pop_data, str) == True:
				sd = "NA"
			else:				
				sd = (variance(pop_data, pop = "no"))**0.5
			sd_array["Population" + str(pop)] = sd
		return sd_array
	elif pop == "pairwise":
		sd_array = {}
		no_pops = int((1+(1+8*len(data))**0.5) / 2)   # this formula (positive quadratic solution) gives the number of populations from the dimensions (length) of the data. We do this so that we don't have to add another argument in the function
		for i in range(1,no_pops+1,1):
			for j in range(i+1,no_pops+1,1):		
				pop_data = data["Population" + str(i) + "_" + str(j)]		
				sd = (variance(pop_data, pop = "no"))**0.5
				sd_array["Population" + str(i) + "_" + str(j)] = sd
		return sd_array		

## Calculate Quartiles ##
def quartiles(data, pop = "yes"):
	if pop == "no":	
		quartiles = numpy.percentile(data, numpy.arange(25,100,25))   # to get the 1st, 2nd and 3rd quartiles
		return quartiles
	elif pop == "yes":
		quartiles_array = {}
		for pop in range(1,len(data)+1,1):
			pop_data = data["Population" + str(pop)]
			if isinstance(pop_data, str) == True:
				quartiles = "NA"
			else:	
				quartiles = numpy.percentile(pop_data, numpy.arange(25,100,25))
			quartiles_array["Population" + str(pop)] = quartiles	
		return quartiles_array	
	elif pop == "pairwise":
		quartiles_array = {}
		no_pops = int((1+(1+8*len(data))**0.5) / 2)   # this formula (positive quadratic solution) gives the number of populations from the dimensions (length) of the data. We do this so that we don't have to add another argument in the function
		for i in range(1,no_pops+1,1):
			for j in range(i+1,no_pops+1,1):
				pop_data = data["Population" + str(i) + "_" + str(j)]	
				quartiles = numpy.percentile(pop_data, numpy.arange(25,100,25))
				quartiles_array["Population" + str(i) + "_" + str(j)] = quartiles	
		return quartiles_array	

## Calculate Outliers (based on IQR method) ##		
def outliers(data,pop = "yes", outlier_coeff = 1.5, output = "count"):
	if pop == "no":	
		quartiles_ = quartiles(data,pop = "no")
		IQR = quartiles_[2] - quartiles_[0]  # inter-quartile range
		lower_bound = quartiles_[0] - (outlier_coeff * IQR)
		upper_bound = quartiles_[2] + (outlier_coeff * IQR)
		outliers = {}	
		for count,a in enumerate(data):
			if a > upper_bound or a < lower_bound:
				outliers["locus_" + str(count)] = a
		outliers_count = len(outliers)
		if output == "count":
			return outliers_count
		elif output == "list":
			return outliers
	elif pop == "yes":
		outliers_array = {}
		outliers_count = {}
		for pop in range(1,len(data)+1,1):		
			pop_data = data["Population" + str(pop)]
			if isinstance(pop_data, str) == True:
				outliers = "NA"
			else:				
				quartiles_ = quartiles(pop_data,pop = "no")
				IQR = quartiles_[2] - quartiles_[0]  # inter-quartile range
				lower_bound = quartiles_[0] - (outlier_coeff * IQR)
				upper_bound = quartiles_[2] + (outlier_coeff * IQR)
				outliers = {}	
				for count,a in enumerate(pop_data):
					if a > upper_bound or a < lower_bound:
						outliers["locus_" + str(count)] = a
			outliers_array["Population" + str(pop)] = outliers
			outliers_count["Population" + str(pop)] = len(outliers)
		if output == "count":
			return outliers_count
		elif output == "list":
			return outliers_array
	elif pop == "pairwise":
		outliers_array = {}
		outliers_count = {}		
		no_pops = int((1+(1+8*len(data))**0.5) / 2)   # this formula (positive quadratic solution) gives the number of populations from the dimensions (length) of the data. We do this so that we don't have to add another argument in the function
		for i in range(1,no_pops+1,1):
			for j in range(i+1,no_pops+1,1):		
				pop_data = data["Population" + str(i) + "_" + str(j)]			
				quartiles_ = quartiles(pop_data,pop = "no")
				IQR = quartiles_[2] - quartiles_[0]  # inter-quartile range
				lower_bound = quartiles_[0] - (outlier_coeff * IQR)
				upper_bound = quartiles_[2] + (outlier_coeff * IQR)
				outliers = {}	
				for count,a in enumerate(pop_data):
					if a > upper_bound or a < lower_bound:
						outliers["locus_" + str(count)] = a
				outliers_array["Population" + str(i) + "_" + str(j)] = outliers
				outliers_count["Population" + str(i) + "_" + str(j)] = len(outliers)				
		if output == "count":
			return outliers_count
		elif output == "list":
			return outliers_array

## Calculate Confidence Intervals ##
def confidence_intervals(data, pop = "yes", confidence = 0.95):
	if pop == "no":	
		ci = scipy.stats.norm.interval(confidence, loc=numpy.mean(data), scale=scipy.stats.sem(data))
		return ci
	elif pop == "yes":
		ci_array = {}		
		for pop in range(1,len(data)+1,1):		
			pop_data = data["Population" + str(pop)]
			if isinstance(pop_data, str) == True:
				norm_sem = "NA"
			else:					
				norm_sem = scipy.stats.norm.interval(confidence, loc=numpy.mean(pop_data), scale=scipy.stats.sem(pop_data))
			ci_array["Population" + str(pop)] = norm_sem		
		return ci_array		
	elif pop == "pairwise":
		ci_array = {}
		no_pops = int((1+(1+8*len(data))**0.5) / 2)   # this formula (positive quadratic solution) gives the number of populations from the dimensions (length) of the data. We do this so that we don't have to add another argument in the function
		for i in range(1,no_pops+1,1):
			for j in range(i+1,no_pops+1,1):		
				pop_data = data["Population" + str(i) + "_" + str(j)]				
				norm_sem = scipy.stats.norm.interval(confidence, loc=numpy.mean(pop_data), scale=scipy.stats.sem(pop_data))
				ci_array["Population" + str(i) + "_" + str(j)] = norm_sem		
		return ci_array	
		

#%%

############################################## FUNCTIONS (POPULATION STATISTICS) ##############################################


### Return the indexes of sites which are invariant, following/due to the stochastic pooling processes (recall that stochasiticy from the pooling simulation can lead to some originally variant sites being returned as invariant (across all populations) ###
# input: raw_allele_freqs array \ output: list of invariant site indexes  #
# Example usage: index_leftover_invariant_sites = index_remaining_invariant_sites (raw_allele_freqs, no_iterations, pop_size)
def index_remaining_invariant_sites (raw_allele_freqs, no_iterations, pop_size):
	indexes_to_remove_array = {}
	for iteration in range(1,no_iterations+1,1):	
		indexes_to_remove = []
		for index, site in enumerate(range(len(raw_allele_freqs["Allele_freqs"]["Iteration_" + str(iteration)]["Population1"]))):  # We take Population 1 here because all population within an iteration anyway have the same length, so doesn't matter
			sum_site = 0
			for pop in range(1,len(pop_size)+1,1):
				site_AF = raw_allele_freqs["Allele_freqs"]["Iteration_" + str(iteration)]["Population" + str(pop)][site]
				sum_site += site_AF
			if sum_site == 0:
				indexes_to_remove.append(site)
		indexes_to_remove_array["Iteration_" + str(iteration)] = indexes_to_remove
	return indexes_to_remove_array


### DISCARD SITES that don't meet the conditions in the lsd_high function (e.g. minimum allele frequency) and output retained sites ###
# input: raw_allele_freqs array, index_leftover_invariant_sites \ output: retained_allele_freqs array  #
# Example usage: retained_allele_freqs = discard_sites (raw_allele_freqs, index_leftover_invariant_sites, no_iterations, pop_size)
def discard_sites (raw_allele_freqs, index_leftover_invariant_sites, no_iterations, pop_size):
	results = {}   # to be consistent with lsd_high output. For consistency in downstream functions	
	retained_allele_freqs_master_array = {}
	for iteration in range(1,no_iterations+1,1):	
		retained_allele_freqs_pop_array = {}
		for pop in range(1,len(pop_size)+1,1):
			retained_allele_freqs = []		
			for index, site in enumerate(range(len(raw_allele_freqs["Allele_freqs"]["Iteration_" + str(iteration)]["Population" + str(pop)]))):
				all_indexes_to_remove = sorted(list(set((raw_allele_freqs["index_discarded_sites"]["Iteration_" + str(iteration)]) + (index_leftover_invariant_sites["Iteration_" + str(iteration)]))))
				if index not in all_indexes_to_remove:
					retained_sites = raw_allele_freqs["Allele_freqs"]["Iteration_" + str(iteration)]["Population" + str(pop)][site]
					retained_allele_freqs.append(retained_sites)
			retained_allele_freqs_pop_array["Population" + str(pop)] = retained_allele_freqs
		retained_allele_freqs_master_array["Iteration_" + str(iteration)] = retained_allele_freqs_pop_array
	results["Allele_freqs"] = retained_allele_freqs_master_array
	return results


### Count the number of discarded sites ###
# input: raw_allele_freqs array, index_leftover_invariant_sites \ output: count of discarded sites  #
# Example usage: discarded_sites = no_discarded_sites (raw_allele_freqs, index_leftover_invariant_sites, no_iterations, pop_size)
def no_discarded_sites (raw_allele_freqs, index_leftover_invariant_sites, no_iterations, pop_size):
	mean_no_discarded_sites_array = {}	
	list_no_discarded_sites = []
	for iteration in range(1,no_iterations+1,1):
		all_indexes_removed = list(set((raw_allele_freqs["index_discarded_sites"]["Iteration_" + str(iteration)]) + (index_leftover_invariant_sites["Iteration_" + str(iteration)])))
		len_discarded_sites = len(all_indexes_removed)
		list_no_discarded_sites.append(len_discarded_sites)	
	for pop in range(1,len(pop_size)+1,1):
		mean_no_discarded_sites_array["Population" + str(pop)] = list_no_discarded_sites
	return mean_no_discarded_sites_array
	

### Calculate nucleotide diversity PER SITE, for each locus, Pi ###
# To get total nucleotide diversity, just multiply nucleotide diversity per site with seq_length.
# input: allele_freqs array \ output: locus pi array #
# Example usage: locus_pi = pi_calc(retained_allele_freqs, no_iterations, pop_size, seq_length)
# For biallelic SNPs: nucleotide diversity Pi per site = SUM [ p * q * 2 * (n/(n-1)) ] / seq_length; where n/n-1 is the correction factor, p and q are allele frequencies, and q = 1-p (Nei & Li (1979); Gillespie (2004) Population Genetics: A concise guide. 2nd Ed., p. 45 top; Hamilton (2009) Population Genetics. pg 250 - though this book doesn't have the 2x coefficient of Gillespie - double check)
def pi_calc(allele_freqs, no_iterations, pop_size, seq_length):
	list_pi_array = {}	
	for pop in range(1,len(pop_size)+1,1):
		list_pi_iterations = []
		for iteration in range(1,no_iterations+1,1):	
			sum_pi = 0		
			for site in range(len(allele_freqs["Allele_freqs"]["Iteration_" + str(iteration)]["Population" + str(pop)])):
				derived_allele_freq = allele_freqs["Allele_freqs"]["Iteration_" + str(iteration)]["Population" + str(pop)][site]
				pi_per_site = 2 * derived_allele_freq * (1 - derived_allele_freq) * (pop_size[pop-1] / (pop_size[pop-1] - 1)) # where n in the original formulation is the number of alleles or in our case the number of samples (chromosomes) in a population
				sum_pi += pi_per_site
			mean_pi_per_site = sum_pi / seq_length # we divide by entire sequence length and not just by number of segregating sites, to get unbiased estimate of pi
			list_pi_iterations.append(mean_pi_per_site)
		list_pi_array["Population" + str(pop)] = list_pi_iterations			
	return list_pi_array

	
### Calculate nucleotide diversity PER SITE, for all sites, Pi ###
# input: allele_freqs array \ output: site pi array #
# Example usage: locus_site_pi = site_pi_calc(retained_allele_freqs, no_iterations, pop_size, output = "locus")
# For biallelic SNPs: nucleotide diversity Pi per site = SUM [ p * q * 2 * (n/(n-1)) ] / seq_length; where n/n-1 is the correction factor, p and q are allele frequencies, and q = 1-p (Nei & Li (1979); Gillespie (2004) Population Genetics: A concise guide. 2nd Ed., p. 45 top; Hamilton (2009) Population Genetics. pg 250 - though this book doesn't have the 2x coefficient of Gillespie - double check)
def site_pi_calc(allele_freqs, no_iterations, pop_size, output_ = "locus"):
	pi_site_array = {}
	for pop in range(1,len(pop_size)+1,1):
		list_pi_iterations = []
		dict_pi_iterations = {}
		for iteration in range(1,no_iterations+1,1):	
			list_pi_per_site = []
			for site in range(len(allele_freqs["Allele_freqs"]["Iteration_" + str(iteration)]["Population" + str(pop)])):
				derived_allele_freq = allele_freqs["Allele_freqs"]["Iteration_" + str(iteration)]["Population" + str(pop)][site]
				pi_per_site = 2 * derived_allele_freq * (1 - derived_allele_freq) * (pop_size[pop-1] / (pop_size[pop-1] - 1)) # where n in the original formulation is the number of alleles or in our case the number of samples (chromosomes) in a population
				list_pi_per_site.append(pi_per_site)
			if output_ == "mean":	
				list_pi_iterations.append(list_pi_per_site)
			elif output_ == "locus":
				dict_pi_iterations["Iteration_" + str(iteration)] = list_pi_per_site
		if output_ == "mean":		
			mean_pi_across_iterations = [float(sum(col))/len(col) for col in zip(*list_pi_iterations)]    # This calculates the pi per site averaged across all iterations. see https://docs.python.org/3/tutorial/controlflow.html#tut-unpacking-arguments for how * works 
			pi_site_array["Population" + str(pop)] = mean_pi_across_iterations
		elif output_ == "locus":
			pi_site_array["Population" + str(pop)] = dict_pi_iterations
	return pi_site_array
	
	
### Plot pi across all sites ###
# First identify WHAT you want to plot. E.g. all sites, retained sites, private sites, etc?
# input: site_pi \ output: plot of pi per site #
# Example usage: site_pi_mean = site_pi_calc(retained_allele_freqs, no_iterations, pop_size, output = "mean")
# Example usage: Pi_plot = plot_pi (site_pi, pop_size, file_format = "png")
def plot_pi (site_pi, pop_size, file_format):
	figures = []
	for pop in range(1,len(pop_size)+1,1):
		pi_site_pop = site_pi["Population" + str(pop)]
		figures.append(plt.figure())		
		plt.plot(pi_site_pop)
		my_xticks = list(range(len(pi_site_pop)))
		plt.xticks(numpy.arange(min(my_xticks), max(my_xticks)+1, 50))
		plt.ylabel("Site nucleotide diversity, pi")
		plt.xlabel("Position")
		plt.title("Pi per site; Population  " + str(pop))
		plt.savefig("Pi_Population" + str(pop) + "." + str(file_format), bbox_inches='tight', transparent = True)
	return plt.show()


### Calculate TOTAL number of segregating sites, S (INCLUDING erroneously called invariant sites), across all populations ###
# input: retained_allele_freqs array \ output: locus S array # 
# Example usage: locus_S = S_calc (retained_allele_freqs, pop_size, no_iterations)
def S_calc (retained_allele_freqs, pop_size, no_iterations):
	seg_sites_array = {}
	for pop in range(1,len(pop_size)+1,1):
		list_S_iterations = []
		for iteration in range(1,no_iterations+1,1):	
			seg_sites = len(retained_allele_freqs["Allele_freqs"]["Iteration_" + str(iteration)]["Population" + str(pop)])
			list_S_iterations.append(seg_sites)
		seg_sites_array["Population" + str(pop)] = list_S_iterations
	return seg_sites_array


### Calculate number of population-specific segregating sites, SPop. This is NOT the same as private segregating sites! ###
# input: allele_freqs array \ output: locus SPop array #
# Example usage: locus_SPop = SPop_calc (retained_allele_freqs, pop_size, no_iterations)
def SPop_calc (allele_freqs, pop_size, no_iterations):
	seg_sites_array = {}
	for pop in range(1,len(pop_size)+1,1):
		list_S_iterations = []
		for iteration in range(1,no_iterations+1,1):	
			seg_sites = numpy.count_nonzero(allele_freqs["Allele_freqs"]["Iteration_" + str(iteration)]["Population" + str(pop)])
			list_S_iterations.append(seg_sites)
		seg_sites_array["Population" + str(pop)] = list_S_iterations
	return seg_sites_array


### Calculate number of PRIVATE segregating sites, SP ###
# input: allele_freqs array \ output: locus SP array #
# Example usage: locus_SP = SP_calc (retained_allele_freqs, pop_size, no_iterations)
def SP_calc (allele_freqs, pop_size, no_iterations):
	raw_private_S_array = {}	
	private_S_array={}		
	for iteration in range(1,no_iterations+1,1):		
		# We first convert the allele freqs which are in dict format into list of list format; where each list corresponds to the allele freqs of one population.	
		dictlist_keys =[]	
		dictlist_values=[]
		for key, value in allele_freqs["Allele_freqs"]["Iteration_" + str(iteration)].items():
			dictlist_values.append(value)
			dictlist_keys.append(key)	
		# Note that since dicts are inherently unordered, we need to order the list of lists according to the population numbers (keys), as follows:
		dictlist_values_sorted = []
		for i in range(len(pop_size)):
			pop_index = dictlist_keys.index("Population"+str(i+1))
			dictlist_values_sorted.append(dictlist_values[pop_index])	
		# We then zip the lists together, and find count the number of private segregating sites for each population
		private_S = {}	
		for pop in range(len(pop_size)):
			count = 0
			for site in zip(*dictlist_values_sorted):       # Rememeber, if you have a list of iterables, you can zip the iterables together using zip(*iterables)
		#		print(str(site[0])+", "+str(site[1])+", "+str(site[2]))
				#if (site[pop] != 0 and site[pop] == sum(site)) or (site[pop] == 0 and site.count(0) == 1) or (site[pop] == 1 and site[pop] == (sum(site) + 1)):       # OLD, DEPRECIATED: the 1st condition denotes sites (with allele freq) that are not 0 in the target population but are 0 in all other populations; the 2nd condition denotes sites (with allele freq) that are zero in the target population but non-zero in all other populations; and the 3rd conditions denotes sites (with allele freqs) that are 1 in the target population but not 1 in all other populations. Together,this encapsulates all conditions for a private segregating sites. 
				if (site[pop] != 0 and site[pop] == sum(site)) or (site[pop] != 1 and (site[pop] + len(pop_size) - 1) == sum(site)):       # NEW REVISED: the 1st condition denotes sites (with allele freq) that are not 0 in the target population but are 0 in all other populations; the 2nd condition denotes sites (with allele freq) that are not 1 in the target population but are 1 in all other populations.
		#			print(str(site[0])+", "+str(site[1])+", "+str(site[2])+ "; is private for pop " + str(pop+1))     # To check
					count+=1
			private_S["Population" + str(pop+1)] = count
		raw_private_S_array["Iteration_" + str(iteration)] = private_S			
	for pop in range(1,len(pop_size)+1,1):
		poplist_private_S=[]
		for iteration in range(1,no_iterations+1,1):	
			element_private_S = raw_private_S_array["Iteration_" + str(iteration)]["Population" + str(pop)]
			poplist_private_S.append(element_private_S)
		private_S_array["Population" + str(pop)] = poplist_private_S
	return private_S_array
	

### Calculate number of segregating sites, SRAW, as outputted from ms (raw) output ### 
# input: filename \ output: locus SRAW array #
# Example usage: locus_SRAW = SRAW_calc(filename, no_iterations, pop_size)
def SRAW_calc (filename, no_iterations, pop_size, start_line=6, lines_between_endblock_startblock = 4):
	list_segsites_array = {}
	list_no_segsites = []	
	lines_between_startblock_startblock = sum(pop_size) + lines_between_endblock_startblock   # Open the file, define structure of (iteration) blocks
	with open(filename, encoding = "utf-8") as f:
		msms_raw = f.read().splitlines()				
		for count,iteration in enumerate(range(0,(no_iterations*lines_between_startblock_startblock),lines_between_startblock_startblock)):   # This first loop iterates over the iteration blocks
			index = []   # We make an index (list) of the start positions (line number) of the population blocks within each iteration block in the file, including the final end position. We need this to iterate over the populations within each iteration block.
			for j in range(len(pop_size)+1): # +1 since we also need the position of the final end position of the last population
				index_element = start_line + iteration + sum(pop_size[:j])
				index.append(index_element)		
			no_segsites = int(msms_raw[iteration + start_line - 2] [10:])   # We acquire the number of segregating sites directly from the ms output file:
			list_no_segsites.append(no_segsites)		
	for pop in range(1,len(pop_size)+1,1):
		list_segsites_array["Population" + str(pop)] = list_no_segsites		
	return list_segsites_array	


### Calculate number of segregating sites, SI (EXCLUDING invariant sites) ###
# input: filename, raw_allele_freqs \ output: locus SI array #
# Example usage: locus_SI = SI_calc (filename, raw_allele_freqs, no_iterations, pop_size)
def SI_calc (filename, raw_allele_freqs, no_iterations, pop_size):
	SI_array = {}
	locus_SRAW = SRAW_calc(filename, no_iterations, pop_size)
	discarded_sites = no_discarded_sites (raw_allele_freqs, no_iterations, pop_size)
	for pop in range(1,len(pop_size)+1,1):
		SI = [x-y for x,y in zip(locus_SRAW["Population" + str(pop)], discarded_sites["Population" + str(pop)])]
		SI_array["Population" + str(pop)] = SI
	return SI_array


### Calculate the harmonic number, An ###
# Method 1 (using the digamma function) is the the fastest and is exact, however requires loading additional modules
# Example usage: An = digamma_H(n-1)
def digamma_H(n):
    return digamma(n + 1) + euler_gamma
			
# Method 2 is more straightforward to understand plus also exact (gives output in terms of fractions), but slightly slower:
# Example usage: An = H(n-1)
def H(n):
	return sum(Fraction(1, d) for d in range(1, n + 1))
	
### Calculate the 2nd power harmonic number ###
# Example usage: An2 = H2(n-1)
def H2(n):
	return sum(Fraction(1, d**2) for d in range(1, n + 1))
	

### Calculate Watterson's Theta ###
# input: locus S, SI, SPop or SP array \ output: locus Watterson's theta array #
# Example usage: locus_WattersonTheta = WattersonTheta_calc(locus_S, pop_size)
# Watterson's theta is given by the formula: Watterson_theta = S / An, where S is the number of segregating sites and An is the (n-1)th harmonic number (n is the number of samples (chromosomes) in a population). See http://stacSoverflow.com/questions/404346/python-program-to-calculate-harmonic-series for different ways to calculate the harmonic number (depending on desired levels of accuracy and speed)
def WattersonTheta_calc (S, pop_size):
	Watterson_theta_array = {}
	for pop in range(1,len(pop_size)+1,1):
		Watterson_theta = [x / digamma_H(pop_size[pop-1] - 1) for x in S["Population" + str(pop)]]
		Watterson_theta_array["Population" + str(pop)] = Watterson_theta		
	return Watterson_theta_array
	
	
### Calculate Tajima's D ###
# Let's define some identities to below to make things simpler (Source: Hamilton 2009. Population Genetics; https://en.wikipedia.org/wiki/Tajima%27s_D):
"""
e1 = (c1 / An)
e2 = (c2 / ((An**2) + An2))
c1 = (b1 - (1/An))
c2 = (b2 - ((n+2) / (An*n)) + (An2 / An**2))
b1 = ((n+1) / (3*(n-1)))
b2 = ((2*((n**2)+n+3)) / (9 * n * (n-1)))
An = (digamma_H(n-1))                             # this is the (n-1)th harmonic number
An = (H(n-1))                                     # this is the (n-1)th harmonic number (alternative method)
An2 = (H2(n-1))	                                 # this is the (n-1)th 2nd-order harmonic number
S = number of segregating sites                   # this is taken from the output of mean_S 			

Var_d = ((e1*S) + (e2*S*(S-1)))                   # this if the variance of Tajima's d
SD_d = (Var_d**0.5)                               # this is the standard deviation of Tajima's d

tajima_d = nucleotide_distance - theta_W          # Tajima's d
tajima_D = tajima_d / SD_d                        # Tajima's D (normalised Tajima's d)
"""
				
# We can write the above equations succinctly as:
"""				
An = (digamma_H(n-1))                               # this is the (n-1)th harmonic number
An = (H(n-1))                                       # this is the (n-1)th harmonic number (alternative method)
An2 = (H2(n-1))	                                   # this is the (n-1)th 2nd-order harmonic number	
SD_d = ((((((n+1) / (3*(n-1))) - (1/An)) / An) * S) + (((((2*((n**2)+n+3)) / (9 * n * (n-1))) - ((n+2) / (An*n)) + (An2 / An**2)) / ((An**2) + An2)) * S * (S-1)))**0.5
tajima_d = nucleotide_distance - theta_W
tajima_D = tajima_d / SD_d

"""

### Using the above equations, we calculate Tajima's D ###
# input: locus S, SI, SPop or SP array, locus Watterson theta array, locus pi array \ output: locus Tajima's D array #
# Example usage: locus_TajimaD = TajimaD_calc (locus_S, locus_WattersonTheta, locus_pi, pop_size, seq_length)
def TajimaD_calc (segsites, WattersonTheta, pi, pop_size, seq_length):
	Tajima_D_array = {}
	for pop in range(1,len(pop_size)+1,1):	
		theta_W = WattersonTheta["Population" + str(pop)]
		nucleotide_distance = [p * seq_length for p in pi["Population" + str(pop)]]   # we use total sequence pi here, not pi per site 
		S = segsites["Population" + str(pop)]
		n = pop_size[pop-1]		
		tajima_d = [x-y for x,y in zip(nucleotide_distance, theta_W)]
		if n > 3:     # Test of Tajima's D requires at least 4 samples per population. Samples of < 3 return a 0 value for SD_d.
			An = (H(n-1))
			An2 = (H2(n-1))
			SD_d = [((((((n+1) / (3*(n-1))) - (1/An)) / An) * s) + (((((2*((n**2)+n+3)) / (9 * n * (n-1))) - ((n+2) / (An*n)) + (An2 / An**2)) / ((An**2) + An2)) * s * (s-1)))**0.5  for s in S] 
			tajima_D = [float(t/sd) if s>0 else 0 for t, sd, s in zip(tajima_d, SD_d, S)]
			Tajima_D_array["Population" + str(pop)] = tajima_D
		else:
			Tajima_D_array["Population" + str(pop)] = "Not enough samples in population. Test of Tajima's D requires at least 4 samples per population." 
	return Tajima_D_array


##################### NOTE: FOR THE CALCULATION OF HS, HT & FST BELOW, ALWAYS USE RETAINED_ALLELE_FREQS! It doesn't make sense to use populationspecific_allele_freqs as this has different lengths! #####################

### Calculate the average expected heterozygosity of subpopulations, Hs ###
# input: locus pi array \ output: locus Hs #
# Example usage: locus_Hs = Hs_calc (pi = locus_pi, compare_pops = [1,3])
def Hs_calc (pi, compare_pops):   # mean_pi -> pi
	list_pi_pop = []
	for pop in compare_pops:
		pi_pop = pi["Population" + str(pop)]
		list_pi_pop.append(pi_pop)
	Hs = numpy.mean(list_pi_pop,axis=0)
	return Hs	


### Calculate the expected heterozygosity of the total population, Ht ###
# input: retained_allele_freqs array  \ output: locus Ht #
# Example usage: locus_Ht = Ht_calc (retained_allele_freqs, no_iterations, pop_size, seq_length, compare_pops = [1,3])
def Ht_calc (retained_allele_freqs, no_iterations, pop_size, seq_length, compare_pops):
	list_mean_pi_iterations = []
	for iteration in range(1,no_iterations+1,1):	
		counts_list = []
		for pop in compare_pops:
			allele_freq_pop = retained_allele_freqs["Allele_freqs"]["Iteration_" + str(iteration)]["Population" + str(pop)]
			count_pop = [x * pop_size[pop-1] for x in allele_freq_pop]   # original definition
			counts_list.append(count_pop)
		counts_array = numpy.array(counts_list)
		if len(counts_array[0]) == len(counts_array[1]): 
			counts_total = numpy.sum(counts_array, axis = 0)
			total_no_samples = sum ([pop_size[x-1] for x in compare_pops])  # we calculate the total number of samples within the populations defined in compare_pops
			allele_freq_total = [x / total_no_samples for x in counts_total]
		else:
			print("Error, the populations have different numbers of sites; homologous (one-to-one) site comparison not possible!")
	# now calculate pi
		sum_pi = 0		
		for site in range(len(allele_freq_total)):
			derived_allele_freq = allele_freq_total[site]
			pi_per_site = 2 * derived_allele_freq * (1 - derived_allele_freq) * (total_no_samples / (total_no_samples - 1)) # where n in the original formulation is the number of alleles or in our case the number of samples (chromosomes) in a population
			sum_pi += pi_per_site
		mean_pi_per_site = sum_pi / seq_length # we divide by entire sequence length and not just by number of invariant sites, to get unbiased estimate of pi
		list_mean_pi_iterations.append(mean_pi_per_site)
	return list_mean_pi_iterations


### Calculate Fst (CLASSICAL, UNWEIGHTED): Fst = (Ht - Hs) / Ht ### DOUBLE-CHECK, POTENTIALLY REQUIRES REVISION
## Note, this (classical) Fst assumes equal population sizes (?).
# Note: the mean of locus_Fst =/= mean Fst (because of the difference in the calculation of mean vs locus-by-locus Fst = (Ht - Hs) / Ht)
# input: Hs, Ht \ output: mean or locus-by-locus Fst #
# Example usage: locus_Fst = Fst_unweighted_locus_calc (locus_pi, locus_S, retained_allele_freqs, no_iterations, pop_size, seq_length, span = "segsites")
def Fst_unweighted_locus_calc (locus_pi, locus_S, retained_allele_freqs, no_iterations, pop_size, seq_length, span = "segsites"):
	Fst_dict = {}                                    # this list will be of length no_pops choose 2 once filled                          
	for i in range(1,len(pop_size)+1,1):
		for j in range(i+1,len(pop_size)+1,1):
			Hs = Hs_calc (locus_pi, compare_pops = [i,j])
			Ht = Ht_calc (retained_allele_freqs, no_iterations, pop_size, seq_length, compare_pops = [i,j])
			Fst_segsites = [float((y-x)/y) for x,y in zip(Hs,Ht)]   # this is the Fst across all segregating sites
			if span == "segsites":
				Fst_dict["Population" + str(i) + "_" + str(j)] = Fst_segsites   # note the filling order: pop1pop2->pop1pop3...->pop1popn->pop2pop3->pop2pop4->...pop2popn->...popn-1popn				
			elif span == "allsites":
				no_segsites = locus_S["Population1"]  # this is the number of segregating sites. remember, all populations within an iteration have the same length
				Fst_all = [f * (s / seq_length) for f, s in zip(Fst_segsites, no_segsites)]    # this is the Fst across ALL sites
				Fst_dict["Population" + str(i) + "_" + str(j)] = Fst_all   # note the filling order: pop1pop2->pop1pop3...->pop1popn->pop2pop3->pop2pop4->...pop2popn->...popn-1popn
		## To output the results to a triangular matrix 		
		#	Fst_matrix = numpy.zeros((len(pop_size),len(pop_size)))      # empty matrix of dimension no_pops X no_pops
		#	indices = numpy.tril_indices(len(pop_size),-1)         # define matrix filling indices (offset the diagonal by -1 (since we want the main diagonal to comprise of 0's))
		#	Fst_matrix[indices] = Fst_list	                # fills the matrix with the Fst values accordingly
	return Fst_dict


### Calculate Fst (Weir & Cockenheim 1984) ###
# Let's define some identities to below to make things simpler (Source: Weir & Cockenheim 1984 pg 1359-1360):
"""
n_i = sample size of population i
r = no_pops = len(pop_size) (in our case, since we're always dealing with pairwise comparisons, r = 2)
h_i = average heterozygosity (heterozygote frequency or alternatively pi) for allele A in population i
h_hat = average heterozygosity of allele A across populations (samples) = SUM(over i)n_i*h_i / r*n_hat
n_hat = average sample size = sum(pop_size) / len(pop_size) = sum(pop_size) / 2
n_c = ((r * n_hat) - (sum([x**2 for x in pop_size]) / (r*n_hat))) / (r - 1)
p_i = average frequency of allele A in population i
p_hat = average frequency of allele A across populations (samples) = SUM(over i)n_i*p_i / r*n_hat
s_squared = sample variance of allele A frequencies over populations = SUM(over i) n_i*(p_i - p_hat)**2 / (r - 1)*n_hat

a = component of variance between subpops
b = component of variance between individuals within subpops
c = component of variance between gametes within individuals

a = (n_hat/n_c) * (s_square - ((1/(n_hat - 1)) * ((p_hat*(1-p_hat)) - (((r-1)/r)*s_square) - (0.25*h_hat))))
b = (n_hat/(n_hat - 1)) * ((p_hat*(1-p_hat)) - (((r-1)/r)*s_square) - ((((2*n_hat) - 1)/(4*n_hat))*h_hat))
c = 0.5*h_hat
theta = a / (a + b + c)   # theta = Fst
"""

### Calculate Fst PER SITE (Weir-Cockenheim 1984) ###
# This code is modified from Mathias Scharmann 2016 (ms2stats.stats.py), who in turn modified from Eva Chan 2008 (calc_wcFstats.R), and is based in the equations above
# Negative estimates are returned as 0.0, because that is what they mean.
# We use pi in the place of heterozygosity, in light that nucleotide diversity pi and expected heterozygosity are similar, and because ms outputs haploid samples (not diploids)
# input: allele frequencies, pi, pop_size \ output: site Fst #
# Example usage: Fst_WC_site = Fst_WC_site_calc (p_i, p_j, n_i, n_j, h_i, h_j)
def Fst_WC_site_calc (p_i, p_j, n_i, n_j, h_i, h_j):	
	# allele freq of allele A in pop i
	p_i = float(p_i)
	p_j = float(p_j)
	## sample size in pop i
	n_i = float(n_i)
	n_j = float(n_j)
	## number of populations; here of course 2
	r = 2.0 
	## observed proportion of individuals heterozygous for allele A in pop i
	h_i = float(h_i)
	h_j = float(h_j)	
	## average sample size
	n = (n_i + n_j) / r  
	## scv = squared coefficient of variation of sample sizes = square of (standard deviation / mean)
#	scv = ( ( math.sqrt( (1/r)*( (n_i-n)**2 + (n_j-n)**2 ) ) )/n )**2
#	scv = ( ( (n_i**2) + (n_j**2) ) - (n*n*r) ) / ( (n*n) * (r-1) ) # Eva Chan 2008 calc_wcFstats.R 	
#	n_c = n*(1-(scv/r)) 	
	n_c = ((r*n) - ( ((n_i*n_i)/(r*n)) + ((n_j*n_j)/(r*n)) ) ) / (r - 1.0) # Eva Chan 2008 calc_wcFstats.R	
	## average sample frequency of allele A
	p = (n_i*p_i)/(r*n) + (n_j*p_j)/(r*n)
	## sample variance of allele A frequnecies over populations
	s_2 = (n_i*((p_i-p)**2)) / ((r-1.0)*n) + (n_j*((p_j-p)**2)) / ((r-1.0)*n)		
	## the average heterozygote frequency for allele A
	h = (n_i*h_i)/(r*n) + (n_j*h_j)/(r*n)
	# now the "three quantities"
	a = (n/n_c) * ( s_2 - ((1.0/(n-1.0))*((p*(1.0-p)) - (((r-1.0)/r)*s_2) - ((1.0/4.0)*h))) )
	b = (n/(n-1.0)) * ((p*(1.0-p)) - (((r-1.0)/r)*s_2) - ((((2.0*n)-1.0)/(4.0*n))*h))
	c = h/2.0 
	# finally, the Fst estimator theta; it is undefined for (a + b + c = 0) caused by fixation of the same allele in all samples
	if (a + b + c) == 0.0:
#		theta_hat = float("-nan")
		theta_hat = 0.0   # since theta is 0.0 if same allele in all samples
	else:
		theta_hat = a / (a + b + c)	
	return theta_hat
	

### Calculate Fst PER LOCUS (Weir-Cockenheim 1984) ###
# Returns Weir & Cockerham's 1984 Fst estimator theta_hat averaged over the number of SNPs per locus
# Negative estimates are returned as 0.0, because that is what they mean.
# Example usage: locus_Fst_WC = Fst_WC_locus_calc(retained_allele_freqs, locus_site_pi, pop_size, no_iterations, seq_length, span = "segsites")
def Fst_WC_locus_calc (retained_allele_freqs, locus_site_pi, pop_size, no_iterations, seq_length, span = "segsites"):	
	Fst_dict = {} 	
	for i in range(1,len(pop_size)+1,1):
		for j in range(i+1,len(pop_size)+1,1):
			Fst_WC_locus_array = []
			for iteration in range(1,no_iterations+1,1):	 
				Fst_sum = 0.0
				for site in range(len(retained_allele_freqs["Allele_freqs"]["Iteration_" + str(iteration)]["Population1"])):   # All populations within an iteration have the same length, hence we can use any here
					p_A = retained_allele_freqs["Allele_freqs"]["Iteration_" + str(iteration)]["Population" + str(i)][site]
					p_B = retained_allele_freqs["Allele_freqs"]["Iteration_" + str(iteration)]["Population" + str(j)][site]
					n_A = pop_size[i-1]
					n_B = pop_size[j-1]
					h_A = locus_site_pi["Population" + str(i)]["Iteration_" + str(iteration)][site]
					h_B = locus_site_pi["Population" + str(j)]["Iteration_" + str(iteration)][site]	
					Fst_WC_site = Fst_WC_site_calc (p_A, p_B, n_A, n_B, h_A, h_B)
					Fst_sum += Fst_WC_site
				if span == "segsites" and len(retained_allele_freqs["Allele_freqs"]["Iteration_" + str(iteration)]["Population1"]) > 0:
					Fst_locus = Fst_sum / len(retained_allele_freqs["Allele_freqs"]["Iteration_" + str(iteration)]["Population1"])	
				elif span == "segsites" and len(retained_allele_freqs["Allele_freqs"]["Iteration_" + str(iteration)]["Population1"]) == 0:   # in the case that all variant sites are removed through filtering, to avoid zero division error
					Fst_locus = 0
				elif span == "allsites":	
					Fst_locus = Fst_sum / seq_length		# we divide by entire sequence length and not just by number of segregating sites, to get across all sites (within locus) estimate of Fst		
				Fst_WC_locus_array.append([Fst_locus if Fst_locus >= 0.0 else 0.0])			
			Fst_dict["Population" + str(i) + "_" + str(j)] = [item for sublist in Fst_WC_locus_array for item in sublist]    # note the filling order: pop1pop2->pop1pop3...->pop1popn->pop2pop3->pop2pop4->...pop2popn->...popn-1popn		
	return Fst_dict


### Calculate Dxy PER LOCUS ###
# Formula taken from Nei & Li 1979 (between eqs. 24 and 25), Nei 1987 eq 10.20, Cruickshank & Hahn 2014 (Box 1)
# for a single biallelic SNP (1,2) in two pops it simplifies to: Dxy = pop1_freq1 * ( 1.0 - pop2_freq1 ) + pop2_freq1 * ( 1.0 - pop1_freq1 )
# Dxy can take values between 0 and 2		
# Example usage: locus_Dxy = Dxy_calc (retained_allele_freqs, pop_size, no_iterations, seq_length, span = "segsites")
#def Dxy_calc (retained_allele_freqs, pop_size, no_iterations, seq_length, span = "segsites"):
def Dxy_calc (retained_allele_freqs, pop_size, no_iterations, seq_length):	
	Dxy_dict = {} 	
	for i in range(1,len(pop_size)+1,1):
		for j in range(i+1,len(pop_size)+1,1):
			Dxy_array = []
			for iteration in range(1,no_iterations+1,1):	 
				Dxy_sum = 0.0
				for site in range(len(retained_allele_freqs["Allele_freqs"]["Iteration_" + str(iteration)]["Population1"])):   # All populations within an iteration have the same length, hence we can use any here					
					p_i = retained_allele_freqs["Allele_freqs"]["Iteration_" + str(iteration)]["Population" + str(i)][site]
					p_j = retained_allele_freqs["Allele_freqs"]["Iteration_" + str(iteration)]["Population" + str(j)][site]
					Dxy = (p_i * (1.0 - p_j)) + (p_j * (1.0 - p_i))
					Dxy_sum += Dxy
				#if span == "segsites" and len(retained_allele_freqs["Allele_freqs"]["Iteration_" + str(iteration)]["Population1"]) > 0:
				#	Dxy_locus = Dxy_sum / len(retained_allele_freqs["Allele_freqs"]["Iteration_" + str(iteration)]["Population1"])	
				#elif span == "segsites" and len(retained_allele_freqs["Allele_freqs"]["Iteration_" + str(iteration)]["Population1"]) == 0:   # in the case that all variant sites are removed through filtering, to avoid zero division error
				#	Dxy_locus = 0
				#elif span == "allsites":
				# Recall that Dxy is formally calculated over all sites, not just segregating sites.
				Dxy_locus = Dxy_sum / seq_length
				Dxy_array.append([Dxy_locus if Dxy_locus >= 0.0 else 0.0])			
			Dxy_dict["Population" + str(i) + "_" + str(j)] = [item for sublist in Dxy_array for item in sublist]    # note the filling order: pop1pop2->pop1pop3...->pop1popn->pop2pop3->pop2pop4->...pop2popn->...popn-1popn		
	return Dxy_dict	
	

### Calculate site frequency spectrum (SFS) from allele frequencies ###
# input: retained_allele_freqs array \ output: SFS array #
# Example usage: SFS = SFS_calc(retained_allele_freqs, pop_size, no_iterations, polarised)
def SFS_calc (allele_freqs, pop_size, no_iterations, polarised):
	SFS_master_array = {}
	for iteration in range(1,no_iterations+1,1):			
		SFS_pop_array = {}		
		for pop in range(1,len(pop_size)+1,1):		
			SFS_vector_raw = []
			for c in range(1,pop_size[pop-1],1):   # We go from 1 to n-1, since we consider only variant sites
				allele_frequencies = allele_freqs["Allele_freqs"]["Iteration_" + str(iteration)]["Population" + str(pop)]
				counts = [round(x * (pop_size[pop-1])) for x in allele_frequencies] # round to bin fractions to nearest integer
				SFS_count = counts.count(c)
				SFS_vector_raw.append(SFS_count)
			if polarised == False:   # We fold, by mirroring and summing
				folded_full_length = list(numpy.array(SFS_vector_raw) + numpy.array(SFS_vector_raw[::-1])) 
				SFS_vector = folded_full_length[:int(len(SFS_vector_raw)/2)] # Then we cut in half. Recall, int() rounds down a fraction
				if (len(SFS_vector_raw)/2).is_integer() == False:   # In the case of odd-sized populations, we append the middle (odd) element (the last element in the folded SFS)
					SFS_vector.append(SFS_vector_raw[int(len(SFS_vector_raw)/2)])   
			if polarised == True:				
				SFS_vector = SFS_vector_raw
			SFS_pop_array["Population" + str(pop)] = SFS_vector		
		SFS_master_array["Iteration_" + str(iteration)] = SFS_pop_array
	return SFS_master_array	
	
	
### Calculate mean SFS across all loci/iterations ###
# input: SFS array \ output: mean SFS array #
# Example usage: mean_SFS = mean_SFS_calc(SFS, pop_size, no_iterations)
def mean_SFS_calc (SFS, pop_size, no_iterations):		
	SFS_mean_array = {}			
	for pop in range(1,len(pop_size)+1,1):
		list_SFSiterations = []
		for iteration in range(1,no_iterations+1,1):
			list_SFSiterations.append(SFS["Iteration_" + str(iteration)]["Population" + str(pop)])
		mean_across_iterations = numpy.mean(list_SFSiterations, axis = 0)
		SFS_mean_array["Population" + str(pop)] = mean_across_iterations
	return SFS_mean_array


### Plot SFS ###
# input: mean SFS array \ output: plot of mean SFS #
# Example usage: SFS_plot = plot_SFS(mean_SFS, pop_size, file_format = "png")
# Consider if you want to set a min and max possible value for the counts. Min being 1 and max being no. of samples - 1 (since we can't have fixed segregating sites). For it theoretically shouldn't be able to take these values, though in practice, it is possible for the pooled data/simulation to take these values.	
def plot_SFS (mean_SFS, pop_size, file_format):
	figures = []	
	for pop in range(1,len(pop_size)+1,1):
		SFS_pop = mean_SFS["Population" + str(pop)].tolist()
		figures.append(plt.figure())		
		x_axis = range(1,len(SFS_pop)+1,1)
		y_axis = numpy.arange(len(SFS_pop))
		plt.bar(y_axis, SFS_pop, width = 0.9, align = "center")
		plt.xticks(y_axis, x_axis)
		plt.ylabel("Number of SNPS")
		plt.xlabel("Number of variants at site")
		plt.title("SFS Population " + str(pop))
		plt.savefig("SFS_" + str(pop) + "." + str(file_format), bbox_inches='tight', transparent = True)
	return plt.show()


############################################## THIS IS THE MAIN SUMMARY STATISTICS CALCULATOR ##############################################
# input: raw_allele_freqs \ output: sumstats_output.txt, fst_output.txt, SFS_plots.png/pdf #
# Example usage: sumstats_calc (raw_allele_freqs, pop_size, seq_length, plotSFS = True)
def sumstats_calc (raw_allele_freqs, pop_size, seq_length, span, polarised, output = "full", plotSFS = False, plotPi = False):
#def sumstats_calc (raw_allele_freqs, pop_size, seq_length, span, output = "full", plotSFS = False, plotPi = False):     # We remove/hash out estimates of SI (and thus the requirement to define a filename here) because SI cannot be calculated for observed data.

#	startTime = datetime.now()
	print(version)
	
	# We define the no_iterations
	no_iterations = redefine_no_iterations(raw_allele_freqs)
	
	# We calculate the population genetic summary statistics
	index_leftover_invariant_sites = index_remaining_invariant_sites (raw_allele_freqs, no_iterations, pop_size)
	retained_allele_freqs = discard_sites (raw_allele_freqs, index_leftover_invariant_sites, no_iterations, pop_size)
	discarded_sites = no_discarded_sites (raw_allele_freqs, index_leftover_invariant_sites, no_iterations, pop_size)
	locus_pi = pi_calc(retained_allele_freqs, no_iterations, pop_size, seq_length)
	locus_site_pi = site_pi_calc(retained_allele_freqs, no_iterations, pop_size, output_ = "locus")
#	locus_SI = SI_calc (filename, raw_allele_freqs, no_iterations, pop_size)
	locus_S = S_calc (retained_allele_freqs, pop_size, no_iterations)
	locus_SPop = SPop_calc (retained_allele_freqs, pop_size, no_iterations)
	locus_SP = SP_calc (retained_allele_freqs, pop_size, no_iterations)
	locus_WattersonTheta = WattersonTheta_calc(locus_SPop, pop_size)	
	locus_TajimaD = TajimaD_calc (locus_SPop, locus_WattersonTheta, locus_pi, pop_size, seq_length)
#	locus_Fst = Fst_unweighted_locus_calc (locus_pi, locus_S, retained_allele_freqs, no_iterations, pop_size, seq_length, span)
	locus_Fst_WC = Fst_WC_locus_calc(retained_allele_freqs, locus_site_pi, pop_size, no_iterations, seq_length, span)
	locus_Dxy = Dxy_calc (retained_allele_freqs, pop_size, no_iterations, seq_length)
	SFS = SFS_calc(retained_allele_freqs, pop_size, no_iterations, polarised)	
	mean_SFS = mean_SFS_calc(SFS, pop_size, no_iterations)
	
	
	# We calculate the mean, variance, quartiles, confidence intervals and outliers of the population genetic summary statistics     # We hash out calculations for SD because the observed data is a single iteration (thus SD between iterations doesn't exist)
	mean_pi = mean(locus_pi)
#	sd_pi = standard_deviation(locus_pi)
	mean_S = mean(locus_S)
#	sd_S = standard_deviation(locus_S)
	mean_SPop = mean(locus_SPop)
#	sd_SPop = standard_deviation(locus_SPop)	
	mean_SP = mean(locus_SP)
#	sd_SP = standard_deviation(locus_SP)
#	mean_SI = mean(locus_SI)
#	sd_SI = standard_deviation(locus_SI)
	mean_discarded_sites = mean(discarded_sites)
#	sd_discarded_sites = standard_deviation(discarded_sites)	
	mean_WattersonTheta = mean(locus_WattersonTheta)
#	sd_WattersonTheta = standard_deviation(locus_WattersonTheta)
	mean_TajimaD = mean(locus_TajimaD)
#	sd_TajimaD = standard_deviation(locus_TajimaD)
#	mean_Fst_unweighted = mean(locus_Fst, pop = "pairwise")
#	sd_Fst_unweighted = standard_deviation(locus_Fst, pop = "pairwise")
	mean_Fst_WC = mean(locus_Fst_WC, pop = "pairwise")
#	sd_Fst_WC = standard_deviation(locus_Fst_WC, pop = "pairwise")
	if no_iterations > 1:
		outliers_Fst_WC = outliers(locus_Fst_WC, pop = "pairwise")
	mean_Dxy = mean(locus_Dxy, pop = "pairwise")
#	sd_Dxy = standard_deviation(locus_Dxy, pop = "pairwise")	
		
#	endTime = datetime.now()
#	total_runtime = "total_runtime:" + str(endTime - startTime)
	
	if output == "full":   # Outputs full results, in a convenient, human-readable format
		# Define header and row headers:
		sumstats_header = ["no_pops:" + str(len(pop_size)), "no_loci:" + str(no_iterations), "pop_size:" + str(pop_size), "estimated_total_sequence_length:" + str(seq_length)]
#		sumstats_header = ["input_file:" + str(filename), "no_pops:" + str(len(pop_size)), "no_loci:" + str(no_iterations), "pop_size:" + str(pop_size), "estimated_total_sequence_length:" + str(seq_length), str(startTime), total_runtime]
		pops = ["Population"]
		mean_pi_output = ["Mean_pi"]
#		sd_pi_output = ["SD_pi"]		
		mean_S_output = ["Mean_S"]
#		sd_S_output = ["SD_S"]
		mean_SPop_output = ["Mean_SPop"]
#		sd_SPop_output = ["SD_SPop"]
		mean_SP_output = ["Mean_Private_S"]
#		sd_SP_output = ["SD_Private_S"]	
#		mean_SI_output = ["Mean_S_w/o_invariants"]
#		sd_SI_output = ["SD_S_w/o_invariants"]
		mean_discarded_sites_output = ["Mean_discarded_sites"]
#		sd_discarded_sites_output = ["SD_discarded_sites"]	
		mean_WattersonTheta_output = ["Mean_WattersonTheta"]
#		sd_WattersonTheta_output = ["SD_WattersonTheta"]
		mean_TajimaD_output = ["Mean_TajimasD"]
#		sd_TajimaD_output = ["SD_TajimasD"]
#		mean_singletons_output = ["Mean_singletons"]
		mean_doubletons_output = ["Mean_doubletons"]	
		mean_tripletons_output = ["Mean_tripletons"]	
		mean_SFS_output = ["Mean_SFS"]
	
		# Define the population (calculated summary statistic data) entries for the output file
		for pop in range(1,len(pop_size)+1,1):
			pops.append("Pop" + str(pop))		
			mean_pi_string = str(mean_pi["Population" + str(pop)])
			mean_pi_output.append(mean_pi_string)		
#			sd_pi_string = str(sd_pi["Population" + str(pop)])
#			sd_pi_output.append(sd_pi_string)	
					
			mean_S_string = str(mean_S["Population" + str(pop)])
			mean_S_output.append(mean_S_string)		
#			sd_S_string = str(sd_S["Population" + str(pop)])
#			sd_S_output.append(sd_S_string)
			
			mean_SPop_string = str(mean_SPop["Population" + str(pop)])
			mean_SPop_output.append(mean_SPop_string)		
#			sd_SPop_string = str(sd_SPop["Population" + str(pop)])
#			sd_SPop_output.append(sd_SPop_string)	
	
			mean_SP_string = str(mean_SP["Population" + str(pop)])
			mean_SP_output.append(mean_SP_string)		
#			sd_SP_string = str(sd_SP["Population" + str(pop)])
#			sd_SP_output.append(sd_SP_string)	
			
#			mean_SI_string = str(mean_SI["Population" + str(pop)])
#			mean_SI_output.append(mean_SI_string)		
#			sd_SI_string = str(sd_SI["Population" + str(pop)])
#			sd_SI_output.append(sd_SI_string)		
			
			mean_discarded_sites_string = str(mean_discarded_sites["Population" + str(pop)])
			mean_discarded_sites_output.append(mean_discarded_sites_string)		
#			sd_discarded_sites_string = str(sd_discarded_sites["Population" + str(pop)])
#			sd_discarded_sites_output.append(sd_discarded_sites_string)			
			
			mean_WattersonTheta_string = str(mean_WattersonTheta["Population" + str(pop)])
			mean_WattersonTheta_output.append(mean_WattersonTheta_string)		
#			sd_WattersonTheta_string = str(sd_WattersonTheta["Population" + str(pop)])
#			sd_WattersonTheta_output.append(sd_WattersonTheta_string)	
	
			mean_TajimaD_string = str(mean_TajimaD["Population" + str(pop)])
			mean_TajimaD_output.append(mean_TajimaD_string)		
#			sd_TajimaD_string = str(sd_TajimaD["Population" + str(pop)])
#			sd_TajimaD_output.append(sd_TajimaD_string)	
			
			mean_SFS_string = [float(x) for x in mean_SFS["Population" + str(pop)]]
			mean_SFS_output.append(str(mean_SFS_string))
#			mean_singletons_output.append(str(mean_SFS_string[0]))
			mean_doubletons_output.append(str(mean_SFS_string[1]))		
			mean_tripletons_output.append(str(mean_SFS_string[2]))
		
		# To write sumstats output file
		with open("summary_stats_temp.txt", "w") as sumstatsfile:
			header_line =  "\t".join(sumstats_header) + "\n" + "\n"
			header_line_pops = "\t".join(pops) + "\n"
			line_mean_pi = "\t".join(mean_pi_output) + "\n"	
#			line_sd_pi = "\t".join(sd_pi_output) + "\n"		
			line_mean_S = "\t".join(mean_S_output) + "\n"
#			line_sd_S = "\t".join(sd_S_output) + "\n"
			line_mean_SPop = "\t".join(mean_SPop_output) + "\n"
#			line_sd_SPop = "\t".join(sd_SPop_output) + "\n"	
			line_mean_SP = "\t".join(mean_SP_output) + "\n"
#			line_sd_SP = "\t".join(sd_SP_output) + "\n"		
#			line_mean_SI = "\t".join(mean_SI_output) + "\n"
#			line_sd_SI = "\t".join(sd_SI_output) + "\n"		
			line_mean_discarded_sites = "\t".join(mean_discarded_sites_output) + "\n"
#			line_sd_discarded_sites = "\t".join(sd_discarded_sites_output) + "\n"			
			line_mean_WattersonTheta = "\t".join(mean_WattersonTheta_output) + "\n"
#			line_sd_WattersonTheta = "\t".join(sd_WattersonTheta_output) + "\n"		
			line_mean_TajimaD = "\t".join(mean_TajimaD_output) + "\n"
#			line_sd_TajimaD = "\t".join(sd_TajimaD_output) + "\n"		
#			line_mean_singletons = "\t".join(mean_singletons_output) + "\n"		
			line_mean_doubletons = "\t".join(mean_doubletons_output) + "\n"		
			line_mean_tripletons = "\t".join(mean_tripletons_output) + "\n"		
			line_mean_SFS = "\t".join(mean_SFS_output) + "\n"	
#			line_mean_Fst_unweighted = "Mean_Fst_unweighted\t" + str(mean_Fst_unweighted) + "\n"
#			line_sd_Fst_unweighted = "SD_Fst_unweighted\t" + str(sd_Fst_unweighted) + "\n"		
			line_mean_Fst_WC = "Mean_Fst_WC\t" + str(mean_Fst_WC) + "\n"		
#			line_sd_Fst_WC = "SD_Fst_WC\t" + str(sd_Fst_WC) + "\n"		
			if no_iterations > 1:
				line_outliers_Fst_WC = "Outliers_Fst_WC\t" + str(outliers_Fst_WC) + "\n"
			line_mean_Dxy = "Mean_Dxy\t" + str(mean_Dxy) + "\n"		
#			line_sd_Dxy = "SD_Dxy\t" + str(sd_Dxy) + "\n"
			
			if no_iterations == 1:
				sumstatsfile.writelines([header_line, header_line_pops, line_mean_pi, line_mean_S, line_mean_SPop, line_mean_SP, line_mean_discarded_sites, line_mean_WattersonTheta, line_mean_TajimaD, line_mean_doubletons, line_mean_tripletons, line_mean_SFS, line_mean_Fst_WC, line_mean_Dxy])
			elif no_iterations > 1: 
				sumstatsfile.writelines([header_line, header_line_pops, line_mean_pi, line_mean_S, line_mean_SPop, line_mean_SP, line_mean_discarded_sites, line_mean_WattersonTheta, line_mean_TajimaD, line_mean_doubletons, line_mean_tripletons, line_mean_SFS, line_mean_Fst_WC, line_outliers_Fst_WC, line_mean_Dxy])
#			sumstatsfile.writelines([header_line, header_line_pops, line_mean_pi, line_sd_pi, line_mean_S, line_sd_S, line_mean_SP, line_sd_SP, line_mean_SI, line_sd_SI, line_mean_discarded_sites, line_sd_discarded_sites, line_mean_WattersonTheta, line_sd_WattersonTheta, line_mean_TajimaD, line_sd_TajimaD, line_mean_singletons, line_mean_doubletons, line_mean_tripletons, line_mean_SFS, line_mean_Fst_unweighted, line_sd_Fst_unweighted, line_mean_Fst_WC, line_sd_Fst_WC, line_outliers_Fst_WC, line_mean_Dxy, line_sd_Dxy])
	
	elif output == "ABC":   # Outputs results in the format required for ABCToolbox, i.e. 2 rows with first row comprising of summary statistics headers and second row comprising of summary statistics values.
		# Define lists for headers and relevant summary statistics
		header_sumstats = ["Pi", "WattersonTheta", "Private_S", "TajimasD"]
#		header_sumstats = ["Mean_pi", "SD_pi", "Mean_WattersonTheta", "SD_WattersonTheta", "Mean_Private_S", "SD_Private_S", "Mean_TajimasD", "SD_TajimasD"]
		results_sumtats = [mean_pi, mean_WattersonTheta, mean_SP, mean_TajimaD]
		header_SFS = ["Doubletons", "Tripletons"]
		if no_iterations == 1:
			header_Fst_Dxy = ["Fst_WC", "Dxy"]
			results_Fst_Dxy = [mean_Fst_WC, mean_Dxy]
		elif no_iterations > 1: 
			header_Fst_Dxy = ["Fst_WC", "Outliers_Fst_WC", "Dxy"]			
			results_Fst_Dxy = [mean_Fst_WC, outliers_Fst_WC, mean_Dxy]
		
		# Make two lists that contain, respectively, all headers and all summary statistic values
		sumstat_header = ["S"]
		sumstat_results = [mean_S["Population1"]]          # Just take one value, since this is the same for all populations (the total number of segregating sites) 
		for header, results in zip(header_sumstats, results_sumtats):
			for pop in range(1,len(pop_size)+1,1):
				pop_header = header + "_P" + str(pop)
				sumstat_header.append(pop_header)
				pop_result = results["Population" + str(pop)]
				sumstat_results.append(pop_result)
		for header, i in zip(header_SFS, range(1,3)):	 # note this calculates until tripletons (modified: singletons now excluded because they are too sensitive to error). This requires that the population has at least 3 samples (individuals) per population, if data is polarised, or 2*3 samples per population, if data is unpolarised. Or generally, n samples (polarised data) or 2n samples (polarised data) per population where n is the variant class frequency desired.
			for pop in range(1,len(pop_size)+1,1):
				SFS_header = header + "_P" + str(pop)
				sumstat_header.append(SFS_header)
				SFS_entry = float(mean_SFS["Population" + str(pop)][i])
				sumstat_results.append(SFS_entry)
		for header, results in zip (header_Fst_Dxy, results_Fst_Dxy):		
			for i in range(1,len(pop_size)+1,1):
				for j in range(i+1,len(pop_size)+1,1):		
					poppair_header = header + "_P" + str(i) + "_" + str(j)
					sumstat_header.append(poppair_header)
					poppair_entry = results["Population" + str(i) + "_" + str(j)]
					sumstat_results.append(poppair_entry)

		# To write sumstats output file
		with open("summary_stats_temp.txt", "w") as sumstatsfile:
			if len(sumstat_header) == len(sumstat_results):
				sumstat_header = "\t".join(sumstat_header) + "\n"
				sumstat_results = '\t'.join(map(str,sumstat_results))					
				sumstatsfile.writelines([sumstat_header, sumstat_results])					
			else:
				error_warning = "Number of headers does not match number of summary statistis, please recheck script."
				sumstatsfile.writelines([error_warning])

	elif output == "ABC-":   # Outputs results in the format required for ABCToolbox, i.e. 2 rows with first row comprising of summary statistics headers and second row comprising of summary statistics values. This format outputs a reduced number of summary statistics, by collapsing P1, P2 and P3 stats into a mean PLow stat, and P4, P5 and P6 stats into a mean PHigh stat.
		# First we reduce the many summary statistics. We do this by collapsing P1, P2 and P3 stats into a mean PLow stat, and P4, P5 and P6 stats into a mean PHigh stat. Thus we calculate the means of the high and low classes.
		mean_pi_list = []
		mean_WattersonTheta_list = []
		mean_SP_list = []
		mean_TajimaD_list = []
		for pop in range(1,len(pop_size)+1,1):
			mean_pi_list.append(mean_pi["Population" + str(pop)])
			mean_WattersonTheta_list.append(mean_WattersonTheta["Population" + str(pop)])
			mean_SP_list.append(mean_SP["Population" + str(pop)])
			mean_TajimaD_list.append(mean_TajimaD["Population" + str(pop)])
		mean_pi_high = numpy.mean(mean_pi_list[:3],axis=0)
		mean_pi_low = numpy.mean(mean_pi_list[3:],axis=0)
		mean_WattersonTheta_high = numpy.mean(mean_WattersonTheta_list[:3],axis=0)
		mean_WattersonTheta_low = numpy.mean(mean_WattersonTheta_list[3:],axis=0)				
		mean_SP_high = numpy.mean(mean_SP_list[:3],axis=0)
		mean_SP_low = numpy.mean(mean_SP_list[3:],axis=0) 
		mean_TajimaD_high = numpy.mean(mean_TajimaD_list[:3],axis=0)
		mean_TajimaD_low = numpy.mean(mean_TajimaD_list[3:],axis=0)
		# We do the same for the pairwise sumstats, namely Fst and Dxy		
		mean_Fst_WC_within_list = [mean_Fst_WC["Population1_2"], mean_Fst_WC["Population1_3"], mean_Fst_WC["Population2_3"], mean_Fst_WC["Population4_5"], mean_Fst_WC["Population4_6"], mean_Fst_WC["Population5_6"]]
		mean_Fst__WC_between_list = [mean_Fst_WC["Population1_4"], mean_Fst_WC["Population1_5"], mean_Fst_WC["Population1_6"], mean_Fst_WC["Population2_4"], mean_Fst_WC["Population2_5"], mean_Fst_WC["Population2_6"], mean_Fst_WC["Population3_4"], mean_Fst_WC["Population3_5"], mean_Fst_WC["Population3_6"]]
		mean_Fst_WC_within = numpy.mean(mean_Fst_WC_within_list,axis=0)
		mean_Fst_WC_between = numpy.mean(mean_Fst__WC_between_list,axis=0)	
		mean_Dxy_within_list = [mean_Dxy["Population1_2"], mean_Dxy["Population1_3"], mean_Dxy["Population2_3"], mean_Dxy["Population4_5"], mean_Dxy["Population4_6"], mean_Dxy["Population5_6"]]
		mean_Dxy_between_list = [mean_Dxy["Population1_4"], mean_Dxy["Population1_5"], mean_Dxy["Population1_6"], mean_Dxy["Population2_4"], mean_Dxy["Population2_5"], mean_Dxy["Population2_6"], mean_Dxy["Population3_4"], mean_Dxy["Population3_5"], mean_Dxy["Population3_6"]]
		mean_Dxy_within = numpy.mean(mean_Dxy_within_list,axis=0)
		mean_Dxy_between = numpy.mean(mean_Dxy_between_list,axis=0)
		# We do the same for the SFS categories
		mean_SFS_high = []
		mean_SFS_low = []
		for i in range(1,3):	 # note this calculates until tripletons (modified: singletons now excluded because they are too sensitive to error). This requires that the population has at least 3 samples (individuals) per population, if data is polarised, or 2*3 samples per population, if data is unpolarised. Or generally, n samples (polarised data) or 2n samples (polarised data) per population where n is the variant class frequency desired.
			mean_SFS_high_list = []
			for pop in range(1,int((len(pop_size)/2)+1),1):
				SFS_entry = float(mean_SFS["Population" + str(pop)][i])
				mean_SFS_high_list.append(SFS_entry)
			mean_SFS_high_entry = numpy.mean(mean_SFS_high_list, axis=0)
			mean_SFS_high.append(mean_SFS_high_entry)
			mean_SFS_low_list = []	
			for pop in range(int((len(pop_size)/2)+1),int(len(pop_size)+1),1):
				SFS_entry = float(mean_SFS["Population" + str(pop)][i])
				mean_SFS_low_list.append(SFS_entry)
			mean_SFS_low_entry = numpy.mean(mean_SFS_low_list, axis=0)
			mean_SFS_low.append(mean_SFS_low_entry)
		mean_doubleton_high = mean_SFS_high[0]
		mean_doubleton_low = mean_SFS_low[0]
		mean_tripleton_high = mean_SFS_high[1]
		mean_tripleton_low = mean_SFS_low[1]

		# Make two lists that contain, respectively, all headers and all summary statistic values 
		# Remember for S, we can take any population value, since this is the same for all populations (the total number of segregating sites)
		sumstat_header = ["S", "Pi_Hi", "Pi_Lo", "ThetaW_Hi", "ThetaW_Lo", "PrivS_Hi", "PrivS_Lo", "TajD_Hi", "TajD_Lo", "Doubletons_Hi", "Doubletons_Lo", "Tripletons_Hi", "Tripletons_Lo", "Fst_within", "Fst_between", "Dxy_within", "Dxy_between"]
		sumstat_results = [mean_S["Population1"], mean_pi_high, mean_pi_low, mean_WattersonTheta_high, mean_WattersonTheta_low, mean_SP_high, mean_SP_low, mean_TajimaD_high, mean_TajimaD_low, mean_doubleton_high, mean_doubleton_low, mean_tripleton_high, mean_tripleton_low, mean_Fst_WC_within, mean_Fst_WC_between, mean_Dxy_within, mean_Dxy_between]

		# To write sumstats output file
		with open("summary_stats_temp.txt", "w") as sumstatsfile:
			if len(sumstat_header) == len(sumstat_results):
				sumstat_header = "\t".join(sumstat_header) + "\n"
				sumstat_results = '\t'.join(map(str,sumstat_results))					
				sumstatsfile.writelines([sumstat_header, sumstat_results])					
			else:
				error_warning = "Number of headers does not match number of summary statistics, please recheck script."
				sumstatsfile.writelines([error_warning])	
			
#	# To write Fst output file (use numpy.loadtxt to easily load and read this output file back into python)
#	header_Fst = "\t".join(pops[1:])
#	numpy.savetxt('fst_output_temp.txt', Fst, fmt='%1.8f', delimiter = "\t", header = header_Fst)   # fmt option here outputs values to 8 decimal figures
	
	# To generate SFS figures in png or pdf format
	if plotSFS == True:
		plot_SFS(mean_SFS, pop_size, file_format = "png")
	if plotPi == True:
		site_pi = site_pi_calc(retained_allele_freqs, no_iterations, pop_size)
		plot_pi (site_pi, pop_size, file_format = "png")
		

#%%