#!/usr/bin/python
# This first line (the shebang line) tells the script to use the interpreter specified by the path

# sumstatscalc version
version = "lsd_high version: 22_03_2019"

### Load dependencies and set working directory/path ###

#from __future__ import division # to allow fractional division with Python 2.X. Not necessary for Python 3.X
import numpy
from scipy.stats import binom
import argparse
import os
import sumstatscalc
from datetime import datetime
#import csv
#import glob

#os.chdir("/Users/luqman/Documents/Hirzi/ETHPHD/Analysis/PoolSimulator")

#%%

######################  LSD-HIGH FUNCTION AND WRAPPER FUNCTION  ######################

# This script contains the main function (lines 168-658), and the wrapper function (lines 659-687) that runs the lsd_high_simulate function and summary statistics calculator function (imported), given a list of supplied arguments (on the command line).

########## NOTES AND EXAMPLE USAGE OF MAIN FUNCTION ##########

### EXAMPLE USAGE (for argument and function options, see lines 535-549):
## Example usage using command line input (simulating pooling process, sequence errors and invariant sites): 
# (single iteration)   python lsd_high.py output_1.txt -d 40 -d 40 -d 40 -d 40 -d 40 -d 40 -l 12500 -p -i --error_method 4 --error_rate 0.001 --minallelecount 1 --mindepth 10 --maxdepth 100 --sampler nbinom -c Dsyl_coverage_nbinomDist.txt -o output_1.sumstats -f ABC
# (100 iterations)     python lsd_high.py output.txt -n 100 -d 40 -d 40 -d 40 -d 40 -d 40 -d 40 -l 12500 -p -i --error_method 4 --error_rate 0.001 --minallelecount 1 --mindepth 10 --maxdepth 100 --sampler nbinom -o output.sumstats -f ABC
## Example usage using command line input (individual data (aka not pooled), no simulation of invariant sites and errors (valid only for high coverage)): 
# (single iteration)   python lsd_high.py output_1.txt -d 40 -d 40 -d 40 -d 40 -d 40 -d 40 -l 12500 --mindepth 10 --maxdepth 100 --sampler nbinom -o output_1_test4.sumstats -f ABC 
## Example usage using separate parameter file input (DEPRECIATED): 
#                      python lsd_high.py dianthus_test1.params

### NOTES
# Note 0: This script was originally written for pooled data, and later extended to work with individual data. As such, certain functions need to be updated to ensure full functionality with individual data.
# Note 1: E.g. currently simulation of error is only implemented for pooled data. For non-pooled  (i.e. individual) data, invariant sites should be set to false (since we simulate invariant errors according to a pooling scheme). Thus, high-coverage is a requirement for using individual data with lsd_high. UPDATE NEEDED TO ALSO SIMULATE PER-INDIVIDUAL SEQUENCING ERRORS!
# Note 2: Number of independent iterations/simulations can be interpreted as number of independent loci.
# Note 3: Here, the minimum allele count is defined as follows: If allele frequency is greater than minimum allele count, then we keep. If equal to or below minimum allele count, we discard. E.g. a minimum allele count = 2 means any alleles with counts of 2 or below are discarded (need a minimum of 3 to keep).
# Note 4: The error rate (currently only implemented for pooled data) in this script defines the rate at which an individual base read is read incorrectly during the random sampling process of the pooling process (random sampling n coverage times with replacement), per site. It is NOT exactly analogous/comparable to genotyping or SNP call error rate (it is much lower).
# Note 5: This error simulation operates in addition to the sampling bias that will result from stochasticity from the random sampling.
# Note 6: msms, ms, cosi2, etc (i.e. most coalescent simulators) work with haploid samples. To use with diploids, simply take 2N individuals.
# Note 7: At present, we are simulating without recombination. We make the assumption that within simulated loci, there is no recombination; hoever between loci, there is free recombination. This is a simplified approach and should work well assuming you use relatively short loci.
# Note 8: Ghost pops (pops of sample size 0) need/should not be included in the above command line/params file! (since no samples are outputted by coalescent simulator). Thus leave them out.

### MISC
## We define and calculate sequence length: based off the equation theta_parameter = 4 * Ne * Mu, where Mu = Mu_site * seq_length (ONLY used for simulations generated using msms)
#def sequence_length (Ne, Mu_site, theta_parameter):
#	const = 1 / (4 * Ne * Mu_site)
#	seq_length = round(const * theta_parameter)
#	return seq_length

#%%

######################  PARSE ARGUMENTS  ######################
# This is the main (master or wrapper) function that runs the lsd_high_simulate function (defined above) and summary statistics calculator function (imported), given a list of supplied arguments (on the command line). Use this to run analyses #

## We parse arguments from the command line using argparse 
parser = argparse.ArgumentParser(description="LSD-high sequence simulator and summary statistics calculator (version 0.1, Luqman 2020)", add_help=True)

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
	raw_no_iterations = int(arg_dict["raw_no_iterations"])
	pop_size = [int(x) for x in arg_dict["pop_size"].replace(" ","").split(",")]
	seq_length = int(arg_dict["seq_length"])
	pooled = bool(arg_dict["pooled"] == "true")
	segsites = arg_dict["segsites"]
	invariant_sites = bool(arg_dict["invariant_sites"] == "true")
	span =  arg_dict["across"]	
	error_method = arg_dict["error_method"]
	error_rate = float(arg_dict["error_rate"])
	min_allele_count = int(arg_dict["min_allele_count"])
	plotSFS = bool(arg_dict["plotSFS"] == "true")
	plotPi = bool(arg_dict["plotPi"] == "true")
	output_filename = arg_dict["output"]
	output_format = arg_dict["output_format"]
"""

# 2. Or we parse arguments directly from the command line:	
#else:
parser.add_argument("filename", action="store", help="The input file, which is a ms format output file.")
parser.add_argument("-n", action="store", dest="nrep", default=1, type=int, help="The number of independent iterations (alternatively, can be interpreted as the number of loci simulated).")
parser.add_argument("-d", action="append", required=True, dest="demesizes", help="The sample sizes of the different populations.")
parser.add_argument("-l", action="store", dest="seqlength", default=0, type=int, help="The length of the simulated region (i.e. sequence length). (required for simulations generated with COSI2).")
parser.add_argument("-p", action="store_true", default=False, dest="pooled", help="To pool sequences.")
parser.add_argument("-i", action="store_true", default=False, dest="invariants", help="To estimate and account for invariants sites.")
parser.add_argument("-a", action="store", dest="across", default="segsites", help="Define whether to calculate Fst and Dxy summary statistics across only segregating sites (segsites) or across all sites (allsites).")
parser.add_argument("--polarised", action="store_true", default=False, dest="polarised", help="Define whether to consider the simulated sequences as polarised (ancestral:0 - derived:1) or unpolarised (major:most frequent allele - minor: 2nd most frequent allele).")
parser.add_argument("--error_method", action="store", dest="error_method", default="4", help="Error calculation method (five options: 0, 1, 2, 3 and 4). Error method 0 doesn't account for error, only for stochastic sampling bias from the pooling process. Error methods 1, 2, 3 and 4 account for stochastic sampling bias from the pooling process in addition to randomly including error at each read per site at (probability) error_rate. Error method 1: If the error leads to a change to alleles 0 or 1, the read with the error is retained and appended to the total counts. If the error leads to inobservable alleles 2 and 3, the reads are discarded (no sites discarded in this method). Error method 2: All reads are retained, regardless of whether the error leads to alleles 0, 1, 2 or 3. If allele 2 and/or 3 is above the minimum allele count, we consider it a triallelic site and discard the site. If they are equal to or below the minimum allele count, and allele 1 is above the minimum allele count, we consider the site biallelic (for alleles 0 and 1) , and retain the original coverage (no reads are discarded in this method). Error method 3: All reads are retained, regardless of whether the error leads to alleles 0, 1, 2 or 3. This method counts the two most frequent alleles (the ancestral + the most frequent minor allele between allele 1, allele 2 and allele 3), given that the 1) most frequent minor allele is > the minumum allele count, and 2) the the most frequent minor allele is NOT equal to the second most frequent minor allele (otherwise site is considered triallelic and discarded). If these two points are satisfied, the two least frequent alleles are ignored and we consider the site biallelic, whilst retaining the original coverage (no reads are discarded here). Error method 4: An approximation of error method 3, which results in approximately 5x shorter runtimes. The approximation assumes that all retained erroneous sites have (min_allele_count + 1) number of errors. Sites with less number of errors are filtered by the minimum allele criterion, sites with more errors are exponentially much less likely, and to the closest approximation, results in negligable number of entries.")
parser.add_argument("--error_rate", action="store", dest="error_rate", default=0.0, type=float, help="Error rate in sampling/calling SNPs during the simulated pooling process.")
parser.add_argument("--minallelecount", action="store", dest="min_allele_count", default=2, type=int, help="The minimum allele count to call a SNP; here a minimum allele count of 2 means that alleles that occur with frequencies 2 or below are discarded. Only required if using Error Methods 2, 3 or 4.")
parser.add_argument("--mindepth", action="store", dest="min_cov", default=0, type=int, help="Minimum depth filter")
parser.add_argument("--maxdepth", action="store", dest="max_cov", default=999, type=int, help="Maximum depth filter")
parser.add_argument("--plotSFS", action="store_true", default=False, dest="plotSFS", help="To output plots of site frequency spectrum (one for each population).")
parser.add_argument("--plotPi", action="store_true", default=False, dest="plotPi", help="To output plots of nucleotide diversity, pi, across site positions (one for each population).")
parser.add_argument("--sampler", action="store", dest="cov_sampler", default="norm", help="Define distribution for which to simulate coverage.")
parser.add_argument("-c", action="store", default=None, dest="obs_cov", help="A file containing per population coverage distribution parameters, to inform coverage simulation.")
parser.add_argument("-o", action="store", default="summary_stats_temp.txt", dest="output", help="Prefix for output files.")
parser.add_argument("-f", action="store", default="full", dest="output_format", help="Format of output file. 'Full' outputs a convenient, human-readable summary output while 'ABC' outputs the file in the format required for ABCToolbox.")
args = parser.parse_args()

##  We associate the parsed arguments with our global variables (defined in the functions below)
filename = args.filename
raw_no_iterations = args.nrep
pop_size = [int(x) for x in args.demesizes]
pooled = args.pooled
invariant_sites = args.invariants
span = args.across
polarised = args.polarised
error_method = args.error_method
error_rate = args.error_rate
min_allele_count = args.min_allele_count
min_cov = args.min_cov
max_cov = args.max_cov
plotSFS = args.plotSFS
plotPi = args.plotPi
cov_sampler = args.cov_sampler
obs_cov = args.obs_cov
output_filename = args.output
output_format = args.output_format
seq_length = args.seqlength

##  Print selected arguments and program version:
print(version)
print("Filename is: " + str(filename))
print("Number of iterations (loci) is: " +str(raw_no_iterations))
print("List of respective population sample sizes are: " + str(pop_size))
print("The length of the simulated sequences are: " + str(seq_length))
print("The sequences are pooled: " + str(pooled))
print("The simulated sequences are polarised: " + str(polarised))
#print("Fst statistics are calculated across " + str(span))
print("Invariant sites are estimated and accounted for: " + str(invariant_sites))
print("Error method: " + str(error_method))
print("The sequencing error rate is: " + str(error_rate))
print("The minimum allele count (only used in error methods 2,3 and 4) is: " + str(min_allele_count))
print("Minimum depth filter: " + str(min_cov))
print("Maximum depth filter: " + str(max_cov))
if cov_sampler == "norm":
	print("Coverage will be simulated according to a normal distribution")
elif cov_sampler == "nbinom":
	print("Coverage will be simulated according to a negative binomial distribution")
if obs_cov == None:
	print("Empirical estimates of coverage distribution parameters not provided. Coverage will be simulated with default parameter settings (mean of 30 and standard deviation or dispersion parameter of 15).")
elif obs_cov != None:	
	print("Empirical estimates of per population coverage distribution parameters is provided and will be used.")
print("Plots of SFS are produced: " + str(plotSFS))
print("Plots of pi are produced: " + str(plotPi))
print("Prefix for output files is: " + str(output_filename))
print("Output format is: " + str(output_format))

#%%

###################### DEFINE SUPPORTING FUNCTIONS ######################

# We define a function to read in user given estimates of coverage (distribution parameters). If a normal distribution is defined, the first and second rows represent mean and sd respectively, if the negative binomial is defined, the first and second rows correspond to mean and dispersal (size) parameters respectively.
# If obs_cov is not defined (given by user), the coverage distribution will be that of a normal or negative binomial distribution with mean 30 and sd/size of 15 for all populations.
def read_obs_cov_file (obs_cov, pop_size, mean = 30, sd = 15):
	if obs_cov != None:
		with open(obs_cov) as f:
			obs_cov_raw = f.read().splitlines()	
		obs_cov_lines = [line.split("\t") for line in obs_cov_raw]
		mean_coverage_observed = [float(x) for x in obs_cov_lines[0]]
		sd_coverage_observed = [float(x) for x in obs_cov_lines[1]]
	elif obs_cov == None:	
		mean_coverage_observed = [mean] * len(pop_size)
		sd_coverage_observed = [sd] * len(pop_size)
	obs_cov_results = [mean_coverage_observed, sd_coverage_observed]
	return obs_cov_results

# We define the coverage sampler, from which samples are drawn from the observed coverage distribution.
# To find the most suitable coverage sampler, first fit the observed data to possible distributions (e.g. negative binomial, binomial, normal, etc). Use Fit_dist.R (which makes use of the fitdistrplus package) for this.
# For this data, negative binomial appears to be the best fit.
if cov_sampler == "norm":
	## To sample coverage from a NORMAL DISTRIBUTION.
	def cov_distr (mean_cov, sd_cov, min_depth = 1):
		coverage = 0
		while coverage < min_depth:
			coverage = round(numpy.random.normal(mean_cov,sd_cov))
		return coverage
elif cov_sampler == "nbinom":	
	## To sample coverage from a NEGATIVE BINOMIAL DISTRIBUTION. Since the Fit_dist.R parametises in terms of paramters mu and size, while numpy parametises in terms of n and p, we first have to convert these parameters accordingly.
	def convert_params(mu, size):
	    """
	    Convert mean/dispersion parameterization of a negative binomial to the ones scipy/numpy supports
	    See https://en.wikipedia.org/wiki/Negative_binomial_distribution#Alternative_formulations
		   https://stat.ethz.ch/R-manual/R-patched/library/stats/html/NegBinomial.html
		   https://stackoverflow.com/questions/40846992/alternative-parametrization-of-the-negative-binomial-in-scipy
		   https://stats.stackexchange.com/questions/260580/negative-binomial-distribution-with-python-scipy-stats
	    """
	    r = size
	    var = mu + 1 / r * mu ** 2
	    p = (var - mu) / var
	    return r, 1 - p	
	def cov_distr (mu, size, min_depth = 1):
		coverage = 0
		while coverage < min_depth:
			coverage = numpy.random.negative_binomial(*convert_params(mu, size))
		return coverage

# Error method 1 randomly includes error at each read per site at error_rate. If the error leads to a change to alleles 0 or 1, the read with the error is retained and appended to the total counts. If the error leads to inobservable alleles 2 and 3, the reads are discarded.
# This function below defines the sampling process for error method 1. The filtering process is included in the main script.
def error_method1 (column_random_sample, error_rate):	
	column_random_sample_with_errors = []	
	for sample in column_random_sample:
		if numpy.random.random() < error_rate:
			call_with_errors = int(((int(sample)) + (numpy.random.choice([1,2,3], size=1, replace=True))) % 4)
			if call_with_errors <= 1:
				column_random_sample_with_errors.append(str(call_with_errors))  # we discard error calls to alleles 2 and 3; which if/when occurs will consequently reduce effective coverage to ((total coverage) - (# of discards))
		else:
			column_random_sample_with_errors.append(sample)
	return column_random_sample_with_errors

# Error method 2 randomly includes error at each read per site at error_rate. Here, all reads are retained, regardless of whether the error leads to alleles 0, 1, 2 or 3. If allele 2 or 3 is above minimum allele count, we consider it a triallelic site and discard the site. If they are equal to or below the minimum allele count, and allele 1 is above the minimum allele count, we consider the site biallelic (alleles 0 and 1), and retain the original coverage (no reads are discarded here).
# This function below defines the sampling process for error methods 2, 3 and 4. The filtering process is included in the main script.
def error_method2_3_4 (column_random_sample, error_rate):					
	column_random_sample_with_errors = []
	for sample in column_random_sample:
		if numpy.random.random() < error_rate:
			call_with_errors = int(((int(sample)) + (numpy.random.choice([1,2,3], size=1, replace=True))) % 4)
			column_random_sample_with_errors.append(str(call_with_errors))   # we're keeping all alleles here
		else:
			column_random_sample_with_errors.append(sample)
	return column_random_sample_with_errors

##  We define the function below to produce an empty (null) sumstats output file, necessary in the case that the msms output file contains no segsites (data). In such a case, only 6 lines will be written. So we condition on this (on line 181)
def create_null_obs_file(pop_size, raw_no_iterations, output):
	# Define lists for headers and relevant summary statistics
	if output == "ABC":
		header_sumstats = ["Pi", "WattersonTheta", "Private_S", "TajimasD"]
		header_SFS = ["Doubletons", "Tripletons"]
		if raw_no_iterations > 1:
			header_Fst_Dxy = ["Fst_WC", "Outliers_Fst_WC", "Dxy"]
		elif raw_no_iterations == 1:
			header_Fst_Dxy = ["Fst_WC", "Dxy"]	
		# Make two lists that contain, respectively, all headers and all summary statistic values
		sumstat_header = ["S"]
		sumstat_results = [0]          # Just take one value, since this is the same for all populations (the total number of segregating sites) 
		for header in header_sumstats:
			for pop in range(1,len(pop_size)+1,1):
				pop_header = header + "_P" + str(pop)
				sumstat_header.append(pop_header)
				pop_result = 0
				sumstat_results.append(pop_result)
		for header, i in zip(header_SFS, range(1,3)):	 # note this calculates until tripletons. This requires that the population has at least 3+1 samples (individuals) per population. Or generally, n+1 samples per population where n is the variant class frequency desired.
			for pop in range(1,len(pop_size)+1,1):
				SFS_header = header + "_P" + str(pop)
				sumstat_header.append(SFS_header)
				SFS_entry = 0
				sumstat_results.append(SFS_entry)
		for header in header_Fst_Dxy:		
			for i in range(1,len(pop_size)+1,1):
				for j in range(i+1,len(pop_size)+1,1):		
					poppair_header = header + "_P" + str(i) + "_" + str(j)
					sumstat_header.append(poppair_header)
					poppair_entry = 0
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
	elif output == "ABC-":
		sumstat_header = ["S", "Pi_Hi", "Pi_Lo", "ThetaW_Hi", "ThetaW_Lo", "PrivS_Hi", "PrivS_Lo", "TajD_Hi", "TajD_Lo", "Doubletons_Hi", "Doubletons_Lo", "Tripletons_Hi", "Tripletons_Lo", "Fst_within", "Fst_between", "Dxy_within", "Dxy_between"]
		sumstat_results = ["0"]*len(sumstat_header)	
		with open("summary_stats_temp.txt", "w") as sumstatsfile:
			sumstat_header = "\t".join(sumstat_header) + "\n"
			sumstat_results = '\t'.join(sumstat_results)					
			sumstatsfile.writelines([sumstat_header, sumstat_results])	

########## DEFINE MAIN FUNCTION ##########

# This function converts the output of ms or msms into individual or pooled sequence format, while incorporating error, and calculates/outputs allele frequencies #

### NOTES ###
## General ##
# Note 0: We allow 5 different methods for error simulation. Error method 0 doesn't account for error, only for stochastic sampling bias from the pooling process. Error methods 1, 2, 3 and 4 account for stochastic sampling bias from the pooling process in addition to randomly including error at each read per site at (probability) error_rate. Error method 1: If the error leads to a change to alleles 0 or 1, the read with the error is retained and appended to the total counts. If the error leads to inobservable alleles 2 and 3, the reads are discarded (no sites discarded in this method). Error method 2: All reads are retained, regardless of whether the error leads to alleles 0, 1, 2 or 3. If allele 2 and/or 3 is above the minimum allele count, we consider it a triallelic site and discard the site. If they are equal to or below the minimum allele count, and allele 1 is above the minimum allele count, we consider the site biallelic (for alleles 0 and 1) , and retain the original coverage (no reads are discarded in this method). Error method 3: All reads are retained, regardless of whether the error leads to alleles 0, 1, 2 or 3. This method counts the two most frequent alleles (the ancestral + the most frequent minor allele between allele 1, allele 2 and allele 3), given that the 1) most frequent minor allele is > the minumum allele count, and 2) the the most frequent minor allele is NOT equal to the second most frequent minor allele (otherwise site is considered triallelic and discarded). If these two points are satisfied, the two least frequent alleles are ignored and we consider the site biallelic, whilst retaining the original coverage (no reads are discarded here). Error method 4: An approximation of error method 3, which results in approximately 5x shorter runtimes. The approximation assumes that all retained erroneous sites have (min_allele_count + 1) number of errors. Sites with less number of errors are filtered by the minimum allele criterion, sites with more errors are exponentially much less likely, and to the closest approximation, results in negligable number of entries.")
# Note 1: When normalising counts, we inevitably end up with fractions, which we round (since we want counts which are whole numbers). Beware that rounding always results in a loss of information, and thus in downstream calculations, always use allele frequencies or the non-normalised values for the calculation.
# Note 2: Number of independent iterations/simulations can be interpreted as number of independent loci.
# Note 3: Here, the minimum allele count is defined as follows: If allele frequency is greater than minimum allele count, than we keep. If equal to or below minimum allele count, we discard. E.g. a minimum allele count = 2 means any alleles with counts of 2 or below are discarded (need a minimum of 3 to keep).
# Note 4: The error rate (currently only implemented for pooled data) in this script defines the rate at which an individual base read is read incorrectly during the random sampling process of the pooling process (random sampling n coverage times with replacement), per site. It is NOT exactly analogous/comparable to genotyping or SNP call error rate (it is much lower).
# Note 5: This error simulation operates in addition to the sampling bias that will result from stochasticity from the random sampling.
# Note 7: msms, ms, cosi2, etc (i.e. most coalescent simulators) work with haploid samples. To use with diploids, simply take 2N individuals.
# Note 8: At present, we are simulating without recombination. We make the assumption that within simulated loci, there is no recombination; hoever between loci, there is free recombination. This is a simplified approach and should work well assuming you use relatively short loci.
# NOTE 9: Approximation assumption of error method 4; setting the number of errors per erroneous site (which is now defined as a site with > min_allele_count number of errors) to min_allele_count + 1; which is the case the large majority of the time. I.e. the approximation assumes that all retained erroneous sites have (min_allele_count + 1) number of errors. Sites with less number of errors are filtered by the minimum allele criterion, sites with more errors are exponentially much less likely, and to the closest approximation, results in negligable number of entries.
# 	     This approximation results in a >4x speed increase in performance.
"""
		# To see this, we can plot minimum number of errors per site vs number of erroneous sites:
		import matplotlib.pyplot as plt
		import numpy
		from scipy.stats import binom
		import seaborn as sns		
		mean_cov = 45
		error_rate = 0.001
		no_invariant_sites = 25000
		min_no_of_errors_persite = [1,2,3,4,5]		
		#ax = sns.barplot(min_allele_count, numpy.random.binomial(no_invariant_sites, (1 - (binom.cdf([x-1 for x in min_allele_count],mean_cov,error_rate)))))
		ax = sns.barplot(min_no_of_errors_persite, [numpy.log10(i) for i in numpy.random.binomial(no_invariant_sites, (1 - (binom.cdf([x-1 for x in min_no_of_errors_persite],mean_cov,error_rate))))])
		ax.set(xlabel='variant frequency', ylabel='log10(# of erroneous sites)')
		#
"""

## Regarding function usage @@
# Note 1: To run with error_method = "2" or "3", the min_allele_count parameter must be defined.
# Note 2: Error method 0 should only be used with invariant sites turned off (doesn't make sense otherwise).
# Note 3: To avoid these notations, make the function more user friendly, include error message if input parameters incorrectly defined (use try and except, and flow control). Or even better, design it in such a way so as to try to prevent it happening at all.
# Note 4: Convert numbers to floats, when dealing with fractions.

## Example usage:
# pooled_allele_freqs = lsd_high_simulate("output3.txt", raw_no_iterations=100, pop_size=[40,40,40,40,40,40], seq_length = 11500, pooled = True, polarised = True, error_method = "4", invariant_sites = True, error_rate = 0.001, min_allele_count = 1, min_cov = 10, max_cov = 100)
# pooled_allele_freqs = lsd_high_simulate("output2.txt", raw_no_iterations=1, pop_size=[40,40,40,40,40,40], seq_length = 11500, pooled = True, polarised = True, error_method = "4", invariant_sites = True, error_rate = 0.001, min_allele_count = 1, min_cov = 10, max_cov = 100)
# pooled_allele_freqs = lsd_high_simulate("test_short_n1.txt", raw_no_iterations=1, pop_size=[4,6], seq_length = 12500, pooled = True, polarised = True, error_method = "4", invariant_sites = True, error_rate = 0.001, min_allele_count = 1, min_cov = 10, max_cov = 100)
# pooled_allele_freqs = lsd_high_simulate("test_medium.txt", raw_no_iterations=5, pop_size=[4,5,3], seq_length = 12500, pooled = True, polarised = True, error_method = "4", invariant_sites = True, error_rate = 0.001, min_allele_count = 1, min_cov = 10, max_cov = 100)
# pooled_allele_freqs = lsd_high_simulate("test_short2.txt", raw_no_iterations=5, pop_size=[4,6], seq_length = 12500, pooled = True, polarised = True, error_method = "4", invariant_sites = True, error_rate = 0.001, min_allele_count = 1, min_cov = 10, max_cov = 100)
# pooled_allele_freqs = lsd_high_simulate("output_1.txt", raw_no_iterations=1, pop_size=[40,40,40,40,40,40], seq_length = 12500, pooled = True, polarised = True, error_method = "4", invariant_sites = True, error_rate = 0.001, min_allele_count = 2, min_cov = 20, max_cov = 80)

## To explore results: pooled_allele_freqs["Iteration_X"]["PopulationY"]; where X goes from 1 to n number of iterations, and Y from 1 to p number of population.

def lsd_high_simulate (filename, raw_no_iterations, pop_size, obs_cov, seq_length, pooled = False, polarised = False, error_method = "4", invariant_sites = False, error_rate = 0, start_line=6, lines_between_endblock_startblock = 4, min_allele_count = None, min_cov = None, max_cov = None):   ## Let's perhaps  make the function read the parameter info directly from the msms output file?

	### Open the file, define structure of (iteration) blocks
	lines_between_startblock_startblock = sum(pop_size) + lines_between_endblock_startblock
	
	with open(filename, encoding = "utf-8") as f:
		msms_raw = f.read().splitlines()
	
		### We first need to remove fixed loci (because 1. they are uninformative, and 2. they mess up (make irregular) the format (number of lines between blocks) of the ms output file
		# First, we identify the lines which correspond to fixed loci
		fixed_loci = 'segsites: 0'    # We look for fixed loci	
		fixed_loci_indices_list = []	
		for num, line in enumerate(msms_raw, 1):
			if fixed_loci in line:
				fixed_loci_indices_list.extend([num-3, num-2, num-1])	
		fixed_loci_indices = list(set(fixed_loci_indices_list))
		# We then delete the lines corresponding to the fixed loci, making sure to maintain the correct interblock distance (so here we need to remove 3 lines), and to delete in reverse order (to maintain indexing)
		if len(fixed_loci_indices) > 0:
			for i in sorted(fixed_loci_indices, reverse=True):
				del msms_raw[i]
		
		# And then we have to redefine the new (trimmed) number of iterations:
		no_iterations = int(raw_no_iterations - (len(fixed_loci_indices_list)/3))

		# We define empty dicts, where the results (e.g. allele frequencies) will be written to.
		results = {}
		allelefreqs_master_array = {}
		index_discarded_sites_master_array = {}
		
		# We define the coverage distribution from which to sample from.
		observed_coverages = read_obs_cov_file(obs_cov = obs_cov, pop_size = pop_size)
		mean_coverages = observed_coverages[0]
		sd_coverages = observed_coverages[1]
		
		### This first loop iterates over the iteration blocks
#		for count, (iteration, loci_coverages) in enumerate(zip(range(0,(no_iterations*lines_between_startblock_startblock),lines_between_startblock_startblock), coverages)):	 # Should you later wish to allow for different coverage values for different iteratations (or loci). Will require redefining the 'coverages' below. Not sure how reasonable this action is...	
		for count,iteration in enumerate(range(0,(no_iterations*lines_between_startblock_startblock),lines_between_startblock_startblock)):
			
			# We note down the raw number of segregating sites (from the msms output)
			range_segsites = list(range(int(msms_raw[iteration + start_line - 2] [10:])))
			
			
			########## SIMULATING ERROR IN INVARIANT SITES ##########
				
			if invariant_sites == True and seq_length > len(range_segsites):
				
				# We then calculate counts from erroneously called invariant sites (called as polymorphic sites).			
				no_invariant_sites = seq_length - len(range_segsites)

				no_erroneous_polymorphic_sites_pop_array = []				
				for pop in range(len(pop_size)):		
					# We define the population specific mean and standard deviation of the coverage distribution from which the coverage values will be drawn (based on observed coverage distribution). We also define the sample size of the population.
					mean_cov = mean_coverages[pop]
					sd_cov = sd_coverages[pop]

					# We estimate the population specific number of erroneous polymorphic sites, by first calculating the probability of getting at least min_allele_count errors in a site (of cov coverage; hence trials = coverage, prob = error_rate), and then using this probabilty to draw from a binomial distribution (trials = number of invariant sites).
					prob_atleast_min_allele_count_error_in_coverage_trials_per_site = 1 - (binom.cdf(min_allele_count,round(mean_cov),error_rate))         # Remember: in binomial disributions, "at least" n is the complement of "at most" n-1, the latter which is calculated via the cumulative distribution function (cdf). binom.cdf(1,mean_cov,error_rate) is the probabilty of at least one error in n=coverage trials per site (Binomial distribution). See http://mathbits.com/MathBits/TISection/Statistics2/binomialAtMost.htm for calculation
					no_erroneous_polymorphic_sites = numpy.random.binomial(no_invariant_sites, prob_atleast_min_allele_count_error_in_coverage_trials_per_site)				
					no_erroneous_polymorphic_sites_pop_array.append(no_erroneous_polymorphic_sites)
				
				# We then produce a population array which lists the error alleles of the erroneous polymorphic sites.
				list_of_error_samples_pop_array = []
				for pop, n in enumerate(no_erroneous_polymorphic_sites_pop_array):
					list_of_error_samples = []
					for erroneous_site in range(n):
						if error_method == "1" or error_method == "2" or error_method == "3":
							mean_cov = mean_coverages[pop]
							sd_cov = sd_coverages[pop]							
							coverage_erroneous_site = cov_distr (mean_cov, sd_cov)	
							no_of_errors = 0
							while no_of_errors <= min_allele_count:  
								no_of_errors = numpy.random.binomial(coverage_erroneous_site,error_rate)
							list_of_error_samples_persite = (numpy.random.choice(["1","2","3"],no_of_errors,replace=True)).tolist()  # The large majority (9X%) of sites will have only 1 error (check distribution) - these will be removed by the minimum allele criterion. The remaining erroneous sites will mostly be of length 2, and this drops exponentially (what power?) to 3 and so on.
						# See NOTE 1 for approximation and performance increase for method below (approx 4-5x faster)
						# Error method 4 approximates the step in the lines above, by setting the number of errors per erroneous site (which is now defined as a site with > min_allele_count number of errors) to min_allele_count + 1; which is the case the large majority of the time. I.e. the approximation assumes that all retained erroneous sites have (min_allele_count + 1) number of errors. Sites with less number of errors are filtered by the minimum allele criterion, sites with more errors are exponentially much less likely, and to the closest approximation, results in negligable number of entries.
						elif error_method == "4":
							list_of_error_samples_persite = (numpy.random.choice(["1","2","3"],min_allele_count+1,replace=True)).tolist()  # The large majority (9X%) of sites will have only 1 error (check distribution) - these will be removed by the minimum allele criterion. The remaining erroneous sites will mostly be of length 2, and this drops exponentially (what power?) to 3 and so on.
						list_of_error_samples.append(list_of_error_samples_persite) 
					list_of_error_samples_pop_array.append(list_of_error_samples)	
				
				# We then produce a list of indexes for which all the erroneous polymorphic sites will lie on. Here, we assume that the errors are randomly distributed along the length of invariant sites. We capture this by giving each erroneous polymorphic site a randomly drawn index (numpy.random.choice(no_invariant_sites, len(error_samples), replace = False)). We do this because we want to take into acount the possibility that erroneous polymorphic sites can, by chance, land on the same positional index (i.e. site) across populations
				error_samples_pop_dict={}
				indexes_pop_dict = {}
				all_indexes = []
				for pop in range(len(pop_size)):
					error_samples = list_of_error_samples_pop_array[pop]
					indexes = sorted(numpy.random.choice(no_invariant_sites, len(error_samples), replace = False).tolist()) # True or false for replacement?
					all_indexes.append(indexes)	
					error_samples_pop_dict["Population"+str(pop+1)] = error_samples
					indexes_pop_dict ["Population"+str(pop+1)] = indexes						
				all_indexes = sorted(list(set([item for sublist in all_indexes for item in sublist])))
					
				# We then produce the population specific lists that detail the error samples per site, positioned along the previously defined all_indexes
				invariant_samples_pop_dict = []
				for pop in range(len(pop_size)):
					indexes = indexes_pop_dict["Population"+str(pop+1)]
					error_samples = error_samples_pop_dict["Population"+str(pop+1)]
					# We note down the sites (indexes) that are not erroneously called in this population; we will convert these to zeroes later. We don't do it now because converting to zeroes now may cause errors should the 0 positional index be found in the population indexes.
					index_of_zero_indexes = []
					indexes_zero_sites = sorted([x for x in all_indexes if x not in indexes])
					for ind in indexes_zero_sites:
						if ind in all_indexes:
							index_of_zero_indexes.append(all_indexes.index(ind))	
					# We then convert the sites in the list which are erroneously called in this population to the error samples for the site, by index	
					error_samples_pop = [x for x in all_indexes]
					for ind, error_sample in zip(indexes, error_samples):
						if ind in all_indexes:
							error_samples_pop[error_samples_pop.index(ind)] = error_sample
					# Finally, we convert the rest (the non-erroneously called sites) in this population to zero
					for ind in index_of_zero_indexes:
						error_samples_pop[ind] = ["0"]
					invariant_samples_pop_dict.append(error_samples_pop)
				invariant_samples_pop_array = numpy.asarray(invariant_samples_pop_dict)
					
				# We convert the list of error samples to list of error counts. We produce a list of indexes for discarded sites; these are sites that will be discarded due to not satisfying the filtering requirements (which is identical to the pooling error filters)						
				index_discarded_sites_invariant_errors_pop_array = []
				erroneous_counts_pop_dict = {}
				effective_coverages_invariants_pop_dict = {}
				for pop in range(len(pop_size)):
					erroneous_counts = []
					index_discarded_sites_invariant_errors = []						
					effective_coverages_invariants = []
					mean_cov = mean_coverages[pop]
					sd_cov = sd_coverages[pop]											
					for site in range(len(all_indexes)):	
						# We estimate the population specific coverage per site
						coverage_erroneous_site = cov_distr (mean_cov, sd_cov)		
						list_of_error_samples_persite = (invariant_samples_pop_array[pop][site])   
						if list_of_error_samples_persite == ['0']:
							erroneous_counts.append(0) 
							effective_coverages_invariants.append(coverage_erroneous_site) 
						else:
							erroneous_column_count1 = list_of_error_samples_persite.count("1")
							erroneous_column_count2 = list_of_error_samples_persite.count("2")
							erroneous_column_count3 = list_of_error_samples_persite.count("3")
							# Same for polarised and unpolarised options, since in these error cases, the 0 state will always be the major allele
							# Error method 1 randomly includes error at each read per site at error_rate. If the error leads to a change to alleles 0 or 1, the read with the error is retained and appended to the total counts. If the error leads to inobservable alleles 2 and 3, the reads are discarded.	
							if error_method == "1":
								erroneous_counts.append(erroneous_column_count1) 
								discarded_coverage = erroneous_column_count3 + erroneous_column_count2  # effective coverage is reduced as alleles 2 and 3 are discarded	
								effective_coverages_invariants.append((coverage_erroneous_site - discarded_coverage)) 
							# Error method 2 randomly includes error at each read per site at error_rate. Here, all reads are retained, regardless of whether the error leads to alleles 0, 1, 2 or 3. If allele 2 or 3 is above minimum allele count, we consider it a triallelic site and discard the site. If they are equal to or below the minimum allele count, and allele 1 is above the minimum allele count, we consider the site biallelic (alleles 0 and 1), and retain the original coverage (no reads are discarded here).
							elif error_method == "2": 
								if (erroneous_column_count2 > min_allele_count) or (erroneous_column_count3 > min_allele_count) or (0 < erroneous_column_count1 <= min_allele_count):  # if pooled_column_count2 or 3 is > minimum allele count, we discard site.		
									index_discarded_sites_invariant_errors.append(site)
								erroneous_counts.append(erroneous_column_count1)
								effective_coverages_invariants.append(coverage_erroneous_site) 
							# Error methods 3 and 4 count the two most frequent alleles (the ancestral + the most common minor allele between allele 1, allele 2 and allele 3), given that the 1) most frequent minor allele is > the minumum allele count, and 2) the the most frequent minor allele is NOT equal to the second most frequent minor allele (i.e. otherwise we consider the site is triallelic and discard site). If these 2 points are satisfied, we  ignore the two least frequent alleles and consider the site biallelic, retaining the original coverage (no reads are discarded here). Here, all reads are retained, regardless of whether the error leads to alleles 0, 1, 2 or 3.
							elif error_method == "3" or error_method == "4": 
								erroneous_column_count_sorted = sorted([erroneous_column_count1, erroneous_column_count2, erroneous_column_count3], reverse=True)
								if (erroneous_column_count_sorted[0] > 0 and erroneous_column_count_sorted[0] < coverage_erroneous_site and erroneous_column_count_sorted[0] == erroneous_column_count_sorted[1]) or (0 < erroneous_column_count_sorted[0] <= min_allele_count):   # Here, we are retaining the the 2 most frequent alleles (the ancestral allele 0 by default + the most frequent out of the 3 derived alleles). If there are >2 most frequent alleles (e.g. the are two second most frequent alleles), we consider the site triallelic and discard site. If most frequent minor allele is > minimum allele count, we keep. Succinctly, if minor allele is > min allele count and is not equal in frequency to 2nd/3rd most frequent minor alleles, we consider the site to be biallelic and append.                                                                                           
									index_discarded_sites_invariant_errors.append(site)
								erroneous_counts.append(erroneous_column_count_sorted[0])	
								effective_coverages_invariants.append(coverage_erroneous_site)   # remember, in error_method 3 and 4, we do not discard reads, even if erroneously called, as we want to avoid biases in allele frequency that are contingent upon coverage							
					erroneous_counts_pop_dict["Population"+str(pop+1)] = erroneous_counts
					effective_coverages_invariants_pop_dict["Population"+str(pop+1)] = effective_coverages_invariants
					index_discarded_sites_invariant_errors_pop_array.append(index_discarded_sites_invariant_errors)
				index_discarded_sites_invariant_errors_set = sorted(list(set([item for sublist in index_discarded_sites_invariant_errors_pop_array for item in sublist])))

				# Here we discard the sites that failed to satisfy the filtering criteria.
				retained_effective_coverages_invariants_pop_dict = {}
				retained_erroneous_invariants_counts_pop_dict = {}
				for pop in range(len(pop_size)):
					retained_erroneous_counts = []
					retained_effective_coverages = []
					for index, (site, cov) in enumerate(zip(erroneous_counts_pop_dict["Population"+str(pop+1)], effective_coverages_invariants_pop_dict["Population"+str(pop+1)])):
						if index not in index_discarded_sites_invariant_errors_set:
							retained_erroneous_counts.append(site)
							retained_effective_coverages.append(cov)								
					retained_erroneous_invariants_counts_pop_dict["Population" + str(pop+1)] = retained_erroneous_counts	
					retained_effective_coverages_invariants_pop_dict["Population" + str(pop+1)] = retained_effective_coverages
			
			
			########## SIMULATING ERROR AND POOLING IN SEG SITES ##########
			
			# We make an index (list) of the start positions (line number) of the population blocks within each iteration block in the file, including the final end position. We need this to iterate over the populations within each iteration block.  
			index = []
			for j in range(len(pop_size)+1): # +1 since we also need the position of the final end position of the last population
				index_element = start_line + iteration + sum(pop_size[:j])
				index.append(index_element)
		
			# For unpolarised data, we need to determine the minor alleles (meta-population wide)	
			if polarised == False:
				msmsraw_block = []	
				for line in msms_raw[index[0]:index[-1]:1]:
					msmsraw_block.append(list(line))
				all_pops_block = numpy.asarray(msmsraw_block)					
				minor_alleles = []
				for site in range_segsites:    # In this for loop, we're iterating base by base, i.e. column by column
					column_nested_list = (all_pops_block[:,[site]]).tolist()    # to make a list (of list) of the column values
					column_list = [item for sublist in column_nested_list for item in sublist]    # to collapse this list of lists
					minor_allele = min(set(column_list), key=column_list.count) # to return the minor allele. If frequency of both alleles are equal, this returns the lower value (in this case "0")				
					minor_alleles.append(minor_allele)
		
			# We define empty dicts where the results from each iteration will be written to, as well as the number of segregating sites per iteration.
			pop_array = {}	
			allele_freqs = {}
			index_discarded_sites_pop_array = []
					
			# We then loop over the 'index' defined above to generate a character list array for each population block.	
			for pop in range(len(pop_size)):
				list_msmsraw = []	
				for line in msms_raw[index[pop]:index[pop+1]:1]:
					char_list=list(line)
					list_msmsraw.append(char_list)

				# We add the character list array to the previously defined open dict. We then define the population specific mean and standard deviation of the coverage distribution from which the coverage values will be drawn (based on observed coverage distribution). We also define the sample size of the population.
				pop_array["Population" + str(pop+1)] = numpy.asarray(list_msmsraw)
				mean_cov = mean_coverages[pop]
				sd_cov = sd_coverages[pop]
				pop_no_samples = pop_size[pop]

				# We define empty lists where the (temporary) population results will be written to. These results will be appended to the dicts defined above at the end, according to population number, and then iteration number.		
				rawcounts = []
				pooledcounts = []
				effective_coverages = []
				index_discarded_sites = []

				## This is the loop that converts the data to pooled data, as well as perform the allele frequency counts
				for site in range_segsites:    # In this for loop, we're iterating base by base, i.e. column by column
					column_nested_list = (pop_array["Population" + str(pop+1)][:,[site]]).tolist()    # to make a list (of list) of the column values
					column_list = [item for sublist in column_nested_list for item in sublist]    # to collapse this list of lists
					column_count = column_list.count("1")    # counts number of derived alleles at one site (per column) across all samples of a population
					rawcounts.append(column_count)

					coverage = cov_distr (mean_cov, sd_cov)   # We define the coverage for this site, which will be used below
	
					# This is the code to produce a pooled sample (count)					
					if pooled == True:     
						column_random_sample = numpy.random.choice(column_list,coverage,replace=True) # random sample from site (column) 'coverage' times

						## These lines then add / simulate errors in the pooled samples ###
						# Error method 0 doesn't simulate error. 
						if error_method == "0":     							
							if polarised == False:
								pooled_column_count_minor_derived = list(column_random_sample).count(minor_allele)   # to count the 1 sites for the pooled data
							elif polarised == True:
								pooled_column_count_minor_derived = list(column_random_sample).count("1")   # to count the derived sites for the pooled data		
							pooledcounts.append(pooled_column_count_minor_derived)
							effective_coverages.append(coverage)	
						
						# Error method 1 randomly includes error at each read per site at error_rate. If the error leads to a change to alleles 0 or 1, the read with the error is retained and appended to the total counts. If the error leads to inobservable alleles 2 and 3, the reads are discarded.
						elif error_method == "1": 
							column_random_sample_with_errors = error_method1 (column_random_sample, error_rate)
							reduced_coverage = len(column_random_sample_with_errors)   # as coverage is lost when we discard calls to alleles 2 and 3
							if polarised == False:						
								pooled_column_count_minor_derived = column_random_sample_with_errors.count(minor_allele)	# to count the 1 sites for the pooled data
							elif polarised == True:
								pooled_column_count_minor_derived = column_random_sample_with_errors.count("1")   # to count the derived sites for the pooled data	
							pooledcounts.append(pooled_column_count_minor_derived)
							effective_coverages.append(reduced_coverage)	
	
						# Error method 2 randomly includes error at each read per site at error_rate. Here, all reads are retained, regardless of whether the error leads to alleles 0, 1, 2 or 3. If allele 2 or 3 is above minimum allele count, we consider it a triallelic site and discard the site. If they are equal to or below the minimum allele count, and allele 1 is above the minimum allele count, we consider the site biallelic (alleles 0 and 1), and retain the original coverage (no reads are discarded here).
						elif error_method == "2":
							column_random_sample_with_errors = error_method2_3_4 (column_random_sample, error_rate)
							pooled_column_count2 = column_random_sample_with_errors.count("2")   # to count the sites (allele 2) for the pooled data						
							pooled_column_count3 = column_random_sample_with_errors.count("3")   # to count the sites (allele 3) for the pooled data												
							if polarised == False:				
								pooled_column_count_minor_derived = column_random_sample_with_errors.count(minor_allele)   # to count the 1 sites for the pooled data
								pooled_column_count_sorted = sorted([pooled_column_count_minor_derived, pooled_column_count2, pooled_column_count3], reverse=True)
								if (pooled_column_count_sorted[1] > min_allele_count) or (pooled_column_count_sorted[2] > min_allele_count) or (0 < pooled_column_count_sorted[0] <= min_allele_count):  # if pooled_column_count2 or 3 is > minimum allele count, we discard site.		
									index_discarded_sites.append(site)
								pooledcounts.append(pooled_column_count_sorted[0])														
							elif polarised == True:	
								pooled_column_count_minor_derived = column_random_sample_with_errors.count("1")   # to count the derived sites (allele 1) for the pooled data					
								if (pooled_column_count2 > min_allele_count) or (pooled_column_count3 > min_allele_count) or (0 < pooled_column_count_minor_derived <= min_allele_count):  # if pooled_column_count2 or 3 is > minimum allele count, we discard site.		
									index_discarded_sites.append(site)
								pooledcounts.append(pooled_column_count_minor_derived)														
							effective_coverages.append(coverage)
								
						# Error method 3 and 4 randomly includes error at each read per site at error_rate. Here, all reads are retained, regardless of whether the error leads to alleles 0, 1, 2 or 3. Error methods 3 and 4 count the two most frequent alleles (the ancestral + the most common minor allele between allele 1, allele 2 and allele 3), given that the 1) most frequent minor allele is > the minumum allele count, and 2) the the most frequent minor allele is NOT equal to the second most frequent minor allele (i.e. otherwise we consider the site is triallelic and discard site). If these 2 points are satisfied, we  ignore the two least frequent alleles and consider the site biallelic, retaining the original coverage (no reads are discarded here).						
						elif error_method == "3" or error_method == "4":     
							column_random_sample_with_errors = error_method2_3_4 (column_random_sample, error_rate)
							pooled_column_count2 = column_random_sample_with_errors.count("2")   # to count the sites (allele 2) for the pooled data						
							pooled_column_count3 = column_random_sample_with_errors.count("3")   # to count the sites (allele 3) for the pooled data		
							if polarised == False:
								pooled_column_count_minor_derived = column_random_sample_with_errors.count(minor_allele)   # to count the 1 sites for the pooled data	
							elif polarised == True:
								pooled_column_count_minor_derived = column_random_sample_with_errors.count("1")   # to count the derived sites (allele 1) for the pooled data				
							pooled_column_count_sorted = sorted([pooled_column_count_minor_derived, pooled_column_count2, pooled_column_count3], reverse=True)	
							if (pooled_column_count_sorted[0] > 0 and pooled_column_count_sorted[0] < coverage and pooled_column_count_sorted[0] == pooled_column_count_sorted[1]) or (0 < pooled_column_count_sorted[0] <= min_allele_count):   # Here, we are retaining the the 2 most frequent alleles (the ancestral allele 0 by default + the most frequent out of the 3 derived alleles). If there are >2 most frequent alleles (e.g. the are two second most frequent alleles), we consider the site triallelic and discard site. If most frequent minor allele is > minimum allele count, we keep. Succinctly, if minor allele is > min allele count and is not equal in frequency to 2nd/3rd most frequent minor alleles, we consider the site to be biallelic and append.                                                                                           
								#if pooled_column_count_sorted[1] <= min_allele_count and pooled_column_count_sorted[2] <= min_allele_count:    # These shaded out lines add the additional rule that if the 2nd and 3rd most frequent minor alleles is > minimum allele count, we discard site.
								index_discarded_sites.append(site)								
							pooledcounts.append(pooled_column_count_sorted[0])														
							effective_coverages.append(coverage)	

				# We calculate the pooled counts and effective coverages for the case where invariant sites are considered.
				if invariant_sites == True and seq_length > len(range_segsites):
					pooledcounts_with_falsepositives = pooledcounts + retained_erroneous_invariants_counts_pop_dict["Population"+str(pop+1)]
					effective_coverages_with_falsepositives = effective_coverages + retained_effective_coverages_invariants_pop_dict["Population"+str(pop+1)]		

				# We append the single population results to the overall population dicts and calculate the derived allele frequencies		
				if pooled == True:
					if invariant_sites == True and seq_length > len(range_segsites):
						# Add sites which don't satisfy min/max depth criteria to discarded sites
						index_discarded_sites_failed_coverage = [i for i in range(len(effective_coverages_with_falsepositives)) if effective_coverages_with_falsepositives[i] < min_cov or effective_coverages_with_falsepositives[i] > max_cov]
						index_discarded_sites_concat = sorted(list(set(index_discarded_sites + index_discarded_sites_failed_coverage)))
						index_discarded_sites_pop_array.append(index_discarded_sites_concat)
						index_discarded_sites_set = sorted(list(set([item for sublist in index_discarded_sites_pop_array for item in sublist])))
						allele_freqs["Population" + str(pop+1)] = [x/y for x,y in zip(pooledcounts_with_falsepositives, effective_coverages_with_falsepositives)] # this is equivalent to [x / pop_no_samples for x in normalised_pooledcounts], apart from the rounding performed in the normalised_pooledcounts, hence we use this (non-rounded exact values)
					else:
						# Add sites which don't satisfy min/max depth criteria to discarded sites
						index_discarded_sites_failed_coverage = [i for i in range(len(effective_coverages)) if effective_coverages[i] < min_cov or effective_coverages[i] > max_cov]
						index_discarded_sites_concat = sorted(list(set(index_discarded_sites + index_discarded_sites_failed_coverage)))
						index_discarded_sites_pop_array.append(index_discarded_sites_concat)
						index_discarded_sites_set = sorted(list(set([item for sublist in index_discarded_sites_pop_array for item in sublist])))							
						allele_freqs["Population" + str(pop+1)] = [x/y for x,y in zip(pooledcounts, effective_coverages)] # this is equivalent to [x / pop_no_samples for x in normalised_pooledcounts], apart from the rounding performed in the normalised_pooledcounts, hence we use this (non-rounded exact values)
				else:			
					allele_freqs["Population" + str(pop+1)] = [x / pop_no_samples for x in rawcounts]
					index_discarded_sites_set = []										

			# We append the population dicts to the master (iteration) dicts and return the output of the function				
			allelefreqs_master_array["Iteration_" + str(count+1)] = allele_freqs
			index_discarded_sites_master_array["Iteration_" + str(count+1)] = index_discarded_sites_set
			results["Allele_freqs"] = allelefreqs_master_array
			results["index_discarded_sites"] = index_discarded_sites_master_array
			
		return results

#%%

###################### (WRAPPER FUNCTION AND EXECUTION) LSD-HIGH AND SUMSTATS CALCULATOR ######################

##  Finally, we run the lsd_high_simulate and sumstats_calc functions:
if __name__ == "__main__":

	startTime = datetime.now()	
	
	file_no_lines = sum(1 for line in open(filename))

	if file_no_lines > 7:
		raw_allele_freqs = lsd_high_simulate(filename = filename, raw_no_iterations = raw_no_iterations, pop_size = pop_size, obs_cov = obs_cov, seq_length = seq_length, pooled = pooled, polarised = polarised, error_method = error_method, invariant_sites = invariant_sites, error_rate = error_rate, min_allele_count = min_allele_count, min_cov = min_cov, max_cov = max_cov)
		
		sumstatscalc.sumstats_calc (raw_allele_freqs, pop_size = pop_size, seq_length = seq_length, span = span, polarised = polarised, output = output_format, plotSFS = plotSFS, plotPi = plotPi)
	
		os.rename("summary_stats_temp.txt", output_filename)	
		
	else:
		create_null_obs_file(pop_size, raw_no_iterations, output = output_format)
		
		os.rename("summary_stats_temp.txt", output_filename)			
		
	endTime = datetime.now()
	total_runtime = "Total runtime:" + str(endTime - startTime)

	print(total_runtime)

# End of function	
#%%