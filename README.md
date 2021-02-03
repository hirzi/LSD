# Identifying Loci under Selection via explicit Demographic models (LSD)


  This repository contains a suite of scripts for performing LSD genome scans based on explicit demographic models (Luqman et al. 2021, https://www.biorxiv.org/content/10.1101/2020.07.20.211581v2.full.pdf). The current implementation estimates demographic parameters via an Approximate Bayesian Computation (ABC) framework, and works in two steps. First, neutral demographic parameters are estimated (see requirements for LSD point (ii)). Second, per-locus parameter estimates (e.g. for a sliding window across the chromosome) are compared to the neutral estimates, to identify selected loci. 

  As LSD is an ABC approach, it relies on simulations to estimate the posterior distribution of model parameters. The current implementation takes ms-format coalescent samples as input for simulated data and mpileup format (e.g. from BAM files) as input for observed data. A large range of modern coalescent simulators (or those that approximate the coalescent) output ms-format data including e.g. ms (Hudson, 2002), msHOT (Hellenthal & Stephens, 2007), msms (Ewing & Hermisson, 2010), msprime (Kelleher & Etheridge, 2015), MaCS (Chen, Marjoram, & Wall, 2009), cosi2 (Shlyakhter et al., 2014) and SCRM (Staab et al., 2015). 

  The processing, format and final output of observed genetic data will often differ from that of raw coalescent simulations, given that observed genetic data may be subject to various pre-sequencing (e.g. pooling), sequencing (e.g. sequencing errors, stochastic sampling of reads) and post-sequencing (e.g. filters) events that perturb and reformat the data from the original source. We thus provide programs that interface with coalescent simulators to replicate observed sequencing pipelines and generate simulated sequencing data. LSD-High can accommodate and simulate both individual and pooled data and assumes mid to high coverage (>10x) data, while LSD-Low accepts individual data and can additionally accommodate low coverage (>2x) data by utilising genotype likelihoods via msToGLF and ANGSD (Korneliussen, Albrechtsen, & Nielsen, 2014). These programs then calculate a suite of summary statistics for the simulated and observed data. Summary statistics currently implemented include the number of segregating sites (S), private S, nucleotide diversity (pi), Watterson’s theta estimator, Tajima’s D, relative divergence (FST), absolute divergence (DXY), and site frequencies, though in principle any summary statistic can be included with appropriate additions or modifications to the programs’ scripts. 

  ABC is currently implemented via ABCtoolbox (Wegmann et al., 2010). LSD-Scan takes the output of the ABC parameter estimates and estimates the departure of the inferred posteriors from neutral expectations. If conditioning the detection of selected (and linked) loci on multiple (joint) parameters, LSD also outputs the directionality of the deviation in the joint posterior with respect to the marginal parameters, and represents these as colours in the genome scan Manhattan plot.

  We note that LSD is NOT a program, rather it is an analytical framework for identifying loci under selection via deviations in demographic parameters. As such, it is not constrained to any particular program. Rather, we envision a custom and modular implementation that may interface with any appropriate combination of coalescent simulator, summary statistics calculator and ABC program. Here, we simply propose an example implementation utilising msms, LSD-High/Low and ABCtoolbox. As currently implemented, LSD takes whole genome sequence (WGS) data and conditions the inference of selection on genomic windows (regions).

<img src="https://github.com/hirzi/LSD/blob/master/ABC_demography_pipeline_general_landscape3.png" width="800"> 

# Requirements for LSD

  Firstly, a demographic model needs to be defined. Definition and choice of the demographic model should i) be informed by knowledge of the study system, ii) be motivated by the model’s capacity to provide a useful approximation of a biological process of interest, and iii) be sufficiently simple to remain computationally tractable. Additionally, given that we condition the inference of selection on demographic parameters, the model should be formulated according to whether deviation in N<sub>E</sub> or M<sub>E</sub> is desired for the inference of selection. Finally, the model should be able to sufficiently describe the neutral genetic variation of the system. This can be validated by demonstrating that the observed data can be accurately and sufficiently captured by the simulated data (e.g. in terms of summary statistics)

  To estimate neutral demographic parameters in the first step, we need an a priori set of regions that we believe reflect neutral evolution. Such a set may be informed by the particular structural or functional class the sites belong to and may e.g. consist of genomic regions not linked to structural annotations. Alternatively, given that LSD is very robust to mis-specification of the neutral set, we may rely on the whole genome or a random subset of the genome to reflect neutral diversity.
  
  LSD, as currently implemented (based on ABC), is computationally demanding and requires a fair amount of computing resources. This currently limits its use to those with access to computer clusters. 

===================================================================

# Installation

As currently implemented, LSD interfaces with ms-output format coalescent simulators and ABCtoolbox. In this example, we’ll use msms as the coalescent simulator, though any of the previously described simulators may work.
 
msms is available for download at: https://www.mabs.at/ewing/msms/download.shtml

msms follows ms syntax. The authors of the program have written a convenient accessory program called ‘PopPlanner’ that facilitates visualizing ms and msms command lines into demographic models: https://www.mabs.at/ewing/msms/popplanner.shtml
 
ABCtoolbox is available for download at: https://bitbucket.org/wegmannlab/abctoolbox/wiki/Home

LSD accessory scripts are coded in Python 3 (https://www.python.org/downloads/) and R (available at https://www.r-project.org/).

Other programs that may be needed include SAMtools (http://www.htslib.org/download/) for manipulating observed data (SAM/BAM files) and ANGSD (http://www.popgen.dk/angsd/index.php/Installation) for handling low-coverage data (via genotype likelihoods) in LSD-Low.

===================================================================

# Instructions

Being reliant on ABC for parameter estimation, LSD requires summary statistics to be calculated for 1) the observed sequenced data and 2) simulated data. 

#  i) Calculate summary statistics for observed data
To calculate observed summary statistics, we supply the command with a text file containing a list of mpileup files.
For the neutral regions (1st step), this list of files will comprise of neutral or genome wide regions e.g.:

	python lsd_high_sumstats_calculator_OBS.py extractedNeutralRegions_filelist.txt -d 40 -d 40 -q 0 -m 2 -o 2popModel_observedNeutralSumStats -f ABC -r single --startPos 1 --endPos 99999 --mindepth 10 --maxdepth 500 --windowSize 5000 –pooled

For the genome scan (2nd step), we supply a list of genome or chromosome-wide sequences, and may add an additional argument (--windowStep 1000) should we wish to modulate the step-size of the sliding window e.g.:

	python lsd_high_sumstats_calculator_OBS.py genomes_filelist.txt -d 40 -d 40 -q 0 -m 2 -o 2popModel_observedSumStats -f ABC -r single --startPos 1 --endPos 99999 --mindepth 10 --maxdepth 500 --windowSize 5000 --windowStep 1000 –pooled

In these commands, we have applied the same filtering regime as in the simulated data.

For more information on the available options for LSD-High, you can run:

	python lsd_high_sumstats_calculator_OBS.py -h

To calculate observed summary statistics under low-coverage, we use ANGSD. A convenient wrapper function for this is not written yet, but one can modify LSD-Low.sh for use with observed data, to ensure that the calculation of same set of summary statistics and formatting of the output file (NEED TO UPDATE).

#  ii) Generate simulated data

   a) Generating coalescent simulations
   
   We first generate coalescent samples under a defined demographic model (see LSD requirements (1). E.g. let us assume we have 2 populations inhabiting contrasting environments, with each population comprising 20 diploid individuals each. The msms command line to generate coalescent samples for this demographic model would be:  
	
	msms 80 1 -t 10 -I 2 40 40 -n 1 1 -n 2 1 -m 1 2 M_12 -m 2 1 M_21
	
   where M is the scaled migration rate Nm (demographic parameter) that we condition the detection of selection on. Note that in msms, M as well as most other parameters are scaled to a fixed, global N<sub>E</sub> (=10,000). msms assumes haploid number of samples, hence we provide 40 haploid individuals per population here.

If parameters are not known with confidence, as will usually be the case with empirical systems, we can define the effective population sizes (and other demographic parameters, e.g. theta) as variables. Note that msms defines N in fractions of the global N<sub>E</sub> (= 10,000).

	msms 80 1 -t theta -I 2 40 40 -n 1 fraction_N1 -n 2 fraction_N2 -m 1 2 M_12 -m 2 1 M_21
	
To explore parameter space (for parameter estimation), we want these variables to be drawn from large, prior ranges. We will do this by embedding the msms command, as well as the following LSD-High command, in ABCtoolbox. 

   b) Calculating simulated summary statistics
   	
   To replicate observed sequencing pipelines, generate appropriate simulated sequencing data, and calculate a suite of summary statistics for ABC, we use LSD-High or LSD-Low. Similar to msms, LSD-High assumes haploid sample numbers. Given the simulated coalescent sample, we can generate summary statistics by e.g.:
	
	python lsd_hi.py msms_output -d 40 -d 40 -l 5000 -f ABC

   Or if we want to simulate errors (at a certain error rate), filtering, pooled samples, and a specific coverage distrubtion, we can do e.g.:
	
	python lsd_hi.py msms_output -d 40 -d 40 -l 5000 -p -i --error_method 4 --error_rate 0.001 --minallelecount 2 --mindepth 10 --maxdepth 500 --sampler nbinom -c covDist_moments.txt -f ABC
	
   where we sample according a coverage distribution fitted to the empirical coverage distribution, whose moments are described here in covDist_moments.txt (OPTIONAL). 

To acquire this fitted distribution, we first need to acquire the empirical coverage distribution. Per-site coverage of empirical data can be acquired e.g. via SAMtools

	samtools depth bam.file > bam.depthPerSite

We can then explore which theoretical distribution that best captures this observed coverage distribution. 

	# Load these libraries
	library(fitdistrplus)
	library(logspline)
	library(MASS)

	# Read in data
	raw_data <- read.table(“bam.depthPerSite”, header = FALSE, sep = "\t", row.names = NULL)

	# For testing, you may want to subsample data (since there can be millions of rows)
	raw_data <- raw_data[sample(nrow(raw_data), 100000), ]

	# Column 3 lists the coverages for the population.
	pop_data <- raw_data$V3

	# Consider if you want to apply a max depth filter
	max_depth <- 500
	pop_data_filtered <- pop_data[pop_data <= max_depth]

	# Evaluate and plot which theoretical distribution best fits data
	# Full data
	descdist(pop_data, discrete = TRUE, boot = 100)
	# Filtered data
	descdist(pop_data_filtered, discrete = TRUE, boot = 100)

	# Fit the distribution to the data and plot
	# Full data
	pop_data_nbinom <- fitdist(pop_data, "nbinom")
	plot(pop_data_nbinom)
	# Filtered data
	pop_data_filtered_nbinom <- fitdist(pop_data_filtered, "nbinom")
	plot(pop_data_filtered_nbinom)

	# Output the summary
	summary(pop_data_nbinom)
	summary(pop_data_filtered_nbinom)
	# This summary contains the moments of the distribution (e.g. mean and s.d. for normal distributions; mean and size (dispersal for negative binomial distributions). 

We specify the moments of the fitted distribution to a file; with columns representing individuals or pooled populations, the first row representing the mean of the distribution and the second row the standard deviation or dispersal. Note that sequencing data is typically best fit by a negative binomial distribution.

For more information on the available options for LSD-High, you can run:

	python lsd_hi.py -h
	
If the observed data is of low-coverage (< 10x), it is better to work with genotype likelihoods than with called genotypes, to ensure that the uncertainties in the genotypes are propagated and treated fairly. We can simulate low-coverage data from ms/msms via LSD-Low, which is a bash wrapper function for ANGSD and MsToGLF (http://www.popgen.dk/angsd/index.php/MsToGlf). MsToGLF and hence LSD-Low assume diploid sample size. Similar to LSD-High, LSD-Low allows the user to replicate (define) the depth and error rate. Additionally, the reference fasta index (of the observed data) must be supplied.

	lsd_low.sh -f msms_output -p 20,20 -l 5000 -d 5 -e 0.001 -r ref.fai -w working_dir -o output 2>&1 | tee -a log.file

c) Efficiently generating simulations with ABCtoolbox

Steps ii.a) and ii.b), that is the generation of simulated summary statistics, can be embedded and performed efficiently under ABCtoolbox. See: https://bitbucket.org/wegmannlab/abctoolbox/wiki/simulation/Performing%20Simulations%20with%20ABCtoolbox. Running these two steps under ABCtoolbox confers the convenient ability to draw variables (e.g. M and N) from defined prior ranges and thus automate the process of generating simulated data.

Example input file for generating simulations with ABCtoolbox:

	//	ABCtoolbox input file
	//	*********************

	// To set up ABCtoolbox to perform simulations
	task	simulate

	// Define the type of sampler to be used (standard, MCMC or PMC)
	samplerType	standard
	// samplerType	MCMC

	// Define .est file, which defines the priors
	estName	ABCSampler.priors

	// Define the file which contains the observed sumstats. This is to ensure that the observed data (summary statistics) and the simulated data (summary statistics) are in the same format
	obsName	/path/2popModel_observedSumStats

	//	Output file name
	outName	2pop_simpleModel_4params_simulatedSumStats

	//	Number of simulations to perform
	numSims	5000

	//	Name of the simulation program. Must be an executable.
	simProgram	/path/msms

	//	This is the msms command line we use
	//	We replace parameter (argument) values with tags defined in .est file, and define under simArgs (removing "msms")
	simArgs 80 no_loci -t theta -I 2 40 40 -n 1 fraction_N1 -n 2 fraction_N2 -m 1 2 M_12 -m 2 1 M_21

	//	This redirects the standard output to a file
	simOutputRedirection SIMDATANAME

	//	Name of the program calculating summary statistic
	//	LSD-High
	sumStatProgram /path/lsd_hi.py
	//	LSD-Low
	//	sumStatProgram /path/lsd_low.sh

	//	Arguments to be passed to the program calculating summary statistics.
	//	LSD-High
	sumStatArgs	SIMDATANAME -d 40 -d 40 -l 5000 -p -i --error_method 4 --error_rate 0.001 --minallelecount 2 --mindepth 10 --maxdepth 500 --sampler nbinom -c /path/covDist_moments.txt  -f ABC
	//	LSD-Low
	//sumStatArgs -f SIMDATANAME -p 20,20 -l 5000 -d 50 -e 0.001 -r /path/ref.fai -w . -o summary_stats_temp.txt

	//	Indicate name of summary statistics output
	sumStatName summary_stats_temp.txt

	//	Verbose output
	verbose

	//	Define name of log file
	logFile	2pop_simpleModel_4params_simulatedSumStats.log

	//Additional argument tags (required for ABC-MCMC)
	//numCaliSims 100
	//thresholdProp 0.1
	//rangeProp 1


Example priors file for generating simulations with ABCtoolbox:

	// Example ABCtoolbox priors and rules file
	// *********************

	// #### Example 4 parameter, 2 population model ####

	// msms command line in inut file:
	//	simArgs 80 no_loci -t theta -I 2 40 40 -n 1 fraction_N1 -n 2 fraction_N2 -m 1 2 M_12 -m 2 1 M_21

	[PARAMETERS]

	//	#isInt? #name	#dist.#min	#max	

	// 	Migration rates 
	//	In msms, M_i_j represents the fraction of subpopulation i that is made up of migrants from subpopulation j in forward time. Hence pastward we have the rate that a lineage moves from deme i to j as M_i_j.
	0	log_M_12	unif	-4	4 output
	0	log_M_21	unif	-4	4 output

	//	Effective population sizes 
	//	Subpopulations are defined as fractions relative to N_effective 
	0	log_N1	unif	2	7 output
	0	log_N2	unif	2	7 output

	//	Fixed parameters
	//	Here, we'll want to keep sequence length fixed (so that ABCtoolbox can conveniently vary theta for us) and bounded to a certain maximum length (since we're operating under the assumption of no recombination within-locus and free recombination between-loci)
	1	sequence_length	fixed	5000	hide
	//	We define the mutation rate at a reasonable value (e.g. 1x10-8).
	0	mutation_rate	fixed	0.00000001	hide
	//	Define the number of independent loci (iterations)
	1	no_loci	fixed	1	output


	[RULES]

	// E.g. continent population sizes should be larger than island population sizes


	[COMPLEX PARAMETERS]

	//	Various
	0	theta = (4 * 5000 * mutation_rate * sequence_length)	output

	// We output population size as absolute rather than ratio.
	0	fraction_N1 = pow10(log_N1 - 4)	output
	0	fraction_N2 = pow10(log_N2 - 4)	output

	0	M_12 = pow10(log_M_12) output
	0	M_21 = pow10(log_M_21) output

Once we have defined this input (let's call this ABCSampler.input) and priors file (let's call this ABCSampler.priors), we can run generate the simulated data with ABCtoolbox as such:

	ABCtoolbox ABCSampler.input

#  iii) Remove correlation between summary statistics
To account for potential correlation between summary statistics and to retain only their informative components, we apply a Partial Least Squares transformation. We can calculate PLS coefficients via find_pls.r. We want to find the minimum number of PLS components that explains the majority of the signal. Hence, our strategy is two run this in two steps: 

a) run for # PLS components = # of summary statistics. find_pls.r will output a plot which helps determine what the optimum number of PLS components is. 

<img src="https://github.com/hirzi/LSD/blob/master/PLS_example.png" width="400">

b) Re-run find_pls.r with the optimum number of PLS components. Be sure to modify the following lines in this script depending on the format of your summary statistics file.

		# Define working directory
		directory<-"/path/"

		# Define number of PLS components
		numComp<-5

		# Define the starting column for the summary statistics
		firstStat<-13

		# Define the columns for the free (i.e. non-fixed) parameters
		p<-c(2,3,4,5)

Observed and simulated summary statistics can then be transformed into PLS components via defining such a parameter file (let's name it as transformPLS.input):
  
	task	transform
	linearComb  PLSdef.txt
	input	raw_sumstats.txt
	output	PLStranformed_sumstats.txt
	boxcox
	numLinearComb	5

And running:
	
	ABCtoolbox transformPLS.input

#  iv) Validation of simulations
Before advancing to parameter estimation, we should first make sure that our simulated summary statistics efficiently captures that of the observed data. 

To do this, we can simply plot the simulated and observed summary statistics in summary statistic or PLS space, to assess overlap. 

<img src="https://github.com/hirzi/LSD/blob/master/Sims_vs_Obs_examplePlot.png" width="500">

	# Import libraries
	library(ggplot2)
	library(gridExtra)

	# Import data
	ABC_rej <- read.delim("./PLS_transformed_simulatedSumStats.txt.txt")
	ABC_obs<-read.delim("./PLS_transformed_observedSumStats.txt")

	# Select column with PLS LinearCombination_n
	ABC_rej<-ABC_rej[c(1:5000),c(12:ncol(ABC_rej))]
	ABC_obs<-ABC_obs[,c(2:ncol(ABC_obs))]

	# Number of PLS components (to plot!)
	num_PLS <- 4

	# For storing ggplot objects in a list using a loop: https://stackoverflow.com/questions/31993704/storing-ggplot-objects-in-a-list-from-within-loop-in-r
	plot_list <- list()
	for (i in seq(1, (num_PLS/2))) {
	  local({
	    i <- i
	    plot_list[[i]] <<- ggplot(data = ABC_rej, aes(x=ABC_rej[,2*i-1], y=ABC_rej[,2*i]) ) +
	      geom_hex(bins = 35) + scale_fill_gradientn(colours=c("gray85","gray15"),name = "sim count",na.value=NA) +
	      geom_hex(data = ABC_obs, bins = 70, aes(x=ABC_obs[,2*i-1], y=ABC_obs[,2*i], alpha=..count..), fill="red") +
	      theme_bw() + xlab(paste0("PLS ",2*i-2)) + ylab(paste0("PLS ",2*i-1)) + labs("asd")
	  })
	}

	grid.arrange(grobs = plot_list, ncol=3, top = "Overlap of simulated & observed PLS-transformed summary statistics")


#  v) ABC parameter estimation
  Perform demographic parameter estimation via ABCtoolbox. See: https://bitbucket.org/wegmannlab/abctoolbox/wiki/estimation/parameter_estimation. In our example, we seek to obtain the joint posterior of reciprocal migration rates between the 2 populations. The ABCtoolbox parameter file should reflect this accordingly.

	//----------------------------------------------------------------------
	//ABCtoolbox input file for parameter estimation
	//----------------------------------------------------------------------

	//	To estimate parameters
	task estimate

	//	Define estimation method
	//	Recall that when setting independentReplicates, the obs file should contain multiple entries (one per line), where each entry is the sumstats calculated from one (neutral) region. The algorithm then consider each region as an independent replicate, on which it subsequently calculates the posteriors.
	//estimationType standard
	//estimationType independentReplicates

	//	Observed data
	obsName	/path/example_PLS5.obs

	//	Simulated data
	simName /path/PLS_transformed_simulatedSumStats.txt

	//	Specifies the columns containing the parameters of interest in the file containing the summary statistics of the simulated data, i.e. its assigned values are numbers indicating the position of the respective columns in the file.
	params 2-3

	//	Specify the output file prefix
	outputPrefix ABC_estimation_2pop_simpleModel_

	//	Rejection settings

	//	Specifies the number of simulations in the file containing the summary statistics of the simulated data to be taken into account.
	maxReadSims 10000000

	//	Specifies the number of simulations closest to the observed data to be selected from the simulations.
	numRetained 2500

	//	Calculates the tukey depth P-value. This calculates the Tukey depth (the minimum number of sample points on one side of a hyperplane through the point, i.e. a measure of how centered a point is in an overall cloud of points) of the observed data and contrasts it with the Tukey depth of the distribution of all retained simulation points (hence argument should be equal or less than numRetained), to produce a p-value. If the observed point is close to the center of the retained simulations, we expect that most 
	//tukeyPValue 500

	//	Calculates the marginal density P-value. Similar in approach to the above, this tag calculates the P-value for the marginal density of the observed datapoint by doing so for the observed datapoint and the retained simulations (distribution)
	//marDensPValue 500

	//	If the parameter writeRetained is defined and set to 1, ABCestimator writes two files: one containing the parameter and statistics of the retained simulations and one with the smoothed parameter distribution of the retained simula- tions (see
	writeRetained 1

	//	To remove highly correlated statistics
	pruneCorrelatedStats

	//	Specifies whether (1) or not (0) the statistics are standardized before the distance is calculated.
	standardizeStats 1

	//	Posterior estimation settings
	//	Since ABCestimator standardizes the parameters internally to the range [0, 1] diracPeakWidth, the same diracPeakWidth value is used for all parameters. Too small values of diracPeakWidth will result in wiggly posterior curves, too large values might unduly smear out the curves. The best advice is to run the calculations with several choices for diracPeakWidth. The choice of diracPeakWidth depends on the number of retained simulations: the larger the number of retained parameter values, the sharper the smaller diracPeakWidth can be chosen in order to still get a rather smooth result. If the parameter diracPeakWidth is not defined, ABCestimator uses as value of 0.001, unless the parameter numRetained is defined. In this case ABCestimator sets σk = 1/N, where N is the number of simulations to retain, as proposed by Leuenberger and Wegmann (2009).
	//diracPeakWidth 0.02

	//	ABCestimator calculates the density of the marginal posteriors on a number of equally spaced points along the range of every parameter. The number of such points is specified with the parameter posteriorDensityPoints with default value 100.
	posteriorDensityPoints 100

	//	Should you wish to estimate joint posteriors
	jointPosteriors log_M_12,log_M_21
	// While we don't need the output of jointPosteriorDensityPoints since we're outputting the jointSamplesMCMC, this nonetheless needs to be set to a minimum of 2!
	jointPosteriorDensityPoints 33

	// For reference see: https://bitbucket.org/wegmannlab/abctoolbox/src/master/
	//jointSamplesMCMC 10000
	//sampleMCMCStart jointmode
	//sampleMCMCBurnin 100
	// You'll want to achieve an acceptance rate of circa 0.33 in an MCMC (check the log file for acceptance rate figures). The sampleMCMCRangeProp parameter allows you to tweak this (the higher the sampleMCMCRangeProp, the lower the acceptance rate)
	//sampleMCMCRangeProp 2
	//sampleMCMCSampling 5

	//	For cross-validation of parameter estimates. The pseudo-observed data can either be chosen among the retained simulations (retainedValidation) or among all simulations (randomValidation). The number of simulations to be used as pseudo-observed data is assigned to either one of the argument-tags.
	//randomValidation 1000

	//output settings
	logFile ABC_estimation_2pop_simpleModel.log

	verbose

Once we have defined this parameter file (let's call this file ABCEstimate.input), we can perform the parameter estimation as such:

	ABCtoolbox ABCEstimate.input

#  vi) Estimate neutral demographic parameters
  As LSD first requires an estimate of neutral demographic parameters, this must first be estimated. Following parameter estimation for the putative neutral regions (or the genome-wide representation), we combine the neutral window posterior density distributions to inform of the global (neutral) parameter estimates.

<img src="https://github.com/hirzi/LSD/blob/master/Combine_posteriors_github.png" width="625">

To acquire an estimate of the neutral (global) posteriors, we run can do as follows, taking the output of the previous step (v).

	# Import libraries
	library(MASS)
	library(RColorBrewer)

	# Define input variables
	nparams <- 2 # number of parameters

	# Directory and prefix
	master_dir <- "/Users/luqman/Desktop/ABC Simulation Study/Robustness Analyses/"
	results_dir <- "Estimation_results_2pop_simpleModel_poolSim_t10_neutral_FIXED_without0segsiteloci_PLS15_m5m5_St1_SiPop1_0.1_SiPop2_0.1_SPop1_200_SPop2_200_Chr_withSFC_jointGrid33"
	prefix <- "ABC_estimation_2pop_simpleModelmodel0_MarginalPosteriorDensities_Obs"
	setwd(paste0(master_dir,results_dir))

	####### Plot combined posteriors and all observations (loop over fraction of loci under selection and over models #######

	# Define total number of loci and fraction of loci which are under selection and neutrality. 
	total_num_loci <- 1000

	# Find product of probability densities (take log and sum)
	prod <- matrix(0, ncol=nparams, nrow=100)
	for (i in 0:(total_num_loci-1)) { # recall that ABCEstimator results are 0-indexed so no of loci - 1
	  ABC_GLM<-read.delim(paste(prefix, i, ".txt", sep = "")) 
	  for(p in 1:nparams){
	    prod[,p] <- prod[,p] + log(ABC_GLM[,2*p+1]);
	  }
	}

	# Normalise product of probability densities
	for(p in 1:nparams){
	  # Normalise and plot product of densities
	  prod[,p] <- prod[,p] - max(prod[,p]);
	  prod[,p] <- exp(prod[,p])

	}

	# Get values at peaks of product of probability densities (values)
	combined_estimate_num <- vector()
	for(p in 1:nparams){
	  comb_estimate_param_num <- (ABC_GLM[,2*p][match(max(exp(prod[,p])),exp(prod[,p]))])
	  combined_estimate_num[[p]] <- comb_estimate_param_num
	}

	# Get names and values at peaks of product of probability densities (names and values)
	combined_estimate <- list()
	for(p in 1:nparams){
	  comb_estimate_param <- (paste(colnames(ABC_GLM)[2*p], round(ABC_GLM[,2*p][match(max(exp(prod[,p])),exp(prod[,p]))],3), sep = ": "))
	  combined_estimate[[p]] <- comb_estimate_param
	}

	# Plot
	par(mfrow=c(1,2))
	for(p in 1:nparams){
	  # Normalise and plot product of densities
	  plot(ABC_GLM[,2*p], prod[,p], type='l', lty=2, col='darkred', main=names(ABC_GLM[2*p]), ylim = c(0,(max(ABC_GLM[[2*p + 1]])*2)), lwd = 2, xlab = combined_estimate[[p]]);
	  # Plot all windows
	  for (i in 0:(total_num_loci-1)) {
	    ABC_GLM<-read.delim(paste(prefix, i, ".txt", sep = ""))
	    #normalised_density <- ABC_GLM[[2*p + 1]] / max(ABC_GLM[[2*p + 1]]) # if you want to normalise the height for plotting
	    lines(ABC_GLM[[2*p]],ABC_GLM[[2*p + 1]],type='l',col='dodgerblue', lwd = 0.1)
	    #lines(ABC_GLM[[2*p]],normalised_density,type='l',col='dodgerblue', lwd = 0.1)
	  }
	}

	# Print result
	print(combined_estimate)

#  vii) Calculate deviation from neutral expectations
  Following estimation of neutral posteriors, we may then calculate the departure of window parameter estimates from neutral expectations, via lsd_scan.R (NEED TO UPDATE).

#  viii) Visualise results	
  To visualise the results, we run lsd_blotter.R (NEED TO UPDATE), which outputs a Manhattan plot of loci under selection and, if conditioned on joint (e.g. reciprocal migration) parameters, the asymmetry of the joint posterior for each loci.

<img src="https://github.com/hirzi/LSD/blob/master/Example_ManhattanPlot.png" width="1000"> 

===================================================================

# Citation

Hirzi Luqman, Alex Widmer, Simone Fior, Daniel Wegmann; bioRxiv 2020.07.20.211581; doi: https://doi.org/10.1101/2020.07.20.211581

===================================================================

# Notes

NOTE: This page is currently a work in progress and will be continuously updated, e.g. with scripts and further instructions, in the coming weeks (February 2021).
