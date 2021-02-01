# Identifying Loci under Selection via explicit Demographic models (LSD)

  This repository contains a suite of scripts for performing LSD genome scans based on explicit demographic models (Luqman et al. 2021, https://www.biorxiv.org/content/10.1101/2020.07.20.211581v2.full.pdf). The current implementation estimates demographic parameters via an Approximate Bayesian Computation (ABC) framework, and works in two steps. First, neutral demographic parameters are estimated (see requirements for LSD point (ii)). Second, per-locus parameter estimates (e.g. for a sliding window across the chromosome) are compared to the neutral estimates, to identify selected loci. 

  As LSD is an ABC approach, it relies on simulations to estimate the posterior distribution of model parameters. The current implementation takes ms-format coalescent samples as input for simulated data and mpileup format (e.g. from BAM files) as input for observed data. A large range of modern coalescent simulators (or those that approximate the coalescent) output ms-format data including e.g. ms (Hudson, 2002), msHOT (Hellenthal & Stephens, 2007), msms (Ewing & Hermisson, 2010), msprime (Kelleher & Etheridge, 2015), MaCS (Chen, Marjoram, & Wall, 2009), cosi2 (Shlyakhter et al., 2014) and SCRM (Staab et al., 2015). 

  The processing, format and final output of observed genetic data will often differ from that of raw coalescent simulations, given that observed genetic data may be subject to various pre-sequencing (e.g. pooling), sequencing (e.g. sequencing errors, stochastic sampling of reads) and post-sequencing (e.g. filters) events that perturb and reformat the data from the original source. We thus provide programs that interface with coalescent simulators to replicate observed sequencing pipelines and generate simulated sequencing data. LSD-High can accommodate and simulate both individual and pooled data and assumes mid to high coverage (>10x) data, while LSD-Low accepts individual data and can additionally accommodate low coverage (>2x) data by utilising genotype likelihoods via msToGLF and ANGSD (Korneliussen, Albrechtsen, & Nielsen, 2014). These programs then calculate a suite of summary statistics for the simulated and observed data. Summary statistics currently implemented include the number of segregating sites (S), private S, nucleotide diversity (pi), Watterson’s theta estimator, Tajima’s D, relative divergence (FST), absolute divergence (DXY), and site frequencies, though in principle any summary statistic can be included with appropriate additions or modifications to the programs’ scripts. 

  ABC is currently implemented via ABCtoolbox (Wegmann et al., 2010). LSD-blotter takes the output of the ABC parameter estimates and estimates the departure of the inferred posteriors from neutral expectations. If conditioning the detection of selected (and linked) loci on multiple (joint) parameters, LSD also outputs the directionality of the deviation in the joint posterior with respect to the marginal parameters, and represents these as colours in the genome scan Manhattan plot.

  We note that LSD is NOT a program, rather it is an analytical framework for identifying loci under selection via deviations in demographic parameters. As such, it is not constrained to any particular program. Rather, we envision a custom and modular implementation that may interface with any appropriate combination of coalescent simulator, summary statistics calculator and ABC program. Here, we simply propose an example implementation utilising msms, LSD-High/Low and ABCtoolbox. As currently implemented, LSD takes WGS data and bases its inference of selection on genomic windows (regions), rather than SNPs. That said, it can be extended to work on SNPs if implemented with a SNP-based summary statistics calculator and appropriate formatting of the coalescent simulations, though this would entail likely entail some scripting/programming.

============================================

# Requirements for LSD

  Firstly, a demographic model needs to be defined. Definition and choice of the demographic model should i) be informed by knowledge of the study system, ii) be motivated by the model’s capacity to provide a useful approximation of a biological process of interest, and iii) be sufficiently simple to remain computationally tractable. Additionally, given that we condition the inference of selection on demographic parameters, the model should be formulated according to whether deviation in N_E or M_E is desired for the inference of selection. Finally, the model should be able to sufficiently describe the neutral genetic variation of the system. This can be validated by demonstrating that the observed data can be accurately and sufficiently captured by the simulated data (e.g. in terms of summary statistics)

  To estimate neutral demographic parameters in the first step, we need an a priori set of regions that we believe reflect neutral evolution. Such a set may be informed by the particular structural or functional class the sites belong to and may e.g. consist of genomic regions not linked to structural annotations. Alternatively, given that LSD is very robust to mis-specification of the neutral set, we may rely on the whole genome or a random subset of the genome to reflect neutral diversity. 

============================================

# Installation

As currently implemented, LSD interfaces with ms-output format coalescent simulators and ABCtoolbox. In this example, we’ll use msms as the coalescent simulator, though any of the previously described simulators may work.
 
msms is available for download at: https://www.mabs.at/ewing/msms/download.shtml

msms follows ms syntax. The authors of the program have written a convenient accessory program called ‘PopPlanner’ that facilitates visualizing ms and msms command lines into demographic models: https://www.mabs.at/ewing/msms/popplanner.shtml
 
ABCtoolbox is available for download at: https://bitbucket.org/wegmannlab/abctoolbox/wiki/Home

============================================

# Instructions

We first generate coalescent samples under a defined demographic model. Being reliant on ABC for parameter estimation, LSD requires summary statistics to be calculated for 1) the observed sequenced data and 2) simulated data. 

#  i) Generate simulated data
   a) Generating coalescent simulations
   
   We first generate coalescent samples under a defined demographic model (see LSD requirements (1). E.g. let’s assume we have 2 populations inhabiting contrasting environments, with each population comprising 20 individuals each. The msms command line to generate coalescent samples for this demographic model would be e.g.:  
	
	msms 80 1 -t 10 -I 2 40 40 -n 1 1 -n 2 1 -m 1 2 M -m 2 1 M
	
   where M is the migration rate (demographic parameter) that we condition the detection of selection on, and should be drawn from a reasonably large prior range.
   Here, we simulate a single locus, and hence assume no recombination between loci and fixed recombination within locus.

   b) Calculating simulated summary statistics
   
   To replicate observed sequencing pipelines, generate appropriate simulated sequencing data, and calculate a suite of summary statistics for ABC, we use LSD-High or LSD-Low. Given the simulated coalescent sample, we can generate summary statistics by e.g.:
	
	python3 lsd_hi.py msms_output -d 40 -d 40 -l 5000 -f ABC

   Or if we want to simulate errors (at a certain error rate), filtering, pooled samples, and a specific coverage distrubtion, we can do e.g.:
	
	python3 lsd_hi.py msms_output -d 40 -d 40 -l 5000 -p -i --error_method 4 --error_rate 0.001 --minallelecount 2 --mindepth 10 --maxdepth 500 --sampler nbinom -c covDist_moments.txt -f ABC
	
   where we sample according a coverage distribution fitted to the empirical coverage distribution, whose moments are described here in covDist_moments.txt. Elaborate…

Steps a) and b), that is the generation of simulated summary statistics, can be embedded and performed efficiently under ABCtoolbox. See: https://bitbucket.org/wegmannlab/abctoolbox/wiki/simulation/Performing%20Simulations%20with%20ABCtoolbox. Provide example.

#  ii) Calculating observed summary statistics
To calculate observed summary statistics, we supply the command with a text file containing a list of mpileup files.
For the neutral regions (1st step), this list of files will comprise of neutral or genome wide regions e.g.:

	python lsd_high_sumstats_calculator_OBS.py extractedNeutralRegions_filelist.txt -d 40 -d 40 -q 0 -m 2 -o 2popModel -f ABC -r single --startPos 1 --endPos 99999 --mindepth 10 --maxdepth 500 --windowSize 5000 –pooled

For the genome scan (2nd step), we supply a list of genome or chromosome-wide sequences, and may add an additional argument (--windowStep 1000) should we wish to modulate the step-size of the sliding window e.g.:

	python lsd_high_sumstats_calculator_OBS.py genomes_filelist.txt -d 40 -d 40 -q 0 -m 2 -o 2popModel -f ABC -r single --startPos 1 --endPos 99999 --mindepth 10 --maxdepth 500 --windowSize 5000 --windowStep 1000 –pooled

In these commands, we have applied the same filtering regime as well as mimicked the observed error rate as well as coverage distribution.

#  iii) Remove correlation between summary statistics
  To account for potential correlation between summary statistics and to retain only their informative components, we apply a Partial Least Squares transformation (Wegmann, Leuenberger, & Excoffier, 2009). We can calculate PLS coefficients via find_pls.r. Observed and simulated summary statistics can then be transformed into PLS components via the ABCTransform scripts.
  
#  iv) Validation of simulations
  Before advancing to parameter estimation, we should first make sure that our simulated summary statistics efficiently captures that of the (neutral or genome-wide) observed data. To do this, we can simply plot the simulated and observed summary statistics in summary statistic or PLS space, to assess overlap (script). 

#  v) ABC parameter estimation
  Perform demographic parameter estimation via ABCtoolbox. See: https://bitbucket.org/wegmannlab/abctoolbox/wiki/estimation/parameter_estimation. In our example, we seek to obtain the joint posterior of reciprocal migration rates between the 2 populations. The ABCtoolbox parameter files should reflect this accordingly.

#  vi) Estimating neutral demographic parameters
  In LSD, steps ii)-v) above are first carried out assuming the observed data to constitute neutral (or genome-wide) regions. Assuming such, step v) generates the parameter posterior distributions for numerous putative neutral windows. To acquire an estimate of the neutral (global) posteriors, we run x.script.

#  vii) Calculating deviation from neutral expectations
  Following estimation of neutral posteriors, we may then calculate the departure of window parameter estimates from neutral expectations. This assumes that step ii) has been performed where the observed data constitutes the full while genome or chromosome data. This can be run via the xxx scripts (this needs to be revised such that departure is from the diagonal defined by the neutral point). 

#  viii) Visualise results	
  To visualise the results, we run scriptx, which outputs a Manhattan plot of loci under selection and, if conditioned on joint (e.g. reciprocal migration) parameters, the asymmetry of the joint posterior for each loci.


============================================

# Citation

============================================

# Notes

NOTE: This page is currently a work in progress and will be continuously updated, e.g. with scripts and further instructions, in the coming weeks (February 2021).
