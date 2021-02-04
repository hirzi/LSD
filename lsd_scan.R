############# This script performs a genome scan based on demographic parameter posteriors of sliding windows in a chromosome. #############
## I.e. it identifies outlier windows based on inferred demographic parameters that deviate from neutrality
## It takes as the input file the output of ABCestimator
## This version performs a genome scan even if (the regression in ABCestimator) fails for some windows (e.g. no data, 0 segsites). It does this by removing these failed windows from the scan.
## Before running, make sure to adjust the variables defined under "Define working directories" and "Define input variables and run parameters" at the top of the script.
## Additionally, lines 65,69,134, 231-232, 249-258 should be modified according to the defined paths and naming scheme.

# Read arguments from command line
#args <- commandArgs(trailingOnly = TRUE)
# To parallelise across CPUs. Assign argument to 1) chromosome number or 2) directory name
#Chr <- as.integer(args[1])
Chr <- 1  # in case assigning argument to directory name, leave this as default = 1. This has no bearing on real chromosome number, just for compatibility/legacy purposes.
#ABC_results_dir <- as.character(args[1])

# Import libraries
library(MASS)
library(stringr)
library(gtools)

############# Define working directories #############

# Working directories (parents)
masterdir <- "/path/"
ABC_parent_dir <- "/path/ABCEstimator_results/"
# Working directories (children)
ABC_results_dir <- "Estimation_results_2popModel/"
# Working directories (results)
#prefix_dir <- ""      # Define when parallelising across chromosomes
workdir <- paste0(ABC_parent_dir,ABC_results_dir)
#workdir <- paste0(ABC_parent_dir,prefix_dir,"_Chr",Chr)
output_dir <- paste0(masterdir,"GenomeScan_results/")
output_suffix <- "2popModel"

############# Define input variables and run parameters #############

# Define chromosome name
#Chr_name = Chr
Chr_name = 6
# Define whether joint posteriors were calculated via grid search "grid" or MCMC sampling "MCMC"
joint.type = "grid"
# Credible interval to define significance (for joint parameters)
credible_interval <- 0.999
#credible_interval_1D <- round(credible_interval^0.5,3) # for equivalence, we may consider taking the square root of the 2D CI for the 1D CI
credible_interval_1D <- 0.995
# Define neutral estimates - taking either 1) the marginal neutral estimates or 2) the joint neutrals from Combine_posteriors_2DJoint.R. The former is preferred because 1) it gives a better signal, 2) it has been validated through simulations with pseudo-observed data, and 3) your method of estimating the neutral estimate via the product of probabilites becomes more assumptive in higher dimensions.
# Marginal neutral estimates
neutral_estimate <- c(2, 2, 5.5, 5.5)
# Joint neutral estimates
neutral_estimate_2D <- c(2, 2)
# List of scaffolds
Chr_scaffold_list <- c("Chr6_2popModel.scaffoldlist") # this is an output of lsd_high, when run outputting for windows
# Number of windows per chromosome. This is the length of Chr6_2popModel.scaffoldlist (i.e. the number of windows)
setwd("masterdir")
raw<-read.delim("Chr6_2popModel.scaffoldlist", header = FALSE)
num_loci_list <- nrow(raw)      
# Import data on linkage group position and chromosome number - if available
#scaffold_pos_raw<-read.delim("anchored_scaffolds.txt", sep = "\t", header = TRUE)

############# Run Genome Scan for nChr chromosomes #############

# Define number of windows or loci in chromosome
num_loci <- num_loci_list[Chr]

# Set working directory to point to marginalPosteriorDensities results
setwd(paste0(workdir,"/marginalPosteriorDensities"))

# In case regression in ABCestimator fails for some loci (which can happen), we'll iterate over all succesfully estimated loci rather than all observations.
ABC_output_files <- mixedsort(list.files())
idx_loci_passed <- grep("MarginalPosteriorDensities", ABC_output_files)
# Check that indexes are sorted as expected
sort_result <- !is.unsorted(idx_loci_passed)
print(paste0("Loci are sorted by observation: ", sort_result))
# Define number of windows or loci in chromosome
#num_loci <- num_loci_list[Chr]
num_loci_passed <- length(idx_loci_passed)

# Acquire marginal output format and number of parameters
for (i in ABC_output_files[idx_loci_passed][1]) {
  ABC_GLM_marginal_init <- read.delim(i, sep = "")
  nparams <- (ncol(ABC_GLM_marginal_init)-1)/2
}

############# Calculate posterior density peaks and marginal (1D) confidence intervals #############
posterior_peaks_matrix <- matrix(0, ncol=nparams, nrow=num_loci_passed)
#CI_matrix <- matrix(0, ncol=nparams, nrow=num_loci_passed)
significance_HDR_matrix <- matrix(0, ncol=nparams, nrow=num_loci_passed)
significance_QI_matrix <- matrix(0, ncol=nparams, nrow=num_loci_passed)
neutral_CI_matrix <- matrix(0, ncol=nparams, nrow=num_loci_passed)
#neutral_density_matrix <- matrix(0, ncol=nparams, nrow=num_loci_passed)

counter = 1
for (i in ABC_output_files[idx_loci_passed]) {
  for(p in 1:nparams){
    ABC_GLM<-read.delim(paste(i, sep = "")) 
    # Find peaks of posterior
    posterior_peak <- ABC_GLM[,2*p][match(max(ABC_GLM[,2*p+1]),ABC_GLM[,2*p+1])] # May also want to consider an alternative estimator than the max
    # Output results to matrix
    posterior_peaks_matrix[counter,p] <- posterior_peak
    # Calculate 1D marginal credible interval (highest density region and quantile interval)
    # Highest density region CI
    const <- sum(ABC_GLM[,2*p+1])
    spxx <- sort(ABC_GLM[,2*p+1], decreasing = TRUE) / const
    HDR_CI <- spxx[which(cumsum(spxx) >= credible_interval_1D)[1]] * const
    # % credible interval at neutral estimate (as a continuous measure of significance)
    neutral_density_index <- which.min(abs(ABC_GLM[,2*p] - neutral_estimate[p]))
    neutral_density <- ABC_GLM[,2*p+1][neutral_density_index]
    neutral_index <- match((neutral_density / const), spxx)
    neutral_CI <- cumsum(spxx)[neutral_index]   # this variable indicates at which HDR credible interval the neutral estimate lies (in fraction)
    neutral_CI_matrix[counter,p] <- neutral_CI
    # Quantile interval
    cpxx <- cumsum(ABC_GLM[,2*p+1]) / sum(ABC_GLM[,2*p+1])
    lower_bound_CI <- ABC_GLM[,2*p][which(cpxx >= (0 + ((1 - credible_interval_1D) / 2)))[1]] 
    upper_bound_CI <- ABC_GLM[,2*p][which(cpxx >= (1 - ((1 - credible_interval_1D) / 2)))[1]-1]
    # Calculate density at neutral estimate
    neutral_closestObs <- which.min(abs(ABC_GLM[[2*p]] - neutral_estimate[p]))
    neutral_est_density <- ABC_GLM[[2*p+1]][neutral_closestObs]
    # Calculate whether neutral point lies within 1D credible interval
    neutrality_inf_HDR <- neutral_est_density > HDR_CI
    neutrality_inf_QI <- (upper_bound_CI > neutral_estimate[p]) & (neutral_estimate[p] > lower_bound_CI)
    # Output results to matrix
    #CI_matrix[counter,p] <- HDR_CI
    #neutral_density_matrix[counter,p] <- neutral_est_density
    significance_HDR_matrix[counter,p] <- neutrality_inf_HDR
    significance_QI_matrix[counter,p] <- neutrality_inf_QI
  }
  counter = counter + 1
}

############# Calculate 2D joint posterior confidence intervals #############
# It may be more relevant/powerful to draw credible intervals on the 2D joint posterior. We do this by performing a 2D KDE.
# Here we assume the case of 2 joint parameters. For >2 joint parameters, modify code accordingly

# Change working directory to point to jointPosterior results
setwd(paste0(workdir,"/jointPosteriors"))

# Initiate empty vectors
significance_KDE2D_vector <- vector()
significance_KDE2D_neutralCI_vector <- vector()
asymmetry_KDE2D_vector <- vector()
counter = 1

# In case regression in ABCestimator fails for some loci (which can happen), we'll iterate over all succesfully estimated loci rather than all observations.
ABC_output_files_joint <- mixedsort(list.files())
if (joint.type == "grid") {
  idx_loci_passed_joint <- grep("jointPosterior", ABC_output_files_joint)
} else if (joint.type == "MCMC") {
  idx_loci_passed_joint <- grep("jointPosteriorSamples", ABC_output_files_joint)
}
# Check that indexes are sorted as expected
sort_result_joint <- !is.unsorted(idx_loci_passed_joint)
print(paste0("Loci are sorted by observation: ", sort_result_joint))
# Define number of windows or loci in chromosome
#num_loci <- num_loci_list[Chr]
num_loci_passed_joint <- length(idx_loci_passed_joint)

# Acquire prior range for joint parameters (here assuming 2 joint parameters) and joint output format
# Also acquire the size of the grid along one axis (e.g. the number of sampled points per parameter). Assuming square grid, we can simply take the square root of the number of samples.
for (i in ABC_output_files_joint[idx_loci_passed_joint][1]) {
  ABC_GLM_joint_init <- read.delim(i, sep = "")
  prior_jointParam1 <- c(round(min(ABC_GLM_joint_init[2])),round(max(ABC_GLM_joint_init[2])))
  prior_jointParam2 <- c(round(min(ABC_GLM_joint_init[3])),round(max(ABC_GLM_joint_init[3])))
  density_points = sqrt(nrow(ABC_GLM_joint_init))
}

for (i in ABC_output_files_joint[idx_loci_passed_joint]) {
  ABC_GLM<-read.delim(i, sep = "")
  # Specify whether joint parameters were calculated via grid search or MCMC sampling
  if (joint.type == "grid") {
    dens <- list()
    # This (x-axis) is refers to log_m_highcontinent_lowcontinent
    dens[[1]] <- seq(prior_jointParam1[1],prior_jointParam1[2], length.out = density_points)
    # This (y-axis) is refers to log_m_lowcontinent_highcontinent
    dens[[2]] <- seq(prior_jointParam2[1],prior_jointParam2[2], length.out = density_points)
    # We transform the density vector into matrix.
    dens[[3]] <- matrix(ABC_GLM$density,nrow=density_points, ncol = density_points)
    names(dens) <- c("x", "y", "z")
  } else if (joint.type == "MCMC") {
    dens <- kde2d(ABC_GLM[,2], ABC_GLM[,3], n = density_points)
  }
  # Calculate and output asymmetry.
  if (joint.type == "grid") {
    # Recall that there are points above the diagonal, below the diagonal and on the diagonal. To reduce degrees of freedom, we first account for the points on the diagonal (i.e. where M_12 = M_21).
    symmetry <- sum((ABC_GLM[,2] == ABC_GLM[,3])*ABC_GLM$density)
    # Then the integral of the area above the diagonal equals one minus the integral of the area below the diagonal. So we can just output one variable.
    asymmetry <- sum((ABC_GLM[,2] < ABC_GLM[,3])*ABC_GLM$density) / (sum(ABC_GLM$density) - symmetry)
  } else if (joint.type == "MCMC") {
    # Recall that there are points above the diagonal, below the diagonal and on the diagonal. To reduce degrees of freedom, we first account for the points on the diagonal (i.e. where M_12 = M_21).
    symmetry <- sum(ABC_GLM[,2] == ABC_GLM[,3])
    # Then the integral of the area above the diagonal equals one minus the integral of the area below the diagonal. So we can just output one variable.
    asymmetry <- sum(ABC_GLM[,2] < ABC_GLM[,3]) / (length(ABC_GLM[,2]) - symmetry)
  }
  # Append asymmetry statistic to results vector
  asymmetry_KDE2D_vector[counter] <- asymmetry
  # Output intersection point of neutral estimate and 2D joint distribution
  neutral_contour_x <- which.min(abs(dens$x - neutral_estimate_2D[1]))
  neutral_contour_y <- which.min(abs(dens$y - neutral_estimate_2D[2]))
  neutral_contour_xy <- dens$z[neutral_contour_x, neutral_contour_y]
  # Calculate credible interval (here: highest density interval)
  const2d <- sum(dens$z)
  spxx2d <- sort(dens$z, decreasing = TRUE) / const2d
  crit <- spxx2d[which(cumsum(spxx2d) >= credible_interval)[1]] * const2d
  # Append significance to results vector
  neutrality_inf_KDE2D <- neutral_contour_xy > crit
  significance_KDE2D_vector[counter] <- neutrality_inf_KDE2D
  # Additionally, we indicate at which % credible interval the neutral estimate lies (i.e. the area above the horizontal threshold defined by the neutral estimate)
  neutral_index2d <- match((neutral_contour_xy / const2d), spxx2d)
  neutral_CI_2d <- cumsum(spxx2d)[neutral_index2d]   # this variable indicates at which HDR credible interval the neutral estimate lies (in fraction)
  significance_KDE2D_neutralCI_vector[counter] <- neutral_CI_2d
  counter = counter + 1
}
significance_KDE2D_vector <- as.numeric(significance_KDE2D_vector)

############# For ploting #############
# We want to colour the points based on the credible interval overlap. We will do this by subsetting on matrix significance_matrix, which is currently defined as 1 = TRUE and 0 = FALSE). Since R indexes from 1 (not from 0), we add 1 to all elements of the matrix.
significance_HDR_matrix <- significance_HDR_matrix + 1
significance_QI_matrix <- significance_QI_matrix + 1
significance_KDE2D_vector <- significance_KDE2D_vector + 1

# To color points (windows) by their significance values (% credible interval where the neutral estimate lies), we first rank the variable of interest for colour assignment
significance_KDE2D_neutralCI_df <- as.data.frame(significance_KDE2D_neutralCI_vector)
significance_KDE2D_neutralCI_df$order = findInterval(significance_KDE2D_neutralCI_df$significance_KDE2D_neutralCI_vector, sort(significance_KDE2D_neutralCI_df$significance_KDE2D_neutralCI_vector))

# So that the lower values occupy more space, we transform the data as such: -log10(1-CI)
neutral_logCI_matrix <- -log10(1 - neutral_CI_matrix)
significance_KDE2D_neutralCI_df[,3] <- -log10(1 - significance_KDE2D_neutralCI_df[,1])

############# Make and output genomeScan table #############

# Read in scaffold list information. To extract only relevant information from scaffold/window name, change separator accordingly (depending on naming scheme of scaffold list)
#setwd(masterdir)
setwd(paste0(masterdir,"Chromosomes/"))
window_pos<-read.delim(Chr_scaffold_list[Chr], sep = "l", header = FALSE) # depending on naming scheme of scaffold list, change separator accordingly

# Remove scaffolds if no posterior was calculated (i.e. if regression failed during ABCestimator)
first_passed_obs.raw <- ABC_output_files[idx_loci_passed[1]]
first_passed_obs.temp <- str_split_fixed(first_passed_obs.raw, "Obs", 2)
first_passed_obs <- str_split_fixed(first_passed_obs.temp[2], ".txt", 2)
first_passed_obs <- as.numeric(first_passed_obs[1])
idx_loci_passed_MarginalPosteriorDensities <- idx_loci_passed - idx_loci_passed[1] + (1 + first_passed_obs)
first_passed_obs_joint.raw <- ABC_output_files_joint[idx_loci_passed_joint[1]]
first_passed_obs_joint.temp <- str_split_fixed(first_passed_obs_joint.raw, "Obs", 2)
first_passed_obs_joint <- str_split_fixed(first_passed_obs_joint.temp[2], ".txt", 2)
first_passed_obs_joint <- as.numeric(first_passed_obs_joint[1])
idx_loci_passed_jointPosteriors <- idx_loci_passed_joint - idx_loci_passed_joint[1] + (1 + first_passed_obs_joint)
# Merge (find union) of idx_loci_passed (marginal and joint)
idx_loci_passed_combined <- intersect(idx_loci_passed_MarginalPosteriorDensities,idx_loci_passed_jointPosteriors)
window_pos <- window_pos[idx_loci_passed_combined,]

# Select relevant data (select relevant column from scaffold list)
window_pos<-as.vector(window_pos[,3])
# Remove leading and trailing (useless) characters
scaffold_list.temp <- lapply(window_pos, function(x) substring(x, 2, nchar(x)-4))
scaffold_list <- as.data.frame(t(as.data.frame(scaffold_list.temp)))
rownames(scaffold_list)<-c(1:length(scaffold_list[,1]))
# Split into scaffold name and window range
scaffold_df <- scaffold_list
scaffold_df.windows.temp <- as.data.frame(str_split_fixed(scaffold_list$V1, "window_", 2))
scaffold_df.windows <- as.data.frame(str_split_fixed(scaffold_df.windows.temp$V2, "_", 2)) 

# Define dataframe
# Append scaffold name
genomeScan_table <- as.data.frame(scaffold_df$V1)
# Append scaffold window ranges (start and end)
genomeScan_table[,ncol(genomeScan_table)+1] <- scaffold_df.windows$V1
genomeScan_table[,ncol(genomeScan_table)+1] <- scaffold_df.windows$V2
# Append CI values (1D and 2D)
marginalParam1_idx <- as.numeric(which(colnames(ABC_GLM_marginal_init) == colnames(ABC_GLM_joint_init)[2]))
marginalParam2_idx <- as.numeric(which(colnames(ABC_GLM_marginal_init) == colnames(ABC_GLM_joint_init)[3]))
param_indexes <- c(marginalParam1_idx/2,marginalParam2_idx/2)
for (p in param_indexes) {
  genomeScan_table[,ncol(genomeScan_table)+1] <- neutral_logCI_matrix[,p]
}
genomeScan_table[,ncol(genomeScan_table)+1] <- significance_KDE2D_neutralCI_df[,3]
# Append binary significance states based on 95% CI - 2:not significant, 1: significant (1D and 2D)
for (p in param_indexes) {
  genomeScan_table[,ncol(genomeScan_table)+1] <- significance_HDR_matrix[,p]
}
genomeScan_table[,ncol(genomeScan_table)+1] <- significance_KDE2D_vector
# Append migration rate posterior modes 
for (p in param_indexes) {
  genomeScan_table[,ncol(genomeScan_table)+1] <- posterior_peaks_matrix[,p]
}
# Append asymmetry stastistic
genomeScan_table[,ncol(genomeScan_table)+1] <- asymmetry_KDE2D_vector

# Add linkage group position and chromosome number
genomeScan_table[,ncol(genomeScan_table)+1] <- NA      # scaffold name
genomeScan_table[,ncol(genomeScan_table)+1] <- NA      # chromosomal position or position on linkage group
genomeScan_table[,ncol(genomeScan_table)+1] <- paste0("Chr",Chr_name)  # chromosome or linkage group

# Add column headers
marginalParam1_prefix <- substring(colnames(ABC_GLM_marginal_init)[marginalParam1_idx],5)
marginalParam2_prefix <- substring(colnames(ABC_GLM_marginal_init)[marginalParam2_idx],5)
jointParams_prefix <- paste0(marginalParam1_prefix,"_X_",marginalParam2_prefix)

colnames(genomeScan_table) <- c("scaffold", "window_start", "window_end", paste0("-log(1-p_",marginalParam1_prefix,")"), paste0("-log(1-p_",marginalParam2_prefix,")"), paste0("-log(1-p_joint_",jointParams_prefix,")"), paste0("significance_",credible_interval_1D,"_",marginalParam1_prefix), paste0("significance_",credible_interval_1D,"_",marginalParam2_prefix), paste0("significance_",credible_interval,"_joint_",jointParams_prefix), paste0("mode_",marginalParam1_prefix), paste0("mode_",marginalParam2_prefix), "asymmetry", "pos_index", "Chr_position", "Chr")
# Reorder dataframe columns
genomeScan_table <- subset(genomeScan_table, select =  c("scaffold", "Chr", "Chr_position", "window_start", "window_end",  paste0("-log(1-p_",marginalParam1_prefix,")"), paste0("-log(1-p_",marginalParam2_prefix,")"), paste0("-log(1-p_joint_",jointParams_prefix,")"), paste0("significance_",credible_interval_1D,"_",marginalParam1_prefix), paste0("significance_",credible_interval_1D,"_",marginalParam2_prefix), paste0("significance_",credible_interval,"_joint_",jointParams_prefix), paste0("mode_",marginalParam1_prefix), paste0("mode_",marginalParam2_prefix), "asymmetry") )
# Redefine column types. To convert factor to numeric, need to convert to character first (because factors are stored internally as integers with a table to give the factor level labels)
#sapply(genomeScan_table, class)
genomeScan_table$window_start <- as.numeric(as.character(genomeScan_table$window_start))
genomeScan_table$window_end <- as.numeric(as.character(genomeScan_table$window_end))

# In case you'd like to change sort order, e.g. sort by chromosome position (first) and windows (second)
#genomeScan_table <- genomeScan_table[with(genomeScan_table, order(Chr_position, window_start)), ]

# Write out table
setwd(output_dir)
write.table(genomeScan_table, file = paste0("genomeScan_table_Chr",Chr_name,"_",output_suffix,".txt"), row.names = FALSE)
