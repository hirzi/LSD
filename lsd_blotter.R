
############# This script draws Manhatten plot for genome scans, for single chromosome plots #############

# Note 1: The input is a single chromomome genome scan table, in text format.
# Note 2: in the genome scan input table, a significance state of "1" implies significance and of "2" implies non-significant.
# Note 3: Before running, make sure to adjust the variables defined under plotting parameters at the top of the script.
# Note 4: This script additionally colors points above a defined threshold according to the symmetry of the posterior distribution.

# Load libraries
library(RColorBrewer)
library(DescTools)
library(stringr)

# Define directory containing genome scan table files
setwd("/Users/luqman/Downloads/")

### Plotting parameters

# Define chromosome to plot
#iChr <- 6

# Should you wish to plot a horizontal line representing the credible interval, we can set this here. Note, this DOES NOT define the significance states of the input table, rather this is only used for plotting purposes.
credible_interval <- 0.999
#credible_interval_1D <- round(credible_interval^0.5,3)
credible_interval_1D <- 0.995

# Cap maximum likelihood, in case some values are extremely high.
likelihood_max_cap_2D <- 6
likelihood_max_cap_1D <- 6

### Plotting main function
#for (Chr in seq(iChr,iChr)) {

#input_file <- paste0("genomeScan_table_Chr",iChr,suffix,".txt")
#input_file <- "genomeScan_table_Chr6_region_40MB_56MB_2pop_simpleModel_retSims6750.txt"
input_file <- "genomeScan_table_Chr6_region_40MB_56MB_6pop_IMmodel_m_4_m4_1million_retSims5000.txt"
#input_file <- "genomeScan_table_Chr6_region_52MB_54MB_2pop_simpleModel_retSims4000.txt"
#input_file <- "genomeScan_table_Chr6_region_52MB_54MB_2pop_simpleModel_retSims6750.txt"
#output_suffix <- paste0("Chr", Chr,"xxxModel")
#output_suffix <- "Chr6_region_52MB_54MB_6pop_IMmodel_m_4_m4_retSims5000"
#output_suffix <- "Chr6_region_52MB_54MB_2pop_simpleModel_retSims6750"
#output_suffix <- "Chr6_region_40MB_56MB_2pop_simpleModel_retSims6750"
output_suffix <- "Chr6_region_40MB_56MB_6pop_IMmodel_retSims5000"

# Read in genome scan table
genomeScan_table <- read.delim(input_file, sep = " ")
# Check data classes
sapply(genomeScan_table, class)
# Rename headers (in case header special characters modified)
colnames(genomeScan_table) <- c(colnames(genomeScan_table)[1], colnames(genomeScan_table)[2], colnames(genomeScan_table)[3], colnames(genomeScan_table)[4], colnames(genomeScan_table)[5],  paste0(gsub('.{1}$', '', paste0("-log(1-",substring(colnames(genomeScan_table)[6],9))),")"), paste0(gsub('.{1}$', '', paste0("-log(1-",substring(colnames(genomeScan_table)[7],9))),")"), paste0(gsub('.{1}$', '', paste0("-log(1-",substring(colnames(genomeScan_table)[8],9))),")"), colnames(genomeScan_table)[9], colnames(genomeScan_table)[10], colnames(genomeScan_table)[11], colnames(genomeScan_table)[12], colnames(genomeScan_table)[13], colnames(genomeScan_table)[14])

# Cap max -log(1-p) to a reasonable value, e.g. 5
cap_max_prob <- function(x, likelihood_max_cap){
  if(x > likelihood_max_cap){
    x <- likelihood_max_cap - runif(1,0,1) # add jitter to distinguish between capped points
    return(x)
  }
  else {
    return(x)
  }
}
for (prob_column in c(6,7)) {
  genomeScan_table[,prob_column] <- sapply(genomeScan_table[,prob_column], cap_max_prob, likelihood_max_cap=likelihood_max_cap_1D)
}
for (prob_column in c(8)) {
  genomeScan_table[,prob_column] <- sapply(genomeScan_table[,prob_column], cap_max_prob, likelihood_max_cap=likelihood_max_cap_2D)
}

# To color significant values by asymmetry
genomeScan_table[,15] <- (genomeScan_table[,11]== 1)*genomeScan_table$asymmetry
genomeScan_table[,15][genomeScan_table[,15] == 0] <- NA
# Define number of colour breaks and colour palette
col_breaks <- 11
#genomeScan_table[,16] <- as.numeric(cut(genomeScan_table[,15],breaks = col_breaks)) # breaks with min/max value = min/max of data
genomeScan_table[,16] <- round(genomeScan_table[,15]*10) + 1  # breaks with min/max value = 0/1 (fixed interval)
col_pal <- colorRampPalette(c('red3','snow','dodgerblue3'))
#col_pal <- colorRampPalette(c('red','black','blue'))
#cols <- brewer.pal(col_breaks, "RdYlBu")
#col_pal <- colorRampPalette(cols)
genomeScan_table[,17] <- col_pal(col_breaks)[genomeScan_table[,16]]
genomeScan_table[,18] <- genomeScan_table[,17]
genomeScan_table[,18][!is.na(genomeScan_table[,18])] <- "black"

# Let's add position (take midpoint of window)
genomeScan_table$Chr_position <- (genomeScan_table$window_start + genomeScan_table$window_end) / 2

# Manhatten plot based on joint 2D credible intervals
par(mfrow=c(1,1),mar=c(4,4,4,6))
prefix_joint <- str_split(gsub('.{1}$', '', substring(colnames(genomeScan_table)[8],16)),"_X_")[[1]]
#plot(genomeScan_table[,3], genomeScan_table[,8], xlab="Position", ylab="-log10(1-P)", main = paste0("Joint 2D ", prefix_joint[1], " - ", prefix_joint[2], " significance; ", output_suffix), pch=19, cex = 0.3, col = "grey50")
plot(genomeScan_table[,3], genomeScan_table[,8], xlab="Position", ylab="-log10(1-P)", main = "Model M1", pch=19, cex = 0.3, col = "grey50")
points(genomeScan_table[,3], genomeScan_table[,8], pch=21, cex = 1.25, col=genomeScan_table[,18], bg=genomeScan_table[,17])
abline(h = -log10(1-credible_interval), lty = 2, lwd = 1.5, col = "grey30")
ColorLegend("right", col=col_pal(col_breaks*5), labels=sprintf("%.1f",seq(0,1,0.1)), adj = c(0,0.5), cntrlbl = TRUE, cex=0.5, inset = -0.04, width = round(nrow(genomeScan_table)*40), height = min(likelihood_max_cap_2D,max(genomeScan_table[,8]))*1.08)

# Manhatten plot based on 1D credible intervals
prefix_marginals <- c(gsub('.{1}$', '', substring(colnames(genomeScan_table)[6],10)), gsub('.{1}$', '', substring(colnames(genomeScan_table)[7],10)))
par(mfrow=c(1,1),mar=c(4,4,4,4))
for (p in 1:2) {
  # We want to distinguish between windows that satisfy a single parameter CI signficance and those that satisfy the joint (2) parameter CI significance.
  # We do this by adding a third value (3), which signifies windows that satisfy the 1D CI but not the joint 2D CI
  signficance_combined_vector <- vector()
  for (i in seq(1, nrow(genomeScan_table))) {
    if ( (genomeScan_table[,p+8][i] == 1) & (genomeScan_table[,11][i] == 2)) {
      signficance_combined_vector[i] <- 3
    }
    else {
      signficance_combined_vector[i] <- genomeScan_table[,11][i]
    }
  }
  # And plot
  plot(seq(1,nrow(genomeScan_table)), genomeScan_table[,p+5], xlab="Position", ylab="-log10(1-P)", main = paste0(prefix_marginals[p], " (marginal) significance; ", output_suffix), pch=19, cex = 0.3, col=c("red3", "grey35","darkorange")[signficance_combined_vector])
  abline(h = -log10(1-credible_interval_1D), lty = 2, lwd = 1.5, col = "grey30")
  legend("topleft", legend = c(paste0("Not significant(<",credible_interval*100,"^0.5% CI)"), paste0("Significant for this marginal parameter (>", credible_interval*100,"^0.5% CI)"), paste0("Significant for (both) joint parameters (>", credible_interval*100,"% CI)")), col=c("grey35","darkorange", "red3"), pch = 19, bty = "n", pt.cex = 1.25, cex = 0.6, horiz = FALSE, inset = c(0.01, 0.025))
}
#}