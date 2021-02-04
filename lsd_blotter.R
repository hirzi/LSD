
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
setwd("/path/")

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
# Cap maximum and minumum asymmetry, in case some values are extremely high.
asymmetry_min_max_cap_value <- 3

### Plotting main function
#for (Chr in seq(iChr,iChr)) {

#input_file <- paste0("genomeScan_table_Chr",iChr,suffix,".txt")
input_file <- "genomeScan_table_Chr6_2popModel.txt"
output_suffix <- "2popModel"

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

# Let's transform asymmetry values to emphasize extreme values, here via log(a/(1-a)), known as the odds ratio of the posterior
genomeScan_table[,ncol(genomeScan_table)+1] <- log(genomeScan_table$asymmetry / (1 - genomeScan_table$asymmetry))
colnames(genomeScan_table)[ncol(genomeScan_table)-1] <- "asymmetry_original"
colnames(genomeScan_table)[ncol(genomeScan_table)] <- "asymmetry"

# To color significant values by asymmetry
# Let's cap (new) asymmetry values above 3 to 3. Recall log(a/(1-a)) = 3 corresponds to a = 0.95.
cap_min_max_asymmetry <- function(x, asymmetry_min_max_cap){
  if(x > asymmetry_min_max_cap){
    x <- asymmetry_min_max_cap
    return(x)
  }
  else if(x < -asymmetry_min_max_cap){
    x <- -asymmetry_min_max_cap
    return(x)
  }
  else {
    return(x)
  }
}
genomeScan_table$asymmetry <- sapply(genomeScan_table$asymmetry, cap_min_max_asymmetry, asymmetry_min_max_cap=asymmetry_min_max_cap_value)
genomeScan_table[,ncol(genomeScan_table)+1] <- (genomeScan_table[,11]== 1)*genomeScan_table$asymmetry
genomeScan_table[,ncol(genomeScan_table)][genomeScan_table[,ncol(genomeScan_table)] == 0] <- NA
# Define number of colour breaks and colour palette
# Option 1: Equal sized bins across colour scale
# col_breaks <- 11
# ##genomeScan_table[,ncol(genomeScan_table)+1] <- as.numeric(cut(genomeScan_table[,15],breaks = col_breaks)) # breaks with min/max value = min/max of data
# genomeScan_table[,ncol(genomeScan_table)+1] <- round(genomeScan_table[,ncol(genomeScan_table)]*(5/3)) + 6  # breaks with min/max value = 0/1 (fixed interval)
# col_pal <- colorRampPalette(c('red3','white','dodgerblue3'))
# genomeScan_table[,ncol(genomeScan_table)+1] <- col_pal(col_breaks)[genomeScan_table[,ncol(genomeScan_table)]]
# Option 2: Equal sized bins across colour scale, apart from the central bin around asymmetry=0, which is larger
genomeScan_table[,ncol(genomeScan_table)+1] <- round(genomeScan_table[,ncol(genomeScan_table)]*(14/3)) + 15  # breaks with min/max value = 0/1 (fixed interval)
cols_asymmetry <- c("#CD0000", "#D01212", "#D42424", "#D73636", "#DB4848", "#DE5B5B", "#E26D6D", "#E67F7F", "#E99191", "#EDA3A3", "#F0B6B6", "#F4C8C8", "#F7DADA", "#FFFFFF", "#FFFFFF", "#FFFFFF", "#DEEBF7", "#CDE1F4", "#BDD7F0", "#ACCDED", "#9CC3E9", "#8BB9E6", "#7BAFE2", "#6AA5DE", "#5A9BDB", "#4991D7", "#3987D4", "#287DD0", "#1874CD")
genomeScan_table[,ncol(genomeScan_table)+1] <- cols_asymmetry[genomeScan_table[,ncol(genomeScan_table)]]
# Continue
genomeScan_table[,ncol(genomeScan_table)+1] <- genomeScan_table[,ncol(genomeScan_table)]
genomeScan_table[,ncol(genomeScan_table)][!is.na(genomeScan_table[,ncol(genomeScan_table)])] <- "black"

# Let's add position (take midpoint of window)
genomeScan_table$Chr_position <- (genomeScan_table$window_start + genomeScan_table$window_end) / 2

# Manhatten plot based on joint 2D credible intervals
par(mfrow=c(1,1),mar=c(5,5,5,5), bg = 'white')

prefix_joint <- str_split(gsub('.{1}$', '', substring(colnames(genomeScan_table)[8],16)),"_X_")[[1]]
#plot(seq(1,nrow(genomeScan_table)), genomeScan_table[,8], xlab="Position", ylab="-log10(1-P)", main = paste0("Joint 2D ", prefix_joint[1], " - ", prefix_joint[2], " significance; ", output_suffix), pch=19, cex = 0.3)
par(family = "ArialMT")
#plot(genomeScan_table[,3], genomeScan_table[,8], xlab="", ylab="-log10(1-P)", main = output_suffix, pch=19, cex = 0.3, cex.main=1, cex.lab=1, cex.axis=1)
plot(genomeScan_table[,3], genomeScan_table[,8], xlab="Position", ylab="-log10(1-P)", main = output_suffix, pch=19, cex = 0.3, cex.main=1, cex.lab=1, cex.axis=1, col = "grey50")
points(genomeScan_table[,3], genomeScan_table[,8], pch=21, cex = 1.25, col=genomeScan_table[,ncol(genomeScan_table)], bg=genomeScan_table[,ncol(genomeScan_table)-1])
abline(h = -log10(1-credible_interval), lty = 2, lwd = 1.5, col = "grey30")
#ColorLegend("right", col=col_pal(col_breaks*5), labels=sprintf("%.1f",seq(-3,3,0.5)), adj = c(0,0.5), cntrlbl = TRUE, cex=0.5, inset = -0.04, width = round(nrow(genomeScan_table)*40), height = min(likelihood_max_cap_2D,max(genomeScan_table[,8]))*1.08)
ColorLegend("right", col=cols_asymmetry, labels=sprintf("%.1f",seq(-3,3,0.5)), adj = c(0,0.5), cntrlbl = TRUE, cex=0.5, inset = -0.04, width = round(nrow(genomeScan_table)*40), height = min(likelihood_max_cap_2D,max(genomeScan_table[,8]))*1.08)
#ColorLegend(x=54800000, y=8, col=col_pal(col_breaks), labels=sprintf("%.1f",seq(0,1,0.1)), cntrlbl = TRUE, cex=0.5, inset = -0.05, width = round(nrow(genomeScan_table)*40), height = 10)

# # Manhatten plot based on 1D credible intervals
# prefix_marginals <- c(gsub('.{1}$', '', substring(colnames(genomeScan_table)[6],10)), gsub('.{1}$', '', substring(colnames(genomeScan_table)[7],10)))
# par(mfrow=c(1,1),mar=c(4,4,4,4))
# for (p in 1:2) {
#   # We want to distinguish between windows that satisfy a single parameter CI signficance and those that satisfy the joint (2) parameter CI significance.
#   # We do this by adding a third value (3), which signifies windows that satisfy the 1D CI but not the joint 2D CI
#   signficance_combined_vector <- vector()
#   for (i in seq(1, nrow(genomeScan_table))) {
#     if ( (genomeScan_table[,p+8][i] == 1) & (genomeScan_table[,11][i] == 2)) {
#       signficance_combined_vector[i] <- 3
#     }
#     else {
#       signficance_combined_vector[i] <- genomeScan_table[,11][i]
#     }
#   }
#   # And plot
#   plot(seq(1,nrow(genomeScan_table)), genomeScan_table[,p+5], xlab="Position", ylab="-log10(1-P)", main = paste0(prefix_marginals[p], " (marginal) significance; ", output_suffix), pch=19, cex = 0.3, col=c("red3", "grey35","darkorange")[signficance_combined_vector])
#   abline(h = -log10(1-credible_interval_1D), lty = 2, lwd = 1.5, col = "grey30")
#   legend("topleft", legend = c(paste0("Not significant(<",credible_interval*100,"^0.5% CI)"), paste0("Significant for this marginal parameter (>", credible_interval*100,"^0.5% CI)"), paste0("Significant for (both) joint parameters (>", credible_interval*100,"% CI)")), col=c("grey35","darkorange", "red3"), pch = 19, bty = "n", pt.cex = 1.25, cex = 0.6, horiz = FALSE, inset = c(0.01, 0.025))
# }
#}
