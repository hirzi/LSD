#### This script outputs the PLS components of the simulation file in addition to the RMSEP plots which we can refer to when deciding the number of components to retain.
## Example (in terminal): Rscript find_pls.r 2pop_simpleModel_4params_simulatedSumStats.txt
## Lines which may need to be modified are line 20 (# of PLS components), line 22 (working directory), line 32 (positional (column) indicator of first summary statistics entry) and line 36 (positional (column) indicator of free parameters)

# First, we install libraries if needed and load them
need.pkg <- c("MASS", "pls")
inst.pkg <- which(!need.pkg %in% installed.packages())
if (length(inst.pkg) > 0) {cat("\nyou lack some R packages, installing...\n") ; install.packages(need.pkg[inst.pkg], repos="http://cran.rstudio.com/")}
cat("\nloaded packages and versions:\n")
for (i in need.pkg) {
  suppressWarnings(suppressPackageStartupMessages(do.call(library, list(i))))
  cat(i)
  cat(paste0(" v", packageVersion(i)))
  cat("\n")
  rm(i)
}
cat("\n")

#open File
numComp<-5;

directory<-"/path/";

filename<-commandArgs()[length(commandArgs())];
#filename<-"2pop_simpleModel_4params_simulatedSumStats.txt"
print(paste("Reading file '", filename, "'", sep=""));

#read input (summary statistics) file
a<-read.table(paste(directory, filename, sep=""), header=T);

# Define the starting column for the summary statistics
firstStat<-13;
stats<-a[,firstStat:length(a[1,])];

# Define the columns holding the free parameters
p<-c(2,3,4,5);
if(length(p)==1){
  params <- data.frame(x=a[,p]); names(params)<-names(a)[p];
} else {
  params<-a[,p]; 
}
print(names(params));

#standardize the params
for(i in 1:length(params)){params[,i]<-(params[,i]-mean(params[,i]))/sd(params[,i]);}

#force stats in [1,2]
myMax<-c(); myMin<-c(); lambda<-c(); myGM<-c();
for(i in 1:length(stats)){
  myMax<-c(myMax, max(stats[,i]));
  myMin<-c(myMin, min(stats[,i]));
  stats[,i]<-1+(stats[,i]-myMin[i])/(myMax[i]-myMin[i]);
}

#transform statistics via boxcox  
library("MASS");	
for(i in 1:length(stats)){		  
  d<-cbind(stats[,i], params);
  mylm<-lm(as.formula(d), data=d)			
  myboxcox<-boxcox(mylm, lambda=seq(-50, 80, 1/10), plotit=T, interp=T, eps=1/50);	
  lambda<-c(lambda, myboxcox$x[myboxcox$y==max(myboxcox$y)]);			
  print(paste(names(stats)[i], myboxcox$x[myboxcox$y==max(myboxcox$y)]));
  myGM<-c(myGM, exp(mean(log(stats[,i]))));			
}

#standardize the BC-stats
myBCMeans<-c(); myBCSDs<-c();
for(i in 1:length(stats)){
  stats[,i]<-(stats[,i]^lambda[i] - 1)/(lambda[i]*myGM[i]^(lambda[i]-1));	
  myBCSDs<-c(myBCSDs, sd(stats[,i]));
  myBCMeans<-c(myBCMeans, mean(stats[,i]));		
  stats[,i]<-(stats[,i]-myBCMeans[i])/myBCSDs[i];
}

#perform pls
library("pls");
#myPlsr<-plsr(as.matrix(params) ~ as.matrix(stats), scale=F, ncomp=numComp, validation="LOO");
myPlsr<-plsr(as.matrix(params) ~ as.matrix(stats), scale=F, ncomp=numComp);

#write pls to a file
myPlsrDataFrame<-data.frame(comp1=myPlsr$loadings[,1]);
for(i in 2:numComp) { myPlsrDataFrame<-cbind(myPlsrDataFrame, myPlsr$loadings[,i]); } 
write.table(cbind(names(stats), myMax, myMin, lambda, myGM, myBCMeans, myBCSDs, myPlsrDataFrame), file=paste(directory, "PLSdef_", filename, sep=""), col.names=F, row.names=F, sep="\t", quote=F);

#make RMSE plot
pdf(paste(directory, "RMSE_", filename, ".pdf", sep=""));
plot(RMSEP(myPlsr), ylim=c(0,1));
dev.off();
