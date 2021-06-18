## from https://github.com/anbena/ABC-FDSS/blob/master/2pop/div_and_mig.r
## publication: https://www.biorxiv.org/content/10.1101/507897v2.full

#wrapper script to simulate two population split with asymmetrical gene flow from divergence to the present
#Rscript --vanilla ./div_and_mig.r <nchr> <locuslength> <nlociXsim> <rec.rate=1e-8> <nsimulations>
args<-commandArgs(trailingOnly=TRUE)
ms<-"/opt/software/genetics/ms/ms"
cpd<-"./"
mod<-"m2"
nchr<-as.character(args[1])
tgen<-1
mu<-1e-8
recomb<-as.numeric(args[4])#recombination rate
ll<-as.numeric(args[2])#locus length
nsims<-as.numeric(args[5])#number of ABC simulations
nloci<-as.numeric(args[3])#loci to simulate in each sim
out<-paste(mod,"_ll",as.character(ll),"_nl",as.character(nloci),"_r",as.character(recomb),"_nc",nchr,sep="")

#Prior distributions parameters
n1<-sample(500:50000,nsims,replace=T)
n2<-sample(500:50000,nsims,replace=T)
na<-sample(500:50000,nsims,replace=T)
td<-sample(300:20000,nsims,replace=T)
m12<-rexp(nsims,0.1)*4
m21<-rexp(nsims,0.1)*4

#Parameters transformation for ms simulator
tn1<-4*n1*mu*ll
rn2<-n2*4*mu*ll/tn1
rna<-na*4*mu*ll/tn1
tds<-(td/tgen)/(4*n1)
srec<-4*n1*(recomb*(ll-1))

#Output: parameters files
partable<-cbind(n1,n2,na,td,m12,m21)
colnames(partable)<-c("n1","n2","na","td","m12","m21")
partablescaled<-cbind(tn1,rn2,tds,tds,rna,m12,m21,srec)
write.table(partable,paste(out,".param",sep=""),row.names=F,quote=F,sep="\t")
write.table(partablescaled,paste(out,".paramscaled",sep=""),row.names=F,col.names=T,quote=F,sep="\t")

#Output summary statistics: FDSS
i<-1
write(paste(ms," ",as.character(2*as.numeric(nchr))," ",as.character(nloci)," -t ",as.character(tn1[i])," -r ",as.character(srec[i])," ",as.character(ll)," -I 2 ",nchr," ",nchr," -n 2 ",as.character(rn2[i])," -m 1 2 ",as.character(m12[i])," -m 2 1 ",as.character(m21[i])," -ej ",as.character(tds[i])," 2 1 -en ",as.character(tds[i])," 1 ",as.character(rna[i]),"| ",cpd,"compute_ss.py -np 2 -nc ",nchr," -w 30 -b 50 -s > ",out,".tab",sep=""), stderr())
for (i in 1:nsims){
	print(i)
	if (i==1){
		system(paste(ms," ",as.character(2*as.numeric(nchr))," ",as.character(nloci)," -t ",as.character(tn1[i])," -r ",as.character(srec[i])," ",as.character(ll)," -I 2 ",nchr," ",nchr," -n 2 ",as.character(rn2[i])," -m 1 2 ",as.character(m12[i])," -m 2 1 ",as.character(m21[i])," -ej ",as.character(tds[i])," 2 1 -en ",as.character(tds[i])," 1 ",as.character(rna[i]),"| ",cpd,"compute_ss.py -np 2 -nc ",nchr," -w 30 -b 50 -s > ",out,".tab",sep=""))
	}
	else{
		system(paste(ms," ",as.character(2*as.numeric(nchr))," ",as.character(nloci)," -t ",as.character(tn1[i])," -r ",as.character(srec[i])," ",as.character(ll)," -I 2 ",nchr," ",nchr," -n 2 ",as.character(rn2[i])," -m 1 2 ",as.character(m12[i])," -m 2 1 ",as.character(m21[i])," -ej ",as.character(tds[i])," 2 1 -en ",as.character(tds[i])," 1 ",as.character(rna[i]),"| ",cpd,"compute_ss.py -np 2 -nc ",nchr," -w 30 -b 50 -s >> ",out,".tab",sep=""))
	}
}
