#!/usr/bin/env Rscript

###################################################
# Generates 4 graphs, based on 4 statistics files #
# tss.stats.csv, exon.stats.csv, intron.stats.csv #
# and tss.distance.csv.                           #
# Reads statitics files from graphs folder.       #
# Stores graphs in same graphs folder.            #
#                                                 #
# Author: Alexis Blanchet-Cohen                   #
# Modified from script written by Maxime Caron    #
###################################################

#########################################################
# TSS categories statistics (location of binding sites) #
#########################################################

d1<-read.table("graphs/statistics/tss.stats.csv", header=T, sep=",", check.names=F)
slices <- c(d1[,1] + d1[,2], d1[,3], d1[,4], d1[,5], d1[,6], d1[,7])
lbls <- c("gene", names(d1[3:length(d1)]))
pct <- round(slices/sum(slices)*100)
lbls <- paste(lbls, "(",pct, sep="") # add percents to labels
lbls <- paste(lbls,"%)",sep="") # add % to labels 

# PDF file
pdf("graphs/pdf/location_binding_sites.pdf")
	pie(slices,labels=lbls, main=paste("Location analysis of binding sites"))
dev.off()

# This section is to change the location of the labels.
# pdf("graphs/pdf/location_binding_sites2.pdf")
# 	labels <- lbls
# 	lbls[5] <- " "
# 	pie(slices, main=paste("Location analysis of binding sites"), labels=lbls)
# 	text(0.1,-0.9, "gene_desert(2%)")
# dev.off()

# postscript(file ="graphs/eps/location_binding_sites2.eps", onefile=FALSE, horizontal=FALSE)
# 	pie(slices, main=paste("Location analysis of binding sites"), labels=lbls)
# 	text(0.1,-0.9, "gene_desert(2%)")
# dev.off()
       		
###################################################
# Disitribution of peaks within introns and exons #
###################################################

exons<-paste("graphs/statistics/exon.stats.csv")
if(file.info(exons)$size == 0) {
	d1=data.frame(c(0))
} else {
	d1<-read.table(exons, header=F, sep=",", check.names=F)
}
	introns<- paste("graphs/statistics/intron.stats.csv")
if(file.info(introns)$size == 0) {
	d1=data.frame(c(0))
} else {
	d2<-read.table(introns, header=F, sep=",", check.names=F)
}	

# PDF files
pdf("graphs/pdf/distribution_peaks_exons.pdf")
	hist(d1[,1], breaks=length(levels(as.factor(d1[,1]))), xlim=c(0,20), main=paste("Distribution of peaks found within exons"), xlab="Exon", ylab="Number of peaks", col="lightblue")
dev.off()

pdf("graphs/pdf/distribution_peaks_introns.pdf")
	hist(d2[,1], breaks=length(levels(as.factor(d2[,1]))), xlim=c(0,20), main=paste("Distribution of peaks found within introns"), xlab="Intron", ylab="Number of peaks", col="lightblue")
dev.off()


##############################################		
# Distribution of distance of peaks from TSS #
##############################################
distance<-paste("graphs/statistics/tss.distance.csv")
d1<-read.table(distance, header=F, sep=",", check.names=F)
d1<-subset(d1, d1[,1] > -10000 & d1[,1] < 10000)

# PDF file
pdf("graphs/pdf/distribution_peaks_distance_tss.pdf")
	hist(d1[,1], breaks=seq(-10000,10000,1000), main=paste("Distribution of peak distances relative to TSS"), xlab="Distance to TSS (bp)", ylab="Number of peaks", col="lightblue")		
dev.off()
