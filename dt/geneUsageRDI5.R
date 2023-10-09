# Load required packages
library(alakazam)
library(dplyr)
library(scales)
library(shazam)
library(tidyverse)
library(ggthemes)
library(reshape)
library(ggplot2)
library(tidyr)
library(rdi)

infoTable<-read.table("uniqueSamples.txt",sep="\t",header=FALSE)
colnames(infoTable)<-c("patientName")
#tiggerGenotype<-read.table("tiggerOutputPatient/geno_bayesian_104.txt",sep="\t",header=TRUE)

clusterTable<-read.table("clusters.txt",sep="\t",header=TRUE)

clusterTable$patientName<-as.integer(sub(".*-", "", clusterTable$sample.id))
clusterTable$patientName2<-as.integer(sub(".*-", "", clusterTable$sample.id))

infoTable<-merge(infoTable,clusterTable,by="patientName")

changeoTableMaster<-read.table("changeoTableMasterClonePatient_expandedStatus_SHM.tsv",sep="\t",header=TRUE)

infoTable$patient<-infoTable$patientName
changeoTableMaster<-merge(changeoTableMaster,infoTable,by="patient")


time<-"A"
isotype<-"IgM"
#isotype<-"All"
gene<-"v"

changeoTableMasterIsotype<-changeoTableMaster[which(changeoTableMaster$cprimer==isotype & changeoTableMaster$time==time),]
#changeoTableMasterIsotype<-changeoTableMaster
changeoTableMasterIsotype$snp2<-paste0(changeoTableMasterIsotype$pop,"_",changeoTableMasterIsotype$snp_gt_x)

#RDI

dir.create("plots")
folder<-paste0("rdi_",gene)
dir.create(paste0("plots/",folder))


#subset to gene clicque genes (SKIP IF YOU WANT ALL GENES)
# Sample character vector
text <- changeoTableMasterIsotype[,paste0("germline_",gene,"_call")]
# Define multiple patterns
patterns <- c("IGHV3-43","IGHV3-53","IGHV1-69","IGHV3-66","IGHV4-61","IGHV4-59","IGHV3-64","IGHV4-31")
# Use grep to search for patterns
matching_indices <- grep(paste(patterns, collapse = "|"), text)
changeoTableMasterIsotype<-changeoTableMasterIsotype[matching_indices,]


# RDI Global
changeo<-changeoTableMasterIsotype
genes<-changeo[[paste0("germline_",gene,"_call")]]
seqAnnot <- paste0(changeo$sample,"_",changeo$snp2)
cts<-calcVDJcounts(genes, seqAnnot)
d <- calcRDI(cts)
df  <- melt(as.matrix(d))
p    <- t(apply(df[,c(1,2)],1,FUN=sort))
rmv1 <- which(p[,1] == p[,2])
p    <- paste(p[,1],p[,2],sep="|")
rmv2 <- which(duplicated(p))
df   <- df[-c(rmv1,rmv2),] #removing self distances and duplicated distances
RDI<-df

RDI$sampleA<-sapply(as.character(RDI$X1), function(x) unlist(strsplit(x, "_", fixed = TRUE))[1])
RDI$genotype1A<-sapply(as.character(RDI$X1), function(x) unlist(strsplit(x, "_", fixed = TRUE))[2])
RDI$genotype2A<-sapply(as.character(RDI$X1), function(x) unlist(strsplit(x, "_", fixed = TRUE))[3])
RDI$genotypeA<-paste0(RDI$genotype1A,"_",RDI$genotype2A)

RDI$sampleB<-sapply(as.character(RDI$X2), function(x) unlist(strsplit(x, "_", fixed = TRUE))[1])
RDI$genotype1B<-sapply(as.character(RDI$X2), function(x) unlist(strsplit(x, "_", fixed = TRUE))[2])
RDI$genotype2B<-sapply(as.character(RDI$X2), function(x) unlist(strsplit(x, "_", fixed = TRUE))[3])
RDI$genotypeB<-paste0(RDI$genotype1B,"_",RDI$genotype2B)


#RDI<-as.vector(d)
#RDI<-as.data.frame(RDI)
#colnames(RDI)<-c("RDI","group")
colnames(RDI)<-c("X1","X2","RDI","sampleA","genotype1A","genotype2A","genotypeA","sampleB","genotype1B","genotype2B","genotypeB")


write.table(RDI,paste0("plots/",folder,"/RDI_",isotype,"_",gene,"_Global.txt"),sep="\t",row.names = FALSE,quote=FALSE)
#write.table(RDI,paste0("plots/",folder,"/RDI_",isotype,"_",gene,"_Global_CliqueGenes.txt"),sep="\t",row.names = FALSE,quote=FALSE)