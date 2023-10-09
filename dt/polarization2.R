# Load required packages
library(alakazam)
library(dplyr)
library(scales)
library(shazam)
library(tidyverse)
library(ggthemes)
library(reshape)
library(ggplot2)
library(treemap)
library(ggplot2)
library(ggrepel)

# Subset example data
#data(ExampleDb)
#ExampleDb

changoTableMaster<-read.table("changeoTableMasterClonePatient_expandedStatus.tsv",sep="\t",header=TRUE)

dir.create("clones")
plotFolder<-"clones"

#clones<-countClones(changoTableMaster,group=c("patient", "time"), clone="clone_id")

clones<-countClones(changoTableMaster,group=c("sample"), clone="clone_id")


#Polarization
polarizationList<-list()
topClonesList<-list()
for(ii in unique(clones$sample)){
  clonesSample<-clones[which(clones$sample==ii),]
  y<-0
  jj<-1
  while(y<0.2)  { 
    y<-y+as.numeric(clonesSample[jj,"seq_freq"])
    jj<-jj+1 
    #print(jj)
  }
  polarization<-jj-1
  print(ii)
  print(polarization)
  topClonesList[[ii]]<-clonesSample[1:polarization,]
  polarization<-polarization/nrow(clonesSample)
  print(nrow(clonesSample))
  print(polarization)
  polarizationList[[ii]]<-c(ii,polarization)
}

df<-as.data.frame(do.call(rbind, polarizationList))
colnames(df)<-c("sample","polarization")

write.table(df,"clones/polarizationSample.txt",sep="\t",row.names = FALSE, quote=FALSE)

sample_curve <- alphaDiversity(changoTableMaster, group="sample", clone="clone_id",
                               min_q=2, max_q=2, step_q=1,
                               ci=0.95, nboot=100)

diversity<-sample_curve@diversity 

write.table(diversity[which(diversity$q==2),],"clones/diversityQ2Sample.txt",sep="\t",row.names = FALSE, quote=FALSE)

topClonesTumor<-as.data.frame(do.call(rbind, topClonesList))

topClonesList<-list()
for(ii in unique(clonesTissue$patient_id)){
  clonesTissuePatient<-clonesTissue[which(clonesTissue$patient_id==ii),]
  #clonesTissue<-clonesTissue[order(-clonesTissue$seq_freq),]
  print(ii)

  topClonesList[[ii]]<-clonesTissuePatient[1:5,]
}


changoTableMaster2

changoList<-list()
for(ii in unique(changoTableMaster2$patient_id))
{
  changoList[[ii]]<-changoTableMaster2[which(changoTableMaster2$patient_id==ii),]
  changoList[[ii]]
}

changoList<-list()
topClonesList<-list()
for(ii in unique(clonesTissue$patient_id)){
  clonesTissuePatient<-clonesTissue[which(clonesTissue$patient_id==ii),]
  #clonesTissue<-clonesTissue[order(-clonesTissue$seq_freq),]
  changoTissuePatient<-changoTableMaster2[which(changoTableMaster2$patient_id==ii),]
  #print(ii)
  topClones<-clonesTissuePatient[1:3,]
  topClonesList[[ii]]<-topClones
  changoList[[ii]]<-changoTissuePatient[which(changoTissuePatient$clone_id%in%topClones$clone_id),]
}
topClonesTumor<-as.data.frame(do.call(rbind, topClonesList))

changoTableMaster3<-as.data.frame(do.call(rbind, changoList))

#dfNormal<-as.data.frame(do.call(rbind, polarizationList))
#colnames(dfNormal)<-c("patient_id","polarization","tissue")

#dfPBMC<-as.data.frame(do.call(rbind, polarizationList))
#colnames(dfPBMC)<-c("patient_id","polarization","tissue")

#df<-rbind(dfTumor,dfNormal,dfPBMC)

#df<-dfTumor
#df$polarization<-as.numeric(df$polarization)
df<-topClonesTumor

#dfOrig<-df

pdf(paste0(plotFolder,"/clones/polarization_cloneSample.pdf"))
  ggplot(df,aes( tissue, polarization, fill=tissue)) +
  geom_boxplot() +
  geom_line(aes(group=patient_id), position = position_dodge(0.2)) +
  geom_point(aes(fill=tissue,group=patient_id), position = position_dodge(0.2)) +
  theme(legend.position = "none")
dev.off()


#groups<-read.table("../cytofCancer/plotsCD45Patients1_48_6Nov2020_Corrected_1234_StateCluster_25_Kamille105/patientBGroup27April2022.txt",sep="\t",header=TRUE)
#groups<-read.table("../cytofCancer/plotsCD45Patients1_48_6Nov2020_Corrected_1234_StateCluster_25_Kamille118/kmeansAssignmentsAllDiffInTumor2.txt",sep="\t",header=TRUE)
groups<-read.table("../cytofCancer/plotsCD45Patients1_48_6Nov2020_Corrected_1234_StateCluster_25_Kamille118_2/patentGroupKMeans_8Nov2022.txt",sep="\t",header=TRUE)
groups$patient_id<-groups$X
groups$X<-NULL
groups$patient_id<-paste0("Patient",groups$patient_id)
groups$kMeansGroup<-as.factor(groups$kMeansGroup)

#NEEDS ALL=TRUE OR PATIENTS WITHOUT GROUP (WITHOUT TUMOR SAMPLE) WILL BE EXCLUDED
df2<-merge(df,groups,by="patient_id", all=TRUE)
#df2<-df2[!is.na(df2$polarization),]
df2<-df2[!is.na(df2$seq_freq),]
df2<-df2[!is.na(df2$kMeansGroup),]
#df2<-df2[!is.na(df2$AllDiffInTumor),]
#df2$AllDiffInTumor<-as.factor(df2$AllDiffInTumor)

library(ggplot2)
library(ggrepel)


#df3<-df2[which(df2$tissue=="Tumor"),]
pdf(paste0(plotFolder,"/clones/polarization.pdf"))
ggplot(df2, aes(x=kMeansGroup, y=polarization, fill=kMeansGroup)) +
geom_boxplot()+
geom_point()+
#geom_text_repel(aes(label = patient_id)) + 
ylim(0.1, 0.7)
#geom_text_repel(aes(label = df$patient_id)) 
dev.off()

