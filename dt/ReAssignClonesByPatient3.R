library(shazam)
library(ggplot2)

args <- commandArgs(trailingOnly = TRUE)
# read input file
patientName <- args[1]
#patientName<-"B"
#patientName

infoTable<-read.table("/home/datier01/bCellDevelopment/sampleInfo.txt",sep="\t",header=FALSE)
colnames(infoTable)<-c("samplePath","sampleName","patient_time")
infoTable$patientName = substr(infoTable$patient_time, 1, nchar(infoTable$patient_time)-1)
infoTable$time = substr(infoTable$patient_time, nchar(infoTable$patient_time), nchar(infoTable$patient_time))

infoTable$sample<-substr(infoTable$sampleName, 1, nchar(infoTable$sampleName)-3)
infoTable$chain<-substr(infoTable$sampleName,  nchar(infoTable$sampleName)-1, nchar(infoTable$sampleName))

changeoInput<-"/home/datier01/bCellDevelopment/changeos_igblast_PreTigger_6Sept2023"
changeoOutput<-"/home/datier01/bCellDevelopment/changeoPatientCloneOutputs_PreTigger_6Sept2023"
dir.create(changeoOutput)

#IMPORT CHANGEO, WRITE TO FILE
infoTablePatient<-infoTable[which(infoTable$patientName==patientName & infoTable$chain=="HC"),]
changeoTableTemp<-list()
changeoTableTemp[[1]]<-read.table(paste0(changeoInput,"/",infoTablePatient[which(infoTablePatient$time=="1"),"sampleName"],"_atleast-2_db-pass_clone-pass_germ-pass.tsv"),sep="\t",header=TRUE)
changeoTableTemp[[1]]$sample<-paste0(patientName,"1")
changeoTableTemp[[1]]$patient<-patientName
changeoTableTemp[[1]]$time<-"1"
changeoTableTemp[[2]]<-read.table(paste0(changeoInput,"/",infoTablePatient[which(infoTablePatient$time=="2"),"sampleName"],"_atleast-2_db-pass_clone-pass_germ-pass.tsv"),sep="\t",header=TRUE)
changeoTableTemp[[2]]$sample<-paste0(patientName,"2")
changeoTableTemp[[2]]$patient<-patientName
changeoTableTemp[[2]]$time<-"2"
changeoTableTemp[[3]]<-read.table(paste0(changeoInput,"/",infoTablePatient[which(infoTablePatient$time=="3"),"sampleName"],"_atleast-2_db-pass_clone-pass_germ-pass.tsv"),sep="\t",header=TRUE)
changeoTableTemp[[3]]$sample<-paste0(patientName,"3")
changeoTableTemp[[3]]$patient<-patientName
changeoTableTemp[[3]]$time<-"3"
changeoTableTemp[[4]]<-read.table(paste0(changeoInput,"/",infoTablePatient[which(infoTablePatient$time=="4"),"sampleName"],"_atleast-2_db-pass_clone-pass_germ-pass.tsv"),sep="\t",header=TRUE)
changeoTableTemp[[4]]$sample<-paste0(patientName,"4")
changeoTableTemp[[4]]$patient<-patientName
changeoTableTemp[[4]]$time<-"4"
changeoTable<-do.call(rbind, changeoTableTemp)
write.table(changeoTable,paste0(changeoOutput,"/changeoTable",patientName,".tsv"),sep="\t",row.names = FALSE)
#Get the threshold
dist_ham <- distToNearest(changeoTable, sequenceColumn="junction", vCallColumn="v_call", jCallColumn="j_call",model="ham", normalize="len", nproc=1)
output <- findThreshold(dist_ham$dist_nearest, method="density")
threshold <- output@threshold
print(threshold)
#Run define clones
#system("DefineClones.py -d changeoTableTemp.tsv --act set --model ham --norm len --dist 0.16", intern = FALSE,ignore.stdout = FALSE)
system(paste0("DefineClones.py -d ",changeoOutput,"/changeoTable",patientName,".tsv --outdir ",changeoOutput," --act set --model ham --norm len --dist ",threshold), intern = FALSE,ignore.stdout = FALSE)
