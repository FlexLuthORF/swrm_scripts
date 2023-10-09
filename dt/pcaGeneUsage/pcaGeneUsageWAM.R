library(alakazam)
library(tidyr)

changeoTableMaster<-read.table("changeoTableMasterClonePatient_expandedStatus_SHM.tsv",sep="\t",header=TRUE)

infoTable<-read.table("uniqueSamples.txt",sep="\t",header=FALSE)
colnames(infoTable)<-c("patientName")
#tiggerGenotype<-read.table("tiggerOutputPatient/geno_bayesian_104.txt",sep="\t",header=TRUE)

clusterTable<-read.table("clusters.txt",sep="\t",header=TRUE)

clusterTable$patientName<-as.integer(sub(".*-", "", clusterTable$sample.id))
clusterTable$patientName2<-as.integer(sub(".*-", "", clusterTable$sample.id))

infoTable<-merge(infoTable,clusterTable,by="patientName")

infoTable$patient<-infoTable$patientName
changeoTableMaster<-merge(changeoTableMaster,infoTable,by="patient")

isotype<-"IgM"
#isotype<-"All"
gene<-"v"
time<-"A"

#outfolder<-paste0("geneTables_",time)
outfolder<-paste0("cloneWithgeneTables_",time)
dir.create(outfolder)

changeoTableMasterIsotype<-changeoTableMaster[which(changeoTableMaster$cprimer==isotype & changeoTableMaster$time==time),]
changeoTableMasterIsotype$snp2<-paste0(changeoTableMasterIsotype$pop,"_",changeoTableMasterIsotype$snp_gt_x)
#geneCount <- countGenes(changeoTableMasterIsotype, gene=paste0("germline_",gene,"_call"), groups=c("sample","snp2"),  mode="gene", fill=FALSE)
geneCount <- countGenes(changeoTableMasterIsotype, gene=paste0("germline_",gene,"_call"), groups=c("sample","snp2"), clone="clone_id", mode="gene", fill=FALSE)
#geneCount<-geneCount[-grep("*_NA",geneCount$snp2),]
write.table(geneCount,paste0(outfolder,"/",gene,"_",isotype,"_genes_",time,".tsv"),sep="\t",row.names = FALSE, quote=FALSE)




#geneCount<-read.table(paste0("/scratch/scratch/datier01/CW56_Flu/cloneWithgeneTables_",time,"/",gene,"_",isotype,"_genes_",time,".tsv"),header=TRUE)
library(tidyr)
#print(vGene)
geneCount<-as.data.frame(geneCount)
geneCount$sample_snp2<-paste0(geneCount$sample,":",geneCount$snp2)
#geneCount<-geneCount[-grep("*_NA",geneCount$snp2),]
#geneCount$snp2<-as.factor(geneCount$snp2)
#m<-spread(geneCount[,c(3,5,6)], gene, seq_freq, fill=0)
m<-spread(geneCount[,c(3,5,6)], gene, clone_freq, fill=0)
m_sample_snp2<-m$sample_snp2
m$sample_snp2<-NULL
m_sample<-sub(paste0("\\:", ".*"), "", m_sample_snp2)
m_snp2<-sub(paste0(".*", "\\:"), "", m_sample_snp2)
m_snp<-sapply(strsplit(m_snp2, "_"), "[[", 1)

m<-as.data.frame(m)
rownames(m)<-m_sample
m<-m[,which((colSums(m == 0) < 0.5*dim(m)[1]))]

#m<-m[which(m_snp2%in%c("F_A/A","F_A/G","F_G/G")),]
m<-m[which(m_snp%in%c("F","L")),]

library(factoextra)
res.pca <- prcomp(m, scale = FALSE)
folder<-"pca"
pdf(paste0(folder,"/",gene,"_genes_",isotype,"_PCA_V_clonesWithGenes.pdf"))
p<-fviz_pca_ind(res.pca, geom = "point",
             pointsize = 2,
             pointshape = 19,
             col.ind = m_snp[which(m_snp%in%c("F","L"))],  # Individuals color
             #palette = c("blue", "red"),
             addEllipses = FALSE, # Concentration ellipses
             ellipse.type = "confidence",
             ellipse.level=0.95
             #ggtheme=theme(axis.text=element_text(size=15),
             #axis.title=element_text(size=16, face="bold"), geom_point(size =3))
             #repel = TRUE     # Avoid text overlapping
             ) #+ ggtitle(vGene)
print(p)
#pdf(paste0("pca_d_",vGene,"loadings.pdf"))
q<-fviz_pca_var(res.pca,
            col.var = "contrib", # Color by contributions to the PC
            #gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
            repel = TRUE     # Avoid text overlapping
            )
print(q)
#dev.off()

#}
dev.off()