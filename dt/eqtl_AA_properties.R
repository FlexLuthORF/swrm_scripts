library(alakazam)
library(ggplot)
infoTable<-read.table("changeoPathsM.txt",sep="\t",header=FALSE)

extract_sample_name <- function(text) {
  parts <- strsplit(text, "/")[[1]]
  last_part <- tail(parts, 1)
  return(substr(last_part, 1, nchar(last_part) - 15))
}
infoTable$V2 <- sapply(infoTable$V1, extract_sample_name)
infoTable$V2<- gsub("-", ".", infoTable$V2)
infoTable$V2 

colnames(infoTable)<-c("changeoPath","sampleName")

changeoTableList<-list()
for(ii in 1:length(infoTable$sampleName)){
#for(ii in 1:5){
  tryCatch({
  print(infoTable$sampleName[ii])
  changeoTable<-read.table(infoTable$changeoPath[ii],sep="\t",header=TRUE)
  #Needed because 
  changeoTable$sample<-infoTable$sampleName[ii]
  changeoTableList[[ii]]<-changeoTable
  }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
}

changeoTableMaster<-do.call(rbind, changeoTableList)

#CDR3
db_props <- aminoAcidProperties(changeoTableMaster, seq="junction", trim=TRUE, label="cdr3")
colnames(db_props)

db_props2<-db_props[,c("sample","cdr3_aa_length","cdr3_aa_gravy","cdr3_aa_bulk","cdr3_aa_aliphatic","cdr3_aa_polarity","cdr3_aa_charge","cdr3_aa_basic","cdr3_aa_acidic","cdr3_aa_aromatic")]

# Split the dataframe by SampleID
sample_groups <- split(db_props2[, 2:(length(db_props2))], db_props2$sample)

# Calculate average of Value1 and Value2 for each sample
average_values <- lapply(sample_groups, function(sample_data) {
  colMeans(sample_data)
})

# Convert the list of averages to a dataframe
average_df <- data.frame(do.call(rbind, average_values))
du2<-average_df
IDs <- rownames(du2)

dir.create("eqtl_AA")
dir.create("eqtl_AA/AA_Properties")
#######
#Regression ANalysis
### IGHD locus SNVs: association with IGHD gene usage (all genes)
genes = colnames(du2)
genos = read.table("snps.txt", header=T, row.names=1, sep="\t")
g = genos[,IDs]

uni = NULL
for(i in rownames(g)){
  tt = as.data.frame(length(unique(na.omit(unlist(g[i,])))))
  colnames(tt) = "Uni"
  uni = rbind(uni, tt)
}

g$Uni = uni
g$NAs = rowSums(is.na(g))

g2 = subset(g, NAs < 5)
g2 = subset(g2, Uni > 1)
g2 = g2[,1:(length(g2)-2)]
du2 = du2[colnames(g2),]


interims = NULL
pvals = NULL

for (i in rownames(g2)){
  snp = as.data.frame(t(g2[i,]))

  interims = data.frame(row.names=i)
    for (a in genes){
      str = as.data.frame(du2[,a]) ### Use if not including ethnicity
      rownames(str) = as.vector(colnames(g2))
      p = as.data.frame((summary(lm(str[,1] ~ unlist(snp[,i]), na.action = na.exclude))$coefficients[2,4]))
      rownames(p) = i
      colnames(p) = a
      interims = cbind(interims, p)
    }
  pvals = rbind(pvals, interims)
}
write.table(pvals, paste0("eqtl_AA/AA_Properties/reg.pvals_AAProps.txt"), quote=F, col.names=T, row.names=T, sep="\t")

lten = -log10(pvals)
position<-rownames(lten)

pdf(paste0("eqtl_AA/AA_Properties/manahattan.pos3.pdf"), paper = "USr")
for (i in colnames(lten)){
    cols = ifelse(lten[,i] > 2.3, "red",
            ifelse(lten[,i] > 1.3 & lten[,i] < 2.3, "orange", "black"))
  par(mar=c(12,4,12,2))
    plot(lten[,i] ~ position,
        cex = .8,
        pch = 19,
        col = cols,
        ylab = "-log10 P value",
        xlab = "chr7 position",
        main = i,
        ylim=c(0,12))
    options(scipen=8)
    abline(h=c(2,4,6,8,10,12), col="gray", lty=2)
}
dev.off()

topPvals<-pvals[order(pvals$cdr3_aa_gravy),]
topPvals<-topPvals[1:10,]

pdf("gravy_boxplots.pdf")
for(ii in 1:length(topPvals)){
  topsnp<-rownames(topPvals)[ii]
  print(topsnp)
  df1<-as.data.frame(t(genos[topsnp,]))
  df2<-du2[,"cdr3_aa_gravy", drop=FALSE]
  rownames(df2)<-substr(rownames(df2), 1, nchar(rownames(m)) - 5)

  df<-merge(df1,df2,by="row.names")
  df[,2]<-as.factor(df[,2])
  row.names(df)<-df$Row.names
  df$Row.names<-NULL
  colnames(df)<-c("a","b")
  print(df)
  print(ggplot(df, aes(x = a, y = b)) + geom_boxplot() +  labs(title = topsnp, x = "genotype", y = "aa_gravy"))
}
dev.off()

##
lapply(genos,)

throwOuts<-NULL
for(i in rownames(g)){
    if(sum(is.na(genos[i,]))>20){throwOuts[[i]]<-i}
#print(sum(is.na(genos[i,])))
}

thdf<-as.data.frame(throwOuts)[1,]

g<-g[-which(rownames(g)%in%thdf[1,]),]

interims = NULL
pvals = NULL

for (i in rownames(g)){
    #i<-"chr7:168152001_C/T"
  snp = as.data.frame(t(g[i,]))
  print(i)
  print(sum(is.na(genos[i,])))
  interims = data.frame(row.names=i)
    for (a in genes){
      str = as.data.frame(du2[,a]) ### Use if not including ethnicity
      rownames(str) = IDs
      p = as.data.frame((summary(lm(str[,1] ~ unlist(snp[,i]), na.action = na.exclude))$coefficients[2,4]))
      rownames(p) = i
      colnames(p) = a
      interims = cbind(interims, p)
    }
  pvals = rbind(pvals, interims)
} 

