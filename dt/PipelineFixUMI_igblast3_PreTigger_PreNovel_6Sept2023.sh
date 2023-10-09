#!/usr/bin/env bash
# Example pRESTO pipeline for UMI barcoded Illumina Miseq 2x250 data
# Data from Stern, Yaari and Vander Heiden et al, 2014, Sci Trans Med.
#
# Author:  Jason Anthony Vander Heiden, Gur Yaari
# Date:    2016.03.04

S=${1}
#S=S18_HC
# Define run parameters and input files
workingDir=/home/datier01/bCellDevelopment
NPROC=12

#S=P1
#R1_FILE=${workingDir}/testR1.fastq
#R2_FILE=${workingDir}/testR2.fastq
#R1_FILE=${workingDir}/data/CW-${S}*R1_001.fastq
#R2_FILE=${workingDir}/data/CW-${S}*R2_001.fastq
OUTDIR=${workingDir}/fastqs
OUTNAME=${S}
#R1_PRIMERS=${workingDir}/CPRIMERS.fasta
#R2_PRIMERS=${workingDir}/VPRIMERS.fasta
#R1_PRIMERS=${workingDir}/IGGPRIMER.fasta
#R2_PRIMERS=${workingDir}/VPRIMERS.fasta
PIPELINE_LOG=Pipeline.log
ZIP_FILES=false

#module load python-3.5.2
#module load java

#module load gcc/7.3.0
#module load anaconda2-5.2.0
#module load java

#source activate ~/anaconda3/envs/prestoChangoEnv/

# Make output directory and empty log files
mkdir -p $OUTDIR
cd $OUTDIR
echo '' > $PIPELINE_LOG

# Make tab directory
#TABDIR=${S}_tab
#mkdir -p $TABDIR 

# Start
echo "OUTPUT DIRECTORY: ${OUTDIR}"
echo -e "START"
STEP=0

function presto {
#Fastqc
printf "  %2d: %-*s $(date +'%H:%M %D')\n" $((++STEP)) 24 "Fastqc"
fastqc -o . $R1_FILE $R2_FILE

## Remove low quality reads
#printf "  %2d: %-*s $(date +'%H:%M %D')\n" $((++STEP)) 24 "FilterSeq quality"
#FilterSeq.py quality -s $R1_FILE -q 20 \
#	--outname "${OUTNAME}-R1" --log FS1.log --nproc $NPROC --outdir . >> $PIPELINE_LOG
#FilterSeq.py quality -s $R2_FILE -q 20 \
#	--outname "${OUTNAME}-R2" --log FS2.log --nproc $NPROC --outdir . >> $PIPELINE_LOG

#Filter by trimming ends of reads
printf "  %2d: %-*s $(date +'%H:%M %D')\n" $((++STEP)) 24 "FilterSeq trimqual"
FilterSeq.py trimqual -s $R1_FILE -q 20 --outname "${OUTNAME}-R1" --outdir . \
	--log FS1.log --nproc $NPROC >> $PIPELINE_LOG
FilterSeq.py trimqual -s $R2_FILE -q 20 --outname "${OUTNAME}-R2" --outdir . \
	--log FS2.log --nproc $NPROC >> $PIPELINE_LOG

#Filter short reads
printf "  %2d: %-*s $(date +'%H:%M %D')\n" $((++STEP)) 24 "FilterSeq length"
FilterSeq.py length -s ${OUTNAME}-R1_trimqual-pass.fastq -n 125 --outname "${OUTNAME}-R1_trimqual-pass" --outdir . \
	--log FS3.log --nproc $NPROC >> $PIPELINE_LOG
FilterSeq.py length -s ${OUTNAME}-R2_trimqual-pass.fastq -n 125 --outname "${OUTNAME}-R2_trimqual-pass" --outdir . \
	--log FS4.log --nproc $NPROC >> $PIPELINE_LOG

# Identify primers and UIDs
printf "  %2d: %-*s $(date +'%H:%M %D')\n" $((++STEP)) 24 "MaskPrimers align"
MaskPrimers.py align -s "${OUTNAME}-R1_trimqual-pass_length-pass.fastq" -p $R1_PRIMERS \
	--mode cut --maxlen 50 --pf PRIMER --maxerror 0.2 --failed --outname "${OUTNAME}-R1" \
	--log MP1.log --nproc $NPROC >> $PIPELINE_LOG

printf "  %2d: %-*s $(date +'%H:%M %D')\n" $((++STEP)) 24 "MaskPrimers extract"
MaskPrimers.py extract -s "${OUTNAME}-R2_trimqual-pass_length-pass.fastq" \
	--mode cut --start 0 --len 12 --pf BARCODE --failed \
	--outname "${OUTNAME}-R2" --log MP2.log --nproc $NPROC >> $PIPELINE_LOG

#MaskPrimers.py score -s "${OUTNAME}-R2_trimqual-pass_length-pass.fastq" -p $R2_PRIMERS \
#	--mode cut --start 12 --barcode --maxerror 0.0 \
#	--outname "${OUTNAME}-R2" --log MP2.log --nproc $NPROC >> $PIPELINE_LOG

# Assign UID to read 2 sequences
printf "  %2d: %-*s $(date +'%H:%M %D')\n" $((++STEP)) 24 "PairSeq"
PairSeq.py -1 "${OUTNAME}-R1_primers-pass.fastq" -2 "${OUTNAME}-R2_primers-pass.fastq" \
    --2f BARCODE --coord sra >> $PIPELINE_LOG


#Allign reads with same UMI
printf "  %2d: %-*s $(date +'%H:%M %D')\n" $((++STEP)) 24 "Align reads with same UMI"
#AlignSets.py muscle -s "${OUTNAME}-R1_primers-pass_pair-pass.fastq" --bf BARCODE --exec ~/anaconda3/envs/prestoChangoEnv/bin/muscle
#AlignSets.py muscle -s "${OUTNAME}-R1_primers-pass_pair-pass.fastq" --bf BARCODE --exec /home/datier01/share/muscle5.1.linux_intel64
AlignSets.py muscle -s "${OUTNAME}-R1_primers-pass_pair-pass.fastq" --bf BARCODE --exec /home/datier01/share/muscle3.8.31_i86linux64

#AlignSets.py muscle -s "${OUTNAME}-R2_primers-pass_pair-pass.fastq" --bf BARCODE --exec ~/anaconda3/envs/prestoChangoEnv/bin/muscle
#AlignSets.py muscle -s "${OUTNAME}-R2_primers-pass_pair-pass.fastq" --bf BARCODE --exec /home/datier01/share/muscle5.1.linux_intel64
AlignSets.py muscle -s "${OUTNAME}-R2_primers-pass_pair-pass.fastq" --bf BARCODE --exec /home/datier01/share/muscle3.8.31_i86linux64

#Cluster reads with same UMI
printf "  %2d: %-*s $(date +'%H:%M %D')\n" $((++STEP)) 24 "Cluster reads with same UMI"
ClusterSets.py set -s "${OUTNAME}-R1_primers-pass_pair-pass_align-pass.fastq" -f BARCODE -k CLUSTER --exec ~/anaconda3/envs/prestoChangoEnv/bin/usearch
ParseHeaders.py copy -s "${OUTNAME}-R1_primers-pass_pair-pass_align-pass_cluster-pass.fastq" -f BARCODE -k CLUSTER --act cat
ClusterSets.py set -s "${OUTNAME}-R2_primers-pass_pair-pass_align-pass.fastq" -f BARCODE -k CLUSTER --exec ~/anaconda3/envs/prestoChangoEnv/bin/usearch
ParseHeaders.py copy -s "${OUTNAME}-R2_primers-pass_pair-pass_align-pass_cluster-pass.fastq" -f BARCODE -k CLUSTER --act cat

# Build UID consensus sequences
printf "  %2d: %-*s $(date +'%H:%M %D')\n" $((++STEP)) 24 "BuildConsensus"
BuildConsensus.py -s "${OUTNAME}-R1_primers-pass_pair-pass_align-pass_cluster-pass_reheader.fastq" --bf CLUSTER --pf PRIMER \
	--prcons 0.6 --maxerror 0.1 --maxgap 0.5 --failed\
	--outname "${OUTNAME}-R1" --log BC1.log --nproc $NPROC >> $PIPELINE_LOG
BuildConsensus.py -s "${OUTNAME}-R2_primers-pass_pair-pass_align-pass_cluster-pass_reheader.fastq" --bf CLUSTER \
	--maxerror 0.1 --maxgap 0.5 --failed\
	--outname "${OUTNAME}-R2" --log BC2.log --nproc $NPROC >> $PIPELINE_LOG

# Synchronize consensus sequence files
printf "  %2d: %-*s $(date +'%H:%M %D')\n" $((++STEP)) 24 "PairSeq"
PairSeq.py -1 "${OUTNAME}-R1_consensus-pass.fastq" -2 "${OUTNAME}-R2_consensus-pass.fastq" \
    --coord presto >> $PIPELINE_LOG

## Assemble paired ends via mate-pair alignment
#printf "  %2d: %-*s $(date +'%H:%M %D')\n" $((++STEP)) 24 "AssemblePairs sequential"
#AssemblePairs.py sequential -1 "${OUTNAME}-R2_consensus-pass_pair-pass.fastq" \
#	-2 "${OUTNAME}-R1_consensus-pass_pair-pass.fastq" \
#	-r ${workingDir}/IMGT_refseq_18.12.11_v-gene.fasta \
#	 --coord presto --rc tail --scanrev --1f CONSCOUNT PRIMER BARCODE --2f PRCONS CONSCOUNT \
#	--aligner blastn --maxerror 0.10 --maxlen 200 --minlen 15 --minident 1\
#	--outname "${OUTNAME}-C" --log AP.log --nproc $NPROC >> $PIPELINE_LOG

# Assemble paired ends via mate-pair alignment
printf "  %2d: %-*s $(date +'%H:%M %D')\n" $((++STEP)) 24 "AssemblePairs align"
AssemblePairs.py align -1 "${OUTNAME}-R2_consensus-pass_pair-pass.fastq" \
	-2 "${OUTNAME}-R1_consensus-pass_pair-pass.fastq" \
	--coord presto --rc tail --scanrev --1f CONSCOUNT PRIMER BARCODE --2f PRCONS CONSCOUNT \
	--maxerror 0.10 --maxlen 200 --minlen 5 --failed\
	--outname "${OUTNAME}-C" --log AP.log --nproc $NPROC >> $PIPELINE_LOG

## Assemble paired ends via mate-pair alignment
#printf "  %2d: %-*s $(date +'%H:%M %D')\n" $((++STEP)) 24 "AssemblePairs sequential"
#AssemblePairs.py sequential -1 "${OUTNAME}-R2_consensus-pass_pair-pass.fastq" \
#	-2 "${OUTNAME}-R1_consensus-pass_pair-pass.fastq" \
#	-r /home/datier01/share/igblast_PreTigger_PreNovel/fasta/imgt_human_ig_v.fasta \
#	 --coord presto --rc tail --scanrev --1f CONSCOUNT PRIMER BARCODE --2f PRCONS CONSCOUNT \
#	--aligner blastn --maxerror 0.10 --maxlen 200 --minlen 5 --failed\
#	--outname "${OUTNAME}-C" --log AP.log --nproc $NPROC >> $PIPELINE_LOG


#Filter short reads
printf "  %2d: %-*s $(date +'%H:%M %D')\n" $((++STEP)) 24 "FilterSeq length"
FilterSeq.py length -s ${OUTNAME}-C_assemble-pass.fastq -n 400 --outname "${OUTNAME}-C_assemble-pass" --outdir . \
	--log APFL1.log --nproc $NPROC >> $PIPELINE_LOG

## Annotate with internal C-region
#printf "  %2d: %-*s $(date +'%H:%M %D')\n" $((++STEP)) 24 "MaskPrimers align"
#MaskPrimers.py align -s "${OUTNAME}-C_assemble-pass.fastq" -p $CREGION_FILE \
#    --maxlen 100 --maxerror 0.3 --mode tag --revpr --skiprc --pf CREGION \
#    --outname "${OUTNAME}-C" --log MP3.log --nproc $NPROC >> $PIPELINE_LOG

# Rewrite header with minimum of CONSCOUNT
printf "  %2d: %-*s $(date +'%H:%M %D')\n" $((++STEP)) 24 "ParseHeaders collapse"
#ParseHeaders.py collapse -s "${OUTNAME}-C_primers-pass.fastq" \
ParseHeaders.py collapse -s "${OUTNAME}-C_assemble-pass_length-pass.fastq" \
	-f CONSCOUNT --act min > /dev/null
    
# Remove duplicate sequences
printf "  %2d: %-*s $(date +'%H:%M %D')\n" $((++STEP)) 24 "CollapseSeq"
CollapseSeq.py -s "${OUTNAME}-C_assemble-pass_length-pass_reheader.fastq" -n 20 \
	--uf CREGION --cf CONSCOUNT --act sum --inner \
	--outname "${OUTNAME}-C" >> $PIPELINE_LOG

# Filter to sequences with at least 2 supporting reads
printf "  %2d: %-*s $(date +'%H:%M %D')\n" $((++STEP)) 24 "SplitSeq group"
SplitSeq.py group -s "${OUTNAME}-C_collapse-unique.fastq" -f CONSCOUNT --num 2 \
    --outname "${OUTNAME}-C" >> $PIPELINE_LOG

# Create tables of final repertoire files
printf "  %2d: %-*s $(date +'%H:%M %D')\n" $((++STEP)) 24 "ParseHeaders table"
ParseHeaders.py table -s "${OUTNAME}-C_atleast-2.fastq" -f ID PRCONS CONSCOUNT DUPCOUNT \
    >> $PIPELINE_LOG

# Process log files
printf "  %2d: %-*s $(date +'%H:%M %D')\n" $((++STEP)) 24 "ParseLog"
ParseLog.py -l FS[1-2].log -f ID QUALITY > /dev/null &
ParseLog.py -l MP[1-3].log -f ID BARCODE PRIMER ERROR > /dev/null &
ParseLog.py -l BC[1-2].log -f BARCODE SEQCOUNT CONSCOUNT PRCONS PRFREQ ERROR \
	> /dev/null &
ParseLog.py -l AP.log -f ID REFID LENGTH OVERLAP GAP ERROR PVALUE EVALUE1 EVALUE2 IDENTITY FIELDS1 FIELDS2 \
    > /dev/null &
wait

# Zip intermediate and log files
if $ZIP_FILES; then
    LOG_FILES_ZIP=$(ls FS[1-2].log MP[1-3].log BC[1-2].log AP.log)
    tar -zcf LogFiles.tar.gz $LOG_FILES_ZIP
    rm $LOG_FILES_ZIP

    TEMP_FILES_ZIP=$(ls *.fastq | grep -vP "collapse-unique.fastq|atleast-2.fastq")
    tar -zcf TempFiles.tar.gz $TEMP_FILES_ZIP
    rm $TEMP_FILES_ZIP
fi
}

#############################Chango Pipeline##############################
function changeo {

inDir=${workingDir}/changeos_igblast_PreTigger_PreNovel_6Sept2023
mkdir -p $inDir
#inDir=${OUTDIR}

fastQFile=$OUTDIR/${OUTNAME}.fastq
fastAFile=$inDir/${OUTNAME}_atleast-2.fasta
fmt7File=$inDir/${OUTNAME}_atleast-2.ftm7
tabFile=$inDir/${OUTNAME}_atleast-2_db-pass.tsv
cloneFile=$inDir/${OUTNAME}_atleast-2_db-pass_clone-pass.tsv

sharedir=/home/datier01/share
sed -n '1~4s/^@/>/p;2~4p' $fastQFile > $fastAFile

## Download reference databases
#VERSION="1.21.0"
#wget ftp://ftp.ncbi.nih.gov/blast/executables/igblast_PreTigger_PreNovel/release/${VERSION}/ncbi-igblast-${VERSION}-x64-linux.tar.gz
#tar -zxf ncbi-igblast-${VERSION}-x64-linux.tar.gz

#${sharedir}/scripts/fetch_igblastdb.sh -o ${sharedir}/igblast_PreTigger_PreNovel
#cp -r ${sharedir}/ncbi-igblast-${VERSION}/internal_data ${sharedir}/igblast_PreTigger_PreNovel
#cp -r ${sharedir}/ncbi-igblast-${VERSION}/optional_file ${sharedir}/igblast_PreTigger_PreNovel


#Start here for making a custom database starting with an existing one with sequences added 11 Sept 2023
#${sharedir}/scripts/fetch_imgtdb.sh -o ${sharedir}/germlines_PreTigger_PreNovel/imgt
## Build IgBLAST database from IMGT reference sequences
#imgt2igblast.sh -i ${sharedir}/germlines_PreTigger_PreNovel/imgt -o ${sharedir}/igblast_PreTigger_PreNovel
#ALSO USE IF YOU WANT TO APPEND NOVEL GAPPED SEQUENCES TO FASTAS AND WANT TO BUILD A NEW DATABASE NO NEED FOR THE TWO SEGMENTS BELOW

##DONT RUN THIS UNLESS YOU WANT A VERY SPECIFIC DATABASE BUILT, THE LINE ABOVE DOES THIS FOR ALL SPECIES AND LOCI
#mkdir ${sharedir}/igblast_PreTigger_PreNovel/fasta
## V segment database
#edit_imgt_file.pl ${sharedir}/germlines_PreTigger_PreNovel/imgt/human/vdj/imgt_human_IGHV.fasta > ${sharedir}/igblast_PreTigger_PreNovel/fasta/imgt_human_ig_v.fasta
#makeblastdb -parse_seqids -dbtype nucl -in ${sharedir}/igblast_PreTigger_PreNovel/fasta/imgt_human_ig_v.fasta \
#    -out ${sharedir}/igblast_PreTigger_PreNovel/database/imgt_human_ig_v
## D segment database
#edit_imgt_file.pl ${sharedir}/germlines_PreTigger_PreNovel/imgt/human/vdj/imgt_human_IGHD.fasta > ${sharedir}/igblast_PreTigger_PreNovel/fasta/imgt_human_ig_d.fasta
#makeblastdb -parse_seqids -dbtype nucl -in ${sharedir}/igblast_PreTigger_PreNovel/fasta/imgt_human_ig_d.fasta \
#    -out ${sharedir}/igblast_PreTigger_PreNovel/database/imgt_human_ig_d
## J segment database
#edit_imgt_file.pl ${sharedir}/germlines_PreTigger_PreNovel/imgt/human/vdj/imgt_human_IGHJ.fasta > ${sharedir}/igblast_PreTigger_PreNovel/fasta/imgt_human_ig_j.fasta
#makeblastdb -parse_seqids -dbtype nucl -in ${sharedir}/igblast_PreTigger_PreNovel/fasta/imgt_human_ig_j.fasta \
#    -out ${sharedir}/igblast_PreTigger_PreNovel/database/imgt_human_ig_j

##IF YOU ALREADY APPENDED NOVEL NON GAPPED SEQUENCES TO FASTAS AND WANT TO BUILD A NEW DATABASE, RUN THIS!!!!
#makeblastdb -parse_seqids -dbtype nucl -in ${sharedir}/igblast_PreTigger_PreNovel/fasta/imgt_human_ig_v.fasta \
#    -out ${sharedir}/igblast_PreTigger_PreNovel/database/imgt_human_ig_v
#makeblastdb -parse_seqids -dbtype nucl -in ${sharedir}/igblast_PreTigger_PreNovel/fasta/imgt_human_ig_d.fasta \
#    -out ${sharedir}/igblast_PreTigger_PreNovel/database/imgt_human_ig_d
#makeblastdb -parse_seqids -dbtype nucl -in ${sharedir}/igblast_PreTigger_PreNovel/fasta/imgt_human_ig_j.fasta \
#    -out ${sharedir}/igblast_PreTigger_PreNovel/database/imgt_human_ig_j

##IF YOU NEED TO GAP YOUR SEQS (D and J dont have gaps)
##python ${sharedir}/gap_inferred.py -h
#python ${sharedir}/gap_inferred.py ${sharedir}/igblast_PreTigger_PreNovel/fasta/imgt_human_ig_v.fasta ${sharedir}/germlines/imgt/human/vdj/imgt_human_IGHV.fasta ${sharedir}/germlines_PreTigger_PreNovel/imgt/human/vdj/imgt_human_IGHV.fasta
#cp ${sharedir}/igblast_PreTigger_PreNovel/fasta/imgt_human_ig_d.fasta ${sharedir}/germlines_PreTigger_PreNovel/imgt/human/vdj/imgt_human_IGHD.fasta
#cp ${sharedir}/igblast_PreTigger_PreNovel/fasta/imgt_human_ig_j.fasta ${sharedir}/germlines_PreTigger_PreNovel/imgt/human/vdj/imgt_human_IGHJ.fasta

export IGDATA=${sharedir}/igblast_PreTigger_PreNovel
#igblastn\
/home/datier01/share/ncbi-igblast-1.21.0/bin/igblastn \
	-germline_db_V ${sharedir}/igblast_PreTigger_PreNovel/database/imgt_human_ig_v\
	-germline_db_D ${sharedir}/igblast_PreTigger_PreNovel/database/imgt_human_ig_d \
	-germline_db_J ${sharedir}/igblast_PreTigger_PreNovel/database/imgt_human_ig_j \
	-auxiliary_data ${sharedir}/igblast_PreTigger_PreNovel/optional_file/human_gl.aux \
	-domain_system imgt -ig_seqtype Ig -organism human \
	-outfmt '7 std qseq sseq btop' \
	-query $fastAFile \
	-out $fmt7File

MakeDb.py igblast -i ${fmt7File} -s ${fastAFile} \
	-r ${sharedir}/germlines_PreTigger_PreNovel/imgt/human/vdj/imgt_human_IGHV.fasta \
	${sharedir}/germlines_PreTigger_PreNovel/imgt/human/vdj/imgt_human_IGHD.fasta \
	${sharedir}/germlines_PreTigger_PreNovel/imgt/human/vdj/imgt_human_IGHJ.fasta \
	${sharedir}/germlines_PreTigger_PreNovel/imgt/human/vdj/imgt_human_IGKV.fasta \
	${sharedir}/germlines_PreTigger_PreNovel/imgt/human/vdj/imgt_human_IGKJ.fasta \
	${sharedir}/germlines_PreTigger_PreNovel/imgt/human/vdj/imgt_human_IGLV.fasta \
	${sharedir}/germlines_PreTigger_PreNovel/imgt/human/vdj/imgt_human_IGLJ.fasta \
	--extended --failed

DefineClones.py -d ${tabFile} --act set --model ham \
    --norm len --dist 0.16

CreateGermlines.py -d ${cloneFile} -g dmask --cloned \
    	-r ${sharedir}/germlines_PreTigger_PreNovel/imgt/human/vdj/imgt_human_IGHV.fasta \
	${sharedir}/germlines_PreTigger_PreNovel/imgt/human/vdj/imgt_human_IGHD.fasta \
	${sharedir}/germlines_PreTigger_PreNovel/imgt/human/vdj/imgt_human_IGHJ.fasta \
	${sharedir}/germlines_PreTigger_PreNovel/imgt/human/vdj/imgt_human_IGKV.fasta \
	${sharedir}/germlines_PreTigger_PreNovel/imgt/human/vdj/imgt_human_IGKJ.fasta \
	${sharedir}/germlines_PreTigger_PreNovel/imgt/human/vdj/imgt_human_IGLV.fasta \
	${sharedir}/germlines_PreTigger_PreNovel/imgt/human/vdj/imgt_human_IGLJ.fasta \

}

#presto
changeo

#End the script
printf "DONE\n\n"
cd ../
