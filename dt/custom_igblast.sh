#requires a conda environment with presto and changeo installed

sharedir=/home/datier01/share
igblast_dir_name=igblast_PreTigger
germlines_dir_name=germlines_PreTigger_Test
species=human

## Download reference databases
#VERSION="1.21.0"
#wget ftp://ftp.ncbi.nih.gov/blast/executables/${igblast_dir_name}/release/${VERSION}/ncbi-igblast-${VERSION}-x64-linux.tar.gz
#tar -zxf ncbi-igblast-${VERSION}-x64-linux.tar.gz

#Download igblastdb. Copy files from igblast download to igblastdb
#${sharedir}/scripts/fetch_igblastdb.sh -o ${sharedir}/${igblast_dir_name}
#cp -r ${sharedir}/ncbi-igblast-${VERSION}/internal_data ${sharedir}/${igblast_dir_name}
#cp -r ${sharedir}/ncbi-igblast-${VERSION}/optional_file ${sharedir}/${igblast_dir_name}

#Fetch imhgtdb to germlines
#${sharedir}/scripts/fetch_imgtdb.sh -o ${sharedir}/${germlines_dir_name}/imgt
# Build IgBLAST database from IMGT reference sequences
#imgt2igblast.sh -i ${sharedir}/${germlines_dir_name}/imgt -o ${sharedir}/${igblast_dir_name}

#Start here for making a custom database starting with an existing one with novel non gapped sequences added to igblast folder (11 Sept 2023)
# add novel sequences to the end of ${igblast_dir_name}/fasta/imgt_${species}_ig_x.fasta where x is v,d,j
##IF YOU ALREADY APPENDED NOVEL NON GAPPED SEQUENCES TO FASTAS AND WANT TO BUILD A NEW DATABASE, RUN THIS!!!!
makeblastdb -parse_seqids -dbtype nucl -in ${sharedir}/${igblast_dir_name}/fasta/imgt_${species}_ig_v.fasta \
    -out ${sharedir}/${igblast_dir_name}/database/imgt_${species}_ig_v
makeblastdb -parse_seqids -dbtype nucl -in ${sharedir}/${igblast_dir_name}/fasta/imgt_${species}_ig_d.fasta \
    -out ${sharedir}/${igblast_dir_name}/database/imgt_${species}_ig_d
makeblastdb -parse_seqids -dbtype nucl -in ${sharedir}/${igblast_dir_name}/fasta/imgt_${species}_ig_j.fasta \
    -out ${sharedir}/${igblast_dir_name}/database/imgt_${species}_ig_j

##IF YOU NEED GAPS IN YOUR GERMLINE SEQS (D and J dont have gaps)
##python ${sharedir}/gap_inferred.py -h
python ${sharedir}/gap_inferred.py ${sharedir}/${igblast_dir_name}/fasta/imgt_${species}_ig_v.fasta ${sharedir}/germlines/imgt/${species}/vdj/imgt_${species}_IGHV.fasta ${sharedir}/${germlines_dir_name}/imgt/${species}/vdj/imgt_${species}_IGHV.fasta
cp ${sharedir}/${igblast_dir_name}/fasta/imgt_${species}_ig_d.fasta ${sharedir}/${germlines_dir_name}/imgt/${species}/vdj/imgt_${species}_IGHD.fasta
cp ${sharedir}/${igblast_dir_name}/fasta/imgt_${species}_ig_j.fasta ${sharedir}/${germlines_dir_name}/imgt/${species}/vdj/imgt_${species}_IGHJ.fasta
