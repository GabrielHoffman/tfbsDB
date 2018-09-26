# Gabriel Hoffman
# Icahn School of Medicine at Mount Sinai
# September 25, 2018
#
# Create file of TFBS motifs scores across genome
# 1) Obtain motifs in MEME format from Cis-BP v1.02 and HOCOMOCO v11
# 	Databases from MEME website:
# 	http://meme-suite.org/meme-software/Databases/motifs/motif_databases.12.18.tgz
#
# 	HOCOMOCO v11: HUMAN/HOCOMOCOv11_full_HUMAN_mono_meme_format.meme
# 	CIS-BP v1.02: CIS-BP/Homo_sapiens.meme
#
# 2) Get genome by chromosomes
# 	get UCSC hg19 genome masked and divided by chromosomes
# 	http://hgdownload.cse.ucsc.edu/goldenPath/hg19/bigZips/chromFaMasked.tar.gz
# 
# 3) Score each motif
#	For each motif, for each chromosome, use FIMO in MEME suite to score motifs
#	This can be expensive, so write jobs first and then run in parallel
#	Create BED file for each motif:chromosome pair
#	Concatenate and then sort all results
#	Final file is motifs_fimo_search.bed.gz with columns:
#	#chrom, start, end, name, pValue, strand, score, qValue

DIR=/sc/orga/scratch/hoffmg01/tfbsDB

# Get motifs
# might have to download manually with browser
MOTIFDB=$DIR/motif_databases
mkdir -p $MOTIFDB
# wget -O $MOTIFDB/motif_databases.12.18.tgz http://meme-suite.org/meme-software/Databases/motifs/motif_databases.12.18.tgz 
cd $MOTIFDB
tar zxvf $MOTIFDB/motif_databases.12.18.tgz -C ../

# get hg19 FASTA by chromosome
FILE=$MOTIFDB/../genome/chromFa.tar.gz
mkdir -p $(dirname $FILE)
wget -O $FILE http://hgdownload.cse.ucsc.edu/goldenPath/hg19/bigZips/chromFa.tar.gz
tar zxvf $FILE -C $(dirname $FILE)

# combine chromosomes
echo $(seq 1 22) X Y M | tr ' ' '\n' | parallel -P1 cat $MOTIFDB/../genome/chr{}.fa > $MOTIFDB/../genome/hg19.fa


# Create scripts for each motif and chromosome
OUTDIR=$DIR/fimo_analysis
mkdir -p $OUTDIR/src/
OUT_SCRIPT=$OUTDIR/src/find_motifs.sh
rm -f ${OUT_SCRIPT}

FIMO=/hpc/users/hoffmg01/build/meme-5.0.2/src/fimo

# Genome-wide FIMO run
######################
FASTA=$MOTIFDB/../genome/hg19.fa

# HOCOMOCOv 11
OUT=$OUTDIR/HOCOMOCOv11
mkdir -p $OUT
MOTIFS=$MOTIFDB/HUMAN/HOCOMOCOv11_full_HUMAN_mono_meme_format.meme
for MOTIF in $(grep MOTIF $MOTIFS | cut -f2 -d' ')
do
	mkdir -p $OUT/$MOTIF
	echo "$FIMO -o $OUT/$MOTIF/hg19 --motif $MOTIF $MOTIFS $FASTA" >> ${OUT_SCRIPT}
done

# CIS-BP
OUT=$OUTDIR/cis_bp
mkdir -p $OUT
MOTIFS=$MOTIFDB/CIS-BP/Homo_sapiens.meme
for MOTIF in $(grep MOTIF $MOTIFS | cut -f2 -d' ')
do
	mkdir -p $OUT/$MOTIF
	echo "$FIMO -o $OUT/$MOTIF/hg19 --motif $MOTIF $MOTIFS $FASTA" >> ${OUT_SCRIPT}
done

# run jobs for each motif
# cat ${OUT_SCRIPT} | parallel -P8

# write jobs files
###################
LOG=$OUTDIR/logs
mkdir -p $LOG $OUTDIR/jobs

cat ${OUT_SCRIPT} | while read CMD; 
do
	ID=$(echo $CMD | cut -f3 -d' ')
	KEY=$(echo $(basename $(dirname $ID))_$(basename $ID))

	echo '#!/bin/bash' > $OUTDIR/jobs/${KEY}.lsf
	echo "#BSUB -J fimo_${KEY}
	#BSUB -P acc_psychencode
	#BSUB -q alloc
	#BSUB -n 1
	#BSUB -R span[hosts=1]
	#BSUB -W 24:00 
	#BSUB -o $LOG/fimo_${KEY}_%J.stdout
	#BSUB -eo $LOG/fimo_${KEY}_%J.stderr
	#BSUB -L /bin/bash

	module purge

	cd $DIR

	$CMD
	" >> $OUTDIR/jobs/${KEY}.lsf
done

module load lsf

# submit jobs
ls $OUTDIR/jobs/*.lsf | head | parallel -P1 "bsub < {}; sleep 1"




# convert GFF to startch to faster processing
rm -f $OUTDIR/src/convert_gff_to_starch.sh
for GFF in $(find $OUTDIR -name fimo.gff)
do
	STRCH=$(echo $GFF | sed 's/gff$/starch/g')
	echo "cat $GFF | gff2starch -  > $STRCH" >> $OUTDIR/src/convert_gff_to_starch.sh
done

cat $OUTDIR/src/convert_gff_to_starch.sh | parallel -P8

# concatenate startch files
starchcat $OUTDIR/*/*/*/fimo.starch > $OUTDIR/motifs_fimo_search.starch

# convert to bed
unstarch $OUTDIR/motifs_fimo_search.starch | bgzip > $OUTDIR/motifs_fimo_search.bed.gz 
tabix -fp bed $OUTDIR/motifs_fimo_search.bed.gz 











sort -k1,1 -k4,4n

/Users/gabrielhoffman/Downloads/fimo_analysis/cis_bp/M0177_1.02/chr2/fimo.gff

# clean up intermediate files
# find $OUTDIR -name fimo.txt -delete
# find $OUTDIR -name fimo.bed.gz -delete
# find $OUTDIR -name fimo.gff -delete
find $OUTDIR -name fimo.html -delete
find $OUTDIR -name fimo.xml -delete
find $OUTDIR -name cisml.xml -delete

# create BED of results
rm -f $OUTDIR/src/create_beds.sh
for FILE in $(find $OUTDIR -name fimo.txt)
do
	BED=$(echo $FILE | sed 's/txt$/bed/g')
	echo -e "#chrom\tstart\tend\tname\tpValue\tstrand\tscore\tqValue" > $BED
	echo "grep -v '^#' $FILE | awk -v OFS='\t' '{print \$2, \$3, \$4, \$1, \$7, \$5, \$6, \$8}' | sort -k1,1 -k2,2n -k3,3n -u >> $BED; bgzip -f $BED" >> $OUTDIR/src/create_beds.sh
done

cat $OUTDIR/src/create_beds.sh | parallel -P8


# combine beds
# sort using tempory disk space
find $OUTDIR -name fimo.bed.gz | xargs cat | zcat | sort -k1,1 -k2,2n -k3,3n -u | bgzip > $OUTDIR/motifs_fimo_search.bed.gz
tabix -fp bed $OUTDIR/motifs_fimo_search.bed.gz

# combine beds
# probably faster then with sort
FILE_LST=$(for FILE in $(find $OUTDIR -name fimo.bed.gz)
do
	echo "<(bgzip -dc $FILE)"
done | xargs)

echo "bedops --header -u $FILE_LST | bgzip > $OUTDIR/motifs_fimo_search.bed.gz" > $OUTDIR/src/combine.sh 

source $OUTDIR/src/combine.sh 


bgzip -dc $OUTDIR/motifs_fimo_search.bed.gz | starch --header - > $OUTDIR/motifs_fimo_search.starch



echo $OUTDIR/motifs_fimo_search.bed.gz
bgzip -dc $OUTDIR/motifs_fimo_search.bed.gz > test.bed



library(Rsamtools)
library(GenomicRanges)
library(readr)


readMotifBed = function( bedFile, grQuery){
	res <- scanTabix(bedFile, param=grQuery)

	dff = Map(function(elt) {
	    read.csv(textConnection(elt), sep="\t", header=FALSE, stringsAsFactors=FALSE)
	}, res)

	dff2 = Map(function(elt) {
	    read.csv(textConnection(elt), sep=';', header=FALSE, stringsAsFactors=FALSE)
	}, res)

	dff = dff[[names(dff)]][,-c(9:10)]
	dff2 = dff2[[names(dff2)]][,-1]

	colnames(dff) = c("chrom", "start", "end", "name", "score", "strand", "software", "origin")

	dff2a = apply(dff2, 2, function(x) sapply(strsplit(x, '='), function(y) y[2]))
	dff2a = dff2a[,apply(dff2a, 2, function(x) ! all(is.na(x)))]
	dff2a = data.frame(dff2a, stringsAsFactors=FALSE)

	colnames(dff2a) = sapply(dff2[1,1:ncol(dff2a)], function(x) sapply(strsplit(x, '='), identity))[1,]

	dff2a$motif = sapply(strsplit(dff2a$ID, '-'), function(x) x[1])

	# add better TF annotation here
	# parse names of 

	df = cbind(dff, dff2a)

	gr = with(df, GRanges(chrom, IRanges(start, end), 
		strand = strand,
		score = as.numeric(score), 
		pValue = as.numeric(pvalue),
		qValue = as.numeric(qvalue),
		pValue = as.numeric(pvalue),
		motif = motif))
	gr
}


bedFile = "~/Downloads/fimo_analysis/motifs_fimo_search.bed.gz"
grQuery = GRanges("chr1", IRanges(0, 13948))

gr = readMotifBed( bedFile, grQuery )

# need to save this to GitHub
# add filter
# run on minerva
# write gr to bed file to be used by DeepFIGV















library(rtracklayer)
library(GenomicRanges)
library(bedr)
library(Rsamtools)

bedFile = "~/Downloads/fimo_analysis/motifs_fimo_search.bed.gz"
grQuery = GRanges("chr1", IRanges(0, 13948))

res = tabix(grQuery, bedFile, check.chr = FALSE)





library(rtracklayer)
library(GenomicRanges)

bedFile = "~/Downloads/fimo_analysis/motifs_fimo_search.bed.gz"
grQuery = GRanges("chr1", IRanges(0, 13948))

grOut = import.bed( bedFile, which=grQuery)



grOut = import.bed( 'test.bed', which=grQuery)


file = './HOCOMOCOv11/ALX1_HUMAN.H11MO.0.B/chr1/fimo.gff'

gr= import( file )


library(genomation)
file = '~/Downloads/fimo_analysis/motifs_fimo_search.bed.gz'
gr = readBed( file )



library(rtracklayer)
file = '~/Downloads/fimo_analysis/motifs_fimo_search.bed.gz'
gr = import( file )

awk -f bed2gff.awk <(bgzip -dc $OUTDIR/motifs_fimo_search.bed.gz) > ~/Downloads/fimo_analysis/output.gff



library(rtracklayer)
file = '~/Downloads/fimo_analysis/output.gff'
gr = import( file )






library(motifbreakR)
library(MotifDb)
library(SNPlocs.Hsapiens.dbSNP142.GRCh37)
library(BSgenome.Hsapiens.UCSC.hg19)




library(PWMEnrich)
library(PWMEnrich.Dmelanogaster.background)
# load the pre-compiled lognormal background
data(PWMLogn.dm3.MotifDb.Dmel)
# load the stripe2 sequences from a FASTA file for motif enrichment
sequence = readDNAStringSet(system.file(package="PWMEnrich",
  dir="extdata", file="stripe2.fa"))
sequence

res = motifEnrichment(sequence, PWMLogn.dm3.MotifDb.Dmel)
## Calculating motif enrichment scores ...
report = sequenceReport(res, 1)


db = readMotifs( '~/Downloads/HOCOMOCOv11_core_pwms_HUMAN_mono.txt')



