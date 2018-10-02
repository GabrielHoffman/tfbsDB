
# Create file of TFBS motif scores across genome
##### Gabriel Hoffman
##### Icahn School of Medicine at Mount Sinai
##### September 25, 2018

## Install package
```R
devtools::install_github("GabrielHoffman/tfbsDB")
```

## Generate TFBS annotation files
 1) Obtain motifs in MEME format from Cis-BP v1.02 and HOCOMOCO v11

	Databases from MEME website:
	http://meme-suite.org/meme-software/Databases/motifs/motif_databases.12.18.tgz

	- HOCOMOCO v11: HUMAN/HOCOMOCOv11_full_HUMAN_mono_meme_format.meme
    - CIS-BP v1.02: CIS-BP/Homo_sapiens.meme

 2) Get genome by chromosomes

	get UCSC hg19 genome masked and divided by chromosomes
	http://hgdownload.cse.ucsc.edu/goldenPath/hg19/bigZips/chromFaMasked.tar.gz
 
 3) Score each motif

- For each motif, use FIMO in MEME suite to score motifs.
- This can be expensive, so write jobs first and then run in parallel. 
- Create BED file for each motif analysis. 
- Concatenate and then sort all results
- Final file is motifs_fimo_search_small.bed.gz with columns:
	chrom, start, end, name, score, strand

### shell code to score motifs and create single BED file
This is designed to run on the Sinai cluster, but is easy to change
``` bash
module load lsf bedops/2.4.20 tabix

DIR=/sc/orga/scratch/hoffmg01/tfbsDB

# Get motifs
# might have to download manually with browser
MOTIFDB=$DIR/motif_databases
mkdir -p $MOTIFDB
# wget -O $MOTIFDB/motif_databases.12.18.tgz http://meme-suite.org/meme-software/Databases/motifs/motif_databases.12.18.tgz 
cd $MOTIFDB
tar zxvf $MOTIFDB/motif_databases.12.18.tgz -C ../

# get motifs directly from HOCOMOCOv11, 
# version on the MEME website is outdated
cd $DIR
wget http://hocomoco11.autosome.ru/final_bundle/hocomoco11/full/HUMAN/mono/HOCOMOCOv11_full_HUMAN_mono_meme_format.meme 

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

# HOCOMOCO v11
# store --max-stored-scores 10000000 in memory
# retain matches with p < 1e-5 
# --no-qvalue because multiple testing burden is wrong
#	because only a single motif is analysed
OUT=$OUTDIR/HOCOMOCOv11
mkdir -p $OUT
MOTIFS=$DIR/HOCOMOCOv11_full_HUMAN_mono_meme_format.meme
for MOTIF in $(grep MOTIF $MOTIFS | cut -f2 -d' ')
do
	mkdir -p $OUT/$MOTIF
	echo "$FIMO -o $OUT/$MOTIF/hg19 --motif $MOTIF --max-stored-scores 10000000 --thresh 1e-5 $MOTIFS $FASTA" >> ${OUT_SCRIPT}
done

## CIS-BP
# exclude for now
# OUT=$OUTDIR/cis_bp
# mkdir -p $OUT
# MOTIFS=$MOTIFDB/CIS-BP/Homo_sapiens.meme
# for MOTIF in $(grep MOTIF $MOTIFS | cut -f2 -d' ')
# do
# 	mkdir -p $OUT/$MOTIF
# 	echo "$FIMO -o $OUT/$MOTIF/hg19 --motif $MOTIF $MOTIFS $FASTA" >> ${OUT_SCRIPT}
# done

# run jobs for each motif
# cat ${OUT_SCRIPT} | parallel -P8
# cat ${OUT_SCRIPT} | grep HOCOMOCO | parallel -P60

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

# submit jobs
ls $OUTDIR/jobs/*.lsf | parallel -P1 "bsub < {}; sleep 1"

# Identify jobs that crashed
# rerun this manually
comm -13 <(find $OUTDIR -name fimo.gff | cut -f9 -d'/' | sort) <(ls $OUTDIR/jobs/*.lsf | parallel -P1 basename | sed 's/_hg19.lsf//g' | sort)

# convert GFF to starch for faster processing
rm -f $OUTDIR/src/convert_gff_to_starch.sh
rm -f $OUTDIR/src/convert_gff_to_bed.sh
for GFF in $(find $OUTDIR -name fimo.gff)
do
	STRCH=$(echo $GFF | sed 's/gff$/starch/g')
	BED=$(echo $GFF | sed 's/gff$/bed.gz/g')
	echo "cat $GFF | gff2starch - > $STRCH" >> $OUTDIR/src/convert_gff_to_starch.sh
	echo "unstarch $STRCH > $BED" >> $OUTDIR/src/convert_gff_to_bed.sh
done

# convert to starch
cat $OUTDIR/src/convert_gff_to_starch.sh | parallel -P60

# convert startcfh to bed
cat $OUTDIR/src/convert_gff_to_bed.sh | parallel -P60

# merge all bed files
bedops -u $(find . -name fimo.bed) | bgzip > $OUTDIR/motifs_fimo_search.bed.gz  

# extract only information need in R
# chrom, start, end, name, pvalue, strand, qvalue
zcat $OUTDIR/motifs_fimo_search.bed.gz | tr ';' '\t' | awk -vOFS="\t" '{print $1, $2, $3, $4, $13, $6, $15}' | sed 's/pvalue=//g' | perl -p -i -e "s/-\d+-chr\d+//g" | bgzip > $OUTDIR/motifs_fimo_search_small.bed.gz 
tabix -fp bed $OUTDIR/motifs_fimo_search_small.bed.gz 
```
Note that qvalue is evaluated per motif


