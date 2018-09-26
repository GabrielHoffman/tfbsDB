
# Create file of TFBS motif scores across genome
##### Gabriel Hoffman
##### Icahn School of Medicine at Mount Sinai
##### September 25, 2018

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
ls $OUTDIR/jobs/*.lsf | head | parallel -P1 "bsub < {}; sleep 1"

# gzip results
find . -name gff | parallel gzip
find . -name txt | parallel gzip
find . -name xml | parallel gzip
find . -name html | parallel gzip

# convert GFF to startch to faster processing
rm -f $OUTDIR/src/convert_gff_to_starch.sh
for GFF in $(find $OUTDIR -name fimo.gff.gz)
do
	STRCH=$(echo $GFF | sed 's/gff$/starch/g')
	echo "zcat $GFF | gff2starch - > $STRCH" >> $OUTDIR/src/convert_gff_to_starch.sh
done

# combine all starh files
# this is much faster then combining the BED files directly
cat $OUTDIR/src/convert_gff_to_starch.sh | parallel -P60

# concatenate startch files
# faster than default using --bzip2
starchcat --gzip $OUTDIR/*/*/*/fimo.starch > $OUTDIR/motifs_fimo_search.starch

# convert to bed
unstarch $OUTDIR/motifs_fimo_search.starch | bgzip > $OUTDIR/motifs_fimo_search.bed.gz 
tabix -fp bed $OUTDIR/motifs_fimo_search.bed.gz 

# extract only information need in R
zcat $OUTDIR/motifs_fimo_search.bed.gz | tr ';' '\t' | awk -vOFS="\t" '{print $1, $2, $3, $4, $13, $6}' | sed 's/pvalue=//g' | bgzip > $OUTDIR/motifs_fimo_search_small.bed.gz 
tabix -fp bed $OUTDIR/motifs_fimo_search_small.bed.gz  
```

### R code to parse GRanges
``` R
library(rtracklayer)
library(GenomicRanges)

readTFBSdb = function( bedFile, grQuery){
	gr = rtracklayer::import( bedFile, which=grQuery)
	gr$tf = sapply(strsplit(gr$name, '_'), function(x) x[1])
	gr$quality = sapply(strsplit(gr$name, '\\.'), function(x) gsub("^(\\S).*$", "\\1", x[4]))
	gr
}

bedFile = "motifs_fimo_search_small.bed.gz"
grQuery = GRanges("chr9", IRanges(0, 1214000))

gr = readTFBSdb( bedFile, grQuery )
```

module unload R
module load R/3.2.3 htslib/1.9
R

suppressPackageStartupMessages(library(ggbio))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(gridExtra))
suppressPackageStartupMessages(library(GenomicRanges))
suppressPackageStartupMessages(library(rtracklayer))
suppressPackageStartupMessages(library(foreach))

#' Read BED file of TFBS scores into GRange object
#'
#' Read BED file of TFBS scores into GRange object
#'
#' @param bedFile BED file that is bgzipped and tabix'd
#' @param grQuery GRange object of query interval
#' @param maxP maximum p-value of TFBS returned
readTFBSdb = function( bedFile, grQuery, maxP=1e-4){
	gr = rtracklayer::import( bedFile, which=grQuery)
	if( length(gr) > 0){
		gr$tf = sapply(strsplit(gr$name, '_'), function(x) x[1])
		gr$quality = sapply(strsplit(gr$name, '\\.'), function(x) gsub("^(\\S).*$", "\\1", x[4]))
	}
	gr[score(gr) < maxP]
}

#' Merge overlapping GRanges entries
#'
#' Merge adjacent GRanges entries that overlap by at least minfrac
#'
#' @param x GRange object 
#' @param y GRange object 
#' @param minfrac minimum overlap between to windows to merge entries
mergeOverlapping = function(x, y, minfrac=0.05) {
    x <- granges(x)
    y <- granges(y)
    hits <- findOverlaps(x, y)
    xhits <- x[queryHits(hits)]
    yhits <- y[subjectHits(hits)]
    frac <- width(pintersect(xhits, yhits)) / pmin(width(xhits), width(yhits))
    merge <- frac >= minfrac
    c(reduce(c(xhits[merge], yhits[merge])),
      xhits[!merge], yhits[!merge],
      x[-queryHits(hits)], y[-subjectHits(hits)])
}


#' Merge overlapping TFBS
#'
#' Merge adjacent TFBS that overlap by at least minfrac
#'
#' @param gr GRange object of TFBS 
#' @param minfrac minimum overlap between to windows to merge TFBS's
#' @param fxn score of merged interval is the assigned as fxn(score(.)) of overlapping intervals
merge_same_tf = function( gr, minfrac=0.05, fxn=min ){

	m = foreach( tfid = unique(gr$tf), .combine=c ) %do%{
		
		# original queries coresponding to tf
		a = gr[gr$tf==tfid]

		# merge overlaps and get now merged intervals: b
		b = mergeOverlapping( a, a, minfrac)
		b$tf = tfid

		# for each entry in b, get the best score of intervals in a that overlap with it
		ovlp = findOverlaps(b, a)
		score(b) = sapply( unique(ovlp@queryHits), function(x) fxn(score(a)[x]))
		b
	}
	m
}

#' Plot TFBS in window
#'
#' Produce ggbio plot that ca be combined with others
#'
#' @param gr GRange object of TFBS 
#' @param xlim range of window in genome coordinates
#' @param aspect.ratio aspect ratio of plot
#' @param tf_text_size size of TF text label
#' @param merge_tfbs should overlapping TFBS from same TF be combined?
#' @param merge_min_frac minimum overlap between to windows to merge TFBS's
#' @param segmentColor color of TFBS segments
#' @param textColor color of TFBS text labels


plotTFBSdb = function( gr, xlim=c(start(gr), end(gr)), aspect.ratio=0.1, tf_text_size=6, merge_tfbs=TRUE, merge_min_frac=0.05, segmentColor="lightblue", textColor="black", colorByP=FALSE, maxGradientValue=10 ){

	if( length(table(as.character(gr@seqnames))) > 1){
		stop("Only entries on one chromosome are allow")
	}

	# remove enties that are out of plotting range
	gr_window = GRanges(as.character(gr@seqnames[1]), IRanges(xlim[1], xlim[2]))
	ovlp = findOverlaps(gr, gr_window)
	gr = gr[ovlp@queryHits]

	# merge overlapping TFBS
	if( merge_tfbs ){
		gr = merge_same_tf(gr, merge_min_frac)
	}

	# get segment location in order to write text
	fig = autoplot(gr)
	g = ggplot_build(fig)
	g2 = with(g$data[[2]], GRanges(seqnames(gr)@values[1], IRanges(xmin,xmax), y=(ymin+ymax)/2))
	idx = 1:length(gr)
	g2$name = gr$tf[idx]

	if( colorByP ){

		gr$score = pmax(gr$score, 10^-maxGradientValue)

		fig = autoplot(gr, aes(fill=-log10(score))) + theme_bw(20) + scale_x_sequnit("bp") + theme(aspect.ratio=aspect.ratio)  + geom_text(data=as.data.frame(g2), aes(x=(start+end)/2, y=y, label=name), size=tf_text_size, color=textColor ) + scale_x_continuous(lim=xlim, expand=c(0,0)) + scale_color_manual(values=segmentColor) + scale_fill_gradientn(colors=c("white", "red"), limits=c(4, maxGradientValue), name="-log10 p") + theme(legend.position="bottom")
	}else{
		# make plot
		fig = autoplot(gr, aes(color="1", fill="1")) + theme_bw(20) + scale_x_sequnit("bp") + theme(aspect.ratio=aspect.ratio)  + geom_text(data=as.data.frame(g2), aes(x=(start+end)/2, y=y, label=name), size=tf_text_size, color=textColor ) + scale_x_continuous(lim=xlim, expand=c(0,0)) + scale_color_manual(values=segmentColor) + scale_fill_manual(values=segmentColor) + theme(legend.position="none")
	}
	fig 
}

# xlim=c(start(gr), end(gr))
# aspect.ratio=0.1
# tf_text_size=6
# merge_tfbs=TRUE
# merge_min_frac=0.1
# segmentColor="lightblue"
# textColor="black"
# olorByP=FALSE 

# bedFile = "~/Downloads/motifs_fimo_search_small.bed.gz"
bedFile = '/sc/orga/scratch/hoffmg01/tfbsDB/fimo_analysis/motifs_fimo_search_small.bed.gz'
grQuery = GRanges("chr9", IRanges(10000, 10150))

gr = readTFBSdb( bedFile, grQuery )


file = "~/www/dream/test.png"

png( file, width=1300)
plotTFBSdb( gr, xlim=c(10000, 10150), colorByP=FALSE)
plotTFBSdb( gr, xlim=c(10000, 10150), colorByP=TRUE)
dev.off()










