# Gabriel Hoffman
# Icahn School of Medicine at Mount Sinai
# September 25, 2018

.ls.objects <- function (pos = 1, pattern, order.by, decreasing = FALSE, head = FALSE, 
    n = 5) 
{
    napply <- function(names, fn) sapply(names, function(x) fn(get(x, 
        pos = pos)))
    names <- ls(pos = pos, pattern = pattern)
    obj.class <- napply(names, function(x) as.character(class(x))[1])
    obj.mode <- napply(names, mode)
    obj.type <- ifelse(is.na(obj.class), obj.mode, obj.class)
    obj.prettysize <- napply(names, function(x) {
        capture.output(print(object.size(x), units = "auto"))
    })
    obj.size <- napply(names, object.size)
    obj.dim <- t(napply(names, function(x) as.numeric(dim(x))[1:2]))
    vec <- is.na(obj.dim)[, 1] & (obj.type != "function")
    obj.dim[vec, 1] <- napply(names, length)[vec]
    out <- data.frame(obj.type, obj.size, obj.prettysize, obj.dim)
    names(out) <- c("Type", "Size", "PrettySize", "Rows", "Columns")
    if (!missing(order.by)) 
        out <- out[order(out[[order.by]], decreasing = decreasing), 
            ]
    if (head) 
        out <- head(out, n)
    out
}




#' Read BED file of TFBS scores into GRange object
#'
#' Read BED file of TFBS scores into GRange object
#'
#' @param grQuery GRange object of query interval
#' @param bedFile BED file that is bgzipped and tabix'd. 
#' @param maxP maximum p-value of TFBS returned 
#' @param quality each motif has a quality value in c('A', 'B', 'C', 'D').  But default exclude quality D, the lowest quality motifs 
#' @return GRange of TFBS locations where p-values are scores
#' @export
#' @import ggplot2 foreach grDevices graphics utils stats
#' @importFrom rtracklayer import
#' @examples
#' library(GenomicRanges)
#' library(tfbsDB)
#' 
#' # range of TFBS to be loaded from file
#' grQuery = GRanges("chr1", IRanges(10280, 10410))
#' 
#' # example data from package
#' bedFile = file.path(system.file("inst/extdata/", package = "tfbsDB"), "motifs_fimo_exmaple.bed.gz")
#' 
#' # read from file
#' gr = readTFBSdb( grQuery, bedFile )
#' 
#' # show TFBS locations
#' plotTFBSdb( gr )
#' 
#' # color locations by p-value of motif match
#' plotTFBSdb( gr, colorByP=TRUE)
#' 
readTFBSdb = function( grQuery, bedFile, maxP=1e-4, quality=c('A', 'B', 'C')){

    # read file
    gr = rtracklayer::import( bedFile, which=grQuery, extraCols=c(qvalue="numeric"))

    if( length(gr) > 0){
        gr$tf = sapply(strsplit(gr$name, '_'), function(x) x[1])
        gr$quality = sapply(strsplit(gr$name, '\\.'), function(x) gsub("^(\\S).*$", "\\1", x[4]))
    }

    # remove intervals with a p-value > maxP
    gr = gr[score(gr) < maxP]

    # remove intervals based on quality
    gr[gr$quality %in% quality]  
}

#' Read BED file of TFBS scores from JASPAR into GRange object
#'
#' Read BED file of TFBS scores from JASPAR into GRange object
#'
#' @param grQuery GRange object of query interval
#' @param bedFile BED file that is bgzipped and tabix'd. 
#' @param maxP maximum p-value of TFBS returned. score corresponds to p-value: 200 -> 10e-2; 300 -> 10e-3 
#' @return GRange of TFBS locations where scores are transformed p-values
#' @export
#' @import ggplot2 foreach grDevices graphics utils stats
#' @importFrom rtracklayer import
#' @examples
#' library(GenomicRanges)
#' library(tfbsDB)
#' 
#' # range of TFBS to be loaded from file
#' grQuery = GRanges("chr1", IRanges(10280, 10410))
#' 
#' # example data from package
#' bedFile = file.path(system.file("inst/extdata/", package = "tfbsDB"), "JASPAR2018_hg19_small.bed.gz")
#' 
#' # read from file
#' gr = readJaspar( grQuery, bedFile )
#' 
#' # show TFBS locations
#' plotTFBSdb( gr )
#' 
#' # color locations by p-value of motif match
#' plotTFBSdb( gr, colorByP=TRUE)
#' 
readJaspar = function( grQuery, bedFile, maxP=1e-4 ){

    gr = rtracklayer::import( bedFile, which=grQuery)
    gr$tf = gr$name

    # convert score to p-value
    gr$score = 10^(-gr$score/100)

    gr[score(gr) < maxP]
}
  

#' Merge overlapping GRanges entries
#'
#' Merge adjacent GRanges entries that overlap by at least minfrac
#'
#' @param x GRange object 
#' @param y GRange object 
#' @param minfrac minimum overlap between to windows to merge entries
#' @importFrom GenomicRanges granges findOverlaps width pintersect reduce 
#' @importFrom S4Vectors queryHits subjectHits
#' @return non-redundant set of intervals as GRanges
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
#' @importFrom GenomicRanges findOverlaps score
#' @importFrom S4Vectors queryHits subjectHits
#' @return non-redundant set of intervals as GRanges
merge_same_tf = function( gr, minfrac=0.05, fxn=min ){

    tfid = ''
    foreach( tfid = unique(gr$tf), .combine=c ) %do% {
        
        # original queries coresponding to tf
        a = gr[gr$tf==tfid]

        # merge overlaps and get now merged intervals: b
        b = mergeOverlapping( a, a, minfrac)
        b$tf = tfid

        # for each entry in b, get the best score of intervals in a that overlap with it
        ovlp = findOverlaps(b, a)
      
        # set p-value
        b$score = sapply( unique(queryHits(ovlp)), function(qh){
            idx = subjectHits(ovlp)[queryHits(ovlp) == qh]
            fxn(score(a)[idx])
        } )

        # set q-value, if it exists
        if( ! is.null( a$qvalue ) ){
             b$qvalue = sapply( unique(queryHits(ovlp)), function(qh){
                idx = subjectHits(ovlp)[queryHits(ovlp) == qh]
                fxn(a$qvalue[idx])
            } )
        }
        b
    }
}

#' Plot TFBS in window
#'
#' Produce ggbio plot that ca be combined with others
#'
#' @param gr GRange object of TFBS 
#' @param xlim range of window in genome coordinates
#' @param tf_text_size size of TF text label
#' @param merge_tfbs should overlapping TFBS from same TF be combined?
#' @param merge_min_frac minimum overlap between to windows to merge TFBS's
#' @param segmentColor color of TFBS segments
#' @param textColor color of TFBS text labels
#' @param colorByP color segments by -log10(p) of motif match
#' @param gradientRange ranges of values for colors when colorByP is TRUE   
#' @return plot as ggbio object 
#' @export
#' @importFrom ggbio autoplot scale_x_sequnit
#' @importFrom GenomicRanges GRanges seqnames restrict
#' @importFrom BiocGenerics start end
#' @importFrom IRanges IRanges
#' @examples
#' library(GenomicRanges)
#' library(tfbsDB)
#' 
#' # range of TFBS to be loaded from file
#' grQuery = GRanges("chr1", IRanges(10280, 10410))
#' 
#' # example data from package
#' bedFile = file.path(system.file("inst/extdata/", package = "tfbsDB"), "motifs_fimo_exmaple.bed.gz")
#' 
#' # read from file
#' gr = readTFBSdb( grQuery, bedFile )
#' 
#' # show TFBS locations
#' plotTFBSdb( gr )
#' 
#' # color locations by p-value of motif match
#' plotTFBSdb( gr, colorByP=TRUE)
#'
plotTFBSdb = function( gr, xlim=c(min(start(gr)), max(end(gr))), tf_text_size=6, merge_tfbs=TRUE, merge_min_frac=0.05, segmentColor="lightblue", textColor="black", colorByP=FALSE, gradientRange=c(4,10) ){

    if( length(table(as.character(gr@seqnames))) > 1){
        stop("Only entries on one chromosome are allow")
    }

    # remove enties that are out of plotting range
    gr_window = GRanges(as.character(gr@seqnames[1]), IRanges(xlim[1], xlim[2]))
    ovlp = GenomicRanges::findOverlaps(gr, gr_window)
    gr = gr[queryHits(ovlp)]

    # if there are no TFBS to plot
    if( length(gr) == 0 ){
        gr_empty = GRanges(seqnames(gr_window), IRanges(0, 0))
        fig = autoplot(gr_empty) + theme_bw(20) + scale_x_sequnit("bp") #+ scale_x_continuous(limits=xlim) # + theme(aspect.ratio=aspect.ratio) , expand=c(0,0)
    }else{

        # merge overlapping TFBS
        if( merge_tfbs ){
            gr = merge_same_tf(gr, merge_min_frac)
        }

        # truncate intervals to match xlim
        gr = restrict(gr, start=xlim[1], end=xlim[2])

        # get segment location in order to write text
        fig = autoplot(gr)
        g = ggplot_build(fig@ggplot)
        g2 = with(g$data[[2]], GRanges(seqnames(gr)@values[1], IRanges(xmin,xmax), y=(ymin+ymax)/2))
        idx = 1:length(gr)
        g2$name = gr$tf[idx]

        if( colorByP ){

            gr$score = pmax(gr$score, 10^-max(gradientRange))

            fig = autoplot(gr, aes(fill=-log10(score))) + theme_bw(20) + scale_x_sequnit("bp") + geom_text(data=as.data.frame(g2), aes(x=(start+end)/2, y=y, label=name), size=tf_text_size, color=textColor ) + scale_color_manual(values=segmentColor) + scale_fill_gradientn(colors=c("white", "red"), limits=gradientRange, name="-log10 p") + theme(legend.position="bottom") + scale_x_continuous(limits=xlim) 
        }else{
            # make plot
            fig = autoplot(gr, aes(color="1", fill="1")) + theme_bw(20) + scale_x_sequnit("bp") + geom_text(data=as.data.frame(g2), aes(x=(start+end)/2, y=y, label=name), size=tf_text_size, color=textColor ) + scale_color_manual(values=segmentColor) + scale_fill_manual(values=segmentColor) + theme(legend.position="none") + scale_x_continuous(limits=xlim) 
        }
    }
    fig 
}

 # plotTFBSdb( gr,  colorByP=FALSE)



