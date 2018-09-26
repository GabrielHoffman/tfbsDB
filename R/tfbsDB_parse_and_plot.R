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

# @import rtracklayer foreach S4Vectors grDevices graphics utils stats




#' Read BED file of TFBS scores into GRange object
#'
#' Read BED file of TFBS scores into GRange object
#'
#' @param bedFile BED file that is bgzipped and tabix'd
#' @param grQuery GRange object of query interval
#' @param maxP maximum p-value of TFBS returned
#' @export
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

    foreach( tfid = unique(gr$tf), .combine=c ) %do% {
        
        # original queries coresponding to tf
        a = gr[gr$tf==tfid]

        # merge overlaps and get now merged intervals: b
        b = mergeOverlapping( a, a, minfrac)
        b$tf = tfid

        # for each entry in b, get the best score of intervals in a that overlap with it
        ovlp = findOverlaps(b, a)
        score(b) = sapply( unique(queryHits(ovlp)), function(x) fxn(score(a)[x]))
        b
    }
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
#' @param colorByP color segments by -log10(p) of motif match
#' @param gradientRange ranges of values for colors when colorByP is TRUE  
#' @export
plotTFBSdb = function( gr, xlim=c(start(gr), end(gr)), aspect.ratio=0.1, tf_text_size=6, merge_tfbs=TRUE, merge_min_frac=0.05, segmentColor="lightblue", textColor="black", colorByP=FALSE, gradientRange=c(4,10) ){

    if( length(table(as.character(gr@seqnames))) > 1){
        stop("Only entries on one chromosome are allow")
    }

    # remove enties that are out of plotting range
    gr_window = GRanges(as.character(gr@seqnames[1]), IRanges(xlim[1], xlim[2]))
    ovlp = GenomicRanges::findOverlaps(gr, gr_window)
    gr = gr[queryHits(ovlp)]

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

        fig = autoplot(gr, aes(fill=-log10(score))) + theme_bw(20) + scale_x_sequnit("bp") + theme(aspect.ratio=aspect.ratio)  + geom_text(data=as.data.frame(g2), aes(x=(start+end)/2, y=y, label=name), size=tf_text_size, color=textColor ) + scale_x_continuous(lim=xlim, expand=c(0,0)) + scale_color_manual(values=segmentColor) + scale_fill_gradientn(colors=c("white", "red"), limits=gradientRange, name="-log10 p") + theme(legend.position="bottom")
    }else{
        # make plot
        fig = autoplot(gr, aes(color="1", fill="1")) + theme_bw(20) + scale_x_sequnit("bp") + theme(aspect.ratio=aspect.ratio)  + geom_text(data=as.data.frame(g2), aes(x=(start+end)/2, y=y, label=name), size=tf_text_size, color=textColor ) + scale_x_continuous(lim=xlim, expand=c(0,0)) + scale_color_manual(values=segmentColor) + scale_fill_manual(values=segmentColor) + theme(legend.position="none")
    }
    fig 
}
