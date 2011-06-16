
## readers 
.coverage_as_dataframe <- function(lst, range)
{
    positive <- sapply(lst, "[[", "+")
    negative <- sapply(lst, "[[", "-")
    if (is.matrix(positive))
        positive <- rowSums(positive)
    if (is.matrix(negative))
        negative <- rowSums(negative)
    pos <- seq.int(start(range), end(range),
                   length.out=length(positive))
    snames <- c("positive", "negative")
    group <- factor(rep(snames, each=length(positive)),
                    levels=snames)
    data.frame(data=c(positive, -negative),
               group=group,
               pos=pos)
}

.fine_coverage_reader <- function(x)
{   ## x: a Snapshot instance
    ## create a plot of coverage, separate lines for each file
    rng <- vrange(x)
    wd <- width(rng)
    cvg <- function(aln) {
        if (identical(0L, length(aln))) 
            numeric(wd)
        else
            as.numeric(coverage(aln, shift=-start(rng)+1, width=wd))
            #as.numeric(coverage(aln, width=wd))
    }
    lst <- lapply(as.list(files(x)), function(fl) {
        aln <- readBamGappedAlignments(fl, which=rng)
        seqlevels(aln) <- seqlevels(rng)
        list(`+`=cvg(aln[strand(aln)=="+"]),
             `-`=cvg(aln[strand(aln)=="-"]))
    })
    .coverage_as_dataframe(lst, rng)
}

.coarse_coverage_reader <- function(x)
{
    nbins <- 10000L
    rng <- vrange(x)
    wd <- width(rng)
    lst <- lapply(as.list(files(x)), function(fl) {
        param <- ScanBamParam(which=rng, what=c("pos", "strand"))
        starts <- scanBam(fl, param=param)[[1]]
        bins <- with(starts, {
            lapply(split(pos, strand)[1:2], function(elt) {
                if (length(elt)) cut(elt, breaks=nbins, labels=FALSE)
                else integer()
            })
        })
        lapply(bins, tabulate, nbins)
    })
    .coverage_as_dataframe(lst, rng)
}

.multifine_coverage_reader <- function(x)
{  
    rng <- vrange(x)
    wd <- width(rng)
    cvg <- function(aln) {
        if (!length(aln))
            numeric(wd)
        else
            as.numeric(coverage(aln, shift=-start(rng)+1, width=wd))
    }
    lst <- lapply(as.list(files(x)), function(fl) {
        aln <- readBamGappedAlignments(fl, which=rng)
        seqlevels(aln) <- seqlevels(rng)
        list(`+`=cvg(aln[strand(aln)=="+"]),
             `-`=cvg(aln[strand(aln)=="-"]))
    })

    positive <- sapply(lst, "[[", "+")
    negative <- sapply(lst, "[[", "-")
    
    pos <- rep(seq.int(start(rng), end(rng)), length(2*length(lst)))

    snames <- c("positive", "negative")
    strand <- rep(snames, each=length(positive))
    fnames <- names(lst)
    file <- rep(do.call(cbind, lapply(names(lst), rep, wd)), 2)
    group <-
        factor(paste(strand, file, sep=": "))
#               levels=c(paste(snames[1], fnames, sep=": "),
#                        paste(snames[2], fnames, sep=": "))

               #levels=paste(rep(fnames, each=2), snames))
    data.frame(data=c(positive, -negative),
               pos=pos, group=group)
}


## viewers
.coverage_viewer <- function(sp)
{
    lty <- rep(seq_len(length(levels(sp$group)) / 2), times=2)
    ## sp: data.frame with "data", "group", and "pos" column
    col <- c("#66C2A5", "#FC8D62")
    cv <- xyplot(data ~ pos, data=sp, group=group, type="s",
                 col=col, 
                 ylab="Coverage", xlab="Coordinate",
                 key=list(space="top", column=2, cex=0.8,
                          lines=list(lty=lty, col=col),
                          text=list(levels(sp$group))),
                 panel=function(...) {
                     panel.xyplot(...)
                     panel.grid(h=0, v=-1)
                     panel.abline(a=0, b=0, alpha=0.3)
                 })

    ## annotation here
    SpTrellis(trellis=cv)   
}

.multicoverage_viewer <- function(sp, ...)
{
    lty <- rep(seq_len(length(levels(sp$group)) / 2), times=2)
    # lty is wrong
    f <- unique(unlist(lapply(strsplit(l, ": "), "[[",2)))
    col <- c(rep("#66C2A5", length(f)), rep("#FC8D62", length(f)))
    cv <- xyplot(data ~ pos, data=sp, group=group, type="s",
                 col=col, lty=lty,
                 ylab="Coverage", xlab="Coordinate",
                 key=list(space="top", column=2, cex=0.8,
                          lines=list(lty=lty, col=col),
                          text=list(levels(sp$group))),
                 panel=function(...) {
                     panel.xyplot(...)
                     panel.grid(h=0, v=-1)
                 })
    SpTrellis(trellis=sv)
}
