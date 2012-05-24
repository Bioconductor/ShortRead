.bind <- function(lst, elt)
{
    do.call(rbind,
    subListExtract(lst, elt, keep.names=FALSE))
}

.qa_alphabetFrequency <-
    function(object, ..., collapse=FALSE, baseOnly=FALSE)
{
    ## avoiding integer overflow in Biostrings::alphabetFrequency
    if (!collapse) {
        msg <- "'collapse' must be TRUE for '.qa_alphabetFrequency'"
        .throw(SRError("InternalError", msg))
    }
    alf <- alphabetByCycle(object)
    mode(alf) <- "numeric"
    alf <- rowSums(alf)
    if (baseOnly && is(object, "DNAStringSet")) {
        bases <- names(Biostrings:::xscodes(object, baseOnly=baseOnly))
        idx <- names(alf) %in% bases
        alf <- c(alf[idx], other=sum(alf[!idx]))
    }
    alf
}

## qa summary

.qa_sampleKey <- function(qa)
    ## use numbers to represent samples, updating all elements of 
{
    value <-rownames(qa[["readCounts"]])
    kv <- data.frame(Key=factor(seq_along(value),
                       levels=seq_along(value)),
                     row.names=value)
    lst <- c(list(keyValue=kv), Map(function(elt, nm, kv) {
        switch(nm,
               readCounts=,
               baseCalls=,
               adapterContamination={
                   rownames(elt) <- kv[rownames(elt), "Key"]
                   elt
               },
               readQualityScore=,
               baseQuality=,
               alignQuality=,
               frequentSequences=,
               sequenceDistribution={
                   elt$lane <- kv[elt$lane, "Key"]
                   elt
               },
               depthOfCoverage={
                   elt$Lane <- kv[elt$Lane, "Key"]
                   elt
               },
               perCycle=,
               perTile={
                   Map(function(elt, nm, kv) {
                       elt$lane <- kv[elt$lane, "Key"]
                       elt
                   }, elt, names(elt), MoreArgs=list(kv))
               },{
                   msg <- sprintf("unhandled QA element '%s'", nm)
                   .throw(SRError("InternalError", msg))
               })
    }, .srlist(qa), names(qa), MoreArgs=list(kv)))
    initialize(qa, .srlist=lst)
}

.qa_qdensity <-
    function(quality)
{
    qscore <- alphabetScore(quality) / width(quality)
    if (length(qscore) >= 2) {
        density(qscore)
    } else {
        pseudo <- list(x=NA, y=NA, bw=Inf, n=0)
        class(pseudo) <- "density"
        pseudo
    }
}

.qa_perCycleBaseCall <-
    function(abc, lane)
{
    if (missing(abc) || dim(abc)[[3]] == 0) {
        df <- data.frame(Cycle=integer(0), Base=factor(),
                         Count=integer(0), lane=character(0))
        return(df)
    }
    abc <- apply(abc, c(1, 3), sum)
    df <- data.frame(Cycle=as.integer(colnames(abc)[col(abc)]),
                     Base=factor(rownames(abc)[row(abc)]),
                     Count=as.vector(abc),
                     lane=lane, row.names=NULL)
    df[df$Count != 0,]
}

.qa_perCycleQuality <-
    function(abc, quality, lane)
{
    if (missing(abc) || dim(abc)[[3]] == 0) {
        df <- data.frame(Cycle=integer(0), Quality=numeric(0),
                         Score=numeric(0), Count=integer(0),
                         lane=character(0))
        return(df)
    }
    abc <- apply(abc, 2:3, sum)
    q <- factor(rownames(abc)[row(abc)], levels=rownames(abc))
    q0 <- as(do.call(class(quality), list(rownames(abc))), "matrix")
    df <- data.frame(Cycle=as.integer(colnames(abc)[col(abc)]),
                     Quality=q, Score=as.integer(q0)[q],
                     Count=as.vector(abc),
                     lane=lane, row.names=NULL)
    df[df$Count != 0, ]
}

.qa_depthOfCoverage <-
    function(aln, lane)
{
    idx <- !is.na(position(aln)) & occurrenceFilter(withSread=FALSE)(aln)
    if (0L == sum(idx)) {
        df <- data.frame(Coverage=character(0), Count=integer(0),
                         CumulativePpn=integer(0), Lane=character(0),
                         row.names=NULL)
        return(df)
    }
    aln <- aln[idx]
    cv <- coverage(aln)
    cvg <- Filter(function(x) length(x) > 0, cv)
    if (0L == length(cvg)) {
        df <- data.frame(Coverage=character(0), Count=integer(0),
                         CumulativePpn=integer(0), Lane=character(0),
                         row.names=NULL)
        return(df)
    }

    ## Each chromosome
    count <-
        do.call(rbind, lapply(seq_len(length(cvg)), function(i, cvg) {
        x <- cvg[[i]]
        res <- tapply(runLength(x), runValue(x), sum)
        data.frame(Coverage=as.numeric(names(res)),
                   Count=as.numeric(res),
                   Seqname=names(cvg)[[i]], row.names=NULL)
        }, cvg))

    ## Entire lane, non-zero coverage
    count <- count[count$Coverage != 0,]
    res <- tapply(as.numeric(count$Count), count$Coverage, sum)
    count <- as.vector(res)
    data.frame(Coverage=as.numeric(names(res)), Count= count,
               CumulativePpn=cumsum(count) / sum(count),
               Lane=lane, row.names=NULL)
}

.qa_adapterContamination <-
    function(aln, lane, ..., Lpattern="", Rpattern="",
             max.Lmismatch=.1, max.Rmismatch=.2, min.trim=9L)
{
    if (missing(Lpattern) && missing(Rpattern)) {
        df <- data.frame(contamination="Not run", row.names=lane)
        return(df)
    }
    trim <- trimLRPatterns(Lpattern, Rpattern, subject=sread(aln),
                           max.Lmismatch=max.Lmismatch, 
                           max.Rmismatch=max.Rmismatch, 
                           ranges=TRUE)
    ac <- sum(width(trim) < (width(aln) - min.trim)) / length(trim)
    data.frame(contamination=ac, row.names=lane)
}

## report-generation

.dnaCol <-                 # brewer.pal(6, "Paired")[c(2, 4, 3, 1, 6)]
    c("#1F78B4", "#33A02C", "#B2DF8A", "#A6CEE3", "#E31A1C")

.ppnCount <- function(m)
{
    if(is.null(m) || is.factor(m[,-1])) {
        "Not available."
    } else {
        ## scale subsequent columns to be proportions of first column
        m[,-1] <-
            round(1000 * m[,-1] / ifelse(is.na(m[,1]), 1, m[,1])) / 1000
        m
    }
}

.df2a <- function(df, fmt="%.3g")
{
    a <-
        if (nrow(df) == 1)
            as.data.frame(lapply(df, sprintf, fmt=fmt))
        else
            sapply(df, sprintf, fmt=fmt)
    row.names(a) <- rownames(df)
    a
}

.plotReadCount <-
    function(qa, ...)
{
    df <- qa[["readCounts"]]
    df1 <- data.frame(Count=unlist(df),
                      Sample=factor(rownames(df), levels=rownames(df)),
                      Census=factor(names(df)[col(df)],
                        levels=names(df)))
    col <- .dnaCol[c(1, 4, 2)]
    dotplot(Sample~Count, group=df1$Census, df1,
            type="b", pch=20, col=col,
            key=list(space="top", lines=list(col=rev(col)),
              text=list(rev(names(df))), columns=ncol(df)))
}

.plotNucleotideCount <-
    function(qa, ...)
{
    df <- qa[["baseCalls"]]
    alph <- df / rowSums(df)
    df1 <- data.frame(Frequency=unlist(alph),
                      Sample=factor(rownames(alph), levels=rownames(alph)),
                      Nucleotide=factor(names(alph)[col(alph)],
                        levels=c("A", "C", "G", "T", "N")))
    dotplot(Sample~Frequency, group=df1$Nucleotide, df1,
            type="b", pch=20,  col=.dnaCol,
            key=list(space="top", lines=list(col=.dnaCol),
              text=list(lab=names(df)), columns=ncol(df)))
}

.plotReadQuality <- function(df, ..., strip=FALSE)
{
    xmin <- min(df$quality)
    ymax <- max(df$density)
    xyplot(density~quality|lane, df,
           type="l",
           xlab="Average (calibrated) base quality",
           ylab="Proportion of reads",
           aspect=2,
           panel=function(..., subscripts) {
               lbl <- as.character(unique(df$lane[subscripts]))
               ltext(xmin, ymax, lbl, adj=c(0, 1))
               panel.xyplot(...)
           },
           strip=FALSE)
}

.plotReadOccurrences <- function(df, ..., strip=FALSE)
{
    df <- local({
        nOccur <- tapply(df$nOccurrences, df$lane, c)
        cumulative <- tapply(df$nOccurrences * df$nReads, df$lane, function(elt) {
            cs <- cumsum(elt)
            (cs-cs[1] + 1) / (diff(range(cs)) + 1L)
        })
        lane <- tapply(df$lane, df$lane, c)
        data.frame(nOccurrences=unlist(nOccur),
                   cumulative=unlist(cumulative),
                   lane=unlist(lane),
                   row.names=NULL)
    })
    xmax <- log10(max(df$nOccurrences))
    xyplot(cumulative~log10(nOccurrences)|factor(lane), df,
           xlab=expression(paste(
               "Number of occurrences of each sequence (",
               log[10], ")", sep="")),
           ylab="Cumulative proportion of reads",
           aspect=2, panel=function(x, y, ..., subscripts, type) {
               lbl <- unique(df$lane[subscripts])
               ltext(xmax, .05, lbl, adj=c(1, 0))
               type <-
                   if (1L == length(x)) "p"
                   else "l"
               panel.xyplot(x, y, ..., type=type)
           }, ..., strip=strip)
}

.freqSequences <- function(qa, read, n=20)
{
    cnt <- qa[["readCounts"]]
    df <- qa[["frequentSequences"]]
    df1 <- df[df$type==read,]
    df1[["ppn"]] <- df1[["count"]] / cnt[df1[["lane"]], read]
    df <- head(df1[order(df1$count, decreasing=TRUE),
                   c("sequence", "count", "lane")], n)
    rownames(df) <- NULL
    df
}

.plotAlignQuality <- function(df)
{
    xyplot(count~score|lane, df,
           type="l",
           prepanel=function(x, y, ...) {
               list(ylim=c(0, 1))
           },
           panel=function(x, y, ...) {
               panel.xyplot(x, y/max(y), ...)
           },
           xlab="Alignment quality score",
           ylab="Number of alignments, relative to lane maximum",
           aspect=2)
}

.atQuantile <- function(x, breaks)
{
    at <- unique(quantile(x, breaks))
    if (length(at)==1)
        at <- at * c(.9, 1.1)
    at
}

.colorkeyNames <- function(at, fmt)
{
    paste(names(at), " (", sprintf(fmt, at), ")", sep="")
}

.tileGeometry <- function(tileIndicies)
{
    n <- as.character(max(tileIndicies))
    switch(n,
           "68"=c(8, 4),
           "100"=c(50, 2),
           "120"=c(60, 2),
           "300"=c(100, 3),
           {
               warning(n, " tiles; ",
                       "assuming lane geometry with 50 tiles / row",
                       call.=FALSE)
               c(50, ceiling(as.integer(n) / 50))
           })
}

.plotTileLocalCoords <- function(tile, nrow)
{
    if (nrow == 8) {
        ## HiSeq
        row <- tile %% 20
        col <- floor(tile / 20) %% 4 + 1L
    } else {
        row <- 1 + (tile - 1) %% nrow
        col <- 1 + floor((tile -1) / nrow)
        row[col%%2==0] <- nrow + 1 - row[col%%2==0]
    }
    list(row=as.integer(row), col=as.factor(col))
}

.plotTileCounts <-
    function(df, nrow=.tileGeometry(df$tile)[[1]])
{
    df <- df[df$count != 0,]
    xy <- .plotTileLocalCoords(df$tile, nrow)
    df[,names(xy)] <- xy
    at <- .atQuantile(df$count, seq(0, 1, .1))
    levelplot(cut(count, at)~col*row|lane, df,
              main="Read count (percentile rank)",
              xlab="Tile x-coordinate",
              ylab="Tile y-coordinate",
              cuts=length(at)-2,
              colorkey=list(labels=.colorkeyNames(at, "%d")),
              aspect=2)
}

.plotTileQualityScore <-
    function(df, nrow=.tileGeometry(df$tile)[[1]])
{
    df <- df[!is.na(df$score),]
    xy <- .plotTileLocalCoords(df$tile, nrow)
    df[,names(xy)] <- xy
    at <- .atQuantile(df$score, seq(0, 1, .1))
    levelplot(cut(score, at)~col*row|lane, df,
              main="Read quality (percentile rank)",
              xlab="Tile x-coordinate",
              ylab="Tile y-coordinate",
              cuts=length(at)-2,
              colorkey=list(labels=.colorkeyNames(at, "%.2f")),
              aspect=2)
}

.plotCycleBaseCall <- function(df, ..., strip=FALSE) {
    col <- .dnaCol[1:4]
    df <- df[df$Base != "N" & df$Base != "-",]
    df$Base <- factor(df$Base)
    xmax <- max(df$Cycle)
    ymax <- log10(max(df$Count))
    xyplot(log10(Count)~as.integer(Cycle)|lane,
           group=factor(Base),
           df[order(df$lane, df$Base, df$Cycle),],
           panel=function(..., subscripts) {
               lbl <- as.character(unique(df$lane[subscripts]))
               ltext(xmax, ymax, lbl, adj=c(1, 1))
               panel.xyplot(..., subscripts=subscripts)
           },
           type="l", col=col,
           key=list(space="top", lines=list(col=col, lwd=2),
             text=list(lab=levels(df$Base)),
             columns=length(levels(df$Base))),
           xlab="Cycle", aspect=2, strip=strip, ...)
}

.plotCycleQuality <-
    function(df, ..., strip=FALSE, strip.left=FALSE)
{
    calc_means <- function(x, y, z)
        rowsum(y * z, x) / rowsum(z, x)

    calc_quantile <- function(x, y, z, q=c(.25, .5, .75))
        by(list(y, z), x, function(x) {
            scoreRle <- Rle(x[[1]], x[[2]])
            quantile(scoreRle, q)
        })

    Lane  <- df$lane
    pal <- c("#66C2A5", "#FC8D62") # brewer.pal(3, "Set2")[1:2]
    lvlPal <- c("#F5F5F5", "black" )
    rng <- range(df$Count)
    at <- seq(rng[1], rng[2], length.out=512)
    np <- length(unique(Lane))
    nrow <- ceiling(np / 4)
    layout <- c(ceiling(np/nrow), nrow)
    ymin <- min(df$Score)

    xyplot(Score ~ Cycle | Lane, df,
           panel=function(x, y, ..., subscripts) {
               z <- df$Count[subscripts]
               mean <- calc_means(x, y, z)
               qtiles <- calc_quantile(x, y, z)
               sxi <- sort(unique(x))
               panel.levelplot(x, y, z, subscripts=TRUE, at=at,
                               col.regions=colorRampPalette(lvlPal))
               llines(sxi, mean, type="l", col=pal[[1]], lwd=1)
               llines(sxi, sapply(qtiles, "[[", 1),
                      type="l", col=pal[[2]], lwd=1, lty=3)
               llines(sxi, sapply(qtiles, "[[", 2),
                      type="l", col=pal[[2]], lwd=1)
               llines(sxi, sapply(qtiles, "[[", 3),
                      type="l", col=pal[[2]], lwd=1, lty=3)
               lbl <- as.character(unique(df$lane[subscripts]))
               ltext(1, ymin, lbl, adj=c(0, 0))
           }, ..., ylab="Quality Score", layout=layout,
           strip=strip, strip.left=strip.left)
}

.plotMultipleAlignmentCount <-
    function(df, ...)
{
    xyplot(log10(Count)~log10(Matches + 1) | lane, df,
           xlab="log10(Number of matches + 1)", aspect=2, ...)
}


.plotDepthOfCoverage <-
    function(df, ..., strip=FALSE)
{
    if (is.null(df))
        return(NULL)
    xmin <- log(min(df$Coverage))
    ymax <- max(df$CumulativePpn)
    xyplot(CumulativePpn~Coverage | Lane, df, type="b", pch=20,
           scales=list(x=list(log=TRUE)),
           ylab="Cumulative Proportion of Nucleotides", aspect=2,
           panel=function(..., subscripts) {
               lbl <- as.character(unique(df$Lane[subscripts]))
               ltext(xmin, ymax, lbl, adj=c(0, 1))
               panel.xyplot(..., subscripts=subscripts)
           }, ..., strip=strip)
}

