.bind <- function(lst, elt)
{
    do.call(rbind,
    subListExtract(lst, elt, keep.names=FALSE))
}

## qa summary

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
    q0 <- 1 + 32 * is(quality, "SFastqQuality")
    df <- data.frame(Cycle=as.integer(colnames(abc)[col(abc)]),
                     Quality=q,
                     Score=as.numeric(q)-q0,
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
    with(count[count$Coverage!=0,], {
        res <- tapply(as.numeric(Count), Coverage, sum)
        count <- as.vector(res)
        data.frame(Coverage=as.numeric(names(res)), Count= count,
                   CumulativePpn=cumsum(count) / sum(count),
                   Lane=lane, row.names=NULL)
    })
}

.qa_adapterContamination <-
    function(aln, lane, ..., Lpattern="", Rpattern="",
             max.Lmismatch=.1, max.Rmismatch=.2)
{
    if (missing(Lpattern) && missing(Rpattern)) {
        df <- data.frame(lane=lane, contamination="Not run", row.names=NULL)
        return(df)
    }
    trim <- trimLRPatterns(Lpattern, Rpattern, subject=sread(aln),
                           max.Lmismatch=max.Lmismatch, 
                           max.Rmismatch=max.Rmismatch, 
                           ranges=TRUE)
    ac <- sum(width(trim) < (width(aln) - 9)) / length(trim)
    data.frame(lane=lane, contamination=ac, row.names=NULL)
}

## report-generation

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

.laneLbl <- function(lane) sub("s_(.*)_.*", "\\1", lane)

.plotReadQuality <- function(df)
{
    df$lane <- .laneLbl(df$lane)
    xyplot(density~quality|lane, df,
           type="l",
           xlab="Average (calibrated) base quality",
           ylab="Proportion of reads",
           aspect=2)
}

.plotReadOccurrences <- function(df, ...)
{
    df$lane <- .laneLbl(df$lane)
    df <- with(df, {
        nOccur <- tapply(nOccurrences, lane, c)
        cumulative <- tapply(nOccurrences*nReads, lane, function(elt) {
            cs <- cumsum(elt)
            (cs-cs[1] + 1) / (diff(range(cs)) + 1L)
        })
        lane <- tapply(lane, lane, c)
        data.frame(nOccurrences=unlist(nOccur),
                   cumulative=unlist(cumulative),
                   lane=unlist(lane),
                   row.names=NULL)
    })
    xyplot(cumulative~log10(nOccurrences)|factor(lane), df,
           xlab=expression(paste(
               "Number of occurrences of each sequence (",
               log[10], ")", sep="")),
           ylab="Cumulative proportion of reads",
           aspect=2, panel=function(x, y, ..., type) {
               type <-
                   if (1L == length(x)) "p"
                   else "l"
               panel.xyplot(x, y, ..., type=type)
           }, ...)
}

.freqSequences <- function(qa, read, n=20)
{
    cnt <- qa[["readCounts"]]
    df <- qa[["frequentSequences"]]
    df1 <- df[df$type==read,]
    df1[["ppn"]] <- df1[["count"]] / cnt[df1[["lane"]], read]
    head(df1[order(df1$count, decreasing=TRUE),
             c("sequence", "count", "lane")], n)
}

.plotAlignQuality <- function(df)
{
    df$lane <- .laneLbl(df$lane)
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
    df$lane <- .laneLbl(df$lane)
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
    df$lane <- .laneLbl(df$lane)
    levelplot(cut(score, at)~col*row|lane, df,
              main="Read quality (percentile rank)",
              xlab="Tile x-coordinate",
              ylab="Tile y-coordinate",
              cuts=length(at)-2,
              colorkey=list(labels=.colorkeyNames(at, "%.2f")),
              aspect=2)
}

.plotCycleBaseCall <- function(df) {
    col <- rep(c("red", "blue"), 2)
    lty <- rep(1:2, each=2)
    df <- df[df$Base != "N" & df$Base != "-",]
    df$lane <- .laneLbl(df$lane)
    df$Base <- factor(df$Base)
    xyplot(log10(Count)~as.integer(Cycle)|lane,
           group=factor(Base),
           df[with(df, order(lane, Base, Cycle)),],
           type="l", col=col, lty=lty,
           key=list(space="top",
             lines=list(col=col, lty=lty),
             text=list(lab=levels(df$Base)),
             columns=length(levels(df$Base))),
           xlab="Cycle",
           aspect=2)
}

.plotCycleQuality <-
    function(df, ..., strip=FALSE, strip.left=TRUE)
{
    calc_means <- function(x, y, z)
        rowsum(y * z, x) / rowsum(z, x)

    calc_quantile <- function(x, y, z, q=c(.25, .5, .75))
        by(list(y, z), x, function(x) {
            scoreRle <- Rle(x[[1]], x[[2]])
            quantile(scoreRle, q)
        })

    Lane  <- .laneLbl(df$lane)
    pal <- c("#66C2A5", "#FC8D62") # brewer.pal(3, "Set2")[1:2]
    lvlPal <- c("#F5F5F5", "black" )
    rng <- range(df$Count)
    at <- seq(rng[1], rng[2], length.out=512)
    np <- length(unique(Lane))
    nrow <- ceiling(np / 3)
    layout <- c(ceiling(np/nrow), nrow)

    xyplot(Score ~ Cycle | Lane, df,
           panel=function(x, y, z, ..., groups, subscripts) {
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
           }, ..., layout=layout, ylab="Quality Score",
           strip=strip, strip.left=strip.left)
}

.plotMultipleAlignmentCount <-
    function(df, ...)
{
    df$lane <- .laneLbl(df$lane)
    xyplot(log10(Count)~log10(Matches + 1) | lane, df,
           xlab="log10(Number of matches + 1)", aspect=2, ...)
}


.plotDepthOfCoverage <-
    function(df, ...)
{
    if (is.null(df))
        return(NULL)
    xyplot(CumulativePpn~Coverage | Lane, df, type="b", pch=20,
           scales=list(x=list(log=TRUE)),
           ylab="Cumulative Proportion of Nucleotides", aspect=2, ...)
}

