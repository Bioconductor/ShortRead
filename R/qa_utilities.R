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
                   Lane=lane,
                   row.names=NULL)
    })
}

.qa_adapterContamination <-
    function(aln, lane, ..., Lpattern = "", Rpattern = "")
{
    if (missing(Lpattern) && missing(Rpattern)) {
        df <- data.frame(lane=lane, contamination="Not run", row.names=NULL)
        return(df)
    }
    Lmismatch <-
        ifelse(!missing(Lpattern), ceiling(0.1 * length(Lpattern)), 0)
    Rmismatch <- 
        ifelse(!missing(Rpattern), ceiling(0.2 * length(Rpattern)), 0)
    trim <- trimLRPatterns(Lpattern, Rpattern, subject=sread(aln),
                           max.Lmismatch=Lmismatch,
                           max.Rmismatch=Rmismatch, ranges=TRUE)
    ac <- length(which(width(trim) != width(aln))) / length(id(aln))
    data.frame(lane=lane, contamination=ac, row.names=NULL)
}

## report-generation

.ppnCount <- function(m)
{
    ## if adapterContamination not run
    if(is.factor(m[,-1]))
        return(m)
    else {
    ## scale subsequent columns to be proportions of
    ## first column
        m[,-1] <- round(1000 * m[,-1] / ifelse(is.na(m[,1]), 1, m[,1])) / 1000
        return(m)
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
            (cs-cs[1] + 1) / max(1, diff(range(cs)))
        })
        lane <- tapply(lane, lane, c)
        data.frame(nOccurrences=unlist(nOccur),
                   cumulative=unlist(cumulative),
                   lane=unlist(lane))
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

.plotTileLocalCoords <- function(tile, nrow)
{
    row <- 1 + (tile - 1) %% nrow
    col <- 1 + floor((tile -1) / nrow)
    row[col%%2==0] <- nrow + 1 - row[col%%2==0]
    list(row=as.integer(row), col=as.factor(col))
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

.plotTileCounts <-
    function(df, nrow=.tileGeometry(df$tile)[[1]])
{
    df <- df[!is.na(df$count),]
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
    df <- df[df$Base != "N",]
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

.fivenum <-
    function(rle)
{
    n <- length(rle)
    n4 <- floor((n + 3)/2)/2
    d <- c(1, n4, (n + 1)/2, n + 1 - n4, n)
    0.5 * as.numeric(rle[floor(d)] + rle[ceiling(d)])
}

.boxplot.stats <-
    function(score, count, coef=1.5)
{
    x <- Rle(score, count)
    stats <- .fivenum(x)
    iqr <- diff(stats[c(2, 4)])
    out <- x < (stats[2L] - coef * iqr) | x > (stats[4L] + coef * iqr)
    list(stats = stats, out = runValue(x[out]))
}

.plotCycleQuality <- function(df)
{
    cycleStats <- with(df, {
        tapply(seq_len(nrow(df)), list(lane, Cycle), function(i)
        {
            ns <- .boxplot.stats(Score[i], Count[i])
            data.frame(Score=c(ns$stats, unique(ns$out)), Cycle=Cycle[i][1],
                       Lane=lane[i][1])
        })
    })
    newScore <- do.call(rbind,cycleStats)
    if(length(unique(newScore$Lane)) > 1)
        playout = c(2, ceiling(length(unique(newScore$Lane))/2))
    else playout = c(1,1)
    xyplot(Score ~ Cycle | Lane, data = newScore,
        layout = playout,
        xlab = "Cycle", ylab = "Quality Score",
        panel = function(x, y, ...)
        {
            newy <- split(y, x)
            newx <- sort(unique(x))
            for(i in seq_along(newy))
            {
                elt <- newy[[i]]
                llines(x = newx[i], y = c(elt[1], elt[2]), lty = 2)
                llines(x = newx[i], y = c(elt[4], elt[5]), lty = 2)
                llines(x = newx[i], y = c(elt[2], elt[4]), lty = 1, lwd = 2)
                lpoints(x = newx[i], y = elt[3], pch = "-", col = "black",
                cex = 1.5)
                if(length(elt) > 5)
                    lpoints(x = newx[i], y = elt[-(1:5)], pch = ".")
            }
        }
    )
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
    xyplot(CumulativePpn~Coverage | Lane, df, type="b", pch=20,
           scales=list(x=list(log=TRUE)),
           ylab="Cumulative Proportion of Nucleotides", aspect=2, ...)
}

