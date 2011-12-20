## QAData

QAData <-
    function(seq=ShortReadQ(), filter=logical(length(seq)), ...)
{
    .QAData$new(seq=seq, filter=filter, ...)
}

setMethod(.filter, "QAData", function(object, useFilter, ...) {
    if (useFilter)
        object$seq[!object$filter]
    else
        object$seq
})

.filterUpdate <-
    function(object, add, value)
{
    if (add)
        object$filter <- object$filter | value
    object
}

## QASource

QAFastqSource <-
    function(con=character(), n=1e6, readerBlockSize=1e8, ...)
{
    new("QAFastqSource", con=as.character(con),
        n=mkScalar(as.integer(n)),
        readerBlockSize=mkScalar(as.integer(readerBlockSize)),
        ...)
}

setMethod(show, "QAFastqSource", function(object) {
    callNextMethod()
    cat("n: ", object@n, ";",
        " readerBlockSize: ", object@readerBlockSize, "\n", sep="")
})

## QASummary

.QASummary <- 
    function (class, useFilter = TRUE, addFilter = TRUE, ..., html) 
{
    if (missing(html)) 
        html <- file.path(system.file("template", package = "ShortRead"), 
                          sprintf("%s.html", class))
    if (!is.na(html) && (!file.exists(html) || !nzchar(html))) 
        .throw(SRError("UserArgumentMismatch",
                       "'html' file does not exist:\n    %s",
                       html))
    new(class, useFilter = mkScalar(as.logical(useFilter)),
        addFilter = mkScalar(as.logical(addFilter)), 
        html = mkScalar(html), ...)
}

.QASummaryFactory <-
    function(summaryName)
{
    function (useFilter = TRUE, addFilter = TRUE, ...) 
        .QASummary(summaryName, useFilter = useFilter,
                   addFilter = addFilter, ...)
}

setMethod(show, "QASummary", function(object) {
    callNextMethod()
    cat("useFilter: ", object@useFilter, "; ",
        "addFilter: ", object@addFilter, "\n", sep="")
})

QASources <- .QASummaryFactory("QASources")

QAFiltered <- .QASummaryFactory("QAFiltered")

QANucleotideUse <- .QASummaryFactory("QANucleotideUse")

QAQualityUse <- .QASummaryFactory("QAQualityUse")

QASequenceUse <- .QASummaryFactory("QASequenceUse")

QAReadQuality <- .QASummaryFactory("QAReadQuality")

QAAdapterContamination <-
    function (useFilter=TRUE, addFilter=TRUE,
              Lpattern = NA_character_, Rpattern = NA_character_, 
              max.Lmismatch = 0.1, max.Rmismatch = 0.2, min.trim = 9L,
              ...)
{
    fmt <- "QAAdapterContamination not a DNA sequence\n  %s=\"%s\""
    if (!is.na(Lpattern)) 
        tryCatch(DNAString(Lpattern), error = function(e) {
            .throw(SRError("UserArgumentMismatch", fmt, "Lpattern", 
                Lpattern))
        })
    if (!is.na(Rpattern)) 
        tryCatch(DNAString(Rpattern), error = function(e) {
            .throw(SRError("UserArgumentMismatch", fmt, "Rpattern", 
                Rpattern))
        })
    .QASummary("QAAdapterContamination",
               useFilter=useFilter, addFilter=addFilter,
               Lpattern = mkScalar(toupper(as.character(Lpattern))), 
               Rpattern = mkScalar(toupper(as.character(Rpattern))), 
               max.Lmismatch = mkScalar(as.numeric(max.Lmismatch)), 
               max.Rmismatch = mkScalar(as.numeric(max.Rmismatch)), 
               min.trim = mkScalar(as.integer(min.trim)), ...)
}

setMethod(show, "QAAdapterContamination", function(object) {
    callNextMethod()
    cat("Lpattern:", object@Lpattern, "\n")
    cat("Rpattern:", object@Rpattern, "\n")
    cat("max.Lmismatch: ", object@max.Lmismatch, "; ",
        "max.Rmismatch: ", object@max.Rmismatch, "; ",
        "min.trim: ", object@min.trim, "\n", sep="")
})

QAFrequentSequence <-
    function (useFilter = TRUE, addFilter = TRUE,
              n = NA_integer_, k = NA_integer_,
              reportSequences = FALSE, ...) 
{
    .QASummary("QAFrequentSequence",
               addFilter = addFilter, useFilter = useFilter, 
               n = mkScalar(as.integer(n)), k = mkScalar(as.integer(k)), 
               reportSequences = mkScalar(as.logical(reportSequences)), 
               ...)
}

setMethod(show, "QAFrequentSequence", function(object) {
    callNextMethod()
    if (!is.na(object@n))
        cat("n: ", object@n, "; ", sep="")
    else
        cat("k: ", object@k, "; ", sep="")
    cat("reportSequences:", object@reportSequences, "\n")
})

QANucleotideByCycle <- .QASummaryFactory("QANucleotideByCycle")

QAQualityByCycle <- .QASummaryFactory("QAQualityByCycle")

## QACollate

setMethod(QACollate, "missing",
          function(src, ...)
{
    QACollate(QAFastqSource(), ...)
})

setMethod(QACollate, "QAFastqSource",
          function(src, ...)
{
    if (1L == length(list(...)) && is(..1, "QACollate"))
        renew(..1, src=src)
    else
        new("QACollate", src=src, SimpleList(...))
})

setMethod(show, "QACollate", function(object) {
    callNextMethod()
    cat("source:", class(object@src),
        "of length", length(object@src@con), "\n")
    elts <- paste(sapply(object, class), collapse = " ")
    txt <- paste(strwrap(sprintf("elements: %s", elts), exdent = 2), 
                 collapse = "\n  ")
    cat(txt, "\n")
})

## QA

QA <- 
    function (sources, filtered, ...) 
{
    new("QA", sources = sources, filtered = filtered, ...)
}

## .clone

setMethod(.clone, "QAData",
          function (object, ...) 
{
    .QAData$new(seq = object$seq, filter = object$filter, ...)
})

setMethod(.clone, "QASource",
          function (object, ...) 
{
    object@data <- .clone(object@data)
    object
})

## values

setMethod(values, "QASummary",
          function(x, ...)
{
    x@values
})

setReplaceMethod("values", c("QASummary", "DataFrame"),
                 function (x, ..., value) 
{
    x@values <- value
    x
})

## rbind

setMethod(rbind, "QASummary",
          function(..., deparse.level=1)
{
    class <- class(..1)
    values <- do.call(rbind, lapply(list(...), values))
    renew(..1, values = values)
})

## qa2

setMethod(qa2, "FastqSampler",
          function(object, state, ..., verbose=FALSE)
{
    if (verbose) message("qa2,FastqSampler-method")
    state$seq <- yield(object)
    state$filter <- rep(FALSE, length(state$seq))
    DataFrame(SourceN=object$status()[["total"]],
              SampleN=length(state$seq))
})

setMethod(qa2, "QAFastqSource",
          function(object, state, ..., verbose=FALSE) 
{
    if (verbose) message("qa2,QAFastqSource-method")
    if (1 != length(object@con))
        .throw(SRError("InternalError",
                       "'QAFastqSource' source length != 1"))
    df <- qa2(FastqSampler(object@con, object@n,
                           object@readerBlockSize),
              object@data, verbose=verbose)
    values <-
        cbind(df, DataFrame(AccessTimestamp=date(),
                            FileName=basename(object@con)))
                            ## Path=dirname(path(object@con))))
    metadata(values) <- list(NumberOfRecords=length(object@data$seq))
    renew(object, values=values)
})

setMethod(qa2, "QAAdapterContamination",
          function(object, state, ..., verbose=FALSE)
{
    if (verbose) message("qa2,QAAdapterContamination-method")
    obj <- .filter(state@data, object@useFilter)

    Lpattern <-
        if (is.na(object@Lpattern)) "" else object@Lpattern
    Rpattern <-
        if (is.na(object@Rpattern)) "" else object@Rpattern

    trim <- trimLRPatterns(Lpattern, Rpattern, sread(obj),
                           object@max.Lmismatch,
                           object@max.Rmismatch,
                           ranges=TRUE)
    filt <- width(trim) < (width(obj) - object@min.trim)

    .filterUpdate(state@data, object@addFilter, filt)
    values <- DataFrame(Contaminants=sum(filt))
    metadata(values) <- list(NumberOfRecords=length(filt))
    renew(object, values=values)
})

setMethod(qa2, "QANucleotideUse",
          function(object, state, ..., verbose=FALSE)
{
    if (verbose) message("qa2,QANucleotideUse-method")
    obj <- .filter(state@data, object@useFilter)
    alf <- .qa_alphabetFrequency(sread(obj), baseOnly=TRUE,
                                 collapse=TRUE)
    values <- DataFrame(Nucleotide=factor(sub("other", "N", names(alf)),
                          levels=c("A", "C", "G", "T", "N")),
                         Count=as.vector(alf))
    metadata(values) <- list(NumberOfRecords=length(obj))
    renew(object, values=values)
})

setMethod(qa2, "QAQualityUse",
          function(object, state, ..., verbose=FALSE)
{
    if (verbose) message("qa2,QAQualityUse-method")
    obj <- .filter(state@data, object@useFilter)
    alf <- .qa_alphabetFrequency(quality(obj), collapse=TRUE)
    alf <- alf[alf != 0]
    quality <- factor(names(alf), levels=alphabet(quality(obj)))
    values <- DataFrame(Quality=quality, Count=as.vector(alf))
    metadata(values) <- list(NumberOfRecords=length(obj))
    renew(object, values=values)
})

setMethod(qa2, "QASequenceUse",
          function(object, state, ..., verbose=FALSE)
{
    if (verbose) message("qa2,QASequenceUse-method")
    obj <- .filter(state@data, object@useFilter)
    t <- tabulate(tabulate(srrank(sread(obj))))
    values <- DataFrame(Occurrences=seq_along(t)[t!=0],
                        Reads=t[t!=0])
    metadata(values) <- list(NumberOfRecords=length(obj))
    renew(object, values=values)
})

setMethod(qa2, "QAFrequentSequence",
          function(object, state, ..., verbose=FALSE)
{
    if (verbose) message("qa2,QAFrequentSequence-method")

    if (is.finite(object@n)) {
        n <- thresh <- object@n
    } else {
        n <- 10L
        thresh <- object@k
    }

    obj <- .filter(state@data, object@useFilter)
    r <- srrank(sread(obj))
    t <- tabulate(r)
    ttop <- head(order(t, decreasing=TRUE), n)
    topCount <-
        setNames(t[ttop], as.character(sread(obj)[match(ttop, r)]))

    filt <- if (is.finite(object@n)) {
        r %in% ttop
    } else r %in% which(t >= thresh)
    .filterUpdate(state@data, object@addFilter, filt)

    values <- DataFrame(Threshold=thresh, Count=sum(filt),
                        TopCount=IntegerList(topCount))
    metadata(values) <- list(NumberOfRecords=length(obj))
    renew(object, values=values)
})

setMethod(qa2, "QAReadQuality",
          function(object, state, ..., verbose=FALSE)
{
    if (verbose) message("qa2,QAReadQuality-method")
    obj <- .filter(state@data, object@useFilter)
    dens <- .qa_qdensity(quality(obj))
    values <- DataFrame(Score=dens$x, Density=dens$y)
    metadata(values) <- list(NumberOfRecords=length(obj))
    renew(object, values=values)
})

setMethod(qa2, "QANucleotideByCycle",
          function(object, state, ..., verbose=FALSE)
{
    if (verbose) message("qa2,QANucleotideByCycle-method")
    obj <- .filter(state@data, object@useFilter)
    abc <- alphabetByCycle(sread(obj))
    values <- DataFrame(Cycle=seq_len(ncol(abc))[col(abc)],
                         Base=factor(rownames(abc)[row(abc)]),
                         Count=as.vector(abc), row.names=NULL)
    metadata(values) <- list(NumberOfRecords=length(obj))
    renew(object, values=values[values$Count != 0,])
})

setMethod(qa2, "QAQualityByCycle",
          function(object, state, ..., verbose=FALSE)
{
    if (verbose) message("qa2,QAQualityByCycle-method")
    obj <- .filter(state@data, object@useFilter)
    abc <- alphabetByCycle(quality(obj))
    q <- factor(rownames(abc)[row(abc)], levels = rownames(abc))
    q0 <- 1 + 32 * is(quality, "SFastqQuality")
    values <- DataFrame(Cycle=seq_len(ncol(abc))[col(abc)],
                        Quality=q, Score=as.numeric(q) - q0,
                        Count=as.vector(abc), row.names=NULL)
    metadata(values) <- list(NumberOfRecords=length(obj))
    renew(object, values=values[values$Count != 0,])
})

.qa2_do_collate1 <-
    function(object, state, ..., verbose=FALSE)
{
    if (verbose) message("qa2,QACollate1-method")
    src <- .clone(object@src)
    srcelt <- qa2(src, verbose=verbose) # side effect -- populate seq
    elts <- endoapply(as(object, "SimpleList"), qa2, src, ...,
                      verbose=verbose)
    names(elts) <- sapply(object, class)
    renew(object, elts, src=renew(object@src, values=values(srcelt)))
}

setMethod(qa2, "QACollate",
    function(object, state, ..., verbose=FALSE)
{
    if (verbose) message("qa2,QACollate-method")

    qas <- srapply(seq_along(object@src@con), function(i, object, ...) {
        object@src@con <- object@src@con[i]
        .qa2_do_collate1(object, ...)
    }, object, ..., verbose=verbose)

    ## collapse summary
    srcValues <- do.call(rbind, lapply(qas, function(elt) {
        values(elt@src)
    }))
    srcValues[["Id"]] <- seq_along(qas)
    ncol <- ncol(srcValues)
    srcValues <- srcValues[, c(ncol, seq_len(ncol - 1L))]

    ## collect NumberOfRecords
    filtered <- as(t(sapply(qas, function(lst) {
        sapply(lst, function(elt) {
            metadata(values(elt))[["NumberOfRecords"]]
        })
    })), "DataFrame")
    filtered[["Id"]] <- srcValues[["Id"]]
    ncol <- ncol(filtered)
    filtered <- filtered[,c(ncol, seq_len(ncol - 1L))]

    ## add Id
    qas <- Map(function(elt, id) endoapply(elt, function(elt, id) {
        df <- values(elt)
        df[["Id"]] <- id
        ncol <- ncol(df)
        rotate <- c(ncol(df), seq_len(ncol - 1L))
        values(elt) <- df[,rotate]
        elt
    }, id), qas, srcValues[["Id"]])

    ## collapse
    values <- do.call(mapply, c(function(...) {
        do.call(rbind, list(...))
    }, qas))

    QA(QASources(values=srcValues), QAFiltered(values=filtered),
       do.call(SimpleList, values))
})

## report

.hwrite <-
    function(df)
{
    hwrite(as(df, "data.frame"), border=0)
}

setMethod(report, "QASources",
          function(x, ..., dest=tempfile(), type="html")
{
    plt <- dotplot(Id ~ SourceN, as(values(x), "data.frame"),
                   type = "b", pch = 20, col = .dnaCol)
    list(SAMPLE_KEY=.hwrite(values(x)),
         PPN_COUNT=.html_img(dest, "readCounts", plt))
})

setMethod(report, "QAFiltered",
          function(x, ..., dest=tempfile(), type="html")
{
    list(FILTERED=.hwrite(values(x)))
})

setMethod(report, "QAAdapterContamination",
          function(x, ..., dest=tempfile(), type="html")
{
    list(ADAPTER_CONTAMINATION=.hwrite(values(x)))
})

setMethod(report, "QANucleotideUse",
          function(x, ..., dest=tempfile(), type="html")
{
    plt <- dotplot(Id ~ Count, group = Nucleotide,
                   as(values(x), "data.frame"),
                   type = "b", pch = 20, col = .dnaCol,
                   key = list(space = "top", lines = list(col = .dnaCol), 
                     text = list(lab = levels(values(x)[["Nucleotide"]])),
                       columns = 5L))
    list(BASE_CALL_COUNT=.html_img(dest, "baseCalls", plt))
})

setMethod(report, "QAQualityUse",
          function(x, ..., dest=tempfile(), type="html")
{
    df <- as(values(x), "data.frame")
    id <- df[["Id"]]
    q <- df[["Quality"]]
    q <- factor(q, levels=levels(q)[min(as.integer(q)):max(as.integer(q))]) 
    df[["Quality"]] <- q
    df <- df[order(df$Id, df$Quality),]
    df[["Proportion"]] <-
        with(df, unlist(Map("/",
                            lapply(split(Count, Id), cumsum),
                            lapply(split(Count, Id), sum)),
                        use.names=FALSE))
    col <- colorRampPalette(brewer.pal(9, "RdYlBu"))(length(levels(q)))
    plt <- dotplot(Id ~ Proportion, group=Quality, df,
                   type = "b", pch = 20, col = col,
                   xlab="Cummulative Proportion",
                   key = list(space = "top",
                     lines = list(col = col, size=3L), 
                     text = list(lab = levels(df[["Quality"]])),
                     columns = 10L, cex=.6))
    list(QUALITY_SCORE_COUNT=.html_img(dest, "qualityCalls", plt))
})

setMethod(report, "QAReadQuality",
    function(x, ..., dest=tempfile(), type="html")
{
    df <- as(values(x), "data.frame")
    xmin <- min(df$Score)
    ymax <- max(df$Density)
    plt <- xyplot(Density ~ Score | Id, df,
                  type = "l", xlab = "Average (calibrated) base quality", 
                  ylab = "Proportion of reads", aspect = 2,
                  panel = function(..., subscripts) {
                      panel.xyplot(...)
                      lbl <- as.character(unique(df$Id[subscripts]))
                      ltext(xmin, ymax, lbl, adj = c(0, 1))
                  }, strip = FALSE)
    list(READ_QUALITY_FIGURE=.html_img(dest, "readQuality", plt))
})

setMethod(report, "QASequenceUse",
          function(x, ..., dest=tempfile(), type="html")
{
    
    df <- with(values(x), {
        nOccur <- tapply(Occurrences, Id, c)
        cumulative <- tapply(Occurrences * Reads, Id, function(elt) {
            cs <- cumsum(elt)
            (cs - cs[1] + 1)/(diff(range(cs)) + 1L)
        })
        id <- tapply(Id, Id, c)
        data.frame(Occurrences = unlist(nOccur),
                   Cumulative = unlist(cumulative), 
                   Id = unlist(id), row.names = NULL)
    })
    xmax <- log10(max(df$Occurrences))
    plt <- xyplot(Cumulative ~ log10(Occurrences) | factor(Id), df,
           xlab = expression(paste("Number of occurrences of each sequence (",
               log[10], ")", sep = "")),
           ylab = "Cumulative proportion of reads", 
           aspect = 2, panel = function(x, y, ..., subscripts, type) {
               lbl <- unique(df$Id[subscripts])
               ltext(xmax, 0.05, lbl, adj = c(1, 0))
               type <- if (1L == length(x)) "p" else "l"
               panel.xyplot(x, y, ..., type = type)
           }, strip = FALSE)
    list(SEQUENCE_USE=.html_img(dest, "sequenceUse", plt))
})

setMethod(report, "QAFrequentSequence",
          function(x, ..., dest=tempfile(), type="html")
{
    thresholdLabel <- if (is.finite(x@n)) "n" else "k"
    threshold <- as.character(if (is.finite(x@n)) x@n else x@k)
    freqseq <- if (x@reportSequences) {
        seqdf <- lapply(with(values(x), { #with() gets wrong env for .hwrite
            lapply(split(TopCount, Id), function(elt) {
                data.frame(Sequence=names(elt[[1]]),
                           Count=unname(elt[[1]]))
            })
        }), .hwrite)
        paste(Map(function(id, seq) {
            sprintf("<p>Id: %s</p>%s", id, seq)
        }, names(seqdf), seqdf), collapse="\n")
    } else ""

    df <- values(x)[, c("Id", "Count")]
    list(THRESHOLD_LABEL=thresholdLabel, THRESHOLD=threshold,
         FREQUENT_SEQUENCE_COUNT=.hwrite(df),
         FREQUENT_SEQUENCES=freqseq)
})

setMethod(report, "QANucleotideByCycle",
          function(x, ..., dest=tempfile(), type="html")
{
    df <- as(values(x), "data.frame")
    df <- df[df$Base != "N" & df$Base != "-", ]
    df$Base <- factor(df$Base)
    xmax <- max(df$Cycle)
    ymax <- log10(max(df$Count))
    plt <-
        xyplot(log10(Count) ~ as.integer(Cycle) | Id,
               group = factor(Base),
               df[order(df$Id, df$Base, df$Cycle),],
               panel = function(..., subscripts) {
                   lbl <- as.character(unique(df$Id[subscripts]))
                   ltext(xmax, ymax, lbl, adj = c(1, 1))
                   panel.xyplot(..., subscripts = subscripts)
               },
               type = "l", col = .dnaCol[1:4],
               key = list(
                 space = "top", lines = list(col = .dnaCol[1:4]),
                 text = list(lab = levels(df$Base)),
                 columns = length(levels(df$Base))),
               xlab = "Cycle", aspect = 2, strip = FALSE)
    list(CYCLE_BASE_CALL=.html_img(dest, "cycleBaseCall", plt))
})

setMethod(report, "QAQualityByCycle",
          function(x, ..., dest=tempfile(), type="html")
{
    calc_means <- function(x, y, z)
        rowsum(y * z, x)/rowsum(z, x)
    calc_quantile <- function(x, y, z, q = c(0.25, 0.5, 0.75))
        by(list(y, z), x, function(x) {
            scoreRle <- Rle(x[[1]], x[[2]])
            quantile(scoreRle, q)
        })
    df <- as(values(x), "data.frame")
    Id <- df$Id
    pal <- c("#66C2A5", "#FC8D62")
    lvlPal <- c("#F5F5F5", "black")
    rng <- range(df$Count)
    at <- seq(rng[1], rng[2], length.out = 512)
    np <- length(unique(Id))
    nrow <- ceiling(np/4)
    layout <- c(ceiling(np/nrow), nrow)
    ymin <- min(df$Score)
    plt <- xyplot(Score ~ Cycle | Id, df,
           panel = function(x, y, ..., subscripts) {
               z <- df$Count[subscripts]
               mean <- calc_means(x, y, z)
               qtiles <- calc_quantile(x, y, z)
               sxi <- sort(unique(x))
               panel.levelplot(x, y, z, subscripts = TRUE, at = at, 
                               col.regions = colorRampPalette(lvlPal))
               llines(sxi, mean, type = "l", col = pal[[1]], lwd = 1)
               llines(sxi, sapply(qtiles, "[[", 1), type = "l",
                      col = pal[[2]], lwd = 1, lty = 3)
               llines(sxi, sapply(qtiles, "[[", 2), type = "l",
                      col = pal[[2]], lwd = 1)
               llines(sxi, sapply(qtiles, "[[", 3), type = "l",
                      col = pal[[2]], lwd = 1, lty = 3)
               lbl <- as.character(unique(df$Id[subscripts]))
               ltext(1, ymin, lbl, adj = c(0, 0))
           }, ylab = "Quality Score", layout = layout, strip = FALSE)
    list(CYCLE_QUALITY=.html_img(dest, "cycleQualityCall", plt))
})

setMethod(report, "QA", function(x, ..., dest=tempfile(), type="html") {
    if (any(type != "html"))
        .throw(SRError("UserArgumentMismatch", "'type' must be 'html'"))
    dir.create(dest, recursive=TRUE)

    sections <- c(system.file("template", "QAHeader.html",
                              package="ShortRead", mustWork=TRUE),
                  x@sources@html, x@filtered@html,
                  sapply(x, slot, "html"),
                  system.file("template", "QAFooter.html",
                              package="ShortRead", mustWork=TRUE))

    values0 <- c(list(report(x@sources, dest=dest),
                      report(x@filtered, dest=dest)),
                 lapply(x, report, dest=dest))
    values <- setNames(unlist(values0, recursive=FALSE, use.names=FALSE),
                       unlist(lapply(values0, names)))

    .report_html_do(dest, sections, values, ...)
})