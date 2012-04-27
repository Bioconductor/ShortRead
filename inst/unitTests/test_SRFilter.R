aln <- local({
    sp <- SolexaPath(system.file("extdata", package="ShortRead"))
    readAligned(sp, "s_2_export.txt")
})

test_srFilter <- function() {
    checkTrue(validObject(srFilter()))
    checkTrue(validObject(srFilter(name="Filter")))
    checkIdentical(name(srFilter(name="Filter")), Biobase::mkScalar("Filter"))
    checkTrue(validObject(srFilter(function(x) {})))
    checkException(srFilter(function(){}), silent=TRUE)
    checkException(srFilter(function(x, ...) {}), silent=TRUE)
}


test_occurrenceFilter <- function()
{
    checkTrue(validObject(occurrenceFilter()))

    aln <- AlignedRead(DNAStringSet(character(2)),
                       chromosome=c("chr1", "chr1"),
                       position=c(1L, 1L),
                       strand=rep(strand("+"), 2))
    checkTrue(all(c(TRUE, FALSE) == occurrenceFilter(withSread=TRUE)(aln)))
    checkTrue(all(c(TRUE, FALSE) == occurrenceFilter(withSread=FALSE)(aln)))
    checkTrue(all(c(TRUE, FALSE) == occurrenceFilter(withSread=NA)(aln)))

    aln <- AlignedRead(DNAStringSet(c("A", "T")),
                       chromosome=c("chr1", "chr1"),
                       position=c(1L, 1L),
                       strand=rep(strand("+"), 2))
    checkTrue(all(c(TRUE, TRUE) == occurrenceFilter(withSread=TRUE)(aln)))
    checkTrue(all(c(TRUE, FALSE) == occurrenceFilter(withSread=FALSE)(aln)))
    checkTrue(all(c(TRUE, TRUE) == occurrenceFilter(withSread=NA)(aln)))

    aln <- AlignedRead(DNAStringSet(character(4)),
                       chromosome=rep(c("chr1", "chr2"), each=2),
                       position=rep(1:2, 2),
                       strand=rep(strand("+"), 4))
    checkTrue(all(occurrenceFilter(withSread=FALSE)(aln)))
    checkTrue(all(occurrenceFilter(withSread=TRUE)(aln)))
    checkTrue(all(c(TRUE, FALSE, FALSE, FALSE) == 
                   occurrenceFilter(withSread=NA)(aln)))

    sp <- SolexaPath(system.file("extdata", package="ShortRead"))
    aln <- readAligned(analysisPath(sp), "s_2_export.txt", "SolexaExport")
    checkIdentical(980L, sum(occurrenceFilter(withSread=NA)(aln)))
    checkIdentical(996L, sum(occurrenceFilter(withSread=TRUE)(aln)))
    df <- data.frame(chromosome(aln), position(aln), strand(aln))
    checkIdentical(sum(!duplicated(df)),
                   sum(occurrenceFilter(withSread=FALSE)(aln)))

    checkIdentical(15L,
                   sum(occurrenceFilter(min=5, max=10, withSread=NA)(aln)))

    checkIdentical(13L,
                   sum(occurrenceFilter(min=3, max=5, withSread=NA)(aln)))

    checkIdentical(8L,
                   sum(occurrenceFilter(min=3, max=5,
                                        duplicates="none",
                                        withSread=NA)(aln)))
}

test_chromosomeFilter <- function() {
    checkTrue(validObject(chromosomeFilter()))
    checkException(chromosomeFilter(c("foo", "bar")), silent=TRUE)

    chr <- "chr5.fa"
    obj <- chromosomeFilter(chr)
    checkIdentical(aln[obj(aln)], aln[grep(chr, chromosome(aln))])
}

test_strandFilter <- function() {
    checkTrue(validObject(strandFilter()))
    checkException(strandFilter(1), silent=TRUE)

    str <- character(0)
    obj <- strandFilter(str)
    checkIdentical(aln[obj(aln)], aln[FALSE])

    str <- "+"
    obj <- strandFilter(str)
    checkIdentical(aln[obj(aln)], aln[strand(aln)=="+" & !is.na(strand(aln))])

    str <- c("+", "-")
    obj <- strandFilter(str)
    checkIdentical(aln[obj(aln)],
                   aln[(strand(aln)=="+" |strand(aln)=="-") &
                       !is.na(strand(aln))])
    str <- c("+", NA)
    obj <- strandFilter(str)
    checkIdentical(aln[obj(aln)],
                   aln[strand(aln)=="+" | is.na(strand(aln))])
    
    obj <- strandFilter(NA_character_)
    checkIdentical(aln[obj(aln)], aln[is.na(strand(aln))])
}

test_nFilter <- function() {
    checkTrue(validObject(nFilter()))
    checkTrue(validObject(nFilter(20)))
    checkException(nFilter("alf"), silent=TRUE)
    checkException(nFilter(1:2), silent=TRUE)

    n <- 0
    checkIdentical(aln[nFilter(n)(aln)], clean(aln))
    n <- 30
    alf <- alphabetFrequency(sread(aln), baseOnly=TRUE)
    checkIdentical(aln[nFilter(n)(aln)], aln[alf[,"other"] <= n])
}

test_polynFilter <- function() {
    checkTrue(validObject(polynFilter()))
    checkTrue(validObject(polynFilter(20)))
    checkTrue(validObject(polynFilter(nuc=c("A", "other"))))
    checkException(polynFilter(1:2), silent=TRUE)
    checkException(polynFilter("x"), silent=TRUE)
    checkException(polynFilter(nuc="Z"), silent=TRUE)

    alf <- alphabetFrequency(sread(aln), baseOnly=TRUE)

    n <- 30
    obj <- polynFilter(n)
    checkIdentical(aln[obj(aln)], aln[apply(alf, 1, max) <= n])

    n <- 30
    obj <- polynFilter(n, c("A", "C", "T", "G"))
    checkIdentical(aln[obj(aln)],
                   aln[apply(alf[,1:4], 1, max) <= n])
}

test_dustyFilter <- function() {
    checkTrue(validObject(dustyFilter()))
    checkTrue(validObject(dustyFilter(20)))

    checkTrue(validObject(lgl0 <- dustyFilter(10L)(aln)))
    checkTrue(validObject(lgl1 <- dustyFilter(10L)(sread(aln))))
    checkIdentical(lgl0, lgl1)

    checkIdentical(lgl0, dustyFilter(10L, 100L)(aln))
    checkIdentical(lgl0, dustyFilter(10L, 100L)(sread(aln)))
}
    

test_srdistanceFilter <- function() {
    checkTrue(validObject(srdistanceFilter()))
    checkTrue(validObject(srdistanceFilter("sdf", 1)))
    checkException(srdistanceFilter(123), silent=TRUE)
    checkException(srdistanceFilter("sdfs", 1:2), silent=TRUE)

    obj <- srdistanceFilter()
    checkIdentical(aln[obj(aln)], aln)

    nr <- c("GATCGGAAGAGCTCGTATGCCGTCTTCTGCTTAGA",
            "GATCGGAAGAGCTCGTATGCCGTCTTCTGCTTGAA")

    obj <- srdistanceFilter(nr[[1]], 1L)
    checkIdentical(aln[obj(aln)],
                   aln[as.character(sread(aln))!=nr[[1]]])
    obj <- srdistanceFilter(nr, 1L)
    checkIdentical(aln[obj(aln)],
                   aln[as.character(sread(aln))!=nr[[1]] &
                       as.character(sread(aln))!=nr[[2]] ])
}

test_alignQualityFilter <- function()
{
    checkTrue(validObject(alignQualityFilter()))
    checkTrue(validObject(alignQualityFilter(70)))
    checkException(alignQualityFilter("foo"), silent=TRUE)
    checkException(alignQualityFilter(threshold=1:2), silent=TRUE)

    checkIdentical(aln[alignQualityFilter()(aln)], aln)

    n <- 70
    obj <- alignQualityFilter(n)
    checkIdentical(aln[obj(aln)],
                   aln[quality(alignQuality(aln))>=70])
}

test_alignDataFilter <- function()
{
    checkTrue(validObject(alignDataFilter()))
    ex <- expression(x>200 & y<600)
    checkTrue(validObject(alignDataFilter(ex)))

    ad <- pData(alignData(aln))
    checkIdentical(aln[alignDataFilter(ex)(aln)],
                   aln[eval(ex, ad)])
}

test_compose <- function() {
    f1 <- chromosomeFilter("chr5.fa")
    f2 <- polynFilter(12)
    checkTrue(validObject(compose()))
    checkTrue(validObject(compose(f1)))
    checkTrue(validObject(compose(f1, f2)))
    obj <- compose(f1, f2)
    checkTrue(validObject(obj))
    checkIdentical(name(obj),
                   Biobase::mkScalar(paste(name(f1), name(f2), sep=" o ")))

    checkException(compose("foo"), silent=TRUE)
    checkException(compose(f1, "foo"), silent=TRUE)

    checkIdentical(aln[compose(f1)(aln)], aln[f1(aln)])
    checkIdentical(aln[obj(aln)], aln[f1(aln) & f2(aln)])
}
