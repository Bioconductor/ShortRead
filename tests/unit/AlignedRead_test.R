sp <- SolexaPath(system.file("extdata", package="ShortRead"))
aln <- readAligned(sp, "s_2_export.txt")

.checkAlignedRead_identical<- function(obs, exp)
    ## can't compare external pointers
{
    checkIdentical(as.character(sread(obs)),
                   as.character(sread(exp)))
    checkIdentical(as.character(quality(quality(obs))),
                   as.character(quality(quality(exp))))
    checkIdentical(as.character(id(obs)), as.character(id(exp)))
    checkIdentical(chromosome(obs), chromosome(exp))
    checkIdentical(strand(obs), strand(exp))
    checkIdentical(alignQuality(obs), alignQuality(exp))
    checkIdentical(alignData(obs), alignData(exp))
}

test_AlignedRead_Bowtie <- function() {
    src <- system.file("extdata", "bowtie", package="ShortRead")
    df <- read.table(file.path(src, "s_1_aligned_bowtie.txt"),
                     fill=TRUE, quote="", sep="\t")
    aln <- readAligned(src, "^s_1_aligned_bowtie.txt$", "Bowtie")
    checkIdentical(nrow(df), length(aln))
    checkIdentical(as.character(df[[2]]), as.character(strand(aln)))
    checkIdentical(as.character(df[[3]]), as.character(chromosome(aln)))
    checkIdentical(df[[4]]+1L, position(aln))
    idx <- strand(aln)=="-"
    s1 <- as.character(df[[5]])
    s1[idx] <- as.character(reverseComplement(DNAStringSet(s1[idx])))
    checkIdentical(s1, as.character(sread(aln)))
    q1 <- as.character(df[[6]])
    q1[idx] <- as.character(reverse(BStringSet(q1[idx])))
    checkIdentical(q1, as.character(quality(quality(aln))))
    checkIdentical(as.character(df[[8]]),
                   as.character(alignData(aln)[["mismatch"]]))
}

test_AlignedRead_readAligned_SolexaExport <- function() {
    obj <- readAligned(analysisPath(sp),
                       pattern="s_2_export.txt", type="SolexaExport")
    checkTrue(validObject(obj))
    checkTrue(is(quality(obj), "SFastqQuality"))
    checkTrue(is(alignQuality(obj), "NumericQuality"))
    checkIdentical(varLabels(alignData(obj)),
                   c("run", "lane", "tile", "x", "y", "filtering"))
}

test_AlignedRead_readAligned_SolexaExport_filter <- function()
{
    chr <- "chr5.fa"
    filt <- chromosomeFilter(chr)
    obs <- readAligned(sp, "s_2_export.txt", filter=filt)
    exp <- aln[grep(chr, chromosome(aln))]
    .checkAlignedRead_identical(obs, exp)
    obs <- readAligned(analysisPath(sp), "s_2_export.txt",
                       "SolexaExport", filter=filt)
    .checkAlignedRead_identical(obs, exp)
}

test_AlignedRead_readAligned_MAQMapview <- function() {
    fl <- system.file("extdata", "maq", package="ShortRead")
    obj <- readAligned(fl, pattern=".*aln.*", type="MAQMapview")
    checkTrue(validObject(obj))
    checkTrue(is(quality(obj), "FastqQuality"))
    checkTrue(is(alignQuality(obj), "NumericQuality"))
    checkIdentical(varLabels(alignData(obj)),
                   c("nMismatchBestHit", "mismatchQuality",
                   "nExactMatch24", "nOneMismatch24"))
    checkIdentical(levels(chromosome(obj)), "ChrA")
    checkTrue(!any(is.na(chromosome(obj))) &&
              !any(is.null(chromosome(obj))))
}

readAligned_maq_consistent <- function() {
    ## FIXME: find adequate data to store in ShortRead pkg
    if (!file.exists("/home/jdavison/sharedrsrc/proj/ycao/data/binary_maps/s_5.map"))
        return(TRUE)

    x <- readAligned("/home/jdavison/sharedrsrc/proj/ycao/data/binary_maps",
                     "s_5.map", "MAQMap")
    y <- readAligned("/home/jdavison/sharedrsrc/proj/ycao/data/text_maps",
                     "s_5.txt", "MAQMapview")

    checkIdentical(length(x), length(y))
    checkIdentical(width(x), width(y))
    checkIdentical(as.character(chromosome(x)), as.character(chromosome(y)))

    ## FIXME: we'd really like chromosome to have identical levels,
    ## but info on levels with no mapped reads is not available in the
    ## text version
    idx <- match(levels(chromosome(y)), levels(chromosome(x)))
    checkTrue(all(!is.na(idx)) && all(diff(idx) > 0))
    checkIdentical(position(x), position(y))
    checkIdentical(strand(x), strand(y))
    checkIdentical(alignQuality(x), alignQuality(y))
    checkIdentical(alignData(x), alignData(y))

    .checkXString <- function(x, y) {
        checkIdentical(as.character(x), as.character(y))
    }
    .checkXString(sread(x), sread(y))
    .checkXString(sread(x), sread(y))
    .checkXString(quality(quality(x)), quality(quality(y)))
    .checkXString(id(x), id(y))
}

test_AlignedRead_readAligned_run_as_factor <- function()
{
    aln <- readAligned("./cases", "^s_2_export_run_as_factor.txt$",
                       "SolexaExport")
    checkIdentical(alignData(aln)[["run"]],
                   factor(rep("genome", length(aln))))
}

test_AlignedRead_readAligned_realign_targetpos <- function()
{
    ## column 4 can be target:pos
    fl <- "s_2_0001_realign_head.txt"
    tbl <- read.table(file.path("./cases", fl), fill=TRUE)
    aln <- readAligned("./cases", fl, "SolexaRealign")
    checkIdentical(as.character(tbl[[1]]), as.character(sread(aln)))
    checkIdentical(tbl[[2]], quality(alignQuality(aln)))
    checkIdentical(tbl[[3]], alignData(aln)[["nMatch"]])
    chr <- sub(":.*", "", tbl[[4]])
    chr[nchar(chr)==0] <- NA
    checkIdentical(factor(chr), chromosome(aln))
    checkIdentical(as.integer(sub(".*:", "", tbl[[4]])), position(aln))
    checkIdentical(ShortRead:::.toStrand_Solexa(tbl[[5]]), strand(aln))
    checkIdentical(tbl[[7]], alignData(aln)[["nextBestAlignQuality"]])
}

test_AlignedRead_readAligned_SolexaResult <- function()
{
    fl <- "s_1_results_head.txt"
    tbl <- read.table(file.path("./cases", fl), fill=TRUE,
                      col.names=paste("V", 1:12, sep=""))
    aln <- readAligned("./cases", fl, "SolexaResult")
    
    checkIdentical(as.character(tbl[[2]]), as.character(sread(aln)))
    chr <- tbl[[7]]
    chr[nchar(chr)==0] <- NA
    checkIdentical(factor(chr), chromosome(aln))
    checkIdentical(tbl[[8]], position(aln))
    checkIdentical(ShortRead:::.toStrand_Solexa(tbl[[9]]), strand(aln))

    ad <- alignData(aln)
    checkIdentical(tbl[[3]], ad[[1]])
    checkIdentical(tbl[[4]], ad[[2]])
    checkIdentical(tbl[[5]], ad[[3]])
    checkIdentical(tbl[[6]], ad[[4]])

    checkIdentical(tbl[[10]], ad[[5]])
    checkIdentical(tbl[[11]], ad[[6]])
    checkIdentical(tbl[[12]], ad[[7]])
}
