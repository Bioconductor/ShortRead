test_coverage_leftmost_plus <- function() 
{
    ## 'leftmost'
    ##                 1    2
    ##        8        7    2
    ##        ++++++++++-----
    ## ....|....|....|....|....|
    readlength <- 10
    aln <- new("AlignedRead",
               chromosome=factor("A"),
               position=8L,
               strand=strand("+"),
               alignQuality=new("MatrixQuality", quality=matrix(0, 1, readlength)),
               alignData=AlignedDataFrame(
                 data=data.frame(row.names=1), metadata=data.frame()),
               quality=SFastqQuality(rep(paste(rep(" ", readlength), collapse=""), 1)),
               sread=DNAStringSet(rep(paste(rep("A", readlength), collapse=""), 1)),
               id=BStringSet(rep("", 1)))
    checkIdentical(c(7L, 15L, 3L),
                   runLength(coverage(aln, 1, 25, extend=5L)[[1]]))
}

test_coverage_leftmost_minus <- function() 
{
    ## ....|....|....|....|
    ##   -----++++++++++
    ##   3    8        1
    ##                 7
    readlength <- 10
    aln <- new("AlignedRead",
               chromosome=factor("A"),
               position=8L,
               strand=strand("-"),
               alignQuality=new("MatrixQuality", quality=matrix(0, 1, readlength)),
               alignData=AlignedDataFrame(
                 data=data.frame(row.names=1), metadata=data.frame()),
               quality=SFastqQuality(rep(paste(rep(" ", readlength), collapse=""), 1)),
               sread=DNAStringSet(rep(paste(rep("A", readlength), collapse=""), 1)),
               id=BStringSet(rep("", 1)))
    checkIdentical(c(2L,15L,3L),
                   runLength(coverage(aln, 1, 20, extend=5L)[[1]]))
}

test_coverage_fiveprime_plus <- function() 
{
    ## 5'
    ##                 1    2
    ##        8        7    2
    ##        ++++++++++-----
    ## ....|....|....|....|....|
    readlength <- 10
    aln <- new("AlignedRead",
               chromosome=factor("A"),
               position=8L,
               strand=strand("+"),
               alignQuality=new("MatrixQuality", quality=matrix(0, 1, readlength)),
               alignData=AlignedDataFrame(
                 data=data.frame(row.names=1), metadata=data.frame()),
               quality=SFastqQuality(rep(paste(rep(" ", readlength), collapse=""), 1)),
               sread=DNAStringSet(rep(paste(rep("A", readlength), collapse=""), 1)),
               id=BStringSet(rep("", 1)))
    checkIdentical(c(7L, 15L, 3L),
                   runLength(coverage(aln, 1, 25, coords="fiveprime", extend=5L)[[1]]))
}

test_coverage_fiveprime_minus <- function()
{
    ## ....|....|....|....|....|
    ##   -----++++++++++
    ##   3    8        1
    ##                 7
    readlength <- 10
    aln <- new("AlignedRead",
               chromosome=factor("A"),
               position=17L,
               strand=strand("-"),
               alignQuality=new("MatrixQuality", quality=matrix(0, 1, readlength)),
               alignData=AlignedDataFrame(
                 data=data.frame(row.names=1), metadata=data.frame()),
               quality=SFastqQuality(rep(paste(rep(" ", readlength), collapse=""), 1)),
               sread=DNAStringSet(rep(paste(rep("A", readlength), collapse=""), 1)),
               id=BStringSet(rep("", 1)))
    checkIdentical(c(2L, 15L, 8L),
                   runLength(coverage(aln, 1, 25, coords="fiveprime", extend=5L)[[1]]))
}
