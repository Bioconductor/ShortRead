.width <- 10
.mkAln <- function(position, width, strand)
{
    n <- length(position)
    AlignedRead(sread=DNAStringSet(rep(polyn("A", width), n)),
                chromosome=rep("ChrA", n),
                position=as.integer(position),
                strand=strand)
}                                    

test_coverage_leftmost_plus <- function() 
{
    ## 'leftmost'
    ##                 1    2
    ##        8        7    2
    ##        ++++++++++-----
    ## ....|....|....|....|....|
    aln <- .mkAln(8L, .width, strand("+"))
    cvg <- coverage(aln, width=c(ChrA=25L), extend=5L)
    checkIdentical(c(7L, 15L, 3L), runLength(cvg[[1]]))
}

test_coverage_leftmost_minus <- function() 
{
    ## ....|....|....|....|
    ##   -----++++++++++
    ##   3    8        1
    ##                 7
    aln <- .mkAln(8L, .width, strand("-"))
    cvg <- coverage(aln, width=c(ChrA=20L), extend=5L)
    checkIdentical(c(2L,15L,3L), runLength(cvg[[1]]))
}

test_coverage_fiveprime_plus <- function() 
{
    ## 5'
    ##                 1    2
    ##        8        7    2
    ##        ++++++++++-----
    ## ....|....|....|....|....|
    aln <- .mkAln(8L, .width, strand("+"))
    cvg <- coverage(aln, width=c(ChrA=25L), coords="fiveprime", extend=5L)
    checkIdentical(c(7L, 15L, 3L), runLength(cvg[[1]]))
}

test_coverage_fiveprime_minus <- function()
{
    ## ....|....|....|....|....|
    ##   -----++++++++++
    ##   3    8        1
    ##                 7
    aln <- .mkAln(17L, .width, strand("-"))
    cvg <- coverage(aln, width=c(ChrA=25L), coords="fiveprime", extend=5L)
    checkIdentical(c(2L, 15L, 8L), runLength(cvg[[1]]))
}

test_coverage_width_names <- function()
{
    aln <- .mkAln(1, 10, strand("+"))
    checkTrue(validObject(coverage(aln)))
    checkTrue(validObject(coverage(aln, width=c(ChrA=20L))))
    ## no names on width
    checkException(coverage(aln, width=100), silent=TRUE)
    ## wrong name on width
    checkException(coverage(aln, width=c(ChrB=100)), silent=TRUE)
    ## extra width element -- ok
    checkTrue(validObject(coverage(aln, width=c(ChrA=20L, ChrB=20L))))
}
