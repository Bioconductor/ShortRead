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
    cvg <- coverage(aln, 1, 25, extend=5L)
    checkIdentical(c(7L, 15L, 3L), runLength(cvg[[1]]))
}

test_coverage_leftmost_minus <- function() 
{
    ## ....|....|....|....|
    ##   -----++++++++++
    ##   3    8        1
    ##                 7
    aln <- .mkAln(8L, .width, strand("-"))
    cvg <- coverage(aln, 1, 20, extend=5L)
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
    cvg <- coverage(aln, 1, 25, coords="fiveprime", extend=5L)
    checkIdentical(c(7L, 15L, 3L), runLength(cvg[[1]]))
}

test_coverage_fiveprime_minus <- function()
{
    ## ....|....|....|....|....|
    ##   -----++++++++++
    ##   3    8        1
    ##                 7
    aln <- .mkAln(17L, .width, strand("-"))
    cvg <- coverage(aln, 1, 25, coords="fiveprime", extend=5L)
    checkIdentical(c(2L, 15L, 8L), runLength(cvg[[1]]))
}
