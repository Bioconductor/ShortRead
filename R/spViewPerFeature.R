.checkClass <- function(x, class, length=NULL)
{
   msg <- paste("'", substitute(x), "' must be object of class ",
                "'", class, "'", sep="")
   fail <- !any(sapply(class, function(c, y) is(y, c), x))
   if (!is.null(length) && length(x) != length) {
       fail=TRUE
       msg <- paste(msg, "of length", length)
   }
   if (fail) stop(msg) else invisible()
   
}

spViewPerFeature <- function(GRL,
                          name, files, #ann.by=c("exon", "transcript"),
                          ignore.strand=FALSE,
                          multi.levels=FALSE,
                          fac=character(0L), ...) 
{   
    .checkClass(GRL, "GRangesList")
    .checkClass(name, "character", 1)
    .checkClass(multi.levels, "logical", 1)
    .checkClass(files, c("character", "BamFileList"))
    .checkClass(ignore.strand, "logical", 1)
    .checkClass(fac, "character")
    
    if (!(name %in% names(GRL)))
        stop(sprintf("name '%s' element does not exist", name))
    
    gr <- GRL[[name]]
    gr <- keepSeqlevels(gr, as.character(seqnames(gr)@values))
    which <- reduce(range(gr))
    annTrack <- gr
    
    if (multi.levels & (length(files)>1)) {
        if (width(which) <= 10000)
            currentFunction="multifine_coverage"
        else
            currentFunction="multicoarse_coverage"
        Snapshot(..., files=files, range=which, annTrack=annTrack, fac=fac,
                 currentFunction=currentFunction, ignore.strand=ignore.strand)
    } else ## sigle file 
        Snapshot(..., files=files, range=which, annTrack=annTrack,
                 ignore.strand=ignore.strand, fac=fac)
}
