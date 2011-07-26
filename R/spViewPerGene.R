spViewPerGene <- function(GRL,
                          name, files, #ann.by=c("exon", "transcript"),
                          multi.lines=FALSE, ...) 
{
    if (!(name %in% names(GRL)))
        stop(sprintf("name '%s' element does not exist", name))
    
    gr <- GRL[[name]]
    gr <- keepSeqlevels(gr, as.character(seqnames(gr)@values))
    which <- reduce(range(gr))
    annTrack <- gr
    
    if (multi.lines & (length(files)>1)) {
        if (width(which) <= 10000)
            currentFunction="multifine_coverage"
        else
            currentFunction="multicoarse_coverage"
        Snapshot(files=files, range=which, annTrack=annTrack,
                 currentFunction=currentFunction, ...)
    } else
        Snapshot(files=files, range=which, annTrack=annTrack, ...)
}
