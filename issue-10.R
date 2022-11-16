n <- 14840633L
m <- 151    #Read length
dna <- paste(rep("A", m), collapse = "")
qual <- paste(rep("J", m), collapse = "")

fastq <- tempfile()
writeLines(paste0(
    "@read", 1:n, "\n",
    dna, "\n",
    "+\n",
    qual
), fastq)

devtools::load_all()
countFastq(fastq)
