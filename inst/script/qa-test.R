## fastq

nLanes <- 50

library(ShortRead)
sp <- SolexaPath(system.file('extdata', package='ShortRead'))
rfq <- readFastq(analysisPath(sp), pattern="s_1_sequence.txt")
nms <- sprintf("file_%d.fastq", seq_len(nLanes))
qas <- lapply(nms, qa, dirPath=rfq)
qa <- do.call(rbind, qas)
res <- browseURL(report(qa))
