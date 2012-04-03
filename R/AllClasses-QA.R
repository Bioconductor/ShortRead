setClass(".QA2",
         representation("VIRTUAL", ".ShortReadBase"))

## data sources

.QAData <-
    setRefClass("QAData",
         fields=list(seq="ShortReadQ", filter="logical"),
         methods=list(show=function() {
             cat(class(.self), " ")
             print(.self$seq)
             cat(sprintf("filter: %d of %d", sum(.self$filter),
                         length(.self$filter)), "\n")
         }))

setClass("QASummary",
         representation("VIRTUAL", ".QA2",
                        addFilter="ScalarLogical",
                        useFilter="ScalarLogical",
                        values="DataFrame",
                        flag="integer",
                        html="ScalarCharacter"),
         prototype=prototype(
           addFilter=mkScalar(TRUE),
           useFilter=mkScalar(TRUE)))

## Sources

setClass("QASource",
         representation("VIRTUAL", "QASummary",
                        metadata="DataFrame", data="QAData",
                        flagNSequencesRange="integer"),
         prototype=prototype(
           flagNSequencesRange=NA_integer_),
         validity=function(object) {
             msg <- NULL
             if (!(is.na(object@flagNSequencesRange) ||
                   2L == length(object@flagNSequencesRange)))
                 msg <- "'flagNSequencesRange' must be integer(2)"
             if (is.null(msg)) TRUE else msg
         })

setClass("QAFastqSource",
         representation("QASource",
                        con="character", n="ScalarInteger",
                        readerBlockSize="ScalarInteger"))

## setclass("QABamSource",
##          representation("QASource", "QASummary", src="BamFile"))

## summaries

setClass("QAFlagged", representation("QASummary"))

setClass("QAFiltered", representation("QASummary"))

setClass("QANucleotideUse", representation("QASummary"))

setClass("QAQualityUse", representation("QASummary"))

setClass("QASequenceUse", representation("QASummary"))

setClass("QAReadQuality",
         representation("QASummary",
                        flagK="ScalarNumeric",
                        flagA="ScalarInteger"))

setClass("QAAdapterContamination",
         representation("QASummary",
                        Lpattern="ScalarCharacter",
                        Rpattern="ScalarCharacter",
                        max.Lmismatch="ScalarNumeric",
                        max.Rmismatch="ScalarNumeric",
                        min.trim="ScalarInteger"))

setClass("QAFrequentSequence",
         representation("QASummary",
                        n="ScalarInteger", a="ScalarInteger",
                        flagK="ScalarNumeric",
                        reportSequences="ScalarLogical"),
         prototype=prototype(n=mkScalar(10L)),
         validity=function(object) {
             msg <- NULL
             if (is.finite(object@n) && is.finite(object@a))
                 msg <- c(msg, "only one of 'n' or 'a' can be defined")
             else if (!is.finite(object@n) && !is.finite(object@a))
                 msg <- c(msg, "one of 'n' or 'a' must be defined")
             if (is.null(msg)) TRUE else paste("\n    ", msg)
         })

setClass("QANucleotideByCycle", representation("QASummary"))

setClass("QAQualityByCycle", representation("QASummary"))

## collation

setClass("QACollate",
         representation(".QA2", "SimpleList", src="QASource"),
         prototype=prototype(
           src=new("QAFastqSource"),
           elementType="QASummary"))

setClass("QA",
         representation(".QA2", "SimpleList",
                        src="QASource",
                        filtered="QAFiltered",
                        flagged="QAFlagged"),
         prototype=prototype(
           src=new("QAFastqSource"),
           elementType="QASummary"))
