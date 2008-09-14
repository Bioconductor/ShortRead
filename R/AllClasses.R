.srValidity <- function(object) TRUE

setGeneric(".srValidity")

## Virtual base classes

setClass(".ShortReadBase")

## .SRUtil: SRError / SRList / SRVector

setClass(".SRUtil",
         representation=representation("VIRTUAL"))

setClass("SRError", contains=".SRUtil",
         representation=representation(
           .type="character",
           .message="character"),
         prototype=prototype(
           .type="Unspecified",
           .message="unknown error"),
         validity=.srValidity)

setClass("SRWarn", contains=".SRUtil",
         representation=representation(
           .type="character",
           .message="character"),
         prototype=prototype(
           .type="Unspecified",
           .message="unknown warning"),
         validity=.srValidity)

setClass("SRList", contains=".SRUtil",
         representation=representation(
           .srlist="list"),
         prototype=prototype(
           .srlist=list()))

setClass("SRVector", contains="SRList",
         representation=representation(
           vclass="character"),
         prototype=prototype(
           vclass=NA_character_),
         validity=.srValidity)

## QualityScore

setClass("QualityScore", contains=".ShortReadBase",
         representation=representation("VIRTUAL"))

setClass("NumericQuality", contains="QualityScore",
         representation=representation(
           quality="numeric"))

NumericQuality <- function(quality=numeric(0)) { # used below
    new("NumericQuality", quality=quality)
}

setClass("IntegerQuality", contains="NumericQuality",
         representation=representation(
           quality="integer"))

setClass("MatrixQuality", contains="QualityScore",
         representation=representation(
           quality="matrix"))

setClass("FastqQuality", contains="QualityScore",
         representation=representation(
           quality="BStringSet"),
         prototype=prototype(
           quality=BStringSet(character(0))))

setClass("SFastqQuality", contains="FastqQuality") # Solexa variant

## ShortRead / ShortReadQ

setClass("ShortRead", contains=".ShortReadBase",
         representation=representation(
           sread="DNAStringSet",
           id="BStringSet"),
         prototype=prototype(
           sread=DNAStringSet(character(0)),
           id=BStringSet(character(0))),
         validity=.srValidity)

setClass("ShortReadQ", contains="ShortRead",
         representation=representation(
           quality="QualityScore"),
         prototype=prototype(
           quality=NumericQuality()),
         validity=.srValidity)

## ExperimentPath (base class for experimental data paths)

setClass("ExperimentPath", contains = c(".ShortReadBase"),
         representation = representation(
           basePath="character"),
         prototype = prototype(
           basePath=NA_character_),
         validity = .srValidity)

## SRSet (base class for datasets)

setClass("SRSet", contains = ".ShortReadBase",
         representation = representation(
           sourcePath="ExperimentPath", # for lazy loading
           readIndex="integer", # for tracking subsets and sorting
           readCount="integer", # counts of reads in each sample
           phenoData="AnnotatedDataFrame", # experimental design
           readData="AnnotatedDataFrame"), # arbitrary read annotations
         prototype = prototype(
           sourcePath=new("ExperimentPath"),
           readIndex=integer(0),
           readCount=integer(0),
           phenoData=new("AnnotatedDataFrame"),
           readData=new("AnnotatedDataFrame")),
         validity = .srValidity)

## AlignedRead: AlignedDataFrame

setClass("AlignedDataFrame", contains="AnnotatedDataFrame",
         prototype=prototype(
           new("AnnotatedDataFrame",
               dimLabels=c("readName", "alignColumn"))),
         validity=.srValidity)

setClass("AlignedRead", contains="ShortReadQ",
         representation=representation(
           chromosome="factor",
           position="integer",
           strand="factor",
           alignQuality="QualityScore",
           alignData="AlignedDataFrame"),
         prototype=prototype(
           strand=factor(levels=c("+", "-")),
           alignQuality=NumericQuality()),
         validity=.srValidity)

## represents an alignment scored against multiple targets
## not sure if this will be useful
## setClass("MultiAlignedRead", contains="AlignedRead",
##          representation=representation(
##            alignQuality="MatrixQuality"),
##          prototype=prototype(
##            alignQuality=MatrixQuality()),
##          validity=.srValidity)

## .Solexa

setClass(".Solexa", contains=".ShortReadBase",
         representation=representation("VIRTUAL"))

setClass("SolexaPath", contains=c("ExperimentPath", ".Solexa"),
         representation=representation(
           dataPath="character",
           scanPath="character",
           imageAnalysisPath="character",
           baseCallPath="character",
           analysisPath="character"),
         prototype=prototype(
           scanPath=NA_character_,
           dataPath=NA_character_,
           imageAnalysisPath=NA_character_,
           baseCallPath=NA_character_,
           analysisPath=NA_character_),
         validity=.srValidity)

setClass("SolexaSet", contains=".Solexa",
         representation=representation(
           solexaPath="SolexaPath",
           laneDescription="AnnotatedDataFrame"),
         prototype=prototype(
           solexaPath=new("SolexaPath"),
           laneDescription=new("AnnotatedDataFrame",
             data=data.frame(1:8)[,FALSE],
             dimLabels=c("laneNames", "laneColumns"))),
         validity=.srValidity)

## setClass("SolexaTile", contains=".Solexa",
##          representation=representation(
##            srfile="character",
##            lane="integer",
##            tile="integer"),
##          validity=.srValidity)

## setClass("SolexaQAVector", contains=".Solexa",
##          prototype=prototype(
##            .class="SolexaTile"),
##          validity=.srValidity)


### .Roche

setClass(".Roche", contains=".ShortReadBase",
         representation=representation("VIRTUAL"))

setClass("RochePath", contains=c("ExperimentPath", ".Roche"),
         representation=representation(
           readPath="character",
           qualPath="character"),
         prototype=prototype(
           readPath=NA_character_,
           qualPath=NA_character_),
         validity=.srValidity)

setClass("RocheSet", contains=c("SRSet", ".Roche"),
         representation=representation(
           sourcePath="RochePath"),
         prototype=prototype(
           sourcePath=new("RochePath")),
         validity=.srValidity)

## QA

setClass(".QA", contains=".ShortReadBase",
         representation=representation("VIRTUAL"))

setClass("SolexaExportQA", contains=c("SRList", ".QA"))
