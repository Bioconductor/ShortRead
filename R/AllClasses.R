## .STRAND_LEVELS needs to be early, to be used in class
## prototypes. C-level code retrieves this value. pileup and
## readAligned,type=MAQMap depend on this ordering
.STRAND_LEVELS <- levels(strand())
.toStrand_Solexa <- function(x)
    factor(.STRAND_LEVELS[match(x, c("F", "R"))],
           levels=.STRAND_LEVELS)

.srValidity <- function(object) TRUE

setGeneric(".srValidity")

## Virtual base classes

setClass(".ShortReadBase")

## .SRUtil: SRError / SRList / SRVector / SRFilter

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

setClass("SRFilter",
         contains=c("function", ".SRUtil"),
         representation=representation(
           name="ScalarCharacter"),
         validity=.srValidity)

setClass("SRFilterResult",
         contains=c("logical", ".SRUtil"),
         representation=representation(
           name="ScalarCharacter",
           stats="data.frame"))

## Intensity

setClass("IntensityMeasure", contains=".ShortReadBase",
         representation=representation("VIRTUAL"))

setClass("IntensityInfo", contains=".ShortReadBase",
         representation=representation("VIRTUAL"))

setClass("Intensity", contains=".ShortReadBase",
         representation=representation(
           .hasMeasurementError="ScalarLogical",
           readInfo="IntensityInfo",
           intensity="IntensityMeasure",
           measurementError="IntensityMeasure",
           "VIRTUAL"),
         prototype=prototype(
           .hasMeasurementError=mkScalar(FALSE)),
         validity=.srValidity)

## Intensity, implementation

setClass("ArrayIntensity",
         contains=c("array", "IntensityMeasure"),
         prototype=prototype(array(0, c(0, 0, 0))))

ArrayIntensity <-
    function(intensity=array(0, c(0, 0, 0)), ...)
{
    new("ArrayIntensity", intensity, ...)
}

setClass("SolexaIntensityInfo",
         ## AnnotatedDataFrame as prototype does not work, r46984.
         ## .init is a work-around to identify user-constructed
         ## objects that should be valid; used in .srValidity-method
         contains=c("AnnotatedDataFrame", "IntensityInfo"),
         representation=representation(.init="ScalarLogical"),
         prototype=prototype(.init=mkScalar(FALSE)),
         validity=.srValidity)

SolexaIntensityInfo <-
  function(lane=integer(0), tile=integer(0)[seq_along(lane)],
           x=integer(0)[seq_along(lane)],
           y=integer(0)[seq_along(lane)])
{
    new("SolexaIntensityInfo",
        data=data.frame(
          lane=lane, tile=tile, x=x, y=y),
        varMetadata=data.frame(
          labelDescription=c(
            "Solexa lane nubmer",
            "Solexa tile nubmer",
            "Tile x coordinate",
            "Tile y coordinate")),
        .init=mkScalar(TRUE))
}

setClass("SolexaIntensity", contains="Intensity",
         prototype=prototype(
           readInfo=SolexaIntensityInfo(),
           intensity=ArrayIntensity(),
           measurementError=ArrayIntensity()),
         validity=.srValidity)

setClass("RtaIntensity", contains="SolexaIntensity")

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
           strand=factor(levels=.STRAND_LEVELS),
           alignQuality=NumericQuality()),
         validity=.srValidity)

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

setClass(".QA", contains=c("SRList", ".ShortReadBase"),
         representation=representation("VIRTUAL"))

setClass("ShortReadQQA", contains=".QA")

setClass("FastqQA", contains="ShortReadQQA") # synonym

setClass("SolexaExportQA", contains=".QA")

setClass("SolexaRealignQA", contains=".QA")

setClass("MAQMapQA", contains=".QA")

setClass("BowtieQA", contains=".QA")

setClass("BAMQA", contains=".QA")
