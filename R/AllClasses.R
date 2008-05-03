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
           quality="BStringSet"),
         prototype=prototype(
           quality=BStringSet(character(0))),
         validity=.srValidity)

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
           alignQuality="numeric",
           alignData="AlignedDataFrame"),
         prototype=prototype(
           strand=factor(levels=c("+", "-"))),
         validity=.srValidity)

## .Solexa

setClass(".Solexa", contains=".ShortReadBase",
         representation=representation("VIRTUAL"))

setClass("SolexaPath", contains=".Solexa",
         representation=representation(
           experimentPath="character",
           dataPath="character",
           scanPath="character",
           imageAnalysisPath="character",
           baseCallPath="character",
           analysisPath="character"),
         prototype=prototype(
           experimentPath=NA_character_,
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
