\name{Snapshot-class}

\Rdversion{1.1}
\docType{class}

\alias{Snapshot-class}
\alias{trellis-class}
% Constructor:
\alias{Snapshot}
\alias{Snapshot,character,GRanges-method}
\alias{Snapshot,character,missing-method}
\alias{Snapshot,BamFileList,GRanges-method}

% Accessors:
\alias{files}
\alias{files,Snapshot-method}
\alias{functions}
\alias{functions,Snapshot-method}
\alias{show,Snapshot-method}
\alias{view}
\alias{view,Snapshot-method}
\alias{vrange}
\alias{vrange,Snapshot-method}
\alias{annTrack}
\alias{annTrack,Snapshot-method}
\alias{fac}
\alias{fac,Snapshot-method}
\alias{getTrellis}
\alias{getTrellis,Snapshot-method}
\alias{ignore.strand}
\alias{ignore.strand,Snapshot-method}

% methods: 
\alias{pan}
\alias{pan,Snapshot-method}
\alias{togglefun}
\alias{togglefun,Snapshot-method}
\alias{togglep}
\alias{togglep,Snapshot-method}
\alias{togglez}
\alias{togglez,Snapshot-method}
\alias{zoom}
\alias{zoom,Snapshot-method}

\title{Class \code{"Snapshot"}}

\description{

  A \code{\linkS4class{Snapshot}}-class to visualize genomic data from 
  BAM files with zoom and pan functionalities.

}

\usage{
Snapshot(files, range, ...)
}

\arguments{

  \item{files}{A character() or \code{BamFileList} specifying the file(s)
    to be visualized.} 

  \item{range}{A \code{\link{GRanges}} object specifing the range to be
    visualized.}

   \item{...}{Additional, optional, arguments to be passed to the
     Snapshot \code{initialize} function. Arguments include:

     \describe{

       \item{functions:}{A \code{\link{SnapshotFunctionList}} of
	 functions, in addition to built-in \sQuote{fine_coverage},
	 \sQuote{coarse_coverage}, \sQuote{multifine_coverage}, to be
	 used for visualization.}

       \item{currentFunction:}{character(1) naming the function, from
	 \code{functinos} to be used for data input and
	 visualization. The default chooses a function based on the
	 scale at which the data is being visualized.}

       \item{annTrack:}{Annotation track. If built-in visulization
         functions are to be used, \code{annTrack} should be a
         \code{GRanges} instance and the first column of its
         elementMeatdata would be used to annotate the range.}

       \item{fac:}{Character(1) indicating which factor used for
          grouping the sample files. The factor should be included in
          the elementMetadata of \code{files}, otherwise ignored. Used
          only to visualize multiple files.
	}
	
       \item{.auto_display:}{logical(1) indicating whether the
	 visualization is to be updated when \code{show} is invoked.}

       \item{.debug}{logical(1) indicating whether debug messages are to
	 be printed.}
   }}
}




\section{Methods}{
  \describe{
    \item{zoom}{\code{signature(x = "Snapshot")}: Zoom (in or out) the
      current plot. }
    
    \item{pan}{\code{signature(x = "Snapshot")}: Pan (right or left)
      the current plot. }
    
    \item{togglefun}{\code{signature(x = "Snapshot")}: Toggle the current
       functions which imported records are to be immediately
       evaluated. Note that the active range will be changed to the
       current active window.}

    \item{togglep}{\code{signature(x = "Snapshot")}: Toggle the panning
      effects.}
    
    \item{togglez}{\code{signature(x = "Snapshot")}: Toggle the zooming 
       effects.} 
  }
}

\section{Accessors}{
   \describe{
     \item{show}{\code{signature(object = "Snapshot")}: Display a
       \code{Snapshot} object. }

     \item{files}{\code{signature(x = "Snapshot")}: Get the \code{files}
       field (object of class \code{BamFileList}) of a \code{Snapshot} object.}

     \item{functions}{\code{signature(x = "Snapshot")}: Get the
       \code{functions} field (object of \code{SnapshotFunctionList}) of a
       \code{Snapshot} object.}

     \item{view}{\code{signature(x = "Snapshot")}: Get the \code{view}
       field (object of \code{SpTrellis}) of a \code{Snapshot} object.}

     \item{vrange}{\code{signature(x = "Snapshot")}: Get the
       \code{.range} field (object of \code{GRanges}) of a \code{Snapshot}
       object. }

     \item{getTrellis}{\code{signature(x = "Snapshot")}: Get the
       \code{trellis} object, a field of the \code{SpTrellis} object.}
  }
}

\section{Fields}{
  \describe{
    \item{\code{.debug}:}{Object of class \code{function} to display
    messages while in debug mode }
    \item{\code{.auto_display}:}{Object of class \code{logical} to
       automatically disply the coverage plot. }
    \item{\code{.range}:}{Object of class \code{GRanges} indicating
      which ranges of records to be imported from BAM fiels. }
    \item{\code{.zin}:}{Object of class \code{logical} indicating
      whether the current zooming effect is zoom in. }
    \item{\code{.pright}:}{Object of class \code{logical} indicating
      whether the current panning effect is right. }
    \item{\code{.data}:}{Object of class \code{data.frame} containing
      coverage a position is represented for each strand and BAM file.}
    \item{\code{.data_dirty}:}{Object of class \code{logical}
      indicating whether to re-evaluate the imported records. }
    \item{\code{.initial_functions}:}{Object of class 
       \code{SnapshotFunctionList} available by the \code{Snapshot} object. }
    \item{\code{.current_function}:}{Object of class \code{character}
      of the function the imported recored are currently evaluated and
      visualized.}
    \item{\code{annTrack}:}{Default to \code{NULL} if not intended to
       visulize the annotation track. If default visulization
       function(s) is intended to be used to plot the annotation,
       \code{annTrack} has to be a \code{GRanges} instance.}   
    \item{\code{functions}:}{Object of class
      \code{SnapshotFunctionList} of customized functions to evaluate
       and visualize the imported records. }
    \item{\code{files}:}{Object of class \code{BamFileList} to be imported. }
    \item{\code{view}:}{Object of class \code{SpTrellis} that is
    essentially a reference class wrapper of \code{Trellis} objects.  }
  }
}

\section{Class-Based Methods}{
  \describe{
    \item{\code{display()}:}{Display the current \code{Snapshot} object. }
    \item{\code{pan()}:}{Pan (right or left) the current plot. }
    \item{\code{zoom()}:}{Zoom (in or out) the current plot. }
    % \item{\code{set_range(range)}:}{ ~~ }
    \item{\code{toggle(zoom, pan, currentFunction)}:}{Toggle zooming,
    panning effects or the currentFuction in which the imported
    records are to be evaluated and visualized.}
    % \item{\code{initialize(..., functions, currentFunction, .range, .auto_display, .debug)}:}{ ~~ }
  }
}


\author{Martin Morgan and Chao-Jen Wong \email{cwon2@fhcrc.org}}

\seealso{\code{\link{SpTrellis}}}

\examples{

## example 1: Importing specific ranges of records

file <- system.file("extdata", "SRR002051.chrI-V.bam",
                    package="yeastNagalakshmi")
which <-  GRanges("chrI", IRanges(1, 2e5))
s <- Snapshot(file, range=which)

## methods
zoom(s) # zoom in
## zoom in to a specific region
zoom(s, range=GRanges("chrI", IRanges(7e4, 7e4+8000)))
pan(s)  # pan right
togglez(s) # change effect of zooming
zoom(s) # zoom out
togglep(s) # chage effect of panning
pan(s)

## accessors
functions(s)
vrange(s)
show(s)
ignore.strand(s)
view(s) ## extract the spTrellis object
getTrellis(s) ## extract the trellis object

## example 2: ignore strand
s <- Snapshot(file, range=which, ignore.strand=TRUE)

##
## example 3: visulizing annotation track
##

library(GenomicFeatures)

getAnnGR <- function(txdb, which) {
    ex <- exonsBy(txdb, by="gene")
    sn <- keepSeqlevels(ex, seqlevels(which))
    r <- range(sn)
    gr <- unlist(r)
    values(gr)[["gene_id"]] <- rep.int(names(r), times=elementLengths(r))
    gr
}

txdbFile <- system.file("extdata", "sacCer2_sgdGene.sqlite",
                    package="yeastNagalakshmi")
# txdb <- makeTranscriptDbFromUCSC(genome="sacCer2",
#                                  tablename="sgdGene")
txdb <- loadFeatures(txdbFile)
which <-  GRanges("chrI", IRanges(1, 2e5))
gr <- getAnnGR(txdb, which)
## note that the first column of the elementMetadata annotates of the
## range of the elements.
gr

s <- Snapshot(file, range=which, annTrack=gr)
annTrack(s)
## zoom in to an interesting region
zoom(s, range=GRanges("chrI", IRanges(7e4, 7e4+8000)))

togglez(s) ## zoom out
zoom(s)

pan(s)

## example 4, 5, 6: muliple BAM files with 'multicoarse_covarage'
## and 'mutifine_coverage' view.

## Resolution does not automatically switch for views of multiple
## files. It is important to note if width(which) < 10,000, use
## multifine_coverage.  Otherwise use multicoarse_coverage
file <- system.file("extdata", "SRR002051.chrI-V.bam",
                    package="yeastNagalakshmi")
which <-  GRanges("chrI", IRanges(1, 2e5))
s <- Snapshot(c(file, file), range=which, 
              currentFunction="multicoarse_coverage")

## grouping files and view by 'multicoarse_coverage'
bfiles <- BamFileList(c(a=file, b=file))
values(bfiles) <- DataFrame(sampleGroup=factor(c("normal", "tumor")))
values(bfiles)
s <- Snapshot(bfiles, range=which, 
              currentFunction="multicoarse_coverage", fac="sampleGroup")

## grouping files and view by 'multifine_coverage'
which <- GRanges("chrI", IRanges(7e4, 7e4+8000))
s <- Snapshot(bfiles, range=which, 
              currentFunction="multifine_coverage", fac="sampleGroup")

}

\keyword{classes}
