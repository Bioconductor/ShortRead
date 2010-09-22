\name{readAligned}

\alias{readAligned}
\alias{readAligned,character-method}

\title{Read aligned reads and their quality scores into R representations}

\description{

  Import files containing aligned reads into an internal representation
  of the alignments, sequences, and quality scores. Methods read all
  files into a single R object.

}
\usage{

readAligned(dirPath, pattern=character(0), ...)

}

\arguments{

  \item{dirPath}{A character vector (or other object; see methods
    defined on this generic) giving the directory path (relative or
    absolute; some methods also accept a character vector of file names)
    of aligned read files to be input.}

  \item{pattern}{The (\code{\link{grep}}-style) pattern describing file
    names to be read. The default (\code{character(0)}) results in
    (attempted) input of all files in the directory.}

  \item{...}{Additional arguments, used by methods. When \code{dirPath}
    is a character vector, the argument \code{type} must be
    provided. Possible values for \code{type} and their meaning are
    described below. Most methods implement \code{filter=srFilter()},
    allowing objects of \code{\linkS4class{SRFilter}} to selectively
    returns aligned reads.}

}

\details{

  There is no standard aligned read file format; methods parse
  particular file types.

  The \code{readAligned,character-method} interprets file types based
  on an additional \code{type} argument. Supported types are:
  \describe{

    \item{\code{type="SolexaExport"}}{

      This type parses \code{.*_export.txt} files following the
      documentation in the Solexa Genome Alignment software manual,
      version 0.3.0. These files consist of the following columns;
      consult Solexa documentation for precise descriptions. If parsed,
      values can be retrieved from \code{\linkS4class{AlignedRead}} as
      follows:

      \describe{
            \item{Machine}{see below}
            \item{Run number}{stored in \code{alignData}}
            \item{Lane}{stored in \code{alignData}}
            \item{Tile}{stored in \code{alignData}}
            \item{X}{stored in \code{alignData}}
            \item{Y}{stored in \code{alignData}}
            \item{Multiplex index}{see below}
            \item{Paired read number}{see below}
            \item{Read}{\code{sread}}
            \item{Quality}{\code{quality}}
            \item{Match chromosome}{\code{chromosome}}
            \item{Match contig}{\code{alignData}}
            \item{Match position}{\code{position}}
            \item{Match strand}{\code{strand}}
            \item{Match description}{Ignored}
            \item{Single-read alignment score}{\code{alignQuality}}
            \item{Paired-read alignment score}{Ignored}
            \item{Partner chromosome}{Ignored}
            \item{Partner contig}{Ignored}
            \item{Partner offset}{Ignored}
            \item{Partner strand}{Ignored}
            \item{Filtering}{\code{alignData}}
      }
	  
	  The following optional arguments, set to \code{FALSE} by default,
	  influence data input

      \describe{
        
        \item{withMultiplexIndex}{When \code{TRUE}, include the
          multiplex index as a column \code{multiplexIndex} in
          \code{alignData}.}

        \item{withPairedReadNumber}{When \code{TRUE}, include the paired
          read number as a column \code{pairedReadNumber} in
          \code{alignData}.}

        \item{withId}{When \code{TRUE}, construct an identifier string
          as
          \sQuote{Machine_Run:Lane:Tile:X:Y#multiplexIndex/pairedReadNumber}. The
          substrings \sQuote{#multiplexIndex} and
          \sQuote{/pairedReadNumber} are not present if
          \code{withMultiplexIndex=FALSE} or
          \code{withPairedReadNumber=FALSE}.}

        \item{withAll}{A convencience which, when \code{TRUE}, sets all
          \code{with*} values to \code{TRUE}.}

      }

      Note that not all paired read columns are interpreted.  Different
      interfaces to reading alignment files are described in
      \code{\linkS4class{SolexaPath}} and
      \code{\linkS4class{SolexaSet}}.

    }

    \item{\code{type="SolexaPrealign"}}{See SolexaRealign}
    \item{\code{type="SolexaAlign"}}{See SolexaRealign}
    \item{\code{type="SolexaRealign"}}{

      These types parse \code{s_L_TTTT_prealign.txt},
      \code{s_L_TTTT_align.txt} or \code{s_L_TTTT_realign.txt} files
      produced by default and eland analyses. From the Solexa
      documentation, \code{align} corresponds to unfiltered first-pass
      alignments, \code{prealign} adjusts alignments for error rates
      (when available), \code{realign} filters alignments to exclude
      clusters failing to pass quality criteria.

      Because base quality scores are not stored with alignments, the
      object returned by \code{readAligned} scores all base qualities as
      \code{-32}.

      If parsed, values can be retrieved from
      \code{\linkS4class{AlignedRead}} as follows:

      \describe{
        \item{Sequence}{stored in \code{sread}}
        \item{Best score}{stored in \code{alignQuality}}
        \item{Number of hits}{stored in \code{alignData}}
        \item{Target position}{stored in \code{position}}
        \item{Strand}{stored in \code{strand}}
        \item{Target sequence}{Ignored; parse using
          \code{\link{readXStringColumns}}}
        \item{Next best score}{stored in \code{alignData}}
      }
    }

    \item{\code{type="SolexaResult"}}{

      This parses \code{s_L_eland_results.txt} files, an intermediate
      format that does not contain read or alignment quality
      scores.

      Because base quality scores are not stored with alignments, the
      object returned by \code{readAligned} scores all base qualities as
      \code{-32}.

      Columns of this file type can be retrieved from
      \code{\linkS4class{AlignedRead}} as follows (description of
      columns is from Table 19, Genome Analyzer Pipeline Software User
      Guide, Revision A, January 2008):

      \describe{
        \item{Id}{Not parsed}
        \item{Sequence}{stored in \code{sread}}
        \item{Type of match code}{Stored in \code{alignData} as
          \code{matchCode}. Codes are (from the Eland manual): NM (no
          match); QC (no match due to quality control failure); RM (no
          match due to repeat masking); U0 (best match was unique and
          exact); U1 (best match was unique, with 1 mismatch); U2 (best
          match was unique, with 2 mismatches); R0 (multiple exact
          matches found); R1 (multiple 1 mismatch matches found, no
          exact matches); R2 (multiple 2 mismatch matches found, no
          exact or 1-mismatch matches).}
        \item{Number of exact matches}{stored in \code{alignData} as
          \code{nExactMatch}}
        \item{Number of 1-error mismatches}{stored in \code{alignData}
          as \code{nOneMismatch}}
        \item{Number of 2-error mismatches}{stored in \code{alignData}
          as \code{nTwoMismatch}}
        \item{Genome file of match}{stored in \code{chromosome}}
        \item{Position}{stored in \code{position}}
        \item{Strand}{(direction of match) stored in \code{strand}}
        \item{\sQuote{N} treatment}{stored in \code{alignData}, as
          \code{NCharacterTreatment}. \sQuote{.} indicates treatment of
          \sQuote{N} was not applicable; \sQuote{D} indicates treatment
          as deletion; \sQuote{|} indicates treatment as insertion}
        \item{Substitution error}{stored in \code{alignData} as
          \code{mismatchDetailOne} and \code{mismatchDetailTwo}. Present
          only for unique inexact matches at one or two
          positions. Position and type of first substitution error,
          e.g., 11A represents 11 matches with 12th base an A in
          reference but not read. The reference manual cited below lists
          only one field (\code{mismatchDetailOne}), but two are present
          in files seen in the wild.}

      }
    }

    \item{\code{type="MAQMap", records=-1L}}{Parse binary \code{map}
      files produced by MAQ. See details in the next section. The
      \code{records} option determines how many lines are read;
      \code{-1L} (the default) means that all records are input.}

    \item{\code{type="MAQMapShort", records=-1L}}{The same as
      \code{type="MAQMap"} but for map files made with Maq prior to
      version 0.7.0. (These files use a different maximum read length
      [64 instead of 128], and are hence incompatible with newer Maq map
      files.)}


    \item{\code{type="MAQMapview"}}{

      Parse alignment files created by MAQ's \sQuote{mapiew} command.
      Interpretation of columns is based on the description in the MAQ
      manual, specifically

      \preformatted{
        ...each line consists of read name, chromosome, position,
        strand, insert size from the outer coordinates of a pair,
        paired flag, mapping quality, single-end mapping quality,
        alternative mapping quality, number of mismatches of the
        best hit, sum of qualities of mismatched bases of the best
        hit, number of 0-mismatch hits of the first 24bp, number
        of 1-mismatch hits of the first 24bp on the reference,
        length of the read, read sequence and its quality.
      }

      The read name, read sequence, and quality are read as
      \code{XStringSet} objects. Chromosome and strand are read as
      \code{factor}s.  Position is \code{numeric}, while mapping quality is
      \code{numeric}. These fields are mapped to their corresponding
      representation in \code{AlignedRead} objects.

      Number of mismatches of the best hit, sum of qualities of mismatched
      bases of the best hit, number of 0-mismatch hits of the first 24bp,
      number of 1-mismatch hits of the first 24bp are represented in the
      \code{AlignedRead} object as components of \code{alignData}.

      Remaining fields are currently ignored.

    }

    \item{\code{type="Bowtie"}}{

      Parse alignment files created with the Bowtie alignment
      algorithm. Parsed columns can be retrieved from
      \code{\linkS4class{AlignedRead}} as follows:

      \describe{
        \item{Identifier}{\code{id}}
        \item{Strand}{\code{strand}}
        \item{Chromosome}{\code{chromosome}}
        \item{Position}{\code{position}; see comment below}
        \item{Read}{\code{sread}; see comment below}
        \item{Read quality}{\code{quality}; see comments below}
        \item{Similar alignments}{\code{alignData}, \sQuote{similar}
        column; Bowtie v. 0.9.9.3 (12 May, 2009) documents this as
        the number of other instances where the same read aligns against the
        same reference characters as were aligned against in this
        alignment. Previous versions marked this as \sQuote{Reserved}}
        \item{Alignment mismatch locations}{\code{alignData}
        \sQuote{mismatch}, column}

      }
	  
      NOTE: the default quality encoding changes to \code{FastqQuality}
      with \pkg{ShortRead} version 1.3.24.

      This method includes the argument \code{qualityType} to specify
      how quality scores are encoded.  Bowtie quality scores are
      \sQuote{Phred}-like by default, with
      \code{qualityType='FastqQuality'}, but can be specified as
      \sQuote{Solexa}-like, with \code{qualityType='SFastqQuality'}.

      Bowtie outputs positions that are 0-offset from the left-most end
      of the \code{+} strand. \code{ShortRead} parses position
      information to be 1-offset from the left-most end of the \code{+}
      strand.

      Bowtie outputs reads aligned to the \code{-} strand as their
      reverse complement, and reverses the quality score string of these
      reads. \code{ShortRead} parses these to their original sequence
      and orientation.

    }

    \item{\code{type="SOAP"}}{

      Parse alignment files created with the SOAP alignment
      algorithm. Parsed columns can be retrieved from
      \code{\linkS4class{AlignedRead}} as follows:

      \describe{
        \item{id}{\code{id}}
        \item{seq}{\code{sread}; see comment below}
        \item{qual}{\code{quality}; see comment below}
        \item{number of hits}{\code{alignData}}
        \item{a/b}{\code{alignData} (\code{pairedEnd})}
        \item{length}{\code{alignData} (\code{alignedLength})}
        \item{+/-}{\code{strand}}
        \item{chr}{\code{chromosome}}
        \item{location}{\code{position}; see comment below}
        \item{types}{\code{alignData} (\code{typeOfHit}: integer
          portion; \code{hitDetail}: text portion)}
      }

      This method includes the argument \code{qualityType} to specify
      how quality scores are encoded.  It is unclear from SOAP
      documentation what the quality score is; the default is
      \sQuote{Solexa}-like, with \code{qualityType='SFastqQuality'}, but
      can be specified as \sQuote{Phred}-like, with
      \code{qualityType='FastqQuality'}.

      SOAP outputs positions that are 1-offset from the left-most end of
      the \code{+} strand. \code{ShortRead} preserves this
      representation.

      SOAP reads aligned to the \code{-} strand are reported by SOAP as
      their reverse complement, with the quality string of these reads
      reversed. \code{ShortRead} parses these to their original sequence
      and orientation.

    }

    \item{\code{type="BAM"}}{

      Parse BAM files produced by samtools and other third party
      programs.  This method includes the argument
      \code{param=ScanBamParam()}. The \code{param} argument can
	  be a single ScanBamParam object or a list of ScanBamParam
	  objects. The number of ScanBamParam objects supplied in the
      \code{param} list argument must match the number of files 
	  identified in \code{dirPath}. If a single ScanBamParam is 
	  provided but there are multiple files in \code{dirPath} the 
	  ScanBamParam object will be recycled and applied to all files.  
	  The \code{which} and \code{flag}
      arguments to \code{ScanBamParam()} can be used to influence which
      reads in the BAM file are parsed; see \code{\link{ScanBamParam}}.
      The following values override user settings (issuing a warning if
      contradictory values are provided):

      \describe{

        \item{\code{simpleCigar=TRUE}}{Reads aligned with indels are
          ignored; this is required for representation in
          \code{AlignedRead}.}

        \item{\code{reverseComplement=TRUE}}{By default, BAM stores
          reads as they are aligned to the reference genome, whereas
          \code{AiignedRead} stores them as they are prior to alignment;
          this flag converts reads from the BAM to \code{AlignedRead}
          format.}

        \item{\code{what=c("qname", "flag", "rname", "strand", "pos",
            "mapq", "seq", "qual")}}{These BAM fields are mapped to
          corresponding fields in \code{AlignedRead}.}

      }

      BAM fields are mapped to \code{AlignedRead} as:

      \describe{

        \item{qname}{\code{id}}
        \item{seq}{\code{sread}}
        \item{qual}{\code{quality}}
        \item{strand}{\code{strand}}
        \item{rname}{\code{chromosome}}
        \item{pos}{\code{position}}
        \item{mapq}{\code{alignQuality}}
        \item{flag}{\code{alignData}}

      }
    }
  }
}

\value{

  A single R object (e.g., \code{\linkS4class{AlignedRead}}) containing
  alignments, sequences and qualities of all files in \code{dirPath}
  matching \code{pattern}. There is no guarantee of order in which files
  are read.

}

\seealso{

  The \code{\linkS4class{AlignedRead}} class.

  Genome Analyzer Pipeline Software User Guide, Revision A, January
  2008.

  The MAQ reference manual,
  \url{http://maq.sourceforge.net/maq-manpage.shtml#5}, 3 May, 2008.

  The Bowtie reference manual, \url{http://bowtie-bio.sourceforge.net},
  28 October, 2008.

  The SOAP reference manual, \url{http://soap.genomics.org.cn/soap1},
  16 December, 2008.

  The BAM file format specification,
  \url{http://samtools.sourceforge.net}.

}

\author{
  Martin Morgan <mtmorgan@fhcrc.org>,
  Simon Anders <anders@ebi.ac.uk> (MAQ map)}

\examples{
sp <- SolexaPath(system.file("extdata", package="ShortRead"))
ap <- analysisPath(sp)
## ELAND_EXTENDED
(aln0 <- readAligned(ap, "s_2_export.txt", "SolexaExport"))
## PhageAlign
(aln1 <- readAligned(ap, "s_5_.*_realign.txt", "SolexaRealign"))

## MAQ
dirPath <- system.file('extdata', 'maq', package='ShortRead')
list.files(dirPath)
## First line
readLines(list.files(dirPath, full.names=TRUE)[[1]], 1)
countLines(dirPath)
## two files collapse into one
(aln2 <- readAligned(dirPath, type="MAQMapview"))

## select only chr1-5.fa, '+' strand
filt <- compose(chromosomeFilter("chr[1-5].fa"),
                strandFilter("+"))
(aln3 <- readAligned(sp, "s_2_export.txt", filter=filt))
}
\keyword{manip}
