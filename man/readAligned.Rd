\name{readAligned}

\alias{readAligned}
\alias{readAligned,character-method}

\title{Read aligned reads and their quality scores into R representations}

\description{

  \code{readAligned} reads all aligned read files in a directory
  \code{dirPath} whose file name matches \code{pattern},
  returning a compact internal representation of the alignments,
  sequences, and quality scores in the files. Methods read all files into a
  single R object; a typical use is to restrict input to a single
  aligned read file.

}
\usage{

readAligned(dirPath, pattern=character(0), ...)

}

\arguments{

  \item{dirPath}{A character vector (or other object; see methods
    defined on this generic) giving the directory path (relative or
    absolute) of aligned read files to be input.}

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
  \itemize{

    \item{\code{type="SolexaExport"}}{

      This type parses \code{.*_export.txt} files following the
      documentation in the Solexa Genome Alignment software manual,
      version 0.3.0. These files consist of the following columns;
      consult Solexa documentation for precise descriptions. If parsed,
      values can be retrieved from \code{\linkS4class{AlignedRead}} as
      follows:

      \describe{
	    \item{Machine}{Ignored}
	    \item{Run number}{stored in \code{alignData}}
	    \item{Lane}{stored in \code{alignData}}
	    \item{Tile}{stored in \code{alignData}}
	    \item{X}{stored in \code{alignData}}
	    \item{Y}{stored in \code{alignData}}
	    \item{Index string}{Ignored}
	    \item{Read number}{Ignored}
	    \item{Read}{\code{sread}}
	    \item{Quality}{\code{quality}}
	    \item{Match chromosome}{\code{chromosome}}
	    \item{Match contig}{Ignored}
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

      Paired read columns are not interpreted.  The resulting
      \code{\linkS4class{AlignedRead}} object does \emph{not} contain a
      meaningful \code{id}; instead, use information from
      \code{alignData} to identify reads.

      Different interfaces to reading alignment files are described in
      \code{\linkS4class{SolexaPath}} and \code{\linkS4class{SolexaSet}}.

    }

    \item{\code{type="SolexaPrealign"}}{See SolexaRealign}
    \item{\code{type="SolexaAlign"}}{See SolexaRealign}
    \item{\code{type="SolexaRealign"}}{

      These types parse \code{s_L_TTTT_prealign.txt},
      \code{s_L_TTTT_align.txt} or \code{s_L_TTTT_realign.txt} files
      produced by default and eland analyses. From the Solexa
      documentation, \code{align} corresponds to unfiltered first-pass
      alignements, \code{prealign} adjusts alignments for error rates
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

    \item{\code{type="MAQMap", records=-1L}}{Parse binary \code{map}
      files produced by MAQ. See details in the next section. The
      \code{records} option determines how many lines are read;
      \code{-1L} (the default) means that all records are input.}

    \item{\code{type="MAQMapShort", records=-1L}}{The same as \code{type="MAQMap"}
      but for map files made with Maq prior to version 0.7.0. (These files
      use a different maximum read length [64 instead of 128], and are hence
      incompatible with newer Maq map files.)}


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
        \item{Identifier}{ignored}
        \item{Strand}{\code{strand}}
        \item{Chromosome}{\code{chromosome}}
        \item{Position}{\code{position}; see comment below}
        \item{Read}{\code{sread}; see comment below}
        \item{Read quality}{\code{quality}; see comments below}
        \item{Bowtie reserved}{ignored}
        \item{Alignment mismatch locations}{\code{alignData}}
      }

      This method includes the argument \code{qualityType} to specify
      how quality scores are encoded.  Bowtie quality scores are
      \sQuote{Solexa}-like by default, with
      \code{qualityType='SFastqQuality'}, but can be specified as
      \sQuote{Phred}-like, with \code{qualityType='FastqQuality'}.

      Bowtie outputs positions that are 0-offset from the left-most end
      of the \code{+} strand. \code{ShortRead} parses position
      information to be 1-offset from the left-most end of the \code{+}
      strand.

      Bowtie outputs reads aligned to the \code{-} strand as their
      reverse complement, and reverses the quality score string of these
      reads. \code{ShortRead} parses these to their original sequence
      and orientation.

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

  A \code{\linkS4class{AlignedRead}} object.

  The MAQ reference manual,
  \url{http://maq.sourceforge.net/maq-manpage.shtml#5}, 3 May, 2008

  The Bowtie reference manual, \url{http://bowtie-bio.sourceforge.net},
  28 October, 2008

}

\author{
  Martin Morgan <mtmorgan@fhcrc.org>,
  Simon Anders <anders@ebi.ac.uk> (MAQ map)}

\examples{
sp <- SolexaPath(system.file("extdata", package="ShortRead"))
ap <- analysisPath(sp)
## ELAND_EXTENDED
readAligned(ap, "s_2_export.txt", "SolexaExport")
## PhageAlign
readAligned(ap, "s_5_.*_realign.txt", "SolexaRealign")

## MAQ
dirPath <- system.file('extdata', 'maq', package='ShortRead')
list.files(dirPath)
## First line
readLines(list.files(dirPath, full.names=TRUE)[[1]], 1)
countLines(dirPath)
## two files collapse into one
readAligned(dirPath, type="MAQMapview")

## select only chr1-5.fa, '+' strand
filt <- compose(chromosomeFilter("chr[1-5].fa"),
                strandFilter("+"))
readAligned(sp, "s_2_export.txt", filter=filt)
}
\keyword{manip}
