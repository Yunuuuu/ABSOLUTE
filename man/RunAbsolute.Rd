\name{RunAbsolute}
\alias{RunAbsolute}
\alias{train.obj}
\alias{chr.arms.dat}
\alias{kDiseaseMap}
\title{
ABSOLUTE: Absolute Quantification of Somatic DNA Alterations in Human Cancer
}
\description{
When DNA is extracted from an admixed population of cancer and normal cells, the information on absolute copy number per cancer cell is lost in the mixing.  The purpose of ABSOLUTE is to re-extract these data from the mixed DNA population.  This process begins by generation of segmented copy number data, which is input to the ABSOLUTE algorithm together with pre-computed models of recurrent cancer karyotypes and, optionally, allelic fraction values for somatic point mutations.  The output of ABSOLUTE then provides re-extracted information on the absolute cellular copy number of local DNA segments and, for point mutations, the number of mutated alleles.
}
\usage{
RunAbsolute(seg.dat.fn, sigma.p, max.sigma.h, min.ploidy, max.ploidy, 
            primary.disease, platform, sample.name, results.dir,
            max.as.seg.count, max.non.clonal, max.neg.genome,
            copy_num_type, maf.fn = NULL, min.mut.af = NULL, 
            output.fn.base = NULL, verbose = FALSE)
}
\arguments{
  \item{seg.dat.fn}{
    A filename which contains the output object from a call to \code{RunHapSeg} from the \code{HAPSEG} package.
  }
  \item{sigma.p}{
    Provisional value of excess sample level variance used for mode search 
  }
  \item{max.sigma.h}{
  Maximum value of excess sample level variance (Eq. 6)
  }
  \item{min.ploidy}{
  Minimum ploidy value to consider. Solutions implying lower ploidy values will be discarded
  }
  \item{max.ploidy}{
  Maximum ploidy value to consider. Solutions implying greater ploidy values will be discarded
  }
  \item{primary.disease}{
    Primary disease of the sample
  }
  \item{platform}{
    The chip type used, supported values are currently \code{SNP_250K_STY} and \code{SNP_6.0}	
  }
  \item{sample.name}{
    The name of the sample, included in output plots
  }
  \item{max.as.seg.count}{
  Maximum number of allelic segments.  Samples with a higher segment count will be flagged as 'failed'
  }
  \item{max.neg.genome}{
     Maximum genome fraction that may be modeled as non-clonal with copy-ratio below that of clonal homozygous deletion.  Solutions implying greater values will be discarded
  }
  \item{max.non.clonal}{
     Maximum genome fraction that may be modeled as non-clonal (subclonal SCNA).  Solutions implying greater values will be discarded
  }
  \item{copy_num_type}{
	The type of copy number to be handled. Either \code{total} or \code{allelic}
  }
  \item{min.mut.af}{
  Minimum mutation allelic fraction.  Mutations with lower allelic fractions will be filtered out before analysis
  }
  \item{maf.fn}{
    Minor allele frequency file, or \code{NULL} if one is not available. This specifies the data for somatic point mutations to be used by ABSOLUTE 
  }
  \item{results.dir}{
    Directory to place output files, defaults to the current directory.
  }
  \item{output.fn.base}{
    Base filename for output files, see \code{value} below
  }
  \item{verbose}{
    If \code{TRUE}, will display a verbose description of the running of this function
  }
}
\value{
  \code{RunAbsolute} simply returns \code{TRUE} on success, but will write a
  combination of files as a side effect, all within \code{results.dir}.

  The main output of \code{RunAbsolute} is a file \code{ofb.ABSOLUTE.RData}, where
  \code{ofb} is the \code{output.fn.base} argument. This \code{RData} file contains
  a list object named \code{seg.dat}.

  The \code{seg.dat} list has the following members:
    \itemize{
      \item{seg.info}{ A \code{matrix} containing the segmentation information}
      \item{matched.name}{ A \code{matrix} name of paired normal SNP array, if used}
      \item{found.matched.normal}{ \code{TRUE} if the matched normal sample was found}
      \item{allele.segs}{ A \code{matrix} containing the genomic segmentation information, including the HSCR values, if available}
      \item{as.seg.dat}{ A \code{matrix} containing the information in allele.segs, but with the two allelic values for each seg as separate rows}
      \item{error.model}{ A \code{list} containing estimated values of sample level error-model properties. These are used to model the uncertainty of the segmental observations}
      \item{array.name}{ The \code{array.name} from the original \code{seg.dat}}
      \item{normal}{ The \code{normal} argument that was passed in}
      \item{sample.name}{ The \code{sample.name} argument that was passed in}
      \item{primary.disease}{ The \code{primary.disease} argument that was passed in}
      \item{group}{ Used to select a karyotype model specific to the disease being analyzed}
      \item{platform}{ The \code{platform} argument that was passed in}
      \item{obs.scna}{ A \code{list} containing the SCNA copy ratios as used for modelling by ABSOLUTE}
      \item{mode.res}{ A \code{list} containing the internal representation of the modeled data for all purity/ploidy modes identified}  
    }

  This file can be read into \R using the \code{\link{load}} function.

  A PDF is also generated withe the name \code{ofb.ABSOLUTE.plot.pdf} where
  \code{ofb} is the \code{output.fn.base} argument. 
}
\author{
Scott L. Carter
}
\seealso{
 \code{\link{CreateReviewObject}}
}
\examples{
}
\keyword{ manip }
