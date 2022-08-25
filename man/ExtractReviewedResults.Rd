\name{ExtractReviewedResults}
\alias{ExtractReviewedResults}
\title{
A function to extracted called ABSOLUTE data
}
\description{
Given a human reviewed ABSOLUTE summary file, will extract the desired calls
}
\usage{
ExtractReviewedResults(reviewed.pp.calls.fn, analyst.id, modes.fn, 
                                  out.dir.base, obj.name, verbose=FALSE)
}
\arguments{
  \item{reviewed.pp.calls.fn}{
  A filename pointing to the output of a \code{CreateReviewObject} call, potentially with
  an additional column added to the left side, with the number of desired calls if one
  chooses to override the ABSOLUTE called solution
} 
  \item{analyst.id}{
  An ID value representing the user who reviewed the data
}
  \item{modes.fn}{
  A filename pointing to the \code{.data.RData} file output by ABSOLUTE
}
\item{out.dir.base}{
  The base directory to write output files. Output will be placed in \code{out.dir.base/reviewed}
}
\item{obj.name}{
  A descriptive name for the group of samples being called
}
\item{verbose}{
  When \code{TRUE} will provide more descriptive output
}
}
\value{
As a side effect, will write a file named \code{NAME}.segtab.txt where \code{NAME} an array name. 
}
\author{
Scott L. Carter
}
\seealso{\code{\link{CreateReviewObject}}}
\examples{
}
\keyword{ utilities }

