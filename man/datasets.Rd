% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/datasets.R
\docType{data}
\name{datasets}
\alias{datasets}
\alias{C_IC_D_06_wID}
\alias{C_IC_N_06_wID}
\alias{CD_06_Mall_wID}
\alias{CN_06_Mall_wID}
\title{Example Data Sets}
\format{
An object of class \code{data.frame} with 251 rows and 2 columns.

An object of class \code{data.frame} with 124 rows and 2 columns.

An object of class \code{data.frame} with 251 rows and 501 columns.

An object of class \code{data.frame} with 124 rows and 501 columns.
}
\usage{
C_IC_D_06_wID

C_IC_N_06_wID

CD_06_Mall_wID

CN_06_Mall_wID
}
\description{
Four Example data sets for use the pTITAN2 package.
}
\details{
Example data sets are from publicly available macroinvertebrate survey data
from California.  The data sets are broken down between the environmental
variable, in this case percent impervious cover, and macroinvertebrate data.
Separate data files are provided for each 'treatment' that is explored. In
this case, the treatments are data from either drought (dry) or normal
precipitation years in the Chaparral region of California.

\code{CN_06_Mall_wID} (Chaparral Region, Treatment = Normal) file contains
raw macroinvertebrate density data for 500 possible macroinvertebrate codes
for each taxonomic level (class, order, family, genus).
\code{CD_06_Mall_wID} has the same information for the dry treatement.

Names for the columns within the Taxonomic files are expected to have the
following structure:
\itemize{
\item{8 characters in length}
\item{characters 1 and 2 denote the class}
\item{characters 3 and 4 denote the order}
\item{characters 5 and 6 denote the family}
\item{characters 7 and 8 denote the genus}
}
If no information at a level exists, use "00" to hold the place. For example:
A code that is 'Bi000000' is the Bivalvia class, while BiVe0000 is the
Bivalvia class, Veneroida order. BiVeSh00 is the Bivalvia class, Veneroida
order, Spheriridae family. BiVeSh01 is a genus within that family.

The raw data files are provided for your use as well.  See example below for
accessing these files.
}
\examples{
head(C_IC_D_06_wID)  # Environemntal Gradient, Dry Treatment
head(C_IC_N_06_wID)  # Environemntal Gradient, Normal Treatment
head(CD_06_Mall_wID) # Taxonomic, Dry Treatment
head(CN_06_Mall_wID) # Taxonomic, Normal Treatment

# Get the paths to the raw data files
list.files(system.file("extdata", package = "pTITAN2"))

}
\keyword{datasets}
