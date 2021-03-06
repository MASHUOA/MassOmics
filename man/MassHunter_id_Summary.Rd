% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/GC-MS functions.R
\name{MassHunter_id_Summary}
\alias{MassHunter_id_Summary}
\title{MassHunter report file (.cefs) parsing and summary}
\usage{
MassHunter_id_Summary(
  workdir = NULL,
  MS.L = NULL,
  MassHunter.report = NULL,
  File.name = "Summary Report.csv",
  MsLibrary = c("NIST", "InHouse"),
  Ret.Time.Filter = 2.5,
  RT.shift.limt = 60,
  mz_L = 38,
  mz_U = 550,
  generate_rt_shift_graph = F,
  generate_rt_shift_graphs = F,
  RTcorrection = F
)
}
\arguments{
\item{workdir}{locate the working Dir}

\item{MS.L}{specify the library file}

\item{MsLibrary}{specify the library file origin, could be "NIST" or "InHouse"}

\item{Ret.Time.Filter}{set the retetion time filter window around the expected retention time (in +/-min).}

\item{RT.shift.limt}{set a threshole in seconds to define a molecule have}

\item{mz_L}{lower mz limit for quantification ion selection}

\item{mz_U}{upper mz limit for quantification ion selection}

\item{generate_rt_shift_graph}{Set TRUE to generate the retetion time shift plots for each molecules that have retention time discrepencies greater than the \code{"RT.shift.limt"}}

\item{RTcorrection}{Set TRUE to enable a retention time correction before multiple peaks resolving}

\item{amdis.report}{specify the identification report file}
}
\value{
None
}
\description{
This is a function to process the AMDIS_report.txt file. specify a in-house library or a subset of NIST library to find the most intense quantification ion within the GC-MS experiment m/z range. The script will invoke selection windows to let the user specify the files location.
}
\examples{
MassHunter_id_Summary()

}
