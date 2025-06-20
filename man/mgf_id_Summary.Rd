% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/GC-MS functions.R
\name{mgf_id_Summary}
\alias{mgf_id_Summary}
\title{mgf report file parsing and summary}
\usage{
mgf_id_Summary(
  workdir = getwd(),
  mgf.report = NULL,
  mgf.peak.rds = NULL,
  RT.shift.limt = 30,
  quantile_LB = 0.25,
  quantile_UB = 0.75,
  generate_rt_shift_graph = T,
  generate_rt_shift_graphs = F,
  Save_prefix = "",
  map_batch = NULL
)
}
\arguments{
\item{workdir}{locate the working Dir}

\item{RT.shift.limt}{set a threshole in seconds to define a molecule have}

\item{generate_rt_shift_graph}{Set TRUE to generate the retetion time shift plots for each molecules that have retention time discrepencies greater than the \code{"RT.shift.limt"}}

\item{MS.L}{specify the library file}

\item{MsLibrary}{specify the library file origin, could be "NIST" or "InHouse"}

\item{amdis.report}{specify the identification report file}

\item{Ret.Time.Filter}{set the retetion time filter window around the expected retention time (in +/-min).}

\item{mz_L}{lower mz limit for quantification ion selection}

\item{mz_U}{upper mz limit for quantification ion selection}

\item{RTcorrection}{Set TRUE to enable a retention time correction before multiple peaks resolving}
}
\value{
None
}
\description{
This is a function to process the AMDIS_report.txt file. specify a in-house library or a subset of NIST library to find the most intense quantification ion within the GC-MS experiment m/z range. The script will invoke selection windows to let the user specify the files location.
}
\examples{
mgf_id_Summary()

}
