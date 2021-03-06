% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/imgExperimentClass.R
\docType{class}
\name{imageExp-class}
\alias{imageExp-class}
\alias{imgExp}
\title{An S4 class representing a single experiment}
\description{
An S4 class representing a single experiment
}
\section{Slots}{

\describe{
\item{\code{expID}}{Experiment ID}

\item{\code{experiment}}{list, each element corresponds to a
single run; each run contains 1 to n folders as list}

\item{\code{data}}{data.frame containing the sinlge measurements;
altered during filtering steps}

\item{\code{dataCC}}{result of the cellCycleFractIntegrDNAInt() function,
which determines the cell cycle}

\item{\code{dataNC}}{result of the ncAreaFilter() function, used 
to filter cells by cell/nuleus area ratio}

\item{\code{dataFB}}{result of the cellCycleFractIntegrDNAInt(), used 
to filter fibroblasts}

\item{\code{dataCutoff}}{TODO}

\item{\code{dataAgg}}{per well aggregated data}

\item{\code{ddataZ}}{z-transformed dataAgg data}

\item{\code{log}}{event-log of the different processing steps}
}}

