% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/imgExperimentClass.R
\name{calc}
\alias{calc}
\title{Determine Cell Cycle, filter for Cell/Nucleus ratio and fibroblasts}
\usage{
calc(
  obj,
  fun = NULL,
  varFB = "Intensity_MedianIntensity_DNA.nucl",
  varCC = "Intensity_IntegratedIntensity_DNA.nucl",
  xMinGlobMaxFB = 0.2,
  xMaxGlobMaxFB = 0.6,
  ...
)
}
\arguments{
\item{obj}{imgExp instance}

\item{fun}{function to calculate, can be cc (cell cycle) or nc
(nucler/cell ratio) or fb (gibrobaslten). defaults to NULL 
(all calculations are performed)}
}
\description{
Calculates derived variables
}
