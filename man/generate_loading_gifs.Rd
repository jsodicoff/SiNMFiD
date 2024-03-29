% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/3d_plotting.R
\name{generate_loading_gifs}
\alias{generate_loading_gifs}
\title{Generate gifs of cell type distributions derived from deconvolution in space}
\usage{
generate_loading_gifs(
  filepath,
  analysis.name,
  spatial.data.name,
  rand.seed = 123,
  mat.use = "proportions",
  cell.types.plot = NULL,
  filter = NULL,
  dims = c(500, 500)
)
}
\arguments{
\item{filepath}{Path to analysis directory}

\item{analysis.name}{String identifying the analysis}

\item{spatial.data.name}{String identifying the spatial sample}

\item{rand.seed}{Integer random seed}

\item{mat.use}{A string, either "raw" or "proportions"
referring to what version of the results to summarize}

\item{cell.types.plot}{A character vector of cell types to plot}

\item{dims}{Integer vector of length 2 corresponding to the width and
height of the RGL window}
}
\value{
nothing
}
\description{
Generate gifs of cell type distributions derived from deconvolution in space
}
