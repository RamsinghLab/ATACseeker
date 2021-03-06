% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plotMT.R
\name{plotMT}
\alias{plotMT}
\title{plot mitochondrial variant calls in the context of heteroplasmy}
\usage{
plotMT(mtCalls, filterVAF = TRUE, rot = -1, title = NULL)
}
\arguments{
\item{mtCalls}{a VRanges object with variant calls and associated VAFs}

\item{filterVAF}{filter low-quality and 0 or 1 VAF calls? (TRUE)}

\item{rot}{plot counterclockwise (-1, default) or clockwise (1)?}

\item{title}{a title for the plot (NULL)}
}
\value{
a ggplot object
}
\description{
use a plot (inspired by the original, which was written by Stephen Turner)
to visualize mitochondrial heteroplasmy from mitochondrial variant calls
}
\details{
by default, calls with a VAF of 0 or 1 are excluded, but this can be changed
}
