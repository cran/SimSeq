\name{SimSeq-package}
\alias{SimSeq-package}
\alias{SimSeq}
\docType{package}
\title{Nonparametric Simulation of RNA-Seq Data}
\description{RNA-seq analysis methods are often derived by relying on hypothetical parametric models for read counts that are not likely to be precisely satisfied in practice.  Methods are often tested by analyzing data that have been simulated according to the assumed model.  This testing strategy can result in an overly optimistic view of the performance of an RNA-seq analysis method.  Rather than generating data from a questionable parametric model, the R package SimSeq provides tools for simulating RNA-seq data by subsampling from existing datasets.  The vector of read counts simulated for a given experimental unit has a joint distribution that closely matches the distribution of actual RNA-seq data.  Users can control the proportion of genes simulated to be differentially expressed (DE) and can adjust the magnitude of differences for DE genes. SimSeq requires a matrix of RNA-seq read counts with large sample sizes in at least two treatment groups. Although a full version of such a dataset is not distributed with the package, numerous candidate datasets are publicly available.}
\details{
\tabular{ll}{
Package: \tab SimSeq\cr
Type: \tab Package\cr
Version: \tab 1.0.0\cr
Date: \tab 2014-02-17\cr
License: \tab GPL-2\cr
}
SimSeq performs data based simulation of RNA-Seq data creating a dataset with a known list of DE and EE genes. The core function that implements of the methodology of SimSeq is the SimData function.
}
\author{
Samuel Benidt

Maintainer: Samuel Benidt <sbenidt@iastate.edu>
}