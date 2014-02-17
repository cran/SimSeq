\name{kidney}
\alias{kidney}
\docType{data}
\title{
Kidney Renal Clear Cell Carcinoma [KIRC] RNA-Seq data
}
\description{
A list of RNA-seq data from The Cancer Genome Atlas, subsetted down to 1000 genes and 69 paired columns of data. Provides a matrix of paired data with rows corresponding to genes and columns corresponding to replicates; replic vector specifies replicates and treatment vector specifies non-tumor and tumor group samples respectively within replicate.
}
\usage{data(kidney)}
\format{
List containing:
\itemize{
\item counts: matrix of RNA-seq data for 1000 sampled genes and 69 paired columns from individuals with Kidney Renal Clear Cell Carcinoma.
\item replic: vector detailing which column in counts matrix belongs to each individual.
\item treatment: vector detailing whether each column in counts matrix is a non-tumor or tumor sample.
}
}
\source{
https://tcga-data.nci.nih.gov/tcga/

The Cancer Genome Atlas Research Network (2013). Comprehensive molecular
characterization of clear cell renal cell carcinoma. Nature, 499(7456), 43-49.
}
\examples{
data(kidney)
}
\keyword{datasets}
