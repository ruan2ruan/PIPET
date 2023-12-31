% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Group_PIPET.R
\name{Group_PIPET}
\alias{Group_PIPET}
\title{Grouping to predict subpopulations in single-cell data}
\usage{
Group_PIPET(
  Seurat_data,
  markers,
  gene_col = NULL,
  class_col = NULL,
  group = NULL,
  rm_NA = TRUE,
  freq_counts = NULL,
  normalize = TRUE,
  scale = TRUE,
  nPerm = 1000,
  distance = "cosine",
  nCores = 4,
  nslices = 4
)
}
\arguments{
\item{Seurat_data}{A Seurat object of single cell data.}

\item{markers}{A data frame of phenotypic information from bulk data.}

\item{gene_col}{A character, the variable name of genes, and must be contained in the markers columns.}

\item{class_col}{A character, the variable name of subclasses, and must be contained in the markers columns.}

\item{group}{A character, name of one metadata column to group cells by (for example, orig.ident).}

\item{rm_NA}{Select Whether to remove NA values. The default value is TRUE.}

\item{freq_counts}{An integer, keep genes expressed in more than a certain number of cells.}

\item{normalize}{Select whether to perform normalization of count data. The default value is TRUE.}

\item{scale}{Select whether to scale and center features in the dataset. The default value is TRUE.}

\item{nPerm}{An integer, number of permutations to do. The default value is 1000.}

\item{distance}{A character, the distance algorithm must be included in "cosine", "pearson", "spearman", "kendall","euclidean","maximum".}

\item{nCores}{Set the number of cores. The default value is 4.}

\item{nslices}{It specifies the number of slices to divide the input data into for parallel processing. The default value is 4.}
}
\value{
This function returns a \code{data.frame} with rows are cells and the columns contain the following attributes:
\item{prediction}{Subpopulation labels corresponding to single cell data determined based on distance.}
\item{dist_}{Distance or similarity for each subclasses (determined by the chosen distance algorithm).}
\item{Pvalue}{A nominal p-value estimated based on a null distribution for the distance  to the feature vectors.}
\item{FDR}{Adjusted p values with false discovery rate.}
}
\description{
Make group predictions for relevant subpopulations in single-cell data from phenotypic information in bulk data.
}
\examples{
# Similar to function PIPET.
}
\references{
Hoshida, Y. (2010). Nearest Template Prediction: A Single-Sample-Based Flexible Class Prediction with Confidence Assessment. PLoS ONE 5, e15543.
}
