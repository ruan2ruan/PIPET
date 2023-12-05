#' @title Predict subpopulations in single-cell data
#' @description Predicting relevant subpopulations in single-cell data from phenotypic information in bulk data.
#'
#' @param SC_data A numeric matrix/data.frame of single-cell data with row features and sample columns or a Seurat object.
#' @param markers A data frame of phenotypic information from bulk data.
#' @param gene_col A character, the variable name of genes, and must be contained in the markers columns.
#' @param class_col A character, the variable name of subclasses, and must be contained in the markers columns.
#' @param rm_NA Select Whether to remove NA values. The default value is TRUE.
#' @param freq_counts An integer, keep genes expressed in more than a certain number of cells.
#' @param normalize Select whether to perform normalization of count data. The default value is TRUE.
#' @param scale Select whether to scale and center features in the dataset. The default value is TRUE.
#' @param nPerm An integer, number of permutations to do. The default value is 1000.
#' @param distance A character, the distance algorithm must be included in "cosine", "pearson", "spearman", "kendall","euclidean","maximum".
#' @param nCores Set the number of cores. The default value is 4.
#'
#' @return  This function returns a \code{data.frame} with rows are cells and the columns contain the following attributes:
#'     \item{prediction}{Subpopulation labels corresponding to single cell data determined based on distance.}
#'     \item{dist_}{Distance or similarity for each subclasses (determined by the chosen distance algorithm).}
#'     \item{Pvalue}{A nominal p-value estimated based on a null distribution for the distance  to the feature vectors.}
#'     \item{FDR}{Adjusted p values with false discovery rate.}
#'
#' @import tibble Seurat parallel stats matrixStats Matrix
#' @export
#'
#' @examples # Please refer to vignette.
#' @references Hoshida, Y. (2010). Nearest Template Prediction: A Single-Sample-Based Flexible Class Prediction with Confidence Assessment. PLoS ONE 5, e15543.
PIPET <- function(SC_data, markers, gene_col=NULL, class_col=NULL,
                  rm_NA=TRUE, freq_counts=NULL, normalize=TRUE, scale=TRUE,
                  nPerm = 1000, distance = "cosine", nCores=4) {
  library(tibble)
  library(Seurat)
  library(parallel)
  library(stats)
  library(matrixStats)
  library(Matrix)

  if(gene_col != "genes" & class_col != "class"){
    markers$genes <- markers[,gene_col]
    markers$class <- markers[,class_col]
  }
  if(gene_col != "genes" & class_col == "class"){
    markers$genes <- markers[,gene_col]
  }
  if(gene_col == "genes" & class_col != "class"){
    markers$class <- markers[,class_col]
  }

  #### Check SC data
  if (class(SC_data) == "Seurat"){
    SC <- as.matrix(SC_data@assays[["RNA"]]@counts)
  } else {
    SC <- as.matrix(SC_data)
  }
  if (has_rownames(SC))
    stop("The row names of the SC matrix must be gene names.")

  #### Check markers
  if (!is.data.frame(markers))
    stop("The data type of markers must be data frame.")
  if (is.null(markers$genes) | is.null(markers$class))
    stop("markers missing gene column or class column.")
  if (length(markers$genes) > length(unique(markers$genes)))
    stop("There cannot be duplicate genes in markers data.")
  if (!is.factor(markers$class))
    markers$class <- as.factor(markers$class)

  #### Preprocess SC data
  ## remove gene rows with missing values
  if(rm_NA){
    index_NA <- stats::complete.cases(SC)
    if (sum(!index_NA) > 0)
      SC <- SC[index_NA,,drop = FALSE]
  }
  ## keep those genes expressed in more than 'freq_counts' cells
  if(!is.null(freq_counts)){
    nonzero <- SC > 0
    keep_genes <- Matrix::rowSums(nonzero) >= freq_counts
    SC <- SC[keep_genes, ]
  }

  #### Preprocess markers
  keep_gene <- markers$genes %in% rownames(SC)
  if (sum(keep_gene) > 0) {
    message(paste0(sum(keep_gene)," genes are in the final feature vector."))
    markers <- markers[keep_gene,]
    message(paste0("The classification of markers is:"));print(table(markers$class))
  } else {
    stop("No overlapping genes, please check the markers data.")
  }

  if (min(table(markers$class ))<2)
    stop("The number of features in a class is less than 2, please rebuild the feature vector.",call. = FALSE)
  if (min(table(markers$class ))<5)
    warning("The number of features in a class is less than 5, the predictions may be unstable.",call.= FALSE)

  #### Preprocessing
  n_cells <- ncol(SC)
  n_levels <- nlevels(markers$class)
  n_markers <- nrow(markers)
  n_SC <- nrow(SC)

  class_names <- levels(markers$class)
  markers$class <- as.numeric(markers$class)

  ## Normalizing and scaling SC data
  SC[is.na(SC)] <- 0
  if(normalize){
    SC <- log1p(sweep(SC,2,Matrix::colSums(SC),FUN = "/")*1e4)
  }
  if(scale){
    SC <- t(scale(t(SC), center=TRUE, scale=TRUE))
  }
  SC[is.na(SC)] <- 0

  if (!is.matrix(SC))
    SC <- as.matrix(SC)

  ## Match vector for SC and markers
  mm <- match(markers$genes, rownames(SC),nomatch = 0)

  if (!all(rownames(SC)[mm] == markers$genes)) {
    stop("No matching genes, please check SC data and markers data.")
  }

  ## Prepare templates
  M_mat <- matrix(rep(markers$class,n_levels), ncol = n_levels) # markers matrix
  for (i in seq_len(n_levels))
    M_mat[,i] <- as.numeric(M_mat[,i] == i)
  if (n_levels == 2) M_mat[M_mat==0] <- -1

  #### Prediction function
  pred_fun <- function(n) {

    if (!distance %in% c("cosine", "pearson", "spearman", "kendall","euclidean","maximum"))
      stop("This distance algorithm is not included.")

    DistToCor  <- function(x) 1-2*x^2
    CorToDist <- function(x) sqrt(1/2*(1-(x)))

    corCosine <- function(x, y) {
      x  <- as.matrix(x);y <- as.matrix(y)
      crossprod(x,y) /
        outer(sqrt(matrixStats::colSums2(x^2)), sqrt(matrixStats::colSums2(y^2)))
    }

    if(distance %in% c("cosine", "pearson", "spearman", "kendall")){
      if (distance == "cosine") {
        corFun <- function(x,y) corCosine(x,y)
      } else {
        corFun <- function(x, y) {
          stats::cor(x, y, method = distance)
        }
      }

      # sample-markers correlations
      cor <- as.vector(corFun(SC[mm,n, drop = FALSE],M_mat))

      cor.perm.max <- matrixStats::rowMaxs(corFun(
        matrix(SC[,n][sample.int(n_SC, n_markers*nPerm, replace=TRUE)],
               ncol = nPerm), M_mat))

      pred <- which.max(cor)
      # estimate p-value
      cor.ranks <- rank(-c(cor[pred],(cor.perm.max)))
      pval <- cor.ranks[1]/length(cor.ranks)

      dist <- CorToDist(cor)
    }

    if(distance %in% c("euclidean","maximum")){
      disFun <- function(x,y) {
        tmp <- dist(t(cbind(x,y)),method = distance)
        tmp1 <- as.vector(tmp)[1:2]
        return(tmp1)
      }

      dist <- as.vector(disFun(SC[mm,n, drop = FALSE],M_mat))
      # sample-markers correlations
      cor <- DistToCor(dist)

      cor.perm.max <- matrixStats::rowMaxs(DistToCor(
        t(apply(matrix(SC[,n][sample.int(n_SC, n_markers*nPerm, replace=TRUE)],
                       ncol = nPerm),2,disFun,y=M_mat))))

      pred <- which.min(dist)
      # estimate p-value
      cor.ranks <- rank(-c(cor[pred],(cor.perm.max)))
      pval <- cor.ranks[1]/length(cor.ranks)
    }

    return(c(
      pred,               # prediction
      dist,               # distance
      pval))              # p-value
  }

  #### Parallel computing
  message(paste0("Calculating......"))

  library(parallel)
  nn <- makeCluster(getOption("cl.cores", nCores))
  clusterEvalQ(nn, { set.seed(123) })
  res <- parLapply(nn,seq_len(n_cells), pred_fun)
  stopCluster(nn)

  res <- data.frame(do.call(rbind,res))

  #### Output
  colnames(res) <- c("prediction",paste0("dist_",class_names),"Pvalue")
  res$prediction <- factor(class_names[res$prediction], levels = class_names)
  rownames(res) <- colnames(SC)
  res$FDR <- stats::p.adjust(res$Pvalue, "fdr")

  message(paste0("Finished!"))
  return(res)
}
