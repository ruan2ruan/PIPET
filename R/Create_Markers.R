#' @title Create markers on bulk data
#' @description Based on the differential genes of bulk data (differential expression analysis was performed by DESeq2),
#'     the feature vector of each subclass is established. The bulk data to be analyzed should have at least two subclasses.
#'     You could customize the levels of subclasses, otherwise the subclasses will be automatically converted into a factor variable.
#'
#' @param bulk_data A matrix of non-negative integers with row features and sample columns.
#' @param colData A data.frame with at least a single column. Rows of colData correspond to columns of bulk_data.
#' @param class_col The variable name of subclasses, and must be contained in the colData columns.
#' @param lg2FC In the DESeq differential expression analysis results, the cutoff value of lg2FC. The default value is 1.
#' @param p.adjust In the DESeq differential expression analysis results, the cutoff value of adjust P. The default value is 0.05.
#' @param show_lg2FC Select whether to show log2 fold changes. The default value is TRUE.
#'
#' @return This function returns a \code{data.frame} with rows are samples and the columns contain the following attributes:
#'     \item{genes}{Differential genes screened out based on lg2FC and p.adjust.}
#'     \item{class}{The subclass of bulk data to which the marker gene belongs.}
#'
#' @import DESeq2 tidyverse
#' @export
#'
#' @examples # Please refer to vignette.
Create_Markers <- function(bulk_data,colData,class_col=NULL,
                           lg2FC=1,p.adjust=0.05,show_lg2FC=TRUE){
  library(DESeq2)

  if(!(class_col %in% colnames(colData)))
    stop("Please check, class_col must be a column in colData.")

  if(is.factor(colData[,class_col])){
    message(paste0("The classification is:"));print(unique(colData[,class_col]))
    message(paste0("If the order of levels is not appropriate, please re-customize ",class_col,"."))
  }else{
    colData[,class_col] <- as.factor(colData[,class_col])
    message(paste0(class_col," has been automatically converted to factor type.",
                   "If the order of levels is not appropriate, please re-customize ",class_col,"."))
  }

  colData$class <- colData[,class_col]
  c <- levels(colData$class)

  if(length(c)<2)
    stop("Please check, there are at least two classes in class_col.")

  if(length(c)==2){
    dds <- DESeqDataSetFromMatrix(countData=bulk_data, colData=colData, design= ~ class)
    dds <- DESeq(dds)
    res <- as.data.frame(results(dds))

    UP <- res %>%
      filter(padj<p.adjust) %>%
      filter(log2FoldChange>=lg2FC)
    DOWN <- res %>%
      filter(padj<p.adjust) %>%
      filter(log2FoldChange<=-lg2FC)

    if(show_lg2FC){
      UP <- UP %>% mutate(class=c[2],genes=rownames(UP))%>%
        select(genes,class,log2FoldChange)
      DOWN <- DOWN %>% mutate(class=c[1],genes=rownames(DOWN)) %>%
        select(genes,class,log2FoldChange)
    }else{
      UP <- UP %>% mutate(class=c[2],genes=rownames(UP))%>%
        select(genes,class)
      DOWN <- DOWN %>% mutate(class=c[1],genes=rownames(DOWN)) %>%
        select(genes,class)
    }

    markers <- rbind(DOWN,UP)

  }else{

    markers <- data.frame(genes=NULL,class=NULL)
    for (i in 1:length(c)) {
      print(paste0("==========",c[i],"=========="))

      colData$class_new <- ifelse(colData$class==c[i],c[i],"Others")
      colData$class_new <- factor(colData$class_new, levels = c("Others",c[i]))

      dds <- DESeqDataSetFromMatrix(countData=bulk_data, colData=colData, design= ~ class_new)
      dds <- DESeq(dds)
      res <- as.data.frame(results(dds))

      UP <- res %>%
        filter(padj<p.adjust) %>%
        filter(log2FoldChange>=lg2FC)

      if(show_lg2FC){
        UP <- UP %>%
          mutate(class=c[i],genes=rownames(UP))%>%
          select(genes,class,log2FoldChange)
      }else{
        UP <- UP %>%
          mutate(class=c[i],genes=rownames(UP))%>%
          select(genes,class)
      }

      markers <- rbind(markers,UP)
    }
  }

  return(markers)
}
