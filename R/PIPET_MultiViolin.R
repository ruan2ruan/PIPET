#' @title Multi-violin plot of single cell data
#' @description Draws a stacked violin plot with multiple modules of single cell data.
#'
#' @param SC A Seurat object of single cell data.
#' @param features A set of features to plot.
#' @param preprocess Select whether to perform normalization and scale of data. The default value is TRUE.
#' @param col A character, the variable name of prediction, and must be contained in the metadata columns.
#' @param palette If a string, will use that named palette. If a number, will index into the list of palettes of appropriate type. The default value is "Set1".
#' @param xlab Set the label name of the x axis.
#' @param ylab Set the label name of the y axis.
#'
#' @return An output of graphic or a list of ggplot objects.
#'
#' @import Seurat ggplot2 cowplot
#' @export
#'
#' @examples # Please refer to vignette.
PIPET_MultiViolin <- function(SC, features, preprocess=TRUE, col="prediction",
                              palette = "Set1",xlab="Gene Expression", ylab=""){
  library(Seurat)
  library(ggplot2)
  library(cowplot)

  if(preprocess){
    ## Normalizing
    SC <- NormalizeData(SC, normalization.method = "LogNormalize", scale.factor = 1e4)
    ## Find highly variable features
    SC <- FindVariableFeatures(SC, selection.method = "vst", nfeatures = 2000)
    ## Scaling the data
    SC <- ScaleData(SC2, features = rownames(SC))
  }

  data <- data.frame(t(data.frame(SC@assays$RNA@data[features,])))
  data$Cell <- colnames(SC)
  data$Var <- SC@meta.data[[col]]

  # Use melt to change data.frame format
  data1 <- reshape2::melt(data, id.vars = c("Cell","Var"), measure.vars = features,
                          variable.name = "Genes", value.name = "Exp")
  head(data1)

  p <- ggplot(data1, aes(Exp, Var, fill = Var)) +
    geom_violin(scale = "width", adjust = 1, trim = TRUE) +
    scale_x_continuous(expand = c(0, 0), labels = function(x)
      c(rep(x = "", times = length(x)-2), x[length(x) - 1], "")) +
    facet_grid(cols = vars(Genes), scales = "free")  +
    theme_cowplot(font_size = 15) +
    scale_fill_brewer(palette = palette)+
    theme(axis.text.x = element_text(size = 13, color = "black",vjust = 0.3),
          axis.text.y = element_text(size = 15, color = "black"),
          axis.title.x = element_text(size = 16),
          axis.title.y = element_text(size = 16),
          legend.text = element_text(size = 12),

          legend.position = "none",
          #              panel.spacing = unit(0, "lines"),
          plot.title = element_text(hjust = 0.5),
          #              panel.background = element_rect(fill = NA, color = "black"),
          strip.background = element_blank(),
          strip.text = element_text(face = "bold"),
          strip.text.x.top = element_text(angle = 90, hjust = 0, vjust = 0.5)) +
    xlab(xlab) + ylab(ylab)
  p

  return(p)
}
