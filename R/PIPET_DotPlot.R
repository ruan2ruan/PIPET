#' @title Dimensional reduction plot of single cell data
#' @description Graphs the output of a dimensional reduction technique on a 2D scatter plot where each point is a cell and it's
#'     positioned based on the cell embeddings determined by the reduction technique. Cells are colored by subpopulations.
#'
#' @param SC A Seurat object of single cell data.
#' @param colors A set of colors to map data values to.
#' @param color_alpha An alpha level in [0,1]. The default value is 1.
#' @param reduction Which dimensionality reduction to use. The default value is "umap".
#' @param group A character, name of one metadata columns to group cells by.
#' @param pt.size Adjust point size for plotting. The default value is 1.2.
#' @param label Whether to label the clusters. The default value is FALSE.
#' @param label.size Sets size of labels. The default value is 5
#' @param title Add a title to the plot.
#'
#' @return An output of graphic or a list of ggplot objects.
#'
#' @import Seurat ggplot2
#' @export
#'
#' @examples # Please refer to vignette.
PIPET_DotPlot <- function(SC, colors, color_alpha=1,
                          reduction = "umap", group,
                          pt.size = 1.2,label = FALSE, label.size = 5,title){

  type_colors <- c(colors[1:length(unique(SC$type))-1],"grey")

  p <- DimPlot(SC, reduction = reduction, group.by = group, pt.size = pt.size,
               label = label, label.size = label.size, seed = 2021)+
    scale_color_manual(values = alpha(type_colors,color_alpha))+
    ggtitle(title)
  p

  return(p)
}
