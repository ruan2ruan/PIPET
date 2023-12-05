#' @title Bar plot of single cell data
#' @description Plot a bar chart and make the height of the bar proportional to the number of cells in each group.
#'
#' @param SC A Seurat object of single cell data.
#' @param group A character, name of one metadata columns to group cells by.
#' @param prediction A character, the variable name of prediction, and must be contained in the metadata columns.
#' @param palette If a string, will use that named palette. If a number, will index into the list of palettes of appropriate type. The default value is "Set1".
#' @param ylim An integer, set limits for the y axis.
#'
#' @return An output of graphic or a list of ggplot objects.
#'
#' @import Seurat ggplot2
#' @export
#'
#' @examples # Please refer to vignette.
PIPET_BarPlot <- function(SC,group,prediction,palette = "Set1",ylim){

  library(Seurat)
  library(ggplot2)

  tmp <- data.frame(Var1=SC@meta.data[[group]],Var2=SC@meta.data[[prediction]])

  if(is.factor(tmp[,"Var1"])){
    message(paste0("The levels of ", group , " are: "));print(levels(tmp[,"Var1"]))
    message(paste0("If the order of levels is not appropriate, please re-customize ",group,"."))
  }else{
    tmp[,"Var1"] <- as.factor(tmp[,"Var1"])
    message(paste0(group," has been automatically converted to factor type.",
                   "If the order of levels is not appropriate, please re-customize ",group,"."))
  }

  p <- ggplot(tmp, aes(x = Var2, group = Var1, fill = Var1)) +
    geom_bar(stat ="count",lwd = 1.5, colour = "white")+
    # scale_fill_manual(values = c("#E7B800","#8ac926"))+
    labs(x = "",y = "Number of cells")+
    # geom_text(aes(label = tmp$Freq),position=position_dodge(width = 0.8),size =4 ,vjust = -0.45)+
    # guides(fill = guide_legend(reverse = F))+
    scale_fill_brewer(palette = palette) +
    theme_classic()+
    scale_y_continuous(expand=c(0,0))+
    coord_cartesian(ylim = c(0, ylim)) +
    theme(axis.text.x = element_text(size = 13, color = "black",vjust = 0.3),
          axis.text.y = element_text(size = 12, color = "black"),
          axis.title.x = element_text(size = 16),
          axis.title.y = element_text(size = 16),
          legend.text = element_text(size = 12),
          legend.title=element_blank())
  p

  return(p)
}
