setwd("D:/Ruan/bulk and single-cell/simulation")

################################## Three Class ########################################
library("splatter")
library("scater")
library("VariantAnnotation")
library("ggplot2")

set.seed(123)
vcf <- mockVCF()
gff <- mockGFF(n.genes = 5000)
params.cond <- newSplatPopParams(eqtl.n = 0.3, 
                                 batchCells = 1000,
                                 similarity.scale = 5,
                                 condition.prob = c(0.4, 0.4, 0.2),
                                 eqtl.condition.specific = 0.5,
                                 cde.facLoc = 0.1, 
                                 cde.facScale = 0.2)

sim.pop.cond <- splatPopSimulate(vcf = vcf, gff = gff, params = params.cond, 
                                 sparsify = FALSE)

sim.pop.cond <- logNormCounts(sim.pop.cond)
sim.pop.cond <- runPCA(sim.pop.cond, ncomponents = 10)
plotPCA(sim.pop.cond, colour_by = "Condition", shape_by = "Sample")

################################################################################
## Prepare Single cell exp and metadata
count0 <- as.data.frame(counts(sim.pop.cond))
metadata0 <- as.data.frame(sim.pop.cond@colData@listData)

rm(gff)
rm(vcf)
rm(sim.pop.cond)
gc()

metadata0$Subtype <- ifelse(metadata0$Condition=="Condition1","Subtype1",
                            ifelse(metadata0$Condition=="Condition2","Subtype2","Subtype3"))
metadata0$Cell <- paste0('Cell',1:5000)
rownames(metadata0) <- metadata0$Cell
colnames(count0) <- metadata0$Cell
table(metadata0$Subtype)

library(caTools)
set.seed(123)
split <- sample.split(metadata0$Subtype, SplitRatio = 0.7)
train <- subset(metadata0, split==TRUE)
test <- subset(metadata0, split==FALSE)
count <- count0[,train$Cell]
metadata <- metadata0[train$Cell,]

## Prepare Bulk data
Subtype1 <- count[,metadata[metadata$Subtype=="Subtype1", ]$Cell]
Subtype2 <- count[,metadata[metadata$Subtype=="Subtype2", ]$Cell]
Subtype3 <- count[,metadata[metadata$Subtype=="Subtype3", ]$Cell]

################################################################################
library(SimBu)
annotation1 <- metadata[metadata$Subtype=="Subtype1",c(1,7)]
colnames(annotation1)[1] <- "ID"
annotation1$cell_type <- "All"
ds1 <- dataset(annotation = annotation1, 
               count_matrix = Matrix::Matrix(as.matrix(Subtype1), sparse = TRUE), 
               filter_genes = FALSE,name = "Cell")
simulation1 <- simulate_bulk(data = ds1,
                             scenario = "random", 
                             scaling_factor = "NONE", 
                             ncells=500, 
                             nsamples = 50)
bulk1 <- as.matrix(SummarizedExperiment::assays(simulation1$bulk)[["bulk_counts"]])
colnames(bulk1) <- paste0('sample',1:50)

annotation2 <- metadata[metadata$Subtype=="Subtype2",c(1,7)]
colnames(annotation2)[1] <- "ID"
annotation2$cell_type <- "All"
ds2 <- dataset(annotation = annotation2, 
               count_matrix = Matrix::Matrix(as.matrix(Subtype2), sparse = TRUE), 
               filter_genes = FALSE,name = "Cell")
simulation2 <- simulate_bulk(data = ds2,
                             scenario = "random", 
                             scaling_factor = "NONE", 
                             ncells=500, 
                             nsamples = 50)
bulk2 <- as.matrix(SummarizedExperiment::assays(simulation2$bulk)[["bulk_counts"]])
colnames(bulk2) <- paste0('sample',51:100)

annotation3 <- metadata[metadata$Subtype=="Subtype3",c(1,7)]
colnames(annotation3)[1] <- "ID"
annotation3$cell_type <- "All"
ds3 <- dataset(annotation = annotation3, 
               count_matrix = Matrix::Matrix(as.matrix(Subtype3), sparse = TRUE), 
               filter_genes = FALSE,name = "Cell")
simulation3 <- simulate_bulk(data = ds3,
                             scenario = "random", 
                             scaling_factor = "NONE", 
                             ncells=500, 
                             nsamples = 50)
bulk3 <- as.matrix(SummarizedExperiment::assays(simulation3$bulk)[["bulk_counts"]])
colnames(bulk3) <- paste0('sample',101:150)

bulk <- cbind(bulk1,bulk2,bulk3)
colData <- data.frame(Sample=c(colnames(bulk1),colnames(bulk2),colnames(bulk3)),
                      Subtype=rep(c("Subtype1","Subtype2","Subtype3"), each = 50))
rownames(colData) <- colData$Sample
colData$Subtype <- factor(colData$Subtype)

################################################################################
colData$Subtype1 <- ifelse(colData$Subtype=="Subtype1","Subtype1","Control")
colData$Subtype2 <- ifelse(colData$Subtype=="Subtype2","Subtype2","Control")
colData$Subtype3 <- ifelse(colData$Subtype=="Subtype3","Subtype3","Control")
colData$Subtype1 <- factor(colData$Subtype1)
colData$Subtype2 <- factor(colData$Subtype2)
colData$Subtype3 <- factor(colData$Subtype3)

## Prepare Bulk data
library("DESeq2")
dds <- DESeqDataSetFromMatrix(countData=bulk, colData=colData, design= ~ Subtype1)
dds <- DESeq(dds)
res <- as.data.frame(results(dds))

library(tidyverse)
UP1 <- res %>%
  filter(padj<0.05) %>%
  filter(log2FoldChange>=1) 
UP1 <- UP1 %>% mutate(class="Subtype1",genes=rownames(UP1))%>%
  select(genes,class)    # 51

dds <- DESeqDataSetFromMatrix(countData=bulk, colData=colData, design= ~ Subtype2)
dds <- DESeq(dds)
res <- as.data.frame(results(dds))
UP2 <- res %>%
  filter(padj<0.05) %>%
  filter(log2FoldChange>=1) 
UP2 <- UP2 %>% mutate(class="Subtype2",genes=rownames(UP2))%>%
  select(genes,class)    # 28

dds <- DESeqDataSetFromMatrix(countData=bulk, colData=colData, design= ~ Subtype3)
dds <- DESeq(dds)
res <- as.data.frame(results(dds))
UP3 <- res %>%
  filter(padj<0.05) %>%
  filter(log2FoldChange>=1) 
UP3 <- UP3 %>% mutate(class="Subtype3",genes=rownames(UP3))%>%
  select(genes,class)    # 36

markers <- rbind(UP1,UP2,UP3)
SC <- count0[,test$Cell]

library(PIPET)
tmp <- PIPET(SC, markers=markers,
             gene_col = "genes",
             class_col = "class")
tmp$Cell <- rownames(tmp)
tmp1 <- merge(tmp,metadata0[test$Cell,],by="Cell")
tmp1$Subtype <- as.factor(tmp1$Subtype)
# tmp1 <- tmp1 %>% filter(FDR<0.05)

################################################################################
library(caret)
confusionMatrix(tmp1$prediction, tmp1$Subtype)

tmp2 <- tmp1 %>% filter(FDR<0.05)
confusionMatrix(tmp2$prediction, tmp2$Subtype)
res_kappa <- 0.938
res_kappa_adj <- 1


subt1 <- as.numeric(tmp1$Subtype)
subt2 <- as.numeric(tmp1$prediction)
# calculate consistency
comb.subt <- data.frame(subt1 = paste0("Class", subt1),
                        subt2 = paste0("Class", subt2),
                        stringsAsFactors = F)
tab_classify <- as.data.frame.array(table(comb.subt$subt1,comb.subt$subt2))

# calculate kappa statistic
x <- table(comb.subt$subt1,comb.subt$subt2)
nr <- nrow(x); nc <- ncol(x); N <- sum(x)
Po <- sum(diag(x))/N; Pe <- sum(rowSums(x) * colSums(x)/N)/N
kappa <- (Po - Pe)/(1 - Pe)
seK0 <- sqrt(Pe/(N * (1 - Pe)))
p.v <- 1 - pnorm(kappa/seK0)
p.lab <- ifelse(p.v < 0.001, "P < 0.001", paste0("P = ", format(round(p.v,3), nsmall = 3)))

# generate consistency table
blue   <- "#cb392b"
lblue  <- "#f48c06"
dwhite <- "#f9c74f"
white  <- "#fffae5"

par(bty="n", mgp = c(2,0.5,0), mar = c(4.1,4.1,4.1,2.1),tcl=-.25, font.main=1)
par(xpd=NA)
plot(c(0,ncol(tab_classify)),c(0,nrow(tab_classify)),
     col = "white",main = "Three Classes",cex.main=1.5,line = 2.5,
     xlab = "",xaxt = "n",
     ylab = "",yaxt = "n")
title(paste0("Kappa = ", format(round(kappa,3), nsmall = 2),
             "\n", p.lab),
      adj = 0, line = 0)

# add y-axis
axis(2, at = 0.5:(nrow(tab_classify)-0.5), labels = FALSE)
text(y = 0.5:(nrow(tab_classify)-0.5),
     par("usr")[1],
     labels = rownames(tab_classify)[nrow(tab_classify):1],
     srt = 0, pos = 2, cex=1, xpd = TRUE)
mtext(paste0("Raw classes"),cex=1.2,  side = 2, line = 3)

# add x-axis
axis(1, at = 0.5:(ncol(tab_classify)-0.5), labels = FALSE)
text(x = 0.5:(ncol(tab_classify)-0.5),
     par("usr")[1] - 0.1,
     labels = colnames(tab_classify),
     srt = 0, pos = 1,cex=1, xpd = TRUE)
mtext(paste0("Predicted classes"),cex=1.2, side = 1, line = 2)

# generate colors
input_matrix <- as.matrix(tab_classify)
mat.max = max(input_matrix)
unq.value <- unique(sort(as.vector(input_matrix)))
rbPal <- colorRampPalette(c(white,dwhite,lblue,blue))
col.vec <- rbPal(max(unq.value) + 1)
col.mat <- matrix(NA,byrow = T,ncol = ncol(input_matrix),nrow = nrow(input_matrix))

# fill matrix
for (i in 1:nrow(col.mat)) {
  for (j in 1:ncol(col.mat)) {
    col.mat[i,j] <- col.vec[input_matrix[i,j] + 1]
  }
}

# generate heatmap
x_size <- ncol(input_matrix)
y_size <- nrow(input_matrix)

my_xleft = rep(c(0:(x_size-1)),each = x_size)
my_xright = my_xleft + 1
my_ybottom = rep(c((y_size-1):0),y_size)
my_ytop = my_ybottom + 1
rect(xleft = my_xleft,
     ybottom = my_ybottom,
     xright = my_xright,
     ytop = my_ytop,
     col=col.mat,
     border = F)

# fill count
text(my_xleft + 0.5,my_ybottom + 0.5,input_matrix, cex = 1)
# 
# # output to pdf
# outFig <- paste0("Simulation_Three_Class.pdf")
# invisible(dev.copy2pdf(file = file.path("E:/bulk and single-cell/", outFig), width = 6, height = 8))

########################## Three Class -- dropout ###############################
library("splatter")
library("scater")
library("VariantAnnotation")
library("ggplot2")

set.seed(123)
vcf <- mockVCF()
gff <- mockGFF(n.genes = 5000)
params.cond <- newSplatPopParams(eqtl.n = 0.3, 
                                 batchCells = 1000,
                                 similarity.scale = 5,
                                 condition.prob = c(0.4, 0.4, 0.2),
                                 eqtl.condition.specific = 0.5,
                                 cde.facLoc = 0.1, 
                                 cde.facScale = 0.2,
                                 dropout.mid =0.5,
                                 dropout.type = "experiment",
                                 dropout.shape=-1)   ## dropout 35%

sim.pop.cond <- splatPopSimulate(vcf = vcf, gff = gff, params = params.cond, 
                                 sparsify = FALSE)
# data.drop1 <- splatter:::splatSimDropout(sim, setParam(params, "dropout.mid", 0))

sim.pop.cond <- logNormCounts(sim.pop.cond)
sim.pop.cond <- runPCA(sim.pop.cond, ncomponents = 10)
plotPCA(sim.pop.cond, colour_by = "Condition", shape_by = "Sample")

################################################################################
## Prepare Single cell exp and metadata
count0 <- as.data.frame(counts(sim.pop.cond))
metadata0 <- as.data.frame(sim.pop.cond@colData@listData)

rm(gff)
rm(vcf)
rm(sim.pop.cond)
gc()

metadata0$Subtype <- ifelse(metadata0$Condition=="Condition1","Subtype1",
                            ifelse(metadata0$Condition=="Condition2","Subtype2","Subtype3"))
metadata0$Cell <- paste0('Cell',1:5000)
rownames(metadata0) <- metadata0$Cell
colnames(count0) <- metadata0$Cell
table(metadata0$Subtype)

library(caTools)
set.seed(123)
split <- sample.split(metadata0$Subtype, SplitRatio = 0.7)
train <- subset(metadata0, split==TRUE)
test <- subset(metadata0, split==FALSE)
count <- count0[,train$Cell]
metadata <- metadata0[train$Cell,]

sum(count==0)/(dim(count)[1]*dim(count)[2])
sum(count0==0)/(dim(count0)[1]*dim(count0)[2])

## Prepare Bulk data
Subtype1 <- count[,metadata[metadata$Subtype=="Subtype1", ]$Cell]
Subtype2 <- count[,metadata[metadata$Subtype=="Subtype2", ]$Cell]
Subtype3 <- count[,metadata[metadata$Subtype=="Subtype3", ]$Cell]

################################################################################
library(SimBu)
annotation1 <- metadata[metadata$Subtype=="Subtype1",c(1,8)]
colnames(annotation1)[1] <- "ID"
annotation1$cell_type <- "All"
ds1 <- dataset(annotation = annotation1, 
               count_matrix = Matrix::Matrix(as.matrix(Subtype1), sparse = TRUE), 
               filter_genes = FALSE,name = "Cell")
simulation1 <- simulate_bulk(data = ds1,
                             scenario = "random", 
                             scaling_factor = "NONE", 
                             ncells=500, 
                             nsamples = 50)
bulk1 <- as.matrix(SummarizedExperiment::assays(simulation1$bulk)[["bulk_counts"]])
colnames(bulk1) <- paste0('sample',1:50)

annotation2 <- metadata[metadata$Subtype=="Subtype2",c(1,8)]
colnames(annotation2)[1] <- "ID"
annotation2$cell_type <- "All"
ds2 <- dataset(annotation = annotation2, 
               count_matrix = Matrix::Matrix(as.matrix(Subtype2), sparse = TRUE), 
               filter_genes = FALSE,name = "Cell")
simulation2 <- simulate_bulk(data = ds2,
                             scenario = "random", 
                             scaling_factor = "NONE", 
                             ncells=500, 
                             nsamples = 50)
bulk2 <- as.matrix(SummarizedExperiment::assays(simulation2$bulk)[["bulk_counts"]])
colnames(bulk2) <- paste0('sample',51:100)

annotation3 <- metadata[metadata$Subtype=="Subtype3",c(1,8)]
colnames(annotation3)[1] <- "ID"
annotation3$cell_type <- "All"
ds3 <- dataset(annotation = annotation3, 
               count_matrix = Matrix::Matrix(as.matrix(Subtype3), sparse = TRUE), 
               filter_genes = FALSE,name = "Cell")
simulation3 <- simulate_bulk(data = ds3,
                             scenario = "random", 
                             scaling_factor = "NONE", 
                             ncells=500, 
                             nsamples = 50)
bulk3 <- as.matrix(SummarizedExperiment::assays(simulation3$bulk)[["bulk_counts"]])
colnames(bulk3) <- paste0('sample',101:150)

bulk <- cbind(bulk1,bulk2,bulk3)
colData <- data.frame(Sample=c(colnames(bulk1),colnames(bulk2),colnames(bulk3)),
                      Subtype=rep(c("Subtype1","Subtype2","Subtype3"), each = 50))
rownames(colData) <- colData$Sample
colData$Subtype <- factor(colData$Subtype)

################################################################################
colData$Subtype1 <- ifelse(colData$Subtype=="Subtype1","Subtype1","Control")
colData$Subtype2 <- ifelse(colData$Subtype=="Subtype2","Subtype2","Control")
colData$Subtype3 <- ifelse(colData$Subtype=="Subtype3","Subtype3","Control")
colData$Subtype1 <- factor(colData$Subtype1)
colData$Subtype2 <- factor(colData$Subtype2)
colData$Subtype3 <- factor(colData$Subtype3)

## Prepare Bulk data
library("DESeq2")
dds <- DESeqDataSetFromMatrix(countData=bulk, colData=colData, design= ~ Subtype1)
dds <- DESeq(dds)
res <- as.data.frame(results(dds))

library(tidyverse)
UP1 <- res %>%
  filter(padj<0.05) %>%
  filter(log2FoldChange>=1) 
UP1 <- UP1 %>% mutate(class="Subtype1",genes=rownames(UP1))%>%
  select(genes,class)  # 85

dds <- DESeqDataSetFromMatrix(countData=bulk, colData=colData, design= ~ Subtype2)
dds <- DESeq(dds)
res <- as.data.frame(results(dds))
UP2 <- res %>%
  filter(padj<0.05) %>%
  filter(log2FoldChange>=1) 
UP2 <- UP2 %>% mutate(class="Subtype2",genes=rownames(UP2))%>%
  select(genes,class)   # 65

dds <- DESeqDataSetFromMatrix(countData=bulk, colData=colData, design= ~ Subtype3)
dds <- DESeq(dds)
res <- as.data.frame(results(dds))
UP3 <- res %>%
  filter(padj<0.05) %>%
  filter(log2FoldChange>=1) 
UP3 <- UP3 %>% mutate(class="Subtype3",genes=rownames(UP3))%>%
  select(genes,class)  # 71

markers <- rbind(UP1,UP2,UP3)
SC <- count0[,test$Cell]

library(PIPET)
tmp <- PIPET(SC, markers=markers,
             gene_col = "genes",
             class_col = "class")
tmp$Cell <- rownames(tmp)
tmp1 <- merge(tmp,metadata0[test$Cell,],by="Cell")
tmp1$Subtype <- as.factor(tmp1$Subtype)
# tmp1 <- tmp1 %>% filter(FDR<0.05)

################################################################################
library(caret)
confusionMatrix(tmp1$prediction, tmp1$Subtype)

tmp2 <- tmp1 %>% filter(FDR<0.05)
confusionMatrix(tmp2$prediction, tmp2$Subtype)
res_kappa[2] <- 0.904
res_kappa_adj[2] <- 0.9965


subt1 <- as.numeric(tmp1$Subtype)
subt2 <- as.numeric(tmp1$prediction)
# calculate consistency
comb.subt <- data.frame(subt1 = paste0("Class", subt1),
                        subt2 = paste0("Class", subt2),
                        stringsAsFactors = F)
tab_classify <- as.data.frame.array(table(comb.subt$subt1,comb.subt$subt2))

# calculate kappa statistic
x <- table(comb.subt$subt1,comb.subt$subt2)
nr <- nrow(x); nc <- ncol(x); N <- sum(x)
Po <- sum(diag(x))/N; Pe <- sum(rowSums(x) * colSums(x)/N)/N
kappa <- (Po - Pe)/(1 - Pe)
seK0 <- sqrt(Pe/(N * (1 - Pe)))
p.v <- 1 - pnorm(kappa/seK0)
p.lab <- ifelse(p.v < 0.001, "P < 0.001", paste0("P = ", format(round(p.v,3), nsmall = 3)))

# generate consistency table
blue   <- "#cb392b"
lblue  <- "#f48c06"
dwhite <- "#f9c74f"
white  <- "#fffae5"

par(bty="n", mgp = c(2,0.5,0), mar = c(4.1,4.1,4.1,2.1),tcl=-.25, font.main=1)
par(xpd=NA)
plot(c(0,ncol(tab_classify)),c(0,nrow(tab_classify)),
     col = "white",main = "Three Classes (35% dropout)",cex.main=1.5,line = 2.5,
     xlab = "",xaxt = "n",
     ylab = "",yaxt = "n")
title(paste0("Kappa = ", format(round(kappa,3), nsmall = 2),
             "\n", p.lab),
      adj = 0, line = 0)

# add y-axis
axis(2, at = 0.5:(nrow(tab_classify)-0.5), labels = FALSE)
text(y = 0.5:(nrow(tab_classify)-0.5),
     par("usr")[1],
     labels = rownames(tab_classify)[nrow(tab_classify):1],
     srt = 0, pos = 2, cex=1, xpd = TRUE)
mtext(paste0("Raw classes"),cex=1.2,  side = 2, line = 3)

# add x-axis
axis(1, at = 0.5:(ncol(tab_classify)-0.5), labels = FALSE)
text(x = 0.5:(ncol(tab_classify)-0.5),
     par("usr")[1] - 0.1,
     labels = colnames(tab_classify),
     srt = 0, pos = 1,cex=1, xpd = TRUE)
mtext(paste0("Predicted classes"),cex=1.2, side = 1, line = 2)

# generate colors
input_matrix <- as.matrix(tab_classify)
mat.max = max(input_matrix)
unq.value <- unique(sort(as.vector(input_matrix)))
rbPal <- colorRampPalette(c(white,dwhite,lblue,blue))
col.vec <- rbPal(max(unq.value) + 1)
col.mat <- matrix(NA,byrow = T,ncol = ncol(input_matrix),nrow = nrow(input_matrix))

# fill matrix
for (i in 1:nrow(col.mat)) {
  for (j in 1:ncol(col.mat)) {
    col.mat[i,j] <- col.vec[input_matrix[i,j] + 1]
  }
}

# generate heatmap
x_size <- ncol(input_matrix)
y_size <- nrow(input_matrix)

my_xleft = rep(c(0:(x_size-1)),each = x_size)
my_xright = my_xleft + 1
my_ybottom = rep(c((y_size-1):0),y_size)
my_ytop = my_ybottom + 1
rect(xleft = my_xleft,
     ybottom = my_ybottom,
     xright = my_xright,
     ytop = my_ytop,
     col=col.mat,
     border = F)

# fill count
text(my_xleft + 0.5,my_ybottom + 0.5,input_matrix, cex = 1)
# 
# # output to pdf
# outFig <- paste0("Simulation_Three_Class.pdf")
# invisible(dev.copy2pdf(file = file.path("E:/bulk and single-cell/", outFig), width = 6, height = 8))

########################## Three Class -- dropout ###############################
library("splatter")
library("scater")
library("VariantAnnotation")
library("ggplot2")

set.seed(123)
vcf <- mockVCF()
gff <- mockGFF(n.genes = 5000)
params.cond <- newSplatPopParams(eqtl.n = 0.3, 
                                 batchCells = 1000,
                                 similarity.scale = 5,
                                 condition.prob = c(0.4, 0.4, 0.2),
                                 eqtl.condition.specific = 0.5,
                                 cde.facLoc = 0.1, 
                                 cde.facScale = 0.2,
                                 dropout.mid =1.4,
                                 dropout.type = "experiment",
                                 dropout.shape=-1)   ## dropout 50%

sim.pop.cond <- splatPopSimulate(vcf = vcf, gff = gff, params = params.cond, 
                                 sparsify = FALSE)
# data.drop1 <- splatter:::splatSimDropout(sim, setParam(params, "dropout.mid", 0))

sim.pop.cond <- logNormCounts(sim.pop.cond)
sim.pop.cond <- runPCA(sim.pop.cond, ncomponents = 10)
plotPCA(sim.pop.cond, colour_by = "Condition", shape_by = "Sample")

################################################################################
## Prepare Single cell exp and metadata
count0 <- as.data.frame(counts(sim.pop.cond))
metadata0 <- as.data.frame(sim.pop.cond@colData@listData)

rm(gff)
rm(vcf)
rm(sim.pop.cond)
gc()

metadata0$Subtype <- ifelse(metadata0$Condition=="Condition1","Subtype1",
                            ifelse(metadata0$Condition=="Condition2","Subtype2","Subtype3"))
metadata0$Cell <- paste0('Cell',1:5000)
rownames(metadata0) <- metadata0$Cell
colnames(count0) <- metadata0$Cell
table(metadata0$Subtype)

library(caTools)
set.seed(123)
split <- sample.split(metadata0$Subtype, SplitRatio = 0.7)
train <- subset(metadata0, split==TRUE)
test <- subset(metadata0, split==FALSE)
count <- count0[,train$Cell]
metadata <- metadata0[train$Cell,]

sum(count==0)/(dim(count)[1]*dim(count)[2])
sum(count0==0)/(dim(count0)[1]*dim(count0)[2])

## Prepare Bulk data
Subtype1 <- count[,metadata[metadata$Subtype=="Subtype1", ]$Cell]
Subtype2 <- count[,metadata[metadata$Subtype=="Subtype2", ]$Cell]
Subtype3 <- count[,metadata[metadata$Subtype=="Subtype3", ]$Cell]

################################################################################
library(SimBu)
annotation1 <- metadata[metadata$Subtype=="Subtype1",c(1,8)]
colnames(annotation1)[1] <- "ID"
annotation1$cell_type <- "All"
ds1 <- dataset(annotation = annotation1, 
               count_matrix = Matrix::Matrix(as.matrix(Subtype1), sparse = TRUE), 
               filter_genes = FALSE,name = "Cell")
simulation1 <- simulate_bulk(data = ds1,
                             scenario = "random", 
                             scaling_factor = "NONE", 
                             ncells=500, 
                             nsamples = 50)
bulk1 <- as.matrix(SummarizedExperiment::assays(simulation1$bulk)[["bulk_counts"]])
colnames(bulk1) <- paste0('sample',1:50)

annotation2 <- metadata[metadata$Subtype=="Subtype2",c(1,8)]
colnames(annotation2)[1] <- "ID"
annotation2$cell_type <- "All"
ds2 <- dataset(annotation = annotation2, 
               count_matrix = Matrix::Matrix(as.matrix(Subtype2), sparse = TRUE), 
               filter_genes = FALSE,name = "Cell")
simulation2 <- simulate_bulk(data = ds2,
                             scenario = "random", 
                             scaling_factor = "NONE", 
                             ncells=500, 
                             nsamples = 50)
bulk2 <- as.matrix(SummarizedExperiment::assays(simulation2$bulk)[["bulk_counts"]])
colnames(bulk2) <- paste0('sample',51:100)

annotation3 <- metadata[metadata$Subtype=="Subtype3",c(1,8)]
colnames(annotation3)[1] <- "ID"
annotation3$cell_type <- "All"
ds3 <- dataset(annotation = annotation3, 
               count_matrix = Matrix::Matrix(as.matrix(Subtype3), sparse = TRUE), 
               filter_genes = FALSE,name = "Cell")
simulation3 <- simulate_bulk(data = ds3,
                             scenario = "random", 
                             scaling_factor = "NONE", 
                             ncells=500, 
                             nsamples = 50)
bulk3 <- as.matrix(SummarizedExperiment::assays(simulation3$bulk)[["bulk_counts"]])
colnames(bulk3) <- paste0('sample',101:150)

bulk <- cbind(bulk1,bulk2,bulk3)
colData <- data.frame(Sample=c(colnames(bulk1),colnames(bulk2),colnames(bulk3)),
                      Subtype=rep(c("Subtype1","Subtype2","Subtype3"), each = 50))
rownames(colData) <- colData$Sample
colData$Subtype <- factor(colData$Subtype)

################################################################################
colData$Subtype1 <- ifelse(colData$Subtype=="Subtype1","Subtype1","Control")
colData$Subtype2 <- ifelse(colData$Subtype=="Subtype2","Subtype2","Control")
colData$Subtype3 <- ifelse(colData$Subtype=="Subtype3","Subtype3","Control")
colData$Subtype1 <- factor(colData$Subtype1)
colData$Subtype2 <- factor(colData$Subtype2)
colData$Subtype3 <- factor(colData$Subtype3)

## Prepare Bulk data
library("DESeq2")
dds <- DESeqDataSetFromMatrix(countData=bulk, colData=colData, design= ~ Subtype1)
dds <- DESeq(dds)
res <- as.data.frame(results(dds))

library(tidyverse)
UP1 <- res %>%
  filter(padj<0.05) %>%
  filter(log2FoldChange>=1) 
UP1 <- UP1 %>% mutate(class="Subtype1",genes=rownames(UP1))%>%
  select(genes,class)  # 126

dds <- DESeqDataSetFromMatrix(countData=bulk, colData=colData, design= ~ Subtype2)
dds <- DESeq(dds)
res <- as.data.frame(results(dds))
UP2 <- res %>%
  filter(padj<0.05) %>%
  filter(log2FoldChange>=1) 
UP2 <- UP2 %>% mutate(class="Subtype2",genes=rownames(UP2))%>%
  select(genes,class)   # 101

dds <- DESeqDataSetFromMatrix(countData=bulk, colData=colData, design= ~ Subtype3)
dds <- DESeq(dds)
res <- as.data.frame(results(dds))
UP3 <- res %>%
  filter(padj<0.05) %>%
  filter(log2FoldChange>=1) 
UP3 <- UP3 %>% mutate(class="Subtype3",genes=rownames(UP3))%>%
  select(genes,class)  # 91

markers <- rbind(UP1,UP2,UP3)
SC <- count0[,test$Cell]

library(PIPET)
tmp <- PIPET(SC, markers=markers,
             gene_col = "genes",
             class_col = "class")
tmp$Cell <- rownames(tmp)
tmp1 <- merge(tmp,metadata0[test$Cell,],by="Cell")
tmp1$Subtype <- as.factor(tmp1$Subtype)
# tmp1 <- tmp1 %>% filter(FDR<0.05)

################################################################################
library(caret)
confusionMatrix(tmp1$prediction, tmp1$Subtype)

tmp2 <- tmp1 %>% filter(FDR<0.05)
confusionMatrix(tmp2$prediction, tmp2$Subtype)
res_kappa[3] <- 0.871
res_kappa_adj[3] <- 0.9952


subt1 <- as.numeric(tmp1$Subtype)
subt2 <- as.numeric(tmp1$prediction)
# calculate consistency
comb.subt <- data.frame(subt1 = paste0("Class", subt1),
                        subt2 = paste0("Class", subt2),
                        stringsAsFactors = F)
tab_classify <- as.data.frame.array(table(comb.subt$subt1,comb.subt$subt2))

# calculate kappa statistic
x <- table(comb.subt$subt1,comb.subt$subt2)
nr <- nrow(x); nc <- ncol(x); N <- sum(x)
Po <- sum(diag(x))/N; Pe <- sum(rowSums(x) * colSums(x)/N)/N
kappa <- (Po - Pe)/(1 - Pe)
seK0 <- sqrt(Pe/(N * (1 - Pe)))
p.v <- 1 - pnorm(kappa/seK0)
p.lab <- ifelse(p.v < 0.001, "P < 0.001", paste0("P = ", format(round(p.v,3), nsmall = 3)))

# generate consistency table
blue   <- "#cb392b"
lblue  <- "#f48c06"
dwhite <- "#f9c74f"
white  <- "#fffae5"

par(bty="n", mgp = c(2,0.5,0), mar = c(4.1,4.1,4.1,2.1),tcl=-.25, font.main=1)
par(xpd=NA)
plot(c(0,ncol(tab_classify)),c(0,nrow(tab_classify)),
     col = "white",main = "Three Classes (50% dropout)",cex.main=1.5,line = 2.5,
     xlab = "",xaxt = "n",
     ylab = "",yaxt = "n")
title(paste0("Kappa = ", format(round(kappa,3), nsmall = 2),
             "\n", p.lab),
      adj = 0, line = 0)

# add y-axis
axis(2, at = 0.5:(nrow(tab_classify)-0.5), labels = FALSE)
text(y = 0.5:(nrow(tab_classify)-0.5),
     par("usr")[1],
     labels = rownames(tab_classify)[nrow(tab_classify):1],
     srt = 0, pos = 2, cex=1, xpd = TRUE)
mtext(paste0("Raw classes"),cex=1.2,  side = 2, line = 3)

# add x-axis
axis(1, at = 0.5:(ncol(tab_classify)-0.5), labels = FALSE)
text(x = 0.5:(ncol(tab_classify)-0.5),
     par("usr")[1] - 0.1,
     labels = colnames(tab_classify),
     srt = 0, pos = 1,cex=1, xpd = TRUE)
mtext(paste0("Predicted classes"),cex=1.2, side = 1, line = 2)

# generate colors
input_matrix <- as.matrix(tab_classify)
mat.max = max(input_matrix)
unq.value <- unique(sort(as.vector(input_matrix)))
rbPal <- colorRampPalette(c(white,dwhite,lblue,blue))
col.vec <- rbPal(max(unq.value) + 1)
col.mat <- matrix(NA,byrow = T,ncol = ncol(input_matrix),nrow = nrow(input_matrix))

# fill matrix
for (i in 1:nrow(col.mat)) {
  for (j in 1:ncol(col.mat)) {
    col.mat[i,j] <- col.vec[input_matrix[i,j] + 1]
  }
}

# generate heatmap
x_size <- ncol(input_matrix)
y_size <- nrow(input_matrix)

my_xleft = rep(c(0:(x_size-1)),each = x_size)
my_xright = my_xleft + 1
my_ybottom = rep(c((y_size-1):0),y_size)
my_ytop = my_ybottom + 1
rect(xleft = my_xleft,
     ybottom = my_ybottom,
     xright = my_xright,
     ytop = my_ytop,
     col=col.mat,
     border = F)

# fill count
text(my_xleft + 0.5,my_ybottom + 0.5,input_matrix, cex = 1)

########################## Three Class -- dropout ###############################
library("splatter")
library("scater")
library("VariantAnnotation")
library("ggplot2")

set.seed(123)
vcf <- mockVCF()
gff <- mockGFF(n.genes = 5000)
params.cond <- newSplatPopParams(eqtl.n = 0.3, 
                                 batchCells = 1000,
                                 similarity.scale = 5,
                                 condition.prob = c(0.4, 0.4, 0.2),
                                 eqtl.condition.specific = 0.5,
                                 cde.facLoc = 0.1, 
                                 cde.facScale = 0.2,
                                 dropout.mid =2.4,
                                 dropout.type = "experiment",
                                 dropout.shape=-1)   ## dropout 65%

sim.pop.cond <- splatPopSimulate(vcf = vcf, gff = gff, params = params.cond, 
                                 sparsify = FALSE)
# data.drop1 <- splatter:::splatSimDropout(sim, setParam(params, "dropout.mid", 0))

sim.pop.cond <- logNormCounts(sim.pop.cond)
sim.pop.cond <- runPCA(sim.pop.cond, ncomponents = 10)
plotPCA(sim.pop.cond, colour_by = "Condition", shape_by = "Sample")

################################################################################
## Prepare Single cell exp and metadata
count0 <- as.data.frame(counts(sim.pop.cond))
metadata0 <- as.data.frame(sim.pop.cond@colData@listData)

rm(gff)
rm(vcf)
rm(sim.pop.cond)
gc()

metadata0$Subtype <- ifelse(metadata0$Condition=="Condition1","Subtype1",
                            ifelse(metadata0$Condition=="Condition2","Subtype2","Subtype3"))
metadata0$Cell <- paste0('Cell',1:5000)
rownames(metadata0) <- metadata0$Cell
colnames(count0) <- metadata0$Cell
table(metadata0$Subtype)

library(caTools)
set.seed(123)
split <- sample.split(metadata0$Subtype, SplitRatio = 0.7)
train <- subset(metadata0, split==TRUE)
test <- subset(metadata0, split==FALSE)
count <- count0[,train$Cell]
metadata <- metadata0[train$Cell,]

sum(count==0)/(dim(count)[1]*dim(count)[2])
sum(count0==0)/(dim(count0)[1]*dim(count0)[2])

## Prepare Bulk data
Subtype1 <- count[,metadata[metadata$Subtype=="Subtype1", ]$Cell]
Subtype2 <- count[,metadata[metadata$Subtype=="Subtype2", ]$Cell]
Subtype3 <- count[,metadata[metadata$Subtype=="Subtype3", ]$Cell]

################################################################################
library(SimBu)
annotation1 <- metadata[metadata$Subtype=="Subtype1",c(1,8)]
colnames(annotation1)[1] <- "ID"
annotation1$cell_type <- "All"
ds1 <- dataset(annotation = annotation1, 
               count_matrix = Matrix::Matrix(as.matrix(Subtype1), sparse = TRUE), 
               filter_genes = FALSE,name = "Cell")
simulation1 <- simulate_bulk(data = ds1,
                             scenario = "random", 
                             scaling_factor = "NONE", 
                             ncells=500, 
                             nsamples = 50)
bulk1 <- as.matrix(SummarizedExperiment::assays(simulation1$bulk)[["bulk_counts"]])
colnames(bulk1) <- paste0('sample',1:50)

annotation2 <- metadata[metadata$Subtype=="Subtype2",c(1,8)]
colnames(annotation2)[1] <- "ID"
annotation2$cell_type <- "All"
ds2 <- dataset(annotation = annotation2, 
               count_matrix = Matrix::Matrix(as.matrix(Subtype2), sparse = TRUE), 
               filter_genes = FALSE,name = "Cell")
simulation2 <- simulate_bulk(data = ds2,
                             scenario = "random", 
                             scaling_factor = "NONE", 
                             ncells=500, 
                             nsamples = 50)
bulk2 <- as.matrix(SummarizedExperiment::assays(simulation2$bulk)[["bulk_counts"]])
colnames(bulk2) <- paste0('sample',51:100)

annotation3 <- metadata[metadata$Subtype=="Subtype3",c(1,8)]
colnames(annotation3)[1] <- "ID"
annotation3$cell_type <- "All"
ds3 <- dataset(annotation = annotation3, 
               count_matrix = Matrix::Matrix(as.matrix(Subtype3), sparse = TRUE), 
               filter_genes = FALSE,name = "Cell")
simulation3 <- simulate_bulk(data = ds3,
                             scenario = "random", 
                             scaling_factor = "NONE", 
                             ncells=500, 
                             nsamples = 50)
bulk3 <- as.matrix(SummarizedExperiment::assays(simulation3$bulk)[["bulk_counts"]])
colnames(bulk3) <- paste0('sample',101:150)

bulk <- cbind(bulk1,bulk2,bulk3)
colData <- data.frame(Sample=c(colnames(bulk1),colnames(bulk2),colnames(bulk3)),
                      Subtype=rep(c("Subtype1","Subtype2","Subtype3"), each = 50))
rownames(colData) <- colData$Sample
colData$Subtype <- factor(colData$Subtype)

################################################################################
colData$Subtype1 <- ifelse(colData$Subtype=="Subtype1","Subtype1","Control")
colData$Subtype2 <- ifelse(colData$Subtype=="Subtype2","Subtype2","Control")
colData$Subtype3 <- ifelse(colData$Subtype=="Subtype3","Subtype3","Control")
colData$Subtype1 <- factor(colData$Subtype1)
colData$Subtype2 <- factor(colData$Subtype2)
colData$Subtype3 <- factor(colData$Subtype3)

## Prepare Bulk data
library("DESeq2")
dds <- DESeqDataSetFromMatrix(countData=bulk, colData=colData, design= ~ Subtype1)
dds <- DESeq(dds)
res <- as.data.frame(results(dds))

library(tidyverse)
UP1 <- res %>%
  filter(padj<0.05) %>%
  filter(log2FoldChange>=1) 
UP1 <- UP1 %>% mutate(class="Subtype1",genes=rownames(UP1))%>%
  select(genes,class)  # 139

dds <- DESeqDataSetFromMatrix(countData=bulk, colData=colData, design= ~ Subtype2)
dds <- DESeq(dds)
res <- as.data.frame(results(dds))
UP2 <- res %>%
  filter(padj<0.05) %>%
  filter(log2FoldChange>=1) 
UP2 <- UP2 %>% mutate(class="Subtype2",genes=rownames(UP2))%>%
  select(genes,class)   # 166

dds <- DESeqDataSetFromMatrix(countData=bulk, colData=colData, design= ~ Subtype3)
dds <- DESeq(dds)
res <- as.data.frame(results(dds))
UP3 <- res %>%
  filter(padj<0.05) %>%
  filter(log2FoldChange>=1) 
UP3 <- UP3 %>% mutate(class="Subtype3",genes=rownames(UP3))%>%
  select(genes,class)  # 151

markers <- rbind(UP1,UP2,UP3)
SC <- count0[,test$Cell]

library(PIPET)
tmp <- PIPET(SC, markers=markers,
             gene_col = "genes",
             class_col = "class")
tmp$Cell <- rownames(tmp)
tmp1 <- merge(tmp,metadata0[test$Cell,],by="Cell")
tmp1$Subtype <- as.factor(tmp1$Subtype)
# tmp1 <- tmp1 %>% filter(FDR<0.05)

################################################################################
library(caret)
confusionMatrix(tmp1$prediction, tmp1$Subtype)

tmp2 <- tmp1 %>% filter(FDR<0.05)
confusionMatrix(tmp2$prediction, tmp2$Subtype)
res_kappa[4] <- 0.829
res_kappa_adj[4] <- 0.9811


subt1 <- as.numeric(tmp1$Subtype)
subt2 <- as.numeric(tmp1$prediction)
# calculate consistency
comb.subt <- data.frame(subt1 = paste0("Class", subt1),
                        subt2 = paste0("Class", subt2),
                        stringsAsFactors = F)
tab_classify <- as.data.frame.array(table(comb.subt$subt1,comb.subt$subt2))

# calculate kappa statistic
x <- table(comb.subt$subt1,comb.subt$subt2)
nr <- nrow(x); nc <- ncol(x); N <- sum(x)
Po <- sum(diag(x))/N; Pe <- sum(rowSums(x) * colSums(x)/N)/N
kappa <- (Po - Pe)/(1 - Pe)
seK0 <- sqrt(Pe/(N * (1 - Pe)))
p.v <- 1 - pnorm(kappa/seK0)
p.lab <- ifelse(p.v < 0.001, "P < 0.001", paste0("P = ", format(round(p.v,3), nsmall = 3)))

# generate consistency table
blue   <- "#cb392b"
lblue  <- "#f48c06"
dwhite <- "#f9c74f"
white  <- "#fffae5"

par(bty="n", mgp = c(2,0.5,0), mar = c(4.1,4.1,4.1,2.1),tcl=-.25, font.main=1)
par(xpd=NA)
plot(c(0,ncol(tab_classify)),c(0,nrow(tab_classify)),
     col = "white",main = "Three Classes (65% dropout)",cex.main=1.5,line = 2.5,
     xlab = "",xaxt = "n",
     ylab = "",yaxt = "n")
title(paste0("Kappa = ", format(round(kappa,3), nsmall = 2),
             "\n", p.lab),
      adj = 0, line = 0)

# add y-axis
axis(2, at = 0.5:(nrow(tab_classify)-0.5), labels = FALSE)
text(y = 0.5:(nrow(tab_classify)-0.5),
     par("usr")[1],
     labels = rownames(tab_classify)[nrow(tab_classify):1],
     srt = 0, pos = 2, cex=1, xpd = TRUE)
mtext(paste0("Raw classes"),cex=1.2,  side = 2, line = 3)

# add x-axis
axis(1, at = 0.5:(ncol(tab_classify)-0.5), labels = FALSE)
text(x = 0.5:(ncol(tab_classify)-0.5),
     par("usr")[1] - 0.1,
     labels = colnames(tab_classify),
     srt = 0, pos = 1,cex=1, xpd = TRUE)
mtext(paste0("Predicted classes"),cex=1.2, side = 1, line = 2)

# generate colors
input_matrix <- as.matrix(tab_classify)
mat.max = max(input_matrix)
unq.value <- unique(sort(as.vector(input_matrix)))
rbPal <- colorRampPalette(c(white,dwhite,lblue,blue))
col.vec <- rbPal(max(unq.value) + 1)
col.mat <- matrix(NA,byrow = T,ncol = ncol(input_matrix),nrow = nrow(input_matrix))

# fill matrix
for (i in 1:nrow(col.mat)) {
  for (j in 1:ncol(col.mat)) {
    col.mat[i,j] <- col.vec[input_matrix[i,j] + 1]
  }
}

# generate heatmap
x_size <- ncol(input_matrix)
y_size <- nrow(input_matrix)

my_xleft = rep(c(0:(x_size-1)),each = x_size)
my_xright = my_xleft + 1
my_ybottom = rep(c((y_size-1):0),y_size)
my_ytop = my_ybottom + 1
rect(xleft = my_xleft,
     ybottom = my_ybottom,
     xright = my_xright,
     ytop = my_ytop,
     col=col.mat,
     border = F)

# fill count
text(my_xleft + 0.5,my_ybottom + 0.5,input_matrix, cex = 1)


data_kappa <- data.frame(kappa=res_kappa,adj_kappa=res_kappa_adj,
                         type=c("3C_no_drops","3C_35_drops","3C_50_drops","3C_65_drops"))
write.csv(data_kappa,file="3C_kappa_results.csv",quote=F,row.names = F)