source("/usr/local/hdd2/data/alejandro/functions.R")
source("/usr/local/hdd2/data/alejandro/Christoph_adrb1/seurat_ADRB1/functions_new.R")
#
# this searches for all input matrices
#
files <- Sys.glob("/usr/local/hdd2/data/alejandro/Christoph_adrb1/ADRB1_blood/outs/raw_feature_bc_matrix/*.mtx.gz")
#files = grep(x=files, pattern="old", value=T, invert=T)
#files = grep(x=files, pattern="sepsis", value=T)

allfiles.raw = list()
allABs.raw = list()
for (file in files)
{
  samplename = str_split(dirname(file), "/")[[1]][8] #estos numeros indican el nombre de la carpeta
  foldername = dirname(file)
  
  print(paste(samplename, foldername))
  
  h5file = Read10X(foldername,unique.features = TRUE)

  if (is.null(names(h5file)))
  {
      print(paste("WITHOUT AB", samplename))
    allfiles.raw[[samplename]] = h5file
  } else {
      print(paste("WITH AB", samplename))
    allfiles.raw[[samplename]] = h5file$`Gene Expression`
    allABs.raw[[samplename]] = h5file$`Antibody Capture`
  }

  print(paste(samplename, nrow(allfiles.raw[[samplename]]), "x", ncol(allfiles.raw[[samplename]]), "genes x cells"))
}

joint.bcs <- intersect(colnames(allfiles.raw$ADRB1_blood), colnames(allABs.raw$ADRB1_blood))
obj.umis <- allfiles.raw$ADRB1_blood[, joint.bcs]
obj.htos <- as.matrix(allABs.raw$ADRB1_blood[, joint.bcs])


# Confirm that the HTO have the correct names
rownames(obj.htos)


#
# here we create a list of seurat object. each entry corresponds to an input matrix from above
#

objlist = list()
for (x in names(allfiles.raw))
{

    matrix = allfiles.raw[[x]]
    
    # this creates a Seurat object from the count matrix. it sets the object's project to x and prepends the sample name to all cells
    # the patternlist.mouse contains patterns for mt and RP-genes
    filteredObj = makeSeuratObj(matrix, x, patternList.mouse)
    
    # this creates log-normalized count matrices in RNA assay
    filteredObj <- NormalizeData(filteredObj, verbose = FALSE)
    # this calculates the most (2000) variable features per data set. variable features are features which show a high variance between all cells of a sample
    filteredObj <- FindVariableFeatures(filteredObj, verbose = FALSE)
    
    objlist[[x]] = filteredObj

    print(x)
    print(filteredObj)
    
    
}

names(objlist)

objlist.raw = objlist

objlist <- lapply(X = objlist.raw, FUN = function(obj) {
  # mt content: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6072887/
  print(paste("Seurat obj project", obj@project.name))
  print(obj)
  obj <- subset(obj, subset = nFeature_RNA > 100 & nFeature_RNA < 6000 & nCount_RNA > 100)
  obj <- subset(obj, subset = percent.mt < 15)
  print(obj)
  
  return(obj)
})



for (name in names(objlist))
{
  p=VlnPlot(objlist[[name]], features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0)
  save_plot(p, paste("QCplots", paste(name, "filtered_violins_qc", sep="_"), sep="/"), fig.width=10, fig.height=6)

  p=VlnPlot(objlist[[name]], features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0, combine=F)
  p[[1]] = p[[1]] + scale_y_continuous(limits = c(0, 1000), breaks = seq(0,1000,100))
  p[[2]] = p[[2]] + scale_y_continuous(limits = c(0, 1000), breaks = seq(0,1000,100))
  p[[3]] = p[[3]] + scale_y_continuous(limits = c(0, 100), breaks = seq(0,100,5))
  p = combine_plot_grid_list(plotlist=p, ncol=3)
  save_plot(p, paste("QCplots", paste(name, "filtered_violins_detail_qc", sep="_"), sep="/"), fig.width=18, fig.height=6)
  
  
  plot1 <- FeatureScatter(objlist[[name]], feature1 = "nCount_RNA", feature2 = "percent.mt")
  plot2 <- FeatureScatter(objlist[[name]], feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
  save_plot(plot1 + plot2, paste("QCplots", paste(name, "filtered_scatter_ncount_mt", sep="_"), sep="/"), fig.width=10, fig.height=6)
  
  plot1 <- FeatureScatter(objlist[[name]], feature1 = "nCount_RNA", feature2 = "percent.rp")
  plot2 <- FeatureScatter(objlist[[name]], feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
  save_plot(plot1 + plot2, paste("QCplots", paste(name, "filtered_scatter_ncount_rp", sep="_"), sep="/"), fig.width=10, fig.height=6)
}

ADRB1_blood= objlist$ADRB1_blood

ADRB1_blood <- FindVariableFeatures(ADRB1_blood, selection.method = "vst", nfeatures = 3000)

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(ADRB1_blood), 10)

# plot variable features with and without labels
plot1 <- VariableFeaturePlot(ADRB1_blood)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
save_plot(plot2, "HVG_ADRB1_blood", fig.width=10, fig.height=6)


ADRB1_blood <- ScaleData(ADRB1_blood, features = VariableFeatures(ADRB1_blood))

colnames(allABs.raw$ADRB1_blood)= paste("ADRB1_blood_", colnames(allABs.raw$ADRB1_blood), sep="")

obj.htos <- allABs.raw$ADRB1_blood[, colnames(ADRB1_blood)]
ADRB1_blood[["HTO"]] <- CreateAssayObject(counts = obj.htos)


ADRB1_blood <- NormalizeData(ADRB1_blood, assay = "HTO", normalization.method = "CLR")
ADRB1_blood <- HTODemux(ADRB1_blood, assay = "HTO", positive.quantile = 0.99)

table(ADRB1_blood$HTO_classification.global)


p= RidgePlot(ADRB1_blood, assay = "HTO", features = rownames(ADRB1_blood[["HTO"]]), ncol = 4)
save_plot(p, "HTOridge_ADRB1_blood", fig.width=30, fig.height=10)

#Singlet analysis
ADRB1_bloodHTO=ADRB1_blood[,ADRB1_blood$HTO_classification.global=="Singlet"]


ADRB1_bloodHTO <- ScaleData(ADRB1_bloodHTO)
ADRB1_bloodHTO= RunPCA(ADRB1_bloodHTO)
ADRB1_bloodHTO= RunUMAP(ADRB1_bloodHTO, dims = 1:50, reduction.key = "UMAP_")
ADRB1_bloodHTO <- FindNeighbors(ADRB1_bloodHTO, dims = 1:50)
ADRB1_bloodHTO <- FindClusters(ADRB1_bloodHTO, resolution = 0.2)


p=DimPlot(ADRB1_bloodHTO, shuffle = T, seed = 1, group.by= "HTO_classification", raster=FALSE)
save_plot(p, "singlets/wnn_ig_dimplot", 12, 8)
p=DimPlot(ADRB1_bloodHTO, group.by= "seurat_clusters", raster=FALSE)
save_plot(p, "singlets/wnn_pca_ig_dimplot", 12, 8)


####
deResTT = makeDEResults(ADRB1_bloodHTO, group.by="seurat_clusters", assay="RNA", test="t")
exprdfTT = getDEXpressionDF(ADRB1_bloodHTO, deResTT, assay="RNA", group.by="seurat_clusters")
write.table(exprdfTT,"singlets/expr_test_clusters_t.tsv", sep="\t", row.names=F, quote = F)
write_xlsx(exprdfTT, "singlets/expr_test_clusters_t.xlsx")


exprdfTT<-read_tsv("singlets/expr_test_clusters_t.tsv")
DefaultAssay(ADRB1_bloodHTO) <- "RNA"
markers.use.tt= subset(exprdfTT , avg_log2FC>0&p_val_adj<0.05&!startsWith(gene, "mt-")&!startsWith(gene, "rp"))
finalMarkers.use.tt = markers.use.tt %>% arrange(p_val_adj, desc(abs(pct.1)*abs(avg_log2FC))) %>% group_by(clusterID) %>% dplyr::slice(1:20)
finalMarkers.use.tt


data_dupli= finalMarkers.use.tt[!duplicated(finalMarkers.use.tt[ , "gene"]), ]


events= data_dupli %>% count(clusterID)
inser=cumsum(events$n)+0.5-events$n
insert=replace(inser, inser==0.5, 0)

xmi<- insert
xmin<- xmi[c(FALSE, TRUE)]
xma<- insert+events$n
xmax<- xma[c(FALSE, TRUE)]
ymi<- 0*events$n
ymin<- ymi[c(FALSE, TRUE)]
yma<- rep(length(events$n)+0.5, each=length(events$n))
ymax<- yma[c(FALSE, TRUE)]



p_dp_genes_idents = DotPlot(ADRB1_bloodHTO, features = data_dupli$gene, assay="RNA", dot.scale = 5, group.by="seurat_clusters")+
    theme(axis.text.x = element_text(angle = 90, hjust = 1.0, vjust = 0.5))+
    annotate("rect", xmin=xmin, xmax=xmax, ymin=ymin , ymax=ymax, alpha=0.2, fill="blue") #rep(c("blue", "grey"), times= length(events$n)/2)
  save_plot(p_dp_genes_idents, "singlets/dotplot_cluster_genes_colored", 45, 8)




#########################Volcanos
###### Separate cells by condition

cellList = colnames(ADRB1_bloodHTO)

featVec <- vector(mode="character", length=length(cellList))
featVec = ADRB1_bloodHTO$HTO_classification


featVec[featVec == "WT1"] = "WT"
featVec[featVec == "WT2"] = "WT"
featVec[featVec == "WT3"] = "WT"
featVec[featVec == "WT4"] = "WT"
featVec[featVec == "KO1"] = "KO"
featVec[featVec == "KO2"] = "KO"
featVec[featVec == "KO3"] = "KO"
featVec[featVec == "KO4"] = "KO"


ADRB1_bloodHTO$condition=featVec


##
save.image("ADRB1_christoph.Rdata")
load("/usr/local/hdd2/data/alejandro/Christoph_adrb1/seurat_ADRB1/ADRB1_christoph.Rdata")

############################################################################################### BHAVISHYA'S CODE
current.cluster.ids <- c(0, 1, 2, 3, 4, 5, 6, 7, 8)
new.cluster.ids <- c("B Cells", "Neutrophils", "Platelets", "T Cells", "Monocytes", "Eosinophils", "Erythrocytes", "NK Cells", "Basophils")
Idents(ADRB1_bloodHTO) <- c("B Cells", "Neutrophils", "Platelets", "T Cells", "Monocytes", "Eosinophils", "Erythrocytes", "NK Cells", "Basophils")

levels(ADRB1_bloodHTO$seurat_clusters) <- c("B Cells", "Neutrophils", "Platelets", "T Cells", "Monocytes", "Eosinophils", "Erythrocytes", "NK Cells", "Basophils")

dot_plot_features = readLines("/usr/local/hdd2/data/alejandro/Christoph_adrb1/seurat_ADRB1/dot_plot_features_christoph.csv")
dot_plot_features = unique(dot_plot_features)
Dot_plot_new = DotPlot(ADRB1_bloodHTO, features = dot_plot_features, assay="RNA", dot.scale = 5, group.by="seurat_clusters", cols = c("blue", "red"))+
    theme(axis.text.x = element_text(angle = 90, hjust = 1.0, vjust = 0.5))
save_plot(Dot_plot_new, "dot_plot_Adrb1_finalx", fig.width=15, fig.height=8)

violin_plot_features = readLines("/usr/local/hdd2/data/alejandro/Christoph_adrb1/seurat_ADRB1/violin_plot_features.csv")

firstup <- function(x) {
  substr(x, 1, 1) <- toupper(substr(x, 1, 1))
  x
}
violin_plot_features = readLines("/usr/local/hdd2/data/alejandro/Christoph_adrb1/seurat_ADRB1/violin_plot_features.csv")
violin_plot_features = tolower(violin_plot_features)
violin_plot_features = firstup(violin_plot_features)
violin_plot_features
Feature_plot= FeaturePlot(ADRB1_bloodHTO, features = violin_plot_features, split.by= "condition")
save_plot(Feature_plot, "Feature_plot", fig.width=12, fig.height=40)
Feature_violin= VlnPlot(ADRB1_bloodHTO, features = violin_plot_features, split.by= "condition", idents = "Neutrophils")
save_plot(Feature_violin, "Feature_violin_plotx", fig.width=20, fig.height=16)

genesetDF = as.data.frame(read.csv("/usr/local/hdd2/data/alejandro/Christoph_adrb1/seurat_ADRB1/sepsis_genes.csv", sep=","))
genesets = list()
identifiers_level_1 <- levels(ADRB1_bloodHTO$seurat_clusters)[ADRB1_bloodHTO$seurat_clusters == 1]
identifiers_level_1 <- as.character(ADRB1_bloodHTO$seurat_clusters[ADRB1_bloodHTO$seurat_clusters == levels(ADRB1_bloodHTO$seurat_clusters)[2]])
print(identifiers_level_1)
obj.integrated.patient_thr = ADRB1_bloodHTO
DefaultAssay(obj.integrated.patient_thr) = "RNA"

for (col in colnames(genesetDF))
{

  genes = genesetDF[, col]
  genes = genes[genes != ""]

  genesets[[col]] = unlist(lapply(genes, FUN=tolower))
  #genesets = firstup(genesets)
  print(paste(col, length(genesets[[col]])))
  print(setdiff(genesets[[col]], rownames(obj.integrated.patient_thr)))
}


genesets$neutrophil_extravasation = firstup(genesets$neutrophil_extravasation)
genesets$azurophilic_granule = firstup(genesets$azurophilic_granule)
genesets$gelatinase_granule = firstup(genesets$gelatinase_granule)
genesets$secretory_vesicle = firstup(genesets$secretory_vesicle)
genesets$specific_granule = firstup(genesets$specific_granule)
genesets$isg = firstup(genesets$isg)
genesets$rps = firstup(genesets$rps)
genesets$rpl = firstup(genesets$rpl)
genesets$NETosis = firstup(genesets$NETosis)
genesets$complement = firstup(genesets$complement)
genesets$transmigration = firstup(genesets$transmigration)
genesets$leukocyte = firstup(genesets$leukocyte)
genesets$Phagocytosis = firstup(genesets$Phagocytosis)

for (namex in names(genesets))
{
  namey = noquote(namex)
  #print(namex)
  genesets$namey = firstup(genesets$namey)
}


Idents(obj.integrated.patient_thr) <- c("B Cells", "Neutrophils", "Platelets", "T Cells", "Monocytes", "Eosinophils", "Erythrocytes", "NK Cells", "Basophils")

obj.integrated.patient_thr$Idents=Idents(obj.integrated.patient_thr)
levels(obj.integrated.patient_thr$seurat_clusters) <- c("B Cells", "Neutrophils", "Platelets", "T Cells", "Monocytes", "Eosinophils", "Erythrocytes", "NK Cells", "Basophils")

for (gsname in names(genesets))
{

  if (!dir.exists(dirname("genesets1")))
  {
    dir.create(dirname("genesets1"), recursive = TRUE)
  }

  print(gsname)
  print(genesets[[gsname]])

  print("missing")
  print(setdiff(genesets[[gsname]], rownames(obj.integrated.patient_thr)))

  DefaultAssay(obj.integrated.patient_thr) = "RNA"

  obj.integrated.patient_thr <- AddModuleScore(obj.integrated.patient_thr,
                    features = list(genesets[[gsname]]),
                    name=gsname)

  baseplotname = paste("/usr/local/hdd2/data/alejandro/Christoph_adrb1/seurat_ADRB1/genesets1/", gsname, sep="")

  ffname = paste(gsname, 1, sep="")
  split.by="condition"
  group.by="Idents"

  splitFeaturePlot(obj.integrated.patient_thr, feature = ffname, split.by=split.by, title=paste("Module Score", gsname), filename=paste(baseplotname, "split", sep="_"), limits=c(-2, 2), low="#713fdf", high="#df713f",  mid="grey")

  p=comparativeVioBoxPlot(obj.integrated.patient_thr, ffname, group.by, split.by)+ylab(paste("Module Score", gsname)) 
  save_plot(p, paste(baseplotname, "violent", sep="_"), fig.height=6, fig.width=20)

  px=comparativeVioBoxPlot(subset(obj.integrated.patient_thr, Idents == "Neutrophils"), ffname, group.by, split.by)+ylab(paste("Module Score", gsname)) 
  save_plot(p, paste(baseplotname, "neutro_violent", sep="_"), fig.height=6, fig.width=20)

  p=VlnPlot(subset(obj.integrated.patient_thr, Idents == "Neutrophils"), ffname, group.by=group.by)+ylab(paste("Module Score", gsname)) 
  save_plot(p, paste(baseplotname, "neutro_violin", sep="_"), fig.height=4, fig.width=10)

  p=comparativeVioBoxPlot(subset(obj.integrated.patient_thr, Idents == "Neutrophils"), ffname, group.by, split.by, onelineLabel=TRUE, yStepIncrease=0.75)+
  ylab(paste("Module Score", gsname))+ylim(c(0,10))
  save_plot(p, paste(baseplotname, "neutro_violent_pooled", sep="_"), fig.height=4, fig.width=10)

  p = FeaturePlot(subset(obj.integrated.patient_thr, Idents == "Neutrophils"), features=ffname )+ggtitle(paste("Module Score", gsname))
  save_plot(p, paste(baseplotname, "neutro_featureplot", sep="_"), fig.height=10, fig.width=12)


}

library(gridExtra)
all_plots <- list()
library(cowplot)

for (gsname in names(genesets))
{

  

  print(gsname)
  print(genesets[[gsname]])

  print("missing")
  print(setdiff(genesets[[gsname]], rownames(obj.integrated.patient_thr)))

  DefaultAssay(obj.integrated.patient_thr) = "RNA"

  obj.integrated.patient_thr <- AddModuleScore(obj.integrated.patient_thr,
                    features = list(genesets[[gsname]]),
                    name=gsname)

  baseplotname = paste("/usr/local/hdd2/data/alejandro/Christoph_adrb1/seurat_ADRB1/genesets1/", gsname, sep="")

  ffname = paste(gsname, 1, sep="")
  split.by="condition"
  group.by="Idents"

  
  p=comparativeVioBoxPlot(subset(obj.integrated.patient_thr, Idents == "Neutrophils"), feature=ffname, group.by=group.by, split.by=split.by)+ylab(paste("Module Score", gsname)) 
  #save_plot(p, paste(baseplotname, "neutro_violent", sep="_"), fig.height=6, fig.width=20)
  all_plots[[gsname]] <- p

}

combined_plot <- plot_grid(plotlist = all_plots, ncol = 6)

# Save the combined plot as an image
ggsave(filename = "/usr/local/hdd2/data/alejandro/Christoph_adrb1/seurat_ADRB1/genesets1/all_plots_new.svg", plot = combined_plot, height = 9, width = 16, units = "in")



rstatix::pairwise_t_test(as.data.frame(ptest),
        as.formula(paste("azurophilic_granule1", " ~ ", "Idents", sep="")), paired = FALSE, 
        p.adjust.method = "BH"
      )







score_analysis = read.csv("/usr/local/hdd2/data/alejandro/Christoph_adrb1/seurat_ADRB1/sepsis_genes.csv")
# Basic function to convert human to mouse gene names
install.packages("biomaRt")
library(biomaRt)
convertHumanGeneList <- function(x){
require("biomaRt")
human = useMart("ensembl", dataset = "hsapiens_gene_ensembl")
mouse = useMart("ensembl", dataset = "mmusculus_gene_ensembl")
genesV2 = getLDS(attributes = c("hgnc_symbol"), filters = "hgnc_symbol", values = x , mart = human, attributesL = c("mgi_symbol"), martL = mouse, uniqueRows=T)
humanx <- unique(genesV2[, 2])
# Print the first 6 genes found to the screen
print(head(humanx))
return(humanx)
}
genes <- convertHumanGeneList(humGenes)

tolower()
#################################################################################################

Feat= FeaturePlot(ADRB1_bloodHTO, features = "Adrb1", split.by= "condition")
save_plot(Feat, "Featureplot_Adrb1", fig.width=12, fig.height=8)

Feat= FeaturePlot(ADRB1_bloodHTO, features = "Adrb2", split.by= "condition")
save_plot(Feat, "Featureplot_Adrb2", fig.width=12, fig.height=8)

Feat= FeaturePlot(ADRB1_bloodHTO, features = "Adrb1", split.by= "condition")
save_plot(Feat, "singlets/FeatPlot_Adrb1", fig.width=12, fig.height=8)

Feat= FeaturePlot(ADRB1_bloodHTO, features = "Itgam")
save_plot(Feat, "ViolinPlot_monocytes", fig.width=12, fig.height=8)


Neut=ADRB1_bloodHTO[,ADRB1_bloodHTO$condition=="WT"]
Neut=Neut[,Neut$seurat_clusters=="1"]

# Extract data for the desired features
features <- c("Adra1a", "Adra1b", "Adra1d", "Adra2a", "Adra2b", "Adra2c", "Adrb1", "Adrb2", "Adrb3")
Neu <- as.data.frame(t(as.matrix(Neut[["RNA"]]@data[features, ])))

# Add cell barcodes (rownames) as a column
Neu_df <- Neu
Neu_df["ref_names"] <- rownames(Neu_df)

# Reshape the data for ggplot2
library(reshape2)
Neu_df <- melt(Neu_df, id.vars = "ref_names")

# Create a violin plot
library(ggplot2)
p <- ggplot(Neu_df, aes(x = variable, y = value, fill = variable)) +
  geom_violin(trim = FALSE) +
  theme_minimal() +
  labs(x = "Genes", y = "Expression", title = "Violin Plot of Adrenoreceptor Genes") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Save the plot
save_plot(p, "Violin_all_adrenoreceptors", fig.width=8, fig.height=6)



Feat= FeaturePlot(ADRB1_bloodHTO, features = "Itgam")
save_plot(Feat, "FeatPlot_monocytes", fig.width=12, fig.height=8)

p=DimPlot(ADRB1_bloodHTO, shuffle = T, seed = 1, group.by= "condition", raster=FALSE)
save_plot(p, "singlets/wnn_ig_dimplot_condition", 12, 8)


####

deGroup = "condition"


allConds = unique(ADRB1_bloodHTO[[deGroup]][[deGroup]])


print("all conds")
print(allConds)

all.cells = list()

for (cond in allConds)
{

  all.cells[[cond]] = cellIDForClusters(ADRB1_bloodHTO, deGroup, c(cond))

}


for (i in 1:(length(allConds)-1))
{

  for (j in (i+1):length(allConds))
  {
    print(paste(i,"<->",j))

    condI = allConds[i]
    condJ = allConds[j]

    condNameI = tolower(str_replace_all(condI, " ", "_"))
    condNameJ = tolower(str_replace_all(condJ, " ", "_"))

    print(paste(condJ, condNameJ))
    print(paste(condI, condNameI))

    deMarkers= de.condi.condj = compareCellsByCluster(ADRB1_bloodHTO, all.cells[[condJ]], all.cells[[condI]], condNameJ, condNameI,
                                                outfolder=paste("singlets", "de", sep="/"), fcCutoff=0.1)


    deMarkersW= dewilcox.condi.condj = compareCellsByCluster(ADRB1_bloodHTO, all.cells[[condJ]], all.cells[[condI]], condNameJ, condNameI,
                                                    outfolder=paste("singlets", "dewilcox", sep="/"), test="wilcox", fcCutoff=0.1)

    compName = paste(condNameJ, condNameI, sep="_")

    saveRDS(de.condi.condj, file=paste("singlets", "de", paste("der_", compName, ".rds", sep=""), sep="/"))
    saveRDS(dewilcox.condi.condj, file=paste("singlets", "dewilcox", paste("der_", compName, ".rds", sep=""), sep="/"))
    
    makeVolcanos(de.condi.condj, paste("DE", condJ, "vs", condI), paste("singlets", "de_volcano", compName, sep="/"),
                turnExpression=F, FCcutoff=0.1, pCutoff = 0.05)

    makeVolcanos(dewilcox.condi.condj, paste("DE", condJ, "vs", condI), paste("singlets", "dewilcox_volcano", compName, sep="/"),
                turnExpression=F, FCcutoff=0.1, pCutoff = 0.05)
    }

}


####### Neutrophils

Neut=ADRB1_bloodHTO[,ADRB1_bloodHTO$seurat_clusters=="1"]

Neut <- ScaleData(Neut)
Neut= RunPCA(Neut)
Neut= RunUMAP(Neut, dims = 1:50, reduction.key = "UMAP_")
Neut <- FindNeighbors(Neut, dims = 1:50)
Neut <- FindClusters(Neut, resolution = 0.2)


p=DimPlot(Neut, group.by= "condition", raster=FALSE)
save_plot(p, "Neutrophils/UMAP_condition", 12, 8)
p=DimPlot(Neut, group.by= "seurat_clusters", raster=FALSE)
save_plot(p, "Neutrophils/UMAP_clusters", 12, 8)
