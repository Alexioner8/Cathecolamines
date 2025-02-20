library(Tempora)
library(Seurat)
library(ggplot2)
library(cowplot)
library(stringr)

library(dplyr)
library(tidyverse)
library(pbapply)

library(clusterProfiler)
library(enrichplot)
library(data.table)
library(writexl)
library(EnhancedVolcano)
library(ggsunburst)
library(ggrepel)
library(RColorBrewer)
library(svglite)
library(ggpubr)
library(gtools)

library("purrr")
library("dplyr")

splitObjListByGroup = function(objlist, group.by)
{
  finalList = list()
  for (objname in names(objlist))
  {

      if (objname %in% names(objlist))
      {
          for (abTag in unique(objlist[[objname]][[group.by]][[group.by]]))
          {
              sampleName = paste(objname, abTag, sep="_")
              print(sampleName)

              selDF = objlist[[objname]][[group.by]]
              selCells = rownames(selDF)[selDF[[group.by]] == abTag]

              finalList[[sampleName]] = subset(objlist[[objname]], cells = selCells)
              finalList[[sampleName]]$orig_project = objname
              finalList[[sampleName]]$library = paste(objname, abTag, sep="_")
              
              print(sampleName)
              print(finalList[[sampleName]])
          }

      } else {
          finalList[[objname]] = objlist[[objname]]

          print(objname)
          print(finalList[[objname]])
      }

  }

  return(finalList)
}


splitObjListByHTO = function(objlist)
{
finalList = list()
for (objname in names(objlist))
{

    if (objname %in% names(objlist))
    {
        for (abTag in unique(objlist[[objname]]$HTO_classification))
        {
            sampleName = paste(objname, abTag, sep="_")
            print(sampleName)

            finalList[[sampleName]] = subset(objlist[[objname]], HTO_classification==abTag)
            finalList[[sampleName]]$orig_project = objname
            print(sampleName)
            print(finalList[[sampleName]])
        }

    } else {
        finalList[[objname]] = objlist[[objname]]

        print(objname)
        print(finalList[[objname]])
    }

}

return(finalList)
}

splitObjListByLibrary = function(objlist)
{
    finalList = list()
    for (objname in names(objlist))
    { 

        xobj = objlist[[objname]]
        xobj$orig_project = objname
        finalList[[objname]] = xobj

    }

    return(finalList)
}



combine_plot_grid = function(...)
{
  inplot = list(...)
  print(inplot)
  
  dataList = list()
  for (i in 1:length(inplot))
  {
    dataList[[length(dataList)+1]] = inplot[[i]]$data
  }
  
  p=cowplot::plot_grid(...)
  
  p$data = dataList
  
  return(p)
}

combine_plot_grid_list = function(plotlist, ...)
{
 
  dataList = list()
  for (i in 1:length(plotlist))
  {
    dataList[[i]] = plotlist[[i]]$data
  }
  
  p=cowplot::plot_grid(plotlist=plotlist, ...)
  
  p$data = dataList
  
  return(p)
}



save_plot = function(plotobj, outname, fig.width, fig.height)
{
  print(paste(outname, fig.width, fig.height))
  
  fname=paste(outname, "png", sep=".")
  print(paste("Saving to file", fname))
  png(filename=fname, width = fig.width, height = fig.height, units = 'in', res = 300)#width = fig.width*100, height=fig.height*100)
  plot(plotobj)
  dev.off()
  
  fname=paste(outname, "pdf", sep=".")
  print(paste("Saving to file", fname))
  pdf(file=fname, width = fig.width, height=fig.height)
  plot(plotobj)
  dev.off()
  

  fname=paste(outname, "svg", sep=".")
  print(paste("Saving to file", fname))
  svglite::svglite(file = fname, width = fig.width, height = fig.height)
  plot(plotobj)
  dev.off()
  
  if (class(plotobj$data) %in% c("list"))
  {
    print("list case")
    for (i in 1:length(plotobj$data))
    {
      fname = paste(outname,i, "data", sep=".")
      print(paste("Saving to file", fname))
      
      if (class(plotobj$data[[i]]) %in% c("list"))
      {
          print("multi list case")
          for (j in 1:length(plotobj$data[[i]]))
          {
              fname = paste(outname,i, j, "data", sep=".")
              print(paste("Saving to file", fname, class(plotobj$data[[i]][[j]])))

              if (class(plotobj$data[[i]][[j]]) %in% c("list", "waiver"))
              {
                next()
              }
              write.table(plotobj$data[[i]][[j]], fname, row.names = TRUE, sep="\t")    

          }
      } else {
          
          tryCatch(write.table(plotobj$data[[i]], fname, row.names = TRUE, sep="\t"), error = function(e) NULL)
      }
      
    }
  } else {
    
      fname = paste(outname,"data", sep=".")
      print(paste("Saving to file", fname))

      write.table(plotobj$data, paste(outname, "data", sep="."), row.names = TRUE, sep="\t")
  }
  
  return(plotobj)
}

write_cell_barcodes = function(scobj, outfolder, valid_orig_idents=NULL)
{
    dir.create(outfolder)

    orig_idents = unique(scobj$orig.ident)

    if (!is.null(valid_orig_idents))
    {
      orig_idents = intersect(valid_orig_idents, orig_idents)
      print("Performing writeout on orig_idents")
      print(orig_idents)
    }

    for (oident in orig_idents)
    {
      print(oident)
      scobjs = subset(scobj, orig.ident == oident)

      allCellNames = colnames(scobjs)

      cellBarcodes = unlist(lapply(strsplit(unlist(lapply(strsplit(colnames(scobjs), "_(?!.*_)", perl = TRUE), function(x){x[2]})), "-"), function(x){x[1]}))

      outfilename = paste(outfolder, paste(oident, "barcodes.tsv", sep="_"), sep="/")
      write.table(cellBarcodes, outfilename, sep="\t", row.names=FALSE, col.names=FALSE, quote=F)
    }

}


toCellPrefix = function(x) {
    datasetprefix = x
    datasetprefix = str_replace_all(datasetprefix, "[./]", "")
    datasetprefix = str_split(datasetprefix, "_")[[1]][1]
    
    return(paste(datasetprefix, sep=""))
}

makesum = function(a, suffix)
{
  out = {}
  out["sum"] = sum(a)

  return(out)
}

patternList.human = list()
patternList.human[["MT"]] = "^MT-"
patternList.human[["RPL"]] = "^RPL"
patternList.human[["RPS"]] = "^RPS"


patternList.mouse = list()
patternList.mouse[["MT"]] = "^mt-"
patternList.mouse[["RPL"]] = "^Rpl"
patternList.mouse[["RPS"]] = "^Rps"


makeSeuratObj = function(matrix, proj, pl)
{
    obj = CreateSeuratObject(matrix, project=proj)
    print("Renaming Cells")
    obj <- RenameCells(obj, add.cell.id=proj)
    
    print(paste("Seurat obj project", obj@project.name))
    
    mtPattern = pl[["MT"]]
    rplPattern = pl[["RPL"]]
    rpsPattern = pl[["RPS"]]
    rpPattern = paste(c(pl[["RPL"]], pl[["RPS"]]), sep="", collapse="|")

    
    selGenes = rownames(obj)[grepl(rownames(obj), pattern=mtPattern)]
    print(paste("Got a total of mt-Genes:", length(selGenes), paste(head(selGenes), collapse=", ")))
    
    selGenes = rownames(obj)[grepl(rownames(obj), pattern=rplPattern)]
    print(paste("Got a total of Rpl-Genes:", length(selGenes), paste(head(selGenes), collapse=", ")))
    
    selGenes = rownames(obj)[grepl(rownames(obj), pattern=rpsPattern)]
    print(paste("Got a total of Rps-Genes:", length(selGenes), paste(head(selGenes), collapse=", ")))
    
    selGenes = rownames(obj)[grepl(rownames(obj), pattern=rpPattern)]
    print(paste("Got a total of Rp-Genes:", length(selGenes), paste(head(selGenes), collapse =", ")))
    
    obj[["percent.mt"]] <- PercentageFeatureSet(obj, pattern = mtPattern)
    obj[["percent.rpl"]] <- PercentageFeatureSet(obj, pattern = rplPattern)
    obj[["percent.rps"]] <- PercentageFeatureSet(obj, pattern = rpsPattern)
    obj[["percent.rp"]] <- PercentageFeatureSet(obj, pattern = rpPattern)
    
    return(obj)
}












convertHumanGeneList <- function(x){
require("biomaRt")
httr::set_config(httr::config(ssl_verifypeer = FALSE))
human = useEnsembl("ensembl", dataset = "hsapiens_gene_ensembl")
mouse = useEnsembl("ensembl", dataset = "mmusculus_gene_ensembl")
genesV2 = getLDS(attributes = c("hgnc_symbol"), filters = "hgnc_symbol", values = x , mart = human, attributesL = c("mgi_symbol"), martL = mouse, uniqueRows=T)
humanx <- unique(genesV2[, 2])
# Print the first 6 genes found to the screen
print(head(humanx))
return(humanx)
}

getCellCountDF = function(scobj, prefix="", group_by="orig.ident", select_by="idents", split_by=NULL, relative=F, outname=NULL,show.percent=F)
{
  
  allClusters = as.character(sort(unique(scobj[[select_by]][,])))
  
  allSplits = NULL
  if (!is.null(split_by))
  {
    allSplits = as.character(sort(unique(scobj[[split_by]][,])))
  }
  
  cellCounts = list()
  
  print("All Clusters")
  print(allClusters)
  print("All Splits")
  print(allSplits)


  for (clusterID in allClusters)
  {
    print(clusterID)
    
    if (!is.null(split_by))
    {
      clusterList = list()
      clusterList[["cluster"]] = clusterID
      cs = scobj[,scobj[[select_by]] == clusterID]
      clusterList[["all"]] = nrow(cs[[group_by]])
      for (split_cat in allSplits)
      {
          print(paste(clusterID, split_cat))
        
          cs = scobj[,scobj[[select_by]] == clusterID & scobj[[split_by]] == split_cat]
          allElems = table(cs[[group_by]])
          cs = NULL
        
          for (grp in names(allElems))
          {
            if (!is.null(prefix))
            {
              ngrp = paste(prefix, paste(split_cat, grp, sep="."), sep=".")
  
            } else {
              ngrp = paste(split_cat, grp, sep=".")
            }
            clusterList[[ngrp]] = allElems[[grp]]
          }
      }
      
    } else {
      
        cs = scobj[,scobj[[select_by]] == clusterID]
      
        allElems = table(cs[[group_by]])

        clusterList = list()
        clusterList[["cluster"]] = clusterID
        clusterList[["all"]] = nrow(cs[[group_by]])
        cs = NULL
        
        for (grp in names(allElems))
        {
          if (!is.null(prefix))
          {
            ngrp = paste(prefix, grp, sep=".")
          } else {
            ngrp = grp
          }
          clusterList[[ngrp]] = allElems[[grp]]
        }
    }
    
    
    
    cellCounts[[clusterID]] = clusterList
    
  }
  df_bysamplerep = cellCounts %>% map(as.data.frame) %>% bind_rows()
  df_bysamplerep[is.na(df_bysamplerep)] <- 0
  
  rownames(df_bysamplerep) = df_bysamplerep$cluster
  df_bysamplerep$cluster = NULL
  
  
  if (relative)
  {
    df_bysamplerep = sweep(df_bysamplerep,2,colSums(df_bysamplerep),"/")
    
    if (show.percent)
    {
      df_bysamplerep = df_bysamplerep*100;
    }
  }
  
  totals=t(colSums(df_bysamplerep))
  totals.df = data.frame(totals)
  rownames(totals.df) = "Total"
  df_bysamplerep=rbind(df_bysamplerep, totals.df)
  
  df_bysamplerep = cbind("cluster"=rownames(df_bysamplerep), df_bysamplerep)
  rownames(df_bysamplerep) = NULL
  
  if (!is.null(outname))
  {
    write.table(df_bysamplerep, file=outname, row.names = F,  quote=FALSE, sep='\t')
    write_xlsx( df_bysamplerep, path = paste(outname, ".xlsx", sep="") )

  }
  
  return(df_bysamplerep)
  #  
}



makeBarCountPlot = function( seuratObj, outpath, select_by, group_by)
{
countByManualCellname = getCellCountDF(seuratObj, prefix="", select_by = select_by, group_by=group_by, relative=TRUE, show.percent=T, outname=paste(outpath, "tsv", sep="."))

cbc_cname = melt(countByManualCellname)
cbc_cname$variable_factors = factor(cbc_cname$variable)#, levels=c("all", ".FIRE_NI_CTRL",".FIRE_NI_KO"))
cbc_cname = cbc_cname[!cbc_cname$cluster=="Total",]

p=DimPlot(seuratObj, group.by = select_by)
g <- ggplot_build(p)
#unique(g$data[[1]]["colour"])
#cellname2color = data.frame(group=g$data[[1]]$group, color=g$data[[1]]$colour)# , cellnames=seuratObj[[group_by]][g$data[[1]]$group])

cbc_cname$cluster = factor(cbc_cname$cluster, levels=levels(p$data[[select_by]]))
cbc_cname$total = 100
cbc_cname$label = paste0(cbc_cname$cluster, " (", round(cbc_cname$value, 2), "% )")

numColumns = length(unique(cbc_cname$variable))

print(cbc_cname)

plot_distribution_celltype = cbc_cname[cbc_cname$variable != "all",] %>%
ggplot(aes(variable_factors, value, fill=cluster)) + 
geom_col(width=0.75)+#scale_fill_manual()+
geom_text_repel(aes(label = paste(round(value, digits = 1))),  size=10, position = position_stack(vjust = .5), direction ="x", segment.size  = 1, segment.color = "black", min.segment.length=0, arrow = arrow(length = unit(0.015, "cm"), angle = 0, type = "closed", ends = "first"), force=2)+
  ylab("Percentage")+xlab("Condition")+
  theme(
        legend.text = element_text(size = 18),
        axis.text.x = element_text(angle = 0, hjust=0.5, vjust=0.5, size=18),
        axis.text.y=element_blank(),
        axis.ticks=element_blank(),
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        panel.background=element_blank(),
        panel.border=element_blank(),
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank())


use_width = numColumns * 3

save_plot(plot_distribution_celltype, outpath, fig.height=16, fig.width=use_width)


}


makesummary_getPop = function(a, suffix)
{
out = {}
out["num"] = length(a)

if (length(a) == 0)
{
    f = c(0,0,0,0,0)
    meanA = 0
} else {
    f = fivenum(a)
    meanA = mean(a)
}

out["min"] = f[1]
out["lower_hinge"] = f[2]
out["median"] = f[3]
out["upper_hinge"] = f[4]
out["max"] = f[5]
out["mean"] = meanA

names(out) = paste(names(out), suffix, sep=".")

return(out)
}

getExprData_getPop = function(markerObj, markerCells, sampleSuffix, slot="data", assay="RNA")
{
expTable = GetAssayData(object = subset(x=markerObj, cells=markerCells), slot = slot, assay=assay)
allgenes = rownames(expTable)
cellnames = colnames(expTable)

expt.r = as(expTable, "TsparseMatrix")
expt.df = data.frame(r = expt.r@i + 1, c = expt.r@j + 1, x = expt.r@x)

DT <- data.table(expt.df)
res = DT[, as.list(makesummary_getPop(x, sampleSuffix)), by = r]
anumCol = paste("anum", sampleSuffix, sep=".")
res[[anumCol]] = length(cellnames)
res$gene = allgenes[res$r]

res = res[,r:=NULL]

return(res)
}

getDEXpressionDF = function ( scdata, markers, assay="SCT", group.by=NULL)
{
outDF = NULL
DefaultAssay(object=scdata) = assay
print(group.by)
if (is.null(group.by))
{
  clusterIDs = as.character(sort(unique(Idents(scdata))))
} else {
  clusterIDs = as.character(sort(unique(scdata[[group.by]][,])))
}
scCells = Idents(scdata)
scCells = names(scCells)
scCells = unlist(as.character(scCells))
for (clusterID in clusterIDs){
    
    print(clusterID)
    
    cellIdents = Idents(scdata)
    
    if (is.null(group.by))
    {
      cellIdents.c = names(cellIdents[cellIdents == clusterID])
    } else {
      cellIdents.c = colnames(scdata[,scdata[[group.by]] == clusterID])
    }
    
    cellIdents.c = unlist(lapply(cellIdents.c, as.character))  
    
    cellIdents.bg = setdiff(unlist(lapply(names(cellIdents), as.character)), cellIdents.c)

    if (length(cellIdents.c) < 3)
    {
      print(paste("Skipping cluster", clusterID, "due to < 3 cells!"))  
      next
    }
    
    expvals = getExprData_getPop(scdata, cellIdents.c, "cluster", assay=assay)
    expvals.bg = getExprData_getPop(scdata, cellIdents.bg, "bg", assay=assay)
    modmarkers = markers[[clusterID]]
    modmarkers$gene = rownames(modmarkers)
    
    markerdf = as.data.frame(modmarkers)
    
    if ((nrow(markerdf) > 0) && (nrow(expvals) > 0))
    {
      expvals = merge(markerdf, expvals, all.x=T, by.x="gene", by.y = "gene")  
    }
    
    if ((nrow(expvals) > 0) && (nrow(expvals.bg) > 0))
    {
      expvals = merge(expvals, expvals.bg, all.x=T, by.x="gene", by.y = "gene")  
    }
    
    expvals = as.data.frame(cbind(clusterID, expvals))
    
    if (!is.data.frame(outDF) || nrow(outDF)==0)
    {
    outDF = expvals
    } else {
    outDF = as.data.frame(rbind(outDF, expvals))
    }
    
}
return(outDF)
}



makeDEResults = function(inobj, group.by=NULL, assay="SCT", test="wilcox")
{
  if (is.null(group.by))
  {
    clusterIDs = as.character(sort(unique(Idents(inobj))))
  } else {
    clusterIDs = as.character(sort(unique(inobj[[group.by]][,])))
  }
  
  retList = list()
  
  for(clusterID in clusterIDs)
  {
  
  
      cellIdents = Idents(inobj)
      
      if (is.null(group.by))
      {
        cellIdents.c = names(cellIdents[cellIdents == clusterID])
  
      } else {
        cellIdents.c = colnames(inobj[,inobj[[group.by]] == clusterID])
      }
      
      cellIdents.c = unlist(lapply(cellIdents.c, as.character))
  
      print(paste("Processing cluster", clusterID, "with a total of", length(cellIdents.c), "cells"))
      print(length(cellIdents.c) < 3)
  
      if (length(cellIdents.c) < 3)
      {
        print(paste("Skipping cluster", clusterID, "due to < 3 cells!"))  
        next
      }
      deMarkers = FindMarkers(inobj, assay=assay, ident.1 = cellIdents.c, test.use=test) 
      retList[[clusterID]] = deMarkers
  
  }
  
  return(retList)
}



get_population_expression_data = function(scobj, group, outname, assay="RNA", slot="counts")
{

exprData = list()

for (cellPop in unique(as.character(unlist(scobj[[group]]))))
{
    varPop = str_to_lower( str_replace_all(
                            str_replace_all(#
                            str_replace_all( cellPop, "\\(|\\)| |,", "_"),
                            "__", "_"),
                            "_$", "")
                        )
    print(paste(cellPop, varPop))
    
    allPopCells = scobj[[group]]
    allPopCells$cnames = rownames(allPopCells)
    cellPopCells = allPopCells[allPopCells[[group]] == cellPop, ]$cnames
    print(paste("Number of cells: ", length(cellPopCells)))

    exprData[[varPop]] = getExprData_getPop(markerObj=scobj, markerCells=cellPopCells, sampleSuffix=varPop, slot=slot, assay=assay)
}


meanExprData = list()

for (name in names(exprData))
{
    
    exprDf = as.data.frame(exprData[[name]])
    subdf = exprDf[ ,c("gene", paste("mean", name, sep=".")) ]

    meanExprData[[name]] = subdf
}

cellnames_manualExprDF = Reduce(function(x,y) merge(x = x, y = y, by = "gene", all.x=T, all.y=T), meanExprData)
cellnames_manualExprDF[is.na(cellnames_manualExprDF)] = 0

write.table(cellnames_manualExprDF, file = paste(outname, ".tsv", sep=""), quote=FALSE, sep = "\t", row.names = F)
write_xlsx( cellnames_manualExprDF, path = paste(outname, ".xlsx", sep="") )

return(cellnames_manualExprDF)


}

cellIDForClusters = function(obj.in, targetVar, clusters)
{

  targetVarDF = as.data.frame(obj.in[[targetVar]])
  print(paste("orig:", length(rownames(targetVarDF))))
  cellNames = rownames(targetVarDF)[targetVarDF[[targetVar]] %in% clusters]

  print(length(cellNames))
  return(cellNames)

}


compareClusters = function(scdata, cellsID1, cellsID2, suffix1, suffix2, prefix="cluster", test="MAST", assay="RNA", outfolder="./", fcCutoff=0.25, all=FALSE)
{
    logfc.threshold = fcCutoff
    
    if (all==TRUE)
    {
    logfc.threshold = 0.01  
    }

    if (!dir.exists(outfolder)){
        dir.create(outfolder)
        print(outfolder)
        print("Dir created!")

    } else {
        print(outfolder)
        print("Dir already exists!")
    }
    
    markers = FindMarkers(scdata, assay=assay, ident.1 = cellsID1, ident.2 = cellsID2, test.use=test, logfc.threshold=logfc.threshold)
    
    outvalues1 = getExprData_getPop(scdata, cellsID1, suffix1, assay=assay)
    outvalues2 = getExprData_getPop(scdata, cellsID2, suffix2, assay=assay) 
    
    
    markers$gene = rownames(markers)
    joinedData = merge(markers, outvalues1, by="gene", all=T)
    joinedData = merge(joinedData, outvalues2, by="gene", all=T)  
    
    joinedData = joinedData[!is.na(joinedData$p_val),]
    
    suffix1=str_replace_all(str_replace_all(suffix1, "\\/", "_"), " ", "_")
    suffix2=str_replace_all(str_replace_all(suffix2, "\\/", "_"), " ", "_")

    
    outfile = paste(outfolder, "/", prefix, ".", suffix1, "_", suffix2, ".tsv", sep="")
    
    message(outfile)
    write.table(joinedData, file=outfile, row.names = F,  quote=FALSE, sep='\t')
    
    outfile = paste(outfolder, "/", prefix, ".", suffix1, "_", suffix2, ".xlsx", sep="")
    
    message(outfile)
    write_xlsx(joinedData, path=outfile)
    
    return(joinedData)
}

compareCellsByCluster = function(inobj, cellsID1, cellsID2, suffix1, suffix2, group.by=NULL, assay="RNA", test="t", outfolder="./", fcCutoff=0.25, all=FALSE)
{
  if (is.null(group.by))
  {
    clusterIDs = as.character(sort(unique(Idents(inobj))))  
  } else {
    clusterIDs = as.character(sort(unique(inobj[[group.by]][,])))
  }
  
  retList = list()
  
  for(clusterID in clusterIDs)
  {
  
      cellIdents = Idents(inobj)
      
      if (is.null(group.by))
      {
        cellIdents.c = names(cellIdents[cellIdents == clusterID])
  
      } else {
        cellIdents.c = colnames(inobj[,inobj[[group.by]] == clusterID])
      }
      
      cellIdents.c = unlist(lapply(cellIdents.c, as.character))
  
      print(paste("Processing cluster", clusterID, "with a total of", length(cellIdents.c), "cells"))
      
      cells1 = intersect(cellsID1, cellIdents.c)
      cells2 = intersect(cellsID2, cellIdents.c)
      
      print(c(length(cells1), length(cells2)))
      
      if (length(cells1) < 3)
      {
        print("Cells1 too few")
        next
      }
      
      if (length(cells2) < 3)
      {
        print("Cells2 too few")
        next
      }

      clusterID_file=str_replace_all(str_replace_all(clusterID, "\\/", "_"), " ", "_")

  
      deMarkers = compareClusters(scdata=inobj,
                                  cellsID1=cells1,
                                  cellsID2=cells2,
                                  prefix= paste("cluster", clusterID_file, sep="_"),
                                  suffix1=suffix1,#paste(suffix1, clusterID, sep="_"),
                                  suffix2=suffix2, #paste(suffix2, clusterID, sep="_"),
                                  test=test, fcCutoff=fcCutoff, assay=assay, outfolder=outfolder)
  
  
      retList[[clusterID]] = deMarkers
  
  }
  
  return(retList)
}


compareCellsByClusternew = function(inobj, cellsID1, cellsID2, suffix1, suffix2, group.by=NULL, assay="RNA", test="t", outfolder="./", fcCutoff=0.25, all=FALSE)
{
  if (is.null(group.by))
  {
    clusterIDs = as.character(c(0,1,3))  
  } else {
    clusterIDs = as.character(c(0,1,3))
  }
  
  retList = list()
  
  for(clusterID in clusterIDs)
  {
  
      cellIdents = Idents(inobj)
      
      if (is.null(group.by))
      {
        cellIdents.c = names(cellIdents[cellIdents == clusterID])
  
      } else {
        cellIdents.c = colnames(inobj[,inobj[[group.by]] == clusterID])
      }
      
      cellIdents.c = unlist(lapply(cellIdents.c, as.character))
  
      print(paste("Processing cluster", clusterID, "with a total of", length(cellIdents.c), "cells"))
      
      cells1 = intersect(cellsID1, cellIdents.c)
      cells2 = intersect(cellsID2, cellIdents.c)
      
      print(c(length(cells1), length(cells2)))
      
      if (length(cells1) < 3)
      {
        print("Cells1 too few")
        next
      }
      
      if (length(cells2) < 3)
      {
        print("Cells2 too few")
        next
      }

      clusterID_file=str_replace_all(str_replace_all(clusterID, "\\/", "_"), " ", "_")

  
      deMarkers = compareClusters(scdata=inobj,
                                  cellsID1=cells1,
                                  cellsID2=cells2,
                                  prefix= paste("cluster", clusterID_file, sep="_"),
                                  suffix1=suffix1,#paste(suffix1, clusterID, sep="_"),
                                  suffix2=suffix2, #paste(suffix2, clusterID, sep="_"),
                                  test=test, fcCutoff=fcCutoff, assay=assay, outfolder=outfolder)
  
  
      retList[[clusterID]] = deMarkers
  
  }
  
  return(retList)
}








makeVolcanos = function(loMG, titlePrefix, outname, restrict_labels=NULL, turnExpression=F, colors = list(neg=list(sig="#448CCA", nosig="#B4D1E9"), pos=list(sig="#F47B78", nosig="#FACAC9")), FCcutoff=0.5, pCutoff = 0.05)
{

    outfolder = dirname(outname)[1]

    if (!dir.exists(outfolder)){
        dir.create(outfolder)
        print(outfolder)
        print("Dir created!")

    } else {
        print(outfolder)
        print("Dir already exists!")
    }
  
  for (cName in names(loMG))
  {
    print(cName)
    
    title=paste(titlePrefix, "(", cName, ")", sep=" ")
    subtitle=""
    
    
    indf = loMG[[cName]]
    
    cName = str_replace_all(str_replace_all(cName, "\\/", "_"), " ", "_")
    popName = str_to_lower( str_replace_all(str_replace_all(str_replace_all( cName, "\\(|\\)| ", "_"), "__", "_"), "_$", "") )

    plotlabels = NULL
    
    if (cName %in% names(restrict_labels))
    {
        cInfoElem = restrict_labels[[cName]]
      
        popName = cInfoElem$fname
        plotlabels = cInfoElem$genes
    }
    

    print(plotlabels)
    
    if (turnExpression)
    {
      indf$avg_log2FC = -indf$avg_log2FC
    }
    
    
    
    keyvals <- ifelse(indf$avg_log2FC > 0,
                      
                      ifelse(indf$p_val_adj < pCutoff & indf$avg_log2FC > FCcutoff, colors$pos$sig, colors$pos$nosig)
                      ,
                      ifelse(indf$avg_log2FC <= 0,
                             ifelse(indf$p_val_adj < pCutoff & indf$avg_log2FC < -FCcutoff, colors$neg$sig, colors$neg$nosig)
                             ,
                             'black'))
    
    
    keyvals[is.na(keyvals)] <- 'black'
    names(keyvals)[keyvals == '#F47B78'] <- 'UpReg sig'
    names(keyvals)[keyvals == '#FACAC9'] <- 'UpReg non-sig'
    names(keyvals)[keyvals == '#448CCA'] <- 'DownReg sig'
    names(keyvals)[keyvals == '#B4D1E9'] <- 'DownReg non-sig'
    
    txtLabelSize = 5
    axisLabelSize=12
    legendLabSize=12
    legendIconSize=5.0
  
    filename=paste(outname, popName, "png", sep=".")
    print(filename)
    png(filename=filename,width = 1200, height = 700)
    p=EnhancedVolcano(indf,
      lab = indf$gene,
      selectLab=plotlabels,
      x = 'avg_log2FC',
      y = 'p_val_adj',
      colCustom = keyvals,
      pCutoff = pCutoff,
      FCcutoff = FCcutoff,
      pointSize = 2.0,
      legendPosition = 'right',
      drawConnectors = TRUE,
      widthConnectors = 0.5,
      title = title,
      subtitle = subtitle,
      axisLabSize=axisLabelSize,
      labSize = txtLabelSize,
      legendLabSize = legendLabSize,
      legendIconSize = legendIconSize,
      ylab = bquote(~-Log[10]~italic(adj~P-Value)),
      gridlines.major = FALSE,
      gridlines.minor = FALSE
     )
    plot(p)
    dev.off()
    
    
        filename=paste(outname, popName, "pdf", sep=".")
    print(filename)
        pdf(filename,width = 12, height = 7)
    p=EnhancedVolcano(indf,
      lab = indf$gene,
      selectLab=plotlabels,
      x = 'avg_log2FC',
      y = 'p_val_adj',
      colCustom = keyvals,
      pCutoff = pCutoff,
      FCcutoff = FCcutoff,
      pointSize = 2.0,
      legendPosition = 'right',
      drawConnectors = TRUE,
      widthConnectors = 0.5,
      title = title,
      subtitle = subtitle,
      axisLabSize=axisLabelSize,
      labSize = txtLabelSize,
      legendLabSize = legendLabSize,
      legendIconSize = legendIconSize,
      ylab = bquote(~-Log[10]~italic(adj~P-Value)),
      gridlines.major = FALSE,
      gridlines.minor = FALSE
     )
    plot(p)
    dev.off()
    
    filename=paste(outname, popName, "svg", sep=".")
    print(filename)
    svglite(file = filename, width = 12, height = 7)
    p=EnhancedVolcano(indf,
      lab = indf$gene,
      selectLab=plotlabels,
      x = 'avg_log2FC',
      y = 'p_val_adj',
      colCustom = keyvals,
      pCutoff = pCutoff,
      FCcutoff = FCcutoff,
      pointSize = 2.0,
      legendPosition = 'right',
      drawConnectors = TRUE,
      widthConnectors = 0.5,
      title = title,
      subtitle = subtitle,
      axisLabSize=axisLabelSize,
      labSize = txtLabelSize,
      legendLabSize = legendLabSize,
      legendIconSize = legendIconSize,
      ylab = bquote(~-Log[10]~italic(adj~P-Value)),
      gridlines.major = FALSE,
      gridlines.minor = FALSE
     )
    plot(p)
    dev.off()
  
    plot(p)
    
  }
  
}

log_both <- function(x){ifelse(x == 0, 0, log(abs(x)) * sign(x))}

log_both_trans <- 
  function(){
    trans_new(name = 'log_both', 
              transform = log_both,
              inverse = log_both) #not clear what `inverse` does
}


cbind.fill <- function(...){
  nm <- list(...) 
  nm <- lapply(nm, as.matrix)
  n <- max(sapply(nm, nrow)) 
  do.call(cbind, lapply(nm, function (x) 
    rbind(x, matrix(, n-nrow(x), ncol(x))))) 
}

make_descr_label = function(plot, descr)
{
  descrLabel <- ggdraw() + draw_label(descr, fontface='bold', angle = 0)
  
  pe = combine_plot_grid_list(plotlist=list("a"=descrLabel, "b"=plot), ncol=1, nrow=2, labels=NULL,rel_heights = c(0.1, 1), align = "h", axis = "l")
  
  return(pe)
}

makeSideBySideDotPlot = function(scobj, plotElems, featureGenes = c(""), group.by="cellnames_manual", col.min = -3, col.max = 3, cols = c("grey", "blue"), title="", scaled=T, scale.by="GROUP", rotate.x=F, abundance.perelem=FALSE)
{

  stopifnot(scale.by %in% c("GROUP", "FEATURE", "ALL"))

  plot_list = list()
  plot_orig_list = list()
  allFoundFeatures = c()
  
  ctFractions = list()
  
  

  scFullTable = table(scobj[[group.by]])
  scFullDf = as.data.frame(scFullTable)

  fillLimitsMax = 0
  elemCount = 0
  for (plotName in names(plotElems))
  {
    print(plotName)
    plotData = plotElems[[plotName]]
    plotCells = plotData$cells
    plotDescr = plotData$label
    
    scobj_subset = subset(scobj, cells=plotCells)
    #print(scobj_subset)
    plotElem_orig = DotPlot(scobj_subset, features=featureGenes, group.by = group.by, col.min = col.min, col.max=col.max, cols = cols)
    #print(plotElem_orig$data)
    #print(plotElem_orig$data$avg.exp)
    scTable = table(scobj_subset[[group.by]])
    scDf = as.data.frame(scTable)

    
    if (abundance.perelem)
    {
      scDf$perc = scDf$Freq / length(plotCells)
    } else {
      scDf$perc = scDf$Freq / sum(scFullDf$Freq) # total abundance!
    }
    #print(scDf)

    fillLimitsMax = max(c(fillLimitsMax,max(scDf$perc)))

    ctFractions[[plotName]] = scDf
    #print(ctFractions)
    
    scobj_subset = NULL
    
    allFoundFeatures = unique(c(allFoundFeatures, as.character(plotElem_orig$data$features.plot)))
    
    plotElem_orig = plotElem_orig + scale_color_gradient(limits=c(col.min, col.max), low = cols[1], high = cols[2])
    plotElem_orig = plotElem_orig + theme(axis.title.x=element_blank(), axis.title.y=element_blank(), legend.position = "none")
    
    plot_orig_list[[plotName]] = plotElem_orig
  }  
  
  if (T)
  {
    
    # initialize Values
    for (plotName in names(plot_orig_list))
    {
        plot_orig_list[[plotName]]$data$avg.exp.scaled2 = 0
        
        plot_orig_list[[plotName]]$data$id = as.character(plot_orig_list[[plotName]]$data$id)
        plot_orig_list[[plotName]]$data$features.plot = as.character(plot_orig_list[[plotName]]$data$features.plot)
        
    }
    
    allIDs = c()
    allFeatures = c()
    for (plotName in names(plot_orig_list))
    {
        allIDs = c(allIDs, plot_orig_list[[plotName]]$data$id)
        allFeatures = c(allFeatures, plot_orig_list[[plotName]]$data$features.plot)
    }
    allIDs = unique(allIDs)
    allFeatures = unique(allFeatures)
    
    #print(allIDs)
    #print(allFeatures)
    

    if (scale.by == "FEATURE")
    {
          # calculate avg.exp.scaled2 for each feature
      for (featureName in allFoundFeatures)
      {
      
        allUnscaledValues = NULL
        for (plotName in names(plot_orig_list))
        {
          pData = plot_orig_list[[plotName]]$data
          

          missingCelltypes = setdiff(allIDs, unique(pData$id))
        
          
          for (celltype in missingCelltypes)
          {
            for (feature in allFeatures)
            {
              dfrow = data.frame(avg.exp=0.0,pct.exp=0.0,features.plot=feature, id=celltype, avg.exp.scaled=0.0, avg.exp.scaled2=0.0)
              pData = rbind.data.frame(pData, dfrow, stringsAsFactors = F)
              
            }
          }
          
          pData$id = factor(pData$id, levels = allIDs)
          pData$features.plot = factor(pData$features.plot, levels=allFeatures)
          

          plot_orig_list[[plotName]]$data = pData
          
          
          if (is.null(allUnscaledValues))
          {
            allUnscaledValues = data.frame(plotName=pData[ pData$features.plot==featureName, ]$avg.exp)
            allUnscaledValues[[plotName]] = pData[ pData$features.plot==featureName, ]$avg.exp  
            allUnscaledValues[["plotName"]] = NULL
            
          } else {
            allUnscaledValues[[plotName]] = pData[ pData$features.plot==featureName, ]$avg.exp  
          }
          
        }
        allUnscaledValues$rnames = as.numeric(rownames(allUnscaledValues))
        
        allUnscaledLong = allUnscaledValues %>% gather(Type, Value, names(plot_orig_list))
        allUnscaledLong$Value = scale(allUnscaledLong$Value)
        
        allScaledValues = allUnscaledLong %>% spread(Type, Value) %>% arrange( order(rnames))
        
        for (plotName in names(plot_orig_list))
        {

          plotElem_orig = plot_orig_list[[plotName]]
          pData = plotElem_orig$data
    
          
          # https://github.com/satijalab/seurat/issues/2798
          plotElem_orig$layers[[1]] <- NULL # remove original geom_point layer where the color scale is hard-coded to use scaled average expression
          
          origData = plotElem_orig$data[plotElem_orig$data$features.plot==featureName, ]

          if (scaled)
          {
            plotElem_orig$data[plotElem_orig$data$features.plot==featureName, "avg.exp.scaled2"] = allScaledValues[,plotName]
          } else {
            plotElem_orig$data[plotElem_orig$data$features.plot==featureName, "avg.exp.scaled2"] = plotElem_orig$data[plotElem_orig$data$features.plot==featureName, "avg.exp"]
          }
          

          plotElem_orig$data$idn=as.numeric(plotElem_orig$data$id)


          plot_orig_list[[plotName]] = plotElem_orig
        }


        
      }
    } else if ((scale.by == "ALL") && (scaled == T)) {

        combinedDataDF = data.frame()

        for (plotName in names(plot_orig_list))
        {

          plotElem_orig = plot_orig_list[[plotName]]
                   
          # https://github.com/satijalab/seurat/issues/2798
          plotElem_orig$layers[[1]] <- NULL # remove original geom_point layer where the color scale is hard-coded to use scaled average expression



          storeDF = plotElem_orig$data
          storeDF["plotpart"] = plotName

          combinedDataDF = rbind(combinedDataDF, storeDF)
        }


        combinedDataDF[, "avg.exp.scaled2"] = scale(combinedDataDF$avg.exp)
        combinedDataDF$id = factor(combinedDataDF$id, levels = mixedsort(as.character(unique(combinedDataDF$id))))


        for (plotName in names(plot_orig_list))
        {

          plotData = combinedDataDF[combinedDataDF$plotpart == plotName,]
          plotElem_orig = plot_orig_list[[plotName]]
          plotElem_orig$data$avg.exp.scaled2 = plotData$avg.exp.scaled2

          plotElem_orig$data$id = plotData$id
          plotElem_orig$data$idn= as.numeric(plotData$id)
          plotElem_orig$data$features.plot = factor(plotElem_orig$data$features.plot, levels=unique(plotElem_orig$data$features.plot))
          
          #print(plotName)
          #print(plotElem_orig$data$id)

          plot_orig_list[[plotName]] = plotElem_orig

        }
       
    } else if ((scale.by == "GROUP") || ((scale.by == "ALL") && (scaled == F))) {

        for (plotName in names(plot_orig_list))
        {

          plotElem_orig = plot_orig_list[[plotName]]
                   
          # https://github.com/satijalab/seurat/issues/2798
          plotElem_orig$layers[[1]] <- NULL # remove original geom_point layer where the color scale is hard-coded to use scaled average expression


          plotElem_orig$data$id = factor(plotElem_orig$data$id, levels = mixedsort(as.character(unique(plotElem_orig$data$id))))
          plotElem_orig$data$idn= as.numeric(plotElem_orig$data$id)
          plotElem_orig$data$features.plot = factor(plotElem_orig$data$features.plot, levels=unique(plotElem_orig$data$features.plot))
          

          if (scaled)
          {
            plotElem_orig$data[, "avg.exp.scaled2"] = scale(plotElem_orig$data$avg.exp)
          } else {
            plotElem_orig$data[, "avg.exp.scaled2"] = plotElem_orig$data$avg.exp
          }

          plot_orig_list[[plotName]] = plotElem_orig
        }
    } else {
      stopifnot(FALSE)
    }


    idDF = data.frame()

    for (plotName in names(plot_orig_list))
    {
      plotElem_orig = plot_orig_list[[plotName]]
      idDF = rbind(idDF, plotElem_orig$data[, c("id", "idn")])
    }

    #print("idDF")
    rownames(idDF) = NULL
    idDF = unique(idDF)

    
    for (plotName in names(plot_orig_list))
    {
      plotElem_orig = plot_orig_list[[plotName]]
      pData = plotElem_orig$data
      
      pData2 = merge(x=pData,y=ctFractions[[plotName]],by.x="id", by.y="Var1",all.x=TRUE)
      plotElem_orig$data <-pData2 %>% mutate(featuren=as.numeric(features.plot), percn=100*perc)
      

      
      plotElem_orig$data[plotElem_orig$data$avg.exp.scaled2>col.max, 'avg.exp.scaled2'] = col.max
      plotElem_orig$data[plotElem_orig$data$avg.exp.scaled2<col.min, 'avg.exp.scaled2'] = col.min
      
      pData = plotElem_orig$data

      if (scaled==T)
      {
        exprTitle = paste("Avg. Expression (scaled by ", scale.by, ")", sep="")
      } else {
        exprTitle = "Avg. Expression"
      }

      #print(plotName)
      #print(pData)
      #print(idDF)
      #print(fillLimitsMax)


      fillLimits = c(0, ceiling(fillLimitsMax*10)*10)
      
      plotElem_orig <- ggplot(pData) +
        scale_x_continuous(breaks=plotElem_orig$data$featuren, labels=plotElem_orig$data$features.plot) +
        scale_y_continuous(breaks=idDF$idn, labels=idDF$id, limits = c(min(idDF$idn)-0.6, max(idDF$idn)+0.6))+
        geom_rect(aes(xmin=featuren-.5, xmax=featuren+.5, ymin = idn-0.5, ymax = idn+0.5, fill=percn), alpha = 0.4, linetype="blank") +
        scale_fill_distiller(palette='Spectral', limits = fillLimits)+
        scale_size_continuous(range = c(0, 10))+
        geom_point(aes(x=featuren, y=idn, colour = avg.exp.scaled2, size = pct.exp)) +
        scale_color_gradient2(limits=c(col.min, col.max), low = cols[1], mid=cols[2], high = cols[3])+
        guides(color=guide_colourbar(title=exprTitle, order = 1),
        size=guide_legend(title="Percent Expressing", order = 2),
        fill=guide_colourbar(title="Cell Abundance (selected cells)", order = 3))
      
      #plotElem_orig = plotElem_orig + scale_color_gradient(limits=c(col.min, col.max), low = cols[1], high = cols[2])
      plotElem_orig = plotElem_orig + theme(axis.title.x=element_blank(), axis.title.y=element_blank(), legend.position = "none", panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
panel.background = element_blank(), axis.line = element_line())

      if (rotate.x)
      {
        plotElem_orig = plotElem_orig + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
      }
              #scale_color_gradient()+

      plot_orig_list[[plotName]] = plotElem_orig
    }
  }
  
    print("Final plot checks")

  
  elemCount = 0
  plot_list=list()
  for (plotName in names(plot_orig_list))
  {
    
    plotElem_orig = plot_orig_list[[plotName]]
    
    elemCount = elemCount + 1
    
    plotElem = plotElem_orig
    
    #+  theme(
              #axis.text.y = element_blank(),
    #          panel.grid.major.x = element_line( size=.1, color="black" ),
              #axis.line.y = element_line(size = 0),
              #axis.ticks.length.y = unit(0, "points")
    #)
  
    print("descr label")
    pe = make_descr_label(plotElem, plotElems[[plotName]]$label)
  
    plot_list[[plotName]] = pe
  
  }


  legendDescr = 'Average Expression'
  if (scaled)
  {
    legendDescr = 'Average Scaled Expression'  
  }

  print("Preparing Legend")
  
  # extract a legend that is laid out horizontally #
  legend_b <- get_legend(
    plotElem_orig + 
      guides(color = guide_colorbar(title = legendDescr, direction="horizontal"))+ #guide_legend(nrow = 1, override.aes = list(size=6)))+
      theme(legend.position = "bottom")
  )
  
  title <- ggdraw() + draw_label(title, fontface='bold')
  print("Preparing plot")

  dataList = list()
  for (i in 1:length(plot_list))
  {
    dataList[[i]] = plot_list[[i]]$data
  }
  ap=cowplot::plot_grid(
    plotlist = plot_list,
    labels = NULL,
    nrow=1,
    align = "v", axis="bt"
  )

  fp = cowplot::plot_grid(title, ap, legend_b, ncol = 1, rel_heights = c(.05, 1, 0.1) )
  fp$data = dataList

  return(fp)
}


#Cellcycle Scoring 
exp.mat <- read.table(file = "/usr/local/hdd2/data/alejandro/nestorawa_forcellcycle_expressionMatrix.txt", header = TRUE,
    as.is = TRUE, row.names = 1)
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes










makeTopDownDotplot = function(scobj, plotElems, featureGenes = c(""), group.by="cellnames_manual", col.min = -3, col.max = 3, cols = c("grey", "blue"), title="", scaled=T, scale.by="GROUP")
{

  stopifnot(scale.by %in% c("GROUP", "ALL"))

  plot_list = list()
  plot_orig_list = list()
  allFoundFeatures = c()
  
  ctFractions = data.frame()

  scFullTable = table(scobj[[group.by]])
  scFullDf = as.data.frame(scFullTable)

  fillLimitsMax = 0
  elemCount = 0
  for (plotName in names(plotElems))
  {
    print(plotName)
    plotData = plotElems[[plotName]]
    plotCells = plotData$cells
    plotDescr = plotData$label
    
    scobj_subset = subset(scobj, cells=plotCells)
    plotElem_orig = DotPlot(scobj_subset, features=featureGenes, group.by = group.by, col.min = col.min, col.max=col.max, cols = cols)
    scTable = table(scobj_subset[[group.by]])
    scDf = as.data.frame(scTable)
    
    scDf$perc = scDf$Freq / sum(scFullDf$Freq) # total abundance!
    scDf$plotpart = plotName

    fillLimitsMax = max(c(fillLimitsMax,max(scDf$perc)))
    ctFractions = rbind(ctFractions, scDf)
    
    scobj_subset = NULL
    
    allFoundFeatures = unique(c(allFoundFeatures, as.character(plotElem_orig$data$features.plot)))
    
    plotElem_orig = plotElem_orig + scale_color_gradient(limits=c(col.min, col.max), low = cols[1], high = cols[2])
    plotElem_orig = plotElem_orig + theme(axis.title.x=element_blank(), axis.title.y=element_blank(), legend.position = "none")
    
    plot_orig_list[[plotName]] = plotElem_orig
  }
  
  
    
    allIDs = c()
    allFeatures = c()
    for (plotName in names(plot_orig_list))
    {
        plot_orig_list[[plotName]]$data$id = as.character(plot_orig_list[[plotName]]$data$id)
        plot_orig_list[[plotName]]$data$features.plot = as.character(plot_orig_list[[plotName]]$data$features.plot)
    

        allIDs = c(allIDs, plot_orig_list[[plotName]]$data$id)
        allFeatures = c(allFeatures, plot_orig_list[[plotName]]$data$features.plot)
    }
    allIDs = unique(allIDs)
    allFeatures = unique(allFeatures)
    
    print(allIDs)
    print(allFeatures)
    
    if ((scale.by == "ALL") && (scaled == T)) {

        combinedDataDF = data.frame()

        for (plotName in names(plot_orig_list))
        {

          plotElem_orig = plot_orig_list[[plotName]]
                   
          # https://github.com/satijalab/seurat/issues/2798
          plotElem_orig$layers[[1]] <- NULL # remove original geom_point layer where the color scale is hard-coded to use scaled average expression



          storeDF = plotElem_orig$data
          storeDF["plotpart"] = plotName

          storeDF$id = factor(storeDF$id, levels = mixedsort(as.character(unique(storeDF$id))))
          storeDF$idn= as.numeric(storeDF$id)

          combinedDataDF = rbind(combinedDataDF, storeDF)
        }

        combinedDataDF[, "avg.exp.scaled2"] = scale(combinedDataDF$avg.exp)
        combinedDataDF$id = factor(combinedDataDF$id, levels = mixedsort(as.character(unique(combinedDataDF$id))))
        combinedDataDF$idn= as.numeric(combinedDataDF$id)
        combinedDataDF$features.plot = factor(combinedDataDF$features.plot, levels=unique(combinedDataDF$features.plot))

       
    } else if ((scale.by == "GROUP") || ((scale.by == "ALL") && (scaled == F))) {
        combinedDataDF = data.frame()

        for (plotName in names(plot_orig_list))
        {

          plotElem_orig = plot_orig_list[[plotName]]
                   
          # https://github.com/satijalab/seurat/issues/2798
          plotElem_orig$layers[[1]] <- NULL # remove original geom_point layer where the color scale is hard-coded to use scaled average expression

          plotElem_orig$data["plotpart"] = plotName

          plotElem_orig$data$id = factor(plotElem_orig$data$id, levels = mixedsort(as.character(unique(plotElem_orig$data$id))))
          plotElem_orig$data$idn= as.numeric(plotElem_orig$data$id)
          plotElem_orig$data$features.plot = factor(plotElem_orig$data$features.plot, levels=unique(plotElem_orig$data$features.plot))
          

          if (scaled)
          {
            plotElem_orig$data[, "avg.exp.scaled2"] = scale(plotElem_orig$data$avg.exp)
          } else {
            plotElem_orig$data[, "avg.exp.scaled2"] = plotElem_orig$data$avg.exp
          }

          combinedDataDF = rbind(combinedDataDF, plotElem_orig)
        }
    } else {
      stopifnot(FALSE)
    }

    print(combinedDataDF)




    pData2 = merge(x=combinedDataDF,y=ctFractions,by.x=c("id", "plotpart"), by.y=c("Var1", "plotpart"),all.x=TRUE)
    combinedDataDF <-pData2 %>% mutate(featuren=as.numeric(features.plot), percn=100*perc)
    
    combinedDataDF[combinedDataDF$avg.exp.scaled2>col.max, 'avg.exp.scaled2'] = col.max
    combinedDataDF[combinedDataDF$avg.exp.scaled2<col.min, 'avg.exp.scaled2'] = col.min


    if (scaled==T)
    {
      exprTitle = paste("Avg. Expression (scaled by ", scale.by, ")", sep="")
    } else {
      exprTitle = "Avg. Expression"
    }

    combinedDataDF$id = paste(combinedDataDF$id, combinedDataDF$plotpart, sep=" ")

    combinedDataDF$id = factor(combinedDataDF$id, levels = mixedsort(as.character(unique(combinedDataDF$id))))
    combinedDataDF$idn= as.numeric(combinedDataDF$id)

    idDF = combinedDataDF[, c("id", "idn")]
    print("idDF")
    rownames(idDF) = NULL
    idDF = unique(idDF)

    print(plotName)
    print(idDF)
    print(fillLimitsMax)
    print(combinedDataDF)







    fillLimits = c(0, ceiling(fillLimitsMax*10)*10)
        
  mainPlot <- ggplot(combinedDataDF) +
    scale_x_continuous(breaks=combinedDataDF$featuren, labels=combinedDataDF$features.plot) +
    scale_y_continuous(breaks=idDF$idn, labels=idDF$id, limits = c(min(idDF$idn)-0.6, max(idDF$idn)+0.6))+
    geom_rect(aes(xmin=featuren-.5, xmax=featuren+.5, ymin = idn-0.5, ymax = idn+0.5, fill=percn), alpha = 0.4, linetype="blank") +
    scale_fill_distiller(palette='Spectral', limits = fillLimits)+
    scale_size_continuous(range = c(0, 10))+
    geom_point(aes(x=featuren, y=idn, colour = avg.exp.scaled2, size = pct.exp)) +
    scale_color_gradient2(limits=c(col.min, col.max), low = cols[1], mid=cols[2], high = cols[3])+
    guides(color=guide_colourbar(title=exprTitle, order = 1),
    size=guide_legend(title="Percent Expressing", order = 2),
    fill=guide_colourbar(title="Cell Abundance (selected cells)", order = 3))
  
  mainPlot = mainPlot + theme(axis.title.x=element_blank(), axis.title.y=element_blank(), legend.position = "none", panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
  panel.background = element_blank(), axis.line = element_line())
  
  legendDescr = 'Average Expression'
  if (scaled)
  {
    legendDescr = 'Average Scaled Expression'  
  }
  
  # extract a legend that is laid out horizontally #
  legend_b <- get_legend(
    mainPlot + 
      guides(color = guide_colorbar(title = legendDescr, direction="horizontal"))+ #guide_legend(nrow = 1, override.aes = list(size=6)))+
      theme(legend.position = "bottom")
  )
  
  title <- ggdraw() + draw_label(title, fontface='bold')

  print("Combining plots")
  fp = combine_plot_grid_list(plotlist=list(title, mainPlot, legend_b), ncol = 1, rel_heights = c(.05, 1, 0.1) )


  return(fp)
}










splitFeaturePlot = function(obj, feature, split.by, title, filename, limits = c(-2, 2), low = "#713fdf", high = "#df713f", mid = "grey")
{
  pds = FeaturePlot(obj, features = feature, reduction = "umap", split.by=split.by, combine=F,min.cutoff=0, max.cutoff=7,order=T)
  pds[[1]] = pds[[1]] + ggtitle(NULL)

  print(paste(length(pds)))
  pdsRange = c(1:(length(pds)-1))

  for (i in pdsRange)
  {
    print(i)
    pds[[i]] = pds[[i]] + scale_color_gradient(limits = c(0,7), low = "lightgrey", high = "blue", guide=FALSE)+ theme(axis.text.x = element_text(face="bold", color="#000000", size=14, angle=0), axis.text.y = element_text(face="bold", color="#000000", size=14, angle=0))+labs(x="", y="")  
  }

  #pds[[1]] = pds[[1]] + scale_color_gradient(limits = c(0,7), low = "lightgrey", high = "blue", guide=FALSE) + theme(axis.text.x = element_text(face="bold", color="#000000", size=14, angle=0), axis.text.y = element_text(face="bold", color="#000000", size=14, angle=0))+labs(x="", y="")
  #pds[[2]] = pds[[2]] + scale_color_gradient(limits = c(0,7), low = "lightgrey", high = "blue", guide=FALSE)+ theme(axis.text.x = element_text(face="bold", color="#000000", size=14, angle=0), axis.text.y = element_text(face="bold", color="#000000", size=14, angle=0))+labs(x="", y="")
  #pds[[3]] = pds[[3]] + scale_color_gradient(limits = c(0,7), low = "lightgrey", high = "blue", guide=FALSE)+ theme(axis.text.x = element_text(face="bold", color="#000000", size=14, angle=0), axis.text.y = element_text(face="bold", color="#000000", size=14, angle=0))+labs(x="", y="")
  pds[[length(pds)]] = pds[[length(pds)]] + scale_color_gradient(limits = c(0,7), low = "lightgrey", high = "blue"             ) + theme(axis.text.x = element_text(face="bold", color="#000000", size=14, angle=0), axis.text.y = element_text(face="bold", color="#000000", size=14, angle=0))+labs(x="", y="")

  prow =combine_plot_grid_list(plotlist=pds, label_x = "a", ncol=length(pds), align="hv")
  # now add the title
  title <- ggdraw() + 
    draw_label(
      title,
      fontface = 'bold',
      x = 0,
      hjust = 0
    ) +
    theme(
      # add margin on the left of the drawing canvas,
      # so title is aligned with left edge of first plot
      plot.margin = margin(0, 0, 0, 7)
    )

  # add the legend to the row we made earlier. Give it one-third of 
  # the width of one plot (via rel_widths).
  fplot = plot_grid(
    title, prow,
    ncol = 1,
    # rel_heights values control vertical title margins
    rel_heights = c(0.1, 1)
  )

  prow = combine_plot_grid_list(plotlist=pds, label_x = "a", ncol=length(pds), align="hv")
  fplot= combine_plot_grid_list(plotlist=list(title, prow), ncol=1, rel_heights=c(0.1, 1))

  if (!is.null(filename))
  {
    save_plot(fplot, filename, fig.width=length(pds) *6 + 2, fig.height=6)

  }
  return(fplot)
}


VlnBoxPlot = function( obj.sc, gene, group.by="idents", split.by=NULL, pt.size=0, assay="RNA", min.threshold=3, col=NULL, per.sample=FALSE)
{

    remainingCells = NULL

    if (is.null(split.by))
    {

        cellGroupDF = obj.sc[[ group.by ]]

        countDF = table(cellGroupDF)
        countDF = countDF[countDF >= min.threshold]

        remainingClusters = rownames(countDF)
        remainingCells = rownames(cellGroupDF)[ cellGroupDF[[group.by]] %in% remainingClusters ]

    } else {

        cellGroupDF = obj.sc[[ c(group.by, split.by) ]]

        countDF = as.data.frame(table(cellGroupDF))
        countDF = countDF[countDF$Freq >= min.threshold,]

        if (per.sample)
        {

          snames = unique(obj.sc[[ split.by ]][[split.by]])
          remainingCells = c()

          for (sname in snames)
          {
            snameDF = countDF[ countDF[split.by] == sname,]
            snameDF = snameDF[snameDF$Freq >= min.threshold,]

            sCells = rownames(obj.sc[[split.by]])[ obj.sc[[split.by]][[split.by]] == sname ]
            gCells = rownames(obj.sc[[group.by]])[ obj.sc[[group.by]][[group.by]] %in% snameDF[[group.by]]]

            remainingCells = c(remainingCells, intersect(sCells, gCells))
          }

        } else {
          remainingClusters = countDF[duplicated(countDF[group.by]),c(group.by)]

          cellGroupDF = obj.sc[[ group.by ]]
          remainingCells = rownames(cellGroupDF)[ cellGroupDF[[group.by]] %in% remainingClusters ]
        }



    }
    

    print(head(remainingCells))
    obj.sc.subset = subset(obj.sc, cells=remainingCells)
    print(obj.sc.subset)


    if (is.null(split.by))
    {

      ucols = NULL
      if (!is.null(col))
      {
        ncolors = length(unique(obj.sc[[group.by]][[group.by]]))
        ucols = rep(col, ncolors)
      }
      

    plots=VlnPlot(obj.sc.subset, gene, group.by=group.by, split.by=split.by, pt.size=pt.size, assay=assay, col=ucols)
    plots = plots + geom_boxplot(color="grey", alpha=0.4) + stat_summary(fun=mean, geom="point", color="black", size=4)

    } else {

    plots=VlnPlot(obj.sc.subset, gene, group.by=group.by, split.by=split.by, pt.size=pt.size, assay=assay, col=col)
    plots = plots + geom_boxplot(color="grey", alpha=0.4, position =position_dodge(width = 0.9)) + stat_summary(fun=mean, geom="point", aes(group=split), position=position_dodge(.9), color="black", size=4)

    }
    

    return(plots)
}


SplitVlnBoxPlot = function( obj.sc, gene, group.by="idents", split.by=NULL, pt.size=0, assay="RNA", min.threshold=3, col=NULL, per.sample=TRUE)
{

    splitValues = unique(obj.sc[[split.by]][[split.by]])

    vplots = list()

    for (sname in splitValues)
    {
      print(sname)
      subsetcellnames = rownames(obj.sc[[split.by]])[ obj.sc[[split.by]][split.by] == sname ]

      if (is.null(col))
      {
        scol = NULL
      } else {
        scol = col[[sname]]
      }

      p=VlnBoxPlot( subset(obj.sc, cells=subsetcellnames), gene, group.by=group.by, split.by=NULL, pt.size=pt.size, assay=assay, min.threshold=min.threshold, col=scol, per.sample=per.sample)
      p = p + ggtitle(paste(gene, "-", sname))
      vplots[[sname]] = p
    }

    ps = combine_plot_grid_list(plotlist=vplots, nrow=1)

    return(ps)
}











######Alejandro's adaptations


compareSMCs = function(inobj, cellsID1, cellsID2, suffix1, suffix2, group.by=NULL, assay="RNA", test="t", outfolder="./", fcCutoff=0.25, all=FALSE)
{
  if (is.null(group.by))
  {
    clusterIDs = as.character(sort(unique(Idents(inobj))))  
  } else {
    clusterIDs = as.character(sort(unique(inobj[[group.by]][,])))
  }
  
  retList = list()
  
  for(clusterID in clusterIDs)
  {
  
      cellIdents = Idents(inobj)
      
      if (is.null(group.by))
      {
        cellIdents.c = names(cellIdents[cellIdents == clusterID])
  
      } else {
        cellIdents.c = colnames(inobj[,inobj[[group.by]] == clusterID])
      }
      
      cellIdents.c = unlist(lapply(cellIdents.c, as.character))
  
      print(paste("Processing cluster", clusterID, "with a total of", length(cellIdents.c), "cells"))
      
      cells1 = intersect(cellsID1, cellIdents.c)
      cells2 = cellsID2 #change this to select all media SMCs
      
      print(c(length(cells1), length(cells2)))
      
      if (length(cells1) < 3)
      {
        print("Cells1 too few")
        next
      }
      
      if (length(cells2) < 3)
      {
        print("Cells2 too few")
        next
      }

      clusterID_file=str_replace_all(str_replace_all(clusterID, "\\/", "_"), " ", "_")

  
      deMarkers = compareClusters(scdata=inobj,
                                  cellsID1=cells1,
                                  cellsID2=cells2,
                                  prefix= paste("cluster", clusterID_file, sep="_"),
                                  suffix1=suffix1,#paste(suffix1, clusterID, sep="_"),
                                  suffix2=suffix2, #paste(suffix2, clusterID, sep="_"),
                                  test=test, fcCutoff=fcCutoff, assay=assay, outfolder=outfolder)
  
  
      retList[[clusterID]] = deMarkers
  
  }
  
  return(retList)
}





compareCellsByOrigin = function(inobj, cellsID1, cellsID2, suffix1, suffix2, group.by=NULL, assay="RNA", test="t", outfolder="./", fcCutoff=0.25, all=FALSE)
{
  if (is.null(group.by))
  {
    clusterIDs = "All_SMCs"
  } else {
    clusterIDs = as.character(sort(unique(inobj[[group.by]][,])))
  }
  
  retList = list()
  
  for(clusterID in clusterIDs)
  {
  
      cellIdents = inobj$orig_project
      
      if (is.null(group.by))
      {
        cellIdents.c = names(cellIdents)
  
      } else {
        cellIdents.c = colnames(inobj[,inobj[[group.by]] == clusterID])
      }
      
      cellIdents.c = unlist(lapply(cellIdents.c, as.character))
  
      print(paste("Processing cluster", clusterID, "with a total of", length(cellIdents.c), "cells"))
      
      cells1 = intersect(cellsID1, cellIdents.c)
      cells2 = intersect(cellsID2, cellIdents.c)
      
      print(c(length(cells1), length(cells2)))
      
      if (length(cells1) < 3)
      {
        print("Cells1 too few")
        next
      }
      
      if (length(cells2) < 3)
      {
        print("Cells2 too few")
        next
      }

      clusterID_file=str_replace_all(str_replace_all(clusterID, "\\/", "_"), " ", "_")

  
      deMarkers = compareClusters(scdata=inobj,
                                  cellsID1=cells1,
                                  cellsID2=cells2,
                                  prefix= NULL,
                                  suffix1=suffix1,#paste(suffix1, clusterID, sep="_"),
                                  suffix2=suffix2, #paste(suffix2, clusterID, sep="_"),
                                  test=test, fcCutoff=fcCutoff, assay=assay, outfolder=outfolder)
  
  
      retList[[clusterID]] = deMarkers
  
  }
  
  return(retList)
}




### Basic function to convert mouse to human gene names
convertMouseGeneList <- function(x){

require("biomaRt")
human = useMart("ensembl", dataset = "hsapiens_gene_ensembl", host = "https://dec2021.archive.ensembl.org/")
mouse = useMart("ensembl", dataset = "mmusculus_gene_ensembl", host = "https://dec2021.archive.ensembl.org/")

genesV2 = getLDS(attributes = c("hgnc_symbol"), filters = "hgnc_symbol", values = x , mart = human, attributesL = c("mgi_symbol"), martL = mouse, uniqueRows=T)
humanx <- unique(genesV2[, 2])

### Print the first 6 genes found to the screen
print(head(humanx))
return(humanx)
}

convertHumanGeneList <- function(x){

require("biomaRt")
human = useMart("ensembl", dataset = "hsapiens_gene_ensembl", host = "https://dec2021.archive.ensembl.org/")
mouse = useMart("ensembl", dataset = "mmusculus_gene_ensembl", host = "https://dec2021.archive.ensembl.org/")

genesV2 = getLDS(attributes = c("mgi_symbol"), filters = "mgi_symbol", values = x , mart = mouse, attributesL = c("hgnc_symbol"), martL = human, uniqueRows=T)
humanx <- unique(genesV2[, 2])

### Print the first 6 genes found to the screen
print(head(humanx))
return(humanx)
}