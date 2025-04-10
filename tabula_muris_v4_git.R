


# turn off diagnostics in RStudio for this file to allow Matrix package to be loaded

# !diagnostics off 


library(Matrix, lib.loc = "/nas/longleaf/home/wrobel/R/local_R_libraries/R-4.2.2")
library(irlba, lib.loc = "/nas/longleaf/home/wrobel/R/local_R_libraries/R-4.2.2")


library(ggplot2)

library(rasterpdf, lib.loc = "/nas/longleaf/home/wrobel/R/local_R_libraries/R-4.2.2")


library(scales)


# for md5sum
library(tools)


library(Seurat, lib.loc = "/proj/tinglab/users/wrobel/R/local_R_libraries/R-4.2.2")



script_start_time <- Sys.time()






# *** START user input ***

split_report <- TRUE




# sample_name <- "TM_facs_Large_Intestine"
# dataset_name <- "TM Large Intestine (facs)"
# plot_title <- "Murine Large Intestine"
# organism_name <- "mouse"
# tissue_name <- "Large Intestine (facs)"
# Robj_file <- "/work/users/w/r/wrobel/tinglab/tabula_muris/Jan2025_download/Robj/facs/facs_Large_Intestine_seurat_tiss.Robj"


sample_name <- "TM_facs_Lung"
dataset_name <- "TM Lung (facs)"
plot_title <- "Murine Lung"
organism_name <- "mouse"
tissue_name <- "Lung (facs)"
Robj_file <- "/work/users/w/r/wrobel/tinglab/tabula_muris/Jan2025_download/Robj/facs/facs_Lung_seurat_tiss.Robj"




gene_type <- "Kat-genes"

gene_list_file <- "/proj/tinglab/users/wrobel/NLR_gene_list/Nlr_gene_list_mouse.txt"


additional_genes <- c(
  
  "Actb",
  # "Myd88",
  # "Tlr2",
  # "Tlr6",
  # "Gsdma",
  # "Gsdmb",
  # "Gsdmc",
  # "Gsdmd",
  # "Gsdme"
  "Il1b",
  "Il18",
  "Ace2",
  
  "Casp1", 
  "Casp4", # Casp11 (synonym) 
  "Tmprss2"
)



plot_genes <- c(
  
  "Aim2",
  "Nlrp3",
  "Nlrp6",
  
  "Il1b",
  "Il18",
  "Ace2",
  
  "Casp1", 
  "Casp4", # Casp11 (synonym) 
  "Tmprss2"
)




plot_genes_Kat <- c(
  
  "Il18",
  "Il1b",
  "Casp4", # Casp11 (synonym)
  "Casp1", 
  "Nlrp6",
  "Nlrp3",
  "Aim2",
  "Tmprss2",
  "Ace2"
)


plot_genes <- plot_genes_Kat


vln_genes <- c(
  
  "Aim2",
  "Nlrp3",
  "Nlrp6"
)





additional_genes <- c(plot_genes, "Actb")




mgi_file <- "/proj/tinglab/users/wrobel/R/common_input/mgi/MGI_EntrezGene_2023-02-10.rpt"
# mgi_file <- "/proj/tinglab/users/wrobel/gene_datasets/mgi/Jan2025/MGI_EntrezGene_downloaded_2025-01-15.rpt"



# *** END user input ***





# direct output to a file
sink(file = paste0(sample_name, "_report.txt"), append = FALSE, split = split_report)



script_start_time_display <- format(script_start_time, "%Y-%m-%d %H:%M:%S %Z")

cat("R script started at:", "\n")
cat(script_start_time_display, "\n")

cat("\n")

session_info <- sessionInfo()

cat(R.version.string, " -- \"", R.version$nickname, "\"", sep = "")
cat("\n")

cat("Platform:", session_info$platform, "\n")

os_type <- .Platform$OS.type
cat("OS type:", os_type, "\n")

cat("Running under:", session_info$running, "\n")

sys_info <- Sys.info()
cat("sysname:", sys_info["sysname"], sys_info["release"], "\n")
#cat("sysname:", sys_info["sysname"], "\n")
#cat("release:", sys_info["release"], "\n")
#cat("version:", sys_info["version"], "\n")
#cat("machine:", sys_info["machine"], "\n")
cat("Hostname:", sys_info["nodename"], "\n")
cat("User:", sys_info["user"], "\n")

cat("\n")
cat("\n")








Robj_input_filename.md5 <- md5sum(Robj_file)

cat("Robj input file:", Robj_file, "\n")
cat("md5:", Robj_input_filename.md5, "\n")
cat("\n")
cat("\n")





load(Robj_file, verbose = TRUE)

cat("\n")
print(tiss)
cat("\n")

cat("Seurat object version:", as.character(tiss@version), "\n")


cat("\n")
cat("updating seurat object . . .", "\n")
cat("\n")

pbmc <- Seurat::UpdateSeuratObject(tiss)
rm(tiss)

cat("Updated seurat object version:", as.character(pbmc@version), "\n\n")



original_metadata <- pbmc@meta.data


write.table(original_metadata, file = paste0(sample_name, "_original_metadata.txt"), quote = FALSE, sep = "\t", row.names = FALSE)

pbmc@meta.data$cluster.ids


unique(pbmc@meta.data$cluster.ids)

unique(pbmc@meta.data$cell_ontology_class)




pdf(paste0(sample_name, "_tSNE_plot_all-seurat-clusters_paper.pdf"), width = 9, height = 7, useDingbats = FALSE)
p <- DimPlot(pbmc, reduction = "tsne", label = TRUE, group.by = "cluster.ids", pt.size = 0.5, label.size = 2.7) +
  theme(legend.text = element_text(size = 7))
print(p)
dev.off()





if (sample_name == "TM_facs_Lung") { 

# Endothelial Cell for all Endothelial
# use cluster id in legend, not in plot

pbmc$Kat_cell_type_1 <- "error"

pbmc$Kat_cell_type_1[pbmc$cluster.ids == "0"] <- "Endothelial Cell"
pbmc$Kat_cell_type_1[pbmc$cluster.ids == "1"] <- "Endothelial Cell"
pbmc$Kat_cell_type_1[pbmc$cluster.ids == "2"] <- "Fibroblast"
pbmc$Kat_cell_type_1[pbmc$cluster.ids == "3"] <- "Endothelial Cell"
pbmc$Kat_cell_type_1[pbmc$cluster.ids == "4"] <- "Endothelial Cell"
pbmc$Kat_cell_type_1[pbmc$cluster.ids == "5"] <- "AT1/AT2 Epithelial Cell"
pbmc$Kat_cell_type_1[pbmc$cluster.ids == "6"] <- "Fibroblast"
pbmc$Kat_cell_type_1[pbmc$cluster.ids == "7"] <- "Monocyte"
pbmc$Kat_cell_type_1[pbmc$cluster.ids == "8"] <- "Dendritic Cell"
pbmc$Kat_cell_type_1[pbmc$cluster.ids == "9"] <- "Fibroblast"
pbmc$Kat_cell_type_1[pbmc$cluster.ids == "10"] <- "Endothelial Cell"
pbmc$Kat_cell_type_1[pbmc$cluster.ids == "11"] <- "Macrophage"
pbmc$Kat_cell_type_1[pbmc$cluster.ids == "12"] <- "Endothelial Cell"
pbmc$Kat_cell_type_1[pbmc$cluster.ids == "13"] <- "B Cell"
pbmc$Kat_cell_type_1[pbmc$cluster.ids == "14"] <- "T Cell"
pbmc$Kat_cell_type_1[pbmc$cluster.ids == "15"] <- "Fibroblast"
pbmc$Kat_cell_type_1[pbmc$cluster.ids == "16"] <- "Smooth Muscle Cell"
pbmc$Kat_cell_type_1[pbmc$cluster.ids == "17"] <- "remove"
pbmc$Kat_cell_type_1[pbmc$cluster.ids == "18"] <- "Endothelial Cell"
pbmc$Kat_cell_type_1[pbmc$cluster.ids == "19"] <- "NK Cell"
pbmc$Kat_cell_type_1[pbmc$cluster.ids == "20"] <- "Granulocyte"
pbmc$Kat_cell_type_1[pbmc$cluster.ids == "21"] <- "Perictye"
pbmc$Kat_cell_type_1[pbmc$cluster.ids == "22"] <- "Ciliated Epithelial Cell"

unique(pbmc$Kat_cell_type_1)



# pbmc <- subset(pbmc, subset = Kat_cell_type_1 != "remove")



pbmc$Kat_cell_type_id_1 <- "error"

pbmc$Kat_cell_type_id_1[pbmc$cluster.ids == "0"] <- "0: Endothelial Cell"
pbmc$Kat_cell_type_id_1[pbmc$cluster.ids == "1"] <- "1: Endothelial Cell"
pbmc$Kat_cell_type_id_1[pbmc$cluster.ids == "2"] <- "2: Fibroblast"
pbmc$Kat_cell_type_id_1[pbmc$cluster.ids == "3"] <- "3: Endothelial Cell"
pbmc$Kat_cell_type_id_1[pbmc$cluster.ids == "4"] <- "4: Endothelial Cell"
pbmc$Kat_cell_type_id_1[pbmc$cluster.ids == "5"] <- "5: AT1/AT2 Epithelial Cell"
pbmc$Kat_cell_type_id_1[pbmc$cluster.ids == "6"] <- "6: Fibroblast"
pbmc$Kat_cell_type_id_1[pbmc$cluster.ids == "7"] <- "7: Monocyte"
pbmc$Kat_cell_type_id_1[pbmc$cluster.ids == "8"] <- "8: Dendritic Cell"
pbmc$Kat_cell_type_id_1[pbmc$cluster.ids == "9"] <- "9: Fibroblast"
pbmc$Kat_cell_type_id_1[pbmc$cluster.ids == "10"] <- "10: Endothelial Cell"
pbmc$Kat_cell_type_id_1[pbmc$cluster.ids == "11"] <- "11: Macrophage"
pbmc$Kat_cell_type_id_1[pbmc$cluster.ids == "12"] <- "12: Endothelial Cell"
pbmc$Kat_cell_type_id_1[pbmc$cluster.ids == "13"] <- "13: B Cell"
pbmc$Kat_cell_type_id_1[pbmc$cluster.ids == "14"] <- "14: T Cell"
pbmc$Kat_cell_type_id_1[pbmc$cluster.ids == "15"] <- "15: Fibroblast"
pbmc$Kat_cell_type_id_1[pbmc$cluster.ids == "16"] <- "16: Smooth Muscle Cell"
# pbmc$Kat_cell_type[pbmc$cluster.ids == "17"] <- "remove"
pbmc$Kat_cell_type_id_1[pbmc$cluster.ids == "18"] <- "18: Endothelial Cell"
pbmc$Kat_cell_type_id_1[pbmc$cluster.ids == "19"] <- "19: NK Cell"
pbmc$Kat_cell_type_id_1[pbmc$cluster.ids == "20"] <- "20: Granulocyte"
pbmc$Kat_cell_type_id_1[pbmc$cluster.ids == "21"] <- "21: Perictye"
pbmc$Kat_cell_type_id_1[pbmc$cluster.ids == "22"] <- "22: Ciliated Epithelial Cell"




pbmc$Kat_cell_type_id_1 <- factor(pbmc$Kat_cell_type_id_1, 
                                levels = c(
                                  
                                  "0: Endothelial Cell",
                                  "1: Endothelial Cell",
                                  "2: Fibroblast",
                                  "3: Endothelial Cell",
                                  "4: Endothelial Cell",
                                  "5: AT1/AT2 Epithelial Cell",
                                  "6: Fibroblast",
                                  "7: Monocyte",
                                  "8: Dendritic Cell",
                                  "9: Fibroblast",
                                  "10: Endothelial Cell",
                                  "11: Macrophage",
                                  "12: Endothelial Cell",
                                  "13: B Cell",
                                  "14: T Cell",
                                  "15: Fibroblast",
                                  "16: Smooth Muscle Cell",

                                  "18: Endothelial Cell",
                                  "19: NK Cell",
                                  "20: Granulocyte",
                                  "21: Perictye",
                                  "22: Ciliated Epithelial Cell"
                                ))



unique(pbmc$Kat_cell_type_id_1)





pbmc$Kat_cell_type <- "error"

pbmc$Kat_cell_type[pbmc$cluster.ids == "0"] <- "Endothelial cell"
pbmc$Kat_cell_type[pbmc$cluster.ids == "1"] <- "Endothelial cell"
pbmc$Kat_cell_type[pbmc$cluster.ids == "2"] <- "Fibroblast"
pbmc$Kat_cell_type[pbmc$cluster.ids == "3"] <- "Endothelial cell"
pbmc$Kat_cell_type[pbmc$cluster.ids == "4"] <- "Endothelial cell"
pbmc$Kat_cell_type[pbmc$cluster.ids == "5"] <- "AT1/AT2 epithelial cell"
pbmc$Kat_cell_type[pbmc$cluster.ids == "6"] <- "Fibroblast"
pbmc$Kat_cell_type[pbmc$cluster.ids == "7"] <- "Monocyte"
pbmc$Kat_cell_type[pbmc$cluster.ids == "8"] <- "Dendritic cell"
pbmc$Kat_cell_type[pbmc$cluster.ids == "9"] <- "Fibroblast"
pbmc$Kat_cell_type[pbmc$cluster.ids == "10"] <- "Endothelial cell"
pbmc$Kat_cell_type[pbmc$cluster.ids == "11"] <- "Macrophage"
pbmc$Kat_cell_type[pbmc$cluster.ids == "12"] <- "Endothelial cell"
pbmc$Kat_cell_type[pbmc$cluster.ids == "13"] <- "B cell"
pbmc$Kat_cell_type[pbmc$cluster.ids == "14"] <- "T cell"
pbmc$Kat_cell_type[pbmc$cluster.ids == "15"] <- "Fibroblast"
pbmc$Kat_cell_type[pbmc$cluster.ids == "16"] <- "Smooth muscle cell"
pbmc$Kat_cell_type[pbmc$cluster.ids == "17"] <- "remove"
pbmc$Kat_cell_type[pbmc$cluster.ids == "18"] <- "Endothelial cell"
pbmc$Kat_cell_type[pbmc$cluster.ids == "19"] <- "NK cell"
pbmc$Kat_cell_type[pbmc$cluster.ids == "20"] <- "Granulocyte"
pbmc$Kat_cell_type[pbmc$cluster.ids == "21"] <- "Perictye"
pbmc$Kat_cell_type[pbmc$cluster.ids == "22"] <- "Ciliated epithelial cell"




pbmc <- subset(pbmc, subset = Kat_cell_type != "remove")


sort(unique(pbmc$Kat_cell_type))


pbmc$cell_type_number <- "error"


pbmc$cell_type_number[pbmc$Kat_cell_type == "AT1/AT2 epithelial cell"] <- "1"   
pbmc$cell_type_number[pbmc$Kat_cell_type == "B cell"] <- "2"                   
pbmc$cell_type_number[pbmc$Kat_cell_type == "Ciliated epithelial cell"] <- "3" 
pbmc$cell_type_number[pbmc$Kat_cell_type == "Dendritic cell"] <- "4"          
pbmc$cell_type_number[pbmc$Kat_cell_type == "Endothelial cell"] <- "5"         
pbmc$cell_type_number[pbmc$Kat_cell_type == "Fibroblast"] <- "6"               
pbmc$cell_type_number[pbmc$Kat_cell_type == "Granulocyte"] <- "7"              
pbmc$cell_type_number[pbmc$Kat_cell_type == "Macrophage"] <- "8"              
pbmc$cell_type_number[pbmc$Kat_cell_type == "Monocyte"] <- "9"                 
pbmc$cell_type_number[pbmc$Kat_cell_type == "NK cell"] <- "10"                  
pbmc$cell_type_number[pbmc$Kat_cell_type == "Perictye"] <- "11"                 
pbmc$cell_type_number[pbmc$Kat_cell_type == "Smooth muscle cell"] <- "12"      
pbmc$cell_type_number[pbmc$Kat_cell_type == "T cell"] <- "13" 



cell_type_number_levels <- as.character(c(1:13))

pbmc$cell_type_number <- factor(pbmc$cell_type_number, levels = cell_type_number_levels)




pbmc$cell_type_legend <- paste(pbmc$cell_type_number, pbmc$Kat_cell_type, sep = ": ")


sort(unique(pbmc$cell_type_legend))




sort(unique(pbmc$cell_type_legend))

cell_type_legend_levels <- c(
  
"1: AT1/AT2 epithelial cell",  
"2: B cell",                  
"3: Ciliated epithelial cell", 
"4: Dendritic cell",           
"5: Endothelial cell",        
"6: Fibroblast",               
"7: Granulocyte",              
"8: Macrophage",              
"9: Monocyte",  
"10: NK cell",                 
"11: Perictye",               
"12: Smooth muscle cell",      
"13: T cell"   
  
)


pbmc$cell_type_legend <- factor(pbmc$cell_type_legend, levels = cell_type_legend_levels)




pbmc$Kat_cell_type_plot <- as.character(pbmc$Kat_cell_type)

sort(unique(pbmc$Kat_cell_type_plot))


pbmc$Kat_cell_type_plot[pbmc$Kat_cell_type_plot == "AT1/AT2 epithelial cell"] <- "AT1/AT2 EC"  
pbmc$Kat_cell_type_plot[pbmc$Kat_cell_type_plot == "Ciliated epithelial cell"] <- "Ciliated EC" 
pbmc$Kat_cell_type_plot[pbmc$Kat_cell_type_plot == "Smooth muscle cell"] <- "SM cell"      


}






if (sample_name == "TM_facs_Large_Intestine") {
  
  unique(pbmc@meta.data$cell_ontology_class)
  
  
  pbmc$Kat_cell_type_1 <- "error"
  
  pbmc$Kat_cell_type_1[pbmc$cell_ontology_class == "epithelial cell of large intestine"] <- "epithelial cell"
  pbmc$Kat_cell_type_1[pbmc$cell_ontology_class == "large intestine goblet cell"] <- "goblet cell"
  pbmc$Kat_cell_type_1[pbmc$cell_ontology_class == "Brush cell of epithelium proper of large intestine"] <- "Brush cell of epithelium proper"
  pbmc$Kat_cell_type_1[pbmc$cell_ontology_class == "enterocyte of epithelium of large intestine"] <- "enterocyte of epithelium"
  pbmc$Kat_cell_type_1[pbmc$cell_ontology_class == "enteroendocrine cell"] <- "enteroendocrine cell"
  
  unique(pbmc$Kat_cell_type_1)
  
  pbmc$Kat_cell_type_id_1 <- pbmc$Kat_cell_type_1
  
  
  
  pbmc$Kat_cell_type <- "error"
  
  pbmc$Kat_cell_type[pbmc$cell_ontology_class == "epithelial cell of large intestine"] <- "Epithelial cell"
  pbmc$Kat_cell_type[pbmc$cell_ontology_class == "large intestine goblet cell"] <- "Goblet cell"
  pbmc$Kat_cell_type[pbmc$cell_ontology_class == "Brush cell of epithelium proper of large intestine"] <- "Brush cell of epithelium"
  pbmc$Kat_cell_type[pbmc$cell_ontology_class == "enterocyte of epithelium of large intestine"] <- "Enterocyte of epithelium"
  pbmc$Kat_cell_type[pbmc$cell_ontology_class == "enteroendocrine cell"] <- "Enteroendocrine cell"
  
  
  sort(unique(pbmc$Kat_cell_type))
  
  
pbmc$cell_type_number <- "error"

pbmc$cell_type_number[pbmc$Kat_cell_type == "Brush cell of epithelium"] <- "1" 
pbmc$cell_type_number[pbmc$Kat_cell_type == "Enterocyte of epithelium"] <- "2" 
pbmc$cell_type_number[pbmc$Kat_cell_type == "Enteroendocrine cell"] <- "3"     
pbmc$cell_type_number[pbmc$Kat_cell_type == "Epithelial cell"] <- "4"         
pbmc$cell_type_number[pbmc$Kat_cell_type == "Goblet cell"] <- "5"


pbmc$cell_type_legend <- paste(pbmc$cell_type_number, pbmc$Kat_cell_type, sep = ": ")


sort(unique(pbmc$cell_type_legend))



pbmc$Kat_cell_type_plot <- as.character(pbmc$Kat_cell_type)

sort(unique(pbmc$Kat_cell_type_plot))

  
  
}







if (length(unique(pbmc$tissue)) > 1) {
  
  stop("tissue list is greater than 1")
  
}












gene_list_df <- read.table(gene_list_file, sep = "\t", header = TRUE, stringsAsFactors = FALSE)

additional_genes_df <- data.frame(gene = additional_genes)










# ***** update gene list *****

# gene_list_df <- rbind(gene_list_df, additional_genes_df)
gene_list_df <- additional_genes_df



gene_list <- gene_list_df$gene

gene_list <- plot_genes





mgi_filename.md5 <- md5sum(mgi_file)

cat("MGI gene annotation file:", mgi_file, "\n")
cat("md5:", mgi_filename.md5, "\n")
cat("\n")
cat("\n")





mgi_all <- read.table(mgi_file, sep = "\t", header = TRUE, quote = "", comment.char = "", stringsAsFactors = FALSE)

select_mgi_columns <- c(
  
  "MGI_Marker_Accession_ID", 
  "Marker_Symbol",           
  "Status",                  
  "Marker_Name",
  "Chromosome",              
  "Type", 
  "Feature_Types",           
  "Genome_Coordinate_Start",
  "Genome_Coordinate_End",   
  "Strand",                 
  "BioTypes" 
)


mgi <- mgi_all[ , select_mgi_columns]

mgi$gene <- mgi$Marker_Symbol
mgi$Marker_Symbol <- NULL

mgi <- mgi[mgi$Status == "O", ]
mgi <- mgi[mgi$Type %in% c("Gene", "Pseudogene"), ]

mgi_chrY.gene <- mgi$gene[mgi$Chromosome == "Y"]

pbmc[["percent.chrY"]] <- PercentageFeatureSet(pbmc, 
     features = mgi_chrY.gene[mgi_chrY.gene %in% row.names(pbmc@assays$RNA@counts)])








pbmc$cell_type <- pbmc$cell_ontology_class
pbmc$cell_type[is.na(pbmc$cell_type)] <- "NA"

pbmc$free_annotation_original <- pbmc$free_annotation
pbmc$free_annotation[is.na(pbmc$free_annotation)] <- "NA"



pbmc$subtissue <- as.character(pbmc$subtissue)

pbmc$subtissue[is.na(pbmc$subtissue)] <- "NA"
pbmc$subtissue[pbmc$subtissue == ""] <- "missing"




CL_cluster_tbl <- as.data.frame(table(pbmc$cell_type, pbmc$cluster.ids))

names(CL_cluster_tbl)[names(CL_cluster_tbl)=="Var1"] <- "CL_cell_type"
names(CL_cluster_tbl)[names(CL_cluster_tbl)=="Var2"] <- "cluster_id"
names(CL_cluster_tbl)[names(CL_cluster_tbl)=="Freq"] <- "cell_count"

CL_cluster_tbl$cluster_id <- as.numeric(as.character(CL_cluster_tbl$cluster_id))

CL_cluster_tbl <- CL_cluster_tbl[CL_cluster_tbl$cell_count != 0, ]
CL_cluster_tbl <- CL_cluster_tbl[order(CL_cluster_tbl$CL_cell_type, CL_cluster_tbl$cluster_id), ]

write.csv(CL_cluster_tbl, file = paste0(sample_name, "_CL_cluster-id_summary.csv"), row.names = FALSE)





skip <- TRUE

if (skip == FALSE) {

Idents(object = pbmc) <- "cell_type"
pbmc.markers <- FindAllMarkers(pbmc, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.csv(pbmc.markers, file = paste0(sample_name, "_CL_markers.csv"), row.names = FALSE)


Idents(object = pbmc) <- "cluster.ids"
pbmc.markers <- FindAllMarkers(pbmc, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.csv(pbmc.markers, file = paste0(sample_name, "_cluster-ids_markers.csv"), row.names = FALSE)

}




pbmc$cluster.ids_num <- as.numeric(pbmc$cluster.ids)


pdf(paste0(sample_name, "_tSNE_plot_cluster-id.pdf"), width = 9, height = 7, useDingbats = FALSE)
p <- DimPlot(pbmc, reduction = "tsne", label = TRUE, group.by = "cluster.ids_num", pt.size = 0.5, label.size = 2.7) +
  theme(legend.text = element_text(size = 7))
print(p)
dev.off()

write.csv(pbmc@meta.data, "test.csv")






# Cell Ontology version
pbmc$tissue_chr <- as.character(pbmc$tissue)

cell_type_df <- as.data.frame(table(pbmc@meta.data[ , c("tissue_chr", "cell_type")]))
names(cell_type_df)[names(cell_type_df)=="tissue_chr"] <- "tissue"

write.table(cell_type_df, file = paste0(sample_name, "_cell-type_table_CL.txt"), quote = FALSE, sep = "\t", row.names = FALSE)



# free annotation version
cell_type_fa_df <- as.data.frame(table(pbmc@meta.data[ , c("tissue_chr", "free_annotation")]))
names(cell_type_fa_df)[names(cell_type_fa_df)=="tissue_chr"] <- "tissue"

write.table(cell_type_fa_df, file = paste0(sample_name, "_cell-type_table_FA.txt"), quote = FALSE, sep = "\t", row.names = FALSE)






pdf(paste0(sample_name, "_tSNE_plot_cell_type.pdf"), width = 9, height = 7, useDingbats = FALSE)
p <- DimPlot(pbmc, reduction = "tsne", label = TRUE, group.by = "cell_type", pt.size = 0.5, label.size = 2.7) +
  theme(legend.text = element_text(size = 7))
print(p)
dev.off()

pdf(paste0(sample_name, "_tSNE_plot_cell_type_no-legend.pdf"), width = 9, height = 7, useDingbats = FALSE)
p <- DimPlot(pbmc, reduction = "tsne", label = TRUE, group.by = "cell_type", pt.size = 1) +
  ggtitle(tissue_name) +
                          NoLegend()
print(p)
dev.off()






pdf(paste0(sample_name, "_tSNE_plot_cell-type_CL.pdf"), width = 9, height = 7, useDingbats = FALSE)
p <- DimPlot(pbmc, reduction = "tsne", label = TRUE, group.by = "cell_type", pt.size = 0.5, label.size = 2.7) +
  theme(legend.text = element_text(size = 7)) +
  theme(plot.subtitle = element_text(hjust = 0.5)) +
  # ggtitle(sample_name) +
  theme(plot.subtitle = element_text(hjust = 0.5)) +
  labs(subtitle = paste0("Cell Ontology", " - ", paste0(tissue_name, " (", organism_name, ")"))) +
  labs(caption = dataset_name)
print(p)
dev.off()


pdf(paste0(sample_name, "_tSNE_plot_cell-type_CL_no-legend.pdf"), width = 9, height = 7, useDingbats = FALSE)
p <- DimPlot(pbmc, reduction = "tsne", label = TRUE, group.by = "cell_type", pt.size = 0.5) +
  NoLegend() +
  theme(plot.subtitle = element_text(hjust = 0.5)) +
  # ggtitle(sample_name) +
  labs(subtitle = paste0("Cell Ontology", " - ", paste0(tissue_name, " (", organism_name, ")"))) +
  labs(caption = dataset_name)
print(p)
dev.off()



pdf(paste0(sample_name, "_tSNE_plot_cell-type_FA.pdf"), width = 9, height = 7, useDingbats = FALSE)
p <- DimPlot(pbmc, reduction = "tsne", label = TRUE, group.by = "free_annotation", pt.size = 0.5, label.size = 2.7) +
  theme(legend.text = element_text(size = 7)) +
  theme(plot.subtitle = element_text(hjust = 0.5)) +
  # ggtitle(sample_name) +
  labs(subtitle = paste0("free annotation", " - ", paste0(tissue_name, " (", organism_name, ")"))) +
  labs(caption = dataset_name)
print(p)
dev.off()

pdf(paste0(sample_name, "_tSNE_plot_cell-type_FA_no-legend.pdf"), width = 9, height = 7, useDingbats = FALSE)
p <- DimPlot(pbmc, reduction = "tsne", label = TRUE, group.by = "free_annotation", pt.size = 0.5) +
  NoLegend() +
  theme(plot.subtitle = element_text(hjust = 0.5)) +
  # ggtitle(sample_name) +
  labs(subtitle = paste0("free annotation", " - ", paste0(tissue_name, " (", organism_name, ")"))) +
  labs(caption = dataset_name)
print(p)
dev.off()



pdf(paste0(sample_name, "_tSNE_plot_cell-type_Kat.pdf"), width = 9, height = 7, useDingbats = FALSE)
p <- DimPlot(pbmc, reduction = "tsne", label = TRUE, group.by = "Kat_cell_type", pt.size = 0.5, label.size = 2.7) +
  theme(legend.text = element_text(size = 7)) +
  theme(plot.subtitle = element_text(hjust = 0.5)) +
  # ggtitle(sample_name) +
  labs(subtitle = paste0("Cell Type", " - ", paste0(tissue_name, " (", organism_name, ")"))) +
  labs(caption = dataset_name)
print(p)
dev.off()









# pdf(paste0(sample_name, "_tSNE_plot_cell-type-id_Kat_paper.pdf"), width = 9, height = 7, useDingbats = FALSE)
# p <- DimPlot(pbmc, reduction = "tsne", label = FALSE, group.by = "Kat_cell_type_id_1", pt.size = 0.5, label.size = 2.7) +
#   theme(legend.text = element_text(size = 7)) +
#   theme(plot.subtitle = element_text(hjust = 0.5)) +
#   # ggtitle(sample_name) +
#   labs(subtitle = paste0("Cluster ID", " - ", paste0(tissue_name, " (", organism_name, ")"))) +
#   labs(caption = dataset_name)
# print(p)
# dev.off()


# need this tSNE plots to identify a groupings of cell types for paper

pdf(paste0(sample_name, "_tSNE_plot_cell-type-id_Kat_label_paper.pdf"), width = 9, height = 7, useDingbats = FALSE)
p <- DimPlot(pbmc, reduction = "tsne", label = TRUE, group.by = "Kat_cell_type_id_1", pt.size = 0.5, label.size = 2.7) +
  theme(legend.text = element_text(size = 7)) +
  theme(plot.subtitle = element_text(hjust = 0.5)) +
  # ggtitle(sample_name) +
  labs(subtitle = paste0("Cluster ID", " - ", paste0(tissue_name, " (", organism_name, ")"))) +
  labs(caption = dataset_name)
print(p)
dev.off()







pdf(paste0(sample_name, "_tSNE_plot_cell-type_number_fig.pdf"), width = 9, height = 7, useDingbats = FALSE)
p <- DimPlot(pbmc, reduction = "tsne", label = TRUE, group.by = "cell_type_number", pt.size = 0.5, label.size = 5) +
  
  ggtitle(plot_title) +
  # theme(legend.text = element_text(size = 7)) +
  # theme(plot.subtitle = element_text(hjust = 0.5)) +
  # labs(subtitle = paste0("Cell Ontology", " - ", paste0(tissue_name, " (", organism_name, ")"))) +
  # labs(caption = dataset_name)
  theme(legend.position = "none",
          plot.title = element_text(size = 16, face = "plain"),
        axis.ticks.length = unit(0.1, "inch"),
        axis.title = element_text(size = 16, face = "plain"),
        axis.text = element_text(size = 16),
        
        )


xlim_all_samples <- layer_scales(p)$x$range$range
ylim_all_samples <- layer_scales(p)$y$range$range


print(p)
dev.off()




# https://github.com/satijalab/seurat/issues/3899


pdf(paste0(sample_name, "_tSNE_plot_cell-type_legend_fig.pdf"), width = 9, height = 7, useDingbats = FALSE)
p <- DimPlot(pbmc, reduction = "tsne", label = TRUE, group.by = "cell_type_legend", pt.size = 0.5, label.size = 2.7) +
  theme(legend.text = element_text(size = 7)) +
  theme(plot.subtitle = element_text(hjust = 0.5)) +
  labs(subtitle = paste0("Cell Ontology", " - ", paste0(tissue_name, " (", organism_name, ")"))) +
  labs(caption = dataset_name) +
  theme(
    legend.text = element_text(size = 16, margin = margin(l = 1, unit = "pt")),
    legend.box.spacing = unit(0, "pt"),
    
    # legend.key.height= unit(2, "cm"),
    
    # legend.key.size = unit(0.5, "cm"),
    
    legend.key.spacing.y = unit(0.1, "cm")

  ) +
  
  guides(color = guide_legend(override.aes = list(size = 10)))




print(p)
dev.off()





pdf(paste0(sample_name, "_tSNE_plot_individual-cell-types.pdf"), width = 9, height = 7, useDingbats = FALSE)


cell_type_names <- sort(unique(pbmc$cell_type_legend))
sample_colors <- hue_pal()(length(cell_type_names))


Idents(object = pbmc) <- "cell_type_legend"

for (current_cell_type in cell_type_names) {
  
cluster_position <- which(cell_type_names == current_cell_type)

current_cells <- WhichCells(object = pbmc, idents = current_cell_type)

p <- DimPlot(pbmc, reduction = "tsne", label = FALSE, group.by = "cell_type_legend", cells = current_cells, pt.size = 0.5, cols = sample_colors[cluster_position]) +
                          xlim(xlim_all_samples[1], xlim_all_samples[2]) +
                          ylim(ylim_all_samples[1], ylim_all_samples[2]) +
                          # NoLegend() +
                          ggtitle(plot_title)
print(p)
  
}

dev.off()


















pdf(paste0(sample_name, "_tSNE_plot_tissue.pdf"), width = 9, height = 7, useDingbats = FALSE)
p <- DimPlot(pbmc, reduction = "tsne", label = TRUE, group.by = "tissue", pt.size = 0.5, label.size = 2.7) +
  theme(legend.text = element_text(size = 7))
print(p)
dev.off()


pdf(paste0(sample_name, "_tSNE_plot_subtissue.pdf"), width = 9, height = 7, useDingbats = FALSE)
p <- DimPlot(pbmc, reduction = "tsne", label = TRUE, group.by = "subtissue", pt.size = 0.5, label.size = 2.7) +
  theme(legend.text = element_text(size = 7))
print(p)
dev.off()



found_genes <- c()
found_genes <- gene_list_df[gene_list_df$gene %in% row.names(GetAssayData(pbmc)), "gene"]

missing_genes_df <- gene_list_df[which(!(gene_list_df$gene %in% found_genes)), ]

write.table(missing_genes_df, file = paste0(sample_name, "_missing_genes.txt"), quote = FALSE, sep = "\t", row.names = FALSE)

















pdf(paste(sample_name, gene_type, "DotPlot_no-scale_x-axis-genes.pdf", sep = "_"), width = 20, height = 7)

p <- DotPlot(pbmc, features = found_genes, group.by = "cell_type", scale = FALSE) + RotatedAxis()
  # scale_x_discrete(labels = found_genes_convert)

print(p)
dev.off()


pdf(paste(sample_name, gene_type, "DotPlot_scaled_x-axis-genes.pdf", sep = "_"), width = 15, height = 5)

p <- DotPlot(pbmc, features = found_genes, group.by = "cell_type", scale = TRUE) + RotatedAxis() +
  ggtitle(tissue_name)

  # scale_x_discrete(labels = found_genes_convert)

print(p)
dev.off()






found_genes_paper <- c()
found_genes_paper <- plot_genes[plot_genes %in% row.names(GetAssayData(pbmc))]

plot_genes_df <- as.data.frame(plot_genes)

missing_genes_paper <- plot_genes[which(!(plot_genes %in% found_genes_paper))]

write.table(missing_genes_paper, file = paste0(sample_name, "_missing_genes_paper-version.txt"), quote = FALSE, sep = "\t", row.names = FALSE)






# dot plots for paper

# https://rpubs.com/eraz0001/enhanced_dotplot

# https://github.com/satijalab/seurat/issues/3914




pdf(paste(sample_name, gene_type, "DotPlot_FA.pdf", sep = "_"), width = 9, height = 7)

# p <- DotPlot(pbmc.temp, features = found_genes_paper, group.by = "free_annotation", scale = TRUE, dot.scale = 1) + coord_flip() +
p <- DotPlot(pbmc, features = found_genes_paper, group.by = "free_annotation", scale = TRUE) + coord_flip() +
  ggtitle(paste0(tissue_name)) +
  xlab("Gene") + 
  ylab("Cell Type (free annotation)") + 
  labs(caption = dataset_name) +
  # scale_x_discrete(labels = found_genes_convert_paper) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

print(p)
dev.off()

dotplot_data <- p$data
write.csv(dotplot_data, file = paste(sample_name, gene_type, "DotPlot_FA_data.csv", sep = "_"), row.names = FALSE)






pdf(paste(sample_name, gene_type, "DotPlot_CL.pdf", sep = "_"), width = 9, height = 7)

p <- DotPlot(pbmc, features = found_genes_paper, group.by = "cell_type", scale = TRUE) + coord_flip() +
  ggtitle(paste0(tissue_name)) +
  xlab("Gene") + 
  ylab("Cell Type (Cell Ontology)") + 
  labs(caption = dataset_name) +
  # scale_x_discrete(labels = found_genes_convert_paper) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

print(p)
dev.off()


dotplot_data <- p$data
write.csv(dotplot_data, file = paste(sample_name, gene_type, "DotPlot_CL_data.csv", sep = "_"), row.names = FALSE)






pdf(paste(sample_name, gene_type, "DotPlot_Kat-cell-type.pdf", sep = "_"), width = 9, height = 7)

p <- DotPlot(pbmc, features = found_genes_paper, group.by = "Kat_cell_type", scale = TRUE) + coord_flip() +
  ggtitle(paste0(tissue_name)) +
  xlab("Gene") + 
  ylab("Cell Type (Kat)") + 
  labs(caption = dataset_name) +
  # scale_x_discrete(labels = found_genes_convert_paper) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

print(p)
dev.off()


dotplot_data <- p$data
write.csv(dotplot_data, file = paste(sample_name, gene_type, "DotPlot_Kat-cell-type_data.csv", sep = "_"), row.names = FALSE)




pdf(paste(sample_name, gene_type, "DotPlot_Kat-cell-type_fig.pdf", sep = "_"), width = 9, height = 7)

p <- DotPlot(pbmc, features = found_genes_paper, group.by = "Kat_cell_type_plot", scale = TRUE) + coord_flip() +
  ggtitle(plot_title) +
  # xlab("Gene") + 
  # ylab("Cell Type (Kat)") + 
  # labs(caption = dataset_name) +
  # scale_x_discrete(labels = found_genes_convert_paper) +
  
  # theme(axis.text.x = element_text(angle = 90, hjust = 0.95, vjust = 0.5),
        theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.text.y = element_text(face = "italic"),
        axis.title = element_blank(),     
        plot.title = element_text(size = 20, face = "plain", hjust = 0.5, vjust = 0),
      axis.ticks.length = unit(0.075, "inch"),
        # axis.title = element_text(size = 16, face = "plain"),
        axis.text = element_text(size = 16)
        
        
        )

print(p)
dev.off()


dotplot_data <- p$data
write.csv(dotplot_data, file = paste(sample_name, gene_type, "DotPlot_Kat-cell-type_data_fig.csv", sep = "_"), row.names = FALSE)


















compartment_list <- unique(pbmc$tissue)
cell_type_list <- unique(pbmc$cell_type)

unique(pbmc$subtissue)









cat("Making violin plots . . .", "\n\n")


summary_df <- data.frame()

compartment_cell_types_df <- data.frame()




# annotation_type <- c("cell_type", "free_annotation", "Kat_cell_type")
annotation_type <- c("Kat_cell_type")


for (current_annotation_type in annotation_type) {
  
  
  message(current_annotation_type)
  
  
  anno_short <- "error"
  if (current_annotation_type == "cell_type") {
    anno_short <- "CL"
  } else if (current_annotation_type == "free_annotation") {
    anno_short <- "FA"
  } else if (current_annotation_type == "Kat_cell_type") {
    anno_short <- "Kat"
  }
  
  
  anno_long <- "error"
  if (current_annotation_type == "cell_type") {
    anno_long <- "Cell Ontology"
  } else if (current_annotation_type == "free_annotation") {
    anno_long <- "free annotation"
  } else if (current_annotation_type == "Kat_cell_type") {
    anno_long <- "Kat"
  }








  
 
  found_genes <- c()
  found_genes <- gene_list_df[gene_list_df$gene %in% row.names(GetAssayData(pbmc)), "gene"]


  
  # raster_pdf(paste0(gene_type, "_", cell_group_name, "_", current_tissue, "_cell-type_violin_plots.pdf"), res = 300)

  # pdf(paste(sample_name, gsub(" ", "_", current_compartment), gene_type, "cell-type_violin_plots.pdf", sep = "_"), useDingbats = FALSE)
  
  pdf(paste0(sample_name, "_", gene_type, "_violin_plots_", anno_short, "_fig.pdf"), useDingbats = FALSE, width = 9, height = 7)
  
  
    
  # pdf(paste0("NLR-genes_", current_dataset, "_", gsub("/", "_", current_cell_type), "_violin_plots.pdf"), useDingbats = FALSE)
  
  
  for (current_gene in found_genes) {
    
    Idents(object = pbmc) <- "compartment"
    
    
    # matching_gene <- gene_list_df[gene_list_hgnc$ensembl_gene_id == current_gene, "gene"]
    
    

    vln_plot_title <- paste0(current_gene, "\n", tissue_name)
    
    
    message(current_gene)
    
     # message(paste0(matching_gene, " (", current_gene, ")"))
    
    # message(vln_plot_title)
    
    
    
    all_genes <- pbmc@assays$RNA@data@Dimnames[[1]]
    
    gene_test <- all_genes[which(all_genes == current_gene)]
    
    if (length(gene_test) != 1) {
      
      stop("current gene length is not 1")
    }
    
    
    
    p <- VlnPlot(object = pbmc, features = current_gene, group.by = current_annotation_type, pt.size = 1, raster = TRUE) +
      ggtitle(vln_plot_title) +
      xlab(paste0("Cell Type (", anno_long, ")")) +
      ylab("Expression Level") +
      labs(caption = dataset_name) +
      theme(legend.position = "none") +
      theme(axis.text.x = element_text(size = rel(0.8)))

    p$layers[[2]]$aes_params$alpha <- 0.2

    print(p)
    
    
    
    
    p <- VlnPlot(object = pbmc, features = current_gene, group.by = "Kat_cell_type_plot", pt.size = 1, raster = TRUE) + 
      ggtitle(current_gene) +
      # xlab(paste0("Cell Type (", anno_long, ")")) +
      ylab(paste(plot_title, "Expression Level", sep = "\n")) + 
      # labs(caption = dataset_name) +
      theme(legend.position = "none", 
         axis.title.x = element_blank(),
      plot.title = element_text(size = 20, face = "italic"),
      axis.ticks.length = unit(0.1, "inch"),
        axis.title = element_text(size = 16, face = "plain"),
        axis.text = element_text(size = 16)
      
      ) 
      
      # theme(axis.text.x = element_text(size = rel(0.8)))
    
    p$layers[[2]]$aes_params$alpha <- 0.2
    
    print(p)
    
    
    
    
    
  

    # p <- VlnPlot(object = pbmc.temp, features = current_gene, group.by = "cell_type", split.by = "subtissue", pt.size = 1) + 
    #   ggtitle(vln_plot_title) +
    #   xlab("Cell Type") + 
    #   ylab("Expression Level") + 
    #   # theme(legend.position = "none") +
    #   theme(axis.text.x = element_text(size = rel(0.8)))
    # 
    # p$layers[[2]]$aes_params$alpha <- 0.2
    # 
    # print(p)
    
    
    
    
  
    
    # cell_type_list <- unique(pbmc.temp$cell_type)
    # 
    # 
    # Idents(pbmc.temp) <- "cell_type"
    
    
    
    # ***** need to change for cell type (CL or FA) used *****
    # cell_type_list <- unique(pbmc.temp$free_annotation)
    cell_type_list <- unique(pbmc@meta.data[ , current_annotation_type])
    
    
    
    Idents(pbmc) <- current_annotation_type
    
  
  for (current_cell_type in cell_type_list) {
  
    
    current_cells <- WhichCells(object = pbmc, idents = current_cell_type)
  
  
            expression_values <- pbmc@assays$RNA@data[current_gene, current_cells]

          if (!is.vector(expression_values)) {
               
               stop("expression values are not a vector")
          }
          

          
          current_summary <- data.frame(gene = current_gene)
          
          # current_summary$ensembl_id <- current_gene
          
          current_summary$tissue <- tissue_name
          # current_summary$compartment <- current_compartment
          current_summary$cell_type <- current_cell_type
          
          current_summary$annotation_type <- anno_long
          
          current_summary$mean <- mean(expression_values)
          current_summary$median <- median(expression_values)
          current_summary$median_greater0 <- median(expression_values[which(expression_values > 0)])
          current_summary$cell_count <- length(expression_values)
          current_summary$pct_cells_0 <- (length(expression_values[which(expression_values == 0)])/length(expression_values)) * 100
          current_summary$pct_cells_greater0 <- (length(expression_values[which(expression_values > 0)])/length(expression_values)) * 100
        
          
          summary_df <- rbind(summary_df, current_summary)
  
  } # END for current cell type
  
  
  
  } # END for current gene
  
  
  dev.off()
  
  
} # END for annotation type
  





  write.csv(summary_df, file = paste(sample_name, gene_type, "summary.csv", sep = "_"), row.names = FALSE)


  # compartment_cell_types_df <- compartment_cell_types_df[order(compartment_cell_types_df$compartment), ]
  # 
  # write.table(compartment_cell_types_df, paste(sample_name, gene_type, "compartment_cell-types.txt", sep = "_"), quote = FALSE, sep = "\t", row.names = FALSE)
  


  












cat("\n")
cat("\n")
cat("------------------------------------------------------", "\n")
cat("\n")
cat("sessionInfo:", "\n")
cat("\n")

sessionInfo_object <- sessionInfo()
print(sessionInfo_object)


cat("\n")
cat("------------------------------------------------------", "\n")
cat("\n")
cat("\n")



script_end_time <- Sys.time()


script_run_time_seconds <- difftime(script_end_time, script_start_time, units = "secs")
script_run_time_minutes <- difftime(script_end_time, script_start_time, units = "mins")
script_run_time_hours <- difftime(script_end_time, script_start_time, units = "hours")

cat("\n")

if (script_run_time_hours > 1) {
  
  cat("R script run time:", script_run_time_hours, "hours", "\n")
  
} else if (script_run_time_minutes > 1) {
  
  cat("R script run time:", script_run_time_minutes, "minutes", "\n")
  
} else {
  
  cat("R script run time:", script_run_time_seconds, "seconds", "\n")
}

cat("\n")



script_end_time_display <- format(script_end_time, "%Y-%m-%d %H:%M:%S %Z")

cat("\n")

cat("R script completed at:", "\n")
cat(script_end_time_display, "\n")



# return output to the terminal 
sink()



