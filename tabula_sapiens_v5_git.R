



library(ggplot2)


library(rasterpdf, lib.loc = "/nas/longleaf/home/wrobel/R/local_R_libraries/R-4.2.2")



library(scales)


# for md5sum
library(tools)


# location of Seurat v5
library(Seurat, lib.loc = "/nas/longleaf/home/wrobel/R/local_R_libraries/R-4.3.1")





script_start_time <- Sys.time()








# *** START user input ***

split_report <- TRUE


# sample_name <- "TS_lung"
# dataset_name <- "TS_lung"
# plot_title <- "Human Lung"
# organism_name <- "human"
# tissue_name <- "lung"
# rds_input_filename <- "/work/users/w/r/wrobel/tinglab/tabula_sapiens/rds/lung/TS_lung_v5.rds"

sample_name <- "TS_large_intestine"
dataset_name <- "TS_large_intestine"
plot_title <- "Human Large Intestine"
organism_name <- "human"
tissue_name <- "large_intestine"
rds_input_filename <- "/work/users/w/r/wrobel/tinglab/tabula_sapiens/rds/large_intestine/TS_large_intestine_v5.rds"















gene_type <- "Kat-genes"

# gene_list_file <- "/proj/tinglab/users/wrobel/NLR_gene_list/NLR_gene_list_v2.txt"








cell_group_name <- "test"



# hgnc_filename <- "/proj/tinglab/users/wrobel/R/common_input/hgnc/hgnc_complete_set_2023-07-01.txt"
hgnc_filename <- "/proj/tinglab/users/wrobel/gene_datasets/hgnc/Oct2024/hgnc_complete_set_2024-10-01.txt"


additional_genes <- c(
  
  "ACTB",
  "IL1B",
  "IL18",
  "ACE2",
  
  "CASP1", 
  "CASP4", 
  "CASP5", 
  "TMPRSS2"
)


plot_genes <- c(
  
  "AIM2",
  "NLRP3",
  "NLRP6",
  
  "IL1B",
  "IL18",
  "ACE2",
  
  "CASP1", 
  "CASP4", 
  "CASP5", 
  "TMPRSS2"
)



plot_genes_Kat <- c(

  "IL18",
  "IL1B",
  "CASP5",
  "CASP4",
  "CASP1",
  "NLRP6",
  "NLRP3",
  "AIM2",
  "TMPRSS2",
  "ACE2"
)


plot_genes <- plot_genes_Kat


VlnPlot_genes <- c(
  
  "NLRP6",
  "NLRP3",
  "AIM2"
)




additional_genes <- c(plot_genes, "ACTB")




# *** END user input ****









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










# gene_list_df <- read.table(gene_list_file, sep = "\t", header = TRUE, stringsAsFactors = FALSE)

additional_genes_df <- data.frame(gene = additional_genes)


# ***** update gene list *****

# gene_list_df <- rbind(gene_list_df, additional_genes_df)
gene_list_df <- additional_genes_df


# gene_list <- gene_list_df$gene

# found_genes <- gene_list[which(gene_list %in% lung@assays$RNA@counts@Dimnames[[1]])]

gene_list <- plot_genes

# gene_list_df <- data.frame(gene = plot_genes)



hgnc_filename.md5 <- md5sum(hgnc_filename)

cat("HUGO gene annotation file:", hgnc_filename, "\n")
cat("md5:", hgnc_filename.md5, "\n")
cat("\n")
cat("\n")


hgnc <- read.table(hgnc_filename, sep = "\t", header = TRUE, quote = "", comment.char = "", stringsAsFactors = FALSE)

select_hgnc_columns <- c(
     
     "hgnc_id", "symbol", "name", "locus_group"
     
     # , "locus_type", "location", "location_sortable", 
     , "gene_group"
     # , "alias_symbol"
     # , "prev_symbol"
     # , "entrez_id"
     , "ensembl_gene_id"
     # , "ucsc_id"
     # , "refseq_accession"
     # , "ccds_id"
     # , "uniprot_ids"
)



hgnc <- hgnc[ , select_hgnc_columns]



gene_list_hgnc <- merge(x = gene_list_df, y = hgnc, by.x = "gene", by.y = "symbol", all.x = TRUE, sort = FALSE)





rds_input_filename.md5 <- md5sum(rds_input_filename)

cat("RDS input file:", rds_input_filename, "\n")
cat("md5:", rds_input_filename.md5, "\n")
cat("\n")
cat("\n")


pbmc <- readRDS(rds_input_filename)


pbmc.original <- pbmc



unique(pbmc@meta.data$cell_type)

unique(pbmc@meta.data$tissue)

unique(pbmc@meta.data$compartment)


Idents(object = pbmc) <- "tissue"




pbmc$tissue2 <- as.character(pbmc$tissue)
pbmc$compartment2 <- as.character(pbmc$compartment)
pbmc$cell_type2 <- as.character(pbmc$cell_type)

cell_type_df <- as.data.frame(table(pbmc@meta.data[ , c("tissue2", "compartment2", "cell_type2")]))

cell_type_df <- cell_type_df[cell_type_df$Freq != 0, ]

write.table(cell_type_df, file = paste0(sample_name, "_cell_type_table_CL.txt"), quote = FALSE, sep = "\t", row.names = FALSE)




# free annotation version
pbmc$free_annotation2 <- as.character(pbmc$free_annotation)

cell_type_fa_df <- as.data.frame(table(pbmc@meta.data[ , c("tissue2", "compartment2", "free_annotation2")]))

cell_type_fa_df <- cell_type_fa_df[cell_type_fa_df$Freq != 0, ]

write.table(cell_type_fa_df, file = paste0(sample_name, "_cell_type_table_FA.txt"), quote = FALSE, sep = "\t", row.names = FALSE)


unique(pbmc$cell_type)


# lung

if (sample_name == "TS_lung") {

Kat_endothelial_group <- c("blood vessel endothelial cell", 
                           "capillary endothelial cell", 
                           "endothelial cell of artery", 
                           "endothelial cell of lymphatic vessel", 
                           "lung microvascular endothelial cell", 
                           "vein endothelial cell")

Kat_monocyte_group <- c("classical monocyte", "intermediate monocyte", "non-classical monocyte")

Kat_epithelial_group <- c("basal cell", 
                          "club cell", 
                          "lung ciliated cell", 
                          "pulmonary ionocyte", 
                          "respiratory goblet cell", 
                          "serous cell of epithelium of bronchus",
                          "type I pneumocyte",
                          "type II pneumocyte")

Kat_T_cell_group <- c("effector CD4-positive, alpha-beta T cell",
                      "CD4-positive, alpha-beta T cell",
                      "effector CD8-positive, alpha-beta T cell",
                      "CD8-positive, alpha-beta T cell")

Kat_smooth_muscle_group <- c("bronchial smooth muscle cell", 
                             "smooth muscle cell", 
                             "vascular associated smooth muscle cell")



pbmc$Kat_cell_type_1 <- as.character(pbmc$cell_type)

pbmc$Kat_cell_type_1[pbmc$Kat_cell_type_1 %in% Kat_endothelial_group] <- "Endothelial cell"
pbmc$Kat_cell_type_1[pbmc$Kat_cell_type_1 %in% Kat_monocyte_group] <- "Monocyte"
pbmc$Kat_cell_type_1[pbmc$Kat_cell_type_1 %in% Kat_epithelial_group] <- "Epithelial cell"
pbmc$Kat_cell_type_1[pbmc$Kat_cell_type_1 %in% Kat_T_cell_group] <- "T cell"
pbmc$Kat_cell_type_1[pbmc$Kat_cell_type_1 %in% Kat_smooth_muscle_group] <- "Smooth Muscle cell"

unique(pbmc$Kat_cell_type_1)



pbmc$Kat_cell_type <- as.character(pbmc$cell_type)


pbmc$Kat_cell_type[pbmc$Kat_cell_type %in% Kat_endothelial_group] <- "Endothelial cell"
pbmc$Kat_cell_type[pbmc$Kat_cell_type %in% Kat_monocyte_group] <- "Monocyte"
pbmc$Kat_cell_type[pbmc$Kat_cell_type %in% Kat_epithelial_group] <- "Epithelial cell"
pbmc$Kat_cell_type[pbmc$Kat_cell_type %in% Kat_T_cell_group] <- "T cell"
pbmc$Kat_cell_type[pbmc$Kat_cell_type %in% Kat_smooth_muscle_group] <- "Smooth muscle cell"

sort(unique(pbmc$Kat_cell_type))



pbmc$Kat_cell_type[pbmc$cell_type == "adventitial cell"] <- "Adventitial cell"            
pbmc$Kat_cell_type[pbmc$cell_type == "basophil"] <- "Basophil"        
pbmc$Kat_cell_type[pbmc$cell_type == "dendritic cell"] <- "Dendritic cell"              
pbmc$Kat_cell_type[pbmc$cell_type == "fibroblast"] <- "Fibroblast"                   
pbmc$Kat_cell_type[pbmc$cell_type == "macrophage"] <- "Macrophage"                  
pbmc$Kat_cell_type[pbmc$cell_type == "mature NK T cell"] <- "Mature NK T cell"             
pbmc$Kat_cell_type[pbmc$cell_type == "mesothelial cell"] <- "Mesothelial cell"             
pbmc$Kat_cell_type[pbmc$cell_type == "myofibroblast cell"] <- "Myofibroblast cell"          
pbmc$Kat_cell_type[pbmc$cell_type == "neutrophil"] <- "Neutrophil"                   
pbmc$Kat_cell_type[pbmc$cell_type == "pericyte"] <- "Pericyte"                     
pbmc$Kat_cell_type[pbmc$cell_type == "plasma cell"] <- "Plasma cell"                 
pbmc$Kat_cell_type[pbmc$cell_type == "plasmacytoid dendritic cell"] <- "Plasmacytoid dendritic cell" 


sort(unique(pbmc$Kat_cell_type))




pbmc$cell_type_number <- "error"


pbmc$cell_type_number[pbmc$Kat_cell_type == "Adventitial cell"] <- "1"            
pbmc$cell_type_number[pbmc$Kat_cell_type == "B cell"] <- "2"                      
pbmc$cell_type_number[pbmc$Kat_cell_type == "Basophil"] <- "3"                    
pbmc$cell_type_number[pbmc$Kat_cell_type == "Dendritic cell"] <- "4"             
pbmc$cell_type_number[pbmc$Kat_cell_type == "Endothelial cell"] <- "5"            
pbmc$cell_type_number[pbmc$Kat_cell_type == "Epithelial cell"] <- "6"             
pbmc$cell_type_number[pbmc$Kat_cell_type == "Fibroblast"] <- "7"                  
pbmc$cell_type_number[pbmc$Kat_cell_type == "Macrophage"] <- "8"                
pbmc$cell_type_number[pbmc$Kat_cell_type == "Mature NK T cell"] <- "9"           
pbmc$cell_type_number[pbmc$Kat_cell_type == "Mesothelial cell"] <- "10"           
pbmc$cell_type_number[pbmc$Kat_cell_type == "Monocyte"] <- "11"                    
pbmc$cell_type_number[pbmc$Kat_cell_type == "Myofibroblast cell"] <- "12"         
pbmc$cell_type_number[pbmc$Kat_cell_type == "Neutrophil"] <- "13"                  
pbmc$cell_type_number[pbmc$Kat_cell_type == "Pericyte"] <- "14"                    
pbmc$cell_type_number[pbmc$Kat_cell_type == "Plasma cell"] <- "15"                 
pbmc$cell_type_number[pbmc$Kat_cell_type == "Plasmacytoid dendritic cell"] <- "16"
pbmc$cell_type_number[pbmc$Kat_cell_type == "Smooth muscle cell"] <- "17"          
pbmc$cell_type_number[pbmc$Kat_cell_type == "T cell"] <- "18" 



cell_type_number_levels <- as.character(c(1:18))

pbmc$cell_type_number <- factor(pbmc$cell_type_number, levels = cell_type_number_levels)




pbmc$cell_type_legend <- paste(pbmc$cell_type_number, pbmc$Kat_cell_type, sep = ": ")




sort(unique(pbmc$cell_type_legend))

cell_type_legend_levels <- c(
  
"1: Adventitial cell",             
"2: B cell",                       
"3: Basophil",                    
"4: Dendritic cell",               
"5: Endothelial cell",             
"6: Epithelial cell",             
"7: Fibroblast",                   
"8: Macrophage",                   
"9: Mature NK T cell", 
"10: Mesothelial cell",            
"11: Monocyte",                   
"12: Myofibroblast cell",          
"13: Neutrophil",                  
"14: Pericyte",                   
"15: Plasma cell",                 
"16: Plasmacytoid dendritic cell", 
"17: Smooth muscle cell",         
"18: T cell" 
  
)


pbmc$cell_type_legend <- factor(pbmc$cell_type_legend, levels = cell_type_legend_levels)




pbmc$Kat_cell_type_plot <- as.character(pbmc$Kat_cell_type)

sort(unique(pbmc$Kat_cell_type_plot))

               
pbmc$Kat_cell_type_plot[pbmc$Kat_cell_type_plot == "Plasmacytoid dendritic cell"] <- "pDC" 
pbmc$Kat_cell_type_plot[pbmc$Kat_cell_type_plot == "Smooth muscle cell"] <- "SM cell"           




}



# large intestine

if (sample_name == "TS_large_intestine") {
  
  unique(pbmc$cell_type)
  
  
  Kat_T_cell_group <- c("CD4-positive, alpha-beta T cell",
                        "CD8-positive, alpha-beta T cell")
  
  
  pbmc$Kat_cell_type_1 <- as.character(pbmc$cell_type)
  
  pbmc$Kat_cell_type_1 [pbmc$Kat_cell_type_1 %in% Kat_T_cell_group] <- "T cell"
  
  pbmc$Kat_cell_type_1 [pbmc$Kat_cell_type_1  == "large intestine goblet cell"] <- "goblet cell"
  pbmc$Kat_cell_type_1 [pbmc$Kat_cell_type_1  == "intestinal crypt stem cell of large intestine"] <- "crypt stem cell"
  pbmc$Kat_cell_type_1 [pbmc$Kat_cell_type_1  == "enterocyte of epithelium of large intestine"] <- "enterocyte of epithelium"
  pbmc$Kat_cell_type_1 [pbmc$Kat_cell_type_1  == "paneth cell of colon"] <- "paneth cell"
  pbmc$Kat_cell_type_1 [pbmc$Kat_cell_type_1  == "transit amplifying cell of colon"] <- "transit amplifying cell"
  pbmc$Kat_cell_type_1 [pbmc$Kat_cell_type_1  == "intestinal enteroendocrine cell"] <- "enteroendocrine cell"
  pbmc$Kat_cell_type_1 [pbmc$Kat_cell_type_1  == "intestinal tuft cell"] <- "tuft cell"
  
  
  
  pbmc$Kat_cell_type <- as.character(pbmc$cell_type)
  
  pbmc$Kat_cell_type[pbmc$Kat_cell_type %in% Kat_T_cell_group] <- "T cell"
  
  pbmc$Kat_cell_type[pbmc$Kat_cell_type == "large intestine goblet cell"] <- "Goblet cell"
  pbmc$Kat_cell_type[pbmc$Kat_cell_type == "intestinal crypt stem cell of large intestine"] <- "Crypt stem cell"
  pbmc$Kat_cell_type[pbmc$Kat_cell_type == "enterocyte of epithelium of large intestine"] <- "Enterocyte of epithelium"
  pbmc$Kat_cell_type[pbmc$Kat_cell_type == "paneth cell of colon"] <- "Paneth cell"
  pbmc$Kat_cell_type[pbmc$Kat_cell_type == "transit amplifying cell of colon"] <- "Transit amplifying cell"
  pbmc$Kat_cell_type[pbmc$Kat_cell_type == "intestinal enteroendocrine cell"] <- "Enteroendocrine cell"
  pbmc$Kat_cell_type[pbmc$Kat_cell_type == "intestinal tuft cell"] <- "Tuft cell"

pbmc$Kat_cell_type[pbmc$Kat_cell_type == "fibroblast"] <- "Fibroblast"                
pbmc$Kat_cell_type[pbmc$Kat_cell_type == "gut endothelial cell"] <- "Gut endothelial cell"       
pbmc$Kat_cell_type[pbmc$Kat_cell_type == "mast cell"] <- "Mast cell"                 
pbmc$Kat_cell_type[pbmc$Kat_cell_type == "monocyte"] <- "Monocyte"                   
pbmc$Kat_cell_type[pbmc$Kat_cell_type == "neutrophil"] <- "Neutrophil"                 
pbmc$Kat_cell_type[pbmc$Kat_cell_type == "plasma cell"] <- "Plasma cell"               


  
sort(unique(pbmc$Kat_cell_type))
  

pbmc$cell_type_number <- "error"

pbmc$cell_type_number[pbmc$Kat_cell_type == "B cell"] <- "1"                   
pbmc$cell_type_number[pbmc$Kat_cell_type == "Crypt stem cell"] <- "2"          
pbmc$cell_type_number[pbmc$Kat_cell_type == "Enterocyte of epithelium"] <- "3" 
pbmc$cell_type_number[pbmc$Kat_cell_type == "Enteroendocrine cell"] <- "4"    
pbmc$cell_type_number[pbmc$Kat_cell_type == "Fibroblast"] <- "5"               
pbmc$cell_type_number[pbmc$Kat_cell_type == "Goblet cell"] <- "6"              
pbmc$cell_type_number[pbmc$Kat_cell_type == "Gut endothelial cell"] <- "7"     
pbmc$cell_type_number[pbmc$Kat_cell_type == "Mast cell"] <- "8"               
pbmc$cell_type_number[pbmc$Kat_cell_type == "Monocyte"] <- "9"                 
pbmc$cell_type_number[pbmc$Kat_cell_type == "Neutrophil"] <- "10"               
pbmc$cell_type_number[pbmc$Kat_cell_type == "Paneth cell"] <- "11"              
pbmc$cell_type_number[pbmc$Kat_cell_type == "Plasma cell"] <- "12"             
pbmc$cell_type_number[pbmc$Kat_cell_type == "T cell"] <- "13"                   
pbmc$cell_type_number[pbmc$Kat_cell_type == "Transit amplifying cell"] <- "14"  
pbmc$cell_type_number[pbmc$Kat_cell_type == "Tuft cell"] <- "15"

cell_type_number_levels <- as.character(c(1:15))

pbmc$cell_type_number <- factor(pbmc$cell_type_number, levels = cell_type_number_levels)


pbmc$cell_type_legend <- paste(pbmc$cell_type_number, pbmc$Kat_cell_type, sep = ": ")

sort(unique(pbmc$cell_type_legend))

cell_type_legend_levels <- c(
  
"1: B cell",                   
"2: Crypt stem cell",         
"3: Enterocyte of epithelium", 
"4: Enteroendocrine cell",     
"5: Fibroblast",               
"6: Goblet cell",             
"7: Gut endothelial cell",     
"8: Mast cell",                
"9: Monocyte",  
"10: Neutrophil",              
"11: Paneth cell",             
"12: Plasma cell",            
"13: T cell",                  
"14: Transit amplifying cell", 
"15: Tuft cell"   
  
)


pbmc$cell_type_legend <- factor(pbmc$cell_type_legend, levels = cell_type_legend_levels)





pbmc$Kat_cell_type_plot <- as.character(pbmc$Kat_cell_type)

sort(unique(pbmc$Kat_cell_type_plot))



  
}
















tissue_list <- unique(pbmc$tissue)

compartment_list <- unique(pbmc$compartment)




# tissue_df <- as.data.frame(tissue_list)
# tissue_df$compartment <- cell_group_name
# names(tissue_df)[names(tissue_df) == "tissue_list"] <- "tissue"
# tissue_df <- tissue_df[ , c("compartment", "tissue")]
# 
# tissue_df$tissue <- as.character(tissue_df$tissue)
# tissue_df <- tissue_df[order(tissue_df$tissue), ]
# 
# 
# write.table(tissue_df, file = paste0(gene_type, "_tissue_list.txt"), quote = FALSE, sep = "\t", row.names = FALSE)










Idents(pbmc) <- "cell_type"



pdf(paste0(sample_name, "_UMAP_plot_cell-type_CL.pdf"), width = 9, height = 7, useDingbats = FALSE)
p <- DimPlot(pbmc, reduction = "umap", label = TRUE, group.by = "cell_type", pt.size = 0.5, label.size = 2.7) +
  theme(legend.text = element_text(size = 7)) +
  theme(plot.subtitle = element_text(hjust = 0.5)) +
  # ggtitle(sample_name) +
  theme(plot.subtitle = element_text(hjust = 0.5)) +
  labs(subtitle = paste0("Cell Ontology", " - ", paste0(tissue_name, " (", organism_name, ")"))) +
  labs(caption = dataset_name)
print(p)
dev.off()


pdf(paste0(sample_name, "_UMAP_plot_cell-type_CL_no-legend.pdf"), width = 9, height = 7, useDingbats = FALSE)
p <- DimPlot(pbmc, reduction = "umap", label = TRUE, group.by = "cell_type", pt.size = 0.5) +
                          NoLegend() +
  theme(plot.subtitle = element_text(hjust = 0.5)) +
  # ggtitle(sample_name) +
  labs(subtitle = paste0("Cell Ontology", " - ", paste0(tissue_name, " (", organism_name, ")"))) +
  labs(caption = dataset_name)
print(p)
dev.off()



pdf(paste0(sample_name, "_UMAP_plot_cell-type_FA.pdf"), width = 9, height = 7, useDingbats = FALSE)
p <- DimPlot(pbmc, reduction = "umap", label = TRUE, group.by = "free_annotation", pt.size = 0.5, label.size = 2.7) +
  theme(legend.text = element_text(size = 7)) +
  theme(plot.subtitle = element_text(hjust = 0.5)) +
  # ggtitle(sample_name) +
  labs(subtitle = paste0("free annotation", " - ", paste0(tissue_name, " (", organism_name, ")"))) +
  labs(caption = dataset_name)
print(p)
dev.off()

pdf(paste0(sample_name, "_UMAP_plot_cell-type_FA_no-legend.pdf"), width = 9, height = 7, useDingbats = FALSE)
p <- DimPlot(pbmc, reduction = "umap", label = TRUE, group.by = "free_annotation", pt.size = 0.5) +
  NoLegend() +
  theme(plot.subtitle = element_text(hjust = 0.5)) +
  # ggtitle(sample_name) +
  labs(subtitle = paste0("free annotation", " - ", paste0(tissue_name, " (", organism_name, ")"))) +
  labs(caption = dataset_name)
print(p)
dev.off()








pdf(paste0(sample_name, "_UMAP_plot_cell-type_Kat.pdf"), width = 9, height = 7, useDingbats = FALSE)
p <- DimPlot(pbmc, reduction = "umap", label = TRUE, group.by = "Kat_cell_type", pt.size = 0.5, label.size = 2.7) +
  theme(legend.text = element_text(size = 7)) +
  theme(plot.subtitle = element_text(hjust = 0.5)) +
  # ggtitle(sample_name) +
  theme(plot.subtitle = element_text(hjust = 0.5)) +
  labs(subtitle = paste0("Cell Type", " - ", paste0(tissue_name, " (", organism_name, ")"))) +
  labs(caption = dataset_name)
print(p)
dev.off()









pdf(paste0(sample_name, "_UMAP_plot_cell-type_number_fig.pdf"), width = 9, height = 7, useDingbats = FALSE)
p <- DimPlot(pbmc, reduction = "umap", label = TRUE, group.by = "cell_type_number", pt.size = 0.5, label.size = 10) +
  
  ggtitle(plot_title) +
  # theme(legend.text = element_text(size = 7)) +
  # theme(plot.subtitle = element_text(hjust = 0.5)) +
  # labs(subtitle = paste0("Cell Ontology", " - ", paste0(tissue_name, " (", organism_name, ")"))) +
  # labs(caption = dataset_name)
  theme(legend.position = "none",
          plot.title = element_text(size = 24, face = "plain"),
        axis.ticks.length = unit(0.1, "inch"),
        axis.title = element_text(size = 22, face = "plain"),
        axis.text = element_text(size = 22),
        
        )




xlim_all_samples <- layer_scales(p)$x$range$range
ylim_all_samples <- layer_scales(p)$y$range$range


print(p)
dev.off()




# https://github.com/satijalab/seurat/issues/3899


pdf(paste0(sample_name, "_UMAP_plot_cell-type_legend_fig.pdf"), width = 9, height = 7, useDingbats = FALSE)
p <- DimPlot(pbmc, reduction = "umap", label = TRUE, group.by = "cell_type_legend", pt.size = 0.5, label.size = 2.7) +
  theme(legend.text = element_text(size = 7)) +
  theme(plot.subtitle = element_text(hjust = 0.5)) +
  labs(subtitle = paste0("Cell Ontology", " - ", paste0(tissue_name, " (", organism_name, ")"))) +
  labs(caption = dataset_name) +
  theme(
    legend.text = element_text(size = 18, margin = margin(l = 1, unit = "pt")),
    legend.box.spacing = unit(0, "pt"),
    
    # legend.key.height= unit(2, "cm"),
    
    # legend.key.size = unit(0.5, "cm"),
    
    legend.key.spacing.y = unit(0.1, "cm")

  ) +
  
  guides(color = guide_legend(override.aes = list(size = 10)))




print(p)
dev.off()






pdf(paste0(sample_name, "_UMAP_plot_individual-cell-types.pdf"), width = 9, height = 7, useDingbats = FALSE)

cell_type_names <- sort(unique(pbmc$cell_type_legend))
sample_colors <- hue_pal()(length(cell_type_names))


Idents(object = pbmc) <- "cell_type_legend"

for (current_cell_type in cell_type_names) {
  
cluster_position <- which(cell_type_names == current_cell_type)

current_cells <- WhichCells(object = pbmc, idents = current_cell_type)

p <- DimPlot(pbmc, reduction = "umap", label = FALSE, group.by = "cell_type_legend", cells = current_cells, pt.size = 0.5, cols = sample_colors[cluster_position]) +
                          xlim(xlim_all_samples[1], xlim_all_samples[2]) +
                          ylim(ylim_all_samples[1], ylim_all_samples[2]) +
                          # NoLegend() +
                          ggtitle(plot_title)
print(p)
  
}

dev.off()











pdf(paste0(sample_name, "_UMAP_plot_compartment.pdf"), width = 9, height = 7, useDingbats = FALSE)
p <- DimPlot(pbmc, reduction = "umap", label = TRUE, group.by = "compartment", pt.size = 0.5, label.size = 2.7) +
  theme(legend.text = element_text(size = 7)) +
  ggtitle(sample_name) +
  labs(caption = dataset_name)
print(p)
dev.off()

pdf(paste0(sample_name, "_UMAP_plot_compartment_no-legend.pdf"), width = 9, height = 7, useDingbats = FALSE)
p <- DimPlot(pbmc, reduction = "umap", label = TRUE, group.by = "compartment", pt.size = 0.5) +
                          NoLegend() +
 ggtitle(sample_name) +
 labs(caption = dataset_name)
print(p)
dev.off()





  
found_genes <- c()
found_genes <- gene_list_hgnc$ensembl_gene_id[which(gene_list_hgnc$ensembl_gene_id %in% pbmc@assays$RNA$data@Dimnames[[1]])]
  
missing_genes_df <- gene_list_hgnc[which(!gene_list_hgnc$ensembl_gene_id %in% found_genes), ]

write.table(missing_genes_df, file = paste0(sample_name, "_missing_genes.txt"), quote = FALSE, sep = "\t", row.names = FALSE)







summary_df <- data.frame()

compartment_cell_types_df <- data.frame()

compartment_list <- "all"


for (current_compartment in compartment_list) {

  message(current_compartment)
  
  Idents(object = pbmc) <- "compartment"
  
  # current_cells <- WhichCells(object = pbmc, idents = current_tissue)

  
  
  # ***** need to change if separating by compartments
  # pbmc.temp <- subset(pbmc, subset = compartment == current_compartment)
  pbmc.temp <- pbmc
  
  
  
  # disease_list <- unique(lung.temp_dataset$disease2)
  
  # rm(lung.temp_dataset)
  # gc()
  
  
  current_compartment_cell_types_df <- as.data.frame(as.character(unique(pbmc.temp$cell_type)))
  
  names(current_compartment_cell_types_df)[1] <- "cell_type"
  
  current_compartment_cell_types_df$compartment <- current_compartment
  current_compartment_cell_types_df$tissue <- tissue_name
  
  current_compartment_cell_types_df <- current_compartment_cell_types_df[ c("compartment", "tissue", "cell_type")]
  
  current_compartment_cell_types_df <- current_compartment_cell_types_df[order(current_compartment_cell_types_df$cell_type), ]
  
  compartment_cell_types_df <- rbind(compartment_cell_types_df, current_compartment_cell_types_df)
  
  
  
  

  
  
  found_genes <- c()
  found_genes <- gene_list_hgnc$ensembl_gene_id[which(gene_list_hgnc$ensembl_gene_id %in% pbmc.temp@assays$RNA@data@Dimnames[[1]])]
    found_genes_df <- gene_list_hgnc[which(gene_list_hgnc$ensembl_gene_id %in% found_genes), ]

found_genes_convert <- found_genes_df[ , c("gene", "ensembl_gene_id")]
row.names(found_genes_convert) <- found_genes_convert$ensembl_gene_id
found_genes_convert$ensembl_gene_id <- NULL
found_genes_convert <- as.vector(found_genes_convert)

found_genes_convert <- setNames(found_genes_df$gene, found_genes_df$ensembl_gene_id)






cat("Making dot plots . . .", "\n\n")





pdf(paste(sample_name, gsub(" ", "_", current_compartment), gene_type, "DotPlot_no-scale_x-axis-genes.pdf", sep = "_"), width = 15, height = 10)

p <- DotPlot(pbmc.temp, features = found_genes, group.by = "cell_type", scale = FALSE) + RotatedAxis() +
  ggtitle(paste0(current_compartment, " (", tissue_name, ")")) +
  scale_x_discrete(labels = found_genes_convert)
  
print(p)
dev.off()



pdf(paste(sample_name, gsub(" ", "_", current_compartment), gene_type, "DotPlot_scaled_x-axis-genes.pdf", sep = "_"), width = 15, height = 7)

p <- DotPlot(pbmc.temp, features = found_genes, group.by = "cell_type", scale = TRUE) + RotatedAxis() +
  ggtitle(paste0(current_compartment, " (", tissue_name, ")")) +
  scale_x_discrete(labels = found_genes_convert)
  
print(p)
dev.off()

  


found_genes_df_paper <- found_genes_df[found_genes_df$gene %in% plot_genes, ]

found_genes_paper <- found_genes_df_paper$ensembl_gene_id


found_genes_convert_paper <- found_genes_convert[found_genes_paper]




# dot plots for paper

# https://rpubs.com/eraz0001/enhanced_dotplot

# https://github.com/satijalab/seurat/issues/3914



pdf(paste(sample_name, gsub(" ", "_", current_compartment), gene_type, "DotPlot_FA.pdf", sep = "_"), width = 9, height = 7)

# p <- DotPlot(pbmc.temp, features = found_genes_paper, group.by = "free_annotation", scale = TRUE, dot.scale = 1) + coord_flip() +
p <- DotPlot(pbmc.temp, features = found_genes_paper, group.by = "free_annotation", scale = TRUE) + coord_flip() +
  ggtitle(paste0(current_compartment, " (", tissue_name, ")")) +
  xlab("Gene") + 
  ylab("Cell Type (free annotation)") + 
  labs(caption = dataset_name) +
  scale_x_discrete(labels = found_genes_convert_paper) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

print(p)
dev.off()

dotplot_data <- p$data
write.csv(dotplot_data, file = paste(sample_name, gsub(" ", "_", current_compartment), gene_type, "DotPlot_FA_data.csv", sep = "_"), row.names = FALSE)






pdf(paste(sample_name, gsub(" ", "_", current_compartment), gene_type, "DotPlot_CL.pdf", sep = "_"), width = 9, height = 7)

p <- DotPlot(pbmc.temp, features = found_genes_paper, group.by = "cell_type", scale = TRUE) + coord_flip() +
  ggtitle(paste0(current_compartment, " (", tissue_name, ")")) +
  xlab("Gene") + 
  ylab("Cell Type (Cell Ontology)") + 
  labs(caption = dataset_name) +
  scale_x_discrete(labels = found_genes_convert_paper) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

print(p)
dev.off()


dotplot_data <- p$data
write.csv(dotplot_data, file = paste(sample_name, gsub(" ", "_", current_compartment), gene_type, "DotPlot_CL_data.csv", sep = "_"), row.names = FALSE)






pdf(paste(sample_name, gsub(" ", "_", current_compartment), gene_type, "DotPlot_Kat-cell-type.pdf", sep = "_"), width = 11, height = 7)

p <- DotPlot(pbmc.temp, features = found_genes_paper, group.by = "Kat_cell_type", scale = TRUE) + coord_flip() +
  ggtitle(paste0(current_compartment, " (", tissue_name, ")")) +
  xlab("Gene") + 
  ylab("Cell Type (Kat)") + 
  labs(caption = dataset_name) +
  scale_x_discrete(labels = found_genes_convert_paper) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

print(p)
dev.off()


dotplot_data <- p$data
write.csv(dotplot_data, file = paste(sample_name, gsub(" ", "_", current_compartment), gene_type, "DotPlot_Kat-cell-type_data.csv", sep = "_"), row.names = FALSE)







pdf(paste(sample_name, gsub(" ", "_", current_compartment), gene_type, "DotPlot_Kat-cell-type_fig.pdf", sep = "_"), width = 12, height = 7)

p <- DotPlot(pbmc.temp, features = found_genes_paper, group.by = "Kat_cell_type_plot", scale = TRUE) + coord_flip() +
  ggtitle(plot_title) +
  scale_x_discrete(labels = found_genes_convert_paper) +
  
    theme(axis.text.x = element_text(size = 21, angle = 45, hjust = 0.95, vjust = 1),
        axis.text.y = element_text(size = 21, face = "italic"),
        axis.title = element_blank(),     
        plot.title = element_text(size = 24, face = "plain", hjust = 0.5, vjust = 0),
      axis.ticks.length = unit(0.075, "inch"),
        # axis.title = element_text(size = 16, face = "plain"),
      legend.text = element_text(size = 20),
      legend.title = element_text(size = 20)


        )

print(p)
dev.off()


dotplot_data <- p$data
write.csv(dotplot_data, file = paste(sample_name, gsub(" ", "_", current_compartment), gene_type, "DotPlot_Kat-cell-type_data_fig.csv", sep = "_"), row.names = FALSE)








cat("Making violin plots . . .", "\n\n")


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
  


  # raster_pdf(paste0(gene_type, "_", cell_group_name, "_", current_tissue, "_cell-type_violin_plots.pdf"), res = 300)

  pdf(paste0(sample_name, "_", gsub(" ", "_", current_compartment), "_", gene_type, "_violin_plots_", anno_short, "_fig.pdf"), useDingbats = FALSE, width = 7, height = 6)
  
    
  # pdf(paste0("NLR-genes_", current_dataset, "_", gsub("/", "_", current_cell_type), "_violin_plots.pdf"), useDingbats = FALSE)
  
  
  for (current_gene in found_genes) {
    
    Idents(object = pbmc) <- "compartment"
    
    
    matching_gene <- gene_list_hgnc[gene_list_hgnc$ensembl_gene_id == current_gene, "gene"]
    
    
    vln_plot_title = paste0(matching_gene, " (", current_gene, ")\n", paste0(current_compartment, " (", tissue_name, ")"))
    
   
    # message(current_gene)
    
     message(paste0(matching_gene, " (", current_gene, ")"))
    
    # message(vln_plot_title)
    
    
    
    all_genes <- pbmc.temp@assays$RNA@data@Dimnames[[1]]
    
    gene_test <- all_genes[which(all_genes == current_gene)]
    
    if (length(gene_test) != 1) {
      
      stop("current gene length is not 1")
    }
    
    
    
    # p <- VlnPlot(object = pbmc.temp, features = current_gene, group.by = "cell_type", pt.size = 1, raster = TRUE) + 
      p <- VlnPlot(object = pbmc.temp, features = current_gene, group.by = current_annotation_type, pt.size = 1, raster = TRUE) + 
      ggtitle(vln_plot_title) +
      # xlab("Cell Type (free annotation") + 
      xlab(paste0("Cell Type (", anno_long, ")")) +
      ylab("Expression Level") + 
      labs(caption = dataset_name) +
      theme(legend.position = "none") +
      theme(axis.text.x = element_text(size = rel(0.8)))
    
    p$layers[[2]]$aes_params$alpha <- 0.2
    
    print(p)
    

    
    
    
      p <- VlnPlot(object = pbmc.temp, features = current_gene, group.by = "Kat_cell_type_plot", pt.size = 1, raster = TRUE) +

      ggtitle(plot_title) +
      
      labs(subtitle = matching_gene) +

      # ylab(paste(plot_title, "Expression Level", sep = "\n")) + 
      ylab("Expression Level") + 

      theme(legend.position = "none", 
         axis.title.x = element_blank(),
      plot.title = element_text(size = 20, face = "plain", vjust = 0),
      axis.ticks.length = unit(0.1, "inch"),
        axis.title.y = element_text(size = 16, face = "plain"),
        axis.text.x = element_text(size = 16),
      axis.text.y = element_text(size = 16),
      plot.subtitle = element_text(size = 20, face = "italic", hjust = 0.5, vjust = 0)
      
      ) 
    
    p$layers[[2]]$aes_params$alpha <- 0.2
    
    print(p)
    
    
    
    
    
    
    
    
    
  
    
    # cell_type_list <- unique(pbmc.temp$cell_type)
    # Idents(pbmc.temp) <- "cell_type"
    
    
    # ***** need to change for cell type (CL or FA) used *****
    # cell_type_list <- unique(pbmc.temp$free_annotation)
    cell_type_list <- unique(pbmc.temp@meta.data[ , current_annotation_type])
    
    
   
    
    
    # Idents(pbmc.temp) <- "free_annotation"
    Idents(pbmc.temp) <- current_annotation_type
    
  
  for (current_cell_type in cell_type_list) {
  
    
    current_cells <- WhichCells(object = pbmc.temp, idents = current_cell_type)
  
  
            expression_values <- pbmc.temp@assays$RNA@data[current_gene, current_cells]

          if (!is.vector(expression_values)) {
               
               stop("expression values are not a vector")
          }
          

          
          current_summary <- data.frame(gene = matching_gene)
          
          current_summary$ensembl_id <- current_gene
          
          current_summary$tissue <- tissue_name
          current_summary$compartment <- current_compartment
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
  
  
  
  
  
  
  
  
  
  
  
  
  
  rm(pbmc.temp)
  gc()
  
  

} # END current compartment



  write.csv(summary_df, file = paste(sample_name, gene_type, "compartment_summary.csv", sep = "_"), row.names = FALSE)


  compartment_cell_types_df <- compartment_cell_types_df[order(compartment_cell_types_df$compartment), ]
  
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



