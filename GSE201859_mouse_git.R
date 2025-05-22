



library(ggplot2)
library(ggrepel)


library(rasterpdf, lib.loc = "/nas/longleaf/home/wrobel/R/local_R_libraries/R-4.2.2")


library(scales)


# for md5sum
library(tools)


# For a (much!) faster implementation of the Wilcoxon Rank Sum Test,
# (default method for FindMarkers) please install the presto package
library(presto, lib.loc = "/nas/longleaf/home/wrobel/R/local_R_libraries/R-4.3.1")



# location of Seurat v5
library(Seurat, lib.loc = "/nas/longleaf/home/wrobel/R/local_R_libraries/R-4.3.1")





script_start_time <- Sys.time()




# *** START user input ***

split_report <- TRUE




sample_name <- "GSE201859_mouse_si"
dataset_name <- "GSE201859_mouse_si"
plot_title <- "Murine Small Intestine"
organism_name <- "mouse"
tissue_name <- "small intestine"
compartment_name <- "small intestine"
rds_input_filename <- "/work/users/w/r/wrobel/tinglab/scRNA-Seq/GSE201859/rds/mouse/GSE201859_mouse_si_v5.rds"






# gene_type <- "Nlr-genes"
gene_type <- "Kat-genes"


# gene_list_file <- "/proj/tinglab/users/wrobel/NLR_gene_list/Nlr_gene_list_mouse.txt"


additional_genes <- c(
  
  "Actb",
  # "Myd88",
  # "Tlr2",
  # "Tlr6",
  # "Gsdma",
  # "Gsdmb",
  # "Gsdbc",
  # "Gsdmd",
  # "Gsdme",

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








# gene_list_df <- read.table(gene_list_file, sep = "\t", header = TRUE, stringsAsFactors = FALSE)
# 
# plot_genes <- gene_list_df$gene

paper_genes <- c("Meis2", "Plb1", "Fabp6")

additional_genes <- c(plot_genes, "Actb", paper_genes)



cell_group_name <- "test"



# mgi_file <- "/proj/tinglab/users/wrobel/R/common_input/mgi/MRK_ENSEMBL_2024-05-01.rpt"
mgi_file <- "/proj/tinglab/users/wrobel/gene_datasets/mgi/MRK_ENSEMBL/2025-03-31/MRK_ENSEMBL_download_2025-03-31.rpt"



padj_value_cutoff <- 0.05


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

gene_list <- gene_list_df$gene

# found_genes <- gene_list[which(gene_list %in% lung@assays$RNA@counts@Dimnames[[1]])]

















mgi_filename.md5 <- md5sum(mgi_file)

cat("MGI gene annotation file:", mgi_file, "\n")
cat("md5:", mgi_filename.md5, "\n")
cat("\n")
cat("\n")




mgi_all <- read.table(mgi_file, sep = "\t", header = FALSE, quote = "", comment.char = "", stringsAsFactors = FALSE)


mgi_header <- c(

"MGI_Marker_Accession_ID",
"Marker_Symbol",	
"Marker_Name",
"cM_Position",	
"Chromosome",	
"Ensembl_Accession_ID",	
"Ensembl_Transcript_ID",	
"Ensembl_Protein_ID",	
"Feature_Types",	
"Genome_Coordinate_Start",	
"Genome_Coordinate_End",	
"Strand",	
"BioTypes"

)

colnames(mgi_all) <- mgi_header




# select_mgi_columns <- c(
#   
#   "MGI_Marker_Accession_ID", 
#   "Marker_Symbol",           
#   "Status",                  
#   "Marker_Name",
#   "Chromosome",              
#   "Type", 
#   "Feature_Types",           
#   "Genome_Coordinate_Start",
#   "Genome_Coordinate_End",   
#   "Strand",                 
#   "BioTypes" 
# )


# mgi <- mgi_all[ , select_mgi_columns]
mgi <- mgi_all

mgi$gene <- mgi$Marker_Symbol
# mgi$Marker_Symbol <- NULL

names(mgi)[names(mgi) == "Ensembl_Accession_ID"] <- "ensembl_gene_id"

mgi_gene_lookup <- mgi[ , c("ensembl_gene_id", "gene")]




# mgi <- mgi[mgi$Status == "O", ]
# mgi <- mgi[mgi$Type %in% c("Gene", "Pseudogene"), ]




gene_list_mgi <- merge(x = gene_list_df, y = mgi, by.x = "gene", by.y = "gene", all.x = TRUE, sort = FALSE)


write.table(gene_list_mgi, file = paste0(sample_name, "_original_gene_list.txt"), quote = FALSE, sep = "\t", row.names = FALSE)









rds_input_filename.md5 <- md5sum(rds_input_filename)

cat("RDS input file:", rds_input_filename, "\n")
cat("md5:", rds_input_filename.md5, "\n")
cat("\n")
cat("\n")



pbmc <- readRDS(rds_input_filename)

original_metadata <- pbmc@meta.data
original_metadata$observation_joinid <- NULL

write.table(original_metadata, file = paste0(sample_name, "_original_metadata.txt"), quote = FALSE, sep = "\t", row.names = FALSE)







# ***** subset by seg group
# pbmc <- subset(pbmc, subset = seg_group == seg_group_name)

# unique(pbmc$seg_group)
# sort(unique(pbmc$segment))




# as.data.frame(table(pbmc$seg_domain, pbmc$seg_classic))

# ***** need to changed based on method *****


# for smart-seq2

# unique(pbmc$tissue)
# unique(pbmc$tissue_original)
# 
# table(pbmc@meta.data[ , c("tissue", "tissue_original")])
# 
# compare_tissue <- as.data.frame(table(pbmc@meta.data[ , c("tissue", "tissue_original")]))
# write.table(compare_tissue, file = paste0(sample_name, "_compare_tissue.txt"), quote = FALSE, sep = "\t", row.names = FALSE)


# for 10x


unique(pbmc$tissue)


table(pbmc@meta.data[ , c("disease", "organism",	"sex",	"tissue",	"development_stage")])

meta_groups <- as.data.frame(table(pbmc@meta.data[ , c("disease", "organism",	"sex",	"tissue",	"development_stage")]))
write.table(meta_groups, file = paste0(sample_name, "_meta_groups_table.txt"), quote = FALSE, sep = "\t", row.names = FALSE)



cell_type_versions <- as.data.frame(table(pbmc@meta.data[ , c("cell.type", "cell_type")]))
write.table(cell_type_versions, file = paste0(sample_name, "_cell_types_versions_table.txt"), quote = FALSE, sep = "\t", row.names = FALSE)






pbmc$cell_type_original <- pbmc$cell.type
# pbmc$cell_type <- pbmc$CellType




unique(pbmc@meta.data$cell.type)
unique(pbmc@meta.data$cell_type)


unique(pbmc@meta.data$tissue)


# pbmc <- subset(x = pbmc, subset = disease == "normal")
# pbmc <- subset(x = pbmc, subset = tissue == tissue_name)

# gc()


filtered_metadata <- pbmc@meta.data
filtered_metadata$observation_joinid <- NULL

write.table(filtered_metadata, file = paste0(sample_name, "_filtered_metadata.txt"), quote = FALSE, sep = "\t", row.names = FALSE)




# pbmc$tissue2 <- as.character(pbmc$tissue)
# pbmc$cell_type2 <- as.character(pbmc$cell_type)
# 
# cell_type_df <- as.data.frame(table(pbmc@meta.data[ , c("tissue2", "cell_type2")]))
# 
# write.table(cell_type_df, file = paste0(sample_name, "_cell_type_table.txt"), quote = FALSE, sep = "\t", row.names = FALSE)









Idents(object = pbmc) <- "tissue"






pbmc$tissue2 <- as.character(pbmc$tissue)
# pbmc$compartment2 <- as.character(pbmc$compartment)
pbmc$cell_type2 <- as.character(pbmc$cell_type)

cell_type_df <- as.data.frame(table(pbmc@meta.data[ , c("tissue2", "cell_type2")]))

cell_type_df <- cell_type_df[cell_type_df$Freq != 0, ]

write.table(cell_type_df, file = paste0(sample_name, "_cell_type_table_CL.txt"), quote = FALSE, sep = "\t", row.names = FALSE)




# free annotation version
pbmc$free_annotation2 <- as.character(pbmc$cell.type)

cell_type_fa_df <- as.data.frame(table(pbmc@meta.data[ , c("tissue2", "free_annotation2")]))

cell_type_fa_df <- cell_type_fa_df[cell_type_fa_df$Freq != 0, ]

write.table(cell_type_fa_df, file = paste0(sample_name, "_cell_type_table_FA.txt"), quote = FALSE, sep = "\t", row.names = FALSE)


unique(pbmc$cell_type)


pbmc$Kat_cell_type_1 <- as.character(pbmc$cell_type)

pbmc$Kat_cell_type_1[pbmc$Kat_cell_type_1 == "small intestine goblet cell"] <- "goblet cell"
pbmc$Kat_cell_type_1[pbmc$Kat_cell_type_1 == "tuft cell of small intestine"] <- "tuft cell"
pbmc$Kat_cell_type_1[pbmc$Kat_cell_type_1 == "enterocyte of epithelium of small intestine"] <- "enterocyte of epithelium"
pbmc$Kat_cell_type_1[pbmc$Kat_cell_type_1 == "transit amplifying cell of small intestine"] <- "transit amplifying cell"
pbmc$Kat_cell_type_1[pbmc$Kat_cell_type_1 == "enteroendocrine cell of small intestine"] <- "enteroendocrine cell"
pbmc$Kat_cell_type_1[pbmc$Kat_cell_type_1 == "paneth cell of epithelium of small intestine"] <- "paneth cell of epithelium"

unique(pbmc$Kat_cell_type_1)



pbmc$Kat_cell_type <- as.character(pbmc$cell_type)

pbmc$Kat_cell_type[pbmc$Kat_cell_type == "small intestine goblet cell"] <- "Goblet cell"
pbmc$Kat_cell_type[pbmc$Kat_cell_type == "tuft cell of small intestine"] <- "Tuft cell"
pbmc$Kat_cell_type[pbmc$Kat_cell_type == "enterocyte of epithelium of small intestine"] <- "Enterocyte"
pbmc$Kat_cell_type[pbmc$Kat_cell_type == "transit amplifying cell of small intestine"] <- "Transit amplifying cell"
pbmc$Kat_cell_type[pbmc$Kat_cell_type == "enteroendocrine cell of small intestine"] <- "Enteroendocrine cell"
pbmc$Kat_cell_type[pbmc$Kat_cell_type == "paneth cell of epithelium of small intestine"] <- "Paneth cell of epithelium"
pbmc$Kat_cell_type[pbmc$Kat_cell_type == "stem cell"] <- "Stem cell"

sort(unique(pbmc$Kat_cell_type))


pbmc$cell_type_number <- "error"

pbmc$cell_type_number[pbmc$Kat_cell_type == "Enterocyte"] <- "1"                
pbmc$cell_type_number[pbmc$Kat_cell_type == "Enteroendocrine cell"] <- "2"      
pbmc$cell_type_number[pbmc$Kat_cell_type == "Goblet cell"] <- "3"               
pbmc$cell_type_number[pbmc$Kat_cell_type == "Paneth cell of epithelium"] <- "4" 
pbmc$cell_type_number[pbmc$Kat_cell_type == "Stem cell"] <- "5"                
pbmc$cell_type_number[pbmc$Kat_cell_type == "Transit amplifying cell"] <- "6"   
pbmc$cell_type_number[pbmc$Kat_cell_type == "Tuft cell"] <- "7"   



# pbmc$cell_type_legend <- "error"
# 
# pbmc$cell_type_legend[pbmc$Kat_cell_type == "Enterocyte"] <- "1: Enterocyte"                
# pbmc$cell_type_legend[pbmc$Kat_cell_type == "Enteroendocrine cell"] <- "2: Enteroendocrine cell"      
# pbmc$cell_type_legend[pbmc$Kat_cell_type == "Goblet cell"] <- "3: Goblet cell"               
# pbmc$cell_type_legend[pbmc$Kat_cell_type == "Paneth cell of epithelium"] <- "4: Paneth cell of epithelium" 
# pbmc$cell_type_legend[pbmc$Kat_cell_type == "Stem cell"] <- "5: Stem cell"                
# pbmc$cell_type_legend[pbmc$Kat_cell_type == "Transit amplifying cell"] <- "6: Transit amplifying cell"   
# pbmc$cell_type_legend[pbmc$Kat_cell_type == "Tuft cell"] <- "7: Tuft cell"   



pbmc$cell_type_legend <- paste(pbmc$cell_type_number, pbmc$Kat_cell_type, sep = ": ")

sort(unique(pbmc$cell_type_legend))












tissue_list <- unique(pbmc$tissue)

tissue_df <- as.data.frame(tissue_list)
tissue_df$compartment <- compartment_name
names(tissue_df)[names(tissue_df) == "tissue_list"] <- "tissue"
tissue_df <- tissue_df[ , c("compartment", "tissue")]

tissue_df$tissue <- as.character(tissue_df$tissue)
tissue_df <- tissue_df[order(tissue_df$tissue), ]


write.table(tissue_df, file = paste0(sample_name, "_tissue_list.txt"), quote = FALSE, sep = "\t", row.names = FALSE)


# compare_cell_types <- as.data.frame(table(pbmc@meta.data[ , c("CellType", "cell_type_original")]))
# write.table(compare_cell_types, file = paste0(sample_name, "_compare_cell_types.txt"), quote = FALSE, sep = "\t", row.names = FALSE)







Idents(pbmc) <- "cell.type"



pdf(paste0(sample_name, "_UMAP_plot_cell-type_FA.pdf"), width = 9, height = 7, useDingbats = FALSE)
p <- DimPlot(pbmc, reduction = "umap", label = TRUE, group.by = "cell.type", pt.size = 0.5, label.size = 2.7) +
  theme(legend.text = element_text(size = 7)) +
  theme(plot.subtitle = element_text(hjust = 0.5)) +
  labs(subtitle = paste0("free annotation", " - ", paste0(tissue_name, " (", organism_name, ")"))) +
  labs(caption = dataset_name)
print(p)
dev.off()

pdf(paste0(sample_name, "_UMAP_plot_cell-type_CL.pdf"), width = 9, height = 7, useDingbats = FALSE)
p <- DimPlot(pbmc, reduction = "umap", label = TRUE, group.by = "cell_type", pt.size = 0.5, label.size = 2.7) +
  theme(legend.text = element_text(size = 7)) +
  theme(plot.subtitle = element_text(hjust = 0.5)) +
  labs(subtitle = paste0("Cell Ontology", " - ", paste0(tissue_name, " (", organism_name, ")"))) +
  labs(caption = dataset_name)
print(p)
dev.off()

pdf(paste0(sample_name, "_UMAP_plot_cell-type_FA_no-legend.pdf"), width = 9, height = 7, useDingbats = FALSE)
p <- DimPlot(pbmc, reduction = "umap", label = TRUE, group.by = "cell.type", pt.size = 1) +
  # ggtitle(sample_name) +
  NoLegend() +
  theme(plot.subtitle = element_text(hjust = 0.5)) +
  labs(subtitle = paste0("free annotation", " - ", paste0(tissue_name, " (", organism_name, ")"))) +
  labs(caption = dataset_name)
print(p)
dev.off()


pdf(paste0(sample_name, "_UMAP_plot_cell-type_CL_no-legend.pdf"), width = 9, height = 7, useDingbats = FALSE)
p <- DimPlot(pbmc, reduction = "umap", label = TRUE, group.by = "cell_type", pt.size = 1) +
  # ggtitle(sample_name) +
  NoLegend() +
  theme(plot.subtitle = element_text(hjust = 0.5)) +
  labs(subtitle = paste0("Cell Ontology", " - ", paste0(tissue_name, " (", organism_name, ")"))) +
  labs(caption = dataset_name)
print(p)
dev.off()




pdf(paste0(sample_name, "_UMAP_plot_cell-type_Kat.pdf"), width = 9, height = 7, useDingbats = FALSE)
p <- DimPlot(pbmc, reduction = "umap", label = TRUE, group.by = "Kat_cell_type", pt.size = 0.5, label.size = 2.7) +
  theme(legend.text = element_text(size = 7)) +
  theme(plot.subtitle = element_text(hjust = 0.5)) +
  labs(subtitle = paste0("Cell Type", " - ", paste0(tissue_name, " (", organism_name, ")"))) +
  labs(caption = dataset_name)
print(p)
dev.off()







pdf(paste0(sample_name, "_UMAP_plot_cell-type_number_fig.pdf"), width = 8, height = 6.25, useDingbats = FALSE)
p <- DimPlot(pbmc, reduction = "umap", label = TRUE, group.by = "cell_type_number", pt.size = 0.5, label.size = 10) +
  
  ggtitle(plot_title) +
  # theme(legend.text = element_text(size = 7)) +
  # theme(plot.subtitle = element_text(hjust = 0.5)) +
  # labs(subtitle = paste0("Cell Ontology", " - ", paste0(tissue_name, " (", organism_name, ")"))) +
  # labs(caption = dataset_name)
  theme(legend.position = "none",
        plot.title = element_text(size = 32, face = "plain"),
        axis.ticks.length = unit(0.1, "inch"),
        axis.title = element_text(size = 28, face = "plain"),
        axis.text = element_text(size = 28),
        
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









pdf(paste0(sample_name, "_UMAP_plot_tissue.pdf"), width = 9, height = 7, useDingbats = FALSE)
p <- DimPlot(pbmc, reduction = "umap", label = TRUE, group.by = "tissue", pt.size = 0.5, label.size = 2.7) +
  theme(legend.text = element_text(size = 7))
print(p)
dev.off()

pdf(paste0(sample_name, "_UMAP_plot_tissue_no-legend.pdf"), width = 9, height = 7, useDingbats = FALSE)
p <- DimPlot(pbmc, reduction = "umap", label = TRUE, group.by = "tissue", pt.size = 0.5, label.size = 2.7) +
                          NoLegend()
print(p)
dev.off()





current_compartment <- compartment_name


  
found_genes <- c()
found_genes <- gene_list_mgi$ensembl_gene_id[which(gene_list_mgi$ensembl_gene_id %in% pbmc@assays$RNA$data@Dimnames[[1]])]
  

found_genes_df <- gene_list_mgi[which(gene_list_mgi$ensembl_gene_id %in% found_genes), ]
write.table(found_genes_df, file = paste0(sample_name, "_found_genes.txt"), quote = FALSE, sep = "\t", row.names = FALSE)


missing_genes_df <- gene_list_mgi[which(!gene_list_mgi$ensembl_gene_id %in% found_genes), ]
write.table(missing_genes_df, file = paste0(sample_name, "_missing_genes.txt"), quote = FALSE, sep = "\t", row.names = FALSE)




found_genes_convert <- found_genes_df[ , c("gene", "ensembl_gene_id")]
row.names(found_genes_convert) <- found_genes_convert$ensembl_gene_id
found_genes_convert$ensembl_gene_id <- NULL
found_genes_convert <- as.vector(found_genes_convert)

found_genes_convert <- setNames(found_genes_df$gene, found_genes_df$ensembl_gene_id)


pdf(paste(sample_name, gene_type, "DotPlot_no-scale_x-axis-genes.pdf", sep = "_"), width = 20, height = 7)

p <- DotPlot(pbmc, features = found_genes, scale = FALSE) + RotatedAxis() +
  scale_x_discrete(labels = found_genes_convert)

print(p)
dev.off()


pdf(paste(sample_name, gene_type, "DotPlot_scaled_x-axis-genes.pdf", sep = "_"), width = 15, height = 7)

p <- DotPlot(pbmc, features = found_genes, scale = TRUE) + RotatedAxis() +
  ggtitle(sample_name) +
  scale_x_discrete(labels = found_genes_convert)

print(p)
dev.off()





found_genes_df_paper <- found_genes_df[found_genes_df$gene %in% plot_genes, ]

found_genes_paper <- found_genes_df_paper$ensembl_gene_id


found_genes_convert_paper <- found_genes_convert[found_genes_paper]




# dot plots for paper

# pdf(paste(sample_name, gene_type, "DotPlot_CL.pdf", sep = "_"), width = 9, height = 7)
# 
# p <- DotPlot(pbmc, features = found_genes_paper, group.by = "cell_type", scale = TRUE) + coord_flip() +
#   # ggtitle(paste0(current_compartment, " (", tissue_name, ")")) +
#   ggtitle(paste0(tissue_name, " (", organism_name, ")")) +
#       xlab("Gene") + 
#       ylab("Cell Type (Cell Ontology)") + 
#   labs(caption = dataset_name) +
#   scale_x_discrete(labels = found_genes_convert_paper) +
#     theme(axis.text.x = element_text(angle = 45, hjust = 1))
#   
# print(p)
# dev.off()



pdf(paste(sample_name, gene_type, "DotPlot_FA.pdf", sep = "_"), width = 9, height = 7)

# p <- DotPlot(pbmc, features = found_genes_paper, group.by = "cell_type", scale = TRUE) + coord_flip() +
  p <- DotPlot(pbmc, features = found_genes_paper, group.by = "cell.type", scale = TRUE) + coord_flip() +
  # ggtitle(paste0(current_compartment, " (", tissue_name, ")")) +
  ggtitle(paste0(tissue_name, " (", organism_name, ")")) +
      xlab("Gene") + 
      ylab("Cell Type (free annotation)") + 
    labs(caption = dataset_name) +
  scale_x_discrete(labels = found_genes_convert_paper) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
print(p)
dev.off()

dotplot_data <- p$data
write.csv(dotplot_data, file = paste(sample_name, gene_type, "DotPlot_FA_data.csv", sep = "_"), row.names = FALSE)



pdf(paste(sample_name, gene_type, "DotPlot_CL.pdf", sep = "_"), width = 9, height = 7)

# p <- DotPlot(pbmc, features = found_genes_paper, group.by = "cell_type", scale = TRUE) + coord_flip() +
p <- DotPlot(pbmc, features = found_genes_paper, group.by = "cell_type", scale = TRUE) + coord_flip() +
  # ggtitle(paste0(current_compartment, " (", tissue_name, ")")) +
  ggtitle(paste0(tissue_name, " (", organism_name, ")")) +
  xlab("Gene") + 
  ylab("Cell Type (Cell Ontology)") + 
  labs(caption = dataset_name) +
  scale_x_discrete(labels = found_genes_convert_paper) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

print(p)
dev.off()

dotplot_data <- p$data
write.csv(dotplot_data, file = paste(sample_name, gene_type, "DotPlot_CL_data.csv", sep = "_"), row.names = FALSE)





pdf(paste(sample_name, gene_type, "DotPlot_Kat-cell-type.pdf", sep = "_"), width = 9, height = 7)

# p <- DotPlot(pbmc, features = found_genes_paper, group.by = "cell_type", scale = TRUE) + coord_flip() +
p <- DotPlot(pbmc, features = found_genes_paper, group.by = "Kat_cell_type", scale = TRUE) + coord_flip() +
  # ggtitle(paste0(current_compartment, " (", tissue_name, ")")) +
  ggtitle(paste0(tissue_name, " (", organism_name, ")")) +
  xlab("Gene") + 
  ylab("Cell Type (Kat)") + 
  labs(caption = dataset_name) +
  scale_x_discrete(labels = found_genes_convert_paper) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

print(p)
dev.off()

dotplot_data <- p$data
write.csv(dotplot_data, file = paste(sample_name, gene_type, "DotPlot_Kat-cell-type_data.csv", sep = "_"), row.names = FALSE)







pdf(paste(sample_name, gene_type, "DotPlot_Kat-cell-type_fig.pdf", sep = "_"), width = 9, height = 7)

# p <- DotPlot(pbmc, features = found_genes_paper, group.by = "cell_type", scale = TRUE) + coord_flip() +
p <- DotPlot(pbmc, features = found_genes_paper, group.by = "Kat_cell_type", scale = TRUE) + coord_flip() +
  ggtitle(plot_title) +
  scale_x_discrete(labels = found_genes_convert_paper) +
  
    theme(axis.text.x = element_text(size = 22, angle = 45, hjust = 0.95, vjust = 1),
        axis.text.y = element_text(size = 24, face = "italic"),
        axis.title = element_blank(),     
        plot.title = element_text(size = 28, face = "plain", hjust = 0.5, vjust = 0),
      axis.ticks.length = unit(0.075, "inch"),
        # axis.title = element_text(size = 16, face = "plain"),
      legend.text = element_text(size = 21),
      legend.title = element_text(size = 21)

        
        
        )
        

print(p)
dev.off()

dotplot_data <- p$data
write.csv(dotplot_data, file = paste(sample_name, gene_type, "DotPlot_Kat-cell-type_data_fig.csv", sep = "_"), row.names = FALSE)


















cat("Making violin plots . . .", "\n\n")


summary_df <- data.frame()

tissue_cell_types_df <- data.frame()


for (current_tissue in tissue_list) {

  message(current_tissue)
  
  Idents(object = pbmc) <- "tissue"
  
  # current_cells <- WhichCells(object = pbmc, idents = current_tissue)

  # lung.temp_dataset <- subset(pbmc, cells = current_cells)
  pbmc.temp <- subset(pbmc, idents = current_tissue)
  
  # disease_list <- unique(lung.temp_dataset$disease2)
  
  # rm(lung.temp_dataset)
  # gc()
  
  
  current_tissue_cell_types_df <- as.data.frame(as.character(unique(pbmc.temp$cell.type)))
  
  names(current_tissue_cell_types_df)[1] <- "cell_type"
  
  current_tissue_cell_types_df$compartment <- compartment_name
  current_tissue_cell_types_df$tissue <- current_tissue
  
  current_tissue_cell_types_df <- current_tissue_cell_types_df[ c("compartment", "tissue", "cell_type")]
  
  current_tissue_cell_types_df <- current_tissue_cell_types_df[order(current_tissue_cell_types_df$cell_type), ]
  
  tissue_cell_types_df <- rbind(tissue_cell_types_df, current_tissue_cell_types_df)
  
  
  
  


  found_genes <- c()
  found_genes <- gene_list_mgi$ensembl_gene_id[which(gene_list_mgi$ensembl_gene_id %in% pbmc.temp@assays$RNA@data@Dimnames[[1]])]
  



  
  
  
  # annotation_type <- c("cell_type", "cell.type", "Kat_cell_type")
  annotation_type <- c("Kat_cell_type")
  
  
  for (current_annotation_type in annotation_type) {
    
    
    message(current_annotation_type)
    
    
    anno_short <- "error"
    if (current_annotation_type == "cell_type") {
      anno_short <- "CL"
    } else if (current_annotation_type == "cell.type") {
      anno_short <- "FA"
    } else if (current_annotation_type == "Kat_cell_type") {
      anno_short <- "Kat"
    }
    
    anno_long <- "error"
    if (current_annotation_type == "cell_type") {
      anno_long <- "Cell Ontology"
    } else if (current_annotation_type == "cell.type") {
      anno_long <- "free annotation"
    } else if (current_annotation_type == "Kat_cell_type") {
      anno_long <- "Kat"
    }
    
    
  
  
  
  
  
  
  
  
  
  
  
  
  
  # raster_pdf(paste0(gene_type, "_", cell_group_name, "_", current_tissue, "_cell-type_violin_plots.pdf"), res = 300)

    
    # for Fig 3
    pdf(paste0(sample_name, "_", gene_type, "_violin_plots_", anno_short, "_fig3.pdf"), useDingbats = FALSE, width = 6, height = 7)
    
    # for other figures
    # pdf(paste0(sample_name, "_", gene_type, "_violin_plots_", anno_short, "_fig.pdf"), useDingbats = FALSE, width = 7, height = 6)
  
    
  # pdf(paste0("NLR-genes_", current_dataset, "_", gsub("/", "_", current_cell_type), "_violin_plots.pdf"), useDingbats = FALSE)
  
  
  for (current_gene in found_genes) {
    
    
    
    Idents(object = pbmc.temp) <- "tissue"
    
    
    # matching_gene <- gene_list_df[gene_list_mgi$ensembl_gene_id == current_gene, "gene"]
    
    matching_gene <- gene_list_mgi[which(gene_list_mgi$ensembl_gene_id == current_gene), "gene"]
    
    
    message(paste0(matching_gene, " (", current_gene, ")"))
    
    
    vln_plot_title = paste0(matching_gene, " (", current_gene, ")\n", paste0(current_tissue, " (", organism_name, ")"))
    
    
    # message(current_gene)
    
    
    
    
    all_genes <- pbmc.temp@assays$RNA@data@Dimnames[[1]]
    
    gene_test <- all_genes[which(all_genes == current_gene)]
    
    if (length(gene_test) != 1) {
      
      stop("current gene length is not 1")
    }
    
    
    
    p <- VlnPlot(object = pbmc.temp, features = current_gene, group.by = current_annotation_type, pt.size = 1, raster = TRUE) +
      # p <- VlnPlot(object = pbmc.temp, features = current_gene, group.by = current_annotation_type, pt.size = 1, raster = TRUE, split.by = "seg_domain") +
      # p <- VlnPlot(object = pbmc.temp, features = current_gene, group.by = current_annotation_type, pt.size = 1, raster = TRUE, split.by = "seg_classic") +
        
        
      # p <- VlnPlot(object = pbmc.temp, features = current_gene, group.by = "cell_type", pt.size = 1, raster = TRUE) + 
      ggtitle(vln_plot_title) +
      xlab(paste0("Cell Type (", anno_long, ")")) +
      ylab("Expression Level") + 
      labs(caption = dataset_name) +
      theme(legend.position = "none") +
      theme(axis.text.x = element_text(size = rel(0.8)))
    
    p$layers[[2]]$aes_params$alpha <- 0.2
    
    print(p)
    
    
    
    
    
    
    
    
    
    
    # Violin Plot for Figure 3 - Nlrp6
    p <- VlnPlot(object = pbmc.temp, features = current_gene, group.by = current_annotation_type, pt.size = 1, raster = TRUE) +

      ggtitle(plot_title) +
      
      labs(subtitle = matching_gene) +

      # ylab(paste(plot_title, "Expression Level", sep = "\n")) + 
      ylab("Expression Level") + 

      theme(legend.position = "none", 
         axis.title.x = element_blank(),
      plot.title = element_text(size = 26, face = "plain", vjust = 0),
      axis.ticks.length = unit(0.1, "inch"),
        axis.title.y = element_text(size = 21.5, face = "plain"),
        axis.text.x = element_text(size = 20),
      axis.text.y = element_text(size = 21.5),
      plot.subtitle = element_text(size = 26, face = "italic", hjust = 0.5, vjust = 0)
      
      ) 
    
    p$layers[[2]]$aes_params$alpha <- 0.2
    
    print(p)
    
    
    
    
        p <- VlnPlot(object = pbmc.temp, features = current_gene, group.by = current_annotation_type, pt.size = 1, raster = TRUE) +

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
    
    
    
    
    
    
    
    

  
    # ***** need to change for cell type (CL or FA) used *****
    # cell_type_list <- unique(pbmc.temp@meta.data[ , current_annotation_type])
    cell_type_list <- unique(pbmc.temp@meta.data[ , current_annotation_type])
    
    
    # Idents(pbmc.temp) <- "cell.type"
    Idents(pbmc.temp) <- current_annotation_type
    
  
  for (current_cell_type in cell_type_list) {
  
    
    current_cells <- WhichCells(object = pbmc.temp, idents = current_cell_type)
  
  
            expression_values <- pbmc.temp@assays$RNA@data[current_gene, current_cells]

          if (!is.vector(expression_values)) {
               
               stop("expression values are not a vector")
          }
          

          
          current_summary <- data.frame(gene = matching_gene)
          
          current_summary$ensembl_id <- current_gene
          
          current_summary$tissue <- current_tissue
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
  
  

} # END current tissue


  write.csv(summary_df, file = paste(sample_name, gene_type, "tissue_summary.csv", sep = "_"), row.names = FALSE)


  tissue_cell_types_df <- tissue_cell_types_df[order(tissue_cell_types_df$tissue), ]
  
  # write.table(tissue_cell_types_df, file = paste(sample_name, gene_type, "tissue_cell-types_FA.txt", sep = "_"), quote = FALSE, sep = "\t", row.names = FALSE)
  
  
  
  
  
  
  

  
  
  
  
  
  
  
  



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



