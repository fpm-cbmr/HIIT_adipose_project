#### Load libraries ####
library('PhosR')
library('limma')
library('statmod')
library('dplyr')
library('tidyverse')
library('ggpubr')
library('Biobase')
library("AnnotationDbi")
library("org.Hs.eg.db")
library('ggnewscale')
library('ggridges')
library('enrichplot')
library('ggrepel')
library('dplyr')
library("ggfortify")
library("devtools")
library("ggbreak")
library("patchwork")
library("clusterProfiler")
library("GOSemSim")
library("rlang")
library("ggplot2")
library("gridExtra")
library("ggheatmap")
library("reshape2")
library("corrplot")
library("DEP")
library("ggnewscale")

#### Loading and filtering data ####
#loading clinical data and remove missing values
clincal_data <- read.delim("XXX") # Clinical data is not provided
clincal_data <- clincal_data[-c(which(clincal_data$New_ID == "")),]

#loading adipose proteome from MSfragger/DIA-NN
Adipose_proteome <- read.delim("diann-output.pg_matrix.tsv",
                               header = TRUE)
#Read grouping and sample ID
Name_adipose <- read.delim("File_names.txt", header = F)
Name_adipose <- Name_adipose[-c(23),] ## removing one outlier

#making a protein expression matrix and removing outlier (<1200 IDs)
Exprs_adipose <- Adipose_proteome[,grepl(".mzML",colnames(Adipose_proteome))]
Exprs_adipose <- Exprs_adipose[,-grep("23.mzML", colnames(Exprs_adipose))]

### trimming mzML column names and aligning clinical names with expression matrix
edit_names <- word(colnames(Exprs_adipose),10, sep="_")
edit_names <- as.integer(gsub("\\..*", "", edit_names))
edit_names_clin <- as.integer(word(clincal_data$New_ID,4, sep="_"))
Exprs_adipose <- Exprs_adipose[,order(edit_names)]
clinical_data2 <- clincal_data[order(edit_names_clin),]
colnames(Exprs_adipose) <- Name_adipose$V2
rownames(Exprs_adipose) <- Adipose_proteome$Protein.Ids

#Log2 transform expression data
Exprs_adipose = as.matrix(log2(Exprs_adipose))
rownames(Exprs_adipose) <- Adipose_proteome$Protein.Ids

# Defining variables
sample_name = strsplit(gsub("^_", "", colnames(Exprs_adipose)), "_")
df = S4Vectors::DataFrame(
  Group = sapply(sample_name, "[[", 1),
  Condition = sapply(sample_name, "[[", 2),
  replicate = sapply(sample_name, "[[", 4))
rownames(df) = colnames(Exprs_adipose)

#filtering data for 75% valid values
combined <- paste0(df$Group, df$Condition)
total <- rep("total",91)
Exprs_adipose_75filt <- selectGrps(Exprs_adipose,total, percent = 0.75)
complete_Exprs_adipose <- selectGrps(Exprs_adipose, combined, 1, n=6) # for PCA plot


#### PCA & batch effects ####
## PCA plot for 100% filtering of valid values
t_exprs <- as.data.frame(t(complete_Exprs_adipose))
t_exprs <- cbind(t_exprs, df$Group)
colnames(t_exprs)[1064] <- "Group"
pca_res <- prcomp(t_exprs[,1:1063], scale. = TRUE)

## visualising and defining three batch clusters
PCA_data_frame <- data.frame(pca_res$x)

ggplot(PCA_data_frame, aes(x=PC1,y=PC2)) + geom_point(aes(colour = df$Group))

# Defining the batch clusters
batch1 <- which(PCA_data_frame$PC2 < -5)
batch2 <- which(PCA_data_frame$PC1 < -15)
batch3 <- which(PCA_data_frame$PC2 > -5 & PCA_data_frame$PC1 > 0)

PCA_data_frame$batch <- 1
PCA_data_frame[batch1,]$batch <- 1
PCA_data_frame[batch2,]$batch <- 2
PCA_data_frame[batch3,]$batch <- 3
PCA_data_frame$batch <- factor(PCA_data_frame$batch)
PCA_data_frame <- cbind(PCA_data_frame,clinical_data2)

### PCA on batch-corrected data
corrected_df <- removeBatchEffect(complete_Exprs_adipose, batch = PCA_data_frame$batch)

t_exprs <- as.data.frame(t(corrected_df))
t_exprs <- cbind(t_exprs, df$Group, df$replicate)
colnames(t_exprs)[1064] <- "Group"
colnames(t_exprs)[1065] <- "Replicate"

pca_res <- prcomp(t_exprs[,1:1063], scale. = TRUE)

corrected_PCA_data <- data.frame(pca_res$x)
corrected_PCA_data$batch <- PCA_data_frame$batch
corrected_PCA_data$group <- df$Group
corrected_PCA_data$replicate <- df$replicate
corrected_PCA_data <- cbind(corrected_PCA_data, clinical_data2)
corrected_PCA_data$group <- factor(df$Group, levels = c('T2D', 'Obese', 'Lean'))

theme<-theme(axis.title = element_text(size = 20), axis.line = element_line(colour = "black")
             ,panel.background = element_blank(),
             panel.border=element_blank(), panel.grid.major = element_blank(),
             panel.grid.minor = element_blank(),strip.background=element_blank(),
             axis.text.x=element_text(colour="black", size = 18),
             axis.text.y=element_text(colour="black", size = 18),
             axis.ticks=element_line(colour="black"), legend.title = element_text(size=14),
             legend.text = element_text(size = 12))

percentage <- round(factoextra::get_eig(pca_res)$variance.percent, 2)
percentage <- paste(colnames(corrected_PCA_data), "(", paste( as.character(percentage), "%", ")", sep="") )

#Figure_S1C
ggplot(corrected_PCA_data, aes(x=PC1,y=PC2, colour = group)) + geom_point(aes(),size=3) +
  stat_ellipse(aes(colour = group)) + theme + xlab(percentage[1]) + ylab(percentage[2]) +
  labs(color='Group') + scale_color_manual(values=c("#CC6666", "#9999CC", "#616166"))

#Figure_1E
ggplot(corrected_PCA_data, aes(x=PC2,y=PC3, colour = group)) + geom_point(aes(),size=3) +
  stat_ellipse(aes(colour = group)) + theme + xlab(percentage[2]) + ylab(percentage[3]) +
  labs(color='Group') + scale_color_manual(values=c("#CC6666", "#6c6cba", "#616166"))


#Figure_S1D - PCA on Fat mass
ggplot(corrected_PCA_data, aes(x=PC2,y=PC3, color = FM1)) + geom_point(aes(), size = 3) +
  theme + xlab(percentage[2]) + ylab(percentage[3]) +
  scale_colour_gradient(
    low = "#d9dbfc",
    high = "#282c70",
    space = "Lab",
    na.value = "grey50",
    guide = "colourbar",
    aesthetics = "colour"
  ) + labs(color='Fat mass (kg)')


#### Saving and exporting proteome data ####

#save batch-corrected expression matrix for other analyses
Exprs_adipose_clean <- removeBatchEffect(Exprs_adipose_75filt, batch = PCA_data_frame$batch)
Exprs_adipose_whole_clean <- removeBatchEffect(Exprs_adipose, batch = PCA_data_frame$batch)

#Changing Protein accession ID's to Gene symbols:
rownames(Exprs_adipose_75filt) <- mapIds(org.Hs.eg.db, keys=rownames(Exprs_adipose_75filt), column="SYMBOL", keytype="ACCNUM", multiVals="first")
rownames(Exprs_adipose_clean) <- mapIds(org.Hs.eg.db, keys=rownames(Exprs_adipose_clean), column="SYMBOL", keytype="ACCNUM", multiVals="first")
rownames(Exprs_adipose) <- mapIds(org.Hs.eg.db, keys=rownames(Exprs_adipose), column="SYMBOL", keytype="ACCNUM", multiVals="first")
rownames(Exprs_adipose_whole_clean) <- mapIds(org.Hs.eg.db, keys=rownames(Exprs_adipose_whole_clean), column="SYMBOL", keytype="ACCNUM", multiVals="first")



### Export different processing files for down-stream analysis
Export_Exprs_adipose_75filt <- as.data.frame(Exprs_adipose_75filt)
Export_Exprs_adipose_75filt$ID <- rownames(Export_Exprs_adipose_75filt)
write.table(Export_Exprs_adipose_75filt,"Exprs_adipose_75filt.txt", append = FALSE, sep = "\t", dec = ".",
            row.names = F, col.names = T)

Export_Exprs_adipose_clean <- as.data.frame(Exprs_adipose_clean)
Export_Exprs_adipose_clean$ID <- rownames(Export_Exprs_adipose_clean)
write.table(Export_Exprs_adipose_clean,"Exprs_adipose_clean.txt", append = FALSE, sep = "\t", dec = ".",
            row.names = F, col.names = T)

Export_Exprs_adipose <- as.data.frame(Exprs_adipose)
Export_Exprs_adipose$ID <- rownames(Export_Exprs_adipose)
write.table(Export_Exprs_adipose,"Exprs_adipose.txt", append = FALSE, sep = "\t", dec = ".",
            row.names = F, col.names = T)

Export_Exprs_adipose_whole_clean <- as.data.frame(Exprs_adipose_whole_clean)
Export_Exprs_adipose_whole_clean$ID <- rownames(Export_Exprs_adipose_whole_clean)
write.table(Export_Exprs_adipose_whole_clean,"Exprs_adipose_whole_clean.txt", append = FALSE, sep = "\t", dec = ".",
            row.names = F, col.names = T)



#### QC data Figure_S1A-B ####
# Figure_S1A
Exprs_adipose_clean_QC <- Exprs_adipose_clean[,order(df$Group)]
colnames(Exprs_adipose_clean_QC) <- paste0(combined[order(df$Group)],"_", df$replicate[order(df$Group)])

corr_LeanPre <- cor(Exprs_adipose_clean_QC[,combined == "LeanPre"], method = "pearson", use = "complete.obs")
corr_LeanPre <- corr_LeanPre[upper.tri(corr_LeanPre)]
corr_LeanPost <- cor(Exprs_adipose_clean_QC[,combined == "LeanPost"], method = "pearson", use = "complete.obs")
corr_LeanPost <- corr_LeanPost[upper.tri(corr_LeanPost)]
corr_ObesePre <- cor(Exprs_adipose_clean_QC[,combined == "ObesePre"], method = "pearson", use = "complete.obs")
corr_ObesePre <- corr_ObesePre[upper.tri(corr_ObesePre)]
corr_ObesePost <- cor(Exprs_adipose_clean_QC[,combined == "ObesePost"], method = "pearson", use = "complete.obs")
corr_ObesePost <- corr_ObesePost[upper.tri(corr_ObesePost)]
corr_T2DPre <- cor(Exprs_adipose_clean_QC[,combined == "T2DPre"], method = "pearson", use = "complete.obs")
corr_T2DPre <- corr_T2DPre[upper.tri(corr_T2DPre)]
corr_T2DPost <- cor(Exprs_adipose_clean_QC[,combined == "T2DPost"], method = "pearson", use = "complete.obs")
corr_T2DPost <- corr_T2DPost[upper.tri(corr_T2DPost)]

corr_mat <- rbind(merge("LeanPre",corr_LeanPre), 
                  merge("LeanPost",corr_LeanPost),
                  merge("ObesePre",corr_ObesePre), 
                  merge("ObesePost",corr_ObesePost),
                  merge("T2DPre",corr_T2DPre), 
                  merge("T2DPost",corr_T2DPost))

ggplot(corr_mat,aes(factor(x, levels = unique(combined)),y)) +
  geom_boxplot() +
  theme_minimal() +
  theme(panel.grid = element_blank(),
        axis.line = element_line(), 
        axis.title = element_text(size = 30),
        axis.text = element_text(size = 18),
        axis.text.x = element_text(angle = 45, vjust = 0.5)) +
  labs(x = "Group", y = "Correlation Value (r)") 

median_vec <- c()
median_vec <- append(median_vec, median(corr_LeanPre,na.rm = T))
median_vec <- append(median_vec, median(corr_LeanPost,na.rm = T))
median_vec <- append(median_vec, median(corr_ObesePre,na.rm = T))
median_vec <- append(median_vec, median(corr_ObesePost,na.rm = T))
median_vec <- append(median_vec, median(corr_T2DPre,na.rm = T))
median_vec <- append(median_vec, median(corr_T2DPost,na.rm = T))

names(median_vec) <- unique(combined)
median_vec <- append(median_vec, median(c(corr_LeanPre,corr_LeanPost,
                                          corr_ObesePre,corr_ObesePost,
                                          corr_T2DPre,corr_T2DPost),na.rm = T))
names(median_vec)[7] <- "Total"
mean(median_vec) ## median Correlation 0.94%


# Figure_S1B
#2 random samples were used for plotting.
Sample_select = c(60,62)
Exprs_adipose_clean_SA_QC = Exprs_adipose_clean_QC[,Sample_select]
Exprs_adipose_clean_SA_QC <- na.omit(Exprs_adipose_clean_SA_QC)

#Histograms
barplot1 <- ggplot() + geom_histogram(aes(Exprs_adipose_clean_SA_QC[,1], alpha = 1000),position = "stack", color = "black") + 
  theme_minimal() + 
  theme(legend.position = "none",
        panel.grid = element_blank(),
        axis.line = element_line()) +
  theme(panel.grid = element_blank(),
        axis.line = element_line(), 
        axis.title = element_text(size = 30),
        axis.text = element_text(size = 22),
        axis.text.x = element_text(angle = 45, vjust = 0.5)) +
  labs(x = "Intensities", y = "Count")


barplot2 <- ggplot() + geom_histogram(aes(Exprs_adipose_clean_SA_QC[,2], alpha = 1000),position = "stack", color = "black") + 
  theme_minimal() + 
  theme(legend.position = "none",
        panel.grid = element_blank(),
        axis.line = element_line()) +
  theme(panel.grid = element_blank(),
        axis.line = element_line(), 
        axis.title = element_text(size = 30),
        axis.text = element_text(size = 22),
        axis.text.x = element_text(angle = 45, vjust = 0.5)) +
  labs(x = "Intensities", y = "Count")


#QQ-plot
qqplot1 <- ggplot(mapping = aes(sample = Exprs_adipose_clean_SA_QC[,1])) + 
  stat_qq_line(size = 1) +
  stat_qq(aes(alpha = 1000),color = "red") + 
  theme_minimal() + 
  theme(legend.position = "none",
        panel.grid = element_blank(),
        axis.line = element_line()) +
  theme(panel.grid = element_blank(),
        axis.line = element_line(), 
        axis.title = element_text(size = 30),
        axis.text = element_text(size = 22),
        axis.text.x = element_text(angle = 45, vjust = 0.5)) +
  labs(x = "Theoretical", y = "Sample")

qqplot2 <- ggplot(mapping = aes(sample = Exprs_adipose_clean_SA_QC[,2])) + 
  stat_qq_line(size = 1) +
  stat_qq(aes(alpha = 1000),color = "red") + 
  theme_minimal() + 
  theme(legend.position = "none",
        panel.grid = element_blank(),
        axis.line = element_line()) +
  theme(panel.grid = element_blank(),
        axis.line = element_line(), 
        axis.title = element_text(size = 30),
        axis.text = element_text(size = 22),
        axis.text.x = element_text(angle = 45, vjust = 0.5)) +
  labs(x = "Theoretical", y = "Sample")

#### Figure_1C Proteins_identified_plot ####
replicate <- word(Name_adipose$V2,2,4, sep="_")
design_group <- word(Name_adipose$V2,1, sep="_")
design_label <- Name_adipose$V2
data_unique <- make_unique(Adipose_proteome, "Genes", "Protein.Ids", delim = ";")
LFQ_columns <- grep("mzML", colnames(Adipose_proteome))
experimental_design <- data.frame(design_label,design_group,replicate)
colnames(experimental_design) <- c("label","condition","replicate")
colnames(data_unique)[6:96] <- experimental_design$label

data_se_parsed <- make_se_parse(data_unique, LFQ_columns)
data_se_parsed$condition <- design_group
data_se_parsed$replicate <- replicate

ALL_SA_data <- plot_numbers(data_se_parsed, plot = FALSE)
colnames(ALL_SA_data)[5] <- "Group"

ggplot(ALL_SA_data, aes(x=label, y=proteins_in_sample, fill = Group)) + 
  geom_col() + xlab("") +
  geom_hline(yintercept=mean(ALL_SA_data$proteins_in_sample), linetype="solid", 
             color = "black", size=1, alpha = 0.5) +
  ylab("Number of proteins") + theme +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank()) +
  scale_fill_manual(values=c("Dark grey", "#9999CC", "#CC6666"))



#### Figure_1D Ranked_abundance_plot ####
ranked_adipose_data <- apply(Adipose_proteome[,6:ncol(Adipose_proteome)], 1, median, na.rm=T)
names(ranked_adipose_data) <- Adipose_proteome$Genes
ranked_adipose_data <- sort(ranked_adipose_data, decreasing = T)
ranked_adipose_data <- log10(ranked_adipose_data)

ranked_adipose_data <- as.data.frame(cbind(ranked_adipose_data,names(ranked_adipose_data)))

ranked_adipose_data <- as.data.frame(ranked_adipose_data)
rownames(ranked_adipose_data) <- NULL
colnames(ranked_adipose_data) <- c("Log10_intensity","Protein_abundance_rank")

ranked_adipose_data$Log10_intensity <- as.numeric(ranked_adipose_data$Log10_intensity)
ranked_adipose_data$Protein_abundance_rank <- factor(ranked_adipose_data$Protein_abundance_rank)

ranked_adipose_data$rank <- 1:length(ranked_adipose_data$Log10_intensity)

options(ggrepel.max.overlaps = Inf)

ggplot(ranked_adipose_data, aes(y=Log10_intensity, x=sort(rank, decreasing = F))) + 
  geom_point(size = 2, alpha = 0.2, colour = "Dark grey") +
  theme(panel.background = element_blank(),
        axis.title.x = element_text(size=20, colour = "black", margin = margin(t=15,r=0,b=0,l=0)),
        axis.title.y = element_text(size=20, colour = "black", margin = margin(t=0,r=15,b=0,l=0)),
        axis.text = element_text(size = 16)) + scale_x_continuous(breaks = c(0,1000,2000,3000, 4000), limits = c(-50,4200), expand = c(0, 0.5)) +
  scale_y_continuous(expand = c(0, 0.5), limits = c(4,11)) +
  labs(x = "Abundance Rank", y = "Log10 Intensity") +
  theme(axis.line.x = element_line(color="black", size = 0.2),
        axis.line.y = element_line(color="black", size = 0.2)) +
  geom_text_repel(size = 6, aes(label=ifelse(Protein_abundance_rank == "ALB" | Protein_abundance_rank == "HBB"
                                             | Protein_abundance_rank == "HBA1"| Protein_abundance_rank == "FABP4"
                                             | Protein_abundance_rank == "PLIN1" | Protein_abundance_rank == "VCL"
                                             | Protein_abundance_rank == "RAB10" | Protein_abundance_rank == "SLC2A4"
                                             | Protein_abundance_rank == "LEP" | Protein_abundance_rank == "ADIPOQ", as.character(Protein_abundance_rank),''))) +
  geom_point(data = ranked_adipose_data %>% filter(Protein_abundance_rank == "ALB" & rank < 10), color = "red", size = 4, alpha = 1) +
  geom_point(data = ranked_adipose_data %>% filter(Protein_abundance_rank == "HBA1"), color = "red", size = 4, alpha = 1) +
  geom_point(data = ranked_adipose_data %>% filter(Protein_abundance_rank == "HBB"), color = "red", size = 4, alpha = 1) +
  geom_point(data = ranked_adipose_data %>% filter(Protein_abundance_rank == "FABP4"), color = "#7a7af0", size = 4, alpha = 1) +
  geom_point(data = ranked_adipose_data %>% filter(Protein_abundance_rank == "PLIN1"), color = "#7a7af0", size = 4, alpha = 1) +
  geom_point(data = ranked_adipose_data %>% filter(Protein_abundance_rank == "VCL"), color = "#7a7af0", size = 4, alpha = 1) +
  geom_point(data = ranked_adipose_data %>% filter(Protein_abundance_rank == "RAB10"), color = "#7a7af0", size = 4, alpha = 1) +
  geom_point(data = ranked_adipose_data %>% filter(Protein_abundance_rank == "SLC2A4"), color = "#7a7af0", size = 4, alpha = 1) +
  geom_point(data = ranked_adipose_data %>% filter(Protein_abundance_rank == "LEP"), color = "#7a7af0", size = 4, alpha = 1) +
  geom_point(data = ranked_adipose_data %>% filter(Protein_abundance_rank == "ADIPOQ"), color = "#7a7af0", size = 4, alpha = 1)