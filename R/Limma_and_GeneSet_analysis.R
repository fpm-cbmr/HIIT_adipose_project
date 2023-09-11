#### Limma setup ####
Exprs_adipose <- read.delim("Exprs_adipose.txt",header = TRUE)
rownames(Exprs_adipose) <- Exprs_adipose$ID
Exprs_adipose <- Exprs_adipose[,1:91]
Exprs_adipose_75filt <- read.delim("Exprs_adipose_75filt.txt",header = TRUE)
rownames(Exprs_adipose_75filt ) <- Exprs_adipose_75filt $ID
Exprs_adipose_75filt  <- Exprs_adipose_75filt [,1:91]

combined <- factor(combined, levels = c("LeanPre","LeanPost","ObesePre","ObesePost","T2DPre","T2DPost"))

# Adding cluster-batch effect as covariate
design <- model.matrix(~ 0 + combined + PCA_data_frame$batch)
colnames(design)[7:8] <- c("batch2","batch3")

corfit_75filt <- duplicateCorrelation(Exprs_adipose_75filt, design, block=df$replicate)
fit_75filt <- lmFit(Exprs_adipose_75filt,design,block=df$replicate,correlation=corfit_75filt$consensus)
corfit <- duplicateCorrelation(Exprs_adipose, design, block=df$replicate)
fit <- lmFit(Exprs_adipose,design,block=df$replicate,correlation=corfit$consensus)

cm <- makeContrasts(combinedLeanPost - combinedLeanPre,
                    combinedObesePost - combinedObesePre,
                    combinedT2DPost - combinedT2DPre,
                    (combinedLeanPost - combinedLeanPre + combinedObesePost - combinedObesePre + combinedT2DPost - combinedT2DPre)/6,
                    (combinedT2DPost - combinedT2DPre) - (combinedObesePost - combinedObesePre),
                    (combinedObesePost - combinedObesePre) - (combinedLeanPost - combinedLeanPre),
                    (combinedT2DPost - combinedT2DPre) - (combinedLeanPost - combinedLeanPre),
                    combinedObesePre - combinedLeanPre,
                    combinedT2DPre - combinedLeanPre,
                    combinedT2DPre - combinedObesePre, levels=design)

#### Limma for 75% proteome ####

fit2_75filt <- eBayes(contrasts.fit(fit_75filt, cm))

Lean_training_75filt <- topTable(fit2_75filt,
                                 coef="combinedLeanPost - combinedLeanPre",
                                 number = Inf, sort.by = "none")
Obese_training_75filt <- topTable(fit2_75filt,
                                  coef="combinedObesePost - combinedObesePre",
                                  number = Inf, sort.by = "none")
T2D_training_75filt <- topTable(fit2_75filt,
                                coef="combinedT2DPost - combinedT2DPre",
                                number = Inf, sort.by = "none")
MainEffect_training_75filt <- topTable(fit2_75filt,
                                       coef="(combinedLeanPost - combinedLeanPre + combinedObesePost - combinedObesePre + combinedT2DPost - combinedT2DPre)/6",
                                       number = Inf, sort.by = "none")
T2D_Obese_training_interaction_75filt <- topTable(fit2_75filt,
                                                  coef="(combinedT2DPost - combinedT2DPre) - (combinedObesePost - combinedObesePre)",
                                                  number = Inf, sort.by = "none")
Obese_Lean_training_interaction_75filt <- topTable(fit2_75filt,
                                                   coef="(combinedObesePost - combinedObesePre) - (combinedLeanPost - combinedLeanPre)",
                                                   number = Inf, sort.by = "none")
T2D_Lean_training_interaction_75filt <- topTable(fit2_75filt,
                                                 coef="(combinedT2DPost - combinedT2DPre) - (combinedLeanPost - combinedLeanPre)",
                                                 number = Inf, sort.by = "none")
ObesePre_LeanPre_75filt <- topTable(fit2_75filt,
                                    coef="combinedObesePre - combinedLeanPre",
                                    number = Inf, sort.by = "none")
T2DPre_LeanPre_75filt <- topTable(fit2_75filt,
                                  coef="combinedT2DPre - combinedLeanPre",
                                  number = Inf, sort.by = "none")
T2DPre_ObesePre_75filt <- topTable(fit2_75filt,
                                   coef="combinedT2DPre - combinedObesePre",
                                   number = Inf, sort.by = "none")
Main_effect_disease_75filt <- topTable(fit2_75filt,
                                       coef=8:10,
                                       number = Inf, sort.by = "none")

#### Limma for full proteome ####
fit2 <- eBayes(contrasts.fit(fit, cm))
Lean_training <- topTable(fit2,
                          coef="combinedLeanPost - combinedLeanPre",
                          number = Inf, sort.by = "none")
Obese_training <- topTable(fit2,
                           coef="combinedObesePost - combinedObesePre",
                           number = Inf, sort.by = "none")
T2D_training <- topTable(fit2,
                         coef="combinedT2DPost - combinedT2DPre",
                         number = Inf, sort.by = "none")
MainEffect_training <- topTable(fit2,
                                coef="(combinedLeanPost - combinedLeanPre + combinedObesePost - combinedObesePre + combinedT2DPost - combinedT2DPre)/6",
                                number = Inf, sort.by = "none")
T2D_Obese_training_interaction <- topTable(fit2,
                                           coef="(combinedT2DPost - combinedT2DPre) - (combinedObesePost - combinedObesePre)",
                                           number = Inf, sort.by = "none")
Obese_Lean_training_interaction <- topTable(fit2,
                                            coef="(combinedObesePost - combinedObesePre) - (combinedLeanPost - combinedLeanPre)",
                                            number = Inf, sort.by = "none")
T2D_Lean_training_interaction <- topTable(fit2,
                                          coef="(combinedT2DPost - combinedT2DPre) - (combinedLeanPost - combinedLeanPre)",
                                          number = Inf, sort.by = "none")
ObesePre_LeanPre <- topTable(fit2,
                             coef="combinedObesePre - combinedLeanPre",
                             number = Inf, sort.by = "none")
T2DPre_LeanPre <- topTable(fit2,
                           coef="combinedT2DPre - combinedLeanPre",
                           number = Inf, sort.by = "none")
T2DPre_ObesePre <- topTable(fit2,
                            coef="combinedT2DPre - combinedObesePre",
                            number = Inf, sort.by = "none")
Main_effect_disease <- topTable(fit2,
                                coef=8:10,
                                number = Inf, sort.by = "none")



### export
empt_matrix <- matrix(ncol=12,nrow=nrow(Exprs_adipose))
for(i in 1:nrow(Exprs_adipose)){
  empt_matrix[i,1] <- mean(na.omit(Exprs_adipose[i,which(df$Group == "Lean" & df$Condition == "Pre")]))
  empt_matrix[i,2] <- sd(na.omit(Exprs_adipose[i,which(df$Group == "Lean" & df$Condition == "Pre")]))/sqrt(length(na.omit((Exprs_adipose[i,which(df$Group == "Lean" & df$Condition == "Pre")]))))
  empt_matrix[i,3] <- mean(na.omit(Exprs_adipose[i,which(df$Group == "Lean" & df$Condition == "Post")]))
  empt_matrix[i,4] <- sd(na.omit(Exprs_adipose[i,which(df$Group == "Lean" & df$Condition == "Post")]))/sqrt(length(na.omit((Exprs_adipose[i,which(df$Group == "Lean" & df$Condition == "Post")]))))
  empt_matrix[i,5] <- mean(na.omit(Exprs_adipose[i,which(df$Group == "Obese" & df$Condition == "Pre")]))
  empt_matrix[i,6] <- sd(na.omit(Exprs_adipose[i,which(df$Group == "Obese" & df$Condition == "Pre")]))/sqrt(length(na.omit((Exprs_adipose[i,which(df$Group == "Obese" & df$Condition == "Pre")]))))
  empt_matrix[i,7] <- mean(na.omit(Exprs_adipose[i,which(df$Group == "Obese" & df$Condition == "Post")]))
  empt_matrix[i,8] <- sd(na.omit(Exprs_adipose[i,which(df$Group == "Obese" & df$Condition == "Post")]))/sqrt(length(na.omit((Exprs_adipose[i,which(df$Group == "Obese" & df$Condition == "Post")]))))
  empt_matrix[i,9] <- mean(na.omit(Exprs_adipose[i,which(df$Group == "T2D" & df$Condition == "Pre")]))
  empt_matrix[i,10] <- sd(na.omit(Exprs_adipose[i,which(df$Group == "T2D" & df$Condition == "Pre")]))/sqrt(length(na.omit((Exprs_adipose[i,which(df$Group == "T2D" & df$Condition == "Pre")]))))
  empt_matrix[i,11] <- mean(na.omit(Exprs_adipose[i,which(df$Group == "T2D" & df$Condition == "Post")]))
  empt_matrix[i,12] <- sd(na.omit(Exprs_adipose[i,which(df$Group == "T2D" & df$Condition == "Post")]))/sqrt(length(na.omit((Exprs_adipose[i,which(df$Group == "T2D" & df$Condition == "Post")]))))
}
colnames(empt_matrix) <- c("Lean_Pre_mean","Lean_Pre_SEM","Lean_Post_mean","Lean_Post_SEM",
                           "Obese_Pre_mean","Obese_Pre_SEM","Obese_Post_mean","Obese_Post_SEM",
                           "T2D_Pre_mean","T2D_Pre_SEM","T2D_Post_mean","T2D_Post_SEM")

Export_matrix_full_prot <- cbind(T2DPre_LeanPre[,c(1,2,5)],T2DPre_ObesePre[,c(2,5)],ObesePre_LeanPre[,c(2,5)],Lean_training[,c(2,5)],Obese_training[,c(2,5)],T2D_training[,c(2,5)],empt_matrix) ## full 3773 protein matrix

empt_matrix <- matrix(ncol=12,nrow=nrow(Exprs_adipose_75filt))
for(i in 1:nrow(Exprs_adipose_75filt)){
  empt_matrix[i,1] <- mean(na.omit(Exprs_adipose_75filt[i,which(df$Group == "Lean" & df$Condition == "Pre")]))
  empt_matrix[i,2] <- sd(na.omit(Exprs_adipose_75filt[i,which(df$Group == "Lean" & df$Condition == "Pre")]))/sqrt(length(na.omit((Exprs_adipose_75filt[i,which(df$Group == "Lean" & df$Condition == "Pre")]))))
  empt_matrix[i,3] <- mean(na.omit(Exprs_adipose_75filt[i,which(df$Group == "Lean" & df$Condition == "Post")]))
  empt_matrix[i,4] <- sd(na.omit(Exprs_adipose_75filt[i,which(df$Group == "Lean" & df$Condition == "Post")]))/sqrt(length(na.omit((Exprs_adipose_75filt[i,which(df$Group == "Lean" & df$Condition == "Post")]))))
  empt_matrix[i,5] <- mean(na.omit(Exprs_adipose_75filt[i,which(df$Group == "Obese" & df$Condition == "Pre")]))
  empt_matrix[i,6] <- sd(na.omit(Exprs_adipose_75filt[i,which(df$Group == "Obese" & df$Condition == "Pre")]))/sqrt(length(na.omit((Exprs_adipose_75filt[i,which(df$Group == "Obese" & df$Condition == "Pre")]))))
  empt_matrix[i,7] <- mean(na.omit(Exprs_adipose_75filt[i,which(df$Group == "Obese" & df$Condition == "Post")]))
  empt_matrix[i,8] <- sd(na.omit(Exprs_adipose_75filt[i,which(df$Group == "Obese" & df$Condition == "Post")]))/sqrt(length(na.omit((Exprs_adipose_75filt[i,which(df$Group == "Obese" & df$Condition == "Post")]))))
  empt_matrix[i,9] <- mean(na.omit(Exprs_adipose_75filt[i,which(df$Group == "T2D" & df$Condition == "Pre")]))
  empt_matrix[i,10] <- sd(na.omit(Exprs_adipose_75filt[i,which(df$Group == "T2D" & df$Condition == "Pre")]))/sqrt(length(na.omit((Exprs_adipose_75filt[i,which(df$Group == "T2D" & df$Condition == "Pre")]))))
  empt_matrix[i,11] <- mean(na.omit(Exprs_adipose_75filt[i,which(df$Group == "T2D" & df$Condition == "Post")]))
  empt_matrix[i,12] <- sd(na.omit(Exprs_adipose_75filt[i,which(df$Group == "T2D" & df$Condition == "Post")]))/sqrt(length(na.omit((Exprs_adipose_75filt[i,which(df$Group == "T2D" & df$Condition == "Post")]))))
}
colnames(empt_matrix) <- c("Lean_Pre_mean","Lean_Pre_SEM","Lean_Post_mean","Lean_Post_SEM",
                           "Obese_Pre_mean","Obese_Pre_SEM","Obese_Post_mean","Obese_Post_SEM",
                           "T2D_Pre_mean","T2D_Pre_SEM","T2D_Post_mean","T2D_Post_SEM")

Export_matrix_75filt_prot <- cbind(T2DPre_LeanPre_75filt[,c(1,2,5,6)],T2DPre_ObesePre_75filt[,c(2,5,6)],ObesePre_LeanPre_75filt[,c(2,5,6)],Lean_training_75filt[,c(2,5,6)],Obese_training_75filt[,c(2,5,6)],T2D_training_75filt[,c(2,5,6)],empt_matrix) ### 2016 protein matrix



#### volcano plots ####
theme_vol <-theme(axis.title = element_text(size = 16), axis.line = element_line(colour = "black", size = 2)
                  ,panel.background = element_blank(),
                  panel.border= element_blank(), panel.grid.major = element_blank(),
                  panel.grid.minor = element_blank(),strip.background=element_blank(),
                  axis.text.x=element_text(colour="black", size = 16),
                  axis.text.y=element_text(colour="black", size = 16),
                  axis.ticks=element_line(colour="black"), legend.position="none")


#Obese vs Lean Figure_S2B
ObesePre_LeanPre_75filt$diffexpressed <- "NO"
ObesePre_LeanPre_75filt$diffexpressed[ObesePre_LeanPre_75filt$logFC > 0 & ObesePre_LeanPre_75filt$P.Val < 0.01] <- "UP"
ObesePre_LeanPre_75filt$diffexpressed[ObesePre_LeanPre_75filt$logFC < 0 & ObesePre_LeanPre_75filt$P.Val < 0.01] <- "DOWN"
ObesePre_LeanPre_75filt$delabel <- NA
ObesePre_LeanPre_75filt$delabel[ObesePre_LeanPre_75filt$diffexpressed != "NO"] <- ObesePre_LeanPre_75filt$ID[ObesePre_LeanPre_75filt$diffexpressed != "NO"]
plot_ObesePre_LeanPre_75filt <- ggplot(data=ObesePre_LeanPre_75filt, aes(x=logFC, y=-log10(P.Value), col=diffexpressed, label = delabel)) + geom_point(alpha = 0.4, size = 3)+ xlim(-1, 2) + scale_color_manual(values=c("black", "lightgrey", "#7a7af0")) + ggtitle("Obese vs Lean") + theme_vol + geom_text_repel()
plot_ObesePre_LeanPre_75filt +  scale_size_manual(values = c("DOWN" = 4, "NO" = 2, "UP" = 4)) + xlab("LogFC (Obese - Lean)") + ylab("P.value (-Log10)") +
  geom_hline(yintercept=2, linetype="dashed",
             color = "black", size=1, alpha = 0.3)

#T2D vs Lean Figure_2A
T2DPre_LeanPre_75filt$diffexpressed <- "NO"
T2DPre_LeanPre_75filt$diffexpressed[T2DPre_LeanPre_75filt$logFC > 0 & T2DPre_LeanPre_75filt$adj.P.Val < 0.1] <- "UP"
T2DPre_LeanPre_75filt$diffexpressed[T2DPre_LeanPre_75filt$logFC < 0 & T2DPre_LeanPre_75filt$adj.P.Val < 0.1] <- "DOWN"
T2DPre_LeanPre_75filt$delabel <- NA
T2DPre_LeanPre_75filt$delabel[T2DPre_LeanPre_75filt$diffexpressed != "NO"] <- T2DPre_LeanPre_75filt$ID[T2DPre_LeanPre_75filt$diffexpressed != "NO"]
plot_T2DPre_LeanPre_75filt <- ggplot(data=T2DPre_LeanPre_75filt, aes(x=logFC, y=-log10(P.Value), col=diffexpressed, label = delabel)) + geom_point(alpha = 0.4, size = 3)+ xlim(-1, 2) + scale_color_manual(values=c("black", "light grey", "darkred")) + ggtitle("T2D vs Lean") + theme_vol + geom_text_repel(max.overlaps = Inf)
plot_T2DPre_LeanPre_75filt + scale_size_manual(values = c("DOWN" = 4, "NO" = 2, "UP" = 4)) + xlab("LogFC (T2D - Lean)") + ylab("P.value (-Log10)")


#T2D vs OBESE Figure_S2A
T2DPre_ObesePre_75filt$diffexpressed <- "NO"
T2DPre_ObesePre_75filt$diffexpressed[T2DPre_ObesePre_75filt$logFC > 0 & T2DPre_ObesePre_75filt$P.Value < 0.01] <- "UP"
T2DPre_ObesePre_75filt$diffexpressed[T2DPre_ObesePre_75filt$logFC < 0 & T2DPre_ObesePre_75filt$P.Value < 0.01] <- "DOWN"
T2DPre_ObesePre_75filt$delabel <- NA
T2DPre_ObesePre_75filt$delabel[T2DPre_ObesePre_75filt$diffexpressed != "NO"] <- T2DPre_ObesePre_75filt$ID[T2DPre_ObesePre_75filt$diffexpressed != "NO"]
plot_T2DPre_ObesePre_75filt <- ggplot(data=T2DPre_ObesePre_75filt, aes(x=logFC, y=-log10(P.Value), col=diffexpressed, label = delabel)) + geom_point(alpha = 0.4, size = 3)+ xlim(-1, 2) + scale_color_manual(values=c("#7a7af0", "light grey", "darkred")) + ggtitle("T2D vs Obese") + theme_vol + geom_text_repel()
plot_T2DPre_ObesePre_75filt +  scale_size_manual(values = c("DOWN" = 4, "NO" = 2, "UP" = 4)) + xlab("LogFC (T2D - Obese)") + ylab("P.value (-Log10)") +
  geom_hline(yintercept=2, linetype="dashed",
             color = "black", size=1, alpha = 0.3)

#### Training volcano effects
#Figure_4A_Lean
Lean_training_75filt$diffexpressed <- "NO"
Lean_training_75filt$diffexpressed[Lean_training_75filt$logFC > 0 & Lean_training_75filt$adj.P.Val < 0.1] <- "UP"
Lean_training_75filt$diffexpressed[Lean_training_75filt$logFC < 0 & Lean_training_75filt$adj.P.Val < 0.1] <- "DOWN"
Lean_training_75filt$delabel <- NA
Lean_training_75filt$delabel[Lean_training_75filt$diffexpressed != "NO"] <- Lean_training_75filt$ID[Lean_training_75filt$diffexpressed != "NO"]
plot_Lean_training_75filt <- ggplot(data=Lean_training_75filt, aes(x=logFC, y=-log10(P.Value), col=diffexpressed, label = delabel)) + geom_point(alpha = 0.4, size = 3)+ xlim(-1.5, 1.5) + scale_color_manual(values=c("#2e136e","lightgrey", "#2e136e")) + ggtitle("Lean") + theme_vol + geom_text_repel()
plot_Lean_training_75filt + xlab("LogFC (Post - Pre)") + ylab("P.value (-Log10)")

#Figure_4A_Obese
Obese_training_75filt$diffexpressed <- "NO"
Obese_training_75filt$diffexpressed[Obese_training_75filt$logFC > 0 & Obese_training_75filt$adj.P.Val < 0.05] <- "UP"
Obese_training_75filt$diffexpressed[Obese_training_75filt$logFC < 0 & Obese_training_75filt$adj.P.Val < 0.05] <- "DOWN"
Obese_training_75filt$delabel <- NA
Obese_training_75filt$delabel[Obese_training_75filt$diffexpressed != "NO"] <- Obese_training_75filt$ID[Obese_training_75filt$diffexpressed != "NO"]
plot_Obese_training_75filt <- ggplot(data=Obese_training_75filt, aes(x=logFC, y=-log10(P.Value), col=diffexpressed, label = delabel)) + geom_point(alpha = 0.4, size = 3)+ xlim(-1.5, 1.5) + scale_color_manual(values=c("lightgrey", "#2e136e")) + ggtitle("Obese") + theme_vol + geom_text_repel()
plot_Obese_training_75filt + xlab("LogFC (Post - Pre)") + ylab("P.value (-Log10)")

#Figure_4A_T2D
T2D_training_75filt$diffexpressed <- "NO"
T2D_training_75filt$diffexpressed[T2D_training_75filt$logFC > 0 & T2D_training_75filt$adj.P.Val < 0.05] <- "UP"
T2D_training_75filt$diffexpressed[T2D_training_75filt$logFC < 0 & T2D_training_75filt$adj.P.Val < 0.05] <- "DOWN"
T2D_training_75filt$delabel <- NA
T2D_training_75filt$delabel[T2D_training_75filt$diffexpressed != "NO"] <- T2D_training_75filt$ID[T2D_training_75filt$diffexpressed != "NO"]
plot_T2D_training_75filt <- ggplot(data=T2D_training_75filt, aes(x=logFC, y=-log10(P.Value), col=diffexpressed, label = delabel)) + geom_point(alpha = 0.4, size = 3)+ xlim(-1.5, 1.5) + scale_color_manual(values=c("lightgrey", "#2e136e")) + ggtitle("T2D") + theme_vol + geom_text_repel()
plot_T2D_training_75filt + xlab("LogFC (Post - Pre)") + ylab("P.value (-Log10)")




#### Profile plot Figure_2B #####
signif_proteins <- T2DPre_LeanPre_75filt$ID[which(T2DPre_LeanPre_75filt$adj.P.Val < 0.1)]
signif_proteins <- rownames(T2DPre_LeanPre_75filt)[which(T2DPre_LeanPre_75filt$adj.P.Val < 0.1)]
Signi_exprs_list <- Exprs_adipose_75filt[which(rownames(Exprs_adipose_75filt) %in% signif_proteins),]
#extracting only pre values
Signi_exprs_list_pre <- Signi_exprs_list[,-c(which(grepl("Post",colnames(Signi_exprs_list))))]

#calculating average protein abundance per group
Lean_average <- apply(Signi_exprs_list_pre[,which(grepl("Lean",colnames(Signi_exprs_list_pre)))],1,mean, na.rm=TRUE)
Obese_average <- apply(Signi_exprs_list_pre[,which(grepl("Obese",colnames(Signi_exprs_list_pre)))],1,mean, na.rm=TRUE)
T2D_average <- apply(Signi_exprs_list_pre[,which(grepl("T2D",colnames(Signi_exprs_list_pre)))],1,mean, na.rm=TRUE)


log_fc_data <- data.frame(cbind(0,Obese_average-Lean_average,T2D_average-Lean_average, abs(T2D_average-Obese_average)))
colnames(log_fc_data) <- c("Lean","Obese","T2D","T2D_Obese")
log_fc_data$Gene <- rownames(log_fc_data)
log_fc_data$Obese_Logfc <- abs(log_fc_data$Obese) 
log_fc_data$T2D_Logfc <- abs(log_fc_data$T2D) 
long_logfc_data <- gather(log_fc_data,ID,Intensity,-Gene, -Obese_Logfc, -T2D_Logfc, -T2D_Obese)

#Color to highlight difference in profile
long_logfc_data$label <- 2
long_logfc_data$label[which(grepl("ANXA3",long_logfc_data$Gene))] <- 1
long_logfc_data$label[which(grepl("ITGAV",long_logfc_data$Gene))] <- 3
long_logfc_data$label[which(grepl("ALDH7A1",long_logfc_data$Gene))] <- 3
long_logfc_data$label[which(grepl("HTRA1",long_logfc_data$Gene))] <- 1

ggplot(long_logfc_data, aes(x=factor(ID),y=Intensity, colour=factor(label))) + geom_point(aes(group=factor(Gene)),size=0.7) + 
  geom_line(aes(group=factor(Gene)),alpha=0.5) + 
  theme_minimal() + theme(axis.line = element_line(linewidth = 0.1, colour = "black"),
                          plot.title = element_text(hjust = 0.5, size = 18), axis.text = element_text(size=16, hjust = 1),
                          axis.title = element_text(size = 18), legend.text = element_text(size=14),
                          legend.title = element_text(size=16)) + xlab("") + ylab("Log2 logFC") +
  ylim(-0.6,1.8) +
  scale_y_break(c(0.8, 1.2), scale=0.1, space = 0.05) + theme(legend.position = "none",panel.spacing = unit(0.35, "cm")) +
  scale_y_continuous(breaks = c(-0.4,0,0.4,0.8,1.4,1.6)) +
  scale_color_manual(values = c("darkred", "grey","blue"))







#### Boxplot pre-comparison Figure_2C & S2C #####
Baseline_samples <- Exprs_adipose_75filt[,grepl("Pre", colnames(Exprs_adipose_75filt))]

HTRA1_barplots <- as.data.frame(Baseline_samples[which(rownames(Baseline_samples) == "HTRA1"),])
HTRA1_barplots$Sample_ID <- rownames(HTRA1_barplots)
colnames(HTRA1_barplots)[1] <- "Intensity"
HTRA1_barplots$Group <- factor(word(HTRA1_barplots$Sample_ID,1, sep="_"), levels = c("Lean","Obese","T2D"))

HTRA1_plot <- ggplot(HTRA1_barplots, aes(y=Intensity, x=Group, fill=Group)) + geom_boxplot(alpha = .2) +
  scale_fill_manual(values=c("darkgrey", "#9999CC", "#CC6666")) +
  geom_point(position = position_jitterdodge(),size=2) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5, size=26),axis.text.x = element_text(angle = 45, size = 28, hjust = 1),
        axis.text.y = element_text(size = 24, hjust = 1),legend.position = "none",
        axis.title.y = element_text(size = 24)) +
  xlab("") + ylab("Intensity (Log2)") + ggtitle("HTRA1")
HTRA1_plot

ITGAV_barplots <- as.data.frame(Baseline_samples[which(rownames(Baseline_samples) == "ITGAV"),])
ITGAV_barplots$Sample_ID <- rownames(ITGAV_barplots)
colnames(ITGAV_barplots)[1] <- "Intensity"
ITGAV_barplots$Group <- factor(word(ITGAV_barplots$Sample_ID,1, sep="_"), levels = c("Lean","Obese","T2D"))

ITGAV_plot <- ggplot(ITGAV_barplots, aes(y=Intensity, x=Group, fill=Group)) + geom_boxplot(alpha = .2) +
  scale_fill_manual(values=c("darkgrey", "#9999CC", "#CC6666")) +
  geom_point(position = position_jitterdodge(),size=2) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5, size=26),axis.text.x = element_text(angle = 45, size = 28, hjust = 1),
        axis.text.y = element_text(size = 24, hjust = 1),legend.position = "none",
        axis.title.y = element_text(size = 24)) + scale_y_continuous(breaks = c(19,19.5,20,20.5,21,21.5,22)) +
  xlab("") + ylab("Intensity (Log2)") + ggtitle("ITGAV")
ITGAV_plot

ALDH7A1_barplots <- as.data.frame(Baseline_samples[which(rownames(Baseline_samples) == "ALDH7A1"),])
ALDH7A1_barplots$Sample_ID <- rownames(ALDH7A1_barplots)
colnames(ALDH7A1_barplots)[1] <- "Intensity"
ALDH7A1_barplots$Group <- factor(word(ALDH7A1_barplots$Sample_ID,1, sep="_"), levels = c("Lean","Obese","T2D"))

ALDH7A1_plot <- ggplot(ALDH7A1_barplots, aes(y=Intensity, x=Group, fill=Group)) + geom_boxplot(alpha = .2) +
  scale_fill_manual(values=c("darkgrey", "#9999CC", "#CC6666")) +
  geom_point(position = position_jitterdodge(),size=2) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5, size=26),axis.text.x = element_text(angle = 45, size = 28, hjust = 1),
        axis.text.y = element_text(size = 24, hjust = 1),legend.position = "none",
        axis.title.y = element_text(size = 24)) +
  xlab("") + ylab("Intensity (Log2)") + ggtitle("ALDH7A1")
ALDH7A1_plot

ANXA3_barplots <- as.data.frame(Baseline_samples[which(rownames(Baseline_samples) == "ANXA3"),])
ANXA3_barplots$Sample_ID <- rownames(ANXA3_barplots)
colnames(ANXA3_barplots)[1] <- "Intensity"
ANXA3_barplots$Group <- factor(word(ANXA3_barplots$Sample_ID,1, sep="_"), levels = c("Lean","Obese","T2D"))

ANXA3_plot <- ggplot(ANXA3_barplots, aes(y=Intensity, x=Group, fill=Group)) + geom_boxplot(alpha = .2) +
  scale_fill_manual(values=c("darkgrey", "#9999CC", "#CC6666")) +
  geom_point(position = position_jitterdodge(),size=2) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5, size=26),axis.text.x = element_text(angle = 45, size = 28, hjust = 1),
        axis.text.y = element_text(size = 24, hjust = 1),legend.position = "none",
        axis.title.y = element_text(size = 24)) +
  xlab("") + ylab("Intensity (Log2)") + ggtitle("ANXA3")
ANXA3_plot


#Figure_S2C
LEP_barplots <- as.data.frame(Baseline_samples[which(rownames(Baseline_samples) == "LEP"),])
LEP_barplots$Sample_ID <- rownames(LEP_barplots)
colnames(LEP_barplots)[1] <- "Intensity"
LEP_barplots$Group <- factor(word(LEP_barplots$Sample_ID,1, sep="_"), levels = c("Lean","Obese","T2D"))

LEP_plot <- ggplot(LEP_barplots, aes(y=Intensity, x=Group, fill=Group)) + geom_boxplot(alpha = .2) +
  scale_fill_manual(values=c("darkgrey", "#9999CC", "#CC6666")) +
  geom_point(position = position_jitterdodge(),size=2) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5, size=26),axis.text.x = element_text(angle = 45, size = 28, hjust = 1),
        axis.text.y = element_text(size = 24, hjust = 1),legend.position = "none",
        axis.title.y = element_text(size = 24)) +
  xlab("") + ylab("Intensity (Log2)") + ggtitle("Leptin")
LEP_plot

# Figure_S2E,F,H,I
WB_HTRA1_data <- read.delim("WB_HTRA1_data.txt")
WB_ITGAV_data <- read.delim("WB_ITGAV_data.txt")

pre_WB_HTRA1 <- WB_HTRA1_data[WB_HTRA1_data$Training == "Pre",]
pre_WB_ITGAV <- WB_ITGAV_data[WB_ITGAV_data$Training == "Pre",]

ggplot(pre_WB_HTRA1, aes(y=Intensity, x=Group, fill=Group)) + geom_boxplot(alpha = .2) +
  scale_fill_manual(values=c("darkgrey", "#9999CC", "#CC6666")) +
  geom_point(position = position_jitterdodge(),size=1) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5, size=26),axis.text.x = element_text(angle = 45, size = 28, hjust = 1),
        axis.text.y = element_text(size = 24, hjust = 1),legend.position = "none",
        axis.title.y = element_text(size = 24)) +
  xlab("") + ylab("Intensity (Log2)") + ggtitle("HTRA1") + ylim(0,8.5)

ggplot(pre_WB_ITGAV, aes(y=Intensity, x=Group, fill=Group)) + geom_boxplot(alpha = .2) +
  scale_fill_manual(values=c("darkgrey", "#9999CC", "#CC6666")) +
  geom_point(position = position_jitterdodge(),size=1) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5, size=26),axis.text.x = element_text(angle = 45, size = 28, hjust = 1),
        axis.text.y = element_text(size = 24, hjust = 1),legend.position = "none",
        axis.title.y = element_text(size = 24)) +
  xlab("") + ylab("Intensity (Log2)") + ggtitle("ITGAV") + ylim(0,4)


cor.test(WB_HTRA1_data$Intensity,WB_HTRA1_data$Prot_intensity, method="kendall")
cor.test(WB_ITGAV_data$Intensity,WB_ITGAV_data$Prot_intensity, method="kendall")

ggplot(WB_HTRA1_data, aes(x=log2(Intensity), y=Prot_intensity)) + geom_point(aes(color = Group)) +
  geom_smooth(method = "lm", se = TRUE, colour = "darkblue", size = 0.75, alpha=0.15, fill="#808096") +
  theme_minimal() + theme(axis.line = element_line(linewidth = 0.1, colour = "black"),
                          plot.title = element_text(hjust = 0.5, size = 18), axis.text = element_text(size=16, hjust = 1),
                          axis.title = element_text(size = 18), legend.text = element_text(size=14),
                          legend.title = element_text(size=16)
  ) + scale_color_manual(values=c("darkgrey", "#7a7af0", "#CC6666")) +
  ggtitle("HTRA1") +
  xlab("Western blot (Log2[AU])") + ylab("Proteomics (Log2)") +
  annotate(geom="text",x=0,y=23, label="r = 0.48", color = "black",size=6) +
  annotate(geom = "text",x=0,y=22.7,label="p = 2.6e-11",size=6)


ggplot(WB_ITGAV_data, aes(x=log2(Intensity), y=Prot_intensity)) + geom_point(aes(color = Group)) +
  geom_smooth(method = "lm", se = TRUE, colour = "darkblue", size = 0.75, alpha=0.15, fill="#808096") +
  theme_minimal() + theme(axis.line = element_line(linewidth = 0.1, colour = "black"),
                          plot.title = element_text(hjust = 0.5, size = 18), axis.text = element_text(size=16, hjust = 1),
                          axis.title = element_text(size = 18), legend.text = element_text(size=14),
                          legend.title = element_text(size=16)
  ) + scale_color_manual(values=c("darkgrey", "#7a7af0", "#CC6666")) +
  ggtitle("ITGAV") +
  xlab("Western blot (Log2[AU])") + ylab("Proteomics (Log2)") +
  annotate(geom="text",x=-1,y=21.5, label="r = 0.25", color = "black",size=6) +
  annotate(geom = "text",x=-1,y=21.2,label="p = 5.5e-04",size=6)




#### Figure_3 ####
# Clinical data will not be provided
clinical_data2_pre <- clinical_data2[which(grepl("Pre",clinical_data2$New_ID)),]
signif_proteins <- T2DPre_LeanPre_75filt$ID[which(T2DPre_LeanPre_75filt$adj.P.Val < 0.1)]
signif_proteins <- rownames(T2DPre_LeanPre_75filt)[which(T2DPre_LeanPre_75filt$adj.P.Val < 0.1)]

Signi_exprs_list <- Exprs_adipose_clean[which(rownames(Exprs_adipose_clean) %in% signif_proteins),]
Signi_exprs_list_pre <- Signi_exprs_list[,-c(which(grepl("Post",colnames(Signi_exprs_list))))]

half_data <- data.frame(cbind(rownames(Signi_exprs_list_pre),Signi_exprs_list_pre))
colnames(half_data)[1] <- "Gene"

GeneList <- unique(half_data$Gene)
final_list_pre <- data.frame()
comb_final_cor <- data.frame(nrow=48)
half_data <- gather(half_data, ID, Intensity, -Gene)
half_data$Intensity <- as.numeric(half_data$Intensity)

sub_clinical_data_pre <- clinical_data2_pre[,c("GIR1","HbA1c1","FM1","TG1","VO2max1")]
for(k in 1:length(sub_clinical_data_pre)){
  if(shapiro.test(log2(sub_clinical_data_pre[,k]))[2]<0.05){
    for(i in 1:length(GeneList)){
      first_list <- cor.test(log2(sub_clinical_data_pre[,k]),as.numeric(half_data$Intensity[which(half_data$Gene == GeneList[i])]),method = "kendall")
      final_list_pre[i,1] <- first_list$p.value
      final_list_pre[i,2] <- first_list$estimate
      final_list_pre[i,3] <- GeneList[i]
    }}
  else{
    for(i in 1:length(GeneList)){
      first_list <- cor.test(log2(sub_clinical_data_pre[,k]),as.numeric(half_data$Intensity[which(half_data$Gene == GeneList[i])]),method = "pearson")
      final_list_pre[i,1] <- first_list$p.value
      final_list_pre[i,2] <- first_list$estimate
      final_list_pre[i,3] <- GeneList[i]
    }}
  final_list_pre[,4] <- p.adjust(final_list_pre[,1])
  final_list_pre <- final_list_pre[c(3,1,4,2)]
  colnames(final_list_pre)[1] <- c("Gene")
  colnames(final_list_pre)[2:4] <- paste0(colnames(sub_clinical_data_pre)[k],c("p.value","adj.p.val","coef"))
  comb_final_cor <- cbind(comb_final_cor,final_list_pre)
  
}

## Export for Table_S2E
comb_final_cor <- comb_final_cor[,-c(1,6,10,14,18)]
write.table(comb_final_cor,"Baseline_correlation_results.txt", append = FALSE, sep = "\t", dec = ".",
            row.names = F, col.names = TRUE)


### Figure_3A
comb_final_cor <- read.delim("Baseline_correlation_results.txt")
p_val_data <- comb_final_cor[,which(grepl("value|Gene",colnames(comb_final_cor)))]

long_p_val_data <- gather(p_val_data, Parameter,p.val,-Gene)
long_p_val_data$p.val <- -log10(as.numeric(long_p_val_data$p.val))
long_p_val_data$label <- NA
long_p_val_data$sign <- "NO"
long_p_val_data$sign[which(long_p_val_data$p.val> -log10(0.0012))] <- "YES"
long_p_val_data$label[long_p_val_data$sign != "NO"] <- long_p_val_data$Gene[long_p_val_data$sign != "NO"]


ggplot(long_p_val_data,aes(x=factor(Parameter),y=p.val)) + geom_point(aes(colour=p.val), size=2) + 
  geom_hline(yintercept=-log10(0.0012), linetype="solid", 
             color = "black", linewidth=0.5) + geom_text_repel(aes(label=label)) + theme_minimal() +
  xlab("") + ylab("-Log10(p.value)") + scale_color_gradient(high = "#CC6666", low = "darkgrey") +
  theme(axis.line = element_line(linewidth = 0.1, colour = "black"),
        plot.title = element_text(hjust = 0.5, size = 18), axis.text = element_text(size=16, hjust = 1),
        axis.title = element_text(size = 18), legend.text = element_text(size=14),
        legend.title = element_text(size=16))



#### Figure_3B-D, Individual clinical data is not provided ####
Baseline_samples <- Exprs_adipose_clean[,grepl("Pre", colnames(Exprs_adipose_clean))]

#ANXA3 correlation - Figure_3B
ANXA3_barplots <- as.data.frame(Baseline_samples[which(rownames(Baseline_samples) == "ANXA3"),])
ANXA3_barplots$Sample_ID <- rownames(ANXA3_barplots)
colnames(ANXA3_barplots)[1] <- "Intensity"
ANXA3_barplots$Group <- factor(word(ANXA3_barplots$Sample_ID,1, sep="_"), levels = c("Lean","Obese","T2D"))

ANXA3_cor_matrix <- cbind(ANXA3_barplots,clinical_data2_pre$GIR1)
colnames(ANXA3_cor_matrix)[4] <- "GIR"

cor.test(ANXA3_barplots$Intensity,log2(clinical_data2_pre$GIR1), method = "kendall")

ggplot(ANXA3_cor_matrix, aes(x=GIR, y=Intensity)) + geom_point(aes(color = Group)) +
  geom_smooth(method = "lm", se = TRUE, colour = "darkblue", size = 0.75, alpha=0.15, fill="#808096") +
  theme_minimal() + theme(axis.line = element_line(linewidth = 0.1, colour = "black"),
                          plot.title = element_text(hjust = 0.5, size = 18), axis.text = element_text(size=16, hjust = 1),
                          axis.title = element_text(size = 18), legend.text = element_text(size=14),
                          legend.title = element_text(size=16)
  ) + scale_color_manual(values=c("darkgrey", "#7a7af0", "#CC6666")) +
  ggtitle("ANXA3") +
  xlab("Glucose infusion rate (mg/min/m2)") + ylab("Log2 intensity") +
  annotate(geom="text",x=400,y=21.6, label="r = 0.45", color = "black",size=6) +
  annotate(geom = "text",x=400,y=21.3,label="p = 7.63e-06",size=6)




#CD14 correlation - Figure_3C
CD14_barplots <- as.data.frame(Baseline_samples[which(rownames(Baseline_samples) == "CD14"),])
CD14_barplots$Sample_ID <- rownames(CD14_barplots)
colnames(CD14_barplots)[1] <- "Intensity"
CD14_barplots$Group <- factor(word(CD14_barplots$Sample_ID,1, sep="_"), levels = c("Lean","Obese","T2D"))

CD14_cor_matrix <- cbind(CD14_barplots,clinical_data2_pre$GIR1)
colnames(CD14_cor_matrix)[4] <- "GIR"

cor.test(CD14_barplots$Intensity,log2(clinical_data2_pre$GIR1), method = "kendall")

ggplot(CD14_cor_matrix, aes(x=GIR, y=Intensity)) + geom_point(aes(color = Group)) +
  geom_smooth(method = "lm", se = TRUE, colour = "darkblue", size = 0.75, alpha=0.15, fill="#808096") +
  theme_minimal() + theme(axis.line = element_line(linewidth = 0.1, colour = "black"),
                          plot.title = element_text(hjust = 0.5, size = 18), axis.text = element_text(size=16, hjust = 1),
                          axis.title = element_text(size = 18), legend.text = element_text(size=14),
                          legend.title = element_text(size=16)
  ) + scale_color_manual(values=c("darkgrey", "#7a7af0", "#CC6666")) +
  ggtitle("CD14") +
  xlab("Glucose infusion rate (mg/min/m2)") + ylab("Log2 intensity") +
  annotate(geom="text",x=460,y=21, label="r = -0.39", color = "black",size=6) +
  annotate(geom = "text",x=460,y=20.7,label="p = 1.1e-04",size=6)


#MIF correlation - Figure_3D
MIF_barplots <- as.data.frame(Baseline_samples[which(rownames(Baseline_samples) == "MIF"),])
MIF_barplots$Sample_ID <- rownames(MIF_barplots)
colnames(MIF_barplots)[1] <- "Intensity"
MIF_barplots$Group <- factor(word(MIF_barplots$Sample_ID,1, sep="_"), levels = c("Lean","Obese","T2D"))

MIF_cor_matrix <- cbind(MIF_barplots,clinical_data2_pre$HbA1c1)
colnames(MIF_cor_matrix)[4] <- "HbA1c1"
cor.test(MIF_barplots$Intensity,log2(clinical_data2_pre$HbA1c1), method = "kendall")

ggplot(MIF_cor_matrix, aes(x=HbA1c1, y=Intensity)) + geom_point(aes(color = Group)) +
  geom_smooth(method = "lm", se = TRUE, colour = "darkblue", size = 0.75, alpha=0.15, fill="#808096") +
  theme_minimal() + theme(axis.line = element_line(linewidth = 0.1, colour = "black"),
                          plot.title = element_text(hjust = 0.5, size = 18), axis.text = element_text(size=16, hjust = 1),
                          axis.title = element_text(size = 18), legend.text = element_text(size=14),
                          legend.title = element_text(size=16)
  ) + scale_color_manual(values=c("darkgrey", "#7a7af0", "#CC6666")) +
  ggtitle("MIF") + xlab("HbA1c (mmol/mol)") + ylab("Log2 intensity") +
  annotate(geom="text",x=70,y=23.9, label="r = 0.36", color = "black",size=6) +
  annotate(geom = "text",x=70,y=23.7,label="p = 5.3e-04",size=6)



#Figure_S2D
Exprs_adipose_whole_clean <- read.delim("Exprs_adipose_whole_clean.txt",header = TRUE)
rownames(Exprs_adipose_whole_clean) <- Exprs_adipose_whole_clean$ID
Exprs_adipose_whole_clean <- as.matrix(Exprs_adipose_whole_clean[,-c(92)])
Baseline_samples_whole_proteome <- Exprs_adipose_whole_clean[,grepl("Pre", colnames(Exprs_adipose_whole_clean))]

LEP_barplots <- as.data.frame(Baseline_samples_whole_proteome[which(rownames(Baseline_samples_whole_proteome) == "LEP"),])
LEP_barplots$Sample_ID <- rownames(LEP_barplots)
colnames(LEP_barplots)[1] <- "Intensity"
LEP_barplots$Group <- factor(word(LEP_barplots$Sample_ID,1, sep="_"), levels = c("Lean","Obese","T2D"))

LEP_cor_matrix <- cbind(LEP_barplots,clinical_data2_pre$FM1)
colnames(LEP_cor_matrix)[4] <- "FM"

cor.test(LEP_cor_matrix$Intensity,log2(clinical_data2_pre$FM), method = "pearson")
ggplot(LEP_cor_matrix, aes(x=FM, y=Intensity)) + geom_point(aes(color = Group)) +
  geom_smooth(method = "lm", se = TRUE, colour = "darkblue", size = 0.75, alpha=0.15, fill="#808096") +
  theme_minimal() + theme(axis.line = element_line(linewidth = 0.1, colour = "black"),
                          plot.title = element_text(hjust = 0.5, size = 18), axis.text = element_text(size=16, hjust = 1),
                          axis.title = element_text(size = 18), legend.text = element_text(size=14),
                          legend.title = element_text(size=16)
  ) + scale_color_manual(values=c("darkgrey", "#7a7af0", "#CC6666")) +
  ggtitle("Leptin") +
  xlab("Fat mass (kg)") + ylab("Log2 intensity") +
  annotate(geom="text",x=43,y=17.9, label="r = 0.59", color = "black",size=6) +
  annotate(geom = "text",x=43,y=17.6,label="p = 1.8e-04",size=6)






#### Boxplots - Figure_4D-E ####
#Ferritin boxplot of pre/post training abundances
#Figure_4D
Ferritin_samples <- as.data.frame(Exprs_adipose_clean[which(rownames(Exprs_adipose_clean) == "FTL"),])
Ferritin_samples$Sample_ID <- rownames(Ferritin_samples)
colnames(Ferritin_samples)[1] <- "Intensity"
Ferritin_samples$Group <- factor(word(Ferritin_samples$Sample_ID,1, sep="_"), levels = c("Lean","Obese","T2D"))
Ferritin_samples$training <- factor(word(Ferritin_samples$Sample_ID,2, sep="_"), levels = c("Pre","Post"))
Ferritin_samples$paired_ID <- factor(word(Ferritin_samples$Sample_ID,4, sep="_"))

ggplot(Ferritin_samples, aes(y=Intensity, x=training, fill=training)) + geom_boxplot(alpha = .2) +
  scale_fill_manual(values=c("darkgrey", "#2e136e")) +
  geom_point(size=0.5) +
  geom_line(aes(group = paired_ID), size= 0.05) +
  facet_wrap(~ Group, switch = "x",nrow=1) +
  scale_x_discrete("") +
  theme_minimal() +
  theme(legend.position = "top",
        plot.title = element_text(hjust = 0.5, size = 18),
        panel.grid = element_blank(),
        legend.title = element_text(size=16),
        legend.text = element_text(size=16),
        strip.text.x = element_text(size = 18),
        axis.title = element_text(size = 18),
        axis.text.x = element_blank(),
        axis.line.x = element_line(size=.5),
        axis.text.y = element_text(size = 18),
        axis.line.y = element_line(size = .5)) +
  xlab("") + ylab("Intensity (Log2)") + ggtitle("FTL")

#Figure_4E
FTH1_samples <- as.data.frame(Exprs_adipose_clean[which(rownames(Exprs_adipose_clean) == "FTH1"),])
FTH1_samples$Sample_ID <- rownames(FTH1_samples)
colnames(FTH1_samples)[1] <- "Intensity"
FTH1_samples$Group <- factor(word(FTH1_samples$Sample_ID,1, sep="_"), levels = c("Lean","Obese","T2D"))
FTH1_samples$training <- factor(word(FTH1_samples$Sample_ID,2, sep="_"), levels = c("Pre","Post"))
FTH1_samples$paired_ID <- factor(word(FTH1_samples$Sample_ID,4, sep="_"))

ggplot(FTH1_samples, aes(y=Intensity, x=training, fill=training)) + geom_boxplot(alpha = .2) +
  scale_fill_manual(values=c("darkgrey", "#2e136e")) +
  geom_point(size=0.5) +
  geom_line(aes(group = paired_ID), size= 0.05) +
  facet_wrap(~ Group, switch = "x",nrow=1) +
  scale_x_discrete("") +
  theme_minimal() +
  theme(legend.position = "top",
        plot.title = element_text(hjust = 0.5, size = 18),
        panel.grid = element_blank(),
        legend.title = element_text(size=16),
        legend.text = element_text(size=16),
        strip.text.x = element_text(size = 18),
        axis.title = element_text(size = 18),
        axis.text.x = element_blank(),
        axis.line.x = element_line(size=.5),
        axis.text.y = element_text(size = 18),
        axis.line.y = element_line(size = .5)) +
  xlab("") + ylab("Intensity (Log2)") + ggtitle("FTH1")


#### Figure_4F-G & Figure_S3_B ####
WB_data <- read.delim("WB_HIIT_paired_data.txt", dec = ",")

long_WB_data <- gather(WB_data, Protein,Intensity,-c(Random_ID,Group,Intervention))
long_WB_data$Intensity <- as.numeric(long_WB_data$Intensity)
long_WB_data$Group <- factor(long_WB_data$Group, levels = c("Lean","Obese","T2D"))
long_WB_data$Intervention<- factor(long_WB_data$Intervention, levels = c("Pre","Post"))
long_WB_data$Random_ID <- factor(long_WB_data$Random_ID)
FTL <- long_WB_data[which(long_WB_data$Protein == "FTL.Normalised"),]


ggplot(FTL, aes(y=Intensity, x=Intervention, fill=Intervention)) + geom_boxplot(alpha = .2) +
  scale_fill_manual(values=c("darkgrey", "#2e136e")) +
  geom_point(size=0.5) + 
  geom_line(aes(group = Random_ID), size= 0.05) +
  facet_wrap(~ Group, switch = "x",nrow=1) +
  scale_x_discrete("") +
  theme_minimal() +
  theme(legend.position = "top",
        plot.title = element_text(hjust = 0.5, size = 18),
        panel.grid = element_blank(),
        legend.title = element_text(size=16),
        legend.text = element_text(size=16),
        strip.text.x = element_text(size = 18),
        axis.title = element_text(size = 18),
        axis.text.x = element_blank(),
        axis.line.x = element_line(size=.5),
        axis.text.y = element_text(size = 18),
        axis.line.y = element_line(size = .5)) + 
  xlab("") + ylab("Fold change from Lean Pre") + ggtitle("FTL")

FTH1 <- long_WB_data[which(long_WB_data$Protein == "FTH1.Normalised"),]

ggplot(FTH1, aes(y=Intensity, x=Intervention, fill=Intervention)) + geom_boxplot(alpha = .2) +
  scale_fill_manual(values=c("darkgrey", "#2e136e")) +
  geom_point(size=0.5) + 
  geom_line(aes(group = Random_ID), size= 0.05) +
  facet_wrap(~ Group, switch = "x",nrow=1) +
  scale_x_discrete("") +
  theme_minimal() +
  theme(legend.position = "top",
        plot.title = element_text(hjust = 0.5, size = 18),
        panel.grid = element_blank(),
        legend.title = element_text(size=16),
        legend.text = element_text(size=16),
        strip.text.x = element_text(size = 18),
        axis.title = element_text(size = 18),
        axis.text.x = element_blank(),
        axis.line.x = element_line(size=.5),
        axis.text.y = element_text(size = 18),
        axis.line.y = element_line(size = .5)) + 
  xlab("") + ylab("Fold change from Lean Pre") + ggtitle("FTH1")

TfR <- long_WB_data[which(long_WB_data$Protein == "TfR.Normalised"),]

ggplot(TfR, aes(y=Intensity, x=Intervention, fill=Intervention)) + geom_boxplot(alpha = .2) +
  scale_fill_manual(values=c("darkgrey", "#2e136e")) +
  geom_point(size=0.5) + 
  geom_line(aes(group = Random_ID), size= 0.05) +
  facet_wrap(~ Group, switch = "x",nrow=1) +
  scale_x_discrete("") +
  theme_minimal() +
  theme(legend.position = "top",
        plot.title = element_text(hjust = 0.5, size = 18),
        panel.grid = element_blank(),
        legend.title = element_text(size=16),
        legend.text = element_text(size=16),
        strip.text.x = element_text(size = 18),
        axis.title = element_text(size = 18),
        axis.text.x = element_blank(),
        axis.line.x = element_line(size=.5),
        axis.text.y = element_text(size = 18),
        axis.line.y = element_line(size = .5)) + 
  xlab("") + ylab("Fold change from Lean Pre") + ggtitle("TfR1")




























#### Figure_5A-B ####
#Figure_5A
FTL_cor <- Exprs_adipose_clean[which(rownames(Exprs_adipose_clean) == "FTL"),]
FTL_cor_matrix <- cbind(FTL_cor,clinical_data2) # Clinical data is not provided

FTL_GIR_cor <- rmcorr(
  ID,
  FTL_cor,
  log2(GIR1),
  FTL_cor_matrix,
  CI.level = 0.95,
  CIs = c("analytic", "bootstrap"),
  nreps = 100,
  bstrap.out = F
)


ggplot(FTL_cor_matrix, aes(x=GIR1, y=FTL_cor)) + geom_point(aes(color=Condition),alpha=0.4, size=3) + geom_line(aes(group = ID),alpha=0.6, size=0.1) +
  theme_minimal() + theme(axis.line = element_line(linewidth = 0.1, colour = "black"),
                          plot.title = element_text(hjust = 0.5), axis.text = element_text(size=12, hjust = 1),
                          axis.title = element_text(size = 12), legend.text = element_text(size=12),
                          legend.title = element_text(size=14)
  ) + scale_color_manual(values=c("#2e136e","darkgrey")) +
  ggtitle("Adipose tissue FTL") +
  xlab("Glucose infusion rate (mg/min/m2)") + ylab("Log2 intensity") + 
  annotate(geom="text",x=600,y=23.3, label="r = 0.69", color = "black",size=4) +
  annotate(geom = "text",x=600,y=23,label="p = 1.81e-07",size=4)



#Figure_5B
FTH1_cor <- Exprs_adipose_clean[which(rownames(Exprs_adipose_clean) == "FTH1"),]
FTH1_cor_matrix <- cbind(FTL_cor,clinical_data2) #Clinical data is not provided

ggplot(FTH1_cor_matrix, aes(x=GIR1, y=FTH1_cor)) + geom_point(aes(color=Condition),alpha=0.4, size=3) + geom_line(aes(group = ID),alpha=0.6, size=0.1) +
  theme_minimal() + theme(axis.line = element_line(linewidth = 0.1, colour = "black"),
                          plot.title = element_text(hjust = 0.5), axis.text = element_text(size=12, hjust = 1),
                          axis.title = element_text(size = 12), legend.text = element_text(size=12),
                          legend.title = element_text(size=14)
  ) + scale_color_manual(values=c("#2e136e","darkgrey")) +
  ggtitle("Adipose tissue FTH1") +
  xlab("Glucose infusion rate (mg/min/m2)") + ylab("Log2 intensity") + 
  annotate(geom="text",x=600,y=22.3, label="r = 0.68", color = "black",size=4) +
  annotate(geom = "text",x=600,y=22,label="p = 4.00e-07",size=4)



#### Figure_5C-E & S3A ####

### updated correlation analysis
setwd("C:/Users/mct330/Desktop/Work/KURT_HIIT_ADIPOSE_STUDY")
library("rmcorr")
library("dplyr")
library("ggpubr")


adipose_data <- read.delim("updated_Repeated_measure_serum_WB_correlation.txt", dec = ",") #This data contains clinical data and is not provided
Before_intervention <- data.frame(adipose_data$ID,adipose_data$Group,adipose_data$Day_2_Transferin_serum,adipose_data$Day_2_Jern_serum,
                                  adipose_data$Day_2_Ferritin_serum,adipose_data$Day_2_FTL_wb,adipose_data$Day_2_FTH1_wb,
                                  adipose_data$Day_2_FTL_proteomics,adipose_data$Day_2_FTH1_proteomics,adipose_data$Day_2_GIR1)
Before_intervention$Training <- "Pre"
After_intervention <- data.frame(adipose_data$ID,adipose_data$Group,adipose_data$Day_4_Transferin_serum,adipose_data$Day_4_Jern_serum,
                                 adipose_data$Day_4_Ferritin_serum,adipose_data$Day_4_FTL_wb,adipose_data$Day_4_FTH1_wb,
                                 adipose_data$Day_4_FTL_proteomics,adipose_data$Day_4_FTH1_proteomics,adipose_data$Day_4_GIR2)

After_intervention$Training <- "Post"

colnames(Before_intervention) <- c("ID","Group","Transferin_serum","Jern_serum","Ferritin_serum","FTL_WB","FTH_WB","FTL_proteomics",
                                   "FTH_proteomics","GIR","Training")
colnames(After_intervention) <- c("ID","Group","Transferin_serum","Jern_serum","Ferritin_serum","FTL_WB","FTH_WB","FTL_proteomics",
                                  "FTH_proteomics","GIR","Training")

combined_intervention <- rbind(Before_intervention, After_intervention)
combined_intervention$Training <- factor(combined_intervention$Training)
combined_intervention[combined_intervention$Group == 1,]$Group <- "Lean" 
combined_intervention[combined_intervention$Group == 2,]$Group <- "Obese" 
combined_intervention[combined_intervention$Group == 3,]$Group <- "T2D" 

ggplot(combined_intervention, aes(x=factor(Training, levels= c("Pre","Post")), y=Jern_serum, fill=Training)) + geom_boxplot(alpha = .2) +
  scale_fill_manual(values=c("darkgrey", "#2e136e"), breaks = c("Pre","Post")) +
  geom_point(size=0.25) + 
  geom_line(aes(group = ID), size= 0.05) +
  facet_wrap(~ Group, switch = "x",nrow=1) +
  scale_x_discrete("") +
  theme_minimal() + 
  theme(legend.position = "top",
        plot.title = element_text(hjust = 0.5, size = 18),
        panel.grid = element_blank(),
        legend.title = element_text(size=14),
        legend.text = element_text(size=14),
        strip.text.x = element_text(size = 12),
        axis.title = element_text(size = 18),
        axis.text.x = element_blank(),
        axis.line.x = element_line(size=.5),
        axis.text.y = element_text(size = 14),
        axis.line.y = element_line(size = .5)) + 
  xlab("") + ylab("(µmol/L)") + ggtitle("Ferritin serum")

ggplot(combined_intervention, aes(x=factor(Training, levels= c("Pre","Post")), y=Transferin_serum, fill=Training)) + geom_boxplot(alpha = .2) +
  scale_fill_manual(values=c("darkgrey", "#2e136e"), breaks = c("Pre","Post")) +
  geom_point(size=0.25) + 
  geom_line(aes(group = ID), size= 0.05) +
  facet_wrap(~ Group, switch = "x",nrow=1) +
  scale_x_discrete("") +
  theme_minimal() + 
  theme(legend.position = "top",
        plot.title = element_text(hjust = 0.5, size = 18),
        panel.grid = element_blank(),
        legend.title = element_text(size=14),
        legend.text = element_text(size=14),
        strip.text.x = element_text(size = 12),
        axis.title = element_text(size = 18),
        axis.text.x = element_blank(),
        axis.line.x = element_line(size=.5),
        axis.text.y = element_text(size = 14),
        axis.line.y = element_line(size = .5)) + 
  xlab("") + ylab("(µmol/L)") + ggtitle("Ferritin serum")


ggplot(combined_intervention, aes(x=factor(Training, levels= c("Pre","Post")), y=Ferritin_serum, fill=Training)) + geom_boxplot(alpha = .2) +
  scale_fill_manual(values=c("darkgrey", "#2e136e"), breaks = c("Pre","Post")) +
  geom_point(size=0.25) + 
  geom_line(aes(group = ID), size= 0.05) +
  facet_wrap(~ Group, switch = "x",nrow=1) +
  scale_x_discrete("") +
  theme_minimal() + 
  theme(legend.position = "top",
        plot.title = element_text(hjust = 0.5, size = 18),
        panel.grid = element_blank(),
        legend.title = element_text(size=14),
        legend.text = element_text(size=14),
        strip.text.x = element_text(size = 12),
        axis.title = element_text(size = 18),
        axis.text.x = element_blank(),
        axis.line.x = element_line(size=.5),
        axis.text.y = element_text(size = 14),
        axis.line.y = element_line(size = .5)) + 
  xlab("") + ylab("(µg/L)") + ggtitle("Ferritin serum")


#Figure_S3A
ggplot(combined_intervention, aes(x=FTL_proteomics, y=log2(FTL_WB))) + geom_point(aes(color=Training),alpha=0.4, size=3) + geom_line(aes(group = ID),alpha=0.6,size=0.1) +
  theme_bw() + theme(plot.title = element_text(hjust = 0.5), axis.text = element_text(size=12, hjust = 1),
                     axis.title = element_text(size = 16), legend.text = element_text(size=12),
                     legend.title = element_text(size=14)
  ) + scale_color_manual(values=c("#2e136e", "darkgrey")) +
  annotate(geom="text",x=3,y=-7.5, label="r = 0.83", color = "black",size=5) +
  annotate(geom = "text",x=3,y=-8,label="p = 5.01e-12",size=5) +
  xlab("FTL proteomics (Log2)") + ylab("FTL Western blot (Log2[AU]) ")


#### Figure_5F ####
new_clinical_data <- combined_intervention[-c(17),] # FTL/FTH1 not detected in this sample
new_clinical_data <- new_clinical_data[which(new_clinical_data$FTL_WB > 0),]
new_clinical_data <- new_clinical_data[order(new_clinical_data$ID),]
new_clinical_data$ID <- sub("-","",new_clinical_data$ID)

prot_exprs_data <- data.frame(t(Exprs_adipose_clean[c(which(rownames(Exprs_adipose_clean) == "FTL"),which(rownames(Exprs_adipose_clean) == "FTH1")),]))
prot_exprs_data$ID <- factor(clinical_data2$ID)
prot_exprs_data <- prot_exprs_data[order(prot_exprs_data$ID),]

Clin_correlation_matrix <- cbind(new_clinical_data,prot_exprs_data[,1:2])
Clin_correlation_matrix <- Clin_correlation_matrix[,-c(6:9)]
Clin_correlation_matrix <- Clin_correlation_matrix[,c(1,2,7,3:6,8,9)]


Long_clin_cor_matrix <- gather(Clin_correlation_matrix, Training, Group, Transferin_serum:FTH1)
colnames(Long_clin_cor_matrix) <- c("ID","Parameter","Value")

variables <- c("Transferin_serum","Jern_serum","Ferritin_serum","GIR","FTL","FTH1")

output_clin <- matrix(nrow=91, ncol = 4)
colnames(output_clin)<- variables[1:4]
for(i in 1:4){
  output_clin[,i] <- log2(as.numeric(Long_clin_cor_matrix[Long_clin_cor_matrix$Parameter == variables[i],]$Value))
}

Clin_correlation_matrix[,c(4:7)] <- output_clin
Long_clin_cor_matrix <- gather(Clin_correlation_matrix, Training, Group, Transferin_serum:FTH1)
colnames(Long_clin_cor_matrix) <- c("ID","Parameter","Value")
Long_clin_cor_matrix$ID <- as.factor(Long_clin_cor_matrix$ID)

empty_list <- matrix(nrow=6, ncol=6)
length(variables)
new_matrix <- matrix()



### for r value
for(i in 1:6){
  test_list <- data.frame(Long_clin_cor_matrix[Long_clin_cor_matrix$Parameter == variables[i],c(1,3)],
                          Long_clin_cor_matrix[Long_clin_cor_matrix$Parameter == variables[1],3])
  empty_list[1,i] <- rmcorr(ID,test_list[,2],test_list[,3],test_list)$r
  test_list <- data.frame(Long_clin_cor_matrix[Long_clin_cor_matrix$Parameter == variables[i],c(1,3)],
                          Long_clin_cor_matrix[Long_clin_cor_matrix$Parameter == variables[2],3])
  empty_list[2,i] <- rmcorr(ID,test_list[,2],test_list[,3],test_list)$r
  test_list <- data.frame(Long_clin_cor_matrix[Long_clin_cor_matrix$Parameter == variables[i],c(1,3)],
                          Long_clin_cor_matrix[Long_clin_cor_matrix$Parameter == variables[3],3])
  empty_list[3,i] <- rmcorr(ID,test_list[,2],test_list[,3],test_list)$r
  test_list <- data.frame(Long_clin_cor_matrix[Long_clin_cor_matrix$Parameter == variables[i],c(1,3)],
                          Long_clin_cor_matrix[Long_clin_cor_matrix$Parameter == variables[4],3])
  empty_list[4,i] <- rmcorr(ID,test_list[,2],test_list[,3],test_list)$r
  test_list <- data.frame(Long_clin_cor_matrix[Long_clin_cor_matrix$Parameter == variables[i],c(1,3)],
                          Long_clin_cor_matrix[Long_clin_cor_matrix$Parameter == variables[5],3])
  empty_list[5,i] <- rmcorr(ID,test_list[,2],test_list[,3],test_list)$r
  test_list <- data.frame(Long_clin_cor_matrix[Long_clin_cor_matrix$Parameter == variables[i],c(1,3)],
                          Long_clin_cor_matrix[Long_clin_cor_matrix$Parameter == variables[6],3])
  empty_list[6,i] <- rmcorr(ID,test_list[,2],test_list[,3],test_list)$r
}

colnames(empty_list) <- variables
rownames(cormat) <- variables
cormat <- round(empty_list,2)

get_upper_tri <- function(cormat){
  cormat[lower.tri(cormat)]<- NA
  return(cormat)
}

reorder_cormat <- function(cormat){
  # Use correlation between variables as distance
  dd <- as.dist((1-cormat)/2)
  hc <- hclust(dd)
  cormat <-cormat[hc$order, hc$order]
}

# Reorder the correlation matrix
cormat <- reorder_cormat(cormat)
upper_tri <- get_upper_tri(cormat)
# Melt the correlation matrix
melted_cormat <- melt(upper_tri, na.rm = TRUE)

ggheatmap <- ggplot(melted_cormat, aes(Var2, Var1, fill = value))+
  geom_tile(color = "white")+
  theme_minimal()+ # minimal theme
  theme(axis.text.x = element_text(angle = 45, vjust = 1,
                                   size = 12, hjust = 1))+
  coord_fixed()

ggheatmap +
  geom_text(aes(Var2, Var1, label = value), color = "black", size = 3) +
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    axis.text.x = element_text(angle = 45, size = 12, hjust = 1),
    axis.text.y = element_text(size = 12, hjust = 1),
    panel.grid.major = element_blank(),
    panel.border = element_blank(),
    panel.background = element_blank(),
    axis.ticks = element_blank(),
    legend.justification = c(1, 0),
    legend.position = c(0.6, 0.7),
    legend.direction = "horizontal")+
  guides(fill = guide_colorbar(barwidth = 7, barheight = 1,
                               title.position = "top", title.hjust = 0.5))+
  scale_y_discrete(position = "right") +
  scale_fill_gradient2(low = "darkblue", high = "darkgray", mid = "white",
                       midpoint = 0, limit = c(-1,1), space = "Lab",
                       name="r correlation")




#### Go enrichment analysis ####
T2DPre_LeanPre_gene_list <- T2DPre_LeanPre_75filt$logFC
names(T2DPre_LeanPre_gene_list) <- T2DPre_LeanPre_75filt$ID
T2DPre_LeanPre_gene_list = sort(T2DPre_LeanPre_gene_list, decreasing = TRUE)

gse <- gseGO(geneList=T2DPre_LeanPre_gene_list,
             ont ="CC",
             keyType = "SYMBOL",
             nPerm = 10000,
             minGSSize = 3,
             maxGSSize = 800,
             pvalueCutoff = 0.1,
             verbose = TRUE,
             OrgDb = org.Hs.eg.db,
             pAdjustMethod = "BH",
             seed = T)

T2D_Lean_CC <- gse@result

gse <- gseGO(geneList=T2DPre_LeanPre_gene_list,
             ont ="BP",
             keyType = "SYMBOL",
             nPerm = 10000,
             minGSSize = 3,
             maxGSSize = 800,
             pvalueCutoff = 0.1,
             verbose = TRUE,
             OrgDb = org.Hs.eg.db,
             pAdjustMethod = "BH",
             seed = T)

T2D_Lean_BP <- gse@result


#### T2D-Obese
T2DPre_ObesePre_gene_list <- T2DPre_ObesePre_75filt$logFC
names(T2DPre_ObesePre_gene_list) <- T2DPre_ObesePre_75filt$ID
T2DPre_ObesePre_gene_list = sort(T2DPre_ObesePre_gene_list, decreasing = TRUE)

gse <- gseGO(geneList=T2DPre_ObesePre_gene_list,
             ont ="CC",
             keyType = "SYMBOL",
             nPerm = 10000,
             minGSSize = 3,
             maxGSSize = 800,
             pvalueCutoff = 0.1,
             verbose = TRUE,
             OrgDb = org.Hs.eg.db,
             pAdjustMethod = "BH",
             seed = T)

T2D_Obese_CC <- gse@result

gse <- gseGO(geneList=T2DPre_ObesePre_list,
             ont ="BP",
             keyType = "SYMBOL",
             nPerm = 10000,
             minGSSize = 3,
             maxGSSize = 800,
             pvalueCutoff = 0.1,
             verbose = TRUE,
             OrgDb = org.Hs.eg.db,
             pAdjustMethod = "BH",
             seed = T)
T2D_Obese_BP <- gse@result


### Obese-Lean
ObesePre_LeanPre_gene_list <- ObesePre_LeanPre_75filt$logFC
names(ObesePre_LeanPre_gene_list) <- ObesePre_LeanPre_75filt$ID
ObesePre_LeanPre_gene_list = sort(ObesePre_LeanPre_gene_list, decreasing = TRUE)

gse <- gseGO(geneList=ObesePre_LeanPre_gene_list,
             ont ="CC",
             keyType = "SYMBOL",
             nPerm = 10000,
             minGSSize = 3,
             maxGSSize = 800,
             pvalueCutoff = 0.1,
             verbose = TRUE,
             OrgDb = org.Hs.eg.db,
             pAdjustMethod = "BH",
             seed = T)

Obese_Lean_CC <- gse@result

gse <- gseGO(geneList=ObesePre_LeanPre_gene_list,
             ont ="BP",
             keyType = "SYMBOL",
             nPerm = 10000,
             minGSSize = 3,
             maxGSSize = 800,
             pvalueCutoff = 0.1,
             verbose = TRUE,
             OrgDb = org.Hs.eg.db,
             pAdjustMethod = "BH",
             seed = T)

Obese_Lean_BP <- gse@result

#### Obesity-driven GO-terms
#Figure_2E
joined_obese_BP <- inner_join(T2D_Lean_BP,Obese_Lean_BP, by = "ID")
joined_obese_BP$GO <- "BP"
joined_obese <- inner_join(T2D_Lean_CC, Obese_Lean_CC, by = "ID")
joined_obese$GO <- "CC"
combined_Obese <- rbind(joined_obese,joined_obese_BP)

comb_obese_selected <- c("cell morphogenesis involved in differentiation","ERK1 and ERK2 cascade",
                         "actin filament","side of membrane","cell morphogenesis",
                         "vasoconstriction","positive regulation of heterotypic cell-cell adhesion",
                         "leading edge membrane","response to wounding","cell-matrix adhesion",
                         "cellular amino acid metabolic process","alpha-amino acid metabolic process",
                         "cell leading edge","platelet aggregation")

ggplot(combined_Obese, aes(x=enrichmentScore.y, y=enrichmentScore.x, label = ifelse(Description.x %in% comb_obese_selected, Description.x,""))) + geom_point(aes(size = setSize.x, colour = p.adjust.x)) +
  geom_text_repel(size=2.5, max.overlaps = Inf) + theme_bw() + facet_grid(cols = vars(GO)) +
  scale_color_gradient(low = "#CC6666", high = "Dark grey") +
  scale_y_break(c(-0.4, 0.25), scale=2.5, space = 0.05) +
  theme(panel.spacing = unit(0.35, "cm")) +
  ylab("GeneRatio (Obese vs Lean)") + xlab("GeneRatio (T2D vs Lean)")

#Export for Table_S2G
write.table(combined_Obese,"Obesity_driven_GOterms.txt", append = FALSE, sep = "\t", dec = ".",
            row.names = F, col.names = TRUE)





#### Diabetes-driven GO-terms
#Figure_2D
joined_T2D_CC <- inner_join(T2D_Lean_CC, T2D_Obese_CC, by = "ID")
joined_T2D_BP <- inner_join(T2D_Lean_BP, T2D_Obese_BP, by = "ID")

joined_T2D_CC$GO <- "CC"
joined_T2D_BP$GO <- "BP"

combined_T2D_CC_BP <- rbind(joined_T2D_CC,joined_T2D_BP)

ggplot(combined_T2D_CC_BP, aes(x=enrichmentScore.y, y=enrichmentScore.x, label = ifelse(Description.x %in% new_selected_BP_CC, Description.x,""))) + geom_point(aes(size = setSize.x, colour = p.adjust.x)) +
  geom_text_repel(size=2.5, max.overlaps = Inf) + theme_bw() + facet_grid(cols = vars(GO)) +
  scale_color_gradient(low = "#CC6666", high = "Dark grey") +
  scale_y_break(c(-0.25, 0.25), scale=0.5, space = 0.05) + theme(panel.spacing = unit(0.35, "cm")) +
  ylab("GeneRatio (T2D vs Obese)") + xlab("GeneRatio (T2D vs Lean)") +
  scale_y_continuous(breaks = c(-1,-0.8,-0.6,-0.4,0.3,0.5,0.7))

new_selected_BP_CC <- c("immune response","adaptive immune response","defense response",
                        "oxidative phosphorylation","fatty acid metabolic process",
                        "mitochondrial respiratory chain complex assembly","electron transport chain",
                        "collagen-containing extracellular matrix","external encapsulating structure",
                        "extracellular matrix",
                        "mitochondrion","mitochondrial matrix","NADH dehydrogenase complex",
                        "mitochondrial inner membrane","fatty acid elongation")

#Export for Table_S2F
write.table(combined_T2D_CC_BP,"Diabetes_driven_GOterms.txt", append = FALSE, sep = "\t", dec = ".",
            row.names = F, col.names = TRUE)







#### Training GO-enrichment Figure_4D-E
Lean_training_gene_list <- Lean_training_75filt$logFC
names(Lean_training_gene_list) <- Lean_training_75filt$ID
Lean_training_gene_list = sort(Lean_training_gene_list, decreasing = TRUE)

gse <- gseGO(geneList=Lean_training_gene_list,
             ont ="CC",
             keyType = "SYMBOL",
             nPerm = 10000,
             minGSSize = 3,
             maxGSSize = 800,
             pvalueCutoff = 0.1,
             verbose = TRUE,
             OrgDb = org.Hs.eg.db,
             pAdjustMethod = "BH",
             seed = T)

Lean_training_GOCC <- gse@result

gse <- gseGO(geneList=Lean_training_gene_list,
             ont ="BP",
             keyType = "SYMBOL",
             nPerm = 10000,
             minGSSize = 3,
             maxGSSize = 800,
             pvalueCutoff = 0.1,
             verbose = TRUE,
             OrgDb = org.Hs.eg.db,
             pAdjustMethod = "BH",
             seed = T)


Lean_training_GOBP <- gse@result
Lean_training_GOCC$group <- "Lean"

### Obese-training effect
Obese_training_gene_list <- Obese_training_75filt$logFC
names(Obese_training_gene_list) <- Obese_training_75filt$ID
Obese_training_list = sort(Obese_training_gene_list, decreasing = TRUE)

gse <- gseGO(geneList=Obese_training_gene_list,
             ont ="CC",
             keyType = "SYMBOL",
             nPerm = 10000,
             minGSSize = 3,
             maxGSSize = 800,
             pvalueCutoff = 0.1,
             verbose = TRUE,
             OrgDb = org.Hs.eg.db,
             pAdjustMethod = "BH",
             seed = T)
Obese_training_GOCC <- gse@result

gse <- gseGO(geneList=Obese_training_gene_list,
             ont ="BP",
             keyType = "SYMBOL",
             nPerm = 10000,
             minGSSize = 3,
             maxGSSize = 800,
             pvalueCutoff = 0.1,
             verbose = TRUE,
             OrgDb = org.Hs.eg.db,
             pAdjustMethod = "BH",
             seed = T)

Obese_training_GOBP <- gse@result
Obese_training_GOCC$group <- "Obese"

Lean_Obese_training_joined_GOCC <- rbind(Lean_training_GOCC, Obese_training_GOCC)

#Export for Table_S2H
write.table(Lean_Obese_training_joined_GOCC,"Lean_Obese_training_GOCCterms.txt", append = FALSE, sep = "\t", dec = ".",
            row.names = F, col.names = TRUE)


#Figure_4B
Lean_training_names_GOCC <- c("lipoprotein particle","extracellular matrix", "endoplasmic reticulum lumen","supramolecular complex",
                              "immunoglobulin complex","intermediate filament cytoskeleton","HFE-transferrin receptor complex",
                              "mitochondrial membrane")

Obese_training_names_GOCC <- c("lipid droplet","lipoprotein particle","mitochondrion","peroxisome","extracellular matrix",
                               "endoplasmic reticulum lumen","endoplasmic reticulum membrane","cortical cytoskeleton")

comb_OB_LEAN_GOBP_names <- c(Lean_training_names_GOCC,Obese_training_names_GOCC)
Lean_Obese_training_joined_GOCC <- Lean_Obese_training_joined_GOCC[which(Lean_Obese_training_joined_GOCC$Description %in% comb_OB_LEAN_GOBP_names),]

ggplot(Lean_Obese_training_joined_GOCC, aes(x=enrichmentScore, y=Description)) + geom_point(aes(size = setSize, colour = p.adjust)) +
  theme_bw() + facet_grid(cols = vars(group)) +
  geom_vline(xintercept = 0, linetype="solid",
             color = "black", linewidth=0.25) +
  scale_color_gradient(low = "darkblue", high = "gray") +
  theme(panel.spacing = unit(0.35, "cm")) +
  ylab("") + xlab("GeneRatio (Post vs Pre training)")

scale_y_break(c(-0.25, 0.25), scale=0.5, space = 0.05) + theme(panel.spacing = unit(0.35, "cm")) +
  ylab("GeneRatio (T2D vs Obese)") + xlab("GeneRatio (T2D vs Lean)") +
  scale_y_continuous(breaks = c(-1,-0.8,-0.6,-0.4,0.3,0.5,0.7))










####T2D training effect
T2D_training_gene_list <- T2D_training_75filt$logFC
names(T2D_training_gene_list) <- T2D_training_75filt$ID
T2D_training_gene_list = sort(T2D_training_gene_list, decreasing = TRUE)

gse <- gseGO(geneList=T2D_training_gene_list,
             ont ="CC",
             keyType = "SYMBOL",
             nPerm = 10000,
             minGSSize = 3,
             maxGSSize = 800,
             pvalueCutoff = 0.1,
             verbose = TRUE,
             OrgDb = org.Hs.eg.db,
             pAdjustMethod = "BH",
             seed = T)
T2D_training_GOCC <- gse@result
T2D_training_GOCC$GO <- "CC"

gse <- gseGO(geneList=T2D_training_gene_list,
             ont ="BP",
             keyType = "SYMBOL",
             nPerm = 10000,
             minGSSize = 3,
             maxGSSize = 800,
             pvalueCutoff = 0.1,
             verbose = TRUE,
             OrgDb = org.Hs.eg.db,
             pAdjustMethod = "BH",
             seed = T)

T2D_training_GOBP <- gse@result
T2D_training_GOBP$GO <- "BP"

T2D_training_GOCC_GOBP <- rbind(T2D_training_GOCC,T2D_training_GOBP)

#Export for Table_S2I
write.table(T2D_training_GOCC_GOBP,"T2D_training_GOCC_GOBPterms.txt", append = FALSE, sep = "\t", dec = ".",
            row.names = F, col.names = TRUE)

#Figure_4C
T2D_training_names_GOCC <- c("external side of plasma membrane","ribosome","vesicle lumen","secretory granule lumen",
                             "immunoglobulin complex","secondary lysosome","ribonucleoprotein complex",
                             "blood microparticle")

T2D_training_names_GOBP <- c("response to stress","immune response","wound healing","sequestering of iron ion",
                             "RNA splicing","translation","blood coagulation","immunoglobulin mediated immune response")

comb_T2D_GOCC_GOBP_names <- c(T2D_training_names_GOCC,T2D_training_names_GOBP)
T2D_training_GOCC_GOBP <- T2D_training_GOCC_GOBP[which(T2D_training_GOCC_GOBP$Description %in% comb_T2D_GOCC_GOBP_names),]
T2D_training_GOCC_GOBP$GO <- factor(T2D_training_GOCC_GOBP$GO, levels = c("CC","BP"))

ggplot(T2D_training_GOCC_GOBP, aes(x=enrichmentScore, y=Description)) + geom_point(aes(size = setSize, colour = p.adjust)) +
  theme_bw() + facet_grid(cols = vars(GO)) +
  geom_vline(xintercept = 0, linetype="solid",
             color = "black", linewidth=0.25) +
  scale_color_gradient(low = "darkblue", high = "gray") +
  theme(panel.spacing = unit(0.35, "cm")) +
  ylab("") + xlab("GeneRatio (Post vs Pre training)")