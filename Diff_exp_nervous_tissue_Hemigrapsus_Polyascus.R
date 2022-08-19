setwd('')

#Reading data
dat <- data.frame(read.csv("data/Polyascus_proteins.csv"))

#select Area 
dat1_pol <- dat[,c(3, 18:29)]
rownames(dat1_pol) <- dat1_pol[,1]
dat1_pol <- dat1_pol[,-1]
head(dat1_pol)

#####################Information about the samples##################
library(readxl)
fact_pol <- data.frame(read_excel("data/Polyascus_fact.xlsx"))
rownames(fact_pol) <- fact_pol[,1]
fact_pol <- fact_pol[, -1]
head(fact_pol)
fact_pol$Differentiation <- as.factor(fact_pol$Differentiation)
fact_pol$Differentiation

#Infection status as a factor
fact_pol$Status <- as.factor(fact_pol$Status)
fact_pol$Status

#Sex as a factor
fact_pol$Sex <- as.factor(fact_pol$Sex)
fact_pol$Sex

colnames(dat1_pol) <- rownames(fact_pol)

###############Filter##############
#as Differentiation
h_fem <- dat1_pol[which(rowMeans(!is.na(dat1_pol[,rownames(subset(fact_pol,Differentiation=="hf"))])) >= 2/3), ]
h_male <- dat1_pol[which(rowMeans(!is.na(dat1_pol[,rownames(subset(fact_pol,Differentiation=="hm"))])) >= 2/3), ]
in_fem <- dat1_pol[which(rowMeans(!is.na(dat1_pol[,rownames(subset(fact_pol,Differentiation=="if"))])) >= 2/3), ]
in_male <- dat1_pol[which(rowMeans(!is.na(dat1_pol[,rownames(subset(fact_pol,Differentiation=="im"))])) >= 2/3), ]

#as infection status
healthy <- dat1_pol[which(rowMeans(!is.na(dat1_pol[,rownames(subset(fact_pol,Status=="healthy"))])) >= 5/6), ]
infected <- dat1_pol[which(rowMeans(!is.na(dat1_pol[,rownames(subset(fact_pol,Status=="infected"))])) >= 5/6), ]

################Venn diagram#########
library(devtools)
library(ggvenn)


vennn_pol <- list(hth_fem = rownames(h_fem), hth_male = rownames(h_male), inf_fem = rownames(in_fem), inf_male = rownames(in_male))

tiff('Polyascus_Venn_diagram.tiff', bg = 'transparent', width = 1200, height = 900)

ggvenn(
  vennn_pol, 
  fill_color = c("#d200ff", "#04fb04", "#fb0432", "#048afb"),
  stroke_size = 0.5, set_name_size = 8, text_size = 7,
)

dev.off()

venn_status <- list(Healthy = rownames(healthy), Infected = rownames(infected))

tiff('Polyascus_inf_Venn_diagram.tiff', bg = 'transparent', width = 1200, height = 900)

ggvenn(
  venn_status, 
  fill_color = c("#50a4dc", "#ffe255"),
  stroke_size = 0.5, set_name_size = 8, text_size = 7,
)

dev.off()

#Removing NA
colSums(is.na(dat1_pol))
dat2_pol <- dat1_pol[which(rowMeans(!is.na(dat1_pol)) > 0.85), ]
mean(complete.cases(dat2_pol))
colSums(is.na(dat2_pol))

#############Imputation#########
library(impute)
tdat_pol <- t(dat2_pol)
pol_knn1 <- impute.knn(tdat_pol, k = 5)
pol_knn <- t(pol_knn1$data)
head(pol_knn)
mean(complete.cases(pol_knn))


#Expression data distribution
library(RColorBrewer)
pal <- brewer.pal(n = 9, name = "Set1")
cols <- pal[fact_pol$Differentiation]
boxplot(pol_knn, outline = FALSE, col = cols, main = "Raw data")
legend("topright", levels(fact_pol$Differentiation), fill = pal, bty = "n", xpd = T)
colSums(pol_knn)

#######Log-transformation#######
dat_log_pol <- log2(pol_knn+1)
mean(complete.cases(dat_log_pol))
boxplot(dat_log_pol, outline = FALSE, col = cols, main = "Log-transformed data")
legend("topright", levels(fact_pol$Differentiation), fill = pal, bty = "n", xpd = T)

##########Quantile normalization##########
library(limma)
dat_norm_pol <- normalizeQuantiles(dat_log_pol)
boxplot(dat_norm_pol, col = cols, main = "Normalized data")
legend("topright", levels(fact_pol$Differentiation), fill = pal, bty = "n", xpd = T)
mean(complete.cases(dat_norm_pol))
colSums(is.na(dat_norm_pol))

###########nMDS##########
library('vegan')

tdat_norm_pol <- t(dat_norm_pol)
pol_ord <- metaMDS(tdat_norm_pol,
                   distance = "euclidean",
                   autotransform = FALSE)

pol_ord$stress

#######nMDS ggplot######
nmds_scrs <- as.data.frame(scores(pol_ord, display = "sites"))
nmds_scrs <- cbind(nmds_scrs, Differentiation = fact_pol$Differentiation)
pal_colour <- c("#d200ff", "#04fb04", "#fb0432", "#048afb")


nmds_plot <- ggplot(nmds_scrs, aes(x = NMDS1, y = NMDS2)) +
  scale_colour_manual(values = pal_colour) +
  geom_polygon(aes(group = factor(Differentiation), colour = factor(Differentiation)), size = 1, fill = 'transparent') +
  geom_point(aes(x = NMDS1, y = NMDS2, colour = factor(Differentiation)), size = 4) +
  coord_fixed() +
  theme_classic() + 
  labs(colour = 'Differentiation') +
  theme(legend.position = "right", legend.text = element_text(size = 12), legend.title = element_text(size = 12), axis.text = element_text(size = 12)) 

tiff('nMDS.tiff', units="in", width=12, height=8, res=600, compression = 'lzw')
nmds_plot
dev.off()

#######PLS-DA#############
library(mixOmics)
plsda_pol <- plsda(tdat_norm_pol, fact_pol$Differentiation, ncomp = 11)
symb <- c(16)

tiff('plsda.tiff', units="in", width=12, height=8, res=600, compression = 'lzw')
plotIndiv(plsda_pol, legend = TRUE, ellipse = TRUE, ind.names = FALSE,
          pch = symb, col = pal_colour, title = 'PLS-DA')
dev.off()

#####################sPLS-DA################
#Model selection
list.keepX <- c(5:10,  seq(20, 100, 10))

opt_val <- tune.splsda(tdat_norm_pol, fact_pol$Differentiation, ncomp = 3, nrepeat = 1000, folds = 3, test.keepX = list.keepX, cpus = 6)
opt_val$choice.keepX
opt_val$choice.ncomp

pal_colour <- c("#d200ff", "#04fb04", "#fb0432", "#048afb")
symb <- c(16)
splsda_pol <- splsda(tdat_norm_pol, fact_pol$Differentiation, ncomp = 2, keepX = c(40, 20))


tiff('sPLS-DA.tiff', units="in", width=11, height=8, res=600, compression = 'lzw')
plotIndiv(splsda_pol, legend = TRUE, ellipse = TRUE, ind.names = FALSE,
          pch = symb, col = pal_colour, title = 'sPLS-DA')
dev.off()

#####Loadings####
library('dplyr')

plsda_load <- splsda_pol$loadings
plsda_mat <- plsda_load$X
plsda_mat <- as.data.frame(plsda_mat)

comp1_plot <- plotLoadings(splsda_pol, comp = 1, method = 'mean', contrib = 'max', ndisplay = 20, 
             title = 'Loading 1 component', size.title = 2, legend = TRUE)  


comp2_plot <- plotLoadings(splsda_pol, comp = 2, method = 'mean', contrib = 'max', ndisplay = 20,
             title = 'Loadings 2 component', size.title = 2, legend = FALSE)  

dev.off()


#####Component 1####
comp1_loadings <- plsda_load$X[,1]
comp1_loadings <- as.data.frame(comp1_loadings)
max_comp1 <- comp1_loadings %>% 
  arrange(desc(abs(comp1_loadings))) %>%
  slice(1:20)
max_comp1$Accession <- row.names(max_comp1)

group_contrib <- comp1_plot$GroupContrib 
max1_contr <- cbind(max_comp1, group_contrib)

descr <- dat[,c(3,46)]

load_filtered <- descr %>%
  filter(descr$Accession %in% max1_contr$Accession)

max_comp1_descr <- left_join(max1_contr, load_filtered, by = 'Accession')

max_comp1_descr <- max_comp1_descr %>%
  rename(Loadings = comp1_loadings)
max_comp1_descr <- max_comp1_descr %>%
  rename(Contributions = group_contrib)

max1_descr <- max_comp1_descr[c(1,2,4:9,11:15,20),]
row.names(max1_descr) <- max1_descr$Description
max1_descr <- max1_descr[c(1:10),c(1,3,4)]

row.names(max1_descr) <- gsub("OS.*.", replacement = "", row.names(max1_descr))
max1_descr$Description <- gsub("OS.*.", replacement = "", max1_descr$Description)

pal_colour <- c("#d200ff", "#04fb04", "#fb0432", "#048afb")
bar_colour <- c("#d200ff", "#04fb04", "#048afb")

comp1_bar <- ggplot(max1_descr, aes(x = reorder(Description, abs(Loadings)), y = Loadings, fill = Contributions)) +
  geom_bar(stat = 'identity', colour = 'black') + xlab('Proteins') + coord_flip() +
  scale_fill_manual(values = bar_colour) + 
  theme_classic() +
  theme(axis.text = element_text(size = 25), axis.title = element_text(size = 30),
        legend.title = element_text(size = 30), legend.text = element_text(size = 25),
        legend.spacing.y = unit(1, 'cm')) +
  guides(fill = guide_legend(byrow = TRUE))
  
tiff('Loadings1.tiff', units="px", width=13200, height=5500, res=600, compression = 'none')
comp1_bar
dev.off()

#####Component2####
comp2_loadings <- plsda_load$X[,2]
comp2_loadings <- as.data.frame(comp2_loadings)
max_comp2 <- comp2_loadings %>% 
  arrange(desc(abs(comp2_loadings))) %>%
  slice(1:20)
max_comp2$Accession <- row.names(max_comp2)

group_contrib2 <- comp2_plot$GroupContrib 
max2_contr <- cbind(max_comp2, group_contrib2)

descr <- dat[,c(3,46)]

load_filtered2 <- descr %>%
  filter(descr$Accession %in% max2_contr$Accession)

max_comp2_descr <- left_join(max2_contr, load_filtered2, by = 'Accession')

max_comp2_descr <- max_comp2_descr %>%
  rename(Loadings = comp2_loadings)
max_comp2_descr <- max_comp2_descr %>%
  rename(Contributions = group_contrib2)


max2_descr <- max_comp2_descr[c(1,4:8,10,14,17:20),]
row.names(max2_descr) <- max2_descr$Description
max2_descr <- max2_descr[c(1:10),c(1,3,4)]

max2_descr$Description <- gsub("OS.*.", replacement = "", max2_descr$Description)
pal_colour <- c("#d200ff", "#04fb04", "#fb0432", "#048afb")

comp2_bar <- ggplot(max2_descr, aes(x = reorder(makeUnique(Description), abs(Loadings)), y = Loadings, fill = Contributions)) +
  geom_bar(stat = 'identity', colour = 'black') + xlab('Proteins') + coord_flip() +
  scale_fill_manual(values = pal_colour) + 
  theme_classic() +
  theme(axis.text = element_text(size = 25), axis.title = element_text(size = 30),
        legend.title = element_text(size = 30), legend.text = element_text(size = 25),
        legend.spacing.y = unit(1, 'cm')) +
  guides(fill = guide_legend(byrow = TRUE))

tiff('Loadings2.tiff', units="px", width=13200, height=5500, res=600, compression = 'none')
comp2_bar
dev.off()


###############Filtered proteins############## 
sign_prot <- rownames(filter(plsda_mat, comp1 != 0 | comp2 != 0))
head(sign_prot)

#########Differential expression analysis############
X_pol <- model.matrix(~ 0 + fact_pol$Differentiation)
X_pol
colnames(X_pol) <- c('hlth_fem', 'hlth_male', 'inf_fem', 'inf_male')

dat_norm_pol_prot <- as.data.frame(dat_norm_pol)

names <- sign_prot

dat_filtered <- dat_norm_pol_prot %>%
  filter(row.names(dat_norm_pol_prot) %in% names)

dat_filtered <- as.matrix(dat_filtered)

fit_filter <- lmFit(dat_filtered, design = X_pol, method = "robust", maxit = 10000)

contrast.matrix <- makeContrasts(hlth_fem-hlth_male, hlth_fem-inf_fem, hlth_fem-inf_male, hlth_male-inf_male, hlth_male-inf_fem, inf_fem-inf_male, levels = X_pol)

fit2_filter <- contrasts.fit(fit_filter, contrast.matrix)

# Empirical Bayes statistics
efit_filter <- eBayes(fit2_filter)

#Differential expression table
topTable(efit_filter)
numGenes_filter <- length(dat_filtered)
full_list_efit_filter <- topTable(efit_filter, number = numGenes_filter)
write.csv(full_list_efit_filter,'Dif_expr_Polyascus_sPLS-DA.csv')
head(full_list_efit_filter)

#Data filtered
library(dplyr)
exp_fil <- data.frame(read.csv("Dif_expr_Polyascus_sPLS-DA.csv"))
exp_fil <- exp_fil %>%
  rename(Accession = X)

descr <- dat[,c(3,46)]

descr_filtered <- descr %>%
  filter(descr$Accession %in% names)

descr_fil_pol <- left_join(descr_filtered, exp_fil, by ='Accession')
descr_fil_pol$Protein <- as.factor(gsub("OS.*.", replacement = "", descr_fil_pol$Description))
write.csv(descr_fil_pol,'Dif_expr_Polyascus_sPLS-Da_description.csv')

#Differential expressed proteins
p_above_fil <- exp_fil$adj.P.Val <= 0.05 
sum(p_above_fil)
p_protein_fil <- exp_fil[c(1:44), ]
accessions_fil <- p_protein_fil$Accession 
write.csv(accessions_fil,'Dif_expr_sPLS-DA_Polyascus_44_prot.csv')
accessions2_fil <- p_protein_fil[, c(1, 11)]
accessions2_fil
acc <- accessions2_fil$Accession

pv_filt <- descr %>%
  filter(descr$Accession %in% accessions2_fil$Accession)

descr_pv_pol <- left_join(pv_filt, accessions2_fil, by ='Accession')

descr_padj_pol <- descr_pv_pol %>% 
  arrange(adj.P.Val)

write.csv(descr_padj_pol,'Dif_expr_descr_sPLS-DA_Polyascus_44.csv')

######Tables with differential expression####
f_dif_all <- descr_fil_pol$adj.P.Val <= 0.05
f_dif_all2 <- descr_fil_pol[f_dif_all,]

for_obt_heatmap <- rownames(dat_filtered) %in% f_dif_all2$Accession
for_obt_heatmap_filt <- dat_filtered[for_obt_heatmap, ]
for_obt_heatmap_filt1 <- as.data.frame(for_obt_heatmap_filt)
for_obt_heatmap_filt1$Accession <- rownames(for_obt_heatmap_filt1)
heatmap_filt <- left_join(for_obt_heatmap_filt1, pv_filt, by ='Accession')

heatmap_descr <- heatmap_filt[!duplicated(heatmap_filt$Description), ]
row.names(heatmap_descr) <- heatmap_descr$Description
hm_descr <- heatmap_descr[c(1:4, 10:12, 14:22, 24:25, 27:28, 30:32, 34:37), 1:12]
row.names(hm_descr) <- gsub("OS.*.", replacement = "", row.names(hm_descr))

m_hm_descr <- as.matrix(hm_descr)

#Table without repeats
table_filt <- left_join(for_obt_heatmap_filt1, descr_padj_pol, by ='Accession')
table_descr <- table_filt[!duplicated(table_filt$Description), ]
row.names(table_descr) <- table_descr$Description

table_prot_descr <- table_descr[c(1:4, 10:12, 14:22, 24:25, 27:28, 30:32, 34:37), 14:15]

write.csv(table_prot_descr,'Dif_prot_final_sPLS-DA_descr_Pp_27.csv')


######Heatmap with accession####
library('ComplexHeatmap')
library('dendextend')

col_dend_acc <- hclust(dist(t(for_obt_heatmap_filt), method = 'euclidean'), method = 'average')
row_dend <- hclust(dist(for_obt_heatmap_filt, method = 'euclidean'), method = 'average')
list_col <- list(Groups = c('hf' = "#d200ff", 'hm' = "#04fb04", 'if' = "#fb0432", 'im' = "#048afb"))
heat_annot <- HeatmapAnnotation(Groups = fact_pol$Differentiation, col = list_col, 
                                annotation_name_gp= gpar(fontsize = 15),
                                annotation_legend_param = list(Groups = list(
                                  title = 'Groups', title_gp = gpar(fontsize = 15),
                                  labels_gp = gpar(fontsize = 12),
                                  at = c("hf", "hm", 'if', 'im'), 
                                  labels = c("Healthy female", "Healthy male",
                                             'Infected female', 'Infected male'))))

heat_acc <- Heatmap(for_obt_heatmap_filt, 
        name = "Differential protein expression", 
        column_title = "Variables", row_title = "Samples",
        row_names_gp = gpar(fontsize = 7),
        cluster_rows = row_dend,
        cluster_columns = color_branches(col_dend_acc, k = 3),
        top_annotation = heat_annot)

######Heatmap with description#####
col_dend_descr <- hclust(dist(t(hm_descr), method = 'euclidean'), method = 'average')

heat_descr <- Heatmap(m_hm_descr,
        name = "Expression", 
        column_title = "Variables", row_title = "Samples",
        rect_gp = gpar(col = "white", lwd = 2),
        row_names_gp = gpar(fontsize = 15),
        row_names_max_width = max_text_width(
          rownames(m_hm_descr), 
          gp = gpar(fontsize = 20)),
        cluster_columns = color_branches(col_dend_descr, k = 2),
        top_annotation = heat_annot,
        heatmap_legend_param = list(labels_gp = gpar(fontsize = 12), 
                                    title_gp = gpar(fontsize = 15)))

dev.off()


tiff('heatmap_descr_gr.tiff', units="px", width=7000, height=4800, res=600, compression = 'none')
draw(heat_descr, heatmap_legend_side = "left", annotation_legend_side = "bottom")


dev.off()



#######Boxplots######
#Hemocyanin subunit 6 OS=Eriocheir sinensis 0.000554452
descr_fil_pol[descr_fil_pol$Accession == "K4EJG5|K4EJG5_ERISI", ]
boxplot(dat_filtered[c("K4EJG5|K4EJG5_ERISI"),] ~ Differentiation, data = fact_pol,
        varwidth = TRUE, log = "y", las = 1)

#FTCD_N domain-containing protein OS=Scylla olivacea 0.000226069
descr_fil_pol[descr_fil_pol$Accession == "A0A0P4WFR3|A0A0P4WFR3_SCYOL", ]
boxplot(dat_filtered[c("A0A0P4WFR3|A0A0P4WFR3_SCYOL"),] ~ Differentiation, data = fact_pol,
        varwidth = TRUE, log = "y", las = 1)

# Proteasome subunit alpha type OS=Scylla olivacea 0.004762967
accessions2_fil[3,]
descr_fil_pol[descr_fil_pol$Accession == "A0A0N7ZDN1|A0A0N7ZDN1_SCYOL", ]
boxplot(dat_filtered[c("A0A0N7ZDN1|A0A0N7ZDN1_SCYOL"),] ~ Differentiation, data = fact_pol,
        varwidth = TRUE, log = "y", las = 1)

#Glucosamine-6-phosphate isomerase 0.005143184
accessions2_fil[5,]
descr_fil_pol[descr_fil_pol$Accession == "A0A5B7DES9|A0A5B7DES9_PORTR", ]
boxplot(dat_filtered[c("A0A5B7DES9|A0A5B7DES9_PORTR"),] ~ Differentiation, data = fact_pol,
        varwidth = TRUE, log = "y", las = 1)

View(descr_padj_pol)

#Glycogen debrancher OS=Scylla olivacea (?)
boxplot(dat_filtered[c("A0A0P4W3T4|A0A0P4W3T4_SCYOL"),] ~ Differentiation, data = fact_pol,
        varwidth = TRUE, log = "y", las = 1)

#Hemocyanin subunit 2
boxplot(dat_filtered[c("sp|C0HLU8|HCY_SCYSE"),] ~ Differentiation, data = fact_pol,
        varwidth = TRUE, log = "y", las = 1)

#Argininosuccinate lyase OS=Penaeus vannamei
boxplot(dat_filtered[c("A0A3R7MNM6|A0A3R7MNM6_PENVA"),] ~ Differentiation, data = fact_pol,
        varwidth = TRUE, log = "y", las = 1)

#UTP--glucose-1-phosphate uridylyltransferase
boxplot(dat_filtered[c("A0A0P4WF84|A0A0P4WF84_SCYOL"),] ~ Differentiation, data = fact_pol,
        varwidth = TRUE, log = "y", las = 1)

#Histone H4 OS=Tigriopus californicus
boxplot(dat_filtered[c("A0A553NPY7|A0A553NPY7_TIGCA"),] ~ Differentiation, data = fact_pol,
        varwidth = TRUE, log = "y", las = 1)

#Calbindin-32-like
boxplot(dat_filtered[c("A0A6A7G3C7|A0A6A7G3C7_9CRUS"),] ~ Differentiation, data = fact_pol,
        varwidth = TRUE, log = "y", las = 1)

#Alpha-2-macroglobulin (Fragment)
boxplot(dat_filtered[c("A0A068J627|A0A068J627_ERISI"),] ~ Differentiation, data = fact_pol,
        varwidth = TRUE, log = "y", las = 1)

#Fructose-bisphosphate aldolase
boxplot(dat_filtered[c("A0A0P4W855|A0A0P4W855_SCYOL"),] ~ Differentiation, data = fact_pol,
        varwidth = TRUE, log = "y", las = 1)

#Paxillin 
boxplot(dat_filtered[c("A0A0P4W3C8|A0A0P4W3C8_SCYOL"),] ~ Differentiation, data = fact_pol,
        varwidth = TRUE, log = "y", las = 1)

#Actin-related protein = supressor of profilin
boxplot(dat_filtered[c("A0A0P4WLJ6|A0A0P4WLJ6_SCYOL"),] ~ Differentiation, data = fact_pol,
        varwidth = TRUE, log = "y", las = 1)

#	Transcriptional activator protein Pur-alpha
boxplot(dat_filtered[c("A0A5B7E8T1|A0A5B7E8T1_PORTR"),] ~ Differentiation, data = fact_pol,
        varwidth = TRUE, log = "y", las = 1)

#Ferritin
boxplot(dat_filtered[c("D6MXQ4|D6MXQ4_ERISI"),] ~ Differentiation, data = fact_pol,
        varwidth = TRUE, log = "y", las = 1)

#Thioredoxin
boxplot(dat_filtered[c("C4N5V0|C4N5V0_ERISI"),] ~ Differentiation, data = fact_pol,
        varwidth = TRUE, log = "y", las = 1)

#Glutathione transferase 
boxplot(dat_filtered[c("A0A4D7AVC7|A0A4D7AVC7_ERISI"),] ~ Differentiation, data = fact_pol,
        varwidth = TRUE, log = "y", las = 1)

#Prostaglandin D synthase OS=Eriocheir sinensis
boxplot(dat_filtered[c("I6LWU2|I6LWU2_ERISI"),] ~ Differentiation, data = fact_pol,
        varwidth = TRUE, log = "y", las = 1)

#14-3-3 protein 
boxplot(dat_filtered[c("A0A385L4G6|A0A385L4G6_PROCL"),] ~ Differentiation, data = fact_pol,
        varwidth = TRUE, log = "y", las = 1)

############sPLS-DA Volcano plot with proteins name############
tiff('fin_Vulcano_des_Polyascus_hf_if.tiff', units="px", width=6600, height=4800, res=600, compression = 'none')

EnhancedVolcano(descr_fil_pol,
                lab = descr_fil_pol$Protein,
                x = 'hlth_fem...inf_fem',
                y = 'adj.P.Val', # ??? ????? adj. p.val? ? ?????? ?? ?????? p. val ?? ????????????
                pCutoff = 0.05,
                xlim = c(-3, 3),
                ylim = c(0, 3),
                FCcutoff = 1, # FCCutoff ????? ??????? ?????, ?? ???? 1 (???????? ? ??? ???? ? ?????). ??? ????, ??? ????? ???????? ???????? ?????. 
                title ="Healthy vs infected  female crabs",
                labSize = 6,
                boxedLabels = F,
                colAlpha = 1)

dev.off()





############Differential expression of all proteins#############
X_pol <- model.matrix(~ 0 + fact_pol$Differentiation)
X_pol
colnames(X_pol) <- c('hlth_fem', 'hlth_male', 'inf_fem', 'inf_male')

fit_pol <- lmFit(dat_norm_pol, design = X_pol, method = "robust", maxit = 10000)

contrast.matrix <- makeContrasts(hlth_fem-hlth_male, hlth_fem-inf_fem, hlth_male-inf_male, hlth_male-inf_fem, inf_fem-inf_male, inf_male-hlth_fem, levels = X_pol)

fit2_pol <- contrasts.fit(fit_pol, contrast.matrix)

# Empirical Bayes statistics
efit_pol <- eBayes(fit2_pol)

#Differential expression table
topTable(efit_pol)
numGenes_pol <- length(dat_norm_pol)
full_list_efit_pol <- topTable(efit_pol, number = numGenes_pol)

write.csv(full_list_efit_pol,'Dif_expr_Polyascus_all_proteins.csv')
head(full_list_efit_pol)



##########Volcano-plots##############
library(EnhancedVolcano)

head(full_list_efit_pol)
tiff('Vulcano_Polyascus_im_hf.tiff', units="in", width=11, height=8, res=600, compression = 'lzw')

EnhancedVolcano(full_list_efit_pol,
                lab = rownames(full_list_efit_pol),
                x = 'hlth_fem...inf_fem',
                y = 'adj.P.Val', # ??? ????? adj. p.val? ? ?????? ?? ?????? p. val ?? ????????????
                pCutoff = 0.05,
                #xlim = c(-8, 10),
                #ylim = c(0, 5),
                FCcutoff = 1, # FCCutoff ????? ??????? ?????, ?? ???? 1 (???????? ? ??? ???? ? ?????). ??? ????, ??? ????? ???????? ???????? ?????. 
                title ="Vulcano",
                labSize = 4.0,
                boxedLabels = F,
                colAlpha = 1)

dev.off()

############Volcano plots with proteins names############
experiment <- data.frame(read.csv("Vul_Dif_expr_Polyascus.csv"))
head(experiment)
descr <- dat[,c(3,46)]
head(descr)
descr_full_pol <- left_join(descr, experiment, by ='Accession')
descr_full_pol$Protein <- as.factor(gsub("OS.*.", replacement = "", descr_full_pol$Description))
write.csv(descr_full_pol,'Dif_expr_Polyascus_description.csv')
head(descr_full_pol)

tiff('Vulcano_Polyascus_hm_im.tiff', units="px", width=6600, height=4800, res=600, compression = 'none')

EnhancedVolcano(descr_full_pol,
                lab = descr_full_pol$Protein,
                x = 'hlth_fem...inf_fem',
                y = 'adj.P.Val', # ??? ????? adj. p.val? ? ?????? ?? ?????? p. val ?? ????????????
                pCutoff = 0.05,
                xlim = c(-3, 3.5),
                ylim = c(0, 3),
                FCcutoff = 1, # FCCutoff ????? ??????? ?????, ?? ???? 1 (???????? ? ??? ???? ? ?????). ??? ????, ??? ????? ???????? ???????? ?????. 
                title ="Healthy vs infected female crabs",
                labSize = 4.7,
                boxedLabels = F,
                colAlpha = 1)

dev.off()

#Differential expressed proteins
p_above_pol <- experiment$adj.P.Val <= 0.05
sum(p_above_pol)
p_protein <- experiment[c(1:26), ]
accessions2 <- p_protein$Accession 
write.csv(accessions2,'Dif_expr_26_prot.csv')
accessions2 <- p_protein[, c(1, 11)]
head(accessions2)

names <- data.frame(read.csv("Dif_expr_26_prot_with_desc.csv"))
names <- names[, -1]
View(names)

######Heatmaps#############
library(gplots)

p_above_pol2 <- experiment[p_above_pol,]
head(p_above_pol2)
head(dat_norm_pol)

for_obt_hm_all <- rownames(dat_norm_pol) %in% p_above_pol2$Accession
for_obt_hm_all_f <- dat_norm_pol[for_obt_hm_all, ]
head(for_obt_hm_all_f)

for_obt_hm_all_f1 <- as.data.frame(for_obt_hm_all_f)
for_obt_hm_all_f1$Accession <- rownames(for_obt_hm_all_f1)

descr_for_all <- names[, -1]
hm_all_filt <- left_join(for_obt_hm_all_f1, descr_for_all, by ='Accession')

View(hm_all_filt)

hm_all_descr <- hm_all_filt[!duplicated(hm_all_filt$Description), ]
row.names(hm_all_descr) <- hm_all_descr$Description
View(hm_all_descr)
hm_descr2 <- hm_all_descr[c(1:6, 10:13, 15:16, 22:23), 1:12]
View(hm_descr2)
row.names(hm_descr2) <- gsub("OS.*.", replacement = "", row.names(hm_descr2))
head(hm_descr)

prot_table <- hm_all_descr[c(1:6, 10:13, 15:16, 22:23), 14:15]
View(prot_table)
write.csv(prot_table,'Dif_prot_padj_Pp.csv')


#С Accession
heatmap_obt <- colorpanel(75, low = '#00d2ff', high = '#ff004e')
heatmap.2(x = as.matrix(for_obt_hm_all_f), col = heatmap_obt, scale = 'none',
          distfun = function(x) dist(x, method = 'euclidean'),
          hclustfun = function(x) hclust(x, method = 'average'),
          key = TRUE, symkey = FALSE, density.info = 'none',
          trace = 'none', cexRow = 1, cexCol = 1, margins = c(6,9),
          keysize = 1)

#with description
heatmap.2(x = as.matrix(hm_descr2), col = heatmap_obt, scale = 'none',
          distfun = function(x) dist(x, method = 'euclidean'),
          hclustfun = function(x) hclust(x, method = 'average'),
          key = TRUE, symkey = FALSE, density.info = 'none',
          trace = 'none', cexRow = 1, cexCol = 1, margins = c(5,20),
          keysize = 0.7)
dev.off()

###############Boxplots for differential expressed proteins############

boxplot(dat_norm_pol[c("A0A2P2I8Q8|A0A2P2I8Q8_9CRUS"),] ~ Differentiation, data = fact_pol,
        varwidth = TRUE, log = "y", las = 1)

tiff('Hemocyanin_boxplot_norm.tiff', units="in", width=11, height=8, res=600, compression = 'lzw')
boxplot(dat_norm_pol[c("K4EJG5|K4EJG5_ERISI"),] ~ Differentiation, data = fact_pol,
        varwidth = TRUE, log = "y", las = 1)
dev.off()

tiff('FTCD_N domain-containing protein_boxplot_norm.tiff', units="in", width=11, height=8, res=600, compression = 'lzw')
boxplot(dat_norm_pol[c("A0A0P4WFR3|A0A0P4WFR3_SCYOL"),] ~ Differentiation, data = fact_pol,
        varwidth = TRUE, log = "y", las = 1)
dev.off()

tiff('Proteasome subunit alpha type OS=Scylla olivacea_boxplot_norm.tiff', units="in", width=11, height=8, res=600, compression = 'lzw')
boxplot(dat_norm_pol[c("A0A0N7ZDN1|A0A0N7ZDN1_SCYOL"),] ~ Differentiation, data = fact_pol,
        varwidth = TRUE, log = "y", las = 1)
dev.off()

tiff('Glucosamine-6-phosphate isomerase_boxplot_norm.tiff', units="in", width=11, height=8, res=600, compression = 'lzw')
boxplot(dat_norm_pol[c("A0A5B7DES9|A0A5B7DES9_PORTR"),] ~ Differentiation, data = fact_pol,
        varwidth = TRUE, log = "y", las = 1)
dev.off()


tiff('Glycogen debrancher OS=Scylla olivacea_boxplot_norm.tiff', units="in", width=11, height=8, res=600, compression = 'lzw') 
boxplot(dat_norm_pol[c("A0A0P4W3T4|A0A0P4W3T4_SCYOL"),] ~ Differentiation, data = fact_pol,
        varwidth = TRUE, log = "y", las = 1)
dev.off()

tiff('Hemocyanin subunit 2 OS=Carcinus aestuarii _boxplot_norm.tiff', units="in", width=11, height=8, res=600, compression = 'lzw') 
boxplot(dat_norm_pol[c("sp|P84293|HCY2_CARAE"),] ~ Differentiation, data = fact_pol,
        varwidth = TRUE, log = "y", las = 1)
dev.off()


tiff('Argininosuccinate lyase OS=Penaeus vannamei_boxplot_norm.tiff', units="in", width=11, height=8, res=600, compression = 'lzw') 
boxplot(dat_norm_pol[c("A0A3R7MNM6|A0A3R7MNM6_PENVA"),] ~ Differentiation, data = fact_pol,
        varwidth = TRUE, log = "y", las = 1)
dev.off()

tiff('UTP--glucose-1-phosphate uridylyltransferase _boxplot_norm.tiff', units="in", width=11, height=8, res=600, compression = 'lzw')
boxplot(dat_norm_pol[c("A0A0P4WF84|A0A0P4WF84_SCYOL"),] ~ Differentiation, data = fact_pol,
        varwidth = TRUE, log = "y", las = 1)
dev.off()

tiff('Uncharacterized protein OS=Scylla olivacea_boxplot_norm.tiff', units="in", width=11, height=8, res=600, compression = 'lzw')
boxplot(dat_norm_pol[c("A0A0P4VS89|A0A0P4VS89_SCYOL"),] ~ Differentiation, data = fact_pol,
        varwidth = TRUE, log = "y", las = 1)
dev.off()

tiff('Histone H4 OS=Tigriopus californicus_boxplot_norm.tiff', units="in", width=11, height=8, res=600, compression = 'lzw')
boxplot(dat_norm_pol[c("A0A553NPY7|A0A553NPY7_TIGCA"),] ~ Differentiation, data = fact_pol,
        varwidth = TRUE, log = "y", las = 1)
dev.off()


