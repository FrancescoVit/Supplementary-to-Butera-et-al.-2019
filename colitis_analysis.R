########################################################################################
########################################################################################

# This is the R script used for microbiota data analysis in the manuscript

# "IL-13 mRNA tissue content identifies two subsets of adult ulcerative colitis patients 
# with different clinical and mucosa-associated microbiota profiles"

# By authors: Alessia Butera, Monica Di Paola, Francesco Vitali, Daniela De Nitto, 
# Francesco Covotta, Francesco Borrini, Roberta Pica, Carlotta De Filippo, Duccio Cavalieri, 
# Alessandro Giuliani, Annamaria Pronio, Monica Boirivan

# Script author: Francesco Vitali

# DISCLAIMER: this is the collection of bash command that were used for analysis. 
# However, for publication purpose, absolute filepath of my local PC were changed with 
# variables (declared at the beginning of the script). As such, they were not tested and 
# are not guaranteed to work

# DISCLAIMER: Use these command and the information contained here at your own risk!
# I'm not responsible for loss of data

########################################################################################
########################################################################################


#########################################################################
## LOAD LIBRARY #########################################################

library(phyloseq)
library(ggplot2)
library(scales)
library(ggsci)
library(ggpubr)
library(microbiome)
library(cowplot)
library(doBy)
library(vegan)
library(plyr)


################################################################################
## DATA LOADING, CLEANING, FIRST VISUALIZATIONS ################################

# Read data in R 
biom_16S_fp <- ## insert here the filepath to the .biom table obtained with QIIME
map_16S_fp <- ## insert here the filepath to the mapping file
tree_16S_fp <- ## insert here the filepath to the phylogenetic tree file obtained with QIIME

# Importing otu table 
otu_biom_16S<-  import_biom(BIOMfilename=biom_16S_fp)
# Importing mapping file
env_16S <- import_qiime_sample_data(map_16S_fp)
# Importing tree file
tree_16S <- read_tree_greengenes(treefile = tree_16S_fp)
# Building phyloseq object
qiime_file_16S_citochine <- merge_phyloseq(otu_biom_16S, env_16S, tree_16S)

# Convertin variables in sample_data()
dfs <- as(sample_data(qiime_file_16S_citochine), "data.frame") 
dfs$IL17.qPCR. <- as.numeric(as.character(dfs$IL17.qPCR.))
dfs$IL13.qPCR. <- as.numeric(as.character(dfs$IL13.qPCR.))
sample_data(qiime_file_16S_citochine) <- dfs

# Renaming phylogenetic rank 
colnames(tax_table(qiime_file_16S_citochine)) <- c(Rank1 = "Kingdom", Rank2 = "Phylum", Rank3 = "Class", Rank4 = "Order", Rank5 = "Family", Rank6 = "Genus", Rank7 = "Species")

# Filtering otu mitochondria and chloroplast OTUs

qiime_file_16S_citochine <- subset_taxa(qiime_file_16S_citochine, Family  != "f__Mitochondria")
qiime_file_16S_citochine <- subset_taxa(qiime_file_16S_citochine, Class  != "c__Chloroplast")
qiime_file_16S_citochine <- subset_taxa(qiime_file_16S_citochine, Kingdom  == "k__Bacteria")

# Visualize distribution of reads in samples

sdt = data.table(as(sample_data(qiime_no_remission), "data.frame"),
                 TotalReads = sample_sums(qiime_no_remission), keep.rownames = TRUE)
setnames(sdt, "rn", "SampleID")
pSeqDepth = ggplot(sdt, aes(TotalReads)) + geom_histogram(col="Black", fill="gray") 
pSeqDepth + theme_bw() + xlab("Numero di sequenze") + ylab("Frequenza") + 
  scale_x_continuous(breaks=pretty_breaks(19), labels = function(x) format(x, scientific = TRUE)) +
  scale_y_continuous(breaks=pretty_breaks(20)) +  theme(axis.text.x = element_text(size=8, angle=90))

# total number of reads
sums <- sample_sums(qiime_no_remission)
sum(sums)

# mean number of reads per category "Pathology"
means <- aggregate(TotalReads ~  Pathology, sdt, mean)
# plot per category
sdt = data.table(as(sample_data(qiime_no_remission), "data.frame"),
                 TotalReads = sample_sums(qiime_no_remission), keep.rownames = TRUE)
setnames(sdt, "rn", "SampleID")
pSeqDepth = ggplot(sdt, aes(y=TotalReads, x=Pathology)) + geom_boxplot(col="Black", fill="gray") 
pSeqDepth + theme_bw() + ylab("Numero di sequenze") + xlab("Tessuto") + 
  scale_y_continuous(breaks=pretty_breaks(15), labels = function(x) format(x, scientific = TRUE)) +  
  theme(axis.text.x = element_text(size=10)) + 
  stat_summary(fun.y=mean, colour="darkred", geom="point", shape=15, size=3,show_guide = FALSE) +
  geom_text(data = means, aes(label = round(TotalReads, digits = 0), y = TotalReads + 5000))

# mean number of reads per category "Tissue"
means <- aggregate(TotalReads ~  Tissue, sdt, mean)
# plot per category
sdt = data.table(as(sample_data(qiime_no_remission), "data.frame"),
                 TotalReads = sample_sums(qiime_no_remission), keep.rownames = TRUE)
p_rich <- ggboxplot(data = sdt, x = "Tissue", y = "TotalReads", add = "jitter", fill = "Tissue", palette ="npg")
p_rich <- p_rich + stat_compare_means(label.y = 350000, method = "anova")
my_comparisons <- list(c("CM","HM"),c("CM","LM"),c("HM","LM"))
p_rich <- p_rich + stat_compare_means(comparisons = my_comparisons, method = "t.test")
p_rich

# Summary of samples data
summary(sample_data(qiime_no_remission)$Tissue)
summary(sample_data(qiime_no_remission)$Pathology)
summary(sample_data(qiime_no_remission)$sesso) # means "Gender" in italian


################################################################################
## DATA TRANSFORMATION #########################################################

## CSS + log transformation, see ref biblio https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4010126/

data.metagenomeSeq = phyloseq_to_metagenomeSeq(qiime_no_remission)
p = cumNormStat(data.metagenomeSeq) #default is 0.5
data.cumnorm = cumNorm(data.metagenomeSeq, p=p)
#data.cumnorm
data.CSS = MRcounts(data.cumnorm, norm=TRUE, log=TRUE)
dim(data.CSS)  # make sure the data are in a correct formal: number of samples in rows
qiime_file_trasf_css <- qiime_no_remission
otu_table(qiime_file_trasf_css) <- otu_table(data.CSS, taxa_are_rows = T)


##########################################################################
## DATA ANALYSIS #########################################################

#### ALPHA DIVERSITY 

rich <- global(qiime_file_trasf_css, index = "Richness")
shan <- global(qiime_file_trasf_css, index = "Shannon")
eve <- global(qiime_file_trasf_css, index = "Pielou")

df <- as(sample_data(qiime_file_trasf_css), "data.frame")

df$rich <- rich$richness_0
df$shan <- shan$diversities_shannon
df$eve <- eve$evenness_pielou

head(df)

# get summary
summaryBy(rich ~ Tissue, data = df,
          FUN = function(x) { c(m = mean(x), s = sd(x)) } )
summaryBy(shan ~ Tissue, data = df,
          FUN = function(x) { c(m = mean(x), s = sd(x)) } )
summaryBy(eve ~ Tissue, data = df,
          FUN = function(x) { c(m = mean(x), s = sd(x)) } )

# Plot data per Pathology variable
variabile_plot = "Pathology"
p_rich <- ggviolin(data = df, x = variabile_plot, y = "rich", add = "boxplot", fill = variabile_plot, palette ="npg")
p_rich <- p_rich + stat_compare_means(label.y = 650)                   # Add global p-value
p_shan <- ggviolin(data = df, x = variabile_plot, y = "shan", add = "boxplot", fill = variabile_plot, palette ="npg")
p_shan <- p_shan + stat_compare_means(label.y = 7)                   # Add global p-value
p_eve <- ggviolin(data = df, x = variabile_plot, y = "eve", add = "boxplot", fill = variabile_plot, palette ="npg")
p_eve <- p_eve + stat_compare_means(label.y = 1)                   # Add global p-value

plot_grid(p_rich, p_shan, p_eve, ncol = 3)

# Plot data per Tissue variable  --> Supplementary figure 3
variabile_plot = "Tissue"
my_comparisons <- list(c("CM","HM"),c("CM","LM"),c("HM","LM"))
p_rich <- ggboxplot(data = df, x = variabile_plot, y = "rich", add = "jitter", fill = variabile_plot, palette ="jama", ylab = "Observed Richness", xlab= FALSE)
p_rich <- p_rich + stat_compare_means(label.y = 650, method = "anova")                   # Add global p-value
p_rich <- p_rich + stat_compare_means(comparisons = my_comparisons, method = "t.test")
p_shan <- ggboxplot(data = df, x = variabile_plot, y = "shan", add = "jitter", fill = variabile_plot, palette ="jama", ylab = "Shannon's Index", xlab= FALSE)
p_shan <- p_shan + stat_compare_means(label.y = 7,method = "anova")                   # Add global p-value
p_shan <- p_shan + stat_compare_means(comparisons = my_comparisons, method = "t.test")
p_eve <- ggboxplot(data = df, x = variabile_plot, y = "eve", add = "jitter", fill = variabile_plot, palette ="jama",  ylab = "Pielou's Evenness", xlab= FALSE)
p_eve <- p_eve + stat_compare_means(label.y = 0.95, method = "anova")                   # Add global p-value
p_eve <- p_eve + stat_compare_means(comparisons = my_comparisons, method = "t.test")

plot_grid(p_rich, p_shan, p_eve, ncol = 3)


#### BETA DIVERSITY --> Figure 5B

set.seed(123456)
PCoA_ordi <- ordinate(qiime_file_trasf_css, "PCoA", distance = "bray")
p_PCoA_ordi <- plot_ordination(physeq = qiime_file_trasf_css, PCoA_ordi, axes = c(1,2))
p <- p_PCoA_ordi + geom_point(size = 4, stroke = 0.5, pch =21, aes(fill = Tissue))
p <- p + theme_bw() + 
  theme( panel.grid.major =  element_blank(),
         panel.grid.minor =  element_blank())
p<- p + ggtitle("PCoA; Bray-Curtis") + scale_fill_jama() + xlab("PCoA 1 (21.5%)") + ylab("PCoA 2 (9.8%)") +
  guides(fill=guide_legend(title="Tissue"))
p$layers <- p$layers[-1] 
p



## ORDINAZIONE

phylo_cloud <- subset_samples(qiime_file_trasf_css, Tissue != "HM")
sample_data(phylo_cloud)$CLUSTER <- as.factor(sample_data(phylo_cloud)$CLUSTER)
sample_data(phylo_cloud)$Tissue <- as.factor(sample_data(phylo_cloud)$Tissue)
PCoA_ordi <- ordinate(phylo_cloud, "PCoA", distance = "bray")
p_PCoA_ordi <- plot_ordination(physeq = phylo_cloud, PCoA_ordi, axes = c(1,2), color= "CLUSTER", shape = "Tissue")
p <- p_PCoA_ordi + geom_point(size = 4)
p <- p + theme_bw() + 
  theme( panel.grid.major =  element_blank(),
         panel.grid.minor =  element_blank())
p<- p + scale_color_aaas()
p

## ORDINATION ANALYSIS ON VARIOUS SAMPLES CLUSTER COMBINATION

phyloseq_LM <- subset_samples(qiime_file_trasf_css, Tissue == "LM")
phyloseq_LM <- subset_samples(phyloseq_LM, CLUSTER != "4")
fact <- as.factor(sample_data(phyloseq_LM)$CLUSTER)
fact <- revalue(fact, c("2"="Low (2+3)", "3"="Low (2+3)", "5"="High (5+6)", "6"="High (5+6)"))
sample_data(phyloseq_LM)$CLUSTER <- fact
sample_data(phyloseq_LM)$CLUSTER <- as.factor(sample_data(phyloseq_LM)$CLUSTER)
sample_data(phyloseq_LM)$Tissue <- as.factor(sample_data(phyloseq_LM)$Tissue)
PCoA_ordi <- ordinate(phyloseq_LM, "PCoA", distance = "bray")
p_PCoA_ordi <- plot_ordination(physeq = phyloseq_LM, PCoA_ordi, axes = c(1,2), color= "CLUSTER", shape = "Tissue")
p <- p_PCoA_ordi + geom_point(aes(fill=CLUSTER),colour="black",pch=21, size=4)
p <- p + theme_bw() + 
  theme( panel.grid.major =  element_blank(),
         panel.grid.minor =  element_blank())
p<- p + scale_color_npg()
p

phyloseq_LM <- subset_samples(qiime_file_trasf_css, Tissue == "LM")
phyloseq_LM <- subset_samples(phyloseq_LM, CLUSTER != "4")
sample_data(phyloseq_LM)$CLUSTER <- as.factor(sample_data(phyloseq_LM)$CLUSTER)
sample_data(phyloseq_LM)$Tissue <- as.factor(sample_data(phyloseq_LM)$Tissue)
PCoA_ordi <- ordinate(phyloseq_LM, "PCoA", distance = "bray")
p_PCoA_ordi <- plot_ordination(physeq = phyloseq_LM, PCoA_ordi, axes = c(1,2), color= "CLUSTER", shape = "Tissue")
p <- p_PCoA_ordi + geom_point(aes(fill=CLUSTER),colour="black",pch=21, size=4)
p <- p + theme_bw() + 
  theme( panel.grid.major =  element_blank(),
         panel.grid.minor =  element_blank())
p<- p + scale_color_npg()
p

phyloseq_LM <- subset_samples(qiime_file_trasf_css, Tissue == "LM")
phyloseq_LM <- subset_samples(phyloseq_LM, CLUSTER != "4")
phyloseq_LM <- subset_samples(phyloseq_LM, CLUSTER == "5" | CLUSTER == "3")
sample_data(phyloseq_LM)$CLUSTER <- as.factor(sample_data(phyloseq_LM)$CLUSTER)
sample_data(phyloseq_LM)$Tissue <- as.factor(sample_data(phyloseq_LM)$Tissue)
PCoA_ordi <- ordinate(phyloseq_LM, "PCoA", distance = "bray")
p_PCoA_ordi <- plot_ordination(physeq = phyloseq_LM, PCoA_ordi, axes = c(1,2), color= "CLUSTER", shape = "Tissue")
p <- p_PCoA_ordi + geom_point(aes(fill=CLUSTER),colour="black",pch=21, size=4)
p <- p + theme_bw() + 
  theme( panel.grid.major =  element_blank(),
         panel.grid.minor =  element_blank())
p<- p + scale_color_npg()
p

phyloseq_LM <- subset_samples(qiime_file_trasf_css, Tissue == "LM")
phyloseq_LM <- subset_samples(phyloseq_LM, CLUSTER != "4")
phyloseq_LM <- subset_samples(phyloseq_LM, CLUSTER == "6" | CLUSTER == "2")
sample_data(phyloseq_LM)$CLUSTER <- as.factor(sample_data(phyloseq_LM)$CLUSTER)
sample_data(phyloseq_LM)$Tissue <- as.factor(sample_data(phyloseq_LM)$Tissue)
PCoA_ordi <- ordinate(phyloseq_LM, "PCoA", distance = "bray")
p_PCoA_ordi <- plot_ordination(physeq = phyloseq_LM, PCoA_ordi, axes = c(1,2), color= "CLUSTER", shape = "Tissue")
p <- p_PCoA_ordi + geom_point(aes(fill=CLUSTER),colour="black",pch=21, size=4)
p <- p + theme_bw() + 
  theme( panel.grid.major =  element_blank(),
         panel.grid.minor =  element_blank())
p<- p + scale_color_npg()
p

phyloseq_LM <- subset_samples(qiime_file_trasf_css, Tissue == "LM")
phyloseq_LM <- subset_samples(phyloseq_LM, CLUSTER != "4")
phyloseq_LM <- subset_samples(phyloseq_LM, CLUSTER == "3" | CLUSTER == "2")
sample_data(phyloseq_LM)$CLUSTER <- as.factor(sample_data(phyloseq_LM)$CLUSTER)
sample_data(phyloseq_LM)$Tissue <- as.factor(sample_data(phyloseq_LM)$Tissue)
PCoA_ordi <- ordinate(phyloseq_LM, "PCoA", distance = "bray")
p_PCoA_ordi <- plot_ordination(physeq = phyloseq_LM, PCoA_ordi, axes = c(1,2), color= "CLUSTER", shape = "Tissue")
p <- p_PCoA_ordi + geom_point(aes(fill=CLUSTER),colour="black",pch=21, size=4)
p <- p + theme_bw() + 
  theme( panel.grid.major =  element_blank(),
         panel.grid.minor =  element_blank())
p<- p + scale_color_npg()
p

phyloseq_LM <- subset_samples(qiime_file_trasf_css, Tissue == "LM")
phyloseq_LM <- subset_samples(phyloseq_LM, CLUSTER != "4")
phyloseq_LM <- subset_samples(phyloseq_LM, CLUSTER == "5" | CLUSTER == "6")
sample_data(phyloseq_LM)$CLUSTER <- as.factor(sample_data(phyloseq_LM)$CLUSTER)
sample_data(phyloseq_LM)$Tissue <- as.factor(sample_data(phyloseq_LM)$Tissue)
PCoA_ordi <- ordinate(phyloseq_LM, "PCoA", distance = "bray")
p_PCoA_ordi <- plot_ordination(physeq = phyloseq_LM, PCoA_ordi, axes = c(1,2), color= "CLUSTER", shape = "Tissue")
p <- p_PCoA_ordi + geom_point(aes(fill=CLUSTER),colour="black",pch=21, size=4)
p <- p + theme_bw() + 
  theme( panel.grid.major =  element_blank(),
         panel.grid.minor =  element_blank())
p<- p + scale_color_npg()
p


#### COMMUNITY COMPOSITION --> Figure 5A

library(pals)
qiime_file_16S_Tissue <- merge_samples(x = qiime_file_trasf_css, group = "Tissue", fun = mean)
qiime_file_16S_Tissue_collapsed <- tax_glom(qiime_file_16S_Tissue, taxrank="Family")
qiime_file_16S_Tissue_collapsed_proportional <- transform_sample_counts(qiime_file_16S_Tissue_collapsed, function(x) 100 * x/sum(x))
qiime_file_16S_Tissue_collapsed_proportional_oneperc <-  filter_taxa(qiime_file_16S_Tissue_collapsed_proportional, function(x) max(x) > 1, TRUE) ## facendo cosi non è ideale
p <- plot_bar(qiime_file_16S_Tissue_collapsed_proportional_oneperc, fill= "Family",  facet_grid = ~Phylum)
p <- p + scale_fill_manual(values=as.character(polychrome(16))) +theme_linedraw() + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) 

p

#### BUILD AND EXPORT OTU TABLES FOR LEFSE --> Figure 5C and D

otu_table_lefse <- otu_table(tax_glom(subset_samples(qiime_file_trasf_css, Tissue == "LM"), taxrank = "Genus"))
tax_table_lefse <- tax_table(tax_glom(subset_samples(qiime_file_trasf_css, Tissue == "LM"), taxrank = "Genus"))
out_lefse <- cbind(otu_table_lefse@.Data, tax_table_lefse@.Data)
#write.table(x = out_lefse, file = "./lefse_final/lefse_LM/input_lefse_LM.txt", sep = "\t")
metadata_lefse <- sample_data(subset_samples(qiime_file_trasf_css, Tissue == "LM"))
#write.table(x = metadata_lefse, file = "./lefse_final/lefse_LM/metadata_lefse_LM.txt", sep = "\t")

#### BUILD SUPPLEMENTARY FIGURE 5

# ORDINATION WITH POINTS COLORED ON Il17
PCoA_ordi <- ordinate(qiime_file_trasf_css, "PCoA", distance = "bray")
p_PCoA_ordi <- plot_ordination(physeq = qiime_file_trasf_css, PCoA_ordi, axes = c(1,2), color = "IL17.qPCR.")
p <- p_PCoA_ordi + geom_point(size = 4, stroke = 0.5, pch =20)
p <- p + theme_bw() + 
  theme( panel.grid.major =  element_blank(),
         panel.grid.minor =  element_blank())
p<- p + scale_color_viridis_c() + xlab("PCoA 1 (21.5%)") + ylab("PCoA 2 (9.8%)") +
  guides(fill=guide_legend(title="IL-17")) +
  annotate("text", x = -0.1, y = 0.25, label = "PERMANOVA p-val = 0.7789")
p$layers <- p$layers[-1] 
p
#PERMANOVA TO TEST
set.seed(12387)
datas <- as(sample_data(qiime_file_trasf_css), "data.frame")
datas$IL17.qPCR.
adonis(phyloseq::distance(subset_samples(qiime_file_trasf_css, IL17.qPCR. != "NA"), method = "bray") ~ IL17.qPCR., data = datas, permutations = 9999)

# ORDINATION WITH POINTS COLORED ON Il13
set.seed(123456)
PCoA_ordi <- ordinate(qiime_file_trasf_css, "PCoA", distance = "bray")
p_PCoA_ordi <- plot_ordination(physeq = qiime_file_trasf_css, PCoA_ordi, axes = c(1,2), color = "IL13.qPCR.")
p <- p_PCoA_ordi + geom_point(size = 4, stroke = 0.5, pch =20)
p <- p + theme_bw() + 
  theme( panel.grid.major =  element_blank(),
         panel.grid.minor =  element_blank())
p<- p + scale_color_viridis_c() + xlab("PCoA 1 (21.5%)") + ylab("PCoA 2 (9.8%)") +
  guides(fill=guide_legend(title="IL-13")) + 
  annotate("text", x = -0.1, y = 0.25, label = "PERMANOVA p-val = 0.6224")
p$layers <- p$layers[-1] 
p

#PERMANOVA TO TEST
set.seed(12387)
datas <- as(sample_data(qiime_file_trasf_css), "data.frame")
adonis(phyloseq::distance(subset_samples(qiime_file_trasf_css, IL13.qPCR. != "NA"), method = "bray") ~ IL13.qPCR., data = datas, permutations = 9999)


## PERMANOVA ANALYSIS ON VARIABLES ---> Supplementary table 4

set.seed(12387)
datas <- as(sample_data(subset_samples(qiime_file_trasf_css, Tissue == "LM" & IL17.qPCR. != "NA")), "data.frame")
subset <- subset_samples(qiime_file_trasf_css, Tissue == "LM" & IL17.qPCR. != "NA")
adonis(phyloseq::distance(subset, "bray") ~ IL17.qPCR., data = datas, permutations = 9999)

set.seed(12387)
datas <- as(sample_data(subset_samples(qiime_file_trasf_css, Tissue == "LM" & IL13.qPCR. != "NA")), "data.frame")
subset <- subset_samples(qiime_file_trasf_css, Tissue == "LM" & IL13.qPCR. != "NA")
adonis(phyloseq::distance(subset, "bray") ~ IL13.qPCR., data = datas, permutations = 9999)

set.seed(12387)
datas <- as(sample_data(subset_samples(qiime_file_trasf_css, Tissue == "LM" & Mayo_endo != "NA")), "data.frame")
subset <- subset_samples(qiime_file_trasf_css, Tissue == "LM" & Mayo_endo != "NA")
adonis(phyloseq::distance(subset, "bray") ~ Mayo_endo, data = datas, permutations = 9999)

set.seed(12387)
datas <- as(sample_data(subset_samples(qiime_file_trasf_css, Tissue == "LM" & Eta_endoscopia != "NA")), "data.frame")
subset <- subset_samples(qiime_file_trasf_css, Tissue == "LM" & Eta_endoscopia != "NA")
adonis(phyloseq::distance(subset, "bray") ~ Eta_endoscopia, data = datas, permutations = 9999)

set.seed(12387)
datas <- as(sample_data(subset_samples(qiime_file_trasf_css, Tissue == "LM" & sesso != "NA")), "data.frame")
subset <- subset_samples(qiime_file_trasf_css, Tissue == "LM" & sesso != "NA")
adonis(phyloseq::distance(subset, "bray") ~ sesso, data = datas, permutations = 9999)

set.seed(12387)
datas <- as(sample_data(subset_samples(qiime_file_trasf_css, Tissue == "LM" & CLUSTER != "NA")), "data.frame")
subset <- subset_samples(qiime_file_trasf_css, Tissue == "LM" & CLUSTER != "NA")
adonis(phyloseq::distance(subset, "bray") ~ as.factor(CLUSTER), data = datas, permutations = 9999)

## FOCUS ON CLUSTER 6 AND CLUSTER 2

CL_6_2_css <- subset_samples(qiime_file_trasf_css,CLUSTER == "2" | CLUSTER == "6")
sample_data(CL_6_2_css)$CLUSTER <- as.factor(sample_data(CL_6_2_css)$CLUSTER )

PCoA_ordi <- ordinate(CL_6_2_css, "PCoA", distance = "bray")
p_PCoA_ordi <- plot_ordination(physeq = CL_6_2_css, PCoA_ordi, axes = c(1,2), color= "CLUSTER")
p <- p_PCoA_ordi + geom_point(size = 2, pch=21, colour = "black")
p <- p + theme_bw() + 
  theme( panel.grid.major =  element_blank(),
         panel.grid.minor =  element_blank())
p<- p + scale_color_manual(values = c("red", "darkgreen"))
p

## PLOT SAMPLES ON IL17/IL13 PLANE ---> supplementary figure 6

df <- as(sample_data(subset_samples(qiime_file_trasf_css, Tissue == "LM")), "data.frame")
df_nona <- df[!is.na(df$IL17.qPCR.),]
df_nona <- df[!is.na(df$IL13.qPCR.),]
df_nona <- df[!is.na(df$CLUSTER),]

df_nona["SP.27.LM",21] <- 5
df_nona["POL.27.LM",21] <- 3
df_nona["POL.62.LM",21] <- 6
p <- ggplot(df_nona, aes(x=IL17.qPCR., y=IL13.qPCR., col= as.factor(CLUSTER), shape= Tissue, label=rownames(df_nona))) + 
  geom_hline(yintercept = mean(df_nona$IL13.qPCR.), linetype="dotted") + 
  geom_vline(xintercept = mean(df_nona$IL17.qPCR.), linetype="dotted") +
  geom_point(size=4) + theme_linedraw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  scale_color_npg()
p 

## BARPLOT FOR EACH CLUSTER ---> supplementary figure 6

phylo_cluster <- subset_samples(qiime_file_trasf_css, Tissue == "LM")

library(pals)
aaa <- subset_samples(phylo_cluster, CLUSTER != "NA")
paletta <- as.character(polychrome(29))
qiime_file_16S_Tissue <- merge_samples(x = aaa, group = "CLUSTER", fun = mean)
qiime_file_16S_Tissue_collapsed <- tax_glom(qiime_file_16S_Tissue, taxrank="Family")
qiime_file_16S_Tissue_collapsed_proportional <- transform_sample_counts(qiime_file_16S_Tissue_collapsed, function(x) 100 * x/sum(x))
qiime_file_16S_Tissue_collapsed_proportional_oneperc <-  filter_taxa(qiime_file_16S_Tissue_collapsed_proportional, function(x) max(x) > 1, TRUE)
p <- plot_bar(qiime_file_16S_Tissue_collapsed_proportional_oneperc, fill= "Family")
p + scale_fill_manual(values=paletta) + theme_linedraw() + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) 

library(pals)
aaa <- subset_samples(phylo_cluster, CLUSTER != "NA")
paletta <- rev(as.character(alphabet(26)))
qiime_file_16S_Tissue <- merge_samples(x = aaa, group = "CLUSTER", fun = mean)
qiime_file_16S_Tissue_collapsed <- tax_glom(qiime_file_16S_Tissue, taxrank="Genus")
qiime_file_16S_Tissue_collapsed_proportional <- transform_sample_counts(qiime_file_16S_Tissue_collapsed, function(x) 100 * x/sum(x))
qiime_file_16S_Tissue_collapsed_proportional_oneperc <-  filter_taxa(qiime_file_16S_Tissue_collapsed_proportional, function(x) max(x) > 1, TRUE)
p <- plot_bar(qiime_file_16S_Tissue_collapsed_proportional_oneperc, fill= "Genus")
p + scale_fill_manual(values=paletta) + theme_linedraw() + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) 


library(pals)
aaa <- subset_samples(phylo_cluster, CLUSTER != "NA")
paletta <- as.character(glasbey(26))
qiime_file_16S_Tissue <- merge_samples(x = aaa, group = "CLUSTER", fun = mean)
qiime_file_16S_Tissue_collapsed <- tax_glom(qiime_file_16S_Tissue, taxrank="Species")
qiime_file_16S_Tissue_collapsed_proportional <- transform_sample_counts(qiime_file_16S_Tissue_collapsed, function(x) 100 * x/sum(x))
qiime_file_16S_Tissue_collapsed_proportional_oneperc <-  filter_taxa(qiime_file_16S_Tissue_collapsed_proportional, function(x) max(x) > 1, TRUE)
p <- plot_bar(qiime_file_16S_Tissue_collapsed_proportional_oneperc, fill= "Species")
p + scale_fill_manual(values=paletta) + theme_linedraw() + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) 


## ORDINATION WITH CLUSTER ---> supplementary figure 6

phylo_cluster_CMLM <- subset_samples(qiime_file_trasf_css, Tissue != "HM")
sample_data(phylo_cluster_CMLM)$CLUSTER <- as.factor(sample_data(phylo_cluster_CMLM)$CLUSTER)

PCoA_ordi <- ordinate(phylo_cluster_CMLM, "PCoA", distance = "bray")
p_PCoA_ordi <- plot_ordination(physeq = phylo_cluster_CMLM, PCoA_ordi, axes = c(1,2), color="CLUSTER", shape = "Tissue")
p <- p_PCoA_ordi + geom_point(size = 3)
p <- p + theme_bw() + 
  theme( panel.grid.major =  element_blank(),
         panel.grid.minor =  element_blank())
p<- p + ggtitle("PCoA; Bray-Curtis") + scale_color_aaas()
p


