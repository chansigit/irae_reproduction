################################################################################
#
#    Analysis
#    PRIME-TR  MD Anderson
#
#
################################################################################

## R package load:
rm(list = ls())
set.seed(123)
library(phyloseq)
library(tidyr)
library(tidyselect)
library(tidyverse)
library(ggplot2)
library(phyloseq)
library(readxl)
library(ComplexHeatmap)
library(pheatmap)
library(microbiome)
library(patchwork)
library(vegan)
library(openxlsx)
library(ggh4x)
library(Cairo)
library(randomcoloR)
library(magrittr)
library(ggpubr)
library(rstatix)
library(lmerTest)
library(lme4)
library(DESeq2)


# data folder
dfol <- 'data'

# color scheme
BaRU <- c('orange', 'royalblue')
RNR <- c('seagreen','gold')
BaRU <- c('darkolivegreen2','deeppink3')

## prepare file for analysis
microbiome_raw = import_biom(file.path(dfol,"WatowichResubmission.03.04.2022.biom"))
colnames(tax_table(microbiome_raw)) <-
  c("Kingdom", "Phylum",  "Class" ,  "Order" ,  "Family" , "Genus")

## Read meta data
temp1 <- read.csv(file.path(dfol,"Yifan.20211215.C3KO.IL6.Abx.level.abundance.metadata.csv")) %>%
  column_to_rownames('SampleID') %>%
  select(Sample.description) %>%
  dplyr::rename(Treatment = Sample.description) %>%
  mutate(Treatment = paste('C3KO',Treatment,sep = '+'))

## Add sample data to phyloseq object
sample_data(microbiome_raw) <- sample_data(temp1)

## compute alpha diversity using the observed abundance data
microbiome_raw <- subset_taxa(microbiome_raw, taxa_sums(microbiome_raw) >0 )
richness_select <- 'InvSimpson'
sample_data(microbiome_raw)$richness <- estimate_richness(microbiome_raw, measures=richness_select)[,richness_select]

## Prepare data for beta diversity analysis and differential abundance analysis
microbiome_raw <- subset_taxa(microbiome_raw, rowSums(otu_table(microbiome_raw) > 0 ) > 2 )


## Data analysis file
save(microbiome_raw, file = file.path(dfol,'processed_data.rds') )




## ---------------  Composition plot --------------- -------------------------
## Make composition plot
##
mbiome_subset <- microbiome_raw
agg_level <- 'Family' # c("Kingdom", "Phylum",  "Class" ,  "Order" ,  "Family" , "Genus",'Species')
mbiome_subset_fam <- microbiome::aggregate_taxa(mbiome_subset, agg_level)
mbiome_subset_fam.rel <- microbiome::transform(mbiome_subset_fam, "compositional")
seltaxa <- apply(otu_table(mbiome_subset_fam.rel) %>% as.data.frame(),2,
                 function(x) taxa_names(mbiome_subset_fam)[order(x, decreasing = T)[1:4]] ) %>%
  as.vector() %>% unique()

mbiome_subset_fam.rel <- subset_taxa(mbiome_subset_fam.rel,
                                     taxa_names(mbiome_subset_fam.rel) %in%
                                       seltaxa)
temp <- taxa_sums(mbiome_subset_fam.rel)
seltaxa <- names(temp)[order(temp, decreasing = T)]
pldata <- otu_table(mbiome_subset_fam.rel) %>% t() %>% as.data.frame() %>%
  cbind(.,Others = 1- sample_sums(otu_table(mbiome_subset_fam.rel))) %>%
  mutate(Sample = rownames(.))

pldata <- pldata %>% cbind(., sample_data(mbiome_subset_fam.rel)) %>%
  gather(key = 'Tax',value = "Abundance", -(Sample:richness)) %>%
  mutate(Abundance  = as.numeric(Abundance))


## implement sample sort and sample.sort = "Firmicutes", otu.sort = "abundance"
set.seed(123)
palette <- distinctColorPalette(length(seltaxa) + 1)[(length(seltaxa) + 1):1]

pldata4 <- pldata %>% mutate(Tax = factor(Tax, levels = c(seltaxa, 'Others') )) %>%
  mutate(Marker  = paste(Treatment, sep = '')) %>%
  arrange(Marker,Tax,-Abundance) %>%
  mutate(Sample = factor(Sample, levels = unique(Sample) ))


comp_plot <- ggplot(pldata4, aes(x = Sample, y = Abundance,
                                 fill = Tax  )) +
  geom_bar(position = "stack", stat = "identity") +
  scale_x_discrete() + # labels = pldata$xlabel, breaks = pldata$Sample
  facet_nested(~ Marker, scales = "free", space = 'free') +
  theme_bw() + #scale_fill_brewer("Family", palette = "Paired") +
  scale_fill_manual(agg_level, values = palette) +
  xlab('Samples') + ylab('Relative abundance') +
  theme( # legend.title = element_text(size = 18),
    panel.spacing=unit(0,"lines"),
    strip.background = element_blank(),
    axis.text.x = element_blank(),
    legend.position = 'top',
    plot.title = element_text(size = 15, hjust = 0.5, color = 'black'),
    axis.text.y = element_text(size = 13, color = 'black'),
    axis.title = element_text(size = 13, color = 'black'),
    strip.text = element_text(size = 11, color = 'black', angle = 0),
    legend.text = element_text(size = 10, color = 'black'),
    legend.background = element_rect(fill="grey95",
                                     size=0.5, linetype="solid"),
    legend.title =  element_text(size = 13, color = 'black',
                                 angle = 90, hjust = 0.5) )





## ---------------  Composition plot --------------- -------------------------
## Make composition plot
##
mbiome_subset <- microbiome_raw
agg_level <- 'Genus' # c("Kingdom", "Phylum",  "Class" ,  "Order" ,  "Family" , "Genus",'Species')
mbiome_subset_fam <- microbiome::aggregate_taxa(mbiome_subset, agg_level)
mbiome_subset_fam.rel <- microbiome::transform(mbiome_subset_fam, "compositional")
seltaxa <- apply(otu_table(mbiome_subset_fam.rel) %>% as.data.frame(),2,
                 function(x) taxa_names(mbiome_subset_fam)[order(x, decreasing = T)[1:4]] ) %>%
  as.vector() %>% unique()

mbiome_subset_fam.rel <- subset_taxa(mbiome_subset_fam.rel,
                                     taxa_names(mbiome_subset_fam.rel) %in%
                                       seltaxa)
temp <- taxa_sums(mbiome_subset_fam.rel)
seltaxa <- names(temp)[order(temp, decreasing = T)]
pldata <- otu_table(mbiome_subset_fam.rel) %>% t() %>% as.data.frame() %>%
  cbind(.,Others = 1- sample_sums(otu_table(mbiome_subset_fam.rel))) %>%
  mutate(Sample = rownames(.))

pldata <- pldata %>% cbind(., sample_data(mbiome_subset_fam.rel)) %>%
  gather(key = 'Tax',value = "Abundance", -(Sample:richness)) %>%
  mutate(Abundance  = as.numeric(Abundance))


## implement sample sort and sample.sort = "Firmicutes", otu.sort = "abundance"
set.seed(123)
palette <- distinctColorPalette(length(seltaxa) + 1)[(length(seltaxa) + 1):1]

pldata4 <- pldata %>% mutate(Tax = factor(Tax, levels = c(seltaxa, 'Others') )) %>%
  mutate(Marker  = paste(Treatment, sep = '')) %>%
  arrange(Marker,Tax,-Abundance) %>%
  mutate(Sample = factor(Sample, levels = unique(Sample) ))


comp_plot1 <- ggplot(pldata4, aes(x = Sample, y = Abundance,
                                  fill = Tax  )) +
  geom_bar(position = "stack", stat = "identity") +
  scale_x_discrete() + # labels = pldata$xlabel, breaks = pldata$Sample
  facet_nested(~ Marker, scales = "free", space = 'free') +
  theme_bw() + #scale_fill_brewer("Family", palette = "Paired") +
  scale_fill_manual(agg_level, values = palette) +
  xlab('Samples') + ylab('Relative abundance') +
  theme( # legend.title = element_text(size = 18),
    panel.spacing=unit(0,"lines"),
    strip.background = element_blank(),
    axis.text.x = element_blank(),
    legend.position = 'top',
    plot.title = element_text(size = 15, hjust = 0.5, color = 'black'),
    axis.text.y = element_text(size = 13, color = 'black'),
    axis.title = element_text(size = 13, color = 'black'),
    strip.text = element_text(size = 11, color = 'black', angle = 0),
    legend.text = element_text(size = 10, color = 'black'),
    legend.background = element_rect(fill="grey95",
                                     size=0.5, linetype="solid"),
    legend.title =  element_text(size = 13, color = 'black',
                                 angle = 90, hjust = 0.5) )

# ggpubr::ggarrange(comp_plot, comp_plot1)

setEPS()
postscript('./plots/resubmission/supFig_comp_f1.eps',  width = 14.7, height = 7.26)
comp_plot
dev.off()
setEPS()
postscript('./plots/resubmission/supFig_comp_f2.eps',  width = 14.7, height = 7.26)
comp_plot1
dev.off()

#




## ---------------  Plot for alpha diversity analysis ---------------
mbiome_subset <- microbiome_raw
metadf <- sample_data(mbiome_subset) %>% data.frame() %>%
  select(c('Treatment',"richness"))

## Box plot across points:
stat.test <- metadf %>%
  t_test( formula = richness ~ Treatment, data = .) %>%
  add_significance() %>%
  add_xy_position(step.increase = .1) %>%
  mutate(p = paste('P-value = ',p))

plot_view1 <- ggplot(metadf, aes(x = Treatment, y = richness)) +
  geom_boxplot(aes(fill = Treatment), alpha = 1, notch = F) +
  # geom_violin(draw_quantiles = c(0.25, 0.5, 0.75)) +
  # geom_line(aes(color = SubjectID, group = SubjectID),size = 1) +
  # geom_jitter(aes(color = Treatment),width = 0.15, size = 3) +
  geom_jitter(width = 0.15, size = 3) +
  xlab('Treatment') +
  stat_pvalue_manual(stat.test, label = "p", tip.length = 0.05, size= 5) +
  ylab('Inverse simpson') + ylim(c(0, 30)) +
  ggtitle('Diversity analysis') +
  # guides(color = guide_legend(override.aes = list(size = 3) ) ) +
  theme_classic() + scale_fill_manual(values = BaRU)+
  scale_color_manual(values = BaRU)+
  theme( panel.spacing=unit(0,"lines"),
         strip.background = element_blank(),
         legend.position = 'none',
         plot.title = element_text(size = 15, hjust = 0.5, color = 'black'),
         axis.text = element_text(size = 14, color = 'black'),
         axis.title = element_text(size = 15, color = 'black'),
         strip.text = element_text(size = 13, color = 'black', angle = 0),
         legend.text = element_text(size = 10, color = 'black'),
         legend.background = element_rect(fill="grey95",
                                          size=0.5, linetype="solid"),
         legend.title =  element_text(size = 13, color = 'black', face = 'bold',
                                      angle = 0, hjust = 0.5) )

plot_view1

setEPS()
postscript('./plots/resubmission/supFig_alpha_d.eps',  width = 18, height = 6)
plot_view1
dev.off()






## --------------- Beta diversity analysis and dispersion test  ---------
## We prepare the data table for plotting


dt_comp <- microbiome::transform(mbiome_subset, "compositional")
dt_comp_c <- ordinate(dt_comp , "PCoA", "bray", weighted=T)
ax_name <- paste(c('PCo1', 'PCo2'), ' [',
                 signif(dt_comp_c$values$Relative_eig * 100,2)[1:2],'%]',
                 sep='')

metadf <- data.frame(sample_data(dt_comp))
dist_bray <- vegdist(otu_table(mbiome_subset)@.Data %>% t(), method = "bray")
bd <- betadisper(dist_bray, metadf$Treatment, type = 'centroid')
bd

set.seed(123)
permAnovaTest <- adonis(
  otu_table(mbiome_subset)@.Data %>% t() ~  Treatment, data = metadf,
  permutations = 999,
  method = "bray"
)
permAnovaTest

pldata <- dt_comp_c$vectors %>% data.frame() %>% cbind(.,metadf) %>%
  cbind(.,bd$centroids[.$Treatment,1:2])
pldata$Name <- substring(rownames(pldata),20)


plot_view2 <- ggplot(pldata , aes(x =  Axis.1, y = Axis.2, color = Treatment)) +
  geom_point(size = 5) +
  # ggrepel::geom_text_repel() +,label = Name
  # geom_text(hjust=0, vjust=0) +
  theme_classic() +
  geom_segment(aes(x = PCoA1, xend=Axis.1, y=PCoA2, yend=Axis.2, color = Treatment),
               arrow = arrow(length = unit(0.25, "cm")),lwd=0.3) +
  scale_color_manual("Treatment", values =  BaRU) +
  ggtitle("Bray-Curtis dissimilarity") +
  xlab(ax_name[1]) + ylab(ax_name[2]) +
  annotate("text", x=0.2, y=-0.4, size =5,
           label= paste("P-value (Treatment) =",signif(permAnovaTest$aov.tab$`Pr(>F)`[1],1))) +
  theme(legend.position = 'top',
        plot.title = element_text(size = 15, hjust = 0.5, color = 'black'),
        axis.text = element_text(size = 13, color = 'black'),
        axis.title = element_text(size = 13, color = 'black'),
        strip.text = element_text(size = 15, color = 'black'),
        legend.text = element_text(size = 11, color = 'black'),
        legend.background = element_rect(fill="grey99",
                                         size=0.5, linetype="solid"),
        legend.title =  element_blank() )

plot_view2

setEPS()
postscript('./plots/resubmission/mainFig_beta_h.eps',  width = 18, height = 6)
plot_view2
dev.off()



## ---- Differential abundance analysis ----
##
set.seed(123)
analysis_data <- mbiome_subset
analysis_data <- microbiome::aggregate_taxa(analysis_data, 'Family')

# deseq diffential abundance analysis
diagdds = phyloseq_to_deseq2(analysis_data, ~ Treatment )
diagdds = DESeq(diagdds, test="LRT", fitType="parametric",
                sfType = c("ratio", "poscounts", "iterate")[3],
                minReplicatesForReplace=Inf, reduced = ~1)
res = results(diagdds, cooksCutoff = FALSE)
# significant taxa table
alpha =  0.1    ## subjected to change
sigtab = res[which((res$pvalue < alpha)), ]
dim(sigtab)
sigtab = cbind(as(sigtab, "data.frame"), as(tax_table(analysis_data)[rownames(sigtab), ], "matrix"))
sigtab$color<-ifelse(sigtab$log2FoldChange >0, "Up",  "Down")
head(sigtab)

write_csv(sigtab, file = './plots/resubmission/main_diffabu1.csv')

lab_names <- sigtab$Family   ## save the label names

pldata <-data.frame(gen=lab_names,value=sigtab$log2FoldChange,
                    col=sigtab$color, index = 1:nrow(sigtab),
                    ID = rownames(sigtab))
x_lable <- pldata$gen[order(pldata$value)]
plot_view41 <- ggplot(data=pldata ,
                      aes(x=reorder(index, value),y=value,fill=col))+
  geom_bar(stat="identity")+ theme_bw()+
  # scale_fill_manual(values = c('red','blue') )+
  scale_fill_manual(NULL,values = BaRU)+
  xlab("Taxonomic identifier (Family)")+
  ylab("log2 \n (fold change)") + coord_flip()+
  scale_x_discrete(labels= x_lable) +
  theme(legend.position =  'top',#c(0.1,0.87),
        legend.background = element_rect(fill="transparent"),
        plot.title = element_text(size = 15, hjust = 0.5),
        plot.subtitle = element_text(size = 15, hjust = 0.5),
        axis.text = element_text(size = 11,color = 'black'),
        axis.title = element_text(size = 13),
        strip.text = element_text(size = 13),
        legend.text = element_text(size = 13),
        legend.title = element_text(size = 13) )
plot_view41


# prepare data for sister box plot
dt_comp <- analysis_data # microbiome::transform(dietStudy_genus, "compositional")
dt_comp <- subset_taxa(dt_comp, taxa_names(dt_comp) %in% rownames(sigtab))
temp <- sample_data(dt_comp)[,c('Treatment')]
temp$SampleID <- rownames(temp)
plmat<- cbind(as(counts(diagdds, normalized=T)[rownames(sigtab), ], "matrix"),
              as(tax_table(dt_comp)[rownames(sigtab), ], "matrix")) %>%
  as.data.frame()  %>%
  gather(  key = "key",value = "value", -(Kingdom:unique)) %>%
  dplyr::left_join(., temp, by = c("key" = "SampleID")) %>%
  mutate(value = as.numeric(value)) %>%
  mutate(Family = factor(Family, levels = pldata$gen[order(pldata$value)]))


plot_view42 <- ggplot(plmat %>% filter(value > 0), aes(x = Family, y = value, fill = Treatment)) +
  geom_boxplot(notch = F,lwd=0.8) +
  # geom_jitter(aes(color = Treatment),shape=1, width = 0.1,stroke = 1, size = 3) +
  scale_y_continuous(trans='log10',n.breaks=4) +
  scale_fill_manual(name = NULL, values = BaRU) +
  scale_color_manual(name = NULL, values = BaRU) +
  theme_bw() +
  ylab( 'Normalized abundance \n (log scale)') +
  xlab('') +
  coord_flip() +
  scale_x_discrete(labels = lab_names)+
  theme(legend.position = 'top',
        plot.title = element_text(size = 15, hjust = 0.5, color = 'black'),
        axis.text.x = element_text(size = 13, hjust = 1.0, vjust = 0.5,
                                   color = 'black'),
        axis.text.y = element_blank(),
        # axis.text.y = element_text(size = 13, color = 'black',),
        axis.title = element_text(size = 13, color = 'black',),
        strip.text = element_text(size = 15, color = 'black'),
        legend.text = element_text(size = 13, color = 'black'),
        legend.title = element_blank() )
plot_view42

setEPS()
postscript('./plots/resubmission/supFig_deseq_e.eps',  width = 10, height = 4)
ggpubr::ggarrange(plot_view41, plot_view42, nrow = 1, widths = c(1.2, 1.),  align = c("none", "h", "v", "hv")[2], common.legend = F)
dev.off()




## ---- Differential abundance analysis ----
##
set.seed(123)
analysis_data <- mbiome_subset
analysis_data <- microbiome::aggregate_taxa(analysis_data, 'Genus')

# deseq diffential abundance analysis
diagdds = phyloseq_to_deseq2(analysis_data, ~ Treatment )
diagdds = DESeq(diagdds, test="LRT", fitType="parametric",
                sfType = c("ratio", "poscounts", "iterate")[2],
                minReplicatesForReplace=Inf, reduced = ~1)
res = results(diagdds, cooksCutoff = FALSE)
# significant taxa table
alpha =  0.1    ## subjected to change
sigtab = res[which((res$pvalue < alpha)), ]
dim(sigtab)
sigtab = cbind(as(sigtab, "data.frame"), as(tax_table(analysis_data)[rownames(sigtab), ], "matrix"))
sigtab$color<-ifelse(sigtab$log2FoldChange >0, "Up",  "Down")
head(sigtab)

write_csv(sigtab, file = './plots/resubmission/main_diffabu2.csv')

lab_names <- sigtab$Genus   ## save the label names

pldata <-data.frame(gen=lab_names,value=sigtab$log2FoldChange,
                    col=sigtab$color, index = 1:nrow(sigtab),
                    ID = rownames(sigtab))
x_lable <- pldata$gen[order(pldata$value)]
plot_view41 <- ggplot(data=pldata ,
                      aes(x=reorder(index, value),y=value,fill=col))+
  geom_bar(stat="identity")+ theme_bw()+
  # scale_fill_manual(values = c('red','blue') )+
  scale_fill_manual(NULL,values = BaRU)+
  xlab("Taxonomic identifier (Genus)")+
  ylab("log2 \n (fold change)") + coord_flip()+
  scale_x_discrete(labels= x_lable) +
  theme(legend.position =  'top',#c(0.1,0.87),
        legend.background = element_rect(fill="transparent"),
        plot.title = element_text(size = 15, hjust = 0.5),
        plot.subtitle = element_text(size = 15, hjust = 0.5),
        axis.text = element_text(size = 11,color = 'black'),
        axis.title = element_text(size = 13),
        strip.text = element_text(size = 13),
        legend.text = element_text(size = 13),
        legend.title = element_text(size = 13) )
plot_view41


# prepare data for sister box plot
dt_comp <- analysis_data # microbiome::transform(dietStudy_genus, "compositional")
dt_comp <- subset_taxa(dt_comp, taxa_names(dt_comp) %in% rownames(sigtab))
temp <- sample_data(dt_comp)[,c('Treatment')]
temp$SampleID <- rownames(temp)
plmat<- cbind(as(counts(diagdds, normalized=T)[rownames(sigtab), ], "matrix"),
              as(tax_table(dt_comp)[rownames(sigtab), ], "matrix")) %>%
  as.data.frame()  %>%
  gather(  key = "key",value = "value", -(Kingdom:unique)) %>%
  dplyr::left_join(., temp, by = c("key" = "SampleID")) %>%
  mutate(value = as.numeric(value)) %>%
  mutate(Genus = factor(Genus, levels = pldata$gen[order(pldata$value)]))


plot_view42 <- ggplot(plmat %>% filter(value > 0), aes(x = Genus, y = value, fill = Treatment)) +
  geom_boxplot(notch = F,lwd=0.8) +
  # geom_jitter(aes(color = Treatment),shape=1, width = 0.1,stroke = 1, size = 3) +
  scale_y_continuous(trans='log10',n.breaks=4) +
  scale_fill_manual(name = NULL, values = BaRU) +
  scale_color_manual(name = NULL, values = BaRU) +
  theme_bw() +
  ylab( 'Normalized abundance \n (log scale)') +
  xlab('') +
  coord_flip() +
  scale_x_discrete(labels = lab_names)+
  theme(legend.position = 'top',
        plot.title = element_text(size = 15, hjust = 0.5, color = 'black'),
        axis.text.x = element_text(size = 13, hjust = 1.0, vjust = 0.5,
                                   color = 'black'),
        axis.text.y = element_blank(),
        # axis.text.y = element_text(size = 13, color = 'black',),
        axis.title = element_text(size = 13, color = 'black',),
        strip.text = element_text(size = 15, color = 'black'),
        legend.text = element_text(size = 13, color = 'black'),
        legend.title = element_blank() )
plot_view42

setEPS()
postscript('./plots/resubmission/mainFig_deseq_i.eps',  width = 10, height = 4)
ggpubr::ggarrange(plot_view41, plot_view42, nrow = 1, widths = c(1.2, 1.),  align = c("none", "h", "v", "hv")[2], common.legend = F)
dev.off()







