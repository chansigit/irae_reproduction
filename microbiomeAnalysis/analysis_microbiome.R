################################################################################
#
#    Analysis
#    PRIME-TR  MD Anderson: Add plot file with the script
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


pe_col <- c('deeppink','darkolivegreen2')
po_col_w <- c('darkolivegreen','black')
po_col_c <- c('darkmagenta','navy')

## Read microbiome file
microbiome_raw <- import_biom(file.path(dfol,"Watowich.irAWs.03.24.2022.biom" ))
colnames(tax_table(microbiome_raw)) <-
  c("Kingdom", "Phylum",  "Class" ,  "Order" ,  "Family" , "Genus")


## prepare file for analysis
temp <- read_xlsx(file.path(dfol,"WatowichS_Metadata_09012020.xlsx")) %>%
  mutate(SampleID2 = SampleID) %>%
  column_to_rownames('SampleID2') %>%
  data.frame() %>%
  .[sample_names(microbiome_raw), ] %>%
  sample_data()

sample_data(microbiome_raw) <- temp
## Get-rid of texa that appears in less than 10% sample
microbiome_raw <- subset_taxa(microbiome_raw, rowSums(otu_table(microbiome_raw) > 0 ) > 8 )



## compute alpha diversity using the observed abundance data
richness_select <- 'InvSimpson'
sample_data(microbiome_raw)$richness <- estimate_richness(microbiome_raw, measures=richness_select)[,richness_select]


## Data analysis file
save(microbiome_raw, file = file.path(dfol,'processed_data2.rds') )





## ---------------  Composition plot --------------- -------------------------
## Make composition plot
## Paired and unpaired post combination of the composition plots
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
seltaxa <- seltaxa[1:min(length(seltaxa),19)]
mbiome_subset_fam.rel <- subset_taxa(mbiome_subset_fam.rel,
                                     taxa_names(mbiome_subset_fam.rel) %in%
                                       seltaxa)
pldata <- otu_table(mbiome_subset_fam.rel) %>% t() %>% as.data.frame() %>%
  cbind(.,Others = 1- sample_sums(otu_table(mbiome_subset_fam.rel))) %>%
  mutate(Sample = rownames(.))

pldata <- pldata %>% cbind(., sample_data(mbiome_subset_fam.rel)) %>%
  gather(key = 'Tax',value = "Abundance", -(Sample:richness)) %>%
  mutate(Abundance  = as.numeric(Abundance))


## implement sample sort and sample.sort = "Firmicutes", otu.sort = "abundance"
set.seed(1234)
palette <- distinctColorPalette(length(seltaxa) + 1)[(length(seltaxa) + 1):1]

pldata4 <- pldata %>% mutate(Tax = factor(Tax, levels = c(seltaxa, 'Others') )) %>%
  mutate(Treatment_Timing = ifelse(Treatment_Timing == 'Post_Treatment',
                                   'Post-treat','Pre-treat')) %>%
  mutate(Status = ifelse(Status == 'Paired',
                         'Paired: Y','Paired: N')) %>%
  mutate(model_var = ifelse(Murine_Model == 'WT', 'WT', 'C3KO' )) %>%
  mutate(Treatment = ifelse(Treatment == 'CTLA4', '[CTLA-4]', '[IgG]' )) %>%
  mutate(Treatment2 = replace(Treatment, Treatment_Timing == 'Pre-treatment', '')) %>%
  mutate(Marker  = paste( model_var,'\n',
                          Treatment2,'\n', Treatment_Timing , sep = ''))  %>%

  mutate(Marker = factor(Marker,
                         levels = unique(Marker) %>% sort(decreasing = T)
  )) %>%
  arrange(Marker,Tax,-Abundance) %>%
  # arrange(Treatment_Timing,Status,Pair_Group) %>%
  mutate(Sample = factor(Sample, levels = unique(Sample) ))



comp_plot1 <- ggplot(pldata4 %>% filter(Status == 'Paired: Y'),
                     aes(x = Sample, y = Abundance, fill = Tax  )) +
  geom_bar(position = "stack", stat = "identity") +
  scale_x_discrete() + # labels = pldata$xlabel, breaks = pldata$Sample
  facet_nested(~  Marker, scales = "free", space = 'free') +
  theme_bw() + #scale_fill_brewer("Family", palette = "Paired") +
  scale_fill_manual(agg_level, values = palette) +
  xlab('Samples (paired)') + ylab('Relative abundance') +
  theme( # legend.title = element_text(size = 18),
    panel.spacing=unit(0.1,"lines"),
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

comp_plot2 <- ggplot(pldata4 %>% filter(Status != 'Paired: Y'),
                     aes(x = Sample, y = Abundance, fill = Tax  )) +
  geom_bar(position = "stack", stat = "identity") +
  scale_x_discrete() + # labels = pldata$xlabel, breaks = pldata$Sample
  facet_nested(~  Marker, scales = "free", space = 'free') +
  theme_bw() + #scale_fill_brewer("Family", palette = "Paired") +
  scale_fill_manual(agg_level, values = palette) +
  xlab('Samples (unpaired)') + ylab('Relative abundance') +
  theme( # legend.title = element_text(size = 18),
    panel.spacing=unit(0.1,"lines"),
    strip.background = element_blank(),
    axis.text.x = element_blank(),
    legend.position = 'top',
    plot.title = element_text(size = 15, hjust = 0.5, color = 'black'),
    # axis.text.y = element_text(size = 13, color = 'black'),
    axis.text.y = element_blank(),
    axis.title.y = element_blank(),
    axis.title.x = element_text(size = 13, color = 'black'),
    strip.text = element_text(size = 11, color = 'black', angle = 0),
    legend.text = element_text(size = 10, color = 'black'),
    legend.background = element_rect(fill="grey95",
                                     size=0.5, linetype="solid"),
    legend.title =  element_text(size = 13, color = 'black',
                                 angle = 90, hjust = 0.5) )

temp <- pldata4 %>% select(Status,Sample) %>% distinct()
table(temp$Status)

setEPS()
postscript('./plots/manuscript/supFig_comp_family_b1.eps',  width = 15.0, height = 7.26)
# comp_plot
ggpubr::ggarrange(comp_plot1,comp_plot2,align = 'hv', common.legend = T, widths = c(54/28,1))
dev.off()





## ---------------  Composition plot (Genus) --------------- -------------------------
## Make composition plot
## Paired and unpaired post combination of the composition plots
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
seltaxa <- seltaxa[1:min(length(seltaxa),19)]
mbiome_subset_fam.rel <- subset_taxa(mbiome_subset_fam.rel,
                                     taxa_names(mbiome_subset_fam.rel) %in%
                                       seltaxa)
pldata <- otu_table(mbiome_subset_fam.rel) %>% t() %>% as.data.frame() %>%
  cbind(.,Others = 1- sample_sums(otu_table(mbiome_subset_fam.rel))) %>%
  mutate(Sample = rownames(.))

pldata <- pldata %>% cbind(., sample_data(mbiome_subset_fam.rel)) %>%
  gather(key = 'Tax',value = "Abundance", -(Sample:richness)) %>%
  mutate(Abundance  = as.numeric(Abundance))


## implement sample sort and sample.sort = "Firmicutes", otu.sort = "abundance"
set.seed(1234)
palette <- distinctColorPalette(length(seltaxa) + 1)[(length(seltaxa) + 1):1]

pldata4 <- pldata %>% mutate(Tax = factor(Tax, levels = c(seltaxa, 'Others') )) %>%
  mutate(Treatment_Timing = ifelse(Treatment_Timing == 'Post_Treatment',
                                   'Post-treat','Pre-treat')) %>%
  mutate(Status = ifelse(Status == 'Paired',
                         'Paired: Y','Paired: N')) %>%
  mutate(model_var = ifelse(Murine_Model == 'WT', 'WT', 'C3KO' )) %>%
  mutate(Treatment = ifelse(Treatment == 'CTLA4', '[CTLA-4]', '[IgG]' )) %>%
  mutate(Treatment2 = replace(Treatment, Treatment_Timing == 'Pre-treatment', '')) %>%
  mutate(Marker  = paste( model_var,'\n',
                          Treatment2,'\n', Treatment_Timing , sep = ''))  %>%

  mutate(Marker = factor(Marker,
                         levels = unique(Marker) %>% sort(decreasing = T)
  )) %>%
  arrange(Marker,Tax,-Abundance) %>%
  # arrange(Treatment_Timing,Status,Pair_Group) %>%
  mutate(Sample = factor(Sample, levels = unique(Sample) ))



comp_plot1 <- ggplot(pldata4 %>% filter(Status == 'Paired: Y'),
                     aes(x = Sample, y = Abundance, fill = Tax  )) +
  geom_bar(position = "stack", stat = "identity") +
  scale_x_discrete() + # labels = pldata$xlabel, breaks = pldata$Sample
  facet_nested(~  Marker, scales = "free", space = 'free') +
  theme_bw() + #scale_fill_brewer("Family", palette = "Paired") +
  scale_fill_manual(agg_level, values = palette) +
  xlab('Samples (paired)') + ylab('Relative abundance') +
  theme( # legend.title = element_text(size = 18),
    panel.spacing=unit(0.1,"lines"),
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

comp_plot2 <- ggplot(pldata4 %>% filter(Status != 'Paired: Y'),
                     aes(x = Sample, y = Abundance, fill = Tax  )) +
  geom_bar(position = "stack", stat = "identity") +
  scale_x_discrete() + # labels = pldata$xlabel, breaks = pldata$Sample
  facet_nested(~  Marker, scales = "free", space = 'free') +
  theme_bw() + #scale_fill_brewer("Family", palette = "Paired") +
  scale_fill_manual(agg_level, values = palette) +
  xlab('Samples (unpaired)') + ylab('Relative abundance') +
  theme( # legend.title = element_text(size = 18),
    panel.spacing=unit(0.1,"lines"),
    strip.background = element_blank(),
    axis.text.x = element_blank(),
    legend.position = 'top',
    plot.title = element_text(size = 15, hjust = 0.5, color = 'black'),
    # axis.text.y = element_text(size = 13, color = 'black'),
    axis.text.y = element_blank(),
    axis.title.y = element_blank(),
    axis.title.x = element_text(size = 13, color = 'black'),
    strip.text = element_text(size = 11, color = 'black', angle = 0),
    legend.text = element_text(size = 10, color = 'black'),
    legend.background = element_rect(fill="grey95",
                                     size=0.5, linetype="solid"),
    legend.title =  element_text(size = 13, color = 'black',
                                 angle = 90, hjust = 0.5) )

temp <- pldata4 %>% select(Status,Sample) %>% distinct()
table(temp$Status)

setEPS()
postscript('./plots/manuscript/supFig_comp_genus_b2.eps',  width = 15.0, height = 7.26)
ggpubr::ggarrange(comp_plot1,comp_plot2,align = 'hv', common.legend = T, widths = c(54/28,1))
dev.off()










## ---------------  Plot for alpha diversity analysis ---------------
mbiome_subset <- microbiome_raw
metadf <- sample_data(mbiome_subset) %>% data.frame() %>%
  mutate(Treatment_Timing = ifelse(Treatment_Timing == 'Post_Treatment',
                                   'Post-treat','Pre-treat')) %>%
  # mutate(Status = ifelse(Status == 'Paired',
  #                        'Paired: Y','Paired: N')) %>%
  mutate(model_var = ifelse(Murine_Model == 'WT', 'WT', 'C3KO' )) %>%
  mutate(Treatment = ifelse(Treatment == 'CTLA4', '[CTLA-4]', '[IgG]' )) %>%
  mutate(Treatment2 = replace(Treatment, Treatment_Timing == 'Pre-treat', '')) %>%
  mutate(Marker  = paste(Treatment_Timing ,'\n',
                         model_var,'\n',Treatment2,'', sep = ''))  %>%
  mutate(Marker = factor(Marker,
                         levels = unique(Marker) %>% sort(decreasing = T)
  )) %>%
  # arrange(Marker,Tax,-Abundance) %>%
  # arrange(Treatment_Timing,Status,Pair_Group) %>%
  # mutate(Sample = factor(Sample, levels = unique(Sample) ))
  select(c('Treatment_Timing','Marker',"richness"))

# facet = pre/post

## Box plot across points:
stat.test <- metadf %>% #filter(Treatment_Timing == 'Pre-treat') %>%
  t_test( formula = richness ~ Marker, data = .) %>%
  add_significance() %>%
  add_xy_position(step.increase = 0.1) %>%
  mutate(p = paste('P-value = ',p))

selrows <- c(with(stat.test, intersect(grep('Pre',group1),grep('Pre',group2) )),
             with(stat.test, intersect(grep('Post',group1),grep('Post',group2) )))

stat.test <- stat.test[selrows,]

stat.test$y.position <- max(metadf$richness) + 4*seq.int(0,nrow(stat.test)-1)

# pe_col <- c('deeppink','darkolivegreen2')
# po_col_w <- c('darkolivegreen','black')
# po_col_c <- c('darkmagenta','navy')
gecol <- c(pe_col[2:1],po_col_w[2:1],po_col_c[2:1])
plot_view1 <- ggplot(metadf, aes(x = Marker, y = richness)) +
  geom_boxplot(aes(fill = Marker), alpha = 1, notch = F, outlier.size = 0) +
  # geom_violin(draw_quantiles = c(0.25, 0.5, 0.75)) +
  # geom_line(aes(color = SubjectID, group = SubjectID),size = 1) +
  # geom_jitter(aes(color = Treatment),width = 0.15, size = 3) +
  geom_jitter(width = 0.15, size = 2) +
  xlab(NULL) +
  stat_pvalue_manual(stat.test, label = "p", tip.length = 0.01, size= 4 ) +
  ylab('Inverse simpson') + #ylim(c(0, 30)) +
  ggtitle('Diversity analysis') +
  # guides(color = guide_legend(override.aes = list(size = 3) ) ) +
  theme_classic() + scale_fill_manual(values = gecol)+
  # scale_color_manual(values = BaRU)+
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
postscript('./plots/manuscript/mainFig_alpha_a.eps',  width = 9.0, height = 9)
plot_view1
dev.off()





## --------------- Beta diversity analysis and dispersion test  ---------
## We prepare the data table for plotting
## Pre sample



dt_comp <- subset_samples(mbiome_subset, Treatment_Timing == 'Pre_Treatment')
# dt_comp <- microbiome::transform(dt_comp, "compositional")
dt_comp_c <- ordinate(dt_comp , "PCoA", "bray", weighted=T)
ax_name <- paste(c('PCo1', 'PCo2'), ' [',
                 signif(dt_comp_c$values$Relative_eig * 100,2)[1:2],'%]',
                 sep='')

metadf <- sample_data(dt_comp) %>% data.frame() %>%
  mutate(model_var = ifelse(Murine_Model == 'WT', 'WT', 'C3KO' )) %>%
  mutate(Marker  = paste(model_var,'', sep = ''))


dist_bray <- vegdist(otu_table(dt_comp)@.Data %>% t(), method = "bray")
bd <- betadisper(dist_bray, metadf$Marker, type = 'centroid')
bd

set.seed(123)
permAnovaTest <- adonis(
  otu_table(dt_comp)@.Data %>% t() ~  Marker, data = metadf,
  permutations = 999,
  method = "bray"
)
permAnovaTest

pldata <- dt_comp_c$vectors %>% data.frame() %>% cbind(.,metadf) %>%
  cbind(.,bd$centroids[.$Marker,1:2])
# pldata$Name <- substring(rownames(pldata),20)



plot_view1 <- ggplot(pldata , aes(x =  Axis.1, y = Axis.2, color = Marker)) +
  geom_point(size = 5) +
  # ggrepel::geom_text_repel() +,label = Name
  # geom_text(hjust=0, vjust=0) +
  theme_classic() +
  geom_segment(aes(x = PCoA1, xend=Axis.1, y=PCoA2, yend=Axis.2, color = Marker),
               arrow = arrow(length = unit(0.25, "cm")),lwd=0.3) +
  scale_color_manual(NULL, values =  pe_col) +
  ggtitle("Pre-treatment") +
  xlab(ax_name[1]) + ylab(ax_name[2]) +
  annotate("text", x=-0.2, y=-0.4, size =5,
           label= paste("P-value =",signif(permAnovaTest$aov.tab$`Pr(>F)`[1],1))) +
  theme(legend.position = 'top',
        plot.title = element_text(size = 15, hjust = 0.5, color = 'black'),
        axis.text = element_text(size = 13, color = 'black'),
        axis.title = element_text(size = 13, color = 'black'),
        strip.text = element_text(size = 15, color = 'black'),
        legend.text = element_text(size = 11, color = 'black'),
        legend.background = element_rect(fill="grey99",
                                         size=0.5, linetype="solid"),
        legend.title =  element_blank() )

plot_view1



## Table for post comparison
dt_comp <- subset_samples(mbiome_subset, Treatment_Timing == 'Post_Treatment')
metadf <- sample_data(dt_comp) %>% data.frame() %>%
  mutate(model_var = ifelse(Murine_Model == 'WT', 'WT', 'C3KO' )) %>%
  mutate(Treatment = ifelse(Treatment == 'CTLA4', '[CTLA-4]', '[IgG]' )) %>%
  mutate(Marker  = paste(model_var,' ',Treatment, sep = ''))
gen_comb <- combn(unique(metadf$Marker),2) %>% rbind(NA)
for (i in 1:ncol(gen_comb)) {
  ind <- metadf$Marker %in% gen_comb[1:2,i]
  temp <- otu_table(dt_comp)@.Data %>% t()
  set.seed(123)
  permAnovaTest <- adonis(
    temp[ind,] ~  Marker, data = metadf[ind,],
    permutations = 999,
    method = "bray"
  )
  gen_comb[3,i] <- permAnovaTest$aov.tab$`Pr(>F)`[1]
}
gen_comb %<>% t()
colnames(gen_comb) <- c('Compare:1',"Compare2",'P-value')
write.csv(gen_comb,'./plots/manuscript/main_postcomp.csv')





## Post sample
dt_comp <- subset_samples(mbiome_subset, Treatment_Timing == 'Post_Treatment')
# dt_comp <- microbiome::transform(dt_comp, "compositional")
dt_comp_c <- ordinate(dt_comp , "PCoA", "bray", weighted=T)
ax_name <- paste(c('PCo1', 'PCo2'), ' [',
                 signif(dt_comp_c$values$Relative_eig * 100,2)[1:2],'%]',
                 sep='')

metadf <- sample_data(dt_comp) %>% data.frame() %>%
  mutate(model_var = ifelse(Murine_Model == 'WT', 'WT', 'C3KO' )) %>%
  mutate(Treatment = ifelse(Treatment == 'CTLA4', '[CTLA-4]', '[IgG]' )) %>%
  mutate(Marker  = paste(model_var,' ',Treatment, sep = ''))


dist_bray <- vegdist(otu_table(dt_comp)@.Data %>% t(), method = "bray")
bd <- betadisper(dist_bray, metadf$Marker, type = 'centroid')
bd

set.seed(123)
permAnovaTest <- adonis(
  otu_table(dt_comp)@.Data %>% t() ~  Marker, data = metadf,
  permutations = 999,
  method = "bray"
)
permAnovaTest

pldata <- dt_comp_c$vectors %>% data.frame() %>% cbind(.,metadf) %>%
  cbind(.,bd$centroids[.$Marker,1:2])
# pldata$Name <- substring(rownames(pldata),20)


# pe_col <- c('deeppink','darkolivegreen2')
# po_col_w <- c('darkolivegreen','black')
# po_col_c <- c('darkmagenta','navy')
gecol <- c(po_col_w[1],po_col_c,po_col_w[2])
gecol <- c(po_col_c,po_col_w)


plot_view2 <- ggplot(pldata , aes(x =  Axis.1, y = Axis.2, color = Marker)) +
  geom_point(size = 5) +
  # ggrepel::geom_text_repel() +,label = Name
  # geom_text(hjust=0, vjust=0) +
  theme_classic() +
  geom_segment(aes(x = PCoA1, xend=Axis.1, y=PCoA2, yend=Axis.2, color = Marker),
               arrow = arrow(length = unit(0.25, "cm")),lwd=0.3) +
  scale_color_manual("Treatment", values =  gecol) +
  ggtitle("Post-treatment") +
  xlab(ax_name[1]) + ylab(ax_name[2]) +
  # annotate("text", x=-0.2, y=-0.4, size =5,
  #          label= paste("P-value =",signif(permAnovaTest$aov.tab$`Pr(>F)`[1],1))) +
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


##


## Post sample: WT



dt_comp <- subset_samples(mbiome_subset, (Treatment_Timing == 'Post_Treatment') &
                            (Murine_Model == 'WT'))
# dt_comp <- microbiome::transform(dt_comp, "compositional")
dt_comp_c <- ordinate(dt_comp , "PCoA", "bray", weighted=T)
ax_name <- paste(c('PCo1', 'PCo2'), ' [',
                 signif(dt_comp_c$values$Relative_eig * 100,2)[1:2],'%]',
                 sep='')

metadf <- sample_data(dt_comp) %>% data.frame() %>%
  mutate(model_var = ifelse(Murine_Model == 'WT', 'WT', 'C3KO' )) %>%
  mutate(Treatment = ifelse(Treatment == 'CTLA4', 'CTLA-4', 'IgG' )) %>%
  mutate(Marker  = paste(Treatment, sep = ''))


dist_bray <- vegdist(otu_table(dt_comp)@.Data %>% t(), method = "bray")
bd <- betadisper(dist_bray, metadf$Marker, type = 'centroid')
bd

set.seed(123)
permAnovaTest <- adonis(
  otu_table(dt_comp)@.Data %>% t() ~  Marker, data = metadf,
  permutations = 999,
  method = "bray"
)
permAnovaTest

pldata <- dt_comp_c$vectors %>% data.frame() %>% cbind(.,metadf) %>%
  cbind(.,bd$centroids[.$Marker,1:2])
# pldata$Name <- substring(rownames(pldata),20)




plot_view3 <- ggplot(pldata , aes(x =  Axis.1, y = Axis.2, color = Marker)) +
  geom_point(size = 5) +
  # ggrepel::geom_text_repel() +,label = Name
  # geom_text(hjust=0, vjust=0) +
  theme_classic() +
  geom_segment(aes(x = PCoA1, xend=Axis.1, y=PCoA2, yend=Axis.2, color = Marker),
               arrow = arrow(length = unit(0.25, "cm")),lwd=0.3) +
  scale_color_manual("Treatment", values =  po_col_w) +
  ggtitle("Post-treatment: WT") +
  xlab(ax_name[1]) + ylab(ax_name[2]) +
  annotate("text", x=-0.0, y=-0.4, size =5,
           label= paste("P-value =",signif(permAnovaTest$aov.tab$`Pr(>F)`[1],1))) +
  theme(legend.position = 'top',
        plot.title = element_text(size = 15, hjust = 0.5, color = 'black'),
        axis.text = element_text(size = 13, color = 'black'),
        axis.title = element_text(size = 13, color = 'black'),
        strip.text = element_text(size = 15, color = 'black'),
        legend.text = element_text(size = 11, color = 'black'),
        legend.background = element_rect(fill="grey99",
                                         size=0.5, linetype="solid"),
        legend.title =  element_blank() )

plot_view3




## Post sample: C3KO



dt_comp <- subset_samples(mbiome_subset, (Treatment_Timing == 'Post_Treatment') &
                            (Murine_Model == 'C3KO'))
# dt_comp <- microbiome::transform(dt_comp, "compositional")
dt_comp_c <- ordinate(dt_comp , "PCoA", "bray", weighted=T)
ax_name <- paste(c('PCo1', 'PCo2'), ' [',
                 signif(dt_comp_c$values$Relative_eig * 100,2)[1:2],'%]',
                 sep='')

metadf <- sample_data(dt_comp) %>% data.frame() %>%
  mutate(model_var = ifelse(Murine_Model == 'WT', 'WT', 'C3KO' )) %>%
  mutate(Treatment = ifelse(Treatment == 'CTLA4', 'CTLA-4', 'IgG' )) %>%
  mutate(Marker  = paste(Treatment, sep = ''))


dist_bray <- vegdist(otu_table(dt_comp)@.Data %>% t(), method = "bray")
bd <- betadisper(dist_bray, metadf$Marker, type = 'centroid')
bd

set.seed(123)
permAnovaTest <- adonis(
  otu_table(dt_comp)@.Data %>% t() ~  Marker, data = metadf,
  permutations = 999,
  method = "bray"
)
permAnovaTest

pldata <- dt_comp_c$vectors %>% data.frame() %>% cbind(.,metadf) %>%
  cbind(.,bd$centroids[.$Marker,1:2])
# pldata$Name <- substring(rownames(pldata),20)




plot_view4 <- ggplot(pldata , aes(x =  Axis.1, y = Axis.2, color = Marker)) +
  geom_point(size = 5) +
  # ggrepel::geom_text_repel() +,label = Name
  # geom_text(hjust=0, vjust=0) +
  theme_classic() +
  geom_segment(aes(x = PCoA1, xend=Axis.1, y=PCoA2, yend=Axis.2, color = Marker),
               arrow = arrow(length = unit(0.25, "cm")),lwd=0.3) +
  scale_color_manual("Treatment", values =  po_col_c) +
  ggtitle("Post-treatment: C3KO") +
  xlab(ax_name[1]) + ylab(ax_name[2]) +
  annotate("text", x=-0.0, y=-0.4, size =5,
           label= paste("P-value =",signif(permAnovaTest$aov.tab$`Pr(>F)`[1],1))) +
  theme(legend.position = 'top',
        plot.title = element_text(size = 15, hjust = 0.5, color = 'black'),
        axis.text = element_text(size = 13, color = 'black'),
        axis.title = element_text(size = 13, color = 'black'),
        strip.text = element_text(size = 15, color = 'black'),
        legend.text = element_text(size = 11, color = 'black'),
        legend.background = element_rect(fill="grey99",
                                         size=0.5, linetype="solid"),
        legend.title =  element_blank() )

plot_view4



setEPS()
postscript('./plots/manuscript/mainFig_beta_b1.eps',  width = 7, height = 6)
plot_view1
dev.off()

setEPS()
postscript('./plots/manuscript/mainFig_beta_b2.eps', width = 7, height = 6)
plot_view2
dev.off()

setEPS()
postscript('./plots/manuscript/supFig_beta_a1.eps',  width = 7, height = 6)
plot_view3
dev.off()

setEPS()
postscript('./plots/manuscript/supFig_beta_a2.eps',  width = 7, height = 6)
plot_view4
dev.off()









## ---------------  Deseq2 plots (Family) --------------- ---------------------

## Pre-treatment comparison
mbiome_subset <- microbiome_raw
an_data <- subset_samples(mbiome_subset, Treatment_Timing == 'Pre_Treatment')
an_data <- microbiome::aggregate_taxa(an_data, 'Family')
an_data <- subset_taxa(an_data, rowSums(otu_table(an_data) > 0) > 3 )
diagdds = phyloseq_to_deseq2(an_data, ~ Murine_Model)

diagdds = DESeq(diagdds, test="LRT", fitType="parametric",
                sfType = c("ratio", "poscounts", "iterate")[3],
                minReplicatesForReplace=Inf, reduced = ~1)
# diagdds = DESeq(diagdds, test="Wald", fitType="parametric",minReplicatesForReplace=Inf)
res = results(diagdds, cooksCutoff = FALSE)
alpha =  0.05
# sigtab = res[which((res$padj < alpha) & (res$baseMean > 150) & (abs(res$log2FoldChange) > 10)), ]
sigtab = res[which((res$pvalue < alpha)), ]
dim(sigtab)
sigtab = cbind(as(sigtab, "data.frame"), as(tax_table(an_data)[rownames(sigtab), ], "matrix"))
head(sigtab)

# addWorksheet(wb,'Figure de1')
# writeData(wb, sheet='Figure de1', sigtab)

lab_names <- sigtab$unique
tind <- which(lab_names == 'Clostridia.vadinBB60.group_unknown')
lab_names[tind] <- 'Clostridia.vadinBB60'
tind <- which(lab_names == 'RF39_unknown')
lab_names[tind] <- 'Oscillospiraceae'


sigtab$color<-ifelse(sigtab$log2FoldChange >0, "Increase",  "Decrease")
bp_data <-data.frame(gen=lab_names,value=sigtab$log2FoldChange,
                     col=sigtab$color, index = 1:nrow(sigtab),
                     ID = rownames(sigtab))
lvl <- bp_data$gen[order(bp_data$value)]
a <- ggplot(data=bp_data ,
            aes(x=reorder(index, value),y=value,fill=col))+
  geom_bar(stat="identity")+ theme_classic()+
  xlab("Taxonomic identifier \n (Family)")+
  ylab("log2(fold change)") + coord_flip()+
  scale_x_discrete(labels= lvl) +
  scale_fill_manual(values = pe_col) +
  ggtitle('Differentially abundant taxa: Pre-treatment \n (WT vs C3KO)') +
  # annotate("text", y= -2, x=7, size = 5, angle = 90,
  #          label= paste("P-value <",0.05 ) ) +
  # annotate("text", x=8, y= 25, size = 7, angle = 270,
  #          label= "Nivo + CT" ) +
  theme(legend.position =  'none',#c(0.1,0.87),
        legend.background = element_rect(fill="transparent"),
        plot.title = element_text(size = 15, hjust = 0.5),
        plot.subtitle = element_text(size = 15, hjust = 0.5),
        axis.text = element_text(size = 13,color = 'black'),
        axis.title = element_text(size = 13),
        axis.line = element_line(colour = 'black', size = 0.2),
        strip.text = element_text(size = 13),
        legend.text = element_text(size = 13),
        legend.title = element_text(size = 13) )
a



dt_comp <- an_data# microbiome::transform(neostar_genus, "compositional")
dt_comp <- subset_taxa(dt_comp, taxa_names(dt_comp) %in% rownames(sigtab))


plmat<- cbind(as(counts(diagdds, normalized=T)[rownames(sigtab), ], "matrix"),
              as(tax_table(dt_comp)[rownames(sigtab), ], "matrix")) %>%
  as.data.frame()

lab_names <- plmat$unique
tind <- which(lab_names == 'Clostridia.vadinBB60.group_unknown')
lab_names[tind] <- 'Clostridia.vadinBB60'
tind <- which(lab_names == 'RF39_unknown')
lab_names[tind] <- 'Oscillospiraceae'

plmat$unique <- lab_names
names(lab_names) <- plmat$ID <- rownames(plmat)

plmat <- plmat[,c(1:27,c(34))] %>% gather(  key = "key",
                                            value = "value", -ID)

plmat <- dplyr::left_join(plmat, sample_data(dt_comp)[,c('Murine_Model','SampleID')],
                          by = c("key" = "SampleID")) %>%
  # mutate(MPR = gsub('Yes ','',paste(MPR, 'MPR'))) %>%
  mutate(value = as.numeric(value))

names(bp_data)[names(bp_data) == 'value'] <- 'lfc_value'
plmat <- dplyr::left_join(plmat, bp_data[,c('lfc_value','ID')],
                          by = c("ID" = "ID"))

plmat$ID <- factor(plmat$ID, levels = bp_data$ID[order(bp_data$lfc_value)])


kk <- ggplot(plmat %>% filter(value >1.5), aes(x = ID, y = value, fill = Murine_Model)) +
  geom_boxplot(notch = F,lwd=0.3) +
  scale_y_continuous(trans='log10') +
  scale_fill_manual(name = NULL, values = pe_col) +
  theme_classic() +
  ylab( 'Observed abundance (log scale)') +
  xlab('') +
  coord_flip() +
  scale_x_discrete(labels = lab_names)+
  theme(legend.position = 'top',
        plot.title = element_text(size = 15, hjust = 0.5, color = 'black'),
        axis.text.x = element_text(size = 13, hjust = 1.0, vjust = 0.5,
                                   color = 'black'),
        axis.line = element_line(colour = 'black', size = 0.2),
        # axis.text.y = element_text(size = 13, color = 'black'),
        axis.text.y = element_blank(),
        axis.title = element_text(size = 13, color = 'black',),
        strip.text = element_text(size = 15, color = 'black'),
        legend.text = element_text(size = 13, color = 'black'),
        legend.title = element_blank() )
kk





setEPS()
postscript('./plots/manuscript/supFig_deseq_c1.eps',  width = 18.7, height = 7.26)
ggpubr::ggarrange(a, kk, nrow = 1,
                  widths = c(1., 1.),
                  align = c("none", "h", "v", "hv")[2],
                  common.legend = F)
dev.off()




## Post-treatment C3KO: CTLA4 IgG
mbiome_subset <- microbiome_raw
an_data <- subset_samples(mbiome_subset, (Treatment_Timing == 'Post_Treatment')&
                            (Murine_Model == 'C3KO'))
an_data <- microbiome::aggregate_taxa(an_data, 'Family')
an_data <- subset_taxa(an_data, rowSums(otu_table(an_data) > 0) > 3 )
diagdds = phyloseq_to_deseq2(an_data, ~ Treatment)

diagdds = DESeq(diagdds, test="LRT", fitType="parametric",
                sfType = c("ratio", "poscounts", "iterate")[3],
                minReplicatesForReplace=Inf, reduced = ~1)
# diagdds = DESeq(diagdds, test="Wald", fitType="parametric",minReplicatesForReplace=Inf)
res = results(diagdds, cooksCutoff = FALSE)
alpha =  0.05
# sigtab = res[which((res$padj < alpha) & (res$baseMean > 150) & (abs(res$log2FoldChange) > 10)), ]
sigtab = res[which((res$pvalue < alpha)), ]
dim(sigtab)
sigtab = cbind(as(sigtab, "data.frame"), as(tax_table(an_data)[rownames(sigtab), ], "matrix"))
head(sigtab)

# addWorksheet(wb,'Figure de3')
# writeData(wb, sheet='Figure de3', sigtab)

lab_names <- sigtab$unique
# tind <- which(lab_names == 'Clostridia.vadinBB60.group_unknown')
# lab_names[tind] <- 'Clostridia.vadinBB60'



sigtab$color<-ifelse(sigtab$log2FoldChange >0, "Increase",  "Decrease")
bp_data <-data.frame(gen=lab_names,value=sigtab$log2FoldChange,
                     col=sigtab$color, index = 1:nrow(sigtab),
                     ID = rownames(sigtab))
lvl <- bp_data$gen[order(bp_data$value)]
a <- ggplot(data=bp_data ,
            aes(x=reorder(index, value),y=value,fill=col))+
  geom_bar(stat="identity")+ theme_classic()+
  xlab("Taxonomic identifier \n (Family)")+
  ylab("log2(fold change)") + coord_flip()+
  scale_x_discrete(labels= lvl) +
  scale_fill_manual(values = po_col_c) +
  ggtitle('Differentially abundant taxa: \n Post-treatment C3KO  (CTLA-4 vs IgG)') +
  # annotate("text", y= -2, x=7, size = 5, angle = 90,
  #          label= paste("P-value <",0.05 ) ) +
  # annotate("text", x=8, y= 25, size = 7, angle = 270,
  #          label= "Nivo + CT" ) +
  theme(legend.position =  'none',#c(0.1,0.87),
        legend.background = element_rect(fill="transparent"),
        plot.title = element_text(size = 15, hjust = 0.5),
        plot.subtitle = element_text(size = 15, hjust = 0.5),
        axis.text = element_text(size = 13,color = 'black'),
        axis.title = element_text(size = 13),
        axis.line = element_line(colour = 'black', size = 0.2),
        strip.text = element_text(size = 13),
        legend.text = element_text(size = 13),
        legend.title = element_text(size = 13) )
a



dt_comp <- an_data# microbiome::transform(neostar_genus, "compositional")
dt_comp <- subset_taxa(dt_comp, taxa_names(dt_comp) %in% rownames(sigtab))


plmat<- cbind(as(counts(diagdds, normalized=T)[rownames(sigtab), ], "matrix"),
              as(tax_table(dt_comp)[rownames(sigtab), ], "matrix")) %>%
  as.data.frame()

lab_names <- plmat$unique
tind <- which(lab_names == 'Clostridia.vadinBB60.group_unknown')
lab_names[tind] <- 'Clostridia.vadinBB60'

plmat$unique <- lab_names
names(lab_names) <- plmat$ID <- rownames(plmat)

plmat <- plmat[,c(1:28,c(35))] %>% gather(  key = "key",
                                            value = "value", -ID)

plmat <- dplyr::left_join(plmat, sample_data(dt_comp)[,c('Treatment','SampleID')],
                          by = c("key" = "SampleID")) %>%
  # mutate(MPR = gsub('Yes ','',paste(MPR, 'MPR'))) %>%
  mutate(value = as.numeric(value))

names(bp_data)[names(bp_data) == 'value'] <- 'lfc_value'
plmat <- dplyr::left_join(plmat, bp_data[,c('lfc_value','ID')],
                          by = c("ID" = "ID"))

plmat$ID <- factor(plmat$ID, levels = bp_data$ID[order(bp_data$lfc_value)])


kk <- ggplot(plmat %>% filter(value >1.5), aes(x = ID, y = value, fill = Treatment)) +
  geom_boxplot(notch = F,lwd=0.3) +
  scale_y_continuous(trans='log10') +
  scale_fill_manual(name = NULL, values = po_col_c) +
  theme_classic() +
  ylab( 'Observed abundance (log scale)') +
  xlab('') +
  coord_flip() +
  scale_x_discrete(labels = lab_names)+
  theme(legend.position = 'top',
        plot.title = element_text(size = 15, hjust = 0.5, color = 'black'),
        axis.text.x = element_text(size = 13, hjust = 1.0, vjust = 0.5,
                                   color = 'black'),
        axis.line = element_line(colour = 'black', size = 0.2),
        # axis.text.y = element_text(size = 13, color = 'black'),
        axis.text.y = element_blank(),
        axis.title = element_text(size = 13, color = 'black',),
        strip.text = element_text(size = 15, color = 'black'),
        legend.text = element_text(size = 13, color = 'black'),
        legend.title = element_blank() )
kk





setEPS()
postscript('./plots/manuscript/subFig_deseq_c4.eps',  width = 18.7, height = 5.26)
ggpubr::ggarrange(a, kk, nrow = 1,
                  widths = c(1., 1.),
                  align = c("none", "h", "v", "hv")[2],
                  common.legend = F)
dev.off()









## Post-treatment comparison: treatment IgG
mbiome_subset <- microbiome_raw
an_data <- subset_samples(mbiome_subset, (Treatment_Timing == 'Post_Treatment')&
                            (Treatment == 'IgG'))
an_data <- microbiome::aggregate_taxa(an_data, 'Family')
an_data <- subset_taxa(an_data, rowSums(otu_table(an_data) > 0) > 3 )
diagdds = phyloseq_to_deseq2(an_data, ~ Murine_Model)

diagdds = DESeq(diagdds, test="LRT", fitType="parametric",
                sfType = c("ratio", "poscounts", "iterate")[3],
                minReplicatesForReplace=Inf, reduced = ~1)
# diagdds = DESeq(diagdds, test="Wald", fitType="parametric",minReplicatesForReplace=Inf)
res = results(diagdds, cooksCutoff = FALSE)
alpha =  0.05
# sigtab = res[which((res$padj < alpha) & (res$baseMean > 150) & (abs(res$log2FoldChange) > 10)), ]
sigtab = res[which((res$pvalue < alpha)), ]
dim(sigtab)
sigtab = cbind(as(sigtab, "data.frame"), as(tax_table(an_data)[rownames(sigtab), ], "matrix"))
head(sigtab)

# addWorksheet(wb,'Figure de4')
# writeData(wb, sheet='Figure de4', sigtab)

lab_names <- sigtab$unique
tind <- which(lab_names == '[Eubacterium].coprostanoligenes.group')
lab_names[tind] <- '[Eubacterium].coprostanoligenes'



sigtab$color<-ifelse(sigtab$log2FoldChange >0, "Increase",  "Decrease")
bp_data <-data.frame(gen=lab_names,value=sigtab$log2FoldChange,
                     col=sigtab$color, index = 1:nrow(sigtab),
                     ID = rownames(sigtab))
lvl <- bp_data$gen[order(bp_data$value)]
a <- ggplot(data=bp_data ,
            aes(x=reorder(index, value),y=value,fill=col))+
  geom_bar(stat="identity")+ theme_classic()+
  xlab("Taxonomic identifier \n (Family)")+
  ylab("log2(fold change)") + coord_flip()+
  scale_x_discrete(labels= lvl) +
  scale_fill_manual(values = c(po_col_w[2],po_col_c[2])[2:1]) +
  ggtitle('Differentially abundant taxa: Post-treatment IgG \n (WT vs C3KO)') +
  # annotate("text", y= -2, x=7, size = 5, angle = 90,
  #          label= paste("P-value <",0.05 ) ) +
  # annotate("text", x=8, y= 25, size = 7, angle = 270,
  #          label= "Nivo + CT" ) +
  theme(legend.position =  'none',#c(0.1,0.87),
        legend.background = element_rect(fill="transparent"),
        plot.title = element_text(size = 15, hjust = 0.5),
        plot.subtitle = element_text(size = 15, hjust = 0.5),
        axis.text = element_text(size = 13,color = 'black'),
        axis.title = element_text(size = 13),
        axis.line = element_line(colour = 'black', size = 0.2),
        strip.text = element_text(size = 13),
        legend.text = element_text(size = 13),
        legend.title = element_text(size = 13) )
a



dt_comp <- an_data# microbiome::transform(neostar_genus, "compositional")
dt_comp <- subset_taxa(dt_comp, taxa_names(dt_comp) %in% rownames(sigtab))


plmat<- cbind(as(counts(diagdds, normalized=T)[rownames(sigtab), ], "matrix"),
              as(tax_table(dt_comp)[rownames(sigtab), ], "matrix")) %>%
  as.data.frame()

lab_names <- plmat$unique
tind <- which(lab_names == '[Eubacterium].coprostanoligenes.group')
lab_names[tind] <- '[Eubacterium].coprostanoligenes'

plmat$unique <- lab_names
names(lab_names) <- plmat$ID <- rownames(plmat)

plmat <- plmat[,c(1:28,c(35))] %>% gather(  key = "key",
                                            value = "value", -ID)

plmat <- dplyr::left_join(plmat, sample_data(dt_comp)[,c('Murine_Model','SampleID')],
                          by = c("key" = "SampleID")) %>%
  # mutate(MPR = gsub('Yes ','',paste(MPR, 'MPR'))) %>%
  mutate(value = as.numeric(value))

names(bp_data)[names(bp_data) == 'value'] <- 'lfc_value'
plmat <- dplyr::left_join(plmat, bp_data[,c('lfc_value','ID')],
                          by = c("ID" = "ID"))

plmat$ID <- factor(plmat$ID, levels = bp_data$ID[order(bp_data$lfc_value)])


kk <- ggplot(plmat %>% filter(value >1.5), aes(x = ID, y = value, fill = Murine_Model)) +
  geom_boxplot(notch = F,lwd=0.3) +
  scale_y_continuous(trans='log10') +
  scale_fill_manual(name = NULL, values = c(po_col_w[2],po_col_c[2])[2:1] ) +
  theme_classic() +
  ylab( 'Observed abundance (log scale)') +
  xlab('') +
  coord_flip() +
  scale_x_discrete(labels = lab_names)+
  theme(legend.position = 'top',
        plot.title = element_text(size = 15, hjust = 0.5, color = 'black'),
        axis.text.x = element_text(size = 13, hjust = 1.0, vjust = 0.5,
                                   color = 'black'),
        axis.line = element_line(colour = 'black', size = 0.2),
        # axis.text.y = element_text(size = 13, color = 'black'),
        axis.text.y = element_blank(),
        axis.title = element_text(size = 13, color = 'black',),
        strip.text = element_text(size = 15, color = 'black'),
        legend.text = element_text(size = 13, color = 'black'),
        legend.title = element_blank() )
kk





setEPS()
postscript('./plots/manuscript/supFig_deseq_c2.eps',  width = 18.7, height = 8.26)
ggpubr::ggarrange(a, kk, nrow = 1,
                  widths = c(1., 1.),
                  align = c("none", "h", "v", "hv")[2],
                  common.legend = F)
dev.off()






##----------- Post-treatment comparison: treatment CTLA4-----------
mbiome_subset <- microbiome_raw
an_data <- subset_samples(mbiome_subset, (Treatment_Timing == 'Post_Treatment')&
                            (Model_Treatment %in% c('WT_CTLA4','C3KO_CTLA4')))
an_data <- microbiome::aggregate_taxa(an_data, 'Family')
an_data <- subset_taxa(an_data, rowSums(otu_table(an_data) > 0) > 3 )
diagdds = phyloseq_to_deseq2(an_data, ~ Model_Treatment)

diagdds = DESeq(diagdds, test="LRT", fitType="parametric",
                sfType = c("ratio", "poscounts", "iterate")[3],
                minReplicatesForReplace=Inf, reduced = ~1)
# diagdds = DESeq(diagdds, test="Wald", fitType="parametric",minReplicatesForReplace=Inf)
res = results(diagdds, cooksCutoff = FALSE)
alpha =  0.05
# sigtab = res[which((res$padj < alpha) & (res$baseMean > 150) & (abs(res$log2FoldChange) > 10)), ]
sigtab = res[which((res$pvalue < alpha)), ]
dim(sigtab)
sigtab = cbind(as(sigtab, "data.frame"), as(tax_table(an_data)[rownames(sigtab), ], "matrix"))
head(sigtab)

# addWorksheet(wb,'Figure de6')
# writeData(wb, sheet='Figure de6', sigtab)


lab_names <- sigtab$unique
tind <- which(lab_names == 'Clostridia.vadinBB60.group_unknown')
lab_names[tind] <- 'Clostridia.vadinBB60'



sigtab$color<-ifelse(sigtab$log2FoldChange >0, "Increase",  "Decrease")
bp_data <-data.frame(gen=lab_names,value=sigtab$log2FoldChange,
                     col=sigtab$color, index = 1:nrow(sigtab),
                     ID = rownames(sigtab))
lvl <- bp_data$gen[order(bp_data$value)]
a <- ggplot(data=bp_data ,
            aes(x=reorder(index, value),y=value,fill=col))+
  geom_bar(stat="identity")+ theme_classic()+
  xlab("Taxonomic identifier \n (Family)")+
  ylab("log2(fold change)") + coord_flip()+
  scale_x_discrete(labels= lvl) +
  scale_fill_manual(values = c(po_col_c[1],po_col_w[1])) +
  ggtitle('Differentially abundant taxa: Post-treatment \n (WT [CTLA-4] vs C3KO [CTLA-4])') +
  theme(legend.position =  'none',#c(0.1,0.87),
        legend.background = element_rect(fill="transparent"),
        plot.title = element_text(size = 15, hjust = 0.5),
        plot.subtitle = element_text(size = 15, hjust = 0.5),
        axis.text = element_text(size = 13,color = 'black'),
        axis.title = element_text(size = 13),
        axis.line = element_line(colour = 'black', size = 0.2),
        strip.text = element_text(size = 13),
        legend.text = element_text(size = 13),
        legend.title = element_text(size = 13) )
a



dt_comp <- an_data# microbiome::transform(neostar_genus, "compositional")
dt_comp <- subset_taxa(dt_comp, taxa_names(dt_comp) %in% rownames(sigtab))


plmat<- cbind(as(counts(diagdds, normalized=T)[rownames(sigtab), ], "matrix"),
              as(tax_table(dt_comp)[rownames(sigtab), ], "matrix")) %>%
  as.data.frame()

lab_names <- plmat$unique
tind <- which(lab_names == 'Clostridia.vadinBB60.group_unknown')
lab_names[tind] <- 'Clostridia.vadinBB60'

plmat$unique <- lab_names
names(lab_names) <- plmat$ID <- rownames(plmat)

plmat <- plmat[,c(1:27,c(34))] %>% gather(  key = "key",
                                            value = "value", -ID)

plmat <- dplyr::left_join(plmat, sample_data(dt_comp)[,c('Model_Treatment','SampleID')],
                          by = c("key" = "SampleID")) %>%
  # mutate(MPR = gsub('Yes ','',paste(MPR, 'MPR'))) %>%
  mutate(value = as.numeric(value))

names(bp_data)[names(bp_data) == 'value'] <- 'lfc_value'
plmat <- dplyr::left_join(plmat, bp_data[,c('lfc_value','ID')],
                          by = c("ID" = "ID"))

plmat$ID <- factor(plmat$ID, levels = bp_data$ID[order(bp_data$lfc_value)])


kk <- ggplot(plmat %>% filter(value >1.5), aes(x = ID, y = value, fill = Model_Treatment)) +
  geom_boxplot(notch = F,lwd=0.3) +
  scale_y_continuous(trans='log10') +
  scale_fill_manual(name = NULL, values = c(po_col_c[1],po_col_w[1])) +
  theme_classic() +
  ylab( 'Observed abundance (log scale)') +
  xlab('') +
  coord_flip() +
  scale_x_discrete(labels = lab_names)+
  theme(legend.position = 'top',
        plot.title = element_text(size = 15, hjust = 0.5, color = 'black'),
        axis.text.x = element_text(size = 13, hjust = 1.0, vjust = 0.5,
                                   color = 'black'),
        axis.line = element_line(colour = 'black', size = 0.2),
        # axis.text.y = element_text(size = 13, color = 'black'),
        axis.text.y = element_blank(),
        axis.title = element_text(size = 13, color = 'black',),
        strip.text = element_text(size = 15, color = 'black'),
        legend.text = element_text(size = 13, color = 'black'),
        legend.title = element_blank() )
kk



setEPS()
postscript('./plots/manuscript/supFig_deseq_c3.eps',  width = 18.7, height = 8)
ggpubr::ggarrange(a, kk, nrow = 1,
                  widths = c(1., 1.),
                  align = c("none", "h", "v", "hv")[2],
                  common.legend = F)
dev.off()













## Pre-treatment comparison
mbiome_subset <- microbiome_raw
an_data <- subset_samples(mbiome_subset, Treatment_Timing == 'Pre_Treatment')
an_data <- microbiome::aggregate_taxa(an_data, 'Genus')
an_data <- subset_taxa(an_data, rowSums(otu_table(an_data) > 0) > 3 )
diagdds = phyloseq_to_deseq2(an_data, ~ Murine_Model)

diagdds = DESeq(diagdds, test="LRT", fitType="parametric",
                sfType = c("ratio", "poscounts", "iterate")[3],
                minReplicatesForReplace=Inf, reduced = ~1)
# diagdds = DESeq(diagdds, test="Wald", fitType="parametric",minReplicatesForReplace=Inf)
res = results(diagdds, cooksCutoff = FALSE)
alpha =  0.05
# sigtab = res[which((res$padj < alpha) & (res$baseMean > 150) & (abs(res$log2FoldChange) > 10)), ]
sigtab = res[which((res$pvalue < alpha)), ]
dim(sigtab)
sigtab = cbind(as(sigtab, "data.frame"), as(tax_table(an_data)[rownames(sigtab), ], "matrix"))
head(sigtab)

# addWorksheet(wb,'Figure de1')
# writeData(wb, sheet='Figure de1', sigtab)

lab_names <- sigtab$unique
tind <- which(lab_names == 'Clostridia.vadinBB60.group_unknown')
lab_names[tind] <- 'Clostridia.vadinBB60'
tind <- which(lab_names == 'RF39_unknown')
lab_names[tind] <- 'Oscillospiraceae'


sigtab$color<-ifelse(sigtab$log2FoldChange >0, "Increase",  "Decrease")
bp_data <-data.frame(gen=lab_names,value=sigtab$log2FoldChange,
                     col=sigtab$color, index = 1:nrow(sigtab),
                     ID = rownames(sigtab))
lvl <- bp_data$gen[order(bp_data$value)]
a <- ggplot(data=bp_data ,
            aes(x=reorder(index, value),y=value,fill=col))+
  geom_bar(stat="identity")+ theme_classic()+
  xlab("Taxonomic identifier \n (Genus)")+
  ylab("log2(fold change)") + coord_flip()+
  scale_x_discrete(labels= lvl) +
  scale_fill_manual(values = pe_col) +
  ggtitle('Differentially abundant taxa: Pre-treatment \n (WT vs C3KO)') +
  # annotate("text", y= -2, x=7, size = 5, angle = 90,
  #          label= paste("P-value <",0.05 ) ) +
  # annotate("text", x=8, y= 25, size = 7, angle = 270,
  #          label= "Nivo + CT" ) +
  theme(legend.position =  'none',#c(0.1,0.87),
        legend.background = element_rect(fill="transparent"),
        plot.title = element_text(size = 15, hjust = 0.5),
        plot.subtitle = element_text(size = 15, hjust = 0.5),
        axis.text = element_text(size = 13,color = 'black'),
        axis.title = element_text(size = 13),
        axis.line = element_line(colour = 'black', size = 0.2),
        strip.text = element_text(size = 13),
        legend.text = element_text(size = 13),
        legend.title = element_text(size = 13) )
a



dt_comp <- an_data# microbiome::transform(neostar_genus, "compositional")
dt_comp <- subset_taxa(dt_comp, taxa_names(dt_comp) %in% rownames(sigtab))


plmat<- cbind(as(counts(diagdds, normalized=T)[rownames(sigtab), ], "matrix"),
              as(tax_table(dt_comp)[rownames(sigtab), ], "matrix")) %>%
  as.data.frame()

lab_names <- plmat$unique
tind <- which(lab_names == 'Clostridia.vadinBB60.group_unknown')
lab_names[tind] <- 'Clostridia.vadinBB60'
tind <- which(lab_names == 'RF39_unknown')
lab_names[tind] <- 'Oscillospiraceae'

plmat$unique <- lab_names
names(lab_names) <- plmat$ID <- rownames(plmat)

plmat <- plmat[,c(1:27,c(35))] %>% gather(  key = "key",
                                            value = "value", -ID)

plmat <- dplyr::left_join(plmat, sample_data(dt_comp)[,c('Murine_Model','SampleID')],
                          by = c("key" = "SampleID")) %>%
  # mutate(MPR = gsub('Yes ','',paste(MPR, 'MPR'))) %>%
  mutate(value = as.numeric(value))

names(bp_data)[names(bp_data) == 'value'] <- 'lfc_value'
plmat <- dplyr::left_join(plmat, bp_data[,c('lfc_value','ID')],
                          by = c("ID" = "ID"))

plmat$ID <- factor(plmat$ID, levels = bp_data$ID[order(bp_data$lfc_value)])


kk <- ggplot(plmat %>% filter(value >1.5), aes(x = ID, y = value, fill = Murine_Model)) +
  geom_boxplot(notch = F,lwd=0.3) +
  scale_y_continuous(trans='log10') +
  scale_fill_manual(name = NULL, values = pe_col) +
  theme_classic() +
  ylab( 'Observed abundance (log scale)') +
  xlab('') +
  coord_flip() +
  scale_x_discrete(labels = lab_names)+
  theme(legend.position = 'top',
        plot.title = element_text(size = 15, hjust = 0.5, color = 'black'),
        axis.text.x = element_text(size = 13, hjust = 1.0, vjust = 0.5,
                                   color = 'black'),
        axis.line = element_line(colour = 'black', size = 0.2),
        # axis.text.y = element_text(size = 13, color = 'black'),
        axis.text.y = element_blank(),
        axis.title = element_text(size = 13, color = 'black',),
        strip.text = element_text(size = 15, color = 'black'),
        legend.text = element_text(size = 13, color = 'black'),
        legend.title = element_blank() )
kk





setEPS()
postscript('./plots/manuscript/mainFig_deseq_c1.eps',  width = 18.7, height = 7.26)
ggpubr::ggarrange(a, kk, nrow = 1,
                  widths = c(1., 1.),
                  align = c("none", "h", "v", "hv")[2],
                  common.legend = F)
dev.off()










##----------- Post-treatment comparison: treatment CTLA4-----------
mbiome_subset <- microbiome_raw
an_data <- subset_samples(mbiome_subset, (Treatment_Timing == 'Post_Treatment')&
                            (Model_Treatment %in% c('WT_CTLA4','C3KO_CTLA4')))
an_data <- microbiome::aggregate_taxa(an_data, 'Genus')
an_data <- subset_taxa(an_data, rowSums(otu_table(an_data) > 0) > 3 )
diagdds = phyloseq_to_deseq2(an_data, ~ Model_Treatment)

diagdds = DESeq(diagdds, test="LRT", fitType="parametric",
                sfType = c("ratio", "poscounts", "iterate")[3],
                minReplicatesForReplace=Inf, reduced = ~1)
# diagdds = DESeq(diagdds, test="Wald", fitType="parametric",minReplicatesForReplace=Inf)
res = results(diagdds, cooksCutoff = FALSE)
alpha =  0.05
# sigtab = res[which((res$padj < alpha) & (res$baseMean > 150) & (abs(res$log2FoldChange) > 10)), ]
sigtab = res[which((res$pvalue < alpha)), ]
dim(sigtab)
sigtab = cbind(as(sigtab, "data.frame"), as(tax_table(an_data)[rownames(sigtab), ], "matrix"))
head(sigtab)

# addWorksheet(wb,'Figure de6')
# writeData(wb, sheet='Figure de6', sigtab)


lab_names <- sigtab$unique
tind <- which(lab_names == 'Clostridia.vadinBB60.group_unknown')
lab_names[tind] <- 'Clostridia.vadinBB60'



sigtab$color<-ifelse(sigtab$log2FoldChange >0, "Increase",  "Decrease")
bp_data <-data.frame(gen=lab_names,value=sigtab$log2FoldChange,
                     col=sigtab$color, index = 1:nrow(sigtab),
                     ID = rownames(sigtab))
lvl <- bp_data$gen[order(bp_data$value)]
a <- ggplot(data=bp_data ,
            aes(x=reorder(index, value),y=value,fill=col))+
  geom_bar(stat="identity")+ theme_classic()+
  xlab("Taxonomic identifier \n (Genus)")+
  ylab("log2(fold change)") + coord_flip()+
  scale_x_discrete(labels= lvl) +
  scale_fill_manual(values = c(po_col_c[1],po_col_w[1])) +
  ggtitle('Differentially abundant taxa: Post-treatment \n (WT [CTLA-4] vs C3KO [CTLA-4])') +
  # annotate("text", y= -2, x=7, size = 5, angle = 90,
  #          label= paste("P-value <",0.05 ) ) +
  # annotate("text", x=8, y= 25, size = 7, angle = 270,
  #          label= "Nivo + CT" ) +
  theme(legend.position =  'none',#c(0.1,0.87),
        legend.background = element_rect(fill="transparent"),
        plot.title = element_text(size = 15, hjust = 0.5),
        plot.subtitle = element_text(size = 15, hjust = 0.5),
        axis.text = element_text(size = 13,color = 'black'),
        axis.title = element_text(size = 13),
        axis.line = element_line(colour = 'black', size = 0.2),
        strip.text = element_text(size = 13),
        legend.text = element_text(size = 13),
        legend.title = element_text(size = 13) )
a



dt_comp <- an_data# microbiome::transform(neostar_genus, "compositional")
dt_comp <- subset_taxa(dt_comp, taxa_names(dt_comp) %in% rownames(sigtab))


plmat<- cbind(as(counts(diagdds, normalized=T)[rownames(sigtab), ], "matrix"),
              as(tax_table(dt_comp)[rownames(sigtab), ], "matrix")) %>%
  as.data.frame()

lab_names <- plmat$unique
tind <- which(lab_names == 'Clostridia.vadinBB60.group_unknown')
lab_names[tind] <- 'Clostridia.vadinBB60'

plmat$unique <- lab_names
names(lab_names) <- plmat$ID <- rownames(plmat)

plmat <- plmat[,c(1:27,c(35))] %>% gather(  key = "key",
                                            value = "value", -ID)

plmat <- dplyr::left_join(plmat, sample_data(dt_comp)[,c('Model_Treatment','SampleID')],
                          by = c("key" = "SampleID")) %>%
  # mutate(MPR = gsub('Yes ','',paste(MPR, 'MPR'))) %>%
  mutate(value = as.numeric(value))

names(bp_data)[names(bp_data) == 'value'] <- 'lfc_value'
plmat <- dplyr::left_join(plmat, bp_data[,c('lfc_value','ID')],
                          by = c("ID" = "ID"))

plmat$ID <- factor(plmat$ID, levels = bp_data$ID[order(bp_data$lfc_value)])


kk <- ggplot(plmat %>% filter(value >1.5), aes(x = ID, y = value, fill = Model_Treatment)) +
  geom_boxplot(notch = F,lwd=0.3) +
  scale_y_continuous(trans='log10') +
  scale_fill_manual(name = NULL, values = c(po_col_c[1],po_col_w[1])) +
  theme_classic() +
  ylab( 'Observed abundance (log scale)') +
  xlab('') +
  coord_flip() +
  scale_x_discrete(labels = lab_names)+
  theme(legend.position = 'top',
        plot.title = element_text(size = 15, hjust = 0.5, color = 'black'),
        axis.text.x = element_text(size = 13, hjust = 1.0, vjust = 0.5,
                                   color = 'black'),
        axis.line = element_line(colour = 'black', size = 0.2),
        # axis.text.y = element_text(size = 13, color = 'black'),
        axis.text.y = element_blank(),
        axis.title = element_text(size = 13, color = 'black',),
        strip.text = element_text(size = 15, color = 'black'),
        legend.text = element_text(size = 13, color = 'black'),
        legend.title = element_blank() )
kk



setEPS()
postscript('./plots/manuscript/mainFig_deseq_c2.eps',  width = 18.7, height = 8)
ggpubr::ggarrange(a, kk, nrow = 1,
                  widths = c(1., 1.),
                  align = c("none", "h", "v", "hv")[2],
                  common.legend = F)
dev.off()













## Post-treatment C3KO: CTLA4 IgG
mbiome_subset <- microbiome_raw
an_data <- subset_samples(mbiome_subset, (Treatment_Timing == 'Post_Treatment')&
                            (Murine_Model == 'C3KO'))
an_data <- microbiome::aggregate_taxa(an_data, 'Genus')
an_data <- subset_taxa(an_data, rowSums(otu_table(an_data) > 0) > 3 )
diagdds = phyloseq_to_deseq2(an_data, ~ Treatment)

diagdds = DESeq(diagdds, test="LRT", fitType="parametric",
                sfType = c("ratio", "poscounts", "iterate")[3],
                minReplicatesForReplace=Inf, reduced = ~1)
# diagdds = DESeq(diagdds, test="Wald", fitType="parametric",minReplicatesForReplace=Inf)
res = results(diagdds, cooksCutoff = FALSE)
alpha =  0.05
# sigtab = res[which((res$padj < alpha) & (res$baseMean > 150) & (abs(res$log2FoldChange) > 10)), ]
sigtab = res[which((res$pvalue < alpha)), ]
dim(sigtab)
sigtab = cbind(as(sigtab, "data.frame"), as(tax_table(an_data)[rownames(sigtab), ], "matrix"))
head(sigtab)

# addWorksheet(wb,'Figure de3')
# writeData(wb, sheet='Figure de3', sigtab)

lab_names <- sigtab$unique
# tind <- which(lab_names == 'Clostridia.vadinBB60.group_unknown')
# lab_names[tind] <- 'Clostridia.vadinBB60'



sigtab$color<-ifelse(sigtab$log2FoldChange >0, "Increase",  "Decrease")
bp_data <-data.frame(gen=lab_names,value=sigtab$log2FoldChange,
                     col=sigtab$color, index = 1:nrow(sigtab),
                     ID = rownames(sigtab))
lvl <- bp_data$gen[order(bp_data$value)]
a <- ggplot(data=bp_data ,
            aes(x=reorder(index, value),y=value,fill=col))+
  geom_bar(stat="identity")+ theme_classic()+
  xlab("Taxonomic identifier \n (Genus)")+
  ylab("log2(fold change)") + coord_flip()+
  scale_x_discrete(labels= lvl) +
  scale_fill_manual(values = po_col_c) +
  ggtitle('Differentially abundant taxa: \n Post-treatment C3KO  (CTLA-4 vs IgG)') +
  # annotate("text", y= -2, x=7, size = 5, angle = 90,
  #          label= paste("P-value <",0.05 ) ) +
  # annotate("text", x=8, y= 25, size = 7, angle = 270,
  #          label= "Nivo + CT" ) +
  theme(legend.position =  'none',#c(0.1,0.87),
        legend.background = element_rect(fill="transparent"),
        plot.title = element_text(size = 15, hjust = 0.5),
        plot.subtitle = element_text(size = 15, hjust = 0.5),
        axis.text = element_text(size = 13,color = 'black'),
        axis.title = element_text(size = 13),
        axis.line = element_line(colour = 'black', size = 0.2),
        strip.text = element_text(size = 13),
        legend.text = element_text(size = 13),
        legend.title = element_text(size = 13) )
a



dt_comp <- an_data# microbiome::transform(neostar_genus, "compositional")
dt_comp <- subset_taxa(dt_comp, taxa_names(dt_comp) %in% rownames(sigtab))


plmat<- cbind(as(counts(diagdds, normalized=T)[rownames(sigtab), ], "matrix"),
              as(tax_table(dt_comp)[rownames(sigtab), ], "matrix")) %>%
  as.data.frame()

lab_names <- plmat$unique
tind <- which(lab_names == 'Clostridia.vadinBB60.group_unknown')
lab_names[tind] <- 'Clostridia.vadinBB60'

plmat$unique <- lab_names
names(lab_names) <- plmat$ID <- rownames(plmat)

plmat <- plmat[,c(1:28,c(36))] %>% gather(  key = "key",
                                            value = "value", -ID)

plmat <- dplyr::left_join(plmat, sample_data(dt_comp)[,c('Treatment','SampleID')],
                          by = c("key" = "SampleID")) %>%
  # mutate(MPR = gsub('Yes ','',paste(MPR, 'MPR'))) %>%
  mutate(value = as.numeric(value))

names(bp_data)[names(bp_data) == 'value'] <- 'lfc_value'
plmat <- dplyr::left_join(plmat, bp_data[,c('lfc_value','ID')],
                          by = c("ID" = "ID"))

plmat$ID <- factor(plmat$ID, levels = bp_data$ID[order(bp_data$lfc_value)])


kk <- ggplot(plmat %>% filter(value >1.5), aes(x = ID, y = value, fill = Treatment)) +
  geom_boxplot(notch = F,lwd=0.3) +
  scale_y_continuous(trans='log10') +
  scale_fill_manual(name = NULL, values = po_col_c) +
  theme_classic() +
  ylab( 'Observed abundance (log scale)') +
  xlab('') +
  coord_flip() +
  scale_x_discrete(labels = lab_names)+
  theme(legend.position = 'top',
        plot.title = element_text(size = 15, hjust = 0.5, color = 'black'),
        axis.text.x = element_text(size = 13, hjust = 1.0, vjust = 0.5,
                                   color = 'black'),
        axis.line = element_line(colour = 'black', size = 0.2),
        # axis.text.y = element_text(size = 13, color = 'black'),
        axis.text.y = element_blank(),
        axis.title = element_text(size = 13, color = 'black',),
        strip.text = element_text(size = 15, color = 'black'),
        legend.text = element_text(size = 13, color = 'black'),
        legend.title = element_blank() )
kk



setEPS()
postscript('./plots/manuscript/mainFig_deseq_c3.eps',  width = 18.7, height = 5.26)
ggpubr::ggarrange(a, kk, nrow = 1,
                  widths = c(1., 1.),
                  align = c("none", "h", "v", "hv")[2],
                  common.legend = F)
dev.off()
