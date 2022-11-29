---
  title: "16S.Gala.PH.2021"
Author: "Rose.KV."
output: html_notebook
---
  #Libraries needed







#Beta Diversity by Timepoints in control Skin

##Weighted Unifrac stats 
Analyzed by Timepoints
```{r}
## PERMANOVA
library(vegan)
physq_sub <- physeq_filtered2 %>% 
  subset_samples(Treatments != "None") %>%
  subset_samples(Treatments != "ASM-treated")  %>%
  subset_samples(Tissue_types != "Pulp")  %>%
  #subset_samples(Timepoints != "Day 49")  %>%
  subset_samples(Timepoints != "Day 7")  %>%
  subset_samples(Treatments == "Control")
# physeq_filtered3 <- physeq_filtered2 %>%
#   subset_samples(Treatments != "None") 

wuinfrac_dist <- phyloseq::distance(physq_sub, method="wunifrac") #RUN this only once because it takes a lot of time

adonis_wunifrac = adonis2(wuinfrac_dist ~ sample_data(physq_sub)$Timepoints)
adonis_wunifrac
## Significant PERMANOVA indicates that centroid (or spatial median) among groups is different and/or with-group dispersion among groups is different
## PERMDISP

# subset phyloseq object as for sample data
wuni_disp <-betadisper(wuinfrac_dist, sample_data(physq_sub)$Timepoints, type=c("median"))
anova(wuni_disp)

## If PERMANOVA and PERMDISP are both significant, you can use plotting to tell if PERMANOVA was significant based on centroid (or spatial median)
plot(wuni_disp)

#?plot.betadisper
## Would look better with higher replication for groups
plot(wuni_disp, label = F)

## Plot with 1 standard deviation ellipses around the group medians
## sample size issue here, but you get the idea
plot(wuni_disp, label = F, hull = F, ellipse = T)

## Within-group dispersion that PERMDISP is testing
boxplot(wuni_disp, las = 2, cex.lab=1.5)
?boxplot

## pairwise p-values
TukeyHSD(wuni_disp)
scores(wuni_disp, 1:4, display = "centroids")
rda(otu_table)
```

##Weighted unifrac beta diversity plot with phyloseq
Used to check the PCoA %
```{r}
beta_wu <- ordinate(physq_sub, "PCoA", "wunifrac")
beta_wu_plot = plot_ordination(physq_sub, beta_wu, type="Timepoints", color="Timepoints", shape="Timepoints", title="PCoA Weighted Unifrac") + stat_ellipse(type = "t", linetype = 3) + stat_ellipse(type = "t") + theme_bw()+labs(colour = "Timepoints") #To add arrows https://neavemj.github.io/posts/coralMicrobiome
beta_wu_plot
```


#Final weighted unifrac plot
```{r}
label_perm <- expression(paste("PERMANOVA, ",R^2 ,"= 0.45, ", paste(italic('p')),"=0.001"))
beta_scatter = as.data.frame(beta_wu[["vectors"]])
beta_meta = merge(beta_scatter,metadata,by = 0,all=F)
pmain_wuF = ggscatter(beta_meta, x = "Axis.1", y = "Axis.2", color = "Timepoints", palette = "aaas",ellipse = TRUE, ellipse.level=.5,mean.point = F, mean.point.size = 5, star.plot = F) +labs(x = "PCoA 1 (40.2%) ", y = "PCoA 2 (19.3%)", colour = "Timepoints", fill = "Timepoints") +annotate("text", x = -0.04, y = -0.06, label = label_perm, colour = "black")
ggtitle("Control: Skin")
pmain_wuF
```



#Beta Diversity by Timepoints in control pulp

##Weighted Unifrac stats 
Analyzed by Timepoints
```{r}
## PERMANOVA
library(vegan)
physq_sub <- physeq_filtered2 %>% 
  subset_samples(Treatments != "None") %>%
  subset_samples(Treatments != "ASM-treated")  %>%
  subset_samples(Tissue_types != "Skin")  %>%
  #subset_samples(Timepoints != "Day 14")  %>%
  # subset_samples(Timepoints != "Day 49")  %>%
  # subset_samples(Timepoints != "Day 25")  %>%
  subset_samples(Timepoints != "Day 7")  %>%
  subset_samples(Treatments == "Control")
# physeq_filtered3 <- physeq_filtered2 %>%
#   subset_samples(Treatments != "None") 

wuinfrac_dist <- phyloseq::distance(physq_sub, method="wunifrac") #RUN this only once because it takes a lot of time

adonis_wunifrac = adonis2(wuinfrac_dist ~ sample_data(physq_sub)$Timepoints)
adonis_wunifrac
## Significant PERMANOVA indicates that centroid (or spatial median) among groups is different and/or with-group dispersion among groups is different
## PERMDISP

# subset phyloseq object as for sample data
wuni_disp <-betadisper(wuinfrac_dist, sample_data(physq_sub)$Timepoints, type=c("median"))
anova(wuni_disp)

## If PERMANOVA and PERMDISP are both significant, you can use plotting to tell if PERMANOVA was significant based on centroid (or spatial median)
plot(wuni_disp)

#?plot.betadisper
## Would look better with higher replication for groups
plot(wuni_disp, label = F)

## Plot with 1 standard deviation ellipses around the group medians
## sample size issue here, but you get the idea
plot(wuni_disp, label = F, hull = F, ellipse = T)

## Within-group dispersion that PERMDISP is testing
boxplot(wuni_disp, las = 2, cex.lab=1.5)
?boxplot

## pairwise p-values
TukeyHSD(wuni_disp)
scores(wuni_disp, 1:4, display = "centroids")
rda(otu_table)