library(breakaway)
library(ggplot2)
library(forcats)
library(vegan)
library(ggpubr)
library(reshape2)
library(viridis)
library(RColorBrewer)
library(rstatix)

### chao1 for alpha diversity analysis at species level
data_count_species = dplyr::select(data_count, contains("s__"))

# start by creating the frequency count table
br_sp = t(data_count_species)
br_sp_freq = build_frequency_count_tables(br_sp)


# Then run chao1 estimation on each columns. Columns are samples, rows are species
vect_br <- c()
for (i in 1:length(colnames(br_sp))) {
  vect_br <- c(vect_br, chao1_bc(br_sp_freq[[colnames(br_sp)[i]]])$estimate )
}

#  data wrangling for graphic representation
br_sp_NDL = vect_br[1:41]
br_sp_MDG = vect_br[42:74]
br_sp_CND = vect_br[75:95]
br_sp_PRU = vect_br[96:100]
br_sp_USA = vect_br[101:132]
br_sp_SPN = vect_br[133:142]
br_sp_FJI = vect_br[143:177]
br_sp_NUN = vect_br[178:456]

br_sp_non_west = c((br_sp_MDG), (br_sp_FJI), (br_sp_PRU))
br_sp_west = c((br_sp_USA), (br_sp_CND), (br_sp_NDL), (br_sp_SPN))

Lifestyle <-c(rep("Nunavik", 279),
              rep("Non-industrialized",73), 
              rep("Industrialized", 104))
break_sp <- c((br_sp_NUN) , (br_sp_non_west), (br_sp_west))

br_data = data.frame(Lifestyle, break_sp)

br_data$Lifestyle <- factor(br_data$Lifestyle, 
                            levels = c("Nunavik", "Non-industrialized", "Industrialized"))

Lifestyle = br_data$Lifestyle

# Wilcoxon-Mann-Whitney test
stat.test <- compare_means(
  break_sp ~ Lifestyle, data = br_data,
  method = "wilcox.test")

ggplot(br_data, aes( y=break_sp, x=br_data$Lifestyle, color = Lifestyle)) + 
  theme_gray() +
  geom_boxplot()+
  theme(text = element_text(size=10), legend.position = "none")+
  xlab(paste("")) +
  ylab(paste("Chao1 richness estimate"))+
  scale_x_discrete(limits=c("Nunavik","Non-industrialized","Industrialized"))+
  stat_pvalue_manual(stat.test,  y.position = c( 180, 190,200),label = "p.signif")+
  ggtitle("") +
  scale_color_manual(values =  c( "#E69F00","#CC79A7", "#56B4E9" ))


# PCoA on Jaccard similarity index at species level
# start by scaling data
s = scale(species_M3_complete, scale = T, center = T)

# apply vegdist to compute similarity matrix on scaled data 
Jaccard_sp = vegdist(s, method = "jaccard", binary = T, diag = T,
                     upper = T, na.rm = T)   

#### pcoa on Jaccard similarity
pcoa_Jaccard_sp = cmdscale(Jaccard_sp, eig = T, k = 4, x.ret = T)

# variation calculation
mds.var.per.sp <- round(abs(pcoa_Jaccard_sp$eig)/sum(abs(pcoa_Jaccard_sp$eig))*100, 1)

# print to check that the calculation was done
mds.var.per.sp

# Format data for ggplot
mds.data.sp <- pcoa_Jaccard_sp$points
colnames(mds.data.sp) <- c("PCoA1", "PCoA2", "PCoA3", "PCoA4")
mds.data.sp = as.data.frame(mds.data.sp)
mds.data.sp$Age <- Metadata_Complete$Age
mds.data.sp$Sex <- Metadata_Complete$Sex
mds.data.sp$Lifestyle = Metadata_Complete$Lifestyle

X= mds.data.sp$PCoA1
Y= mds.data.sp$PCoA2

######

ggplot(mds.data.sp, aes(x=X, y=Y, color=Lifestyle)) +
  theme_gray() +
  geom_point() +
  xlab(paste("PCoA1 - ", mds.var.per.sp[1], "%", sep="")) +
  ylab(paste("PCoA2 - ", mds.var.per.sp[2], "%", sep="")) +
  ggtitle("PCoA on Jaccard distance species level") +
  theme(text = element_text(size=10))+
  scale_color_manual(breaks=c("Nunavik","Non-industrialized",
                              "Industrialized"), values=c("#E69F00","#CC79A7", "#56B4E9" ))+
  stat_conf_ellipse(
    position = "identity",
    na.rm = FALSE,
    inherit.aes = TRUE,
    level = 0.95,
    npoint = 100,
    bary = TRUE
  )

# PERMANOVA calculation
vectlifestyle <- Metadata_Complete$Lifestyle
vectage <- Metadata_Complete$Age
vectsex <- Metadata_Complete$Sex

Permanova_Jaccard_sp = adonis2(Jaccard_sp ~ vectlifestyle + vectage + vectsex , permutations = 999, method = "jaccard")
Permanova_Jaccard_sp

# Top 20 most important species graphic representation
# data wrangling 

Top_20sp_NDL = as.data.frame(Top_20sp[1:41,])
Top_20sp_MDG = as.data.frame(Top_20sp[42:74,])
Top_20sp_CND = as.data.frame(Top_20sp[75:95,])
Top_20sp_PRU = as.data.frame(Top_20sp[96:100,])
Top_20sp_USA = as.data.frame(Top_20sp[101:132,])
Top_20sp_SPN = as.data.frame(Top_20sp[133:142,])
Top_20sp_FJI = as.data.frame(Top_20sp[143:177,])
Top_20sp_Nun = as.data.frame(Top_20sp[178:456,])

top_Nun = Top_20sp_Nun
top_Nun$Lifestyle  = c(rep("Nunavik", 279))
top_Nun = reshape2::melt(top_Nun)
names(top_Nun)[2] <- "Species"
names(top_Nun)[3] <- "Abundance"

Top_20sp_nwest = rbind(Top_20sp_MDG, Top_20sp_PRU, Top_20sp_FJI)
top_nwest = Top_20sp_nwest
top_nwest$Lifestyle  = c(rep("Non-industrialized", 73))
top_nwest = reshape2::melt(top_nwest)
names(top_nwest)[2] <- "Species"
names(top_nwest)[3] <- "Abundance"

Top_20sp_west = rbind(Top_20sp_NDL, Top_20sp_CND, Top_20sp_USA, Top_20sp_SPN)   
top_west =   Top_20sp_west 
top_west$Lifestyle  = c(rep("Industrialized", 104))
top_west = reshape2::melt(top_west)
names(top_west)[2] <- "Species"
names(top_west)[3] <- "Abundance"

top_sp =   rbind(top_Nun, top_nwest, top_west) 

x = top_sp$Lifestyle
y = log(top_sp$Abundance +1)

top_sp$Lifestyle <- factor(top_sp$Lifestyle, levels = c("Nunavik", "Non-industrialized", "Industrialized"))
Lifestyle = top_sp$Lifestyle

# Wilcoxon-Mann_Whitney calculation
stat.test <- top_sp %>%
  group_by(Species) %>%
  wilcox_test(Abundance ~ Lifestyle)

stat.test <- stat.test %>%
  add_xy_position(x = "Lifestyle")

ggplot(top_sp, aes(x, y, color = Lifestyle)) +
  geom_violin(scale = "width", draw_quantiles = c(0.25, 0.5, 0.75), size = 1) +
  theme_gray()  +
  theme( legend.text = element_text(size = 20),legend.title = element_text(size= 20),title = element_text(size= 20),
         text = element_text(size=12),axis.text.x = element_blank(), strip.text = element_text(face = "bold.italic"), 
         strip.background = element_rect(fill="white"), plot.margin = margin(0.1, 0.1, 0.1, 0.1, "cm"))+
  labs(title = "Lifestyle")+
  ylab(paste("Relative abundance log(y+1)"))+
  scale_x_discrete(limits = c("Nunavik", "Non-industrialized", "Industrialized"))+
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.1)))+
  xlab(paste("")) +
  labs(colour = "")+
  ggtitle("Relative abundance - Top 20 most abuandant species in Nunavik metagenomes") +
  stat_pvalue_manual(stat.test, label = "p.adj.signif", tip.length = 0.01, bracket.nudge.y = 0, hide.ns = T,  bracket.size = 0.6, label.size = 6,
                     y.position = c(c (3, 3.5, 4),c(3.5, 4), c(3, 3.5, 4), c(3.5, 4), c(3), 
                                    c(3.5, 4, 4.5), c( 3.5, 4, 4.5), c(3.5, 4, 4.5),c( 4,4.5,5), c(2.5, 3), 
                                    c(2.5), c(2.5, 3, 3.5), c(2.5, 3),c( 3, 3.5, 4), c( 3, 3.5, 4),
                                    c(2.7),c(2.1, 2.6, 3.1),c(2, 2.5),c(2.7,3.2)))+
  facet_wrap(vars(Species) , scales = "free", nrow = 4, ncol = 5)+
  scale_color_manual( values = c( "#E69F00","#CC79A7", "#56B4E9" ))

# RF identified signature species relative abundance graphic representation
# data wrangling 
RF_NDL = sp_gini[1:41,]
RF_MDG = sp_gini[42:74,]
RF_CND = sp_gini[75:95,]
RF_PRU = sp_gini[96:100,]
RF_USA = sp_gini[101:132,]
RF_SPN = sp_gini[133:142,]
RF_FJI = sp_gini[143:177,]

RF_Nun = as.data.frame(sp_gini[178:456,])
RF_Nun$Lifestyle  = c(rep("Nunavik", 279))
RF_Nun = reshape2::melt(RF_Nun)
names(RF_Nun)[2] <- "Species"
names(RF_Nun)[3] <- "Abundance"
test_nun =as.data.frame(sp_gini[178:456,]) 
test_nun = as.data.frame(apply(test_nun,2,median))
names(test_nun)[1] <- "Nunavik"

RF_nwest = as.data.frame(rbind((RF_MDG), (RF_PRU), (RF_FJI)))
RF_nwest$Lifestyle  = c(rep("Non-industrialized", 73))
RF_nwest = reshape2::melt(RF_nwest)
names(RF_nwest)[2] <- "Species"
names(RF_nwest)[3] <- "Abundance"
test_non_west = as.data.frame(rbind((RF_MDG), (RF_PRU), (RF_FJI))) 
test_non_west = as.data.frame(apply(test_non_west,2,median))
names(test_non_west)[1] <- "Non-industrialized"

RF_west = as.data.frame(rbind((RF_NDL), (RF_CND), (RF_USA), (RF_SPN)))    
test_west = as.data.frame(rbind((RF_NDL), (RF_CND), (RF_USA), (RF_SPN))) 
test_west = as.data.frame(apply(test_west,2,median))
names(test_west)[1] <- "Industriel"
RF_west$Lifestyle  = c(rep("Industrialized", 104))
RF_west = reshape2::melt(RF_west)
names(RF_west)[2] <- "Species"
names(RF_west)[3] <- "Abundance"

RF_sp =   rbind(RF_Nun, RF_nwest, RF_west) 
median_gini = cbind((test_nun), (test_non_west),(test_west))

x = RF_sp$Lifestyle
y = log(RF_sp$Abundance+1)

RF_sp$Lifestyle <- factor(RF_sp$Lifestyle, levels = c("Nunavik", "Non-industrialized", "Industrialized"))
Lifestyle = RF_sp$Lifestyle

# Wilcoxon-Mann-Whitney test
stat.test <- RF_sp %>%
  group_by(Species) %>%
  wilcox_test(Abundance ~ Lifestyle)

stat.test <- stat.test %>%
  add_xy_position( x = "Lifestyle")

ggplot(RF_sp, aes(x, y, color = Lifestyle)) +
  geom_violin(scale = "width", draw_quantiles = c(0.25, 0.5, 0.75), size = 1, position = position_dodge(width = 1)) +
  theme_gray() +
  theme( legend.text = element_text(size = 20),legend.title = element_text(size= 20),title = element_text(size= 20),
         text = element_text(size=12),axis.text.x = element_blank(), strip.text = element_text(face = "bold.italic"), 
         strip.background = element_rect(fill="white"), plot.margin = margin(0.1, 0.1, 0.1, 0.1, "cm"))+
  labs(title = "Lifestyle")+
  ylab(paste("Relative abundance log(y+1)"))+
  xlab(paste("")) +
  labs(colour = "")+
  scale_x_discrete(limits = c("Nunavik", "Non-industrialized", "Industrialized"))+
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.1)))+
  ggtitle("Relative abundance - RF identified species in Nunavik metagenomes") +
  stat_pvalue_manual(stat.test, label = "p.adj.signif", tip.length = 0.01, bracket.nudge.y = 0, hide.ns = T,  bracket.size = 0.6, label.size = 6,
                     y.position =   c(c( 1.6, 2), c(2.2,2.7,3), c(3, 3.5,2.5), c(0.35, 0.45), c(1, 1.25,1.8),
                                      c(0.4, 0.65), c(0.55, 0.77), c(3.5, 3.7), c(2, 2.5,3), c( 0.0075, 0.01) ,
                                      c(2.1, 3),c(1.7, 2.2,2.7), c(3.2, 3.5), c( 1, 1.25,0.75), c(1.1, 1.26,0.9), 
                                      c(3.5, 4,3),c(0.7,0.6), c(1.25,1.1), c(0.4,0.5,0.6), c(2.2, 2.7)))+
  facet_wrap(vars(Species) , scales = "free_y", ncol = 5, nrow =4 )+
  scale_color_manual( values = c( "#E69F00","#CC79A7", "#56B4E9" ))



# PCoA on relative abundance of RF identified signature species
# start by scaling data
sc = scale(sp_gini, center = F, scale = T)

# Generate Jaccard similarity matrix
Jaccard_gini = vegdist(sc, method = "jaccard", binary = T, diag = F,
                       upper = F, na.rm = T)

# perform PcoA
pcoa_jacc_gini = cmdscale(Jaccard_gini,  eig = T, k = 4, x.ret = T)

# strucure data 
mds.var.per.RF <- round(abs(pcoa_jacc_gini$eig)/sum(abs(pcoa_jacc_gini$eig))*100, 1)
mds.data.RF <- pcoa_jacc_gini$points
colnames(mds.data.RF) <- c("PCoA1", "PCoA2", "PCoA3", "PCoA4")
mds.data.RF = as.data.frame(mds.data.RF)
Lifestyle = Metadata_Complete$Lifestyle

X= mds.data.RF$PCoA1
Y= mds.data.RF$PCoA2

ggplot(mds.data.RF, aes(x=X, y=Y, color=Lifestyle)) +
  theme_gray() +
  geom_point() +
  theme(text = element_text(size=10))+
  labs(colour = "")+
  xlab(paste("PCoA1 - ", mds.var.per.RF[1], "%", sep="")) +
  ylab(paste("PCoA2 - ", mds.var.per.RF[2], "%", sep="")) +
  ggtitle("PCoA on Jaccard distance - Relative abundance of RF identified species") +
  scale_color_manual(breaks=c("Nunavik","Non-industrialized",
                              "Industrialized"), values=c( "#E69F00","#CC79A7", "#56B4E9" ))+
  stat_conf_ellipse(
    geom = "path",
    position = "identity",
    na.rm = FALSE,
    inherit.aes = TRUE,
    level = 0.95,
    npoint = 100,
    bary = TRUE
  )

### PERMANOVA calculation
permanova_J_gini  = adonis2(Jaccard_gini ~ vectlifestyle + vectage + vectsex , 
                            permutations = 999, method = "jaccard")
permanova_J_gini


# Top 20 most abundant pathways graphic representation


# data wrangling
Top_20p_Nun = Top_20_p[1:279,]
Top_20p_FJI = Top_20_p[280:314,]
Top_20p_USA =  Top_20_p[315:336,]
USA = Top_20_p[348:357,]
Top_20p_USA = rbind(Top_20p_USA ,USA)
Top_20p_SPN = Top_20_p[337:346,]
Top_20p_PRU = Top_20_p[347,]
PRU = Top_20_p[358:361,]
Top_20p_PRU = rbind(Top_20p_PRU, PRU)
Top_20p_MDG = Top_20_p[362:394,]
Top_20p_CND = Top_20_p[395:415,]
Top_20p_NDL = Top_20_p[416:456,]

top_Nun_p = Top_20p_Nun
top_Nun_p$Lifestyle  = c(rep("Nunavik", 279))
top_Nun_p = reshape2::melt(top_Nun_p)
names(top_Nun_p)[2] <- "Pathways"
names(top_Nun_p)[3] <- "Abundance"

Top_20p_nwest = rbind(Top_20p_MDG, Top_20p_PRU, Top_20p_FJI)
top_nwest_p = Top_20p_nwest
top_nwest_p$Lifestyle  = c(rep("Non-industrialized", 73))
top_nwest_p = reshape2::melt(top_nwest_p)
names(top_nwest_p)[2] <- "Pathways"
names(top_nwest_p)[3] <- "Abundance"

Top_20p_west = rbind(Top_20p_NDL, Top_20p_CND, Top_20p_USA, Top_20p_SPN)   
top_west_p =   Top_20p_west 
top_west_p$Lifestyle  = c(rep("Industrialized", 104))
top_west_p = reshape2::melt(top_west_p)
names(top_west_p)[2] <- "Pathways"
names(top_west_p)[3] <- "Abundance"

top_p =   rbind(top_Nun_p, top_nwest_p, top_west_p) 



# grapgic representation 
x = top_p$Lifestyle
y = top_p$Abundance

top_p$Lifestyle <- factor(top_p$Lifestyle, levels = c("Nunavik", "Non-industrialized", "Industrialized"))
Lifestyle = top_p$Lifestyle

# Wilcoxon-Mann-Whitney test 
stat.test <- top_p %>%
  group_by(Pathways) %>%
  wilcox_test(Abundance ~ Lifestyle)

stat.test <- stat.test %>%
  add_xy_position(x = "Lifestyle")

ggplot(top_p, aes(x, y, color = Lifestyle)) +
  geom_violin(scale = "width", draw_quantiles = c(0.25, 0.5, 0.75), size = 1) +
  theme_gray() +
  theme( legend.text = element_text(size = 20),legend.title = element_text(size= 20),title = element_text(size= 20),
         text = element_text(size=12),axis.text.x = element_blank(), strip.text = element_text(face = "bold"), 
         strip.background = element_rect(fill="white"), plot.margin = margin(0.1, 0.1, 0.1, 0.1, "cm"))+
  labs(title = "Lifestyle")+
  ylab(paste("Relative abundance (cpm)"))+
  scale_x_discrete(limits = c("Nunavik", "Non-industrialized", "Industrialized"))+
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.1)))+
  xlab(paste("")) +
  labs(colour = "")+
  ggtitle("Relative abundance (cpm) - Top 20 most abundant pathways in Nunavik metagenomes") +
  stat_pvalue_manual(stat.test, label = "p.adj.signif", tip.length = 0.01, bracket.nudge.y = 0, hide.ns = T, bracket.size = 0.6,label.size = 6,
                     y.position = c(c(1200,1400),c(1250,1450,1650), c(1450, 1650), c(1250,1450,1650),c(1450,1550), 
                                    c(1450,1650), c(1200, 1400),c( 1200,1400 ),c(1200,1400),c(1200,1400),
                                    c( 1300, 1500), c( 1250,1450), c( 1100, 1200), c(1050,1200), c(1100,1300), 
                                    c(1250,1450),c( 1300,1500), c( 1300, 1500), c( 1350,1550) ))+
  facet_wrap(vars(Pathways) , scales = "free_y")+
  scale_color_manual( values = c( "#E69F00","#CC79A7", "#56B4E9" ))

# Rf identified signature pathways relative abundance graphic representation
# data wrangling 
p_RF_FJI = p_RF[280:314,]
p_RF_USA =  p_RF[315:336,]
USA = p_RF[348:357,]
p_RF_USA = rbind(p_RF_USA ,USA)
p_RF_SPN = p_RF[337:346,]
p_RF_PRU = p_RF[347,]
PRU = p_RF[358:361,]
p_RF_PRU = rbind(p_RF_PRU, PRU)
p_RF_MDG = p_RF[362:394,]
p_RF_CND = p_RF[395:415,]
p_RF_NDL = p_RF[416:456,]
p_RF_NUN = p_RF[1:279,]
p_RF_NUN = as.data.frame(p_RF_NUN)
p_RF_NUN$Lifestyle  = c(rep("Nunavik", 279))
p_RF_NUN = reshape2::melt(p_RF_NUN)
names(p_RF_NUN)[2] <- "Pathways"
names(p_RF_NUN)[3] <- "Abundance"

p_RF_nwest = rbind(p_RF_MDG, p_RF_PRU, p_RF_FJI)
p_RF_nwest$Lifestyle  = c(rep("Non-industrialized", 73))
p_RF_nwest = reshape2::melt(p_RF_nwest)
names(p_RF_nwest)[2] <- "Pathways"
names(p_RF_nwest)[3] <- "Abundance"

p_RF_west = rbind(p_RF_NDL, p_RF_CND, p_RF_USA, p_RF_SPN)   
p_RF_west$Lifestyle  = c(rep("Industrialized", 104))
p_RF_west = reshape2::melt(p_RF_west)
names(p_RF_west)[2] <- "Pathways"
names(p_RF_west)[3] <- "Abundance"

RF_p =   rbind(p_RF_NUN, p_RF_nwest, p_RF_west) 

# graphic representation 
x = RF_p$Lifestyle
y = RF_p$Abundance

RF_p$Lifestyle <- factor(RF_p$Lifestyle, levels = c("Nunavik", "Non-industrialized", "Industrialized"))
Lifestyle = RF_p$Lifestyle

# Wilcoxon-Mann-Whitney test
stat.test <- RF_p %>%
  group_by(Pathways) %>%
  wilcox_test(Abundance ~ Lifestyle)

stat.test <- stat.test %>%
  add_xy_position(x = "Lifestyle")

ggplot(RF_p, aes(x, y, color = Lifestyle)) +
  geom_violin(scale = "width", draw_quantiles = c(0.25, 0.5, 0.75), size = 1) +
  theme_gray() +
  theme( legend.text = element_text(size = 20),legend.title = element_text(size= 20),title = element_text(size= 20),
         text = element_text(size=12),axis.text.x = element_blank(), strip.text = element_text(face = "bold"), 
         strip.background = element_rect(fill="white"), plot.margin = margin(0.1, 0.1, 0.1, 0.1, "cm"))+
  labs(title = "Lifestyle")+
  ylab(paste("Relative abundance (cpm)"))+
  scale_x_discrete(limits = c("Nunavik", "Non-industrialized", "Industrialized"))+
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.1)))+
  xlab(paste("")) +
  labs(colour = "")+
  ggtitle("Relative abundance (cpm) - RF identified pathways in Nunavik metagenomes") +
  stat_pvalue_manual(stat.test, label = "p.adj.signif", tip.length = 0.01, bracket.nudge.y = 0, hide.ns = T, bracket.size = 0.6,label.size = 7,
                     y.position = c(c(180,230),c(25,35,45), c(1050, 1150), c(100,125,150),c(120,150, 175), 
                                    c(15,20,25), c(15,20,25),c( 15,20,25 ),c(15,20,25),c(16,21,26),
                                    c( 380, 480,580), c( 550,650), c( 280, 320), c(120,140), c(150,175), 
                                    c(25,30,35), c(200,250),c( 140,170), c(25,30,35), c( 200,250,300)))+
  facet_wrap(vars(Pathways) , scales = "free_y")+
  scale_color_manual( values = c( "#E69F00","#CC79A7", "#56B4E9" ))