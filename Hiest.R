## 1. generate 012 file with plink
## 193
### vcftools --vcf ../193loci_higherFST_NovosNomes-bugios_samp_sng1_SNP.recode.vcf --out 193loci_fst1.vcftools --012
## 79
### vcftools --vcf ../NovosNomes-bugios_allsamples_07FST_40missingInd.recode.vcf --out 79loci_bugios179 --012

# R CODE
# 79 LOCI
library(HIest)
setwd("~/Desktop")

GdataQ <- read.table("./79loci/79loci_fst1.vcftools.012")
PdataQ <- read.table("./79loci/79loci_AlelleDiag.txt", header=T)

HIC_data <- HIC(Gdata193)
write.table(HIC_data, file = "./HIC_data.it1000.txt", quote=F)
HIC_data

# LINHAS: following https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0127575
HybridDetection <- HIest(GdataQ, PdataQ,type="allele.count", method="SANN",
                         iterations=1000,surf=TRUE)
write.csv(HybridDetection,'./HybridDetection_nogrid_methodSANN_1000int.csv',quote=F)
#write.csv(full.HI.long,'full.HI_long.txt')
#2 version is just the redo
# # optional plot
### HybridDetection = read.csv('HybridDetection.csv')
head(HybridEst)

pdf("./HybridDetection_nogrid_methodSANN_1000int.pdf",5,5)
ggplot(Hybrid79) +
  geom_polygon(data = data.frame(x = c(0, 0.5, 1, 0), y = c(0, 1, 0, 0)),
               aes(x, y), fill = "white", color = "black", size = 0.5) +
  geom_point(data = Hybrid79,aes(x = S, y = H), shape = 19, size = 4, color = "red2") +
  labs(x = expression(italic(S)), y = expression(italic(H[I]))) +
  theme_bw() +
  theme(panel.background = element_rect(fill = "white"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  scale_x_continuous(breaks = c(1, 0.8,0.6,0.4,0.2,0)) +
  scale_y_continuous(breaks = c(1, 0.8,0.6,0.4,0.2,0))
dev.off()

# HIclass to calculate likelihoods for each of the six possible genotype classes in two 
# generations of hybridization (each parental species, F1, F2, and backcrosses, toward 
# each parental species)
HybridClasses <- HIclass(GdataQ, PdataQ, type="allele.count") 
write.csv(HybridClasses,'./HybridClasses_nogrid_methodSANN_1000int.csv')

# # compare classification with maximum likelihood estimates
Hybrid.comp <- HItest(HybridClasses,HybridDetection)
write.csv(Hybrid.comp,'./Hybrid.comp_nogrid_methodSANN_1000int.csv')

# To test which estimation is better 
# allowed us to decide whether the simple classification assuming 
# an early hybridization system is acceptable if its AIC (Akaike Information Criterion) 
# is lower than the AIC of the S and HI maximum log-likelihood estimates
HItest<-HItest(HybridClasses, HybridDetection, thresholds = c(2, 1))
write.csv(HItest,'./HItest_nogrid_methodSANN_1000int.csv')

# 192 LOCI
# REF: https://github.com/kyleaoconnell22/sao_tome_caecilians/blob/main/HIest_script_filt.R
# https://onlinelibrary.wiley.com/doi/full/10.1111/mec.15928?casa_token=3J_YRMd9WsgAAAAA%3A6hQhpU9QZNYr2oFE3_1r6XzQ9HPjXTDdjFVI1AJu2QTp3G87JfPOTczkGYcsqzg01DIqvuJNy4msZw
### admixture file 
###########################################################################
###########################################################################
###########################################################################
################ RUN WITH THE RIGHT DIAGNOSIS LOCI ########################
###########################################################################
###########################################################################
###########################################################################

library(HIest)
setwd("~/Desktop")

Gdata193 <- read.table("./192loci_fst1.vcftools.012")
Pdata193 <- read.table("./192loci_AlelleDiag.txt", header=T)

HIC_data <- HIC(Gdata193)
write.table(HIC_data, file = "./193loci/diagnostic_alelle/HIC_data.it1000.txt", quote=F)
HIC_data

# following https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0127575
HybridDetection <- HIest(Gdata193, Pdata193,type="allele.count", method="SANN",
                         iterations=1000,surf=TRUE)

write.csv(HybridDetection,'./HybridDetection_nogrid_methodSANN_1000int.csv',quote=F)
#write.csv(full.HI.long,'full.HI_long.txt')
#2 version is just the redo
# # optional plot
### HybridDetection = read.csv('HybridDetection.csv')
head(HybridEst)

pdf("./HybridDetection_nogrid_methodSANN_1000int.pdf",5,5)
ggplot(HybridDetection) +
  geom_polygon(data = data.frame(x = c(0, 0.5, 1, 0), y = c(0, 1, 0, 0)),
               aes(x, y), fill = "white", color = "black", size = 0.5) +
  geom_point(data = HybridDetection,aes(x = S, y = H), shape = 19, size = 4, color = "red2") +
  labs(x = expression(italic(S)), y = expression(italic(H[I]))) +
  theme_bw() +
  theme(panel.background = element_rect(fill = "white"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  scale_x_continuous(breaks = c(1, 0.8,0.6,0.4,0.2,0)) +
  scale_y_continuous(breaks = c(1, 0.8,0.6,0.4,0.2,0))
dev.off()

# HIclass to calculate likelihoods for each of the six possible genotype classes in two 
# generations of hybridization (each parental species, F1, F2, and backcrosses, toward 
# each parental species)
HybridClasses <- HIclass(Gdata193, Pdata193, type="allele.count") 
write.csv(HybridClasses,'./HybridClasses_nogrid_methodSANN_1000int.csv')

# # compare classification with maximum likelihood estimates
Hybrid.comp <- HItest(HybridClasses,HybridDetection)
write.csv(Hybrid.comp,'./Hybrid.comp_nogrid_methodSANN_1000int.csv')

# To test which estimation is better 
# allowed us to decide whether the simple classification assuming 
# an early hybridization system is acceptable if its AIC (Akaike Information Criterion) 
# is lower than the AIC of the S and HI maximum log-likelihood estimates
HItest<-HItest(HybridClasses, HybridDetection, thresholds = c(2, 1))
write.csv(HItest,'./HItest_nogrid_methodSANN_1000int.csv')

## PLOTS ##
library(HIest)
library(ggplot2)
library(ggthemes)

setwd("./192_HIest/")

# Load data
# change , for . 
# have two columns determined the colours of each category
HIclass192 <- read.csv(file="HItest_nogrid_1000int_inorder_points.csv",header = TRUE,row.names = 1,sep=';')

HIclass79 <- read.csv(file="../79_HIest/HItest_nogrid_1000int_inorder_point.csv",header = TRUE,
                      row.names = 1,sep=';')


# 192 loci
pdf(paste("../../Desktop/192_pures.pdf", sep=""), width=6, height=5)
p <- ggplot(HIclass192,aes(1-S, H)) +
  geom_point(position = position_jitter(seed = 0.8, width = 0.02, height=0.03),
             shape=21, size=4, colour = factor(HIclass192$Fill), fill=factor(HIclass192$Colour)) + 
  scale_color_manual(values = sample_colors) +
  scale_y_continuous(expand = expansion(c(0.08, 0.2)))+
  geom_segment(x = 0, y = 0, xend =0.5 , yend = 1) +
  geom_segment(x = 0.5, y = 1, xend =1 , yend = 0)  +
  geom_segment(x = 0, y = 0, xend =1 , yend = 0) +
  labs(x = "Hybrid index", y = "Interspecific heterozygosity")
p + theme(
  panel.background = element_rect(fill = "white"),
  panel.border = element_blank(),
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank())
dev.off()

# 79 loci
pdf(paste("../../Desktop/79_pures.pdf", sep=""), width=6, height=5)
p <- ggplot(HIclass79, aes(1-S,H)) +
  geom_point(position = position_jitter(seed = 1, width = 0.02, height=0.03), 
             shape=21, size=4, colour = factor(HIclass79$Fill), fill=factor(HIclass79$Colour))+ 
  scale_y_continuous(expand = expansion(c(0.07, 0.07)))+
  geom_segment(x = 0, y = 0, xend =0.5 , yend = 1) +
  geom_segment(x = 0.5, y = 1, xend =1 , yend = 0)  +
  geom_segment(x = 0, y = 0, xend =1 , yend = 0) +
  labs(x = "Hybrid index", y = "Interspecific heterozygosity")
p + theme(
  panel.background = element_rect(fill = "white"),
  panel.border = element_blank(),
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank())
dev.off()
