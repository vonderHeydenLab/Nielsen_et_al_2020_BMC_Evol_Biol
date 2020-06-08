################################################################################
# Scripts to perfrom gene-environment asscociation analyses as described in
# Nielsen et al. (in press) BMC Evolutionary Biology
#
# Adapted and/or written by ES Nielsen, Stellenbosch University, South Africa
#
# Code is provided as is, without support 
################################################################################

################################################# Getting Nei's genetic distance:
library(StAMPP)
setwd("~/Desktop/PhD_stuffies/CH1/min.cov.20.analyses/Stampp")
PA.stampp <- read.csv("~/PA.stampp.csv", header = T, sep = ";")
PA.s.file <- stamppConvert(PA.stampp, type = "r")
PA.d <- stamppNeisD(PA.s.file, pop = TRUE)
write.table(PA.d, file="PA.Nei.d.txt", col.names = F, row.names= F, quote = F)

###################################################################### Plot PCAs:
library(vegan)
library(ggplot2)

# Read allele frequencies
af.data<-read.table("CP.LD.AF.txt",sep="\t",header=T)

#Give sample names
cpnames <- c("JB","YZ","SP","BT","GB","CA","JF","MB","KY","CF","HH","CB","MG")
#PA
panames <- c("PN", "HB", "DB", "LB", "JB", "SP", "BT", "CA", "MB", "KY", "CF", "PA", "HH", "CB")
#SG
sgnames <- c("PN", "HB", "BB", "LB", "JB", "SP", "BT", "CA", "MB", "KY", "CF", "PA", "HH", "HL")
num <- c("1","2","3","4","5","6","7","8","9","10","11","12","13")
num <- c("1","2","3","4","5","6","7","8","9","10","11","12","13","14")
num <- as.numeric(num)

#transform data

## import as Excel, skip the first column which is contig_pos
CP.af.t <- t(CP_pstat_AF)

CP.af.t <- t(af.out.subset[,2:14])

#run pca
pca_data=prcomp(CP.af.t,
                center = TRUE,
                scale. = TRUE) 

pca_data=prcomp(CP.af.t)

#get % variance
pca_data_perc=round(100*pca_data$sdev^2/sum(pca_data$sdev^2),1)

#create dataframe (!!! make sure to change (CP_pstat_AF) to working df !!!!)
df_pca_data = data.frame(PC1 = pca_data$x[,1], PC2 = pca_data$x[,2], sample = colnames(CP_pstat_AF))

#add number to color by
df_pca_data$site <- num

#ggplot pc1 by pc2 (!!! make sure to change ntnames between species !!!)
PCA<-ggplot(df_pca_data, aes(PC1,PC2, label = ntnames))

#fill in labels and colors
PCA + geom_label(aes(fill = df_pca_data$site), colour = "white", fontface = "bold", size=6) + labs(x=paste0("PC1 (",pca_data_perc[1],")"), y=paste0("PC2 (",pca_data_perc[2],")")) + theme(legend.position="bottom")

#plot neutral and outlier snps per species:
library(ggpubr)
ggarrange(cp.n.pca, pa.n.pca, sg.n.pca, cp.o.pca, pa.o.pca, sg.o.pca, 
          labels = c("a)", "b)", "c)", "d)", "e)", "f)"),
          ncol = 3, nrow = 2)


###################################################################  Mantel tests:

library(vegan)
CP_geo_dis_mat <- read_excel("~/Desktop/PhD_stuffies/CH1/PcoA_mantels/CP.geo.dis.mat.xlsx")

CP.geo.log <- log(CP_geo_dis_mat)

CP_pfst_S_fsts <- read_excel("~/Desktop/PhD_stuffies/CH1/PcoA_mantels/CP.pfst.S.fsts.xlsx")

#with vegan
#mantel(CP_pfst_S_fsts, CP.geo.log, method="pearson", permutations=999)

## with ecodist
library(ecodist)

#read log transformed geographic distance matrix:
CP.geo.mat <- data.matrix(CP.geo.log, rownames.force = NA)

#read slatkin's fst matrix:
CP.pfst.mat <- data.matrix(CP_LD_20_fst, rownames.force = NA)

#Run mantel test
CP.pfst.man <- mantel(as.dist(CP.pfst.mat) ~ as.dist(CP.geo.mat), nperm = 999)
CP.pfst.man

# Read environmental variables, and create matrices per variable:
CP.ss.vars <- read.table("CP.5.env.vars.txt", header = T)
Tr <- dist(CP.ss.vars$Trange, method = "euclidean")
SSSm <- dist(CP.ss.vars$SSSmean, method = "euclidean")
SSSr <- dist(CP.ss.vars$SSSrange, method = "euclidean")
SSTm <- dist(CP.ss.vars$SSTmean, method = "euclidean")
SSTr <- dist(CP.ss.vars$SSTrange, method = "euclidean")

#Run partial mantel tests per variable:
p.T <- mantel(as.dist(CP.pfst.mat) ~ Tr + as.dist(CP.geo.mat), nperm=999)
p.T
p.SSm <- mantel(as.dist(CP.pfst.mat) ~ SSSm + as.dist(CP.geo.mat), nperm=999)
p.SSm
p.STm <- mantel(as.dist(CP.pfst.mat) ~ SSTm + as.dist(CP.geo.mat), nperm=999)
p.STm
p.SSr <- mantel(as.dist(CP.pfst.mat) ~ SSSr + as.dist(CP.geo.mat), nperm=999)
p.SSr
p.STr <- mantel(as.dist(CP.pfst.mat) ~ SSTr + as.dist(CP.geo.mat), nperm=999)
p.STr

# Get significance values:
library(qvalue)
p <- mantel_pvals$CP.pf.1
#If a majority of your p-values are highly significant, maybe the best strategy is to fix pi0=1 (this will apply the BH procedure and is conservative
qobj <- qvalue(p, fdr.level=0.05, pi0 = 1)
qobj
