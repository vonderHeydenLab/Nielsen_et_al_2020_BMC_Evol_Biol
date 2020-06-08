################################################################################
# Scripts to perfrom gene-environment asscociation analyses as described in
# Nielsen et al. (in press) BMC Evolutionary Biology
#
# Adapted and/or written by ES Nielsen, Stellenbosch University, South Africa
#
# Code is provided as is, without support 
################################################################################



######################## DOWNLOAD PRED VARIABLES ################################
##################################################################################

LIB <- c("rgbif", "biomod2", "ggplot2", "gridExtra", "knitr", "raster", 
         "ade4", "rworldmap", "cleangeo", "maptools", "rasterVis", "rgdal","rgeos", "sdmpredictors")
#for(i in LIB) { install.packages(i, repos="http://ftp.sun.ac.za/ftp/pub/mirrors/cran.za.r/") ; library(i, character.only=T) }
for(i in LIB) { library(i, character.only=T) }

library(sdmpredictors)

#download landscape features
a.contp <- load_layers( layercodes = c("WC_bio1", "WC_bio5", "WC_bio6", "WC_bio7", "WC_bio12", "WC_bio13", "WC_bio14", "WC_bio15") , equalarea=FALSE, rasterstack=TRUE)

#download marine features
o.contp <- load_layers( layercodes = c("MS_biogeo08_sss_mean_5m", "MS_biogeo09_sss_min_5m", "MS_biogeo10_sss_max_5m", "MS_biogeo11_sss_range_5m", "MS_biogeo13_sst_mean_5m", "MS_biogeo14_sst_min_5m", "MS_biogeo15_sst_max_5m", "MS_biogeo16_sst_range_5m", "BO_damean", "BO_dissox", "BO_ph") , equalarea=FALSE, rasterstack=TRUE)

#resample because air and ocean vars have different resolution/extent
a<-extent(-180, 180, -90, 90) #input layer extent
o<-extent(-180, 180, -90, 90)
extent_list<-list(a, o)
extent_list<-lapply(extent_list, as.matrix)
matrix_extent<-matrix(unlist(extent_list), ncol=length(extent_list))
rownames(matrix_extent)<-c("xmin", "ymin", "xmax", "ymax")
best_extent<-extent(min(matrix_extent[1,]), max(matrix_extent[3,]), min(matrix_extent[2,]), max(matrix_extent[4,]))
ranges<-apply(as.matrix(best_extent), 1, diff)
reso<-res(a.contp) #choose layer you want to keep resolution
nrow_ncol<-ranges/reso
s2<-raster(best_extent, nrows=nrow_ncol[2], ncols=nrow_ncol[1], crs=a.contp@crs) #choose layer crs you want to keep
a.c.r.2 <-resample(a.contp, s2, method="ngb") #resample by s2
o.c.r.2 <-resample(o.contp, s2, method="ngb") #resample by s2
contemp.r.2=stack(a.c.r.2, o.c.r.2) #stack resampled layers

#chl layer is different extent, need to do the same
chl <- load_layers( layercodes = c("BO2_chlomean_ss"))
sc <-extent(-180, 180, -90, 90)
ch <-extent(-180, 180, -79.5, 90 )
extent_list<-list(sc, ch)
extent_list<-lapply(extent_list, as.matrix)
matrix_extent<-matrix(unlist(extent_list), ncol=length(extent_list))
rownames(matrix_extent)<-c("xmin", "ymin", "xmax", "ymax")
best_extent<-extent(min(matrix_extent[1,]), max(matrix_extent[3,]), min(matrix_extent[2,]), max(matrix_extent[4,]))
ranges<-apply(as.matrix(best_extent), 1, diff)
reso<-res(contemp.r.2)
nrow_ncol<-ranges/reso
s2<-raster(best_extent, nrows=nrow_ncol[2], ncols=nrow_ncol[1], crs=contemp.r.2@crs)
sst.cur <-resample(contemp.r.2, s2, method="ngb")
chl <-resample(chl, s2, method="ngb")
environment.contp=stack(sst.cur, chl)

# Extract env values for our species locations:
envdata <- extract(x=contemp.crop, y=cbind(xy$x, xy$y))

# The `dudi.pca` function allows to perform the PCA over the whole study area.
#this requires library(ade4)
# We decide to keep only 2 principal component axes to summarize the whole environmental niche.
pca1 <- dudi.pca(envdata, scannf = F, nf = 2)
round(pca1$eig/sum(pca1$eig)*100, 2)

# A preliminary test is to look for potential outliers in the environmental data.
plot(pca1$li[,1:2]) # PCA scores on first two axes
summary(pca1$li)

#plot a correlation circle
s.corcircle(pca1$co)

##calculate variance inflation factor (VIF)
library(usdm)
SA.ext <- extent(5, 45, -40, -10)
contemp.crop <- crop(environment.contp, SA.ext)
vif(contemp.crop)

#deselect variables with VIF >10, then run again
VARSEL <- c("WC_bio5", "BO2_tempmax_ss", "BO2_curvelmean_ss", "BO2_chlomean_ss")
contemp.4vars <- stack(subset(contemp.crop, VARSEL))
vif(contemp.4p)

#Write as new raster and extract data per sample site
writeRaster(contemp.3p, filename="contemp.3p.tif", format="GTiff", overwrite=TRUE)

envdata.final <- extract(x=contemp.3p, y=cbind(sg_xy$x, sg_xy$y))


########### If you want to combine certain variables into PCs instead of dropping them: 
###########

library(psych)    # Used to investigate correlations among predictors
library(vegan)
library(factoextra)

#read data
library(readxl)
seascape_variable_selection <- read_excel("~/Desktop/PhD_stuffies/CH1/seascape_variable_selection.xlsx", 
                                          sheet = "Sheet5", col_types = c("numeric", 
                                                                          "numeric", "numeric", "numeric", 
                                                                          "numeric", "numeric", "numeric", 
                                                                          "numeric", "numeric", "numeric", 
                                                                          "numeric", "numeric", "numeric", 
                                                                          "numeric"))
View(seascape_variable_selection)

# veiw correlations

pairs.panels(seascape_variable_selection[,1:14], scale=T)

# run PCA on min/max/mean temp/sss/sst
res.pca <- prcomp(seascape_variable_selection[,1:6], scale = TRUE)

#view dimensions
fviz_eig(res.pca)

#get values
eig.val <- get_eigenvalue(res.pca)
eig.val

#get loadings and create new data frame for further analyses
Loadings <- as.data.frame(res.pca$rotation[,1:4])
axes <- predict(res.pca, newdata = seascape_variable_selection)
View(axes) #these will then be used instead of env variables 



################################# POOLFSTAT ######################################
##################################################################################

library(poolfstat)
library("WriteXLS")
setwd("~/Desktop/PhD_stuffies/CH1/poolfstat/nofilter.syncs")

## read sync file

#SG
psizes <- as.numeric(c('80', '80', '80', '80', '80', '80', '80', '80', '78', '80', '80', '80', '80', '60'))
pnames <- as.character(c('PN', 'HB', 'BB', 'LB', 'JB', 'SP', 'BT', 'CA', 'MB', 'KY', 'CF', 'PA', 'HH', 'HL'))

SG.test.pooldata <- popsync2pooldata(sync.file = "~/Desktop/PhD_stuffies/CH1/poolfstat/nofilter.syncs/CP.dn.b.sync.gz", poolsizes = psizes, poolnames = pnames,
                                     min.rc = 4, min.cov.per.pool = 40, max.cov.per.pool = 400,
                                     min.maf = 0.01, noindel = TRUE, nlines.per.readblock = 1e+06)

## global and per SNP fsts
SG.test.fsts <- computeFST(SG.test.pooldata, method = "Anova", snp.index = NA)

# pairwise FSTs
## I don't think I need to caculate for individual SNP, so will only get 2 output matrices
SG.pair.fst <- computePairwiseFSTmatrix(SG.test.pooldata, method = "Anova",
                                        min.cov.per.pool = 40, max.cov.per.pool = 400,
                                        min.maf = 0.01,
                                        output.snp.values = FALSE)

SG.p.fst <- SG.pair.fst$PairwiseFSTmatrix
SG.p.fst <- as.data.frame(SG.p.fst)
WriteXLS(SG.p.fst, "SG.p.fst.0.xls")

## Convert to BAYPASS
#output= genobaypass (allele count), poolsize, snpdet (snp.info matrix)
pooldata2genobaypass(SG.test.pooldata, writing.dir = getwd(), subsamplesize = -1)



################################# BAYPASS #######################################
##################################################################################

source("~/baypass_utils.R")
require(corrplot) ; require(ape)
library(WriteXLS)

##### Run on HPC:
g_baypass -npop 13 -gfile CP.genobaypass -poolsizefile CP.poolsize -d0yij 8 -outprefix CP.BP -npilot 100

### use the outputs of the above in R
##### In R:
# View omega file
omega=as.matrix(read.table("CP.BP.sim_mat_omega.out"))

pop.names=c("JB","YZ","SP","BT","GB","CA","JF","MB","KY",
            "CF","HH","CB","MG")

#Omega output as heatmap
cor.mat=cov2cor(omega)
corrplot(cor.mat,method="color",mar=c(2,1,2,2)+0.1,
         main=expression("Correlation map based on"~hat(Omega)))

#Omega output as tree
bta14.tree=as.phylo(hclust(as.dist(1-cor.mat**2)))
plot(bta14.tree,type="p",
     main=expression("Hier. clust. tree based on"~hat(Omega)~"("*d[ij]*"=1-"*rho[ij]*")"))

#Plot XTX output
CP.0.snp.res=read.table("CP.BP.sim_summary_pi_xtx.out",h=T)
plot(CP.0.snp.res$M_XtX)

pi.beta.coef=read.table("CP.BP.sim_summary_beta_params.out",h=T)$Mean

CP.0.data<-geno2YN("CP.genobaypass")

#Simulate POD XTX
simu.bta<-simulate.baypass(omega.mat=omega,nsnp=1190,sample.size=CP.0.data$NN,
                           beta.pi=pi.beta.coef,pi.maf=0.01,suffix="CP.BP.pod")

###Run on HPC:
g_baypass -npop 13 -gfile CP.BP.pod -efile cp.bp.5.env.txt -poolsizefile CP.poolsize -d0yij 8 -auxmodel -omegafile CP.BP_mat_omega.out -outprefix CP.BP.5.env -npilot 100

##### In R:
#Plot POD omega:
pod.omega=as.matrix(read.table("CP.BP.POD.sim_mat_omega.out"))
plot(pod.omega,omega) ; abline(a=0,b=1)
fmd.dist(pod.omega,omega)

#Plot POD beta coef:
pod.pi.beta.coef=read.table("CP.BP.POD.sim_summary_beta_params.out",h=T)$Mean
plot(pod.pi.beta.coef,pi.beta.coef) ; abline(a=0,b=1)

pod.xtx=read.table("CP.BP.POD.sim_summary_pi_xtx.out",h=T)$M_XtX

#Plot the original XTX values, and add abline to distinguish outliers.  
pod.thresh=quantile(pod.xtx,probs=0.99)
plot(CP.0.snp.res$M_XtX)
abline(h=pod.thresh,lty=2)

#Save XTX values, and then select those with value > pod.thresh as outliers
CP.snp.scores <- as.data.frame(CP.0.snp.res$M_XtX)
write.table(CP.snp.scores, file="CP.snp.scores.txt", sep="\t")

##### Run on HPC for GEAA BayPass analysis:
g_baypass -npop 13 -CP.genobaypass -poolsizefile CP.poolsize -d0yij 8 -outprefix CP.BP -npilot 100

#check the output *_betai to see the Bayesfactor values 

#################################### RDAs ########################################
##################################################################################
library(vegan)
library(usdm)

## Hellinger transformation of allele frequencies
CP.af.h <- decostand(CP_4_40_400_AF, method="hellinger")
CP.af.h.t <- t(CP.af.h)


setwd("~/Desktop/PhD_stuffies/CH1/RDAs")
#snps <- read.table("CP.4.40.400.AF.txt", header = T)
snps <- read.table("CP.pstat.AF.txt", header = T)
snps.t <- t(snps)
snp.mat <- as.matrix(snps.t)
popnames <- rownames(snp.mat)
popnames
snp.hel <- decostand(snp.mat, method = "hellinger")
setwd("~/Desktop/PhD_stuffies/CH1/RDAs/NO_HP")
CP.ss.vars <- read.table("CP.env.vars.txt", header = T)
rownames(CP.ss.vars)
vif(CP.ss.vars)
env.scale <- scale(CP.ss.vars, scale = T, center = T)

#Run and test significance
rda1 <- rda(snp.hel ~., data = as.data.frame(env.scale))
summary(rda1)
RsquareAdj(rda1)
anova(rda1)
anova(rda1, by = "axis") 

#Plot RDA
plot(rda1, type="n", scaling=3)
points(rda1, display="species", pch=20, cex=1, col="gray32", scaling=3) 
#points(rda1, display="sites", pch=21, cex=1.3, col="gray32", scaling=3)
text(rda1, display="sites", pch=21, cex=1.5, col="firebrick3", scaling=3)
text(rda1, scaling=3, display="bp", col="#0868ac", cex=1.5)

# Identify outliers
outliers <- function(x,z){
  lims <- mean(x) + c(-1, 1) * z * sd(x)                   # find +/- z sd from mean loading     
  x[x < lims[1] | x > lims[2]]   # locus names in these tails
}

load.rda <- summary(rda1)$species[,1:3]
load.rda[,1]
load.rda

cand1.3SD <- outliers(load.rda[,1], 3)
cand2.3SD <- outliers(load.rda[,2], 3)

ncand1.3SD <- length(cand1.3SD)
ncand2.3SD <- length(cand2.3SD)
ncand1.3SD
ncand2.3SD

ncand <- ncand1.3SD+ncand2.3SD  
ncand

cand1.3SD.df <- cbind.data.frame(rep(1, times = length(cand1.3SD)), names(cand1.3SD), unname(cand1.3SD)); colnames(cand1.3SD.df) <- c("axis", "snp", "loading")

cand2.3SD.df <- cbind.data.frame(rep(2, times = length(cand2.3SD)), names(cand2.3SD), unname(cand2.3SD)); colnames(cand2.3SD.df) <- c("axis", "snp", "loading")

cand <- rbind(cand1.3SD.df, cand2.3SD.df)
cand$snp <- as.character(cand$snp)

cand.mat <- matrix(nrow=(ncand), ncol=5)  # ncol = number of predictors
#colnames(cand.mat) <- c("Tmax", "SSSmean", "SSTmean")
colnames(cand.mat) <- c("SSSrange", "SSSmean", "SSTrange", "SSTmean", "Trange")

for (i in 1:length(cand$snp)) {
  nam <- cand[i,2]
  snp.gen <- snp.mat[,nam]
  cand.mat[i,] <- apply(env.scale,2,function(x) cor(x,snp.gen))
}


full.cand.df <- cbind(cand, cand.mat)
full.cand.df

cand$snp[duplicated(cand$snp)]  # check for duplicates

full.cand.df <- full.cand.df[!duplicated(full.cand.df$snp),]

for (i in 1:length(full.cand.df$snp)) {
  bar <- full.cand.df[i,]
  full.cand.df[i,9] <- names(which.max(abs(bar[4:8]))) # the 7 is the column to add to table, 3:6 is the columns of env predictors
  full.cand.df[i,10] <- max(abs(bar[4:8]))              # gives the correlation
}

colnames(full.cand.df)[9] <- "predictor"
colnames(full.cand.df)[10] <- "correlation"

full.cand.df

write.table(full.cand.df, file="CP.full.cand.df.txt", sep="\t")

table(full.cand.df$predictor)  
table(full.cand.df$axis)

### Plotting Outliers

sel <- full.cand.df$snp
env <- full.cand.df$predictor
env[env=="Tmax"] <- '#00CED1'
env[env=="SSSmean"] <- '#FF8C00'
env[env=="SSTmean"] <- '#9932CC'

col.pred <- rownames(rda1$CCA$v) # pull the SNP names

for (i in 1:length(sel)) {           # color code candidate SNPs
  foo <- match(sel[i],col.pred)
  col.pred[foo] <- env[i]
}

col.pred[grep("merged_contig",col.pred)] <- '#f1eef6' # non-candidate SNPs
empty <- col.pred
empty[grep("#f1eef6",empty)] <- rgb(0,1,0, alpha=0) # transparent
empty.outline <- ifelse(empty=="#00FF0000","#00FF0000","gray32")
bg <- c('#00CED1','#FF8C00','#9932CC')

plot(rda1, type="n", scaling=3, xlim=c(-1,1), ylim=c(-1,1))
points(rda1, display="species", pch=21, cex=1.5, col="gray32", bg=col.pred, scaling=3)
points(rda1, display="species", pch=21, cex=1.5, col=empty.outline, bg=empty, scaling=3)
text(rda1, scaling=3, display="sites", col="black", cex=1.5)
text(rda1, scaling=3, display="bp", col="#0868ac", cex=1.5)
legend("bottomleft", legend=c("Tmax","SSSmean","SSTmean"), bty="n", col="gray32", pch=21, cex=1.5, pt.bg=bg)



############################### dbRDAs with dbMEMs #################################
####################################################################################
## Create dbMEM
nsites.per.group <- c(14)
sg.mem <- create.dbMEM.model(coord=sg_snp_coords, nsites=nsites.per.group)

#run rda with just dbMEM
DB.rda <- rda(snp.mat ~., data = as.data.frame(cp.db))

#get r2
R2adj.poly=RsquareAdj(DB.rda)$adj.r.squared

#run forward sel
CP.db.fwd=forward.sel(snp.mat, cp.db, adjR2thresh=R2adj.poly)
CP.db.fwd

summary(DB.rda)
RsquareAdj(DB.rda)
anova(DB.rda)
anova(DB.rda, by = "axis") 
# go on and run same scripts as before to identify outliers...


############################### MSOD/MSR outlier detection ############################
#######################################################################################
library(spdep)
library(adespatial)
library(raster)
library(viridis)
library(dplyr)
library(rgdal)
library(readxl)

#CP_mem_sheet <- read_excel("CP.mem.sheet.xlsx") #this needs to have sites be rows, first 2 columns by xy coords, followed by env vars, followed by loci 
CP_mem_sheet <- read_excel("CP.pst.mem.xlsx")

CP_mem_sheet <- read.table("CP.mem.txt", header = T)

dim(CP_mem_sheet)
#below change Loci to be to the max dimension
Coord <- data.matrix(CP_mem_sheet[,2:3])   # The UTM coordinates 
Env   <- as.data.frame(CP_mem_sheet[,4:6])   # Habitat factor 
Loci  <- data.matrix(CP_mem_sheet[,7:1627]) # Genotypes 

nb      <- graph2nb(gabrielneigh(Coord), sym=TRUE)  # Gabriel graph: neighbor definition
listW   <- nb2listw(nb, style="W")                  # Spatial weights matrix
disttri <- nbdists(nb, Coord)                       # Add longlat=T for lat/long coordinates
fdist   <- lapply(disttri, function(x) x^(-1))      # Use inverse distance weights
listW   <- nb2listw(nb, glist=fdist, style="W")     # Revised spatial weights matrix
tmp     <- scores.listw(listW, MEM.autocor = "all") # Eigenanalysis

mem     <- list(vectors = as.matrix(tmp), values = attr(tmp, "values")) # MEM eigenvectors and eigenvalues
mem$values <- mem$values / abs(sum(mem$values))     # Rescale eigenvalues to Moran's I



setwd("~/Desktop/PhD_stuffies/SDMS/SDM.shps")
map1 <- readOGR("continent.shp")

par(mar=c(1, 1, 1, 1))
plot(map1, axes=F, legend=F, box=F,                 # Plot the selection surface
     col=c("gray70","gray98")) 
plot(nb, coords=Coord, col=1, pch=16, cex=0.8, add=T)  # Add the individuals and Gabriel graph


# Correlations between the Loci and MEM axes
# ------------------------------------------
# Calculate R.YV, which contains for each locus the vector of its correlations with all MEM axes. 
R.YV <- cor(Loci, mem$vectors, use="pairwise.complete.obs")    
S    <- apply(R.YV^2, 2, mean)

# Plot power spectra
## number of bars = number of sites-1 
# ------------------

barplot(S, ylim=c(0, 0.12))             # Average power spectrum 

cutoffs <- abs(qnorm(c(0.05, 0.01)/2))

Dev <- sweep(R.YV^2, 2, S, "/") - 1            # Subtract average power spectrum from each locus.
Dev[Dev > 0] <- 0                              # Set positive deviations to zero.
Dev <- apply(Dev, 1, sum)                      # Sum of negative deviations
z <- scale(Dev)                                # Standardize

## To get axes to show up:
graphics.off()

#Then plot
plot(z, ylim=c(-7,5))                          # Plot the z-scores
for(h in 1:length(cutoffs))                 
{
  lines(c(0,1627), rep(cutoffs[h],2), lty=h) #c(0,#ofloci)
  lines(c(0,1627), rep(-cutoffs[h],2), lty=h)
}


cutoff.msod <- cutoffs[2]                              # Just the middle cutoff of 0.01

Candidates.msod <- c(1:551)[abs(z)>cutoff.msod]        # Candidate loci at this cutoff - THIS GOES FROM 1:(# of LOCI)


cutoff.msr <- 0.05    # Set a less stringent cutoff
nPerm <- 999          # Set number of permutations for MSR test (may choose e.g. 499 or 999)

## This only works with 1 ENV predictor at a time!!!

Env   <- as.data.frame(CP_mem_sheet[,6])


Env2 <- mutate_all(Env, function(x) as.numeric(as.character(x)))
R.XV.Env <- cor(Env2, mem$vectors)
R.XV.xcoord <- cor(Coord[,1], mem$vectors)
R.XV.ycoord <- cor(Coord[,2], mem$vectors)

get.pvalue.msr <- function(r.XV=R.XV, r.YV=R.YV, nPerm=999)
{
  R.XV.rand <- matrix(r.XV, nPerm, ncol(r.XV), byrow=TRUE) 
  R.XV.rand <- R.XV.rand * sample(c(-1,1), length(R.XV.rand), replace=TRUE)
  Cor.obs <- abs(as.vector(r.YV %*% t(r.XV)))
  Cor.rand <- abs(r.YV %*% t(R.XV.rand))
  P.values.MSR <- apply((cbind(Cor.obs,Cor.rand) >= Cor.obs), 1, mean)
  P.values.MSR
}


b.Env <- get.pvalue.msr(r.XV=R.XV.Env, r.YV=R.YV[Candidates.msod,], nPerm=nPerm)

#if this says subscript out of bounds, there is probably a problem with numerical values being read as numeric
#this is because Candidates.msod = number of loci, so = dim(Loci)

print(paste("Loci significantly associated with Env: ", names(b.Env)[b.Env < cutoff.msr]))

b.X <- get.pvalue.msr(r.XV=R.XV.xcoord, r.YV=R.YV[Candidates.msod,], nPerm=nPerm)
b.Y <- get.pvalue.msr(r.XV=R.XV.ycoord, r.YV=R.YV[Candidates.msod,], nPerm=nPerm)



############################## LFMM SCRIPTS #####################################
####################################################################################
library(LEA)
setwd("~/Desktop/PhD_stuffies/CH1/LFMM")
#save MEM xlsx without any labels
#save as txt file
#delete .txt from file name, replace with .lfmm of .env

#CP.lfmm <-read.lfmm("~/Desktop/PhD_stuffies/CH1/LFMM/CP.pst.lfmm")
CP.lfmm <-read.lfmm("~/Desktop/PhD_stuffies/CH1/min.cov.20.analyses/LFMM/CP.lfmm")

#should only run 1 env variable at a time
CP.env <- read.env("~/Desktop/PhD_stuffies/CH1/LFMM/CP.tmax.env")
CP.env <- read.env("~/Desktop/PhD_stuffies/CH1/min.cov.20.analyses/LFMM/CP.tmax.env")

write.lfmm(CP.lfmm, "genotypes.lfmm")

write.env(CP.env, "gradients.env")

pc = pca("genotypes.lfmm", scale = TRUE)

CP.tmax = lfmm( "genotypes.lfmm", "gradients.env", K = 2, repetitions = 10, project = "new")

zs = z.scores(CP.tmax, K = 2, d=1)
zs.median = apply(zs, MARGIN = 1, median)
lambda = median(zs.median^2)/qchisq(0.5, df = 1)
lambda
adj.p.values = pchisq(zs.median^2/lambda, df = 1, lower = FALSE)
hist(adj.p.values, col = "red")

#CP
L = 1190

#PA
L = 822

#SG
L = 1658

q = 0.05
w = which(sort(adj.p.values) < q * (1:L)/L)
candidates.bh = order(adj.p.values)[w]
candidates.bh


