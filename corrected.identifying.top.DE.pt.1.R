# Load libraries####
library(edgeR); library(OmicKriging); library(locfit); library(parallel)
load("~/Desktop/identifying.top.DE.subset.RDA.RData")
source("https://bioconductor.org/biocLite.R")
biocLite("qvalue",lib="/home/arfesta/R/x86_64-pc-linux-gnu-library/3.3/")
#Load counts from asreml output####
load("/media/seagate/Adam.RNA/analyses/combined/tip.full.sample.est.rda")
phenos <- read.table(file="/media/seagate/Adam.RNA/analyses/combined/phenos/all.phenos.192.csv", header=T, sep= ",", row.names=1)
phenos <- phenos[order(phenos$animal, decreasing=F),]
phenos$animal[1:48] <- phenos$animal[1:48]-900; phenos$animal[49:192] <- phenos$animal[49:192]-10000
rownames(phenos) <- phenos$animal ;phenos <- phenos[order(phenos$vol, decreasing=F),]
reads.sample.est <- 2^(tip.full.sample.est)
c <- match(rownames(phenos),rownames(reads.sample.est)); reads.sample.est <- reads.sample.est[c,] #set to be same order as phenos object
identical(rownames(phenos),rownames(reads.sample.est))
phenos$group[which(phenos$group==25)] <- 24
phenos$group[which(phenos$group %in% c(26:57))] <- phenos$group[which(phenos$group %in% c(26:57))] -1

#create DGE object with EdgeR####
groups = phenos$group # groups dictate the family biological reps belong to
table(groups)
# Define vector to record which samples have only a single replicate
singles <- c(1,8,40) 

y <- DGEList(counts=t(reads.sample.est),group=groups)
keep <- rowSums(cpm(y)>75) >=1
table(keep)
y <- y[keep,,keep.lib.sizes=F]
y2 <- calcNormFactors(y); 
#plotMDS(y2[,c(1:2,189:192)])
A <- aveLogCPM(y2); sum(which(A < 1)) # returns 4031145
logcpm2 <- cpm(y2, log=TRUE, prior.count=5, lib.size=y2$samples[,2])
identical(colnames(logcpm2),rownames(phenos)) 


# Use edgeR GLM approach to find genes that are differentially expressed####
# in some families relative to the mean value for all families; try using
# only a subset of differentially-expressed genes for predictions

#Start by creating design matrix:
design <- model.matrix(~0+group, data=y$samples)
colnames(design) <- levels(y$samples$group)
design[1:10,1:10]

y3 <- estimateDisp(y2, design)
fit <- glmFit(y3, design)


###Fam 1 vs. all others####
runit <- function(x){
  if(x == 1){
    end <- 54
    data <- glmLRT(fit, contrast=c(1,-1,rep(0,end)))
    data$table
  } else if(x==55){
    end <- 54
    data <- glmLRT(fit, contrast=c(1,rep(0,end),-1))
    data$table
  } else {
    d <- 54-x
    c <- (x-1)+1
    data <- glmLRT(fit, contrast=c(1,rep(0,c),-1,rep(0,d)))
    data$table}
}
t <- proc.time()
v01.56 <- mclapply(1:55,runit,mc.cores = 55)
proc.time() - t

#now for family 2 vs. all others####
runit <- function(x){
  if(x == 1){
    end <- 53
    data <- glmLRT(fit, contrast=c(0,1,-1,rep(0,end)))
    data$table
  } else if(x==54){
    end <- 53
    data <- glmLRT(fit, contrast=c(0,1,rep(0,end),-1))
    data$table
  } else {
    d <- 53-x
    c <- (x-1)+1
    data <- glmLRT(fit, contrast=c(0,1,rep(0,c),-1,rep(0,d)))
    data$table}
}
t <- proc.time()
v02.56 <- mclapply(1:54,runit,mc.cores = 54)
proc.time() - t

##now for family 3 vs. all others####
runit <- function(x){
  a <-2
  d <- 52-x
  c <- (x-1)+1
  end <- 52
  if(x == 1){
    data <- glmLRT(fit, contrast=c(rep(0,a),1,-1,rep(0,end)))
    data$table
  } else if(x==53){
    data <- glmLRT(fit, contrast=c(rep(0,a),1,rep(0,end),-1))
    data$table
  } else {
    data <- glmLRT(fit, contrast=c(rep(0,a),1,rep(0,c),-1,rep(0,d)))
    data$table}
}
t <- proc.time()
v03.56 <- mclapply(1:53,runit,mc.cores = 53)
proc.time() - t

##now for family 4 vs. all others####
runit <- function(x){
  a <- 3
  d <- 51-x
  c <- (x-1)+1
  end <- 51
  r <- end+1
  if(x == 1){
    data <- glmLRT(fit, contrast=c(rep(0,a),1,-1,rep(0,end)))
    data$table
  } else if(x==r){
    data <- glmLRT(fit, contrast=c(rep(0,a),1,rep(0,end),-1))
    data$table
  } else {
    data <- glmLRT(fit, contrast=c(rep(0,a),1,rep(0,c),-1,rep(0,d)))
    data$table}
}
t <- proc.time()
v04.56 <- mclapply(1:52,runit,mc.cores = 52)
proc.time() - t

##now for family 5 vs. all others####
runit <- function(x){
  a <- 4
  d <- 50-x
  c <- (x-1)+1
  end <- 50
  r <- end+1
  if(x == 1){
    data <- glmLRT(fit, contrast=c(rep(0,a),1,-1,rep(0,end)))
    data$table
  } else if(x==r){
    data <- glmLRT(fit, contrast=c(rep(0,a),1,rep(0,end),-1))
    data$table
  } else {
    data <- glmLRT(fit, contrast=c(rep(0,a),1,rep(0,c),-1,rep(0,d)))
    data$table}
}
t <- proc.time()
v05.56 <- mclapply(1:51,runit,mc.cores = 51)
proc.time() - t

#6####
runit <- function(x){
  a <- 5
  end <- 49
  d <- end-x
  c <- (x-1)+1
  r <- end+1
  if(x == 1){
    data <- glmLRT(fit, contrast=c(rep(0,a),1,-1,rep(0,end)))
    data$table
  } else if(x==r){
    data <- glmLRT(fit, contrast=c(rep(0,a),1,rep(0,end),-1))
    data$table
  } else {
    data <- glmLRT(fit, contrast=c(rep(0,a),1,rep(0,c),-1,rep(0,d)))
    data$table}
}
t <- proc.time()
v06.56 <- mclapply(1:50,runit,mc.cores = 50)
proc.time() - t

#7####
runit <- function(x){
  a <- 6
  end <- 48
  d <- end-x
  c <- (x-1)+1
  r <- end+1
  if(x == 1){
    data <- glmLRT(fit, contrast=c(rep(0,a),1,-1,rep(0,end)))
    data$table
  } else if(x==r){
    data <- glmLRT(fit, contrast=c(rep(0,a),1,rep(0,end),-1))
    data$table
  } else {
    data <- glmLRT(fit, contrast=c(rep(0,a),1,rep(0,c),-1,rep(0,d)))
    data$table}
}
t <- proc.time()
v07.56 <- mclapply(1:49,runit,mc.cores = 49)
proc.time() - t

#8####
runit <- function(x){
  a <- 7
  end <- 47
  d <- end-x
  c <- (x-1)+1
  r <- end+1
  if(x == 1){
    data <- glmLRT(fit, contrast=c(rep(0,a),1,-1,rep(0,end)))
    data$table
  } else if(x==r){
    data <- glmLRT(fit, contrast=c(rep(0,a),1,rep(0,end),-1))
    data$table
  } else {
    data <- glmLRT(fit, contrast=c(rep(0,a),1,rep(0,c),-1,rep(0,d)))
    data$table}
}
t <- proc.time()
v08.56 <- mclapply(1:48,runit,mc.cores = 48)
proc.time() - t

#9####
runit <- function(x){
  a <- 8
  end <- 46
  d <- end-x
  c <- (x-1)+1
  r <- end+1
  if(x == 1){
    data <- glmLRT(fit, contrast=c(rep(0,a),1,-1,rep(0,end)))
    data$table
  } else if(x==r){
    data <- glmLRT(fit, contrast=c(rep(0,a),1,rep(0,end),-1))
    data$table
  } else {
    data <- glmLRT(fit, contrast=c(rep(0,a),1,rep(0,c),-1,rep(0,d)))
    data$table}
}
t <- proc.time()
v09.56 <- mclapply(1:47,runit,mc.cores = 47)
proc.time() - t

#10####
runit <- function(x){
  a <- 9
  end <- 45
  d <- end-x
  c <- (x-1)+1
  r <- end+1
  if(x == 1){
    data <- glmLRT(fit, contrast=c(rep(0,a),1,-1,rep(0,end)))
    data$table
  } else if(x==r){
    data <- glmLRT(fit, contrast=c(rep(0,a),1,rep(0,end),-1))
    data$table
  } else {
    data <- glmLRT(fit, contrast=c(rep(0,a),1,rep(0,c),-1,rep(0,d)))
    data$table}
}
t <- proc.time()
v10.56 <- mclapply(1:46,runit,mc.cores = 46)
proc.time() - t

#11####
runit <- function(x){
  a <- 10
  end <- 44
  d <- end-x
  c <- (x-1)+1
  r <- end+1
  if(x == 1){
    data <- glmLRT(fit, contrast=c(rep(0,a),1,-1,rep(0,end)))
    data$table
  } else if(x==r){
    data <- glmLRT(fit, contrast=c(rep(0,a),1,rep(0,end),-1))
    data$table
  } else {
    data <- glmLRT(fit, contrast=c(rep(0,a),1,rep(0,c),-1,rep(0,d)))
    data$table}
}
t <- proc.time()
v11.56 <- mclapply(1:45,runit,mc.cores = 45)
proc.time() - t

##12####
runit <- function(x){
  a <- 11
  end <- 43
  d <- end-x
  c <- (x-1)+1
  r <- end+1
  if(x == 1){
    data <- glmLRT(fit, contrast=c(rep(0,a),1,-1,rep(0,end)))
    data$table
  } else if(x==r){
    data <- glmLRT(fit, contrast=c(rep(0,a),1,rep(0,end),-1))
    data$table
  } else {
    data <- glmLRT(fit, contrast=c(rep(0,a),1,rep(0,c),-1,rep(0,d)))
    data$table}
}
t <- proc.time()
v12.56 <- mclapply(1:44,runit,mc.cores = 44)
proc.time() - t

#13####
runit <- function(x){
  a <- 12
  end <- 42
  d <- end-x
  c <- (x-1)+1
  r <- end+1
  if(x == 1){
    data <- glmLRT(fit, contrast=c(rep(0,a),1,-1,rep(0,end)))
    data$table
  } else if(x==r){
    data <- glmLRT(fit, contrast=c(rep(0,a),1,rep(0,end),-1))
    data$table
  } else {
    data <- glmLRT(fit, contrast=c(rep(0,a),1,rep(0,c),-1,rep(0,d)))
    data$table}
}
t <- proc.time()
v13.56 <- mclapply(1:43,runit,mc.cores = 43)
proc.time() - t

#14####
runit <- function(x){
  a <- 13
  end <- 41
  d <- end-x
  c <- (x-1)+1
  r <- end+1
  if(x == 1){
    data <- glmLRT(fit, contrast=c(rep(0,a),1,-1,rep(0,end)))
    data$table
  } else if(x==r){
    data <- glmLRT(fit, contrast=c(rep(0,a),1,rep(0,end),-1))
    data$table
  } else {
    data <- glmLRT(fit, contrast=c(rep(0,a),1,rep(0,c),-1,rep(0,d)))
    data$table}
}
t <- proc.time()
v14.56 <- mclapply(1:42,runit,mc.cores = 42)
proc.time() - t

#15####
runit <- function(x){
  a <- 14
  end <- 40
  d <- end-x
  c <- (x-1)+1
  r <- end+1
  if(x == 1){
    data <- glmLRT(fit, contrast=c(rep(0,a),1,-1,rep(0,end)))
    data$table
  } else if(x==r){
    data <- glmLRT(fit, contrast=c(rep(0,a),1,rep(0,end),-1))
    data$table
  } else {
    data <- glmLRT(fit, contrast=c(rep(0,a),1,rep(0,c),-1,rep(0,d)))
    data$table}
}
t <- proc.time()
v15.56 <- mclapply(1:41,runit,mc.cores = 41)
proc.time() - t

#16####
runit <- function(x){
  a <- 15
  end <- 39
  d <- end-x
  c <- (x-1)+1
  r <- end+1
  if(x == 1){
    data <- glmLRT(fit, contrast=c(rep(0,a),1,-1,rep(0,end)))
    data$table
  } else if(x==r){
    data <- glmLRT(fit, contrast=c(rep(0,a),1,rep(0,end),-1))
    data$table
  } else {
    data <- glmLRT(fit, contrast=c(rep(0,a),1,rep(0,c),-1,rep(0,d)))
    data$table}
}
t <- proc.time()
v16.56 <- mclapply(1:40,runit,mc.cores = 40)
proc.time() - t

#17####
runit <- function(x){
  a <- 16
  end <- 38
  d <- end-x
  c <- (x-1)+1
  r <- end+1
  if(x == 1){
    data <- glmLRT(fit, contrast=c(rep(0,a),1,-1,rep(0,end)))
    data$table
  } else if(x==r){
    data <- glmLRT(fit, contrast=c(rep(0,a),1,rep(0,end),-1))
    data$table
  } else {
    data <- glmLRT(fit, contrast=c(rep(0,a),1,rep(0,c),-1,rep(0,d)))
    data$table}
}
t <- proc.time()
v17.56 <- mclapply(1:39,runit,mc.cores = 39)
proc.time() - t

#18####
runit <- function(x){
  a <- 17
  end <- 37
  d <- end-x
  c <- (x-1)+1
  r <- end+1
  if(x == 1){
    data <- glmLRT(fit, contrast=c(rep(0,a),1,-1,rep(0,end)))
    data$table
  } else if(x==r){
    data <- glmLRT(fit, contrast=c(rep(0,a),1,rep(0,end),-1))
    data$table
  } else {
    data <- glmLRT(fit, contrast=c(rep(0,a),1,rep(0,c),-1,rep(0,d)))
    data$table}
}
t <- proc.time()
v18.56 <- mclapply(1:38,runit,mc.cores = 38)
proc.time() - t

#19####
runit <- function(x){
  a <- 18
  end <- 36
  d <- end-x
  c <- (x-1)+1
  r <- end+1
  if(x == 1){
    data <- glmLRT(fit, contrast=c(rep(0,a),1,-1,rep(0,end)))
    data$table
  } else if(x==r){
    data <- glmLRT(fit, contrast=c(rep(0,a),1,rep(0,end),-1))
    data$table
  } else {
    data <- glmLRT(fit, contrast=c(rep(0,a),1,rep(0,c),-1,rep(0,d)))
    data$table}
}
t <- proc.time()
v19.56 <- mclapply(1:37,runit,mc.cores = 37)
proc.time() - t

#20####
runit <- function(x){
  a <- 19
  end <- 35
  d <- end-x
  c <- (x-1)+1
  r <- end+1
  if(x == 1){
    data <- glmLRT(fit, contrast=c(rep(0,a),1,-1,rep(0,end)))
    data$table
  } else if(x==r){
    data <- glmLRT(fit, contrast=c(rep(0,a),1,rep(0,end),-1))
    data$table
  } else {
    data <- glmLRT(fit, contrast=c(rep(0,a),1,rep(0,c),-1,rep(0,d)))
    data$table}
}
t <- proc.time()
v20.56 <- mclapply(1:36,runit,mc.cores = 36)
proc.time() - t

#21####
runit <- function(x){
  a <- 20
  end <- 34
  d <- end-x
  c <- (x-1)+1
  r <- end+1
  if(x == 1){
    data <- glmLRT(fit, contrast=c(rep(0,a),1,-1,rep(0,end)))
    data$table
  } else if(x==r){
    data <- glmLRT(fit, contrast=c(rep(0,a),1,rep(0,end),-1))
    data$table
  } else {
    data <- glmLRT(fit, contrast=c(rep(0,a),1,rep(0,c),-1,rep(0,d)))
    data$table}
}
t <- proc.time()
v21.56 <- mclapply(1:35,runit,mc.cores = 35)
proc.time() - t

#22####
runit <- function(x){
  a <- 21
  end <- 33
  d <- end-x
  c <- (x-1)+1
  r <- end+1
  if(x == 1){
    data <- glmLRT(fit, contrast=c(rep(0,a),1,-1,rep(0,end)))
    data$table
  } else if(x==r){
    data <- glmLRT(fit, contrast=c(rep(0,a),1,rep(0,end),-1))
    data$table
  } else {
    data <- glmLRT(fit, contrast=c(rep(0,a),1,rep(0,c),-1,rep(0,d)))
    data$table}
}
t <- proc.time()
v22.56 <- mclapply(1:34,runit,mc.cores = 34)
proc.time() - t

#23####
runit <- function(x){
  a <- 22
  end <- 32
  d <- end-x
  c <- (x-1)+1
  r <- end+1
  if(x == 1){
    data <- glmLRT(fit, contrast=c(rep(0,a),1,-1,rep(0,end)))
    data$table
  } else if(x==r){
    data <- glmLRT(fit, contrast=c(rep(0,a),1,rep(0,end),-1))
    data$table
  } else {
    data <- glmLRT(fit, contrast=c(rep(0,a),1,rep(0,c),-1,rep(0,d)))
    data$table}
}
t <- proc.time()
v23.56 <- mclapply(1:33,runit,mc.cores = 33)
proc.time() - t

#24####
runit <- function(x){
  a <- 23
  end <- 31
  d <- end-x
  c <- (x-1)+1
  r <- end+1
  if(x == 1){
    data <- glmLRT(fit, contrast=c(rep(0,a),1,-1,rep(0,end)))
    data$table
  } else if(x==r){
    data <- glmLRT(fit, contrast=c(rep(0,a),1,rep(0,end),-1))
    data$table
  } else {
    data <- glmLRT(fit, contrast=c(rep(0,a),1,rep(0,c),-1,rep(0,d)))
    data$table}
}
t <- proc.time()
v24.56 <- mclapply(1:32,runit,mc.cores = 32)
proc.time() - t

#25####
runit <- function(x){
  a <- 24
  end <- 30
  d <- end-x
  c <- (x-1)+1
  r <- end+1
  if(x == 1){
    data <- glmLRT(fit, contrast=c(rep(0,a),1,-1,rep(0,end)))
    data$table
  } else if(x==r){
    data <- glmLRT(fit, contrast=c(rep(0,a),1,rep(0,end),-1))
    data$table
  } else {
    data <- glmLRT(fit, contrast=c(rep(0,a),1,rep(0,c),-1,rep(0,d)))
    data$table}
}
t <- proc.time()
v25.56 <- mclapply(1:31,runit,mc.cores = 31)
proc.time() - t



#26####
runit <- function(x){
  a <- 25
  end <- 29
  d <- end-x
  c <- (x-1)+1
  r <- end+1
  if(x == 1){
    data <- glmLRT(fit, contrast=c(rep(0,a),1,-1,rep(0,end)))
    data$table
  } else if(x==r){
    data <- glmLRT(fit, contrast=c(rep(0,a),1,rep(0,end),-1))
    data$table
  } else {
    data <- glmLRT(fit, contrast=c(rep(0,a),1,rep(0,c),-1,rep(0,d)))
    data$table}
}
t <- proc.time()
v26.56 <- mclapply(1:30,runit,mc.cores = 30)
proc.time() - t

#27####
runit <- function(x){
  a <- 26
  end <- 28
  d <- end-x
  c <- (x-1)+1
  r <- end+1
  if(x == 1){
    data <- glmLRT(fit, contrast=c(rep(0,a),1,-1,rep(0,end)))
    data$table
  } else if(x==r){
    data <- glmLRT(fit, contrast=c(rep(0,a),1,rep(0,end),-1))
    data$table
  } else {
    data <- glmLRT(fit, contrast=c(rep(0,a),1,rep(0,c),-1,rep(0,d)))
    data$table}
}
t <- proc.time()
v27.56 <- mclapply(1:29,runit,mc.cores = 29)
proc.time() - t

#28####
runit <- function(x){
  a <- 27
  end <- 27
  d <- end-x
  c <- (x-1)+1
  r <- end+1
  if(x == 1){
    data <- glmLRT(fit, contrast=c(rep(0,a),1,-1,rep(0,end)))
    data$table
  } else if(x==r){
    data <- glmLRT(fit, contrast=c(rep(0,a),1,rep(0,end),-1))
    data$table
  } else {
    data <- glmLRT(fit, contrast=c(rep(0,a),1,rep(0,c),-1,rep(0,d)))
    data$table}
}
t <- proc.time()
v28.56 <- mclapply(1:28,runit,mc.cores = 28)
proc.time() - t

#29####
runit <- function(x){
  a <- 28
  end <- 26
  d <- end-x
  c <- (x-1)+1
  r <- end+1
  if(x == 1){
    data <- glmLRT(fit, contrast=c(rep(0,a),1,-1,rep(0,end)))
    data$table
  } else if(x==r){
    data <- glmLRT(fit, contrast=c(rep(0,a),1,rep(0,end),-1))
    data$table
  } else {
    data <- glmLRT(fit, contrast=c(rep(0,a),1,rep(0,c),-1,rep(0,d)))
    data$table}
}
t <- proc.time()
v29.56 <- mclapply(1:27,runit,mc.cores = 27)
proc.time() - t

#30####
runit <- function(x){
  a <- 29
  end <- 25
  d <- end-x
  c <- (x-1)+1
  r <- end+1
  if(x == 1){
    data <- glmLRT(fit, contrast=c(rep(0,a),1,-1,rep(0,end)))
    data$table
  } else if(x==r){
    data <- glmLRT(fit, contrast=c(rep(0,a),1,rep(0,end),-1))
    data$table
  } else {
    data <- glmLRT(fit, contrast=c(rep(0,a),1,rep(0,c),-1,rep(0,d)))
    data$table}
}
t <- proc.time()
v30.56 <- mclapply(1:26,runit,mc.cores = 26)
proc.time() - t

#31####
runit <- function(x){
  a <- 30
  end <- 24
  d <- end-x
  c <- (x-1)+1
  r <- end+1
  if(x == 1){
    data <- glmLRT(fit, contrast=c(rep(0,a),1,-1,rep(0,end)))
    data$table
  } else if(x==r){
    data <- glmLRT(fit, contrast=c(rep(0,a),1,rep(0,end),-1))
    data$table
  } else {
    data <- glmLRT(fit, contrast=c(rep(0,a),1,rep(0,c),-1,rep(0,d)))
    data$table}
}
t <- proc.time()
v31.56 <- mclapply(1:25,runit,mc.cores = 25)
proc.time() - t

#32####
runit <- function(x){
  a <- 31
  end <- 23
  d <- end-x
  c <- (x-1)+1
  r <- end+1
  if(x == 1){
    data <- glmLRT(fit, contrast=c(rep(0,a),1,-1,rep(0,end)))
    data$table
  } else if(x==r){
    data <- glmLRT(fit, contrast=c(rep(0,a),1,rep(0,end),-1))
    data$table
  } else {
    data <- glmLRT(fit, contrast=c(rep(0,a),1,rep(0,c),-1,rep(0,d)))
    data$table}
}
t <- proc.time()
v32.56 <- mclapply(1:24,runit,mc.cores = 24)
proc.time() - t

#33####
runit <- function(x){
  a <- 32
  end <- 22
  d <- end-x
  c <- (x-1)+1
  r <- end+1
  if(x == 1){
    data <- glmLRT(fit, contrast=c(rep(0,a),1,-1,rep(0,end)))
    data$table
  } else if(x==r){
    data <- glmLRT(fit, contrast=c(rep(0,a),1,rep(0,end),-1))
    data$table
  } else {
    data <- glmLRT(fit, contrast=c(rep(0,a),1,rep(0,c),-1,rep(0,d)))
    data$table}
}
t <- proc.time()
v33.56 <- mclapply(1:23,runit,mc.cores = 23)
proc.time() - t

#34####
runit <- function(x){
  a <- 33
  end <- 21
  d <- end-x
  c <- (x-1)+1
  r <- end+1
  if(x == 1){
    data <- glmLRT(fit, contrast=c(rep(0,a),1,-1,rep(0,end)))
    data$table
  } else if(x==r){
    data <- glmLRT(fit, contrast=c(rep(0,a),1,rep(0,end),-1))
    data$table
  } else {
    data <- glmLRT(fit, contrast=c(rep(0,a),1,rep(0,c),-1,rep(0,d)))
    data$table}
}
t <- proc.time()
v34.56 <- mclapply(1:22,runit,mc.cores = 22)
proc.time() - t

#35####
runit <- function(x){
  a <- 34
  end <- 20
  d <- end-x
  c <- (x-1)+1
  r <- end+1
  if(x == 1){
    data <- glmLRT(fit, contrast=c(rep(0,a),1,-1,rep(0,end)))
    data$table
  } else if(x==r){
    data <- glmLRT(fit, contrast=c(rep(0,a),1,rep(0,end),-1))
    data$table
  } else {
    data <- glmLRT(fit, contrast=c(rep(0,a),1,rep(0,c),-1,rep(0,d)))
    data$table}
}
t <- proc.time()
v35.56 <- mclapply(1:21,runit,mc.cores = 21)
proc.time() - t

#36####
runit <- function(x){
  a <- 35
  end <- 19
  d <- end-x
  c <- (x-1)+1
  r <- end+1
  if(x == 1){
    data <- glmLRT(fit, contrast=c(rep(0,a),1,-1,rep(0,end)))
    data$table
  } else if(x==r){
    data <- glmLRT(fit, contrast=c(rep(0,a),1,rep(0,end),-1))
    data$table
  } else {
    data <- glmLRT(fit, contrast=c(rep(0,a),1,rep(0,c),-1,rep(0,d)))
    data$table}
}
t <- proc.time()
v36.56 <- mclapply(1:20,runit,mc.cores = 20)
proc.time() - t

#37####
runit <- function(x){
  a <- 36
  end <- 18
  d <- end-x
  c <- (x-1)+1
  r <- end+1
  if(x == 1){
    data <- glmLRT(fit, contrast=c(rep(0,a),1,-1,rep(0,end)))
    data$table
  } else if(x==r){
    data <- glmLRT(fit, contrast=c(rep(0,a),1,rep(0,end),-1))
    data$table
  } else {
    data <- glmLRT(fit, contrast=c(rep(0,a),1,rep(0,c),-1,rep(0,d)))
    data$table}
}
t <- proc.time()
v37.56 <- mclapply(1:19,runit,mc.cores = 19)
proc.time() - t

#38####
runit <- function(x){
  a <- 37
  end <- 17
  d <- end-x
  c <- (x-1)+1
  r <- end+1
  if(x == 1){
    data <- glmLRT(fit, contrast=c(rep(0,a),1,-1,rep(0,end)))
    data$table
  } else if(x==r){
    data <- glmLRT(fit, contrast=c(rep(0,a),1,rep(0,end),-1))
    data$table
  } else {
    data <- glmLRT(fit, contrast=c(rep(0,a),1,rep(0,c),-1,rep(0,d)))
    data$table}
}
t <- proc.time()
v38.56 <- mclapply(1:18,runit,mc.cores = 18)
proc.time() - t

#39####
runit <- function(x){
  a <- 38
  end <- 16
  d <- end-x
  c <- (x-1)+1
  r <- end+1
  if(x == 1){
    data <- glmLRT(fit, contrast=c(rep(0,a),1,-1,rep(0,end)))
    data$table
  } else if(x==r){
    data <- glmLRT(fit, contrast=c(rep(0,a),1,rep(0,end),-1))
    data$table
  } else {
    data <- glmLRT(fit, contrast=c(rep(0,a),1,rep(0,c),-1,rep(0,d)))
    data$table}
}
t <- proc.time()
v39.56 <- mclapply(1:17,runit,mc.cores = 17)
proc.time() - t

#40####
runit <- function(x){
  a <- 39
  end <- 15
  d <- end-x
  c <- (x-1)+1
  r <- end+1
  if(x == 1){
    data <- glmLRT(fit, contrast=c(rep(0,a),1,-1,rep(0,end)))
    data$table
  } else if(x==r){
    data <- glmLRT(fit, contrast=c(rep(0,a),1,rep(0,end),-1))
    data$table
  } else {
    data <- glmLRT(fit, contrast=c(rep(0,a),1,rep(0,c),-1,rep(0,d)))
    data$table}
}
t <- proc.time()
v40.56 <- mclapply(1:16,runit,mc.cores = 16)
proc.time() - t

#41####
runit <- function(x){
  a <- 40
  end <- 14
  d <- end-x
  c <- (x-1)+1
  r <- end+1
  if(x == 1){
    data <- glmLRT(fit, contrast=c(rep(0,a),1,-1,rep(0,end)))
    data$table
  } else if(x==r){
    data <- glmLRT(fit, contrast=c(rep(0,a),1,rep(0,end),-1))
    data$table
  } else {
    data <- glmLRT(fit, contrast=c(rep(0,a),1,rep(0,c),-1,rep(0,d)))
    data$table}
}
t <- proc.time()
v41.56 <- mclapply(1:15,runit,mc.cores = 15)
proc.time() - t

#42####
runit <- function(x){
  a <- 41
  end <- 13
  d <- end-x
  c <- (x-1)+1
  r <- end+1
  if(x == 1){
    data <- glmLRT(fit, contrast=c(rep(0,a),1,-1,rep(0,end)))
    data$table
  } else if(x==r){
    data <- glmLRT(fit, contrast=c(rep(0,a),1,rep(0,end),-1))
    data$table
  } else {
    data <- glmLRT(fit, contrast=c(rep(0,a),1,rep(0,c),-1,rep(0,d)))
    data$table}
}
t <- proc.time()
v42.56 <- mclapply(1:14,runit,mc.cores = 14)
proc.time() - t

#43####
runit <- function(x){
  a <- 42
  end <- 12
  d <- end-x
  c <- (x-1)+1
  r <- end+1
  if(x == 1){
    data <- glmLRT(fit, contrast=c(rep(0,a),1,-1,rep(0,end)))
    data$table
  } else if(x==r){
    data <- glmLRT(fit, contrast=c(rep(0,a),1,rep(0,end),-1))
    data$table
  } else {
    data <- glmLRT(fit, contrast=c(rep(0,a),1,rep(0,c),-1,rep(0,d)))
    data$table}
}
t <- proc.time()
v43.56 <- mclapply(1:13,runit,mc.cores = 13)
proc.time() - t

#44####
runit <- function(x){
  a <- 43
  end <- 11
  d <- end-x
  c <- (x-1)+1
  r <- end+1
  if(x == 1){
    data <- glmLRT(fit, contrast=c(rep(0,a),1,-1,rep(0,end)))
    data$table
  } else if(x==r){
    data <- glmLRT(fit, contrast=c(rep(0,a),1,rep(0,end),-1))
    data$table
  } else {
    data <- glmLRT(fit, contrast=c(rep(0,a),1,rep(0,c),-1,rep(0,d)))
    data$table}
}
t <- proc.time()
v44.56 <- mclapply(1:12,runit,mc.cores = 12)
proc.time() - t

#45####
runit <- function(x){
  a <- 44
  end <- 10
  d <- end-x
  c <- (x-1)+1
  r <- end+1
  if(x == 1){
    data <- glmLRT(fit, contrast=c(rep(0,a),1,-1,rep(0,end)))
    data$table
  } else if(x==r){
    data <- glmLRT(fit, contrast=c(rep(0,a),1,rep(0,end),-1))
    data$table
  } else {
    data <- glmLRT(fit, contrast=c(rep(0,a),1,rep(0,c),-1,rep(0,d)))
    data$table}
}
t <- proc.time()
v45.56 <- mclapply(1:11,runit,mc.cores = 11)
proc.time() - t

#46####
runit <- function(x){
  a <- 45
  end <- 9
  d <- end-x
  c <- (x-1)+1
  r <- end+1
  if(x == 1){
    data <- glmLRT(fit, contrast=c(rep(0,a),1,-1,rep(0,end)))
    data$table
  } else if(x==r){
    data <- glmLRT(fit, contrast=c(rep(0,a),1,rep(0,end),-1))
    data$table
  } else {
    data <- glmLRT(fit, contrast=c(rep(0,a),1,rep(0,c),-1,rep(0,d)))
    data$table}
}
t <- proc.time()
v46.56 <- mclapply(1:10,runit,mc.cores = 10)
proc.time() - t

#47####
runit <- function(x){
  a <- 46
  end <- 8
  d <- end-x
  c <- (x-1)+1
  r <- end+1
  if(x == 1){
    data <- glmLRT(fit, contrast=c(rep(0,a),1,-1,rep(0,end)))
    data$table
  } else if(x==r){
    data <- glmLRT(fit, contrast=c(rep(0,a),1,rep(0,end),-1))
    data$table
  } else {
    data <- glmLRT(fit, contrast=c(rep(0,a),1,rep(0,c),-1,rep(0,d)))
    data$table}
}
t <- proc.time()
v47.56 <- mclapply(1:9,runit,mc.cores = 9)
proc.time() - t

#48####
runit <- function(x){
  a <- 47
  end <- 7
  d <- end-x
  c <- (x-1)+1
  r <- end+1
  if(x == 1){
    data <- glmLRT(fit, contrast=c(rep(0,a),1,-1,rep(0,end)))
    data$table
  } else if(x==r){
    data <- glmLRT(fit, contrast=c(rep(0,a),1,rep(0,end),-1))
    data$table
  } else {
    data <- glmLRT(fit, contrast=c(rep(0,a),1,rep(0,c),-1,rep(0,d)))
    data$table}
}
t <- proc.time()
v48.56 <- mclapply(1:8,runit,mc.cores = 8)
proc.time() - t

#49####
runit <- function(x){
  a <- 48
  end <- 6
  d <- end-x
  c <- (x-1)+1
  r <- end+1
  if(x == 1){
    data <- glmLRT(fit, contrast=c(rep(0,a),1,-1,rep(0,end)))
    data$table
  } else if(x==r){
    data <- glmLRT(fit, contrast=c(rep(0,a),1,rep(0,end),-1))
    data$table
  } else {
    data <- glmLRT(fit, contrast=c(rep(0,a),1,rep(0,c),-1,rep(0,d)))
    data$table}
}
t <- proc.time()
v49.56 <- mclapply(1:7,runit,mc.cores = 7)
proc.time() - t

#50####
runit <- function(x){
  a <- 49
  end <- 5
  d <- end-x
  c <- (x-1)+1
  r <- end+1
  if(x == 1){
    data <- glmLRT(fit, contrast=c(rep(0,a),1,-1,rep(0,end)))
    data$table
  } else if(x==r){
    data <- glmLRT(fit, contrast=c(rep(0,a),1,rep(0,end),-1))
    data$table
  } else {
    data <- glmLRT(fit, contrast=c(rep(0,a),1,rep(0,c),-1,rep(0,d)))
    data$table}
}
t <- proc.time()
v50.56 <- mclapply(1:6,runit,mc.cores = 6)
proc.time() - t

#51####
runit <- function(x){
  a <- 50
  end <- 4
  d <- end-x
  c <- (x-1)+1
  r <- end+1
  if(x == 1){
    data <- glmLRT(fit, contrast=c(rep(0,a),1,-1,rep(0,end)))
    data$table
  } else if(x==r){
    data <- glmLRT(fit, contrast=c(rep(0,a),1,rep(0,end),-1))
    data$table
  } else {
    data <- glmLRT(fit, contrast=c(rep(0,a),1,rep(0,c),-1,rep(0,d)))
    data$table}
}
t <- proc.time()
v51.56 <- mclapply(1:5,runit,mc.cores = 5)
proc.time() - t

#52####
runit <- function(x){
  a <- 51
  end <- 3
  d <- end-x
  c <- (x-1)+1
  r <- end+1
  if(x == 1){
    data <- glmLRT(fit, contrast=c(rep(0,a),1,-1,rep(0,end)))
    data$table
  } else if(x==r){
    data <- glmLRT(fit, contrast=c(rep(0,a),1,rep(0,end),-1))
    data$table
  } else {
    data <- glmLRT(fit, contrast=c(rep(0,a),1,rep(0,c),-1,rep(0,d)))
    data$table}
}
t <- proc.time()
v52.56 <- mclapply(1:4,runit,mc.cores = 4)
proc.time() - t

#53####
runit <- function(x){
  a <- 52
  end <- 2
  d <- end-x
  c <- (x-1)+1
  r <- end+1
  if(x == 1){
    data <- glmLRT(fit, contrast=c(rep(0,a),1,-1,rep(0,end)))
    data$table
  } else if(x==r){
    data <- glmLRT(fit, contrast=c(rep(0,a),1,rep(0,end),-1))
    data$table
  } else {
    data <- glmLRT(fit, contrast=c(rep(0,a),1,rep(0,c),-1,rep(0,d)))
    data$table}
}
t <- proc.time()
v53.56 <- mclapply(1:3,runit,mc.cores = 3)
proc.time() - t

#54####
runit <- function(x){
  a <- 53
  end <- 1
  d <- end-x
  c <- (x-1)+1
  r <- end+1
  if(x == 1){
    data <- glmLRT(fit, contrast=c(rep(0,a),1,-1,rep(0,end)))
    data$table
  } else if(x==r){
    data <- glmLRT(fit, contrast=c(rep(0,a),1,rep(0,end),-1))
    data$table
  } else {
    data <- glmLRT(fit, contrast=c(rep(0,a),1,rep(0,c),-1,rep(0,d)))
    data$table}
}
t <- proc.time()
v54.56 <- mclapply(1:2,runit,mc.cores = 2)
proc.time() - t

#55####
runit <- function(x){
  a <- 54
  end <- 0
  d <- end-x
  c <- (x-1)+1
  r <- end+1
  if(x==r){
    data <- glmLRT(fit, contrast=c(rep(0,a),1,-1))
    data$table
  } }
t <- proc.time()
v55.56 <- mclapply(1,runit,mc.cores = 1)
proc.time() - t

#Remove unneccesary objects and save up to this point####
rm(cor.192,cor.192.2,cor.comp,cor.comp.ht,cor.comp.ht2,esults.mtx3,l,log.DE,log.DE2,mat,matr,matr.chi,matr.exp,
   matr.exp.bot,matr.exp.top,matr.exp.top.bot,matt,mean.cor,mean.out,mean.out.1,output,phen,phen.tx,phenos.1,phenos.2,
   phenos.test,phenos1,r.comp.ht,r.comp.ht2,r.DE,result.df3,SNP.224,SNP.240,SNP.cor,SNP.OP.cor,SNP.OP.mean,table,tip.snp.192,
   tip.snp.192.2,tmp,v21.55,count,cross,data,comps.ht,comps,comps.vol,com,col.2,col,chi,c.list,beg,bad.transcripts,ans1,
   ans,all,data1,data2,data3,df,each,exp,exp.bot,exp.top,exp.top.v.bot,expec,f,fam,family,fdr.,fdr.21.55,full.comparisons,
   full.parents,g,good.transcripts,goodvbad.transcripts,grp,grps,i,j,k.test,last,lom,num.tests,obs,obs.bot,obs.top,
   obs.top.v.bot,observ,observ2,other.21,other.fam,par1,par2,pvalues,pvalues1,pvalues2,q,r,ro,row,row.2,rows,
   rows.top,set,singles,some,States,subs,subs.top,t,t.trans,test,test.group,test.new,test.parents,test2,the.end,the.par2,
   the.row,top,top.cors,top.v.worst.transcripts,total,training.de,try,try.1,try.2,tx.fam,u.all,u.try1,u.try2,uni.group,x,y,yo)
rm(matr.obs,matr.obs.bot,matr.obs.top,matr.obs.top.bot,the.one,v3.56,v4.56,v5.56,v6.56,v7.56,v8.56,v9.56)

save.image("~/Desktop/identifying.top.DE.RDA.RData")

#Combine all DE comparisons into one list####
full.comparisons <- list( v01.56,v02.56,v03.56,
                          v04.56 ,v05.56 ,v06.56 , v07.56, v08.56 , v09.56 , v10.56 , v11.56 , v12.56 , v13.56 , v14.56 , v15.56 ,
                          v16.56, v17.56 , v18.56 ,v19.56 , v20.56 ,v21.56 ,v22.56 , v23.56,v24.56 ,
                          v25.56 ,v26.56 ,v27.56 ,v28.56 ,v29.56 ,v30.56,v31.56,v32.56,v33.56,v34.56,v35.56 ,
                          v36.56 ,v37.56 ,v38.56 ,v39.56 , v40.56 ,v41.56 , v42.56 ,v43.56 ,v44.56 ,v45.56 ,
                          v46.56 ,v47.56 ,v48.56 ,v49.56 ,v50.56 ,v51.56 ,v52.56 ,v53.56 ,v54.56,v55.56)

#Calculate the FDR for all 1540 comparisons####
full.comparisons.fdr <-vector("list",55)
for(i in 1:55){
  a <- length(full.comparisons[[i]])
  for(each in 1:a){
    full.comparisons.fdr[[i]][[each]] <- lfdr(p=full.comparisons[[i]][[each]]$PValue)
  }
  print(i)
}

#Subset the rownames for each comparison that meet the .05 threshold
full.comparisons.fdr.05 <- vector("list",55)
for(i in 1:55){
  a <- length(full.comparisons.fdr[[i]])
  for(each in 1:a){
    full.comparisons.fdr.05[[i]][[each]] <-  rownames(full.comparisons[[i]][[each]])[which(full.comparisons.fdr[[i]][[each]] < .05)] 
  }
  print(i)
}


#Identify the total number of signifcant variables which cross that threshold for all 1540 comparisons
total.sig.fdr.05 <- vector("list",55)
for(i in 1:55){
  a <- length(full.comparisons.fdr[[i]])
  for(each in 1:a){
    total.sig.fdr.05[[i]][[each]] <-  length(full.comparisons.fdr.05[[i]][[each]])
  }
  print(i)
}
hist(unlist(total.sig.fdr.05))
#Which comparison has the highest number of DE transcripts
which.max(unlist(total.sig.fdr.05))
#comparison 481 out of 1540

#What is the total number of DE transcripts for that comparison
unlist(total.sig.fdr.05)[481]
#8259

#To make this comparison easier to identify, create a 55x56 matrix which holds number of DE transcripts for all comparisons
total.sig.fdr.05.matrix <- matrix(NA,55,56)
beg <- 1
end <- 55
match <- 2
all <- 55
for(i in 1:55){
  c <- beg
  d <- end
  total.sig.fdr.05.matrix[i,match:56] <- unlist(total.sig.fdr.05)[c:d]
  print(i)
  beg <- d + 1
  end <- d+(all-i)
  match <- match+1
}

##The last two numbers of all comparisons are same, see what went wrong####
data <- glmLRT(fit, contrast=c(0,1,rep(0,2),-1,rep(0,51)))
test <- lfdr(p=data$table$PValue)
length(which(test < .05))

#Issue is that the 2nd comparison for all families was skipped, 
#first subset #1 from all previous comparisons
corrected.full.comparisons <- vector("list",55)
for(each in 1:55){
  corrected.full.comparisons[[each]][[1]] <- full.comparisons[[each]][[1]]
}

#Now for the first set of comparisons subset 2:54 and put into corrected list of 3:55
g<-3
for(each in 2:54){
  corrected.full.comparisons[[1]][[g]] <- full.comparisons[[1]][[each]]
  g <- g+1
}

#Now for comparisons 2:53 subset 
for(i in 2:53){
  a <- length(full.comparisons[[i]]) - 1
  g <- 3
  for(each in 2:a){
    corrected.full.comparisons[[i]][[g]] <- full.comparisons[[i]][[each]]
    g <- g+1
  }
}

corrected.full.comparisons[[54]][[2]] <- full.comparisons[[54]][[2]]

#Now re-do 2nd contrast for all families
data <- glmLRT(fit, contrast=c(1,rep(0,1),-1,rep(0,53)))
corrected.full.comparisons[[1]][[2]] <- data$table
test <- lfdr(p=data$table$PValue)
length(which(test < .05))

data <- glmLRT(fit, contrast=c(rep(0,1),1,rep(0,1),-1,rep(0,52)))
corrected.full.comparisons[[2]][[2]] <- data$table
test <- lfdr(p=data$table$PValue)
length(which(test < .05))

data <- glmLRT(fit, contrast=c(rep(0,2),1,rep(0,1),-1,rep(0,51)))
corrected.full.comparisons[[3]][[2]] <- data$table
test <- lfdr(p=data$table$PValue)
length(which(test < .05))

a <- 3; b <- 50
for(each in 4:53){
  data <- glmLRT(fit, contrast=c(rep(0,a),1,rep(0,1),-1,rep(0,b)))
  corrected.full.comparisons[[each]][[2]] <- data$table
  a <- a+1
  b <- b-1
  print(each)}

##Now that all comparisons are corrected calulate fdr for all####
corrected.full.comparisons.fdr <-vector("list",55)
for(i in 1:55){
  a <- length(corrected.full.comparisons[[i]])
  for(each in 1:a){
    corrected.full.comparisons.fdr[[i]][[each]] <- lfdr(p=corrected.full.comparisons[[i]][[each]]$PValue)
  }
  print(i)
}
rm(cor.192,cor.comp.number,cor.comp.vol,cor.fam,corrected.total.sig.fdr.1,corrected.total.sig.fdr.01,corrected.total.sig.fdr.05.matrix,
   fam.logcpm2,fam.logcpm2.7.10,fam.phenos,log.DE,output,phenos.fam,phenos.test,r.comp.vol,result.df3,results.mtx3,scaled.fam.logcpm2,
   test.logcpm2,tmp,A,a,all,all.cors,all.fam,all.fam.de,all.transcripts,b,bad.transcripts,beg,c,c.fam.cpm,count,cv.folds,d,data,each,
   end,fam,fam.1.de,fam.2.de,fam.3.de,fam.4.de,fam.5.de,fam.6.de,fam.7.de,fam.8.de,fam.cor,fam.cor.01,fam.de,fam.de.count,family,fdr,fit,
   full.comparisons,full.comparisons.fdr.05,full.comparisons.fdr,full.comparisons.fdr.01,full.comparisons.fdr.1,g, good.transcripts,good.v.bad.txpts,
   groups,grp,i,keep,lm,match,med.q1,med.q2,med.q3,med.q4,num.tests,other.fam,out,out2,p.7.10,pc.fam.cpm,pred,qu,s.7.10,set,test,test.group,test.parents,the.par,
   total,training.de,true,uniq.fam,unique.all.fam.de,unique.fam.de.count,v01.56,v02.56,v03.56,v04.56,v05.56,v06.56,v07.56,v08.56,v09.56,v10.56,
   v11.56,v12.56,v13.56,v14.56,v15.56,v16.56,v17.56,v18.56,v19.56,v20.56,v21.56,v22.56,v23.56,v24.56,v25.56,v26.56,v27.56,v28.56,v29.56,v30.56,
   v31.56,v32.56,v33.56,v34.56,v35.56,v36.56,v37.56,v38.56,v39.56,v40.56,v41.56,v42.56,v43.56,v44.56,v45.56,v46.56,v47.56,v48.56,v49.56,v50.56,
   v51.56,v52.56,v53.56,v54.56,v55.56,vol.grps,x,yo,corrected.total.sig.fdr.01.matrix,corrected.total.sig.fdr.1.matrix,fam.cor.1)
rm(corrected.full.comparisons.fdr.01,corrected.full.comparisons.fdr.1,corrected.total.sig.fdr.05,corrected.full.comparisons.fdr.05)
#save.image("~/Desktop/corrected.identifying.top.DE.RDA.pt1.RData"
