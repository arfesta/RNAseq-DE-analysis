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

##Now that the calculate for fdr is correct subset those less than .01 ####
corrected.full.comparisons.fdr.01 <- vector("list",55)
for(i in 1:55){
  a <- length(corrected.full.comparisons.fdr[[i]])
  for(each in 1:a){
    corrected.full.comparisons.fdr.01[[i]][[each]] <-  rownames(corrected.full.comparisons[[i]][[each]])[which(corrected.full.comparisons.fdr[[i]][[each]] < .01)] 
  }
  print(i)
}

corrected.total.sig.fdr.01 <- vector("list",55)
for(i in 1:55){
  a <- length(corrected.full.comparisons.fdr[[i]])
  for(each in 1:a){
    corrected.total.sig.fdr.01[[i]][[each]] <-  length(corrected.full.comparisons.fdr.01[[i]][[each]])
  }
  print(i)
}
hist(unlist(corrected.total.sig.fdr.01))
summary(unlist(corrected.total.sig.fdr.01))
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
 #1.0   134.8   329.5   561.1   751.2  5727.0 
which.min(unlist(corrected.total.sig.fdr.01))

corrected.total.sig.fdr.01.matrix <- matrix(NA,55,56)
beg <- 1
end <- 55
match <- 2
all <- 55
for(i in 1:55){
  c <- beg
  d <- end
  corrected.total.sig.fdr.01.matrix[i,match:56] <- unlist(corrected.total.sig.fdr.01)[c:d]
  print(i)
  beg <- d + 1
  end <- d+(all-i)
  match <- match+1
}

##Now that the calculate for fdr is correct subset those less than .05 ####
corrected.full.comparisons.fdr.05 <- vector("list",55)
for(i in 1:55){
  a <- length(corrected.full.comparisons.fdr[[i]])
  for(each in 1:a){
    corrected.full.comparisons.fdr.05[[i]][[each]] <-  rownames(corrected.full.comparisons[[i]][[each]])[which(corrected.full.comparisons.fdr[[i]][[each]] < .05)] 
  }
  print(i)
}

corrected.total.sig.fdr.05 <- vector("list",55)
for(i in 1:55){
  a <- length(corrected.full.comparisons.fdr[[i]])
  for(each in 1:a){
    corrected.total.sig.fdr.05[[i]][[each]] <-  length(corrected.full.comparisons.fdr.05[[i]][[each]])
  }
  print(i)
}
hist(unlist(corrected.total.sig.fdr.05))
summary(unlist(corrected.total.sig.fdr.05))
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#   7.0   324.5   746.5  1088.0  1522.0  8259.0 

corrected.total.sig.fdr.05.matrix <- matrix(NA,55,56)
beg <- 1
end <- 55
match <- 2
all <- 55
for(i in 1:55){
  c <- beg
  d <- end
  corrected.total.sig.fdr.05.matrix[i,match:56] <- unlist(corrected.total.sig.fdr.05)[c:d]
  print(i)
  beg <- d + 1
  end <- d+(all-i)
  match <- match+1
}

#identical(corrected.total.sig.fdr.05.matrix[1,4:56],total.sig.fdr.05.matrix[1,3:55])
#TRUE

#since true we can remove the uncorrected comparisons
#rm(full.comparisons,full.comparisons.fdr,full.comparisons.fdr.05,total.sig.fdr.05,total.sig.fdr.05.matrix)

##Now that the calculate for fdr is correct subset those less than .1 ####
corrected.full.comparisons.fdr.1 <- vector("list",55)
for(i in 1:55){
  a <- length(corrected.full.comparisons.fdr[[i]])
  for(each in 1:a){
    corrected.full.comparisons.fdr.1[[i]][[each]] <-  rownames(corrected.full.comparisons[[i]][[each]])[which(corrected.full.comparisons.fdr[[i]][[each]] < .1)] 
  }
  print(i)
}

corrected.total.sig.fdr.1 <- vector("list",55)
for(i in 1:55){
  a <- length(corrected.full.comparisons.fdr[[i]])
  for(each in 1:a){
    corrected.total.sig.fdr.1[[i]][[each]] <-  length(corrected.full.comparisons.fdr.1[[i]][[each]])
  }
  print(i)
}
hist(unlist(corrected.total.sig.fdr.1))
summary(unlist(corrected.total.sig.fdr.1))
#     Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#     17.0   505.2  1099.0  1512.0  2127.0  9780.0 

corrected.total.sig.fdr.1.matrix <- matrix(NA,55,56)
beg <- 1
end <- 55
match <- 2
all <- 55
for(i in 1:55){
  c <- beg
  d <- end
  corrected.total.sig.fdr.1.matrix[i,match:56] <- unlist(corrected.total.sig.fdr.1)[c:d]
  print(i)
  beg <- d + 1
  end <- d+(all-i)
  match <- match+1
}

#Create comparison number 1-1540 by DE pair being compared####
cor.comp.number <- matrix(NA,55,56)
a <- 1
b <- 55
c <- 2
g <- 53
for(i in 1:55){
  cor.comp.number[i,c:56] <- a:b
  a <- b+1
  b <- a+g
  c <- c+1
  g <- g-1
  print(i)
}

#get family mean for counts and phenos####
phenos.fam <- aggregate(x=phenos$vol,by=list(phenos$group),mean)
rownames(phenos.fam) <- 1:56
fam.logcpm2 <- aggregate(t(logcpm2),by=list((phenos$group)),mean)
fam.logcpm2 <- fam.logcpm2[,-1]
scaled.fam.logcpm2 <- scale(fam.logcpm2)
identical(as.numeric(colnames(logcpm2)),phenos$animal)

#Create unrelated groups for CV####
cv.folds <- vector("list")
cv.folds[[1]] <- c(1,11,26,35,36,46,54)
cv.folds[[2]] <- c(2,9,12,20,27,37,47)
cv.folds[[3]] <- c(3,13,14,24,28,38,48)
cv.folds[[4]] <- c(4,15,29,39,43,45,49)
cv.folds[[5]] <- c(5,8,16,30,40,50,56)
cv.folds[[6]] <- c(6,17,22,31,33,41,51)
cv.folds[[7]] <- c(7,18,23,32,42,52,55)
cv.folds[[8]] <- c(10,19,21,25,34,44,53)
sort(unlist(cv.folds))
identical(as.numeric(colnames(logcpm2)),phenos$animal)

#use .01 threshold for prediction####
fam.cor.01 <- vector()
for(fam in 1:55){
num.tests <- length(corrected.full.comparisons[[fam]])
for(test.group in 1:num.tests){
test <- unique(c((full.comparisons.fdr.01[[fam]][[test.group]])))
log.DE <- scaled.fam.logcpm2[,match(test,colnames(scaled.fam.logcpm2))]
cor.fam <- cor(t(log.DE))
colnames(cor.fam) <- 1:56; rownames(cor.fam) <- 1:56
all.fam <- 1:56
pred <- vector();true <- vector(); all.cors <- vector()
for(i in 1:8){
  set <- cv.folds[[i]]
tmp <- okriging(idtest=set,idtrain = all.fam[-set],
                corlist=list(cor.fam),H2vec=c(.99),pheno=phenos.fam,phenoname="x")
pred <- c(pred,tmp[,2])
true <- c(true,tmp[,3])}
fam.cor.01 <- c(fam.cor.1,cor(pred,true))}
print(fam)}

hist(fam.cor.01)

sort(fam.cor,decreasing = T)[1:10]
which(fam.cor > .6)
#17vs51 alabama vs. florida
#7vs10 nc vs. sc
#32x53 FL vs. FL

#Look at top DE pair from .01####
test <- unique(c((corrected.full.comparisons.fdr.05[[7]][[3]])))
log.DE <- fam.logcpm2[,match(test,colnames(fam.logcpm2))]
cor.fam <- cor(t(log.DE))
colnames(cor.fam) <- 1:56; rownames(cor.fam) <- 1:56
all.fam <- 1:56
pred <- vector();true <- vector(); all.cors <- vector()
for(i in 1:8){
  set <- cv.folds[[i]]
  tmp <- okriging(idtest=set,idtrain = all.fam[-set],
                  corlist=list(cor.fam),H2vec=c(.99),pheno=phenos.fam,phenoname="x")
  print(cor(tmp[,2],tmp[,3]))
  pred <- c(pred,tmp[,2])
  true <- c(true,tmp[,3])}

plot(pred,true)

##Lets see how many times transcripts shows up for a family####
all.transcripts <- unlist(corrected.full.comparisons.fdr.05)
all.transcripts <- table(all.transcripts)
str(all.transcripts)
all.txpts <- matrix(NA,27949,1)
rownames(all.txpts) <- names(all.transcripts)
all.txpts[,1] <- all.transcripts
boxplot(all.txpts[,1])
test <- unique(c((corrected.full.comparisons.fdr.05[[7]][[3]])),corrected.full.comparisons.fdr.05[[32]][[20]])
hist(all.txpts[which(rownames(all.txpts) %in% test),1])

test <- unique(c((corrected.full.comparisons.fdr.05[[32]][[20]])))
hist(all.txpts[which(rownames(all.txpts) %in% test),1])
length(which(unique(c((corrected.full.comparisons.fdr.05[[32]][[20]]))) %in% unique(c((corrected.full.comparisons.fdr.05[[7]][[3]])))))

test <- rownames(all.txpts)[which(all.txpts[,1] > 450)]

log.DE <- fam.logcpm2[,match(test,colnames(fam.logcpm2))]
cor.fam <- cor(t(log.DE))
colnames(cor.fam) <- 1:56; rownames(cor.fam) <- 1:56
all.fam <- 1:56
pred <- vector();true <- vector(); all.cors <- vector()
for(i in 1:8){
  set <- cv.folds[[i]]
  tmp <- okriging(idtest=set,idtrain = all.fam[-set],
                  corlist=list(cor.fam),H2vec=c(.99),pheno=phenos.fam,phenoname="x")
  print(cor(tmp[,2],tmp[,3]))
  pred <- c(pred,tmp[,2])
  true <- c(true,tmp[,3])}

plot(pred,true)
cor(pred,true)^2

###otherstuff####
all.fam.de <- vector("list",56)
for(fam in 1:55){
  d=fam
fam.de <- vector()
for(i in 1:fam){
  d <- d-1
  if(d > 0){
    fam.de <- c(fam.de,unlist(corrected.full.comparisons.fdr.05[[i]][[d]]))
  } else {
    fam.de <- c(fam.de,unlist(corrected.full.comparisons.fdr.05[[fam]]))
  }
}
all.fam.de[[fam]] <- fam.de}

fam=56
d=56
fam.de <- vector()
for(i in 1:55){
  d <- d-1
  if(d > 0){
    fam.de <- c(fam.de,unlist(corrected.full.comparisons.fdr.05[[i]][[d]]))}}

all.fam.de[[56]] <- fam.de

#all.fam.de contains the list of DE transcripts called for each family####

#lets see how many DE per family
fam.de.count <- vector()
for(i in 1:56){
  fam.de.count <- c(fam.de.count,length(all.fam.de[[i]]))
}

summary(fam.de.count)
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#3747   37390   54480   59860   77120  185500 
boxplot(fam.de.count)
which(fam.de.count > 185000)
#family 33 has most DE transcripts across 55 comparisons
#family 1 has least DE transcripts across 55 comparisons



unique.all.fam.de <- vector("list",56)
for(fam in 1:55){
  d=fam
  fam.de <- vector()
  for(i in 1:fam){
    d <- d-1
    if(d > 0){
      fam.de <- unique(c(fam.de,unlist(corrected.full.comparisons.fdr.05[[i]][[d]])))
    } else {
      fam.de <- unique(c(fam.de,unlist(corrected.full.comparisons.fdr.05[[fam]])))
    }
  }
  unique.all.fam.de[[fam]] <- fam.de}

fam=56
d=56
fam.de <- vector()
for(i in 1:55){
  d <- d-1
  if(d > 0){
    fam.de <- unique(c(fam.de,unlist(corrected.full.comparisons.fdr.05[[i]][[d]])))}}

unique.all.fam.de[[56]] <- fam.de

#lets see how many DE per family
unique.fam.de.count <- vector()
for(i in 1:56){
  unique.fam.de.count <- c(unique.fam.de.count,length(unique.all.fam.de[[i]]))
}

summary(unique.fam.de.count)
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#1600   11610   13480   13620   16560   22220 
boxplot(unique.fam.de.count)
which.max(unique.fam.de.count)
#family 7 has most unique DE transcripts across 55 comparisons
#family 1 has least unique DE transcripts across 55 comparisons

##
for(each in 1:56){
test <- unique.all.fam.de[[each]]
log.DE <- fam.logcpm2[,match(test,colnames(fam.logcpm2))]
cor.fam <- cor(t(log.DE))
colnames(cor.fam) <- 1:56; rownames(cor.fam) <- 1:56
all.fam <- 1:56
pred <- vector();true <- vector(); all.cors <- vector()
for(i in 1:8){
  set <- cv.folds[[i]]
  tmp <- okriging(idtest=set,idtrain = all.fam[-set],
                  corlist=list(cor.fam),H2vec=c(.99),pheno=phenos.fam,phenoname="x")
  #print(cor(tmp[,2],tmp[,3]))
  pred <- c(pred,tmp[,2])
  true <- c(true,tmp[,3])}
print(cor(pred,true))
plot(pred,true)}
#all roughly r=.3



bad.transcripts <- unique(c(corrected.full.comparisons.fdr.05[[1]][[1]],corrected.full.comparisons.fdr.05[[1]][[2]],corrected.full.comparisons.fdr.05[[1]][[3]],
                            corrected.full.comparisons.fdr.05[[2]][[1]],corrected.full.comparisons.fdr.05[[2]][[2]],
                            corrected.full.comparisons.fdr.05[[3]][[1]]))


good.transcripts <- unique(c(corrected.full.comparisons.fdr.05[[55]][[1]],
                             corrected.full.comparisons.fdr.05[[54]][[2]],corrected.full.comparisons.fdr.05[[54]][[1]],
                             corrected.full.comparisons.fdr.05[[53]][[3]],corrected.full.comparisons.fdr.05[[53]][[2]],
                             corrected.full.comparisons.fdr.05[[53]][[1]]))


good.v.bad.txpts <- unique(c(corrected.full.comparisons.fdr.05[[1]][[55]],corrected.full.comparisons.fdr.05[[1]][[54]],corrected.full.comparisons.fdr.05[[1]][[53]],
                             corrected.full.comparisons.fdr.05[[2]][[54]],corrected.full.comparisons.fdr.05[[2]][[53]],corrected.full.comparisons.fdr.05[[2]][[52]],
                             corrected.full.comparisons.fdr.05[[3]][[53]],corrected.full.comparisons.fdr.05[[3]][[52]],corrected.full.comparisons.fdr.05[[3]][[51]]))


length(good.transcripts)#5399
length(bad.transcripts) #5355
length(good.v.bad.txpts) #4779

length(which(good.transcripts %in% bad.transcripts))#1146



test <- good.v.bad.txpts[!(good.v.bad.txpts %in% unique(c(good.transcripts,bad.transcripts)))]
log.DE <- fam.logcpm2[,match(test,colnames(fam.logcpm2))]
cor.fam <- cor(t(log.DE))
colnames(cor.fam) <- 1:56; rownames(cor.fam) <- 1:56
all.fam <- 1:56
pred <- vector();true <- vector(); all.cors <- vector()
for(i in 1:8){
  set <- cv.folds[[i]]
  tmp <- okriging(idtest=set,idtrain = all.fam[-set],
                  corlist=list(cor.fam),H2vec=c(.99),pheno=phenos.fam,phenoname="x")
  #print(cor(tmp[,2],tmp[,3]))
  pred <- c(pred,tmp[,2])
  true <- c(true,tmp[,3])}
print(cor(pred,true))
plot(pred,true)



#AL ;17
#SCxTX; 21
#FLxTX; 13,38,43
#GAxTX; 35
#VA; 1
#GA; 23,33,34,36,38,48
#FL; 14,16,32,39,44,45,46,47,49,50,51,52,53,54,55
#NC; 2,3,4,5,7,15,20
#SC; 6,10,11,12,18,22,24,25,26,2728,29,30,31,37,40,41,42,56
#SCxGA; 9
#SCxSC; 19
#FLxFL; 38
#SCxTX; 21
#VA; 1

AL <- 17
SC.TX <- 21
FL.TX <- c(13,38,43)
VA <- 1
GA <-  c(23,33,34,36,38,48)
FL<- c( 14,16,32,39,44,45,46,47,49,50,51,52,53,54,55)
NC <- c(2,3,4,5,7,15,20)
SC <- c(6,10,11,12,18,22,24,25,26,2728,29,30,31,37,40,41,42,56)
SC.GA <- 9
SC.SC <-  19
FL.FL <-  38
SC.TX <-  21

#AL vs. Diff Regions####

#vs. SC.TX #.139
all.fam[[AL]][[SC.TX]] #.139
ALa <- all.fam[[17]]
mean(ALa[c(13,37,42)]) #.18 vs. FL.TX
all.fam[[AL]][[1]]  # .0945 vs. VA

mean(ALa[c(23,33,34,36,38,48)]) #. 151 vs. GA
mean(ALa[c( 14,16,31,38,43,44,45,46,48,49,50,51,52,53,54)]) #. 151 vs. GA
mean(ALa[c(23,33,34,36,38,48)]) #. 151 vs. GA
mean(ALa[c(23,33,34,36,38,48)]) #. 151 vs. GA


mean(ALa[c(23,33,34,36,38,48)]) #. 151 vs. GA
mean(ALa[c(23,33,34,36,38,48)]) #. 151 vs. GA
mean(ALa[c(23,33,34,36,38,48)]) #. 151 vs. GA
mean(ALa[c(23,33,34,36,38,48)]) #. 151 vs. GAall.fam[[AL]][[SC.TX]] 
all.fam[[AL]][[SC.TX]] 
all.fam[[AL]][[SC.TX]] 
all.fam[[AL]][[SC.TX]] 
all.fam[[AL]][[SC.TX]] 
all.fam[[AL]][[SC.TX]] 
all.fam[[AL]][[SC.TX]] 
all.fam[[AL]][[SC.TX]] 

#est_pita_67707934
#AT3G62550.1  | Symbols:  | Adenine nucleotide alpha hydr...    82   2e-15  1 
#est_Pita_61779901
#AT4G11360.1  | Symbols: RHA1B | RING-H2 finger A1B | chr...    34   0.26   1
#est_pita_51498052
#AT5G14780.1  | Symbols: FDH | formate dehydrogenase | ch...   248   1e-65  2 
#isotig09939
#AT5G39040.1  | Symbols: ALS1, ATTAP2, TAP2 | transporter...   264   2e-70  3 
#est_Pita_66977002
#AT1G44740.1  | Symbols:  | unknown protein; FUNCTIONS IN...    34   0.57   1
#est_pita_66981047
#AT5G07990.1  | Symbols: TT7, CYP75B1, D501 | Cytochrome ...   150   2e-49  3 
#est_pita_67551644
#AT4G22020.1  | Symbols:  | pseudogene, hypothetical prot...    30   0.064  2 


##TRY PCA OR SVA####
dim(fam.logcpm2)
fam.logcpm2.7.10 <- fam.logcpm2[,which(colnames(fam.logcpm2) %in% unique(c(corrected.full.comparisons.fdr.05[[7]][[3]],corrected.full.comparisons.fdr.05[[32]][[20]])))]
fam.logcpm2.7.10 <- fam.logcpm2[,which(colnames(fam.logcpm2) %in% c(corrected.full.comparisons.fdr.05[[32]][[20]]))]

p.7.10 <- prcomp((fam.logcpm2))
print(head(p.7.10))
summary(p.7.10)
plot(p.7.10,type="l")
biplot(p.7.10)
plot(p.7.10$x[-10,1],p.7.10$x[-10,2])
s.7.10 <- svd(x = fam.logcpm2.7.10)
plot(s.7.10$v[,1],s.7.10$v[,2])
plot(s.7.10$u[,1])
plot(s.7.10$d^2/sum(s.7.10$d^2), xlim = c(0, 56), type = "b", pch = 16, xlab = "principal components", 
     ylab = "variance explained")
plot(s.7.10$v[,1],s.7.10$v[,2])
library(randomForest)
randomForest()

out <- apply(fam.logcpm2,MARGIN = 1,FUN = function(x) sum(p.7.10$rotation[,1]*x)) 

pc.fam.cpm <- prcomp(fam.logcpm2.7.10, retx = TRUE, center = TRUE, scale. = TRUE)
uniq.fam <- as.character(unique(phenos$group))
fam.phenos <- data.frame(fam=uniq.fam,ht=phenos$ht[match(uniq.fam,phenos$group)],
                          vol=phenos$vol[match(uniq.fam,phenos$group)])
summary(fam.phenos$ht)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 96.8   105.2   108.6   108.0   111.2   119.2 
summary(fam.phenos$vol)
#    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#    91.5   109.0   119.6   118.8   128.2   145.7 
vol.grps <- list(q1=which(fam.phenos$vol<109.0),
                 q2=intersect(which(fam.phenos$vol>=109.0),which(fam.phenos$vol<119.6)),
                 q3=intersect(which(fam.phenos$vol>=119.6),which(fam.phenos$vol<128.2)),
                 q4=which(fam.phenos$vol>=128.2))
max(pc.fam.cpm$x[,1]) # 33
min(pc.fam.cpm$x[,1]) # -20
max(pc.fam.cpm$x[,2]) # 25
min(pc.fam.cpm$x[,2]) # -14
plot(pc.fam.cpm$x[vol.grps$q1,1],pc.fam.cpm$x[vol.grps$q1,2],col="red",xlim=c(-21,35),
     ylim=c(-15,26))
points(pc.fam.cpm$x[vol.grps$q2,1],pc.fam.cpm$x[vol.grps$q2,2],col="orange")     
points(pc.fam.cpm$x[vol.grps$q3,1],pc.fam.cpm$x[vol.grps$q3,2],col="purple")     
points(pc.fam.cpm$x[vol.grps$q4,1],pc.fam.cpm$x[vol.grps$q4,2],col="blue")     
# No visible pattern of differences. Try calculating median value of PCs for each group,
# plot those to look for differences.  
med.q1 <- apply(pc.fam.cpm$x[vol.grps$q1,],2,median)
med.q2 <- apply(pc.fam.cpm$x[vol.grps$q2,],2,median)
med.q3 <- apply(pc.fam.cpm$x[vol.grps$q3,],2,median)
med.q4 <- apply(pc.fam.cpm$x[vol.grps$q4,],2,median)
max(c(med.q1,med.q2,med.q3,med.q4)) # 30
min(c(med.q1,med.q2,med.q3,med.q4)) # -50
plot(1:56,med.q1,col="red",pch=19,ylim=c(-3,3))
points(1:56,med.q2,col="orange",pch=20)
points(1:56,med.q3,col="purple",pch=16)
points(1:56,med.q4,col="blue",pch=20)
# PC2 and PC4 show scatter of median values - plot those against each other
plot(pc.fam.cpm$x[vol.grps$q1,3],pc.fam.cpm$x[vol.grps$q1,5],col="red",xlim=c(-10,10),
     ylim=c(-10,10),pch=19,xlab="PC2",ylab="PC4")
points(pc.fam.cpm$x[vol.grps$q2,3],pc.fam.cpm$x[vol.grps$q2,5],col="magenta",pch=20)     
points(pc.fam.cpm$x[vol.grps$q3,3],pc.fam.cpm$x[vol.grps$q3,5],col="purple",pch=16)     
points(pc.fam.cpm$x[vol.grps$q4,3],pc.fam.cpm$x[vol.grps$q4,5],col="blue",pch=20)  
# Some separation of groups is visible that accounts for differences in medians,
# but still lots of overlap. Patterns (if any) seem to be too complex to see in
# pairwise analysis of PCs.

