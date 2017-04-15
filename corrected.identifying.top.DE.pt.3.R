

#sample 250 transcripts at random 1540 times and look at outcome####
set.seed(134242)
p <- proc.time()
random.cors.2501 <- vector()
random.txpts.2501 <- vector("list")
for(each in 1:154000){
test <- sample(rownames(all.txpts),250)
random.txpts.2501[[each]] <- test
log.DE <- scaled.fam.logcpm2[,match(test,colnames(scaled.fam.logcpm2))]
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
random.cors.2501 <- c(random.cors.2501,cor(pred,true))
print(each)}
proc.time() - p
summary(random.cors.250)
#Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
#-0.39180  0.06856  0.15780  0.15750  0.24740  0.70050 
 summary(random.cors.2501)
#Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
#-0.39830  0.06881  0.15830  0.15780  0.24720  0.67650 

t.test(random.cors.250,random.cors.2501) # p=.56
hist(random.cors.250,breaks = seq(-.4,.8,.2))
hist(random.cors.2501, breaks=seq(-.4,.8,.2))
hist(fam.cor.05)

#Compare transcript crossover between top .01 fdr pair and top from random
which.max(random.cors.250)
#88574
r.txt.max <- random.txpts.250[[88574]]
txt.max.01 <- comparisons.fdr.01[[10]][[24]]
length(which(r.txt.max %in% txt.max.01))
#2 <- only 3 out of 243 are in the 250

#see how many random transcripts are similar in everything above .5
which(random.cors.250 > .6)
#163  302  371 1218
length(unique(c(random.txpts.250[[163]],random.txpts.250[[302]],random.txpts.250[[371]],random.txpts.250[[1218]])))
#968/1000 are unique. roughly 97%

##WGCNA####
library(WGCNA)
#test <- good.v.bad.txpts[!(good.v.bad.txpts %in% unique(c(good.transcripts,bad.transcripts)))]
#test <- comparisons.fdr.05[[8]][[9]]
test <- unique(c(comparisons.fdr.05[[1]][[47]],comparisons.fdr.05[[7]][[3]],comparisons.fdr.05[[8]][[9]],
                 comparisons.fdr.05[[8]][[22]],comparisons.fdr.05[[8]][[29]],comparisons.fdr.05[[10]][[24]],
                 comparisons.fdr.05[[17]][[34]],comparisons.fdr.05[[18]][[37]],comparisons.fdr.05[[32]][[21]],
                 comparisons.fdr.05[[37]][[2]]))
datExpr0 <- fam.logcpm2[,match(test,colnames(scaled.fam.logcpm2))]
gsg = goodSamplesGenes(datExpr0, verbose = 3)
sampleTree = hclust(dist(datExpr0), method = "average");
# Plot the sample tree: Open a graphic output window of size 12 by 9 inches
# The user should change the dimensions if the window is too large or too small.
sizeGrWindow(12,9)
#pdf(file = "Plots/sampleClustering.pdf", width = 12, height = 9);
par(cex = 0.6);
par(mar = c(0,4,2,0))
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5,
     cex.axis = 1.5, cex.main = 2)

###If we want to remove samples
# Plot a line to show the cut
#abline(h = 20, col = "red");
# Determine cluster under the line
#clust = cutreeStatic(sampleTree, cutHeight =20, minSize = 10)
#table(clust)
# clust 1 contains the samples we want to keep.
#keepSamples = (clust==1)
datExpr = datExpr0
nGenes = ncol(datExpr)
nSamples = nrow(datExpr)

# Re-cluster samples
sampleTree2 = hclust(dist(datExpr), method = "average")
# Convert traits to a color representation: white means low, red means high, grey means missing entry
datTraits <- phenos.fam$x[]
traitColors = numbers2colors(datTraits, signed = FALSE);
# Plot the sample dendrogram and the colors underneath.
plotDendroAndColors(sampleTree2, traitColors,
                    groupLabels = names(datTraits),
                    main = "Sample dendrogram and trait heatmap")



# Choose a set of soft-thresholding powers
powers = c(c(1:10), seq(from = 12, to=20, by=2))
# Call the network topology analysis function
sft = pickSoftThreshold(datExpr, powerVector = powers, verbose = 5)
# Plot the results:
sizeGrWindow(9, 5)
par(mfrow = c(1,2));
cex1 = 0.9;
# Scale-free topology fit index as a function of the soft-thresholding power
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red");
# this line corresponds to using an R^2 cut-off of h
abline(h=0.90,col="red")
# Mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")

net = blockwiseModules(datExpr, power = 4,
                       TOMType = "unsigned", minModuleSize = 30,
                       reassignThreshold = 0, mergeCutHeight = 0.25,
                       numericLabels = TRUE, pamRespectsDendro = FALSE,
                       saveTOMs = F,
                       saveTOMFileBase = "femaleMouseTOM",
                       verbose = 3)

table(net$colors)

# open a graphics window
sizeGrWindow(12, 9)
# Convert labels to colors for plotting
mergedColors = labels2colors(net$colors)
# Plot the dendrogram and the module colors underneath
plotDendroAndColors(net$dendrograms[[1]], mergedColors[net$blockGenes[[1]]],
                    "Module colors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)

moduleLabels = net$colors
moduleColors = labels2colors(net$colors)
MEs = net$MEs;
geneTree = net$dendrograms[[1]];

# Define numbers of genes and samples
nGenes = ncol(datExpr);
nSamples = nrow(datExpr);
# Recalculate MEs with color labels
MEs0 = moduleEigengenes(datExpr, moduleColors)$eigengenes
MEs = orderMEs(MEs0)
moduleTraitCor = cor(MEs, datTraits, use = "p");
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples);


sizeGrWindow(10,6)
# Will display correlations and their p-values
textMatrix = paste(signif(moduleTraitCor, 2), "\n(",
                   signif(moduleTraitPvalue, 1), ")", sep = "");
dim(textMatrix) = dim(moduleTraitCor)
par(mar = c(6, 8.5, 3, 3));
# Display the correlation values within a heatmap plot
labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = names(datTraits),
               yLabels = names(MEs),
               ySymbols = names(MEs),
               colorLabels = FALSE,
               colors = greenWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.5,
               zlim = c(-1,1),
               main = paste("Module-trait relationships"))



# Calculate topological overlap anew: this could be done more efficiently by saving the TOM
# calculated during module detection, but let us do it again here.
dissTOM = 1-TOMsimilarityFromExpr(datExpr, power = 4);
# Transform dissTOM with a power to make moderately strong connections more visible in the heatmap
plotTOM = dissTOM^7;
# Set diagonal to NA for a nicer plot
diag(plotTOM) = NA;
# Call the plot function
sizeGrWindow(9,9)
TOMplot(plotTOM, geneTree, moduleColors, main = "Network heatmap plot, all genes")


# Recalculate module eigengenes
MEs = moduleEigengenes(datExpr, moduleColors)$eigengenes
# Isolate weight from the clinical traits
weight = as.data.frame(phenos.fam$x);
names(weight) = "weight"
# Add the weight to existing module eigengenes
MET = orderMEs(cbind(MEs, weight))
# Plot the relationships among the eigengenes and the trait
sizeGrWindow(5,7.5);
par(cex = 0.9)
plotEigengeneNetworks(MET, "", marDendro = c(0,4,1,2), marHeatmap = c(3,4,1,2), cex.lab = 0.8, xLabelsAngle
                      = 90)
# Plot the dendrogram
sizeGrWindow(6,6);
par(cex = 1.0)
plotEigengeneNetworks(MET, "Eigengene dendrogram", marDendro = c(0,4,2,0),
                      plotHeatmaps = FALSE)
# Plot the heatmap matrix (note: this plot will overwrite the dendrogram plot)
par(cex = 1.0)
plotEigengeneNetworks(MET, "Eigengene adjacency heatmap", marHeatmap = c(3,4,2,2),
                      plotDendrograms = FALSE, xLabelsAngle = 90)


###Step-by-step####
# Choose a set of soft-thresholding powers
powers = c(c(1:10), seq(from = 12, to=20, by=2))
# Call the network topology analysis function
sft = pickSoftThreshold(datExpr, powerVector = powers, verbose = 5)
# Plot the results:
sizeGrWindow(9, 5)
par(mfrow = c(1,2));
cex1 = 0.9;
# Scale-free topology fit index as a function of the soft-thresholding power
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red");
# this line corresponds to using an R^2 cut-off of h
abline(h=0.90,col="red")
# Mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")

softPower = 7;
adjacency = adjacency(datExpr, power = softPower);
# Turn adjacency into topological overlap
TOM = TOMsimilarity(adjacency);
dissTOM = 1-TOM

# Call the hierarchical clustering function
geneTree = hclust(as.dist(dissTOM), method = "average");
# Plot the resulting clustering tree (dendrogram)
sizeGrWindow(12,9)
plot(geneTree, xlab="", sub="", main = "Gene clustering on TOM-based dissimilarity",
     labels = FALSE, hang = 0.04);

# We like large modules, so we set the minimum module size relatively high:
minModuleSize = 30;
# Module identification using dynamic tree cut:
dynamicMods = cutreeDynamic(dendro = geneTree, distM = dissTOM,
                            deepSplit = 2, pamRespectsDendro = FALSE,
                            minClusterSize = minModuleSize);
table(dynamicMods)

# Convert numeric lables into colors
dynamicColors = labels2colors(dynamicMods)
table(dynamicColors)
# Plot the dendrogram and colors underneath
sizeGrWindow(8,6)
plotDendroAndColors(geneTree, dynamicColors, "Dynamic Tree Cut",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05,
                    main = "Gene dendrogram and module colors")

# Calculate eigengenes
MEList = moduleEigengenes(datExpr, colors = dynamicColors)
MEs = MEList$eigengenes
# Calculate dissimilarity of module eigengenes
MEDiss = 1-cor(MEs);
# Cluster module eigengenes
METree = hclust(as.dist(MEDiss), method = "average");
# Plot the result
sizeGrWindow(7, 6)
plot(METree, main = "Clustering of module eigengenes",
     xlab = "", sub = "")

MEDissThres = 0.25
# Plot the cut line into the dendrogram
abline(h=MEDissThres, col = "red")
# Call an automatic merging function
merge = mergeCloseModules(datExpr, dynamicColors, cutHeight = MEDissThres, verbose = 3)
# The merged module colors
mergedColors = merge$colors;
# Eigengenes of the new merged modules:
mergedMEs = merge$newMEs;

sizeGrWindow(12, 9)
#pdf(file = "Plots/geneDendro-3.pdf", wi = 9, he = 6)
plotDendroAndColors(geneTree, cbind(dynamicColors, mergedColors),
                    c("Dynamic Tree Cut", "Merged dynamic"),
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)
#dev.off()

# Define numbers of genes and samples
nGenes = ncol(datExpr);
nSamples = nrow(datExpr);
# Recalculate MEs with color labels
MEs0 = moduleEigengenes(datExpr, moduleColors)$eigengenes
MEs = orderMEs(MEs0)
moduleTraitCor = cor(MEs, datTraits, use = "p");
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples);

sizeGrWindow(10,6)
# Will display correlations and their p-values
textMatrix = paste(signif(moduleTraitCor, 2), "\n(",
                   signif(moduleTraitPvalue, 1), ")", sep = "");
dim(textMatrix) = dim(moduleTraitCor)
par(mar = c(6, 8.5, 3, 3));
# Display the correlation values within a heatmap plot
labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = names(datTraits),
               yLabels = names(MEs),
               ySymbols = names(MEs),
               colorLabels = FALSE,
               colors = greenWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.5,
               zlim = c(-1,1),
               main = paste("Module-trait relationships"))


########F######
for(i in 0:8){
test <-  colnames(datExpr)[which(net$colors %in% c(0,4))]
log.DE <- scaled.fam.logcpm2[,match(test,colnames(scaled.fam.logcpm2))]
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
plot(pred,true) 
print(cor(pred,true))}



test <- t(p.7.10$rotation)
rf <- randomForest(x=t(cor(t(scaled.fam.logcpm2[,test]))),y=phenos.fam$x,importance = T,keep.forest = T,oob.prox = T,do.trace = T)
g <- cbind(p.7.10$x,phenos.fam$x)
rf <- randomForest(V57 ~ .,g,mtry=5)
plot(rf$predicted,phenos.fam$x)
cor(rf$predicted,phenos.fam$x)
