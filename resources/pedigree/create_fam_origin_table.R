library(DESeq2)
load("/mnt/media/disk6/ARF/DiffExpression/ddstxi_bio_0.v2.RData")

library(readxl)
pedigree <-read_xlsx("/mnt/media/disk6/ARF/DiffExpression/RNAseq-DE-analysis/resources/pedigree/pedigree_2018_04_19_20_36_15.xlsx")
#Note *** The pedigree returned from tiproot does not contain any WG families

ddsTxi$fam_id

unique.parents <- strsplit(x = unique(as.character(ddsTxi$fam_id)),split = "x")
head(unique.parents)
#78 Unique Parents..Including "UC" which is a Pollen Mix

# For each of the 70 crosses, loop thru each cross and identify the county and origin of each parent
## If the length is one then the other parent is OP and if it's two then it may be pollen-pool OR FS

County_origin <- lapply(unique.parents,function(each.cross){
  num.fams <- length(each.cross)
  if(num.fams == 1){ 
    pedigree$County[which(as.character(pedigree$Id) %in% each.cross)]
  } else { p1 <- pedigree$County[which(as.character(pedigree$Id) %in% each.cross[1])]
    p2 <- pedigree$County[which(as.character(pedigree$Id) %in% each.cross[2])]
    c(p1,p2)}
})

State_origin <- lapply(unique.parents,function(each.cross){
  num.fams <- length(each.cross)
  if(num.fams == 1){ 
    pedigree$State[which(as.character(pedigree$Id) %in% each.cross)]
  } else { p1 <- pedigree$State[which(as.character(pedigree$Id) %in% each.cross[1])]
  p2 <- pedigree$State[which(as.character(pedigree$Id) %in% each.cross[2])]
  c(p1,p2)}
})


## Create final table

fam.origins <- data.frame("p1"=unlist(lapply(unique.parents,function(each.cross){each.cross[1]})),
           "p2"=unlist(lapply(unique.parents,function(each.cross){each.cross[2]})),
           "p1.county"=unlist(lapply(County_origin,function(each.cross){each.cross[1]})),
           "p2.county"=unlist(lapply(County_origin,function(each.cross){each.cross[2]})),
           "p1.state"=unlist(lapply(State_origin,function(each.cross){each.cross[1]})),
           "p2.state"=unlist(lapply(State_origin,function(each.cross){each.cross[2]})),
           stringsAsFactors = F)
## If county origin for P1 is missing its WG, lets add those
fam.origins[which(is.na(fam.origins$p1.county) == T),c("p1.county","p1.state")] <- "TX"

# # If P2 exists but there is not state information, it's because its WG family.
TX.crosses <- which(is.na(fam.origins$p2) == F & is.na(fam.origins$p2.state) == T)
fam.origins[TX.crosses,c("p2.county","p2.state")] <- "TX"

# Now we can paste all state and county data for crosses
FS.crosses <- which(is.na(fam.origins$p2) == F)
fam.origins[FS.crosses,"p1.county"] <- paste(fam.origins[FS.crosses,"p1.county"],fam.origins[FS.crosses,"p2.county"],sep = "x")
fam.origins[FS.crosses,"p1.state"] <- paste(fam.origins[FS.crosses,"p1.state"],fam.origins[FS.crosses,"p2.state"],sep = "x")
fam.origins[FS.crosses,"p1"] <- paste(fam.origins[FS.crosses,"p1"],fam.origins[FS.crosses,"p2"],sep = "x")

colnames(fam.origins)
fam.origins <- fam.origins[,-c(2,4,6)]
colnames(fam.origins) <- c("family","county","state")

# Save final table
save(x = fam.origins,file = "/mnt/media/disk6/ARF/DiffExpression/RNAseq-DE-analysis/resources/pedigree/fam_origins.Rdata",compress = T)

