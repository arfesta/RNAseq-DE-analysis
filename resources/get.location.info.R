install.packages("ggmap")
library(ggmap)


fam.origins
get.address <- fam.origins[which(fam.origins$state %in% c("NC","SC","GA","FL","AL","VA")),]
addresses <- unique(paste0(get.address$county," County, ",get.address$state))
geo <- vector("list")
for(out in 1:length(addresses)){
geo[[out]] <- geocode(location = addresses[out], output="latlon", source="google")
print(out)
}
geo
out=15

save.image(file = "/mnt/get.location.info.RData",compress=T)


geo.dat <- matrix(,nrow=16,ncol=4)
geo.dat[,1] <- unlist(lapply(geo,function(each.loc) each.loc[[1]]))
geo.dat[,2] <- unlist(lapply(geo,function(each.loc) each.loc[[2]]))
loc.info <- unique(paste0(get.address$county," ",get.address$state))
loc.info <- strsplit(x = loc.info,split = " ")
geo.dat[,3] <- unlist(lapply(loc.info,function(each.loc) each.loc[[1]]))
geo.dat[,4] <- unlist(lapply(loc.info,function(each.loc) each.loc[[2]]))

colnames(geo.dat) <- c("longitude","latitude","county","state")


install.packages("geosphere")
library(geosphere)
geo.dat <- data.frame(geo.dat,stringsAsFactors = F)
geo.dat$latitude <- as.numeric(geo.dat$latitude)
geo.dat$longitude <- as.numeric(geo.dat$longitude)
distGeo(p1 = geo.dat[3,1:2],p2 = geo.dat[2,1:2])

keep.dat <- all.de.matrix[which(fam.origins$state %in% c("NC","SC","GA","FL","AL","VA")),which(fam.origins$state %in% c("NC","SC","GA","FL","AL","VA"))]

distance.matrix <- matrix("NA",nrow=16,ncol =16 )
for(each.county in 1:16){
  these.counties <- c(1:16)[-c(1:each.county)]
  
  for(each.other.county in 1:length(these.counties)){
   distance.matrix[each.county,each.other.county + each.county] <- distGeo(p1 = geo.dat[each.county,1:2],p2 = geo.dat[each.other.county+each.county,1:2])
  }
  
}
colnames(distance.matrix) <- paste0(geo.dat$county,".",geo.dat$state)
rownames(distance.matrix) <- paste0(geo.dat$county,".",geo.dat$state)


diag(distance.matrix) <- 0

test.mat <- (distance.matrix)
lowerTriangle(test.mat) <- upperTriangle(test.mat,byrow = T)
test.mat[lower.tri(test.mat)] <- (test.mat[upper.tri(test.mat)])

these.matches <- paste0(get.address$county,".",get.address$state)
test.mat2 <- test.mat[match(these.matches,colnames(test.mat)),match(these.matches,colnames(test.mat))]

final.distance.matrix <- test.mat2

rm(distance.matrix,geo,geo.dat,get.address)

dim(final.distance.matrix)
dim(keep.dat)
lowerTriangle(keep.dat) <- upperTriangle(keep.dat,byrow = T)
diag(keep.dat) <- 0


distance.rows.1 <- as.numeric(unlist(final.distance.matrix[upper.tri(final.distance.matrix)]))
de.num.rows.1 <- as.numeric(unlist(keep.dat[upper.tri(keep.dat)]))
plot(scale(distance.rows.1[-which(distance.rows.1 == 0)]),de.num.rows.1[-which(distance.rows.1 == 0)])
hist(de.num.rows.1)

plot(scale(as.numeric(final.distance.matrix[1,-1])+1),scale(as.numeric(keep.dat[1,-1])))
