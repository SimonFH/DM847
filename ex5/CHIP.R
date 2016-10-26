args <- commandArgs(trailingOnly = TRUE)
if (length(args)!=1){stop("Script takes arguments: (inputfile path)")}

#rmap <- read.csv("./ex5_readmapping.txt", sep=" ")
rmap <- read.csv(args[1], sep=" ")
colnames(rmap) <- c("start", "stop", "direction")
size <- 1000

#hst <- by(rmap$start, rmap$direction, function(x){table(factor(x, levels=1:size))})
#hst <- by(rmap, rmap$direction, FUN=function(x){apply(x, 1, function(y){ table(factor(seq(x$start,x$stop,1), levels=1:size))})
# use "by" to split by direction, then apply, per row, a table of the sequence on the levels 1-1000. Finally use rowSums on these. 
hst <- by(rmap, rmap$direction, FUN=function(x){ rowSums(apply(x, 1, function(y){table(factor(seq(as.integer(y["start"]),as.integer(y["stop"]),1), levels=1:size))})) })
names(hst) <- c("fwd", "rev")

## should it be reversed????
#hst$rev <- rev(hst$rev)

dns <- data.frame() 
k <- 5
for(i in 1:size){
  left <- max(1, i-k)
  right <- min(size, i+k)
  nvals <- right-left+1
  dns <- rbind(dns,rbind(Map(hst, f=function(x){sum(x[left:right])/nvals})))
  #dns$fwd[i] <- sum(hst$fwd[left:right])/nvals
}


#c
locMax <- function(x,i){
  left <- max(1, i-1)
  right <- min(size, i+1)
  idx <- unlist(x[(left:right)])
  if(length(idx)==3){
    return(idx[1]<=idx[2]&&idx[2]>idx[3])
  }else{
    if(left==i)
      return(idx[1]>idx[2])
    if(right==i)
      return(idx[1]<=idx[2])
  }  
}

# find local max locations
lmidx <- as.data.frame(apply(dns, 2, FUN=function(s){ unlist(Map(f=function(i){ locMax(s,i) }, 1:size)) }) )
rownames(lmidx) <- 1:size

# d - filter all peaks with < 100 reads
lmidx2 <- apply(lmidx, 2, which)


# indices of local maxima smaller than 100
smidx <- Map(names(dns), f=function(x){y <- c(lmidx2[x][[1]]); z <- dns[,x][y]<100;y[z]})

# indices of local maxima larger than or equal to 100
lmidx3 <- Map(names(dns), f=function(x){y <- c(lmidx2[x][[1]]); z <- dns[,x][y]>=100;y[z]})
#names(lmidx3) <- NULL
names(lmidx3[[1]]) <- NULL
names(lmidx3[[2]]) <- NULL

#dns2 <- dns
#for(n in names(smidx)){
#  for(i in smidx[n][[1]]){
#    dns2[i,n] <- 0
#  }
#}


# e)

# Assuming "+" is the coding strand, we look between 120 and 200 bp upstream for a peak on the "-" strand.
res <- NULL
for (x in lmidx3$fwd){
  q <- ((x+100):(x+250))
  z <- q%in%lmidx3$rev
  if(sum(z)>0){
    res <- c(res,q[z][which.max(unlist(dns$rev[q[z]]))])
  }else{
    res <- c(res, 0)
  }
}
result <- data.frame(fwd=lmidx3$fwd, rev=res)
print(result)
write.csv(result, "./bs_pos.txt")



# Plotting 
png(filename="./out.png", width=16, height=10, res=120, units="in")
# draw histogram(ish) and density
par(cex.axis=1)
plot(1:1000, hst$fwd,col="orange", type="l", ylim=c(-(max(hst$fwd)+50), max(hst$fwd)+50), xlab="index", ylab="value")
lines(1:1000, dns$fwd,col="blue", type="l")
lines(1:1000, -hst$rev,col="orange", type="l")
lines(1:1000, Map('-', dns$rev) ,col="blue", type="l")

# draw limits
lines(1:1000, rep(100, 1000), type="l", lty="dashed", col="red")
lines(1:1000, rep(-100, 1000), type="l", lty="dashed", col="red")

#draw lines connecting peaks on opposing strands
apply(result, 1, FUN=function(x){
        lines(x=c(x["fwd"],x["rev"]), y=c(dns$fwd[x["fwd"]], -unlist(dns$rev[x["rev"]])), col="black", lty="dotted", cex=2)
  }
)
#draw peaks
points(lmidx3$fwd, dns$fwd[lmidx3$fwd], pch=20, cex=2, col="black")
points(lmidx3$fwd, dns$fwd[lmidx3$fwd], pch=20, cex=1, col="cyan")
points(lmidx3$rev, -unlist(dns$rev[lmidx3$rev]), pch=20, cex=2, col="black")
points(lmidx3$rev, -unlist(dns$rev[lmidx3$rev]), pch=20, cex=1, col="magenta")

axis(side=1,at=seq(50,1000,50))
for(x in lmidx3$fwd){
  abline(v=x, lty="dotted", col="cyan",lwd=.5)
}
for(x in lmidx3$rev){
  abline(v=x, lty="dotted",  col="magenta", lwd=.5)
}
abline(h=0, lwd=.75)
dev.off()

