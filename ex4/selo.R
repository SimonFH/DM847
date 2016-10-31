args <- commandArgs(trailingOnly = TRUE)
if (length(args)<1){stop("Script takes arguments: (inputfile path) [-t DNA|LETTERS ] [-s IC|RE] [-f CSV|FAS] [-bg path]")}


format <- args[((1:length(args))[args=="-f"]+1)]
format <- toupper(format)
if(length(format)>0 && format == "CSV"){
  foo <- read.csv(args[1], row.names=1)
} else {
  baz <- read.table(args[1])
  foo <- data.frame(do.call('rbind', strsplit(as.character(baz[,1]),'',fixed=TRUE)))
}


type <- args[((1:length(args))[args=="-t"]+1)]
if (length(type)>0 && toupper(type) == "LETTERS"){
  bgdist <- setNames(rep(1/26, 26), LETTERS)
} else {
  bgdist <- c('A'=0.2, 'T'=0.2, 'C'=0.3, 'G'=0.3)
}

# if background distribution supplied, overwrite the above
bgdpath <- args[((1:length(args))[args=="-bg"]+1)]
if(length(bgdpath)>0){
  bgcsv <- read.csv(bgdpath)
  bgdist <- setNames(bgcsv[,3],toupper(bgcsv[,2]))
}



types <- levels(unlist(foo))
names(foo) <- paste("pos", 1:ncol(foo))

# Background distr

# Absolute PWM (count matrix)
aPWM <- sapply(foo, function(x){table(factor(x, types))})
rownames(aPWM)<-toupper(rownames(aPWM))

# Relative PWM
rPWM <- aPWM/nrow(foo)
rownames(rPWM)<-toupper(rownames(rPWM))

# Relative PWM with pseudo counts
# add pseudo counts to absolute
apPWM <- apply(aPWM, 2, FUN=function(x){x + bgdist[names(x)]*length(types)})

# calculate relative with pseudo counts, aka. "weight matrix"
pPWM <- apPWM/(nrow(foo)+length(types))


types <- toupper(types)

# Entropy
entr <- colSums(t(-sapply(types, function(x){rPWM[x,]*log2(rPWM[x,])})), na.rm=T)

# Relative entropy
relEntr <- colSums(t(sapply(types, function(x){rPWM[x,]*log2(rPWM[x,]/bgdist[x])})), na.rm=T)

# Information content
ic <- log2(length(types))-entr

# Contributions for each position, letters sorted by value
bar <- c();
#for (x in 1:ncol(rPWM)){bar[x] <- list(sort(rPWM[,x],decreasing=T))}
for (x in 1:ncol(rPWM)){bar[x] <- list(sort(rPWM[,x],decreasing=F))}
contr <- Map("*", bar, relEntr)

scale <- args[((1:length(args))[args=="-s"]+1)]
if (length(scale)>0 && toupper(scale) == "IC"){
  scale <- "IC"
  ylabel <- "Information content"
  maxValue <- max(ic)
} else {
  scale <- "RE"
  ylabel <- "Relative entropy"
  maxValue <- max(relEntr)
}

# Printing stuff
print("Position count matrix:")
print(aPWM)
write.csv(aPWM, file="./PCM.csv")
print("Position weight matrix (with pseudo counts):")
print(pPWM)
write.csv(pPWM, file="./PWM.csv")

# Plotting stuff

#open device, ie. the file we want to write to
png(filename="./out.png", width=16, height=10, res=120, units="in")

# Define some pretty colours
#cols = sample(colours(), length(types))
cols <- sample(rainbow(26), length(types))
names(cols) <- types

# Make an empty plot, scaled to fit the data.
plot(x=c(1,(ncol(foo)+1)), y=c(0,maxValue),xlab="sequence", ylab=ylabel, pch='')

# Iterate contributions. One element (x), holds a list of values for all letters.
for (i in 1:length(contr)){
  x <- contr[i]
  xdf <- t(as.data.frame(x, col.names=c("value")))
  if(scale == "IC"){
    scaleFactor <- ic[i]/sum(xdf)
  } else {
    scaleFactor <- 1
  }
  rest <- 0
  # Iterate over letters
  for (y in colnames(xdf)){
    val <- xdf[,y] * scaleFactor
    # Only draw if it appears at this position
    if (val > 0){
      xleft <- i
      ybottom <- rest
      xright <- i+1
      ytop <- rest+val
      rect(xleft=xleft, ybottom=ybottom, xright=xright, ytop=ytop, col=cols[y])
      text(x=(xright-.5),y=ytop-strheight(y), labels=y)
      rest <- rest+val
    }
  }
}
# write to file
dev.off()

cons <- paste(unlist(Map(contr, f=function(x){names(which.max(x))})), collapse='')
print("Consensus sequence:")
print(cons)
writeChar(cons, "./consensus.txt")
