
args <- commandArgs(trailingOnly = TRUE)
if (length(args)<2){stop("Script takes arguments: (inputfile path) (motif length) [-f CSV|FAS]")}

format <- args[((1:length(args))[args=="-f"]+1)]
format <- toupper(format)
if(length(format)>0 && format == "CSV"){
  foo <- read.csv(args[1], row.names=1)
} else {
  baz <- read.table(args[1], comment.char='>')
  #data frame holding individual characters
  foo <- data.frame(do.call('rbind', strsplit(as.character(baz[,1]),'',fixed=T)))
}


#length of motif
W <- as.integer(args[2]) 
#number of sequences
nseq <- nrow(foo)


#types, typically acgt
types <- levels(unlist(foo))


# length of sequences
L <- ncol(foo)

#character count
cc <- table(unlist(foo))
#background distribution
bgdist <- cc/sum(cc)

# number of possible positions for the motif
K <- L-W+1

# guess a start position
motif_init <- function(){
  start = runif(1,min=1, max=K)
  print(c("Guessing start position: ", as.integer(start)))
  mCM <- sapply(foo[,(start:(start+W-1))], function(x){table(factor(x, types))})
  #motif probability matrix
  mPM <- mCM/nrow(foo)
  mPM
}
mPM <- motif_init()


# Primary function that calculates probability matrix for the motif and the Z matrix 
updatemodel <- function(dat, W,K,L, bgdist, mPM){
  Pr <- matrix(0,nrow=nseq, ncol=0)
  for (j in (1:K)){
    a <- matrix()
    b <- matrix()
    c <- matrix()
    if(j>1){
      a <- as.matrix(sapply(FUN=cbind, X=Map(dat[,(1:(j-1))], f=function(x){bgdist[as.character(x)]})))
    }
    b <- sapply(FUN=cbind, X=Map((j:(j+W-1)), f=function(k){cbind(mPM[as.character(dat[,k]), k-j+1])}))

    if((j+W)<L){
      c <- as.matrix(sapply(FUN=cbind, X=Map(dat[,((j+W):L)], f=function(x){bgdist[as.character(x)]})))
    }
    if(!is.na(a)[1,1]){
      b <- cbind(a,b)
    }
    if(!is.na(c)[1,1]){
      b <- cbind(b,c)
    }

    if(j==51)write.csv(b, "/tmp/51s.csv")

    Pr <- cbind(Pr, apply(b,1,FUN=prod))
    #Pr <- cbind(Pr,unlist(Map((1:nrow(dat)),f=function(i){Reduce(f='*', x=t(t(b))[i,], init=1)})))
  }

  Z <- Pr/(rowSums(Pr))

  foo2 <- as.data.frame(matrix(NA, 0, W+1))
  for (i in (1:nrow(Z))){
    for (j in (1:ncol(Z))){
      if (Z[i,j]>0){
        a <- dat[i,(j:(j+W-1))]
        colnames(a) <- (1:W)
        foo2 <- rbind(foo2, cbind(a, Z[i,j]))
      }
    }
  }
  
  #Motif probability matrix
  mPM <- t(as.data.frame(Map(types, f=function(x){colSums((foo2[,(1:W)]==x)*foo2[,(W+1)])/nseq})))
  print(mPM)
  return(list(dat, W,K,L,bgdist, mPM, Z))
}


# Iterate until convergence
updateloop <- function(dat, W, K, L, bgdist, mPM){
  i = 1
  res1 <- updatemodel(dat, W, K, L, bgdist, mPM)
  repeat{
    print(i)
    res2 <- updatemodel(dat, W, K, L, bgdist, res1[[6]])
    if(all.equal(res1[[6]], res2[[6]])==T){
      print(paste(c("CONVERGED AFTER  ",i," ITERATIONS!"), collapse=''))
      return(res2)
    }
    i <- i+1

    res1 <- updatemodel(dat, W, K, L, bgdist, res2[[6]])
  }
}
finalres <- updateloop(foo, W, K, L, bgdist, mPM)

Z <- finalres[[7]]
startpos <- unlist(Map((1:nrow(Z)), f=function(x){which.max(Z[x,])[1]}))

print("Start positions:")
print(startpos)

candidates <- matrix(NA, nrow=0, ncol=W)
candidates <- Map((1:nseq), f=function(x){foo[x,((startpos[x]):(startpos[x]+W-1))]})
cand <- matrix(nrow=0, ncol=W)
colnames(cand) <- (1:W)
for(i in (1:nseq)){colnames(candidates[[i]]) <- (1:W)}
for(i in (1:nseq)){ cand <- rbind(cand, candidates[[i]])}
#print("Motif candidates: ")
#print(cand)

# Write candidates matrix to file
write.csv(cand, "./candidates.csv")

# Also write the background distribution
write.csv(bgdist, "./bgdist.csv")

