# FUNCTION DECLARATIONS

stol <- function(s){unlist(strsplit(s,''))}
# order such that EOF character has lowest value
cidx <- function(c){match(toupper(c),c(eofc,LETTERS))}
cmp1 <- function(s1, s2){
  sl1 <- stol(s1)
  sl2 <- stol(s2)
  for(i in (1:min(length(sl1),length(sl2)))){
    idx1 <- cidx(sl1[i])
    idx2 <- cidx(sl2[i])
    if(idx1>idx2){
      return(1)
    }
    if(idx1<idx2){
      return(-1)
    }
  }
  return(0)
}


swap <- function(ord, i,j){
  a <- ord[i]
  ord[i] <- ord[j]
  ord[j] <- a
  return(ord)
}
partition <- function(A, lo, hi, ord){
    pivot <- A[ord[hi]]
    i <- lo
    for (j in (lo:(hi - 1))){
        if (cmp1(A[ord[j]], pivot)<=0){
            ord <- swap(ord, i,j)
            i <- i + 1
        }
    }
    ord <- swap(ord,i,hi)
    return(list(ord,i))
}
qs_lex_order <- function(A,lo, hi, ord){
    if (lo < hi){
        val <- partition(A, lo, hi, ord)
        ord <- val[[1]]
        p <- val[[2]]
        ord <- qs_lex_order(A, lo, p - 1, ord)
        ord <- qs_lex_order(A, p + 1, hi, ord)
    }
    return(ord)
}

# PROGRAM START

eofc <- '$'
args <- commandArgs(trailingOnly = TRUE)
if (length(args)<1){stop("Input format: {(input string)|-d (input file path)} [-c (EOF char.)] ")}
c <- args[((1:length(args))[args=="-c"]+1)]
if(length(c)==1)
  eofc <- as.character(c)

# encode
if(length(args[((1:length(args))[args=="-d"])])==0){
  s <- paste(args[1],eofc, sep="")
  ss <- unlist(strsplit(s, ''))
  l <- length(ss)
  m <- Map((0:(l-1)), f=function(i){c(ss[(i+1):l],ss[(0:i)])})
  m <- sapply(m, FUN=paste, collapse='')
  ord <- qs_lex_order(m, 1, length(m), ord=(1:length(m)))
  lastcol <- t(sapply(strsplit(m[ord],split=''), unlist))[,length(ss)]
  str <- paste(lastcol, collapse='')
  print(str)
  writeChar(str, "./bwt.txt")
} else {
  fp <- args[((1:length(args))[args=="-d"]+1)]
  lc <- readLines(fp,warn=F)
  lcl <- unlist(strsplit(lc, ''))
  l <- length(lcl)
  m <- Map((0:(l-1)), f=function(i){c(lcl[(i+1):l],lcl[(0:i)])})
  ord <- qs_lex_order(lcl, 1, l, ord=(1:l))
  sorted <- Map((0:(l-1)), f=function(i){c(lcl[ord][(i+1):l],lcl[ord][(0:i)])})[[1]]
  counts <- table(lcl)
  i <- 1
  s <- lcl[1]
  while(T){
    ch <- lcl[i]
    cts <- 0;for(k in 1:i){if(lcl[k]==ch)cts <- cts+1}
    #j <- which(sorted==lcl[i])[[ counts[lcl[i]] ]]
    j <- which(sorted==lcl[i])[[ cts ]]
    counts[lcl[i]] <- counts[lcl[i]]-1
    s <- c(lcl[j],s)
    #sorted[j] <- NA
    i <- j
    if(lcl[i]==eofc)break
  }
  s <- s[s!=eofc]
  print("Original string:")
  print(paste(s,collapse=""))
}




