args <- commandArgs(trailingOnly = TRUE)
if (length(args)!=1){stop("Script takes one string as argument.")}
s <- paste(args[1],"$", sep="")
ss <- unlist(strsplit(s, ''))
su <- sapply(1:length(ss), function(x){ss[x:length(ss)]})
qwe <- sapply(X=su, FUN=paste, collapse='')

#Sys.setlocale("LC_COLLATE", "C")
suf <- data.frame(suf=(0:(length(qwe)-1)),S_suf=qwe) 




#fla <- list()
#for (i in (1:length(su))){
#  fla[[i]] <- sapply(su[[i]], function(x){match(toupper(x), c(LETTERS,"$"))})
#}
#
## now sort using custom comparator
#fli <- list();
#ord <- (1:length(fla))
#j <- 1
#for (i in (1:length(fla))){
#  a <- fla[[ord[i]]][j]
#  
#};

stol <- function(s){unlist(strsplit(s,''))}
cidx <- function(c){match(toupper(c),c(LETTERS,"$"))}
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

ord <- qs_lex_order(qwe,1,length(qwe), (1:length(qwe)))

sufs <- data.frame(i=(0:(nrow(suf)-1)), suf[ord,])
#rownames(sufs) <- 1:(length(sufs))
#sufs[order(sufs$suf),]



lcp_cmp <- function(s1, s2, i){all.equal(strsplit(as.character(s1),'')[[1]][(1:i)], strsplit(as.character(s2),'')[[1]][(1:i)])}
lcp_cmp2 <- function(s1,s2, i){all.equal(s1[[1]][(1:i)], s2[[1]][(1:i)])}
lcp_fun <- function(s1,s2){
  val <- TRUE
  i <- 0
  while(val == TRUE){
    i <- i+1
    val <- lcp_cmp2(s1,s2,i)
  }
  i-1
}

sus <- su[ord]
lcp <- list(0)
for (i in 2:nrow(sufs)){ lcp[i] <- lcp_fun(sus[i-1], sus[i])}

#suf2 <- data.frame(sufs, lcp=cbind(lcp))

lcp <- unlist(lcp)
l <- length(lcp)
skp <- rep(NA, length(lcp))
for( i in (1:l)){
  for (j in (i:l)){
    if(lcp[j]<lcp[i]){
      skp[i] <- j-1
      break
    }
  }
  if(is.na(skp[i])){
   skp[i] <- l 
  }
}


res <- data.frame(i=sufs$i, lcp=lcp, skp=skp, s_suf=sufs$S_suf)
print("Enhanced suffix array:")
print(res)
write.csv(res, "./ESA.csv")
