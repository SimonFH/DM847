# MOSTLY VARIABLE INITIALIZATION

args <- commandArgs(trailingOnly = TRUE)

foo <- read.table(file=args[1], skip=7, sep=" ", comment.char="E", row.names=1)
names(foo) <- c("x","y")

l <- nrow(foo)
all <- 1:l


#parametric exponents, alpha influences pheromonevalue, beta influences heuristic value
#evap_rate describes the intrinsic degredation rate of the pheromones
#n_ants is the number of virtual ants
#max_its is the maximum number of iterations allowed
n_ants <- args[((1:length(args))[args=="-n"]+1)]
if(length(n_ants)>0){
  #ok?
  n_ants <- as.integer(n_ants)
}else{
  n_ants <- 10
}
evap_rate <- args[((1:length(args))[args=="-r"]+1)]
if(length(evap_rate)>0){
  #ok?
  evap_rate <- as.numeric(evap_rate)
}else{
  evap_rate <- 0.1
}
alpha <- args[((1:length(args))[args=="-a"]+1)]
if(length(alpha)>0){
  #ok?
  alpha <- as.numeric(alpha)
}else{
  alpha <- 1
}
beta <- args[((1:length(args))[args=="-b"]+1)]
if(length(beta)>0){
  #ok?
  beta <- as.numeric(beta)
}else{
  beta <- 1
}
max_its <- args[((1:length(args))[args=="-m"]+1)]
if(length(max_its)>0){
  #ok?
  max_its <- as.integer(max_its)
}else{
  max_its <- 10
}


#eucledean distance function
euc_dist <- function(p,q){
  sqrt((p$x-q$x)**2+(p$y-q$y)**2)
}

#euclidean distances
dist <- matrix(nrow=l, ncol=l)
#pheromone value
phvals <- matrix(nrow=l, ncol=l)
rownames(phvals) <- 1:l
colnames(phvals) <- 1:l
# only do the lower diagonal
for (u in 1:l){
  for (v in 1:l){
    if(u>v){
      p <- foo[u,]
      q <- foo[v,]
      dist[u,v] <- euc_dist(p,q)
      phvals[u,v] <- 1
    }
  }
}

#heuristic value, based on the euclidean distance. Larger value == closer distance
hvals <- 1-(dist/max(dist, na.rm=T))
rownames(hvals) <- 1:l
colnames(hvals) <- 1:l



# MOSTLY FUNCTION DECLARATION

# initialize n ants in random starting positions
init_ants <- function(n){
  antnames <- paste("Ant", seq(1,n))
  ants <- as.integer(runif(n, 1, l))
  names(ants) <- antnames
  ants
}

# take one step for each ant
#update_ants <- function(ants){lapply(ants, FUN=function(a){
update_ants <- function(ants){
  nants <- list()
  for(i in 1:length(ants)){
    a <- ants[i]
    aa <- a[[1]]
    la <- length(aa)
    # this check shouldnt be necesarry, but just in case you run more times than the number of cities.
    #if (la<l && names(a) == "Ant 1"){
    if (la<=l){
      u <- as.integer(tail(aa,1))
      #probabilities
      p <- rep(0, l)
      #feasible neighborhood
      N <- all[!all%in%a[[1]]]
      if(length(N)>1){
        # values (as by the supplied formula) for all feasible edges that include the node u
        vals <- unlist(lapply(N, FUN=function(v){
                          uu <- max(u,v)
                          vv <- min(u,v)
                          x <- (phvals[uu,vv]**alpha)*(hvals[uu,vv]**beta)
                          names(x) <- v
                          x
                          }))
        sumVals <- sum(vals, na.rm=T)
        if (length(sumVals)==0 || sumVals==0){
          choice <- N[[1]]
          print("WAT3")
          print(vals)
          break;
        }
        p <- vals/sumVals
        # generate the cumulative sum sequence
        cp <<- cumsum(p)

        choice <- NULL 
        while(length(choice)==0){
          #generate random value between 0 and 1
          x <- runif(1,0,1)
          # choose the corresponding city
          choice <- as.integer(head(names(cp[cp>x]), n=1))
        }
      }else if(length(N)==1){
        choice <- N[[1]]
      }else{
        break
      }

      nants[[i]] <- c(aa, choice)
      names(nants)[i] <- names(a)
    }else{
      nants[i] <- a
      names(nants)[i] <- names(a)
    }
  }
  return(nants)
}

release_ants <- function(ants){
    for(i in 1:(l-1)){
    ants <- update_ants(ants)
  }
  ants
}

get_scores <- function(ants){
  return(lapply(ants, FUN=function(a){
           d <- 0
           for (i in 1:(length(a)-1)){
             p <- foo[a[[i]],]
             q <- foo[a[[i+1]],]
             d <- d + euc_dist(p,q)
           }
           # and wrap
           p <- foo[a[[1]],]
           q <- foo[a[[length(a)]],]
           d <- d + euc_dist(p,q)
  }))
}
q <- 100
normalize_scores <- function(score){
  y <- lapply(X=score, FUN=function(x){
           (q/x)}
  )
  return(y)
}


drop_pheromones <- function(ants, scores, phvals){
  #first handle evaporation
  for(u in 1:l){
   for(v in 1:l){
     if(u>v){
      phvals[u,v] <- phvals[u,v]*(1-evap_rate) 
     }
   }
  }

  # then add score
  for (n in names(ants)){
    a <- ants[[n]]
    for (i in 1:(length(a)-1)){
      u <- a[[i]]
      v <- a[[i+1]]
      uu <- max(u,v)
      vv <- min(u,v)
      phvals[uu,vv] <- phvals[uu,vv]+scores[[n]]
    }
    # remember the wrap
    u <- a[[1]]
    v <- a[[length(a)]]
    uu <- max(u,v)
    vv <- min(u,v)
    phvals[uu,vv] <- phvals[uu,vv]+scores[[n]]
  }
  return(phvals)
}


# takes a normalized tour
partial_convergence_check <- function(ants){
  # take the second and last value of each tour
  # if the tours are sorted but, but some tours are reversed, 
  # the values will still occur in all tours at these positions
  subset <- unlist(lapply(ants, FUN=function(x){c(x[[2]],x[[l]])}))
  sum(subset==subset[[1]]) == n_ants
}


# Check ALL variables for 100% convergency
# takes a normalized tour
full_converge_check <- function(ants){
  val <- lapply(1:(n_ants-1), FUN=function(i){
      lapply((i+1):(n_ants), FUN=function(j){
               (all.equal(ants[[i]], ants[[j]])==T) || (all.equal(ants[[i]], rev(ants[[j]]))==T)
      }) 
    }
  )
  prod(unlist(val))==1
}


# PROGRAM EXECUTION

get_normalized_tours <- function(ants){
  tours <- list()
  for(i in 1:length(ants)){
    a <- ants[[i]]
    startpos <- which(a==1)
    if(startpos!=1){
      tours[[i]] <- c(a[(startpos:l)],a[(1:(startpos-1))])
      names(tours)[i] <- names(ants[i])
    }else{
      tours[[i]] <- a
      names(tours)[i] <- names(ants[i])
    }
  }
  return(tours)
}

best <- 99999999999999
i <- 1
while(i != max_its){
  print(paste(c("Iteration ",i),collapse=""))
  ants <- init_ants(n_ants)
  ants <- release_ants(ants)
  scores <- get_scores(ants)
  
  nbest <- min(c(best,unlist(scores)))
  if(best!=nbest){
    best_ind <- which(unlist(scores)==nbest)
    best_val <- ants[[best_ind]]
  }
  best <- nbest
  print(paste(c("Best so far: ",best),collapse=""))
  # check if any differences between ant routes

  ants2 <- get_normalized_tours(ants)
  if(partial_convergence_check(ants2)){
    print("PARTIAL OK")
    print(ants2)
       if(full_converge_check(ants2)){
         print("REACHED CONVERGENCE!")
         break
       }
  }

  scores <- normalize_scores(scores)
  as <- phvals[1:6,1:6]
  as[is.na(as)] <- 0
  print(as)
  print(c(max(phvals[,1],na.rm=T), max(phvals,na.rm=T)))
  phvals <- drop_pheromones(ants, scores, phvals)
  i <- i+1
}
if(i == max_its){
    print("REACHED MAX ITERATIONS")
}

print("ANTS")
ants <- get_normalized_tours(ants)
print(ants)
print("SCORES")
scores <- get_scores(ants)
#print(scores)
#print(paste(c("Best score last round:",min(c(best,unlist(scores)))),collapse=""))
#print(paste(c("Best score overall:",best),collapse=""))



png(filename="./out_offset.png", width=16, height=10, res=120, units="in")
cols <- sample(rainbow(n_ants),n_ants )
xoffset=1
yoffset=1
a <- ants[[1]]
plot(x=foo[a,1], y=foo[a,2], type="l", col=cols[1], ylim=c(min(foo[,2]),max(foo[,2])+xoffset*n_ants), xlim=c(min(foo[,1]),max(foo[,1])+yoffset*n_ants))
for(i in 2:n_ants){
  a <- ants[[i]]
  lines(x=foo[c(a,a[[1]]),1]+i*xoffset, y=foo[c(a,a[[1]]),2]+i*yoffset, type="l", col=cols[i])
}
for(i in 1:l){
  text(foo[i,], labels=paste(i))
}
dev.off()


png(filename="./out.png", width=16, height=10, res=120, units="in")
cols <- sample(rainbow(n_ants),n_ants )
xoffset=0
yoffset=0
plot(x=foo[c(best_val,best_val[[1]]),1], y=foo[c(best_val,best_val[[1]]),2], type="l", col=cols[1], ylim=c(min(foo[,2]),max(foo[,2])+xoffset*n_ants), xlim=c(min(foo[,1]),max(foo[,1])+yoffset*n_ants))

for(i in 1:l){
  text(foo[i,], labels=paste(i))
}
dev.off()

