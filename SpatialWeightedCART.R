
rm(list=ls())

#remotes::install_github("42n4/dmr.util" , force = TRUE)

library(tictoc)
library(geoR)
library(deldir)
library(rpart)
library(MASS)
library(spatstat)
library(entropy)
library(sp)
library(gstat)
library(distances)
library(automap)
library(lsei)
library(ks)
library(testit)
library(KernelKnn)
library(spdep)


library(foreach)
library(doParallel)
library(parallel)

#install.packages("remotes")
#library(remotes)
#remotes::install_github("42n4/dmr.util" , force = TRUE)
library(dmr.util)

#================== Neyman-Scott function =======================


neyman = function(r , sig2 , rho , nc){

  nclust = function(x0 , y0 , radius , nc) {
    return(runifdisc(nc , radius , centre = c(x0 , y0)))
  }
  
  w = owin(c(-2,2) , c(-2,2))
  S1 = rNeymanScott(kappa=2/16, 0, rcluster=nclust, win=w
                    , radius=rho, nc=nc, saveparents=TRUE)
  
  N1 = npoints(S1)
  S2 = rpoispp((200-N1)/16 , win=w)# 
  N2 = npoints(S2)
  
  
  S = coords(S1)
  SS = coords(S2)
  X = c(S$x,SS$x) ; Y = c(S$y,SS$y)
  
  
  e1 = grf(N1+N2 , cov.model = "exponential", cov.pars = c(sig2,r))
  e2 = grf(N1+N2 , cov.model = "exponential", cov.pars = c(sig2,r)) 
  e3 = grf(N1+N2 , cov.model = "exponential", cov.pars = c(sig2,r)) 
  e4 = grf(N1+N2 , cov.model = "exponential", cov.pars = c(sig2,r)) 
  
  X1 = Y - X + e1$data
  X2 = Y + X + e2$data
  X3 = X + e3$data
  X4 = Y + e4$data
  
  
 
  Z = c()
  for(i in 1:(N1+N2)){
    if (((X[i]<Y[i])&(Y[i]<(-X[i]))&(X[i]<(-1))) | (((-X[i])<Y[i])&(Y[i]<X[i])&(X[i]>1))){
      Z[i]=1
    }else  if(((Y[i]<(-X[i]))&(Y[i]<X[i])&(Y[i]<(-1))) | ((Y[i]>(-X[i]))&(Y[i]>X[i])&(Y[i]>1))){
      Z[i]=2
      
    }else  if((Y[i]<(-X[i]))&(X[i]>(-1))&(Y[i]>(-1))){
      Z[i]=3
    }else  if((Y[i]>(-X[i]))&(X[i]<1)&(Y[i]<1)){
      Z[i]=4
    }else {Z[i]=0}}
  
  
  data = data.frame(X,Y,X1,X2,X3,X4,Z)
  data
}


#========================== Considered weights ==============================

#-------------------- Uniform --------------------------

weight.1 <- function(data){ rep(1/nrow(data) , nrow(data))}

#-------------------- Voronoi -------------------------

weight.2 <- function(data){
  if(nrow(data)>1){
    D <- deldir(data$X , data$Y , rw=c(min(data$X),max(data$X),min(data$Y),max(data$Y)))#,plotit=TRUE)
    D$summary$dir.wts # or D[[3]]# dirichlet(voronoi)weights
    ### D$summary$del.wts
  }else{1}
}

#--------------------- Kernel --------------------------

# H=Hpi(x=data[,1:2])  # default
# H=Hns(x=data[,1:2])
# H=Hscv(x=data[,1:2])

weight.3 <- function(data){ 
  wb <- as.matrix(data[,1:2])
  fhat <- 1/kde(x=wb , eval.points=wb,  H=Hpi(x=wb) )$estimate
  fhat/sum(fhat)
}

#--------------- Class-based kernel --------------------

weight.4 <- function(data){ 
  hhh = c()
  for(i in 1:length(cprobs)){
    if(length(which(data[,class]==i))>2){
      wb <- as.matrix(data[which(data[,class]==i),1:2])
      fhat <- 1/kde(x=wb , eval.points=wb, H=Hpi(x=wb) )$estimate
    }else{fhat <- 0}
    hhh[which(data[,class]==i)] <- fhat
  }
  hhh/sum(hhh)
}


#-------------------- Kriging (Case I) ---------------------------

krige.weights1 <- function(data){ 
  warn <<- c()
  coef.krige = matrix(NA,length(clabs),length(data[,class])) ##--> nrow(data)
  if(length(data[,class])<=3){
    for(i in 1:length(clabs)){
      coef.krige[i,]= rep(1/length(data[,class]),length(data[,class]))
    }}else{
      dataa = data
      coordinates(dataa) = ~X+Y
      sp1 = seq(min(dataa$X),max(dataa$X),length=10)
      sp2 = seq(min(dataa$Y),max(dataa$Y),length=10)
      sp.grid = expand.grid(sp1,sp2)# data.frame
      
      distance = distances(data[,1:2]) # matrix
      dist = matrix(NA , nn <- length(data[,1]) , mm <- length(sp.grid[,1]))
      for(j in 1:mm){ dist[,j] = sqrt((data[,1]-sp.grid[j,1])^2+(data[,2]-sp.grid[j,2])^2) } # dist[,j] vector ???
      
      for(i in 1:length(clabs)){   
        if(length(which(data$Z==i))< 5){ # 0 or 1
          coef.krige[i,]=rep(1/length(data$Z),length(data$Z))
          warn[i] <<- FALSE
        }else{
          experimental_variogram = variogram(I(Z==i) ~ 1, dataa)
          
      
          if(all(experimental_variogram$gamma == 0) || is.null(experimental_variogram)){
            coef.krige[i,]=rep(1/length(data$Z),length(data$Z))
            warn[i] <<- FALSE
          }else{
            warn[i] <<- suppressWarnings(has_warning(variogram <- fit.variogram( experimental_variogram, vgm("Exp") ) )) #,"Gau"
            if(warn[i]==TRUE){
              coef.krige[i,]=rep(1/length(data$Z),length(data$Z))
            }else{
              vario.fit1 = variogramLine(variogram, dist_vector = as.matrix(distance))#, covariance=TRUE)
              vario.fit2 = c()
              for(j in 1:nn ){
                vario.vector = (variogramLine(variogram, dist_vector = dist[j,] ))[[2]]#, covariance=TRUE ))[[2]] 
                vario.fit2[j] = sum(vario.vector)/mm #### or nn ??????
              }
              if(all(vario.fit1==0) & all(vario.fit2==0)){
                coef.krige[i,]=rep(1/length(data$Z),length(data$Z))}else{
                  coef.krige[i,] = pnnqp(-2*vario.fit1 , 2*vario.fit2 , k=0, sum=1, tol=1e-20)$x # zeros found exactly
                }
            }
          }
        }
      }
    }
  kw <- c()
  for(i in 1:nrow(data)){
    kw[i] <- coef.krige[data[i,class],i]  # proportions of misclassified points
  }
  kw/sum(kw)
}


#-------------------- Kriging (Case II) ---------------------------

krige.weights2 <- function(data){ 
  warn <<- c()
  coef.krige = matrix(NA,length(clabs),length(data[,class])) 
  if(length(data[,class])<=3){
    for(i in 1:length(clabs)){
      coef.krige[i,]= rep(1/length(data[,class]),length(data[,class]))
    }}else{
      dataa = data
      coordinates(dataa) = ~X+Y
      sp1 = seq(min(dataa$X),max(dataa$X),length=10)
      sp2 = seq(min(dataa$Y),max(dataa$Y),length=10)
      sp.grid = expand.grid(sp1,sp2)# data.frame
      
      distance = distances(data[,1:2]) # matrix
      dist = matrix(NA , nn <- length(data[,1]) , mm <- length(sp.grid[,1]))
      for(j in 1:mm){ dist[,j] = sqrt((data[,1]-sp.grid[j,1])^2+(data[,2]-sp.grid[j,2])^2) } # dist[,j] vector ???
      
      for(i in 1:length(clabs)){   
        if(length(which(data$Z==i))< 5){ # 0 or 1
          coef.krige[i,]=rep(1/length(data$Z),length(data$Z))
          warn[i] <<- FALSE
        }else{
          experimental_variogram = variogram(I(Z==i) ~ 1, dataa)
          
 
          if(all(experimental_variogram$gamma == 0) || is.null(experimental_variogram)){
            coef.krige[i,]=rep(1/length(data$Z),length(data$Z))
            warn[i] <<- FALSE
          }else{
            warn[i] <<- suppressWarnings(has_warning(variogram <- fit.variogram( experimental_variogram, vgm("Exp") ) )) #,"Gau"
            if(warn[i]==TRUE){
              coef.krige[i,]=rep(1/length(data$Z),length(data$Z))
            }else{
              vario.fit1 = variogramLine(variogram, dist_vector = as.matrix(distance))#, covariance=TRUE)
              vario.fit2 = c()
              for(j in 1:nn ){
                vario.vector = (variogramLine(variogram, dist_vector = dist[j,] ))[[2]]#, covariance=TRUE ))[[2]] 
                vario.fit2[j] = sum(vario.vector)/mm #### or nn ??????
              }
              if(all(vario.fit1==0) & all(vario.fit2==0)){
                coef.krige[i,]=rep(1/length(data$Z),length(data$Z))}else{
                  ident = c(rep(1,nn))
                  coef.krige[i,] = t( (vario.fit2 + ident * (as.vector(1-t(ident) %*% solve(vario.fit1) %*% vario.fit2) / 
                                                               as.vector(t(ident) %*% solve(vario.fit1) %*% ident)
                  ) ) ) %*% solve(vario.fit1)
                  m = -((1-t(ident) %*% solve(vario.fit1) %*% vario.fit2)/(t(ident) %*% solve(vario.fit1) %*% ident))
                  
                }
            }
          }
        }
      }
    }
  
  coef.krige[coef.krige<0]=0 

  kw <- c()
  for(i in 1:nrow(data)){
    kw[i] <- coef.krige[data[i,class],i]  # proportions of misclassified points
  }
  kw/sum(kw)
  
}


#================ Impurity functions =================

#---------------------- Entropy ----------------------

E = function(data){
  p = prob.treee(data)
  if(all(is.nan(p))){NaN}else{
    for(i in 1: length(p)){if(p[i]==0) p[i]=1 } ##sapply or ifelse
    -sum(p*log(p,2))
  }
}


#------------------- Spatial Entropy -----------------

diver.fun <- function(data,lambda,betta){
  #  clabs = cprobs = NULL
  #  clabs <<- levels(factor(data[,ncol(data)]))
  if(length(data[,ncol(data)]) > 3 ){
    index <- lapply(as.numeric(clabs) , function(i) which(data[,class]==i)) ## class = ncol(data)-1
    distance = distances(data[,1:2], normalize = "mahalanobize")
    
    diver <- c()
    for(i in as.numeric(clabs)){
      sum.int <- 0 ; sum.ext <- 0
      if (length(index[[i]]) > 1){
        for(j in 1:length(index[[i]])){
          for(k in 1:length(index[[i]])){
            sum.int <- sum.int + distance[index[[i]][j],index[[i]][k]]
          }
        }
        dint <- sum.int / (length(index[[i]])*(length(index[[i]])-1))
      }else{
        dint <- lambda
        }
      
      if(length(index[[i]])== 0 || length(index[[i]])== nrow(data)){
        dext <- betta
        }else{
          
          for(j in 1:length(index[[i]])){
            for(d in as.numeric(clabs)){
              if(i!=d & length(index[[d]])!=0){
                for(k in 1:length(index[[d]])){
                  sum.ext <- sum.ext + distance[index[[i]][j],index[[d]][k]]
                }
              }
            }
          }
          
          dext <- sum.ext / (length(index[[i]])*(nrow(data)-length(index[[i]])))
        }
      diver[i] <- dint / dext
    }
  }else{
    diver <- rep(0.0001,length(cprobs))
  }
  diver
}


lambda <- 0.01
betta <- 100

spat.E <- function(data){
  p <- prob.treee(data)
  if(all(is.nan(p))){NaN}else{
    diver <- diver.fun(data ,lambda ,betta )
    for(i in 1: length(p)) {if(p[i]==0) p[i]=1} ##sapply
    -sum(diver*p*log(p,2))
  }
}



#=============================== Decision tree function ==================================

grow.dectree = function(formula, data, weight, imp, maxprob, minsplit, maxdepth, cp)
{
  
  init = function()
  {
    clabs <<- levels(factor(data[,class]))
   
    treee <<- data.frame(node=1, attribute=NA, value=NA, index=NA, bestIG=NA,
                         class=NA, count=NA, nodemiss=NA, nodeflag=NA,
                         "names<-"(rep(list(NA), length(clabs)),
                                   paste("p", clabs, sep=".")))
    cprobs <<- (ncol(treee)-length(clabs)+1):ncol(treee) # class probability columns
    nodemap <<- rep(1, nrow(data))
    fullnodemap <<- list()
    n <<- 1
    
    if(weight == "krige.weights1"){ weight1 <<- krige.weights1(data)}
    if(weight == "krige.weights2"){ weight1 <<- krige.weights2(data)}    
    if(weight=="uniform"){ weight1 <<- weight.1(data)}
    if(weight=="voronoi"){ weight1 <<- weight.2(data)}
    if(weight=="kernel"){ weight1 <<- weight.3(data)}
    if(weight=="kernelclass"){ weight1 <<- weight.4(data)}
    if(weight=="kernelknn"){ weight1 <<- weight.5(data )}
    
  }
  
  ##
  
  next.node = function(n){
    if (any(opn <- treee$node>n))
      min(treee$node[opn])
    else Inf }
  
  ##
  
  class.info = function(n){
    treee$count[treee$node == n] <<- length(data[nodemap==n, class])   
    treee[treee$node == n, cprobs] <<- prob.treee(ccc[nodemap==n, ])
    treee$class[treee$node==n] <<- which.max(treee[treee$node == n, cprobs])
    treee$nodemiss[treee$node==n] <<- (1-treee[treee$node == n, cprobs[ treee$class[treee$node==n] ]])
    
  }
  
  ##

  stop.criteria = function(n){
    n >= 2^maxdepth || treee$count[treee$node==n] <= minsplit ||
      max(treee[treee$node==n,cprobs]) > maxprob
  }
  
  ##
  
  split1 = function(data , h){
    n0 = length(data[,class])
    
    IG = c() 
    for(i in 1:nrow(data)){
      v1 = split(data, factor(data[,h] <= data[i,h]))# data frame that ordered by X1<X1[i]
      if(length(v1$'TRUE'[,class]) ==0 || length(v1$'FALSE'[,class]) ==0 ){
        IG[i] = Inf }else{
          
          El = imp(v1$'TRUE') #imp 
          Er = imp(v1$'FALSE') #imp
          
          IG[i] = (length(v1$'TRUE'[,class])/n0)*El+(length(v1$'FALSE'[,class])/n0)*Er
        }
    }
    IG
  }
  
  ##
  
  split.select = function(n){
    splits = data.frame()
    for(h in attributes){ 
      ss = split1(ccc[nodemap==n,],h) 
      splits = rbind(splits , data.frame(attribute = h , value = data[nodemap==n,][a<-which.min(ss),h]
                                         , index=ccc[nodemap==n,][a,"indx"], IG = ss[which.min(ss)]))
    } 
    if ((best.eval <- min(splits$IG))<Inf)
      treee[treee$node==n,2:4] <<- splits[which.min(splits$IG),1:3]
    treee$bestIG[treee$node==n] <<- best.eval
    best.eval
  }
  
  ##
  
  split.apply = function(n){
    treee <<- rbind(treee,
                    data.frame(node=(2*n):(2*n+1),
                               attribute=NA, value=NA,  index=NA, bestIG=NA,
                               class=NA, count=NA, nodemiss=NA, nodeflag=NA,
                               "names<-"(rep(list(NA), length(clabs)),
                                         paste("p", clabs, sep="."))))
    
    av = data[which(attributes==treee$attribute[treee$node==n])+2]
    cond = !is.na(av) & av <= as.numeric(treee$value[treee$node==n])
    
    nodemap[nodemap==n & cond] <<- 2*n
    nodemap[nodemap==n & !cond] <<- 2*n+1
    fullnodemap[[n]] <<- nodemap
    fullnodemap
    
  }
  
  ##-------------------------------------------

  treee = n =  nodemap <<- NULL
  class <- y.var(formula)
  attributes <- x.vars(formula, data)
  clabs = cprobs = NULL
  ccc = cbind(data , indx = c(1:nrow(data)))
  
  start_time = Sys.time()
  init()
  while (is.finite(n)){
    class.info(n)
    treee$nodeflag[treee$node==n] <- ifelse(!stop.criteria(n),TRUE,FALSE)
    if (!stop.criteria(n))
      if (split.select(n)<Inf)
        split.apply(n)
    n <- next.node(n)
  }
  
  Time <<- difftime(Sys.time(), start_time, units = "secs")[[1]]
  treee
  

}

#----------------------------------

weighted.p = function(data ){
  sapply(as.numeric(clabs) ,function(i)
    sum(weight1[data$indx]*I(data[,class]==i))/sum(weight1[data$indx]))#treee$cprobs
}


prob.treee = weighted.p
mytree = grow.dectree(formula <- Z~X1+X2+X3+X4 , data <- mydata , weight <- "uniform" , 
                      imp <- E , maxprob <- 0.95 , minsplit <- 10 , maxdepth <- 8 , cp <- 0.01)

#===================== final predict ==========================

predict.tree <- function(tree, data)
{
  #  treee=tree
  clabs = cprobs = NULL
  clabs <- levels(factor(data[,class]))
  cprobs <- (ncol(tree)-length(clabs)+1):ncol(tree) # class probability columns
  
  stop.criteria <- function(n){
    n>=2^maxdepth || tree$count[tree$node==n]<=minsplit ||
      max(tree[tree$node==n,cprobs])>maxprob
  }
  
  
  descend <- function(n)
  {
    if (!is.na(tree$attribute[tree$node==n]) &  tree$nodeflag[tree$node==n]) #& !stop.criteria(n)) # unless reached a leaf
    {
      av <- data[which(attributes==tree$attribute[tree$node==n])+2]
      cond <- !is.na(av) & av<=tree$value[tree$node==n]
      
      nodemap11[nodemap11==n & cond] <<- 2*n
      nodemap11[nodemap11==n & !cond] <<- 2*n+1
      
      descend(2*n)
      descend(2*n+1)
    }
  }
  
  nodemap11 <<- rep(1, nrow(data))
  descend(1)
  tree$class[match(nodemap11, tree$node)]
}


#====================== Miss classification ==========================
  
missclass.tree11 <- function(data, pred, weight)
{
  if(weight == "krige.weights1"){ weight1 <- krige.weights1(data)}
  if(weight == "krige.weights2"){ weight1 <- krige.weights2(data)}  
  if(weight == "uniform"){ weight1 <- weight.1(data)}
  if(weight == "voronoi"){ weight1 <- weight.2(data)}
  if(weight == "kernel"){ weight1 <<- weight.3(data)}
  if(weight == "kernelclass"){ weight1 <- weight.4(data)}
  if(weight == "kernelknn"){ weight1 <- weight.5(data )}
  
  
  sum((data[,class] != pred) * weight1)  # proportions of misclassified points
  
}


#========================= Prune function ===========================

prune.tree <- function(tree, data, cp){
  TND <- tree[!tree$nodeflag,"node"] # terminal nodes
  IND <- tree$node[!tree$node %in% TND][-1] # internal nodes
  
  branch.map <- function(){
    TND <<- tree[!tree$nodeflag,"node"] # terminal nodes
    IND <<- tree$node[!tree$node %in% TND][-1] # internal nodes
    
    fullnodemap.branch <<- list() ; internal <<- list()
    terminal<<-list() ; terminal[[1]] <<- TND
    for(n in tree$node[-1]){
      #  MM <<- data[fullnodemap[[floor(n/2)]]==n,]
      nodemap22 <<- rep(n, nrow(MM<<-data[fullnodemap[[floor(n/2)]]==n,]))
      h <<- 1
      inter <<- c()
      descend(n)
      fullnodemap.branch[[n]] <<- nodemap22 
      terminal[[n]] <<- sort(unique(fullnodemap.branch[[n]]))
      internal[[n]] <<- inter[-1]
    }
  }
  
  # branch.map()  
  ####
  stop.criteria <- function(n){
    n>=2^maxdepth || tree$count[tree$node==n]<=minsplit ||
      max(tree[tree$node==n,cprobs])>maxprob
  }
  
  #### R(t) = r(t)*p(t)
  nodemiss.t <- function(t1,t2){
    tree$nodemiss[tree$node==t1]*(tree$count[tree$node==t1]/tree$count[tree$node==floor(t2)])
  }
  
  #### R(T_t)=sum(r(t)*p(t))
  branchmiss.t <- function(t){
    k <- 0
    for(m in terminal[[t]]){
      k <- k+nodemiss.t(m,t)
    }
    k
  }
  
 
  # Nodemiss for nonterminal nodes ==> internal nodes
  descend <- function(n)
  {
    if( tree$nodeflag[tree$node==n]){
      inter[h] <<- n
      h <<- h+1
    }
    if ( tree$nodeflag[tree$node==n] ) # unless reached a leaf ##!is.na(tree$attribute[tree$node==n]) &
    {
      av <- MM[which(attributes==tree$attribute[tree$node==n])+2]
      cond <- !is.na(av) & av<=tree$value[tree$node==n]
      
      
      nodemap22[nodemap22==n & cond] <<- 2*n
      nodemap22[nodemap22==n & !cond] <<- 2*n+1
      
      descend(2*n)
      descend(2*n+1)
    }
  }
  
  
  ##
  tree.link <- function(IND1 , terminal1){
    l <- c()
    for(i in 1:length(IND1)){
      R.T <- branchmiss.t(IND1[i]) 
      R.t <- tree$nodemiss[tree$node == IND1[i]] 
      
      l[i] <- ( (R.t-R.T)/(length(terminal1[[IND1[i]]])-1) )
    }
    l
  }
  
  
  for(i in 1:(length(TND)-1)){
    if(TND[i]%%2==0 & TND[i+1]==TND[i]+1){
      a <- nodemiss.t(TND[i]/2,TND[i]/2)
      b <- nodemiss.t(TND[i],TND[i]/2)
      c <- nodemiss.t(TND[i+1],TND[i+1]/2)
      
      if( ( a-(b+c) ) == 0 ){
        tree <- tree[!tree$node==TND[i],]
        tree <- tree[!tree$node==TND[i+1],]
        tree$nodeflag[tree$node == TND[i]/2] <- FALSE
        
      }
    }
  }
 
  
  branch.map()
  
  while( cp >= (link <- tree.link(IND , terminal))[w <<- which.min(link)] ){
    for(t in sort(c(internal[[IND[w]]], terminal[[IND[w]]])) ){
      tree <- tree[!tree$node == t,]
    }
    tree$nodeflag[tree$node == IND[w]] <- FALSE
    branch.map()
    if(length(IND)==0){ break}

  }
  
  tree
}

##=========================== pruned tree ==========================

valid.tree <- function(train , test , formula, weight, imp, prob = c("weighted , kriging") ){ 
  library(dmr.util)
  #library(geoR)
  library(deldir)
  library(rpart)
  library(MASS)
  library(spatstat)
  library(entropy)
  library(sp)
  library(gstat)
  library(distances)
  library(automap)
  library(lsei)
  library(ks)
  library(testit)
  library(KernelKnn)
  library(spdep)
  
  if(prob == "weighted"){ prob.treee <<- weighted.p}
  if(prob == "kriging"){  prob.treee <<- krige.p}
  
  data <- train
  
  treee = n =  nodemap <<- NULL
  class <- y.var(formula)
  attributes <- x.vars(formula, train)
  clabs = cprobs = NULL
  ccc <- cbind(train , indx = c(1:nrow(train)))
  
  init()
  
  mytree <<- grow.dectree(formula ,data <- train , weight , imp , maxprob <<- 0.95 , minsplit <<- 10 , maxdepth <<- 8)
  prunedtree <<- prune.tree( mytree, train, cp <- 0.01)
  
  predictt <<- predict.tree(prunedtree, test)
  mis<-missclass.tree11(test,predictt,weight)
  data.frame(mis,Time)
}


##======================================================
##==============================================================================
valid.tree1 <- function(train , test , formula, weight, imp, prob = c("weighted , kriging") ){ 
  library(dmr.util)
  library(geoR)
  library(deldir)
  library(rpart)
  library(MASS)
  library(spatstat)
  library(entropy)
  library(sp)
  library(gstat)
  library(distances)
  library(automap)
  library(lsei)
  library(ks)
  library(testit)
  library(KernelKnn)
  library(spdep)
  
  if(prob == "weighted"){ prob.treee <<- weighted.p}
  if(prob == "kriging"){  prob.treee <<- krige.p}
  
   treee = n =  nodemap <<- NULL
   class <- y.var(formula)
   attributes <- x.vars(formula, data)
   clabs = cprobs = NULL
   ccc = cbind(data , indx = c(1:nrow(data)))
  
   init()
  
  mytree <<- grow.dectree(formula ,data <- train , weight , imp , maxprob <<- 0.95 , minsplit <<- 10 , maxdepth <<- 8)

  
  predictt <<- predict.tree(mytree , test)
  mis <- missclass.tree11(test,predictt,weight)
  data.frame(mis,Time)
}

