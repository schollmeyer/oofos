##########################################
########			     	              ########
#######						                 #######
######     starshaped sets          ######
#######                            #######
########				                  ########
##########################################

#' Compute the canoncical (stylized) betweennes relation between objects of a formal context
#'
#' @description 'stylized_betweeness' computes the canoncical (stylizedb)
#' betweennes relation between objects of a formal context: given three
#' g,h,i of a formal context we say that h lies betweeen g and i iff
#' Psi({g}) cap Psi({i}) subseteq Psi({h}). Additionally we say in a stylized
#' manner that h lies betweeen g and i iff
#' Psi({g}) cap Psi({i}) subseteq Psi({h}) is almost true with the exception of
#' some attributes of Psi({g}) cap Psi({i}) that do not belong to Psi({h}) but
#' that are not so 'important' in the sense that there are not too much objects
#' that have all attributes in (Psi({g}) cap Psi({i})) -  Psi({h}) but that have
#' not all attributes in Psi({h}). Concretely the stylized betweennes is
#' quantified by the maximum of the weights 'attribute_weights' associated
#' to these attributes. (The weights could be defined for example as the cloumn
#' means of the underlying formal context.
#'
#' @param g is the index of object g w.r.t. the underlying context 'context'.
#'
#' @param h is the index of object h w.r.t. the underlying context 'context'.
#'
#' @param i is the index of object i w.r.t. the underlying context 'context'.
#' @param context is the underlying context.
#' @param attribute_weights is the weight vector for the attributes.
#' @return a ternary fuzzy relation given by an array of dimension
#' n x n x n where n is the number of rows of the context. Higher values
#' of an entry correspond to a larger degree of betweenness.
#' @export
stylized_betweeness <- function(g,h,i,context, attribute_weights){

  common_attributes <- which(g==1 & i==1)
  if(length(common_attributes)==0){return(1)}
  ans <- 1-max((1-h[common_attributes])*attribute_weights[common_attributes])
  return(ans)

}





compute_quotient_order <- function(I){

  #computes quotient order from a quasiorder
  index <- !duplicated(I)&!duplicated(t(I));return(I[index,index])}



compute_stylized_betweeness <- function(g,h,i,context, attribute_weights){

  common_attributes <- which(g==1 & i==1)
  if(length(common_attributes)==0){return(1)}
  ans <- 1-max((1-h[common_attributes])*attribute_weights[common_attributes])
  return(ans)

}



cut_incidence=function(I,width,interval=stats::quantile(unique(as.vector(I)),c(0.01,0.95))){
  vc <- width_hopcroft_karp(compute_quotient_order(compute_transitive_hull(I >= max(I))))
  if(vc <= width){return(I >= max(I))}
  #interval <<- interval
  f <- function(C,I){W=width_hopcroft_karp(compute_quotient_order(compute_transitive_hull(I >=C)))$width;return(W-width)}
  ans <- uniroot(f,interval=interval,I=I)
  return(compute_transitive_hull(I>=ans$root))}








starshaped_subgroup_discovery  <- function(stylized_betweenness,objective,vc_dim,params=list(Outputflag=0)){

  if (dim(stylized_betweenness)[1] != dim(stylized_betweenness)[2] | dim(stylized_betweenness)[1] != dim(stylized_betweenness)[3] | dim(stylized_betweenness)[2] != dim(stylized_betweenness)[3]){print("dimension mismatch")}
  m <- nrow(stylized_betweenness)
  model <- list(modelsense="max",obj=objective,lb=rep(0,m),ub=rep(1,m))
  solutions <- list()
  objvals <- rep(0,m)
  stars <- array(0,c(m,m))

  models=list()
  for(k in (1:m)){  ## quantify over all starcenters


    if (vc_dim == Inf) { incidence <- (stylized_betweenness[k,,] >= max(stylized_betweenness[k,,]))*1}
    else{ incidence <- cut_incidence(Z[k,,],vc_dim) }


    model <- model_from_qoset(t(incidence))# Z[k,,]))
    model$obj <- objective
    model$lb <- rep(0,m)
    model$ub <-rep(1,m)
    model$lb[k] <- 1              ## force centerpoint to be in the set
    model$modelsense <- "max"
    b <- gurobi(model,params=params)
    solutions[[k]] <- b
    objvals[k] <- b$objval
    stars[k,] <- b$x

    models[[k]] =model


  }
  i <- which.max(objvals)


  #I <<- stylized_betweeness[i,,]
  if (vc_dim == Inf) { incidence <- (stylized_betweenness[i,,] >= max(stylized_betweenness[k,,]))*1}
  else{incidence <- cut_incidence(stylized_betweenness[i,,],vc_dim)}


  model <- model_from_qoset(t(incidence))
  model$obj <- objective
  model$lb <- rep(0,m)
  model$ub <-rep(1,m)
  model$lb[k] <- 1              ## Sternmittelpunkt drinnen
  model$modelsense <- "max"
  b <- gurobi(model,params=params)
  #solutions[[k]] <- b
  #objvals[k] <- b$objval
  #stars[k,] <- b$x
  #incidence <- cut_incidence(Z[i,,],vc_dim)

  return(list(models=models,obj=objective,solutions=solutions,objvals=objvals,stars=stars,objval=objvals[i],star=stars[i,],center_id =i,fuzzy_incidence=stylized_betweenness[i,,] , incidence = incidence,model=model) )}

plot_stars <- function(starshaped_result,distance_function){

  i <- which(starshaped_result$star==1)

  j <- compute_maximal_elements(starshaped_result$incidence[i,i])


  i[j]
  print(i[j])





}

plot_corder <- function(corder,main=""){
  m <- nrow(corder)
  plot(as.relation(corder[(1:m),(1:m)]),main=main)

}

starshaped_subgroup_discovery_recompute <- function(models,objective){



  ans <- -Inf
  for(k in (1: length(models))){
    M <- models[[k]]
    M$obj <- objective
    ans <- max(ans,gurobi(M,param=list(outputflag=0))$objval)
  }

  return(ans)}


starshaped_subgroup_discovery_h0 <- function(models,params=list(outputflag=0)){

  v <- sample(models[[1]]$obj)

  ans <- -Inf
  for(k in (1: length(models))){
    M <- models[[k]]
    M$obj <- v
    ans <- max(ans,gurobi(M,params=params)$objval)
  }

  return(ans)}



model_from_qoset <- function(Q){## constructs linear program for the optimization over all upsets of a quasiordered set Q

  QQ <- tr(Q)
  m  <- sum(QQ)
  n  <- dim(QQ)[1]
  A  <- array(as.integer(0),c(m,n))
  t  <- 1
  sense <- rep("<=",m)
  for(k in (1:(n-1))){
    for(l in ((k+1):n)){
      if(QQ[k,l]==1 &QQ[l,k]==0){
        A[t,k]=1;A[t,l]=-1
        t=t+1
      }

      if(QQ[k,l]==0 &QQ[l,k]==1){
        A[t,k]=-1;A[t,l]=1
        t=t+1
      }

      if(QQ[k,l]==1 &QQ[l,k]==1){
        A[t,k]=1;A[t,l]=-1
        sense[t]="="
        t=t+1
      }
    }
  }
  t <- t-1
  ans <- list(A=A[(1:t),],rhs=rep(0,t),sense=sense[(1:t)],lb=rep(0,n),ub=rep(1,n))
  return(ans)}



classification.with.stylized.betweeness=function(x.train,y.train,x.test,y.test,stylizedBetweenness=sb1,p,VCDim,params=list(outputFlag=0,presolve=0,threads=1),VCcut=TRUE,interval){

  print(Sys.time())
  m.train <- dim(x.train)[1]
  m.test  <- dim(x.test)[1]
  m       <- m.train+m.test
  labels  <- levels(y.train)
  X       <- rbind(x.train,x.test)
  Y       <- c(y.train,y.test)
  II      <- array(0,c(m,m,m))
  AP     <- array(0,c(m.test,2*m.train));AM=AP
  bnsplus <- bnsminus <- ansplus <- ansminus <- rep(0,m.train)
  ansplussol=list()
  ansminussol=list()
  bnsplussol=list()
  bnsminussol=list()
  sol=list()
  sol2=list()
  t=1
  sols=list()

  for(kkk in (1:length(stylizedBetweenness))){
    kkk<<-kkk
    SB <- (stylizedBetweenness[kkk])
    SB <- SB[[1]]
    #Print(SB)
    B       <- SB(X,p=p)
    FI      <- pmin(B$A,B$B)
    FII <<- FI



    for(k in (1:m)){

      if(VCcut){J=incidence.cut(FI[k,,],width=VCDim[kkk],interval=interval)}#c(sort(unique(as.vector(FI[k,,])))[2],1))}
      else{J=incidence.cut2(FI[k,,],width=VCDim[kkk],interval=c(sort(unique(as.vector(FI[k,,])))[2],1))}
      II[k,,]=pmax(II[k,,],(J))
    }
  }


  II <<- II
  VCDIMS=rep(0,m)
  for(kkk in (1:m)){

    VCDIMS[kkk]=(width.hopcroft.karp(transitive.hull(II[kkk,,]))$width)}
  Print(table(VCDIMS))
  II=CUT(II,EPS,m.train)

  II <<- II
  print(Sys.time())

  for (ind in (1:m.test)){
    index <- c((1:m.train), ind+m.train)
    #print(Sys.time())
    for(k in (1:m.train)){  ## quant ueber alle sternmittelpunkte## hier auch m moeglich
      M <- model_from_qoset(II[k,index,index])
      M$lb[k] <- 1              ## Sternmittelpunkt drinnen

      #1
      v <- rep(0,m.train+1);v[(1:m.train)]=probabilityDifferences(y.train,label=labels[1])
      M$modelsense <- "max"
      M$lb[m.train+1] <- 1   ## testpunkt drinnen
      M$obj <- v  #+    CC*rep(-1,length(v))###fuer  vortrag
      b <- gurobi(M,params=params)
      ansplus[k] <- b$objval
      ansplussol[[k]]=b$x

      sols[[t]]=b$x
      t=t+1

      ##fuer vortrag:

      #if(k==K){
      #b <<- b
      #}


      #2

      v <- rep(0,m.train+1);v[(1:m.train)]=probabilityDifferences(y.train,label=labels[2])
      M$obj <- v
      b <- gurobi(M,params=params)
      ansminus[k] <- b$objval
      ansminussol[[k]]=b$x

      sols[[t]]=b$x
      t=t+1


      #3
      v <- rep(0,m.train+1);v[(1:m.train)]=probabilityDifferences(y.train,label=labels[1])
      M$obj <- v
      M$lb[m.train+1] <- 0;M$ub[m.train+1] <- 0  ### testpunkt nicht drinnen
      M$modelsense <- "min"
      b <- gurobi(M,params=params)
      if(is.null(b$objval)){print("warning: no feasible solution")}

      if(!is.null(b$objval)){bnsplus[k]=-b$objval}
      bnsplussol[[k]]=b$x

      sols[[t]]=b$x
      t=t+1



      #4
      v <- rep(0,m.train+1);v[(1:m.train)]=probabilityDifferences(y.train,label=labels[2])
      M$obj <- v
      b <- gurobi(M,params=params)
      if(is.null(b$objval)){print("warning: no feasible solution")}
      if(!is.null(b$objval)){bnsminus[k]=-b$objval}
      bnsminussol[[k]]=b$x

      sols[[t]]=b$x
      t=t+1
    }

    AP[ind,] <- c(ansplus,bnsplus)
    AM[ind,] <- c(ansminus,bnsminus)

  }
  predictions <- factor(rep(y.train[1],m.test),level=labels)
  unique=rep(0,m.test)
  for(k in (1:m.test)){
    if(AP[k,] %glex% AM[k,]){predictions[k] <- labels[1];sol[[k]]=ansplussol[[k]];sol2[[k]]=bnsplussol[[k]]}
    else{predictions[k] <- labels[2];sol[[k]]=ansminussol[[k]];sol2[[k]]=bnsminussol[[k]]}
    if(max(AP[k,]) > max(AM[k,]) | max(AP[k,]) < max(AM[k,])){unique[k]=1}
  }

  return(list(AP=AP,AM=AM,predictions=predictions,errorProb=mean(y.test!=predictions),unique=unique, unique.prop=mean(unique),betweenProb=betweennessProb(II[(1:m.train),(1:m.train),(1:m.train)],y.train),sol=sol,sol2=sol2,sols=sols))}


cwsb=classification.with.stylized.betweeness


###### stilisierte Zwischenrelationen


sb1=function(XR,p){
  m=dim(XR)[1]
  n=dim(XR)[2]
  A=array(as.integer(0),c(m,m,m))
  B=A
  C=A
  D=as.matrix(dist(XR,method="minkowski",upper=TRUE,p=p))
  MAXD=max(D)#Q=quantile(D,QQ)

  for(K in (1:m)){#(1:(m-2))){
    for(L in (1:m)){#((K+1):(m-1))){
      for(M in (1:m)){#((L+1):m)){
        #A[K,L,M]=max(0,D[L,K]-D[M,K],D[L,M]-D[K,M])
        A[K,L,M]=sum(D[K,M])/(sum(D[K,L]+D[L,M]))
        B[K,L,M]=(D[K,L]<D[K,M])#*D[K,M]/MAXD#<=Q

        #temp=sign((XR[L,]-XR[K,])*(XR[M,]-XR[K,]))
        #temp2=abs(XR[L,]-XR[K,]) < abs(XR[M,]-XR[K,])
        #i=which(temp>=0 &temp2)
        #A[K,L,M]=length(i)
        #B[K,L,M]=sum(abs(temp[i]))

      }}}
  for(K in (1:m)){
    diag(A[K,,])=1
    diag(B[K,,])=1}
  return(list(A=A,B=B))}

