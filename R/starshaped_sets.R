##########################################
########			     	              ########
#######						                 #######
######     starshaped sets          ######
#######                            #######
########				                  ########
##########################################

get_random_ternary_relation <- function(n,prob,reflexive=FALSE,transitive=FALSE){

  result <- (stats::runif(n^3)>=prob)*1
  dim(result) <- rep(n,3)
  for(k in (1:n)){
   if(reflexive){diag(result[k,,])=1}
    if(transitive){result[k,,] <- compute_transitive_hull(result[k,,])}
  }
  return(result)

}

get_betweenness_from_p_order <- function(p_order){
  n_rows <- nrow(p_order)
  result <- array(0,c(rep(n_rows,3)))
  for(k in seq_len(n_rows)){
    for(l in seq_len(n_rows)){
      for(m in seq_len(n_rows)){
        result[k,l,m] <- (p_order[k,l] & p_order[l,m]) |
                          (p_order[m,l] & p_order[l,k])
      }
    }
  }
return(result)

}

check_if_starshaped <- function(set, ternary_relation){
  if(all(set==0)){return(FALSE)}
  n_rows <- nrow(ternary_relation)
  for(k in which(set==1)){

    if(check_if_center_point(k,set,ternary_relation)){return(TRUE)}
  }
return(FALSE)}

check_if_center_point <- function(point, set,ternary_relation){
  n_rows <- nrow(ternary_relation)
  for(k in which(set==1)){

    indexs <- which(ternary_relation[point,,k]==1)

    if(any(set[indexs]==0)){return(FALSE)}

  }

  return(TRUE)


}

#' Compute the canoncical (stylized) betweennes relation between objects of a formal context
#'
#' @description 'compute_stylized_betweenness' computes the canonical (stylized)
#' betweennes relation between objects of a formal context: Given three
#' objects g,h,i of a formal context we say that h lies betweeen g and i iff
#'
#'
#' Psi({g}) cap Psi({i}) subseteq Psi({h}).
#'
#'
#' Additionally we say in a stylized
#' manner that h lies betweeen g and i iff
#' Psi({g}) cap Psi({i}) subseteq Psi({h}) is almost true with the exception of
#' some attributes of Psi({g}) cap Psi({i}) that do not belong to Psi({h}) but
#' that are not so 'important' in the sense that there are not too much objects
#' that have all attributes in (Psi({g}) cap Psi({i})) -  Psi({h}) but that have
#' not all attributes in Psi({h}). Concretely the stylized betweennes is
#' quantified by the maximum of the weights 'attribute_weights' associated
#' to these attributes. (The weights could be defined for example as the cloumn
#' means of the underlying formal context).
#'
#' @param g is the object g represented by a 0-1 vector of attributes.
#'
#' @param h is the object h represented by a 0-1 vector of attributes.
#'
#' @param i is the object i represented by a 0-1 vector of attributes.
#' @param context is the underlying context.
#' @param attribute_weights is the weight vector for the attributes.
#' @return a ternary fuzzy relation given by an array of dimension
#' n x n x n where n is the number of rows of the context. Higher values
#' of an entry correspond to a larger degree of betweenness.
#' @export
compute_stylized_betweenness <- function(g,h,i,context, attribute_weights=colMeans(context)){

  common_attributes <- which(g==1 & i==1)
  if(length(common_attributes)==0){return(1)}
  ans <- 1-max((1-h[common_attributes])*attribute_weights[common_attributes])
  return(ans)

}

get_whole_stylized_betweenness <- function(context, stylized_betweenness_function=compute_stylized_betweenness,...){
  n_rows <- nrow(context)
  result <- array(0,rep(n_rows,3))
  for(k in (1:n_rows)){
    for(l in (1:n_rows)){
      for(m in (1:n_rows)){
        result[k,l,m] <- stylized_betweenness_function(context[k,],context[l,],context[m,],context, ...)
      }
    }
  }
return(result)
}

#' Compute the quotient order of a quasiorder
#'
#' @description 'compute_quotient_order' computes the quotient order of a
#' quasiorder by deleting equivalent 'duplicates' from the given incidence
#' matrix of the quasiorder
#'
#' @param incidence is the incidence matrix of the quasiorder
#'
#' @return the incidence matrix of the quotient order
#'
#' @export
compute_quotient_order <- function(incidence){
  index <- !duplicated(incidence)&!duplicated(t(incidence))
return(I[index,index])}





cut_incidence=function(incidence,width,interval=stats::quantile(unique(as.vector(incidence)),c(0.001,0.995))){
  vc <- compute_width(compute_quotient_order(compute_transitive_hull(incidence >= interval[2])))$width
  if(vc <= width){return(incidence >= interval[2])}
  f <- function(C,incidence){ width_2 <-compute_width(compute_quotient_order(compute_transitive_hull(incidence >=C)))$width;return(width_2-width)}
  ans <- uniroot(f,interval=interval,incidence=incidence)
return(compute_transitive_hull(incidence>=ans$root))}








starshaped_subgroup_discovery  <- function(stylized_betweenness,objective,local_vc_dim,params=list(Outputflag=0)){

  if (dim(stylized_betweenness)[1] != dim(stylized_betweenness)[2] | dim(stylized_betweenness)[1] != dim(stylized_betweenness)[3] | dim(stylized_betweenness)[2] != dim(stylized_betweenness)[3]){print("dimension mismatch")}
  m <- nrow(stylized_betweenness)
  model <- list(modelsense="max",obj=objective,lb=rep(0,m),ub=rep(1,m))
  solutions <- list()
  objvals <- rep(0,m)
  stars <- array(0,c(m,m))

  models=list()
  for(k in (1:m)){  ## quantify over all starcenters


    if (local_vc_dim == Inf) { incidence <- (stylized_betweenness[k,,] >= max(stylized_betweenness[k,,]))*1}
    else{ incidence <- cut_incidence(Z[k,,],local_vc_dim) }


    model <- model_from_qoset(t(incidence))
    model$obj <- objective
    model$lb <- rep(0,m)
    model$ub <-rep(1,m)
    # force centerpoint to be in the set
    model$lb[k] <- 1
    model$modelsense <- "max"
    b <- gurobi::gurobi(model,params=params)
    solutions[[k]] <- b
    objvals[k] <- b$objval
    stars[k,] <- b$x

    models[[k]] =model


  }
  i <- which.max(objvals)



  if (local_vc_dim == Inf) { incidence <- (stylized_betweenness[i,,] >= max(stylized_betweenness[k,,]))*1}
  else{incidence <- cut_incidence(stylized_betweenness[i,,],local_vc_dim)}


  model <- model_from_qoset(t(incidence))
  model$obj <- objective
  model$lb <- rep(0,m)
  model$ub <-rep(1,m)
  # force centerpoint to be in the set
  model$lb[k] <- 1
  model$modelsense <- "max"
  b <- gurobi::gurobi(model,params=params)
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
    ans <- max(ans,gurobi::gurobi(M,param=list(outputflag=0))$objval)
  }

  return(ans)}


starshaped_subgroup_discovery_h0 <- function(models,params=list(outputflag=0)){

  v <- sample(models[[1]]$obj)

  ans <- -Inf
  for(k in (1: length(models))){
    M <- models[[k]]
    M$obj <- v
    ans <- max(ans,gurobi::gurobi(M,params=params)$objval)
  }

  return(ans)}



model_from_qoset <- function(Q){

  ## constructs linear program for the optimization over all upsets of a quasiordered set Q

  QQ <- compute_pseudoreduction(Q)
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


#######
#######
#######

