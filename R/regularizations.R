#Regularisierungsschema regularization_sceme_a: Context -> DIST -> NMDS -> Context -> best_fitting_attributes

regularization_scheme_a <- function(context,distance_function=dist,d,n_attributes,threshold){
  distance <- (as.matrix(distance_function(context)))
  nmds_result <- vegan::metaMDS(distance,k=d,distance="euclidean",trace=trace,try=20,trymax=1000,maxit=10000,autotransform = FALSE)
  new_distance <- as.matrix(distance_function(nmds_result$points))
  new_context <- get_context_from_distance(new_distance,threshold=threshold)$context
  attribute_fits <- MatrixGenerics::rowMins(Rfast::dista(t(new_context),t(context)))
  n_attributes <- min(n_attributes, ncol(new_context))
  return(list(regularized_context=new_context[,(order(attribute_fits))[seq_len(n_attributes)]]))
}

regularization_scheme_b <- function(context,dist_threshold, stress_threshold,distance_function=dist,dimensions,n_attributes,threshold){
  distance <- (as.matrix(distance_function(context)))
  nmds_results <- list()
  residuals <- list()
  new_distances <- list()
  new_contexts <- list()
  attribute_fits <- list()
  new_context <- NULL
  t <- 1
  for(d in (dimensions)){
    nmds_results[[t]] <- vegan::metaMDS(distance,k=d,distance="euclidean",trace=trace,try=20,trymax=1000,maxit=10000)
    residuals[[t]] <- vegan::goodness(nmds_results[[t]],distance)
    new_distances[[t]] <- as.matrix(distance_function((nmds_results[[t]])$points))
    temp <- get_context_from_distance(new_distances[[t]],threshold=threshold, object_stresses=residuals[[t]])
    new_contexts[[t]] <- temp$context
    attribute_fits[[t]] <- MatrixGenerics::rowMins(Rfast::dista(t(new_contexts[[t]]),t(context)))
    new_context <- cbind(new_context,(new_contexts[[t]])[,which(attribute_fits[[t]] <= dist_threshold & temp$halfspace_stresses*dimensions[t] <= stress_threshold)])
    t <- t+1
  }

  #new_context <- get_context_from_distance(new_distance,threshold=threshold)$context
  #attribute_fits <- MatrixGenerics::rowMins(Rfast::dista(t(new_context),t(context)))
  #n_attributes <- min(n_attributes, ncol(new_context))
  return(list(regularized_context=new_context))
}




get_nmds_results <- function(context,dimensions,distance_function=dist,trace=2){
  distance <- (as.matrix(distance_function(context)))
  nmds_results <- list()
  residuals <- list()
  new_distances <- list()
  diameters <- list()

  t <- 1
  for(d in (dimensions)){
    nmds_results[[t]] <- vegan::metaMDS(distance,k=d,distance="euclidean",trace=trace,try=20,trymax=1000,maxit=10000)
    residuals[[t]] <- vegan::goodness(nmds_results[[t]],distance)
    new_distances[[t]] <- as.matrix(distance_function((nmds_results[[t]])$points))
    diameters[[t]] <- max(new_distances[[t]])
    #temp <- get_context_from_distance(new_distances[[t]],threshold=threshold, object_stresses=residuals[[t]])
    #new_contexts[[t]] <- temp$context
    #attribute_fits[[t]] <- MatrixGenerics::rowMins(Rfast::dista(t(new_contexts[[t]]),t(context)))
    #new_context <- cbind(new_context,(new_contexts[[t]])[,which(attribute_fits[[t]] <= dist_threshold & temp$halfspace_stresses*dimensions[t] <= stress_threshold)])
    t <- t+1
  }

  #new_context <- get_context_from_distance(new_distance,threshold=threshold)$context
  #attribute_fits <- MatrixGenerics::rowMins(Rfast::dista(t(new_context),t(context)))
  #n_attributes <- min(n_attributes, ncol(new_context))
  return(list(nmds_results=nmds_results, residuals=residuals, new_distances=new_distances, dimensions=dimensions))
}


check_three_point_condition <- function(dist_mat,eps=10^-6,lambda=1){
  m <- ncol(dist_mat)
  counterexamples <- array(0,c(m,m))
  for( k in (1:m)){
    #print(k)
    for(l in (1:m)){
      # old version
      #counterexamples[k,l] <- lambda*sum(dist_mat[k,] > eps+ (pmax(dist_mat[k,l],dist_mat[l,])))
      counterexamples[k,l] <-  sum(dist_mat[,k] > eps+ (pmax(dist_mat[k,l],dist_mat[,l])))

      #counterexamples[k,l] <- counterexamples[k,l]+ (1-lambda)*sum(pmax(dist_mat[k,] -( eps+ (pmax(dist_mat[k,l],dist_mat[l,]))),0))



      #counterexamples[k,l] <- sum(dist_mat[k,l] > eps+ (pmax(dist_mat[k,],dist_mat[,l])))
    }



  }
  #print(sum(counterexamples>0))
  return(list(counterexamples=counterexamples,result=all(counterexamples==0)))

}

get_context_distances <- function(context){
  n_col <- ncol(context)
  result <- array(0,c(n_col,n_col))
  for(k in seq_len(n_col)){
    for(l in seq_len(n_col)){
      result[k,l] <- sum(context[,k] !=context[,l])
    }
  }
  return(result)}

test <- function(model,context_distances,threshold,K_max){
  I <- context_distances > threshold
  diag(I) <- 1
  graph <- igraph::graph_from_adjacency_matrix(I)
  graph <- igraph::as_graphnel(graph)
  coloring <- RBGL::sequential.vertex.coloring(graph)
  K <- coloring[[1]]
  result <- model
  A <- model$A
  B <- array(0,c(nrow(A),K))
  t <- 1
  genconmax <- list()
  for(k in seq_len(K)){
    indexs <- as.vector(which(coloring[[2]] == k - 1))
    genconmax[[t]] <- list(resvar = ncol(A)+ t,vars = indexs + nrow(model$context))
    t <- t+1
  }

  result$A <- cbind(A,B)
  result$A <- rbind(result$A,matrix(c(rep(0,ncol(A)),rep(1,K)),nrow=1))
  result$rhs <- c(result$rhs,K_max)
  result$sense <- c(result$sense,"<=")
  result$lb <- c(result$lb,rep(0,K))
  result$ub <- c(result$ub,rep(1,K))
  result$obj <- c(result$obj,rep(0,K))
  result$vtypes <- c(result$vtypes,rep("B",K))
  result$genconmax <- genconmax
  result$groups <- coloring[[2]]
  result$n_groups=coloring[[1]]
  result$K_max <- K_max
  result$index_K_max_contsraint <- length(result$rhs)
  return(result)
}
get_context_from_distance <- function(dist_mat,threshold,object_stresses=rep(0,ncol(dist_mat)^2),complemented=FALSE,indexs=NULL,sampling_proportion=1,remove_duplicates=TRUE,set_seed=TRUE,seed=1234567,eps=10^-10,eps2=10^-10,lambda=1,counterexamples = check_three_point_condition(dist_mat,eps=eps2,lambda=lambda)$counterexamples){
  n_rows <- nrow(dist_mat)
  n_rows_sample <- ceiling(sampling_proportion*n_rows)
  if(set_seed){set.seed(seed)}
  sampled_indexs <- sample(seq_len(n_rows),size=n_rows_sample)
  if(sampling_proportion==1){sampled_indexs <- seq_len(n_rows)}
  if(!is.null(indexs)){sampled_indexs <-indexs;n_rows_sample <- length(indexs)}
  context <- array(FALSE,c(n_rows,n_rows_sample*(n_rows_sample-1)))
  ce <- rep(0,n_rows_sample*(n_rows_sample-1))
  halfspace_stresses <- ce
  #counterexamples <- check_three_point_condition(dist_mat[sampled_indexs,sampled_indexs],eps=eps2,lambda=lambda)$counterexamples
  t <- 1
  col_names <- NULL
  for(k in seq_len(n_rows_sample)){
    #print(k)
    for(l in seq_len(n_rows_sample)[-k]){
      if(counterexamples[k,l] <= threshold){
        #print(counterexamples[k,l])
        #print(c(k,l))
        context[,t] <- (dist_mat[,sampled_indexs[k]] > dist_mat[,sampled_indexs[l]] +eps)
        ce[t] <- counterexamples[k,l]
        halfspace_stresses[t] <- max(object_stresses[k],object_stresses[l])
        #col_names <- c(col_names,paste(k,">",l,collaps=""))
        #print(col_names)
        t <- t + 1
      }

    }
  }
  context <- context[,(1:(t-1))]
  colnames(context) <- col_names
  if(complemented){context <- (cbind(context,1-context))}
  if(remove_duplicates){
    non_duplicates <- which(!duplicated(t(context)))
    context <- context[,non_duplicates]
    ce <- ce[non_duplicates]
    halfspace_stresses <- halfspace_stresses[non_duplicates]

  }
  return(list(context=context,counterexamples=ce,halfspace_stresses=halfspace_stresses))}
