extent_opt_c <- function(context, gen_index, v, binary_variables = "afap") {


  m <- dim(context)[1]
  n <- dim(context)[2]
  mask <- rep(0, m)
  mask[gen_index] <- 1
  N <- 5 * (m + n)
  NN <- m * n - sum(context)
  I <- rep(as.integer(0), 5 * NN)
  J <- I
  V <- I

  lb <- rep(0, m + n)
  ub <- rep(1, m + n)

  rhs <- rep(0, N)
  sense <- rep("", N)
  t <- 1
  tt <- 1
  for (k in (1:m)) {
    i <- which(context[k, ] == 0)

    if (length(i) >= 1) {
      L <- length(i)
      I[tt] <- t
      J[tt] <- k ## A[t,k]=L;
      V[tt] <- L
      tt <- tt + 1
      index <- (tt:(tt + L - 1))
      I[index] <- t
      J[index] <- i + m ### A[t,i+m]=1;
      V[index] <- 1

      tt <- tt + L
      rhs[t] <- length(i)
      sense[t] <- "<="
      t <- t + 1
    } else {
      lb[k] <- 1
    }
  }


  for (k in (1:n)) {
    i <- which(X[, k] == 0)

    if (length(i) >= 1) {
      L <- length(i)
      I[tt] <- t
      J[tt] <- k + m # A[t,k+m]=L
      V[tt] <- L
      tt <- tt + 1
      index <- (tt:(tt + L - 1))
      I[index] <- t
      J[index] <- i #  A[t,i]=1;
      V[index] <- 1
      tt <- tt + L
      rhs[t] <- length(i)
      sense[t] <- "<="
      t <- t + 1
    } else {
      lb[k + m] <- 1
    }
  }



  for (k in (1:m)) {
    i <- which(X[k, ] == 0)
    L <- length(i)
    if (length(i) >= 1) {
      L <- length(i)
      I[tt] <- t
      J[tt] <- k # A[t,k]=1;
      V[tt] <- 1

      tt <- tt + 1

      index <- (tt:(tt + L - 1))
      I[index] <- t
      J[index] <- i + m # A[t,i+m]=1;
      V[index] <- 1
      tt <- tt + L

      rhs[t] <- 1
      sense[t] <- ">="
      t <- t + 1
    }
  }




  for (k in (1:n)) {
    j <- which(mask == 1 & context[, k] == 0)


    if (length(j) >= 1) {
      L <- length(j)

      I[tt] <- t
      J[tt] <- k + m # A[t,k+m]=1;
      V[tt] <- 1
      tt <- tt + 1
      index <- (tt:(tt + L - 1))
      I[index] <- t
      J[index] <- j # A[t,j]=1;
      V[index] <- 1
      tt <- tt + L
      sense[t] <- ">="
      rhs[t] <- 1
      t <- t + 1
    } else {
      lb[k + m] <- 1
    }
  }

  tt <- tt - 1
  t <- t - 1

  jj <- which(I != 0 & J != 0)

  ###  setze je nach Methode gewisse Variablen als binaer
  vtypes <- rep("C", m + n)
  if (binary_variables == "afap") {
    if (length(gen_index) <= min(m, n)) {
      vtypes[gen_index] <- "B"
    }
    if (length(gen_index) > min(m, n) & m <= n) {
      vtypes[(1:m)] <- "B"
    }
    if (length(gen_index) > min(m, n) & n <= m) {
      vtypes[-(1:m)] <- "B"
    }
  }

  if (binary_variables == "sd") {
    if (m <= n) {
      vtypes[(1:m)] <- "B"
    } else {
      vtypes[-(1:m)] <- "B"
    }
  }

  if (binary_variables == "allgen") {
    vtypes[gen_index] <- "B"
    vtypes[-(1:m)] <- "B"
  }

  if (binary_variables == "all") {
    vtypes <- rep("B", m + n)
  }

  if (!(binary_variables %in% c("afap", "sd", "allgen", "all"))) {
    print("invalid argument for binary_variabes")
  }
  return(list(A = simple_triplet_matrix(I[jj], J[jj], V[jj], nrow = t,
                                        ncol = m + n), rhs = rhs[(1:t)],
              sense = sense[(1:t)], modelsense = "max", lb = lb, ub = ub,
              obj = c(v, rep(0, n)), ext_obj = v, intent_obj = rep(0, n),
              m = m, n = n, vtypes = vtypes, n_constr = t, context = context))
}






###
###
###
#  subgroup discovery

quality <- function(sdtask, result, NAMES = colnames(sdtask$context)) { ## berechnet Piatetsky-Shapiro-Qualitätsfunktion für bereits geloestes Model (Variable result). Variable sdtask ist erzeugtes Modell aus Funktion subgroup.discovery.fca.milp
  m <- sdtask$m

  idx <- which(result$x[(1:m)] > 0.5)
  jdx <- which(result$x[-(1:m)] > 0.5)
  n0 <- length(which(sdtask$obj[(1:m)] > 0))
  n <- length(idx)
  p <- length(which(sdtask$obj[(1:m)] > 0 & result$x[(1:m)] > 0.5)) / n
  p0 <- length(which(sdtask$obj > 0)) / m
  rho <- sqrt(n / m) * (p - p0) / sqrt((1 - n / m) * p0 * (1 - p0))
  return(list(n = n, n0 = n0, p = p, p0 = p0, ps = n * (p - p0), rho = rho, obj = result$objval, argmax = NAMES[jdx]))
}




subgroup_discovery_fca_milp <- function(dat, target, target.class, nrep, heuristic, remove.full.columns = TRUE, clarify.cols = FALSE, reduce.cols = FALSE, weighted = FALSE, small) { # Fuehrt Subgroup Discovery durch (erzeugt nur MILP, das dann noch mit z.B. gurobi(...) gelöst werden muss
  XX <- dat
  XX[[target]] <- NULL
  print(objects(XX))
  X <- conceptual.scaling(XX)
  print("E")



  print(c("dim dat: ", dim(dat)), quote = FALSE)
  print(c("dim X (conceptual scaling dat):     ", dim(X)), quote = FALSE)

  if (remove.full.columns) {
    X <- remove.full.cols(X)
    print(c("dim X without full columns:         ", dim(X)), quote = FALSE)
  }
  if (clarify.cols) {
    X <- col.clarify(X)
    print(c("dim X without doubled columns:      ", dim(X)), quote = FALSE)
  }
  if (reduce.cols) {
    X <- col.reduce(X)
    print(c("dim final X (clarified and reduced):", dim(X)), quote = FALSE)
  }


  m <- dim(X)[1]
  n <- dim(X)[2]

  M <- dim(XX)[2]

  # A=NULL
  # A=array(0,c(2*n,m+n))
  # t=1
  # T=1
  # for(k in (1:M)){
  # K=length(unique(XX[,k]))
  # #print(K)
  # temp=rep(0,c(m+n))
  # if(class(XX[,k])[1]=="factor"){

  # temp[((t+m):(t+m+K-1))]=1
  # A[t,]=temp
  # T=T+1
  # #A=rbind(A,temp)

  # t=t+K
  # }

  # if(class(XX[,k])[1]=="ordered" | class(XX[,k])[1]=="numeric" | class(XX[,k])[1]=="integer"){

  # for(l in (1:(K-1))){
  # temp[t+m+l-1]=1
  # temp[t+m+l]=-1
  # A[T,]=temp
  # T=T+1
  # #A=rbind(A,temp)

  # temp=rep(0,c(m+n))
  # temp[t+m+l-1+K]=-1
  # temp[t+m+l+K]=1
  # A[T,]=temp
  # T=T+1
  # #=rbind(A,temp)
  # }
  # t=t+2*K
  # }



  # #t=t+K
  # }
  # T=T-1
  # A=A[(1:T),]
  # A=NULL
  v <- (dat[[target]] == target.class) - (dat[[target]] != target.class) * mean(dat[[target]] == target.class) / mean(dat[[target]] != target.class)
  if (weighted) {
    W <- weighted.repr(X, v)
    v <- W$yw

    X <- W$Xw
  } else {
    W <- list(count = rep(1, dim(X)[1]))
  }

  v <<- v
  m <- dim(X)[1]
  n <- dim(X)[2]
  print("building model...")
  if (small) {
    ans <- extent.opt(X, which(v > 0), v)
  } else {
    ans <- extent.opt(X, (1:m), v)
  }
  print("...done")
  # M=dim(A)[1]
  # ans$A=rbind(ans$A,A);ans$rhs=c(ans$rhs,rep(1,M));ans$sense=c(ans$sense,rep("<=",M))
  print(system.time(temp <- heuristic(X, v, nrep = nrep)))
  CUT <- ceiling(temp$objval)
  print(CUT)

  T <- rep(0, n + m)
  T[(which(v > 0))] <- -W$count[which(v > 0)]
  # ans$A=rbind(ans$A,matrix(T,nrow=1))
  # ans$rhs=c(ans$rhs,-CUT)
  # ans$sense=c(ans$sense,"<=")
  ans$start <- c(temp$solution, PSI(temp$solution, context = X))
  ans$context <- X
  ans$NAMES <- conceptual.scaling.dim(XX)$NAMES
  # TEMP=heuristic.implications(X,v,NREP)

  # ans$A=rbind(ans$A,TEMP$A);ans$rhs=c(ans$rhs,TEMP$rhs);ans$sense=c(ans$sense,TEMP$sense)

  # ans$vtypes[which(ans$start>=0.5)]="B"
  # ans$vtypes=rep("B",length(ans$vtypes))
  return(ans)
}
