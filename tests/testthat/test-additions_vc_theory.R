test_that("compute_extent_vc_dimension works", {

  context <- diag(rep(1,10))
  model <- compute_extent_vc_dimension(context)
  result <- gurobi::gurobi(model)
  expect_equal(result$objval,2)

  context <- 1-diag(rep(1,10))
  model <- compute_extent_vc_dimension(context)
  result <- gurobi::gurobi(model)
  expect_equal(result$objval,10)

  context <- upper.tri(diag(1,20))
  model <- compute_extent_vc_dimension(context)
  result <- gurobi::gurobi(model)
  expect_equal(result$objval,1)

  context <- cbind(upper.tri(diag(1,3)),lower.tri(diag(1,3)))
  model <- compute_extent_vc_dimension(context)
  result <- gurobi::gurobi(model)
  expect_equal(result$objval,2)



})



test_that("check_objset_sufg_candidate works", {

# TODO : hierarchical nominaly scaling method: curate, comment ...
###################################################################
###################################################################
get_hierarchical_scaling_vec <- function(dat){

    # computes hierarchical scaling TODO currently only 4 levels

    m <- length(dat)
    w <- rep(0,m*4)
    t <- 1
    for(l in c(1,10,100,1000)){
      for(k in (1:m)){
        w[t] <- dat[k]%/%l
        t <- t+1
      }
    }
    w <- unique(w)
    context <- NULL
    for(W in w){
      temp <- rep(0,m)
      L <-nchar(W)
      for(k in (1:m)){
        temp[k] <- substr(as.character(dat[k]),1,L)==W
      }
      context <- cbind(context,temp)
    }

    colnames(context) <- w

    return(context)}

############################
############################

levels <- 4
dat <- gtools::permutations(levels,levels,repeats.allowed=TRUE)
n_rows <- nrow(dat)
x <-rep(0,n_rows)
for(k in seq_len(n_rows)){
  x[k] <- as.numeric(paste(dat[k,],collapse=""))

}

x <- sample(x)

context <- get_hierarchical_scaling_vec(x)

n_test <- 7
result_1 <- rep(TRUE,n_test*(n_test-1))
result_2 <- rep(TRUE,n_test*(n_test-1))

t <- 1

for(k in (1:n_test)){
  for(l in (1:n_test)[-k]){
    subset <- rep(0,n_rows)
    subset[c(k,l)] <- 1
    result_1[t] <- check_objset_sufg_candidate(subset,context,2)$result
    result_2[t] <- TRUE#(sum(pmin(context[k,],context[l,]))<=2)
    if(result_2[t]==FALSE){print(c(k,l));break}
    t <- t+1

}
}
expect_equal(all(result_1 == result_2),TRUE)
})

