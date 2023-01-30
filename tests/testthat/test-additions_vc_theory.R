test_that("enumerate_ufg_premises works", {
 p5 <- ddandrda::compute_all_partial_orders(5,list=FALSE,complemented=TRUE)
 i <- sample(seq_len(nrow(p5)))
 p5 <- p5[i,]
 n_row_context <-17

 system.time(result_1 <- enumerate_ufg_premises(p5,n_row_context))


 subsets <- gtools::permutations(2,n_row_context,repeats.allowed=TRUE)-1
 ufg_dim <- 5*4/2
 i <- which(rowSums(subsets)<=ufg_dim)
 subsets <- subsets[i,]
 number_ufgs <- 0
 system.time(
 for(k in seq_len(nrow(subsets))){
   if(test_explicitly_ufg_p_order( c(subsets[k,],rep(0,4231-n_row_context)), p5)){
     #print(which(subsets[k,]==1))
     number_ufgs <- number_ufgs + 1
     #print(number_ufgs)
   }
 }
)
expect_equal(number_ufgs,length(result_1))




while(TRUE){
  subset <- rep(0,nrow(p5))
  index <- sample(seq_len(nrow(p5)),size=sample((3:6),size=1))
  subset[index] <- 1
  if(test_explicitly_ufg_p_order(subset,p5)){break}
}

ans <- FALSE
for(l in index){
  subset_new <- subset; subset_new[l] <-0
  if(test_explicitly_ufg_p_order(subset_new,p5)){ans <- TRUE}

}


expect_equal(ans,TRUE)

})



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

test_that("add_ufg_contraints works", {
  example_posets <- compute_example_posets(8)
  context <- NULL
  for(k in seq_len(length(example_posets))){
    temp <- example_posets[[k]]
    if(nrow(temp)==8 & ncol(temp)==8){
      context <- rbind(context, c(as.vector(temp), 1-as.vector(temp)))
    }
  }
    model <- compute_extent_vc_dimension(context)
    result_1 <- gurobi::gurobi(model)$objval
    model_2 <- add_ufg_constraints(model)
    result_2 <- gurobi::gurobi(model_2)$objval
    expect_equal(result_1 >= result_2,TRUE)



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

# TODO: Was soll das

n_test <- 4
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

