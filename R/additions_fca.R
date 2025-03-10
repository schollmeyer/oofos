reduce_columns <- function(context){

  ## removes reducible attributes from a context
  # TODO description etc.
  result <- NULL
  #m=dim(X)[1]
  #n=dim(X)[2]
  for(k in seq_len(ncol(context))){
    pre_intent <- rep(0,ncol(context))
    pre_intent[k] <- 1
    intent <- operator_closure_attr_input(pre_intent,context=context)
    intent[k] <- 0
    extent <- compute_phi(intent,context)
    if(any(extent !=context[,k])){ result <- c(result,k)}
  }
  return(context[,result])}


get_extreme_attributes <- function(intent,context){
  ## extreme points vs basis points TODO
  extent <- compute_phi(intent,context)
  result <- NULL
  for(k in which(intent==1)){
    intent_new <- intent
    intent_new[k] <- 0
    if(any(compute_phi(intent_new,context)!=extent)){
      result <- c(result,k)
    }


  }
  return(result)
}


compute_example_contexts <- function(){
  planets <- fcaR::planets
  vegas <- fcaR::vegas
  mat <- upper.tri(diag(rep(1,20)))
  diag(mat) <- 1
  mat[3,c(4,5,10)] <- 0
  mat[4,c(5,8)] <- 0
  mat[5,6] <- 0
  mat[6,c((7:11),13) ] <- 0
  mat[7,c((8:10), 14)] <- 0
  mat[8, c((9:11),15)] <- 0
  mat[9, c(10,12)] <- 0
  mat[10,c(13,15)] <- 0
  mat[11,c((12:15),18 )]<- 0
  mat[12,(13:15)]<- 0
  mat[13,c((14:15),17 )]<- 0
  mat[14,15]<- 0
  mat[15,16]<- 0
  mat[16, c(17,18)] <- 0
  mat[17,18] <- 0

  colnames(mat) <- c("Pause","Einklang","Sekunde","Terz","Quart",
                     "quintseptfreier Nonakkord","quintfreier Septakkord",
                     "terzseptfreier Nonakkord","terzfreier Septakkord",
                     "Dreiklang","kompakter Septakkord","septfreier Nonakkord",
                     "quintnonfreier Undezakkord", "terzfreier Nonakkord",
                     "quintfreier Nonakkord","nonfreier Undezakkord",
                     "kompakter Nonakkord","septfreier Undezakkord",
                     "kompakter Undezakkord","kompakter Tredezakkord")
  rownames(mat) <- colnames(mat)

  # Literatur: Rudolf Wille
  #  Musiktheorie und Mathematik 1985. Springer
return(list(planets=planets,vegas=vegas,chords=mat))

}


################################################################################
# Formal Concept Analysis and Implications (general)
################################################################################
################################################################################
# FCA - Calculating of all formal concepts based on a formal context given by a
# cross tabel / incidence
################################################################################
compute_phi <- function(subset_attributes, context) {

  # This function is based on the
  # the master thesis '***''

  # computes for a subset of attributes the minimal extent based on the given context

  # Input: subset_attributes (array): set of attributes
  #         context (matrix): formal context which is used to compute the extent

  # Output: subset (array): the smallest extent (set of objects) in the FCA
  #                         based on subset_attributes and the formalc context

  # Determines and subsets the attributes which are selected
  index_attribute <- which(subset_attributes == 1)
  selected_attributes <- as.matrix(context[, index_attribute])
  dim(selected_attributes) <- c(dim(context)[1], length(index_attribute))

  # Counting for each object how many selected attributes hold and choosing the
  # one where all attributes are true
  count_objects_attribute_hold <- rowSums(selected_attributes)
  index_obejct <- which(count_objects_attribute_hold == length(index_attribute))

  # returning a list which represents which objects correspond to the considered
  # attribute set
  extend <- rep(0, dim(context)[1])
  extend[index_obejct] <- 1

  return(extend)
}


compute_psi <- function(subset_objects, context) {

  # This function (with slight modifications) is taken from the R package
  # ddandrda version 0.0.0.9000 ( https://github.com/hannahblo/ddandrda )

  # computes for a subset of objects the minimal intent based on the given context

  # Input: subset_objects (array): set of objects
  #         context (matrix): formal context which is used to compute the intent

  # Output: subset (array): to smallest intent (set of attributes) in the FCA
  #                         based on subset_objects and context

  # Determines and sub-setting the objects which are selected
  index_object <- which(subset_objects == 1)
  selected_objects <- as.matrix(context[index_object, ])
  dim(selected_objects) <- c(length(index_object), dim(context)[2])

  # Counting for each attribute how many selected objects are related and chose
  # the ones where all objects are related
  count_attributes_object_related <- colSums(selected_objects)
  index_attribute <- which(count_attributes_object_related == length(index_object))

  # returning an array which represents the attributes which correspond to the
  # considered object set
  intent <- rep(0, dim(context)[2])
  intent[index_attribute] <- 1
  return(intent)
}


operator_closure_attr_input <- function(subset_attribute, context) {

  # This function (with slight modifications) is taken from the R package
  # ddandrda version 0.0.0.9000 ( https://github.com/hannahblo/ddandrda )


  # Defines the closure operator for computing all intents (attribute)

  # Input: subset_attribute (array): set of attributes
  #         context (matrix): formal context which is used to compute the intent

  # Output: subset (array): to smallest closure in the FCA based on
  #                         subset_attribute and context

  compute_psi(compute_phi(subset_attribute, context), context)
}


operator_closure_obj_input <- function(subset_object, context) {

  # This function (with slight modifications) is taken from the R package
  # ddandrda version 0.0.0.9000 ( https://github.com/hannahblo/ddandrda )


  # Defines the closure operator for computing all extends (objects)

  # Input: subset_object (array): set of objects
  #         context (matrix): formal context which is used to compute the extent

  # Output: subset (array): to smallest closure in the FCA based on
  #                         subset_object and context
  compute_phi(compute_psi(subset_object, context), context)
}


# Auxiliary functions of compute_all_closure, for algorithm-step 2: next closure
adds_element <- function(old_subset, element) {

  # This function (with slight modifications) is taken from the R package
  # ddandrda version 0.0.0.9000 ( https://github.com/hannahblo/ddandrda )

  # Adds a further element to old_subset and deletes all larger elements
  # based on: Granter (2013), Diskrete Mathematik: Geordnete Mengen, Springer Spektur, p.85

  # input: old_subset (array with 0,1 elements): subset to which the element should be added
  #                                             (1 represents element in subset)
  #         element (integer): element (position) which is added

  # output: subset (array with 0,1 elements): subset with added element
  #                                          (1 represents element in subset)

  # if the element is the first, the subset only consists of this element
  if (element == 1) {
    subset <- rep(0, length(old_subset))
    subset[element] <- 1
  } else {
    index_lower_element_index <- rep(0, length(old_subset))
    index_lower_element_index[(1:(element - 1))] <- 1
    # pmin: A and temp are compared by element by element and the minimum is selected
    subset <- pmin(old_subset, index_lower_element_index)
    subset[element] <- 1
  }
  return(subset)
}


# Auxiliary function of compute_all_closure defining order structure given by 'lektisch' order
compare_closures_lower_i <- function(old_closure, new_closure, element) {

  # This function (with slight modifications) is taken from the R package
  # ddandrda version 0.0.0.9000 ( https://github.com/hannahblo/ddandrda )

    # Tests if the old_closure is smaller than  the new_closure within the meaning of
  # 'lektisch' order
  # based on: Granter (2013), Diskrete Mathematik: Geordnete Mengen, Springer Spekturm, p.26 + 84

  # Input: old_closure (array with 0,1 elements): closure (subset, 1= element within closure)
  #         new_closure (array with 0,1 elements): closure (subset, 1= element within closure)
  #         element (integer): element which is used for comparing

  # Output (logical): returns true if old_closure < new_closure
  if (element == 1) {
    return(new_closure[element] == 1 & old_closure[element] == 0)
  } else {
    temp <- rep(0, length(old_closure))
    temp[(1:(element - 1))] <- 1
    return(new_closure[element] == 1 & old_closure[element] == 0 & all(pmin(old_closure, temp) == pmin(new_closure, temp)))
  }
}


# main functions
compute_all_closure <- function(closure_operator, context,
                                number_attributes = NA,
                                already_computed_closures = 1000) {

  # This function (with slight modifications) is taken from the R package
  # ddandrda version 0.0.0.9000 ( https://github.com/hannahblo/ddandrda )

  # Calculation of all sets of the complete lattice.
  # based on: Granter (2013), Diskrete Mathematik: Geordnete Mengen, Springer Spekturm, p.68

  # Input: closure_operator (func): set-operator which computes the smallest closure
  #               based on a subset
  #         context (matrix): formal context which precises the closure_operator
  #         number_attributes (NA or integer): determines the number of attributes
  #         already_computed_closures (int): states the frequency how often the
  #               information 'how many closures are already computed' is printed

  # Output: (array, elements in 0,1): each row states one computed closure
  #                                   (1 = element in closure)



  if (is.na(number_attributes)) {
    number_attributes <- dim(context)[2]
  }

  # Calculating the first lattice set based on the empty set and the used context
  # Ganter, p68 Algorithm: First Closure
  old_closure <- closure_operator(rep(0, number_attributes), context)
  all_closure <- list()
  all_closure[[1]] <- old_closure

  # In this part all further lattice sets are computed
  t <- 2
  not_all_closures_computed <- TRUE
  while (not_all_closures_computed) {
    attributs_selected <- which(old_closure == 1)

    # Determining all the attributes which could be added, hence which are not
    # selected yet
    if (length(attributs_selected) == 0) {
      index <- (1:number_attributes)
    } else {
      index <- (1:number_attributes)[-attributs_selected]
    }

    # Ganter, p.86 Algorithm: Next Closure
    # Going from the larges to the lowest not yet added attribute until the new
    # computed closure is larger (in sense of the 'lektisch' order, see Ganter p. 26)
    for (element in sort(index, decreasing = TRUE)) {
      # Adding the new element with 'adds_element()' and computing the closure
      new_closure <- closure_operator(adds_element(old_closure, element), context)
      # Test if the new closer is larger then the older closure. If yes, go on.
      if (compare_closures_lower_i(old_closure, new_closure, element)) {
        break # break of the for-loop (not while)
      }
    }

    # Saving the new closure and now it takes the place of the old closure.
    old_closure <- all_closure[[t]] <- new_closure

    # Testing if all closures are computed, the last one has all attributes selected
    if (all(new_closure == 1)) {
      not_all_closures_computed <- FALSE
    }
    # Test if print-information on how many closures are already computed
    if (t %% already_computed_closures == 0) {
      cat(t, "many closures were computed.\n")
    }
    # assignment to the new saving space
    t <- t + 1
  }
  # Convert from list to array and return the object
  return(t(array(unlist(all_closure), dim = c(number_attributes, t - 1))))
}


compute_concept_lattice <- function(context, compute_extents = TRUE) {

  # This function (with slight modifications) is taken from the R package
  # ddandrda version 0.0.0.9000 ( https://github.com/hannahblo/ddandrda )


  # computes the formal concept lattice.
  # Therefore, all formal concept which are defined by the formal context are
  # computed.

  # Input: context (matrix): represents the formal context (rows: objects, columns: attributes)
  #         compute.extends (logical): If it is sufficient to only compute the intent

  # Output: (list) intents(array (0,1 values)): each row represents one intent, (1 = attribute is contained)
  #                 extent(array with 0,1 values): each row represents one extent, (1 = attribute is contained)
  #                 concepts(list): corresponding intent and extend are saved together
  #                                   saving not by index, but directly by their names


  result <- list()

  # Calculating all intents using the closure operator property
  result$intents <- compute_all_closure(operator_closure_attr_input, context)

  number_closure <- dim(result$intents)[1]
  number_objects <- dim(context)[1]
  result$concepts <- rep("", number_closure)

  if (compute_extents) {
    result$extents <- matrix(FALSE, ncol = number_objects, nrow = number_closure)
    for (k in (1:number_closure)) {
      # compute the extends based on the intents
      result$extents[k, ] <- compute_phi(result$intents[k, ], context)
      result$concepts[k] <- paste("{",
        paste((rownames(context))[which(result$extents[k, ] == 1)], collapse = ","),
        "}   {",
        paste((colnames(context))[which(result$intents[k, ] == 1)], collapse = ","),
        "}",
        collapse = ""
      )
    }
  } else {
    for (k in (1:number_closure)) {
      result$concepts[k] <- paste("{",
        paste((colnames(context))[which(result$intents[k, ] == 1)], collapse = ","),
        "}",
        collapse = ""
      )
    }
  }

  return(result)
}


################################################################################
# Order Theory - Calculation of a subconceptrelation based on the closure system
# (here: extends)
################################################################################

compute_incidence <- function(extent_list) {

  # This function (with slight modifications) is taken from the R package
  # ddandrda version 0.0.0.9000 ( https://github.com/hannahblo/ddandrda )

  # generates incidence matrix of a given data table (here it's a closure set)
  # Needed to plot "Begriffsverband"

  # Input: extend_list (matrix): rows represent each object, columns each set
  #                              entry 1: object is in set / entry 0: is not in set

  # Output: subconceptrelations (quadratic matrix, logical): the number of rows
  #         (resp. columns) is the number of sets
  #         entry(i,j) TRUE if i is within j (definition: empty set is within everything)

  number_extends <- dim(extent_list)[1]

  # Defining memory space
  ans <- matrix(FALSE, ncol = number_extends, nrow = number_extends)

  for (k in (1:number_extends)) {
    for (l in (1:number_extends)) {
      # If every element in set k is also in set l, we switch this entry to TRUE
      ans[k, l] <- all(extent_list[k, ] <= extent_list[l, ])
    }
  }
  return(ans)
}


compute_random_context <- function(nrow = 20,
                                   ncol = 10,
                                   prob = 0.5, seed = 1234567) {
  if (!is.null(seed)) {
    withr::with_seed(
      seed = seed,
      result <- matrix(stats::runif(nrow * ncol) <= prob,
        nrow = nrow,
        ncol = ncol
      ) * 1
    )
  } else {
    result <- matrix(stats::runif(nrow * ncol) <= prob,
      nrow = nrow,
      ncol = ncol
    ) * 1
  }
  return(result)
}

###

# # This function (with slight modifications) is taken from the R package
# ddandrda version 0.0.0.9000 ( https://github.com/hannahblo/ddandrda )

# TODO : harmonize with current version of ddandrda

#' Plot a binary relation
#' @description 'plot_relation' plots a binary relation(e.g., a partial order)
#' by drawing the Hasse diagramm of its Dedekind–MacNeille completion. Note
#' that the plotted relation needs not to be a partial order, it can be an
#' arbitrary relation.
#'
#' @param incidence The incidence relation of the partial order.
#'
#' @return The drawing of the partial order in the form of the Hasse diagramm of
#' the Dedekind-MacNeille completion of the given partial order.
#' @examples
#'
#' # plot an interval order
#' upper <- (1:10)
#' lower <- upper - runif(10)
#' interval_order <- array(0, c(10, 10))
#' for (k in (1:10)) {
#'   for (l in (1:10)) {
#'     interval_order[k, l] <- (upper[k] <= lower[l])
#'   }
#' }
#' oofos:::plot_relation(interval_order)
#'
#' # plot a quasiorder
#' q_order <- array(1 * (runif(100) >= 0.9), c(10, 10))
#' diag(q_order) <- 1
#' q_order[1, 2] <- q_order[2, 1] <- 1
#' q_order <- oofos:::compute_transitive_hull(q_order)
#' oofos:::plot_relation(q_order)
#'
#' # plot the Chevron
#' chevron <- rbind(
#'   c(1, 1, 0, 0, 0, 0), c(0, 1, 0, 0, 0, 0), c(1, 1, 1, 0, 1, 1),
#'   c(0, 1, 0, 1, 0, 1), c(0, 0, 0, 0, 1, 1), c(0, 0, 0, 0, 0, 1)
#' )
#' oofos:::plot_relation(chevron)
#'
plot_relation <- function(incidence) {

  # This function (with slight modifications) is taken from the R package
  # ddandrda version 0.0.0.9000 ( https://github.com/hannahblo/ddandrda )


  # define the incidence as a formal context for the fcaR package
  fc <- fcaR::FormalContext$new(incidence[, seq_len(nrow(incidence))])
  # compute all concepts
  fc$find_concepts()
  # plot the Hasse graph of the formal concept lattice, i.e., the diagram of the
  # Dedekind–MacNeille completion of the relation incidence (given as a 0-1
  # matrix)
  fc$concepts$plot()
}
