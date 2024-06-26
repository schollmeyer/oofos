% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/additions_order_theory.R
\name{compute_order_dimension}
\alias{compute_order_dimension}
\title{Compute the order dimension of a partially ordered set}
\usage{
compute_order_dimension(poset, poset_context = NULL)
}
\arguments{
\item{poset}{square incidence matrix with 0-1 entries that describes the poset}

\item{poset_context}{For the computation of the order dimension, the context with all partial orders
as intents is needed. More precisely, one needs a formal context where every object is a linear order L and
every attribute is a pair (a,b) and L I (a,b) iff (a,b) in L. Note that the empty intersection of objects
gives the all relation, and the all relation is not a partial order.
If poset_context is not supplied, then it will be computed beforehand}
}
\value{
the order dimension of the poset.
}
\description{
'compute_order-dimension' computes the order dimension of a partially ordered set,
i.e., the minimal k such that the poset can be represented as the intersection of k linear orderstichain.
}
