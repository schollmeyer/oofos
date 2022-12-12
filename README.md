# oofos
The R package oofos (optimization on fanilies of sets) provides functions for the optimization of linear functions on families of sets (e.g., closure systems, rings of sets and local rings of sets).
For this, the language of formal concept analysis is used. For the concrete implementation, a mixed integer formulation (for closure systems) or a linear programming formulation (for (local) rings of sets) is used (cf., https://www.foundstat.statistik.uni-muenchen.de/personen/mitglieder/schollmeyer/research/tr_209.pdf p.19ff for closure systems and rings of sets and https://georgschollmeyer.weebly.com/uploads/1/3/8/4/138402857/starshaped_subgroup_discovery.pdf for local rings of sets).
The actual computation uses the professional software gurobi ( https://www.gurobi.com ) for which free academic licenses ( https://www.gurobi.com/academia/academic-program-and-licenses/ )are available.

# Installation

You can install the development version of oofos from GitHub with:

``` r
# install.packages("devtools")
devtools::install_github("schollmeyer/oofos")
```
## Example

This is a basic example which shows you how to solve a common problem:

``` r
library(oofos)
context <- (c(1,1,0,0,1,0,1,
	      1,0,0,0,0,1,0,
	      0,1,0,0,0,0,0,
	      0,1,0,1,0,0,0,
	      0,1,0,1,0,0,0,
	      0,0,1,1,1,0,0,
	      0,0,0,0,0,1,0,
	      0,0,0,0,0,1,0,
	      0,1,0,1,0,0,0,
	      1,0,1,0,0,1,0), nrow = 10, ncol = 7, byrow = TRUE)
	      
objective <- c(1,2,3,4,5,6,7,8,9,10) - 55
model <- optimize_on_context_extents(context,(1:10),objective)
model$objval
```
