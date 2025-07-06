# region.R
# This file defines a region in which an outbreak may occur.
# Each function returns a fixed set of region parameters.
region <- function(id) {
  pop <- 1e3 # population size
  ok_prob <- runif(1, 0.01, 0.99) # outbreak probability
  return(list(id = id, pop = pop, ok_prob = ok_prob))
}
