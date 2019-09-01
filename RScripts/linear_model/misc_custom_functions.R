#---- Conditional probabilities function ----
cond_prob <- function(x){
  probs <- vector(length = length(x))
  for(i in 2:length(probs)){
    probs[i] = x[i]/x[i - 1]
  }
  probs[1] <- 1
  return(probs)
}

#---- Recursive subtraction function ----
sub_recurse <- function(x){
  diffs <- vector(length = length(x))
  for(i in 2:length(diffs)){
    diffs[i] = x[i] - x[i - 1]
  }
  diffs[1] <- x[1]
  return(diffs)
}
