iitoi <- function(ii,n){
  i <- floor((ii-1)/n) + 1
  return(i)
}

iitoj <- function(ii,n){
  j <- (ii-1)%%n + 1
  return(j)
}
