# This is the first generation function that returns the probabilitiy that two cases are linked by direct transmission 
# given that they are grouped together by phylogenetic criteria.

Pr.TL <- function(x=0.99,    # specificity when infector not in sample 
                  rho,       # proportion of potential infectors sampled
                  N          # total population size
){
  conx1 <- x >= 0 & x <= 1
  conx2 <- rho >= 0 & rho <= 1
  conx3 <- N > 0
  
  if (all(conx1, conx2, conx3)) {
    out <- 1 - (((1 - x) * (1 - rho)) / ((1 / (N - 1)) + ((1 - x) * (1 - rho))))
  } else {
    out <- NA
  }
  return(out)
}

# A second generation function defining Pr(A=A*|L(a,B)) using sensitivity instead of specificity
Pr.TL2 <- function(N,         # final outbreak size
                   rho,       # sampling proportion 
                   phi        # in sample sensitivity
){
  
  cons1 <- N > 0
  cons2 <- rho >= 0 & rho <= 1
  cons3 <- phi >= 0 & phi <= 1
  M <- N * rho
  
  if (all(cons1, cons2, cons3)) {
    
    num <- phi * (rho * (1 / (N - 1)))
    out <-  num / (num + ((1 / (M - 1)) * ((1 - rho) * ((N - 2) / (N -1))) * ((N - 2) / (N -1))  ))
    
  } else {
    out <- NA
  }
  return(out)
}