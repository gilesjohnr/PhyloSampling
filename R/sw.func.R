### phylosampling functions ###

###############################################################################

# A function returning the probabilitiy that two cases are linked by direct transmission 
# given that they are grouped together by linkage criteria
# assuming: single linkage, single transmission, perfect sensitivity

Pr.A1 <- function(chi,    # specificity of the linkage criteria
                  rho,    # sampling proportion
                  M       # number of cases sampled
){
  conx1 <- chi >= 0 & chi <= 1
  conx2 <- rho >=0 & rho <=1
  conx3 <- M > 0
  
  if (all(conx1, conx2, conx3)){
    pr <- rho / (rho + ((1 - chi^(M-1)) * (1 - rho)))
  } else {
    pr <- NA
  }
  return(pr)
}

# A function returning the expected number of pairs observed in the sample 
# assuming: single linkage, single transmission, perfect sensitivity

EO.A1 <- function(chi,    # specificity of the linkage criteria
                  rho,    # sampling proportion
                  M       # number of cases sampled
){
  conx1 <- chi >= 0 & chi <= 1
  conx2 <- rho >=0 & rho <=1
  conx3 <- M > 0
  
  if (all(conx1, conx2, conx3)){
    pair <- (M / 2) * (rho + ((1 - chi^(M-1)) * (1 - rho)))
  } else {
    pair <- NA
  }
  return(pair)
}

###############################################################################

# A function returning the probabilitiy that two cases are linked by direct transmission 
# given that they are grouped together by linkage criteria
# assuming: single linkage, single transmission

Pr.A2 <- function(chi,    # specificity of the linkage criteria
                  rho,    # sampling proportion
                  M,      # number of cases sampled
                  eta     # sensitivity of the linkage criteria
){
  conx1 <- chi >= 0 & chi <= 1
  conx2 <- rho >=0 & rho <=1
  conx3 <- M > 0
  conx4 <- eta >=0 & eta <=1
  
  if (all(conx1, conx2, conx3, conx4)){
    pr <- (eta * rho) / ((eta * rho) + ((1 - chi^(M-2)) * (1 - eta) * rho) + ((1 - chi^(M-1)) * (1 - rho)))
  } else {
    pr <- NA
  }
  return(pr)
}

# A function returning the expected number of pairs observed in the sample 
# assuming: single linkage, single transmission

EO.A2 <- function(chi,    # specificity of the linkage criteria
                  rho,    # sampling proportion
                  M,      # number of cases sampled
                  eta     # sensitivity of the linkage criteria
){
  conx1 <- chi >= 0 & chi <= 1
  conx2 <- rho >=0 & rho <=1
  conx3 <- M > 0
  conx4 <- eta >=0 & eta <=1
  
  if (all(conx1, conx2, conx3, conx4)){
    pair <- (M / 2) * ((eta * rho) + (rho * (1 - eta) * (1 - chi^(M-2))) + ((1 - rho) * (1 - chi^(M-1))))
  } else {
    pair <- NA
  }
  return(pair)
}