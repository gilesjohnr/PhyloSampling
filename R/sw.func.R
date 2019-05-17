### phylosampling functions ###

###############################################################################

# A function returning the probabilitiy that two cases are linked by direct transmission 
# given that they are grouped together by linkage criteria
# assuming: single linkage, single transmission, perfect sensitivity

Pr.A1 <- function(chi,    # specificity of the linkage criteria
                  rho,    # sampling proportion
                  M,      # number of cases sampled
                  eta,  # perfect sensitivity
                  R       # effective reproductive number
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
                  M,      # number of cases sampled
                  eta,  # perfect sensitivity
                  R       # effective reproductive number
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
                  eta,    # sensitivity of the linkage criteria
                  R       # effective reproductive number
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
                  eta,    # sensitivity of the linkage criteria
                  R       # effective reproductive number
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

###############################################################################

# A function returning the probabilitiy that two cases are linked by direct transmission 
# given that they are grouped together by linkage criteria
# assuming: single linkage, number of transmission links in the sample per case is possion-distributed

Pr.A3 <- function(chi,    # specificity of the linkage criteria
                  rho,    # sampling proportion
                  M,      # number of cases sampled
                  eta,    # sensitivity of the linkage criteria
                  R       # effective reproductive number
){
  conx1 <- chi >= 0 & chi <= 1
  conx2 <- rho >=0 & rho <=1
  conx3 <- M > 0
  conx4 <- eta >=0 & eta <=1
  conx5 <- R >=0
  
  if (all(conx1, conx2, conx3, conx4, conx5)){
    pr <- (1 - exp(-rho * eta * (R+1))) / 
      (1 - ((chi^(M-2)) * exp(rho * (R+1) * ((1-eta)/chi - 1))))
  } else {
    pr <- NA
  }
  return(pr)
}

# A function returning the expected number of pairs observed in the sample 
# assuming: single linkage, number of transmission links in the sample per case is possion-distributed

EO.A3 <- function(chi,    # specificity of the linkage criteria
                  rho,    # sampling proportion
                  M,      # number of cases sampled
                  eta,    # sensitivity of the linkage criteria
                  R       # effective reproductive number
){
  conx1 <- chi >= 0 & chi <= 1
  conx2 <- rho >=0 & rho <=1
  conx3 <- M > 0
  conx4 <- eta >=0 & eta <=1
  conx5 <- R >=0
  
  if (all(conx1, conx2, conx3, conx4, conx5)){
    pair <- (M * rho * (R+1) * eta * (1 - ((chi^(M-2))) * exp(rho * (R+1) * ((1-eta)/chi -1)))) / 
      (2 * (1 - exp(-rho * (R+1) * eta)))
  } else {
    pair <- NA
  }
  return(pair)
}