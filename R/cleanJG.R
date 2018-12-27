
# A function returning the probabilitiy that two cases are linked by direct transmission 
# given that they are grouped together by phylogenetic criteria 

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

# Check
Pr.TL(x=0.99, rho=0.99, N=1001)
Pr.TL(x=0.99, rho=0.9, N=11)

# Explore
rho.vals <- seq(0, 0.99, length.out=50)
N.vals <- seq(1, 1e+06, length.out=50)

m <- matrix(NA, nrow=length(rho.vals), ncol=length(N.vals))

for (i in seq_along(rho.vals)) {
     for(j in seq_along(N.vals)) {
          
          m[i,j] <- Pr.TL(x=0.99, rho=rho.vals[i], N=N.vals[j])
          
     }
}

wireframe(m)

rho.vals <- seq(0, 0.99, length.out=50)
N.vals <- seq(1000, 1e+05, length.out=50)

m <- mesh(rho.vals, N.vals)

Pr.vals <- with(m, Pr.TL(x=0.99, rho=x, N=y))

par(mfrow=c(1,3))
persp3D(x=rho.vals, y=N.vals, z=Pr.vals,
        theta = 210,
        phi = 30,
        colkey=F,
        border='black',
        col=ramp.col (col = c("red", "darkcyan", "cyan", "lightcyan", "greenyellow"), n = 100, alpha = 1),
        ticktype="detailed",
        zlab='Pr( T(A,B) | L(A,B) )',
        xlab='rho',
        ylab='N',
        main='x = 0.99')

Pr.vals <- with(m, Pr.TL(x=0.9, rho=x, N=y))

persp3D(x=rho.vals, y=N.vals, z=Pr.vals,
        theta = 210,
        phi = 30,
        colkey=F,
        border='black',
        col=ramp.col (col = c("red", "darkcyan", "cyan", "lightcyan", "greenyellow"), n = 100, alpha = 1),
        ticktype="detailed",
        zlab='Pr( T(A,B) | L(A,B) )',
        xlab='rho',
        ylab='N',
        main='x = 0.9')

Pr.vals <- with(m, Pr.TL(x=0.8, rho=x, N=y))

persp3D(x=rho.vals, y=N.vals, z=Pr.vals,
        theta = 210,
        phi = 30,
        colkey=F,
        border='black',
        col=ramp.col (col = c("red", "darkcyan", "cyan", "lightcyan", "greenyellow"), n = 100, alpha = 1),
        ticktype="detailed",
        zlab='Pr( T(A,B) | L(A,B) )',
        xlab='rho',
        ylab='N',
        main='x = 0.8')

#################################################
# A function to return the conditional probability of sampling A given that B
# is already in the sample (for small sample sizes)

Pr.AB <- function(N,  # final outbreak size
                  M   # sample size
) {
  
  (M/N) * ((M-1)/(N-1))
  
}

# Defining probability of transmission given linkage using sensitivity instead of specificity
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

# Some checks
Pr.TL2(N=100000, rho=0.99, phi=1)
Pr.TL2(N=100000, rho=0.8, phi=1)
Pr.TL2(N=100000, rho=0.4, phi=1)
Pr.TL2(N=100000, rho=0.2, phi=1)

Pr.TL2(N=1001, rho=0.99, phi=1)
Pr.TL2(N=1001, rho=0.8, phi=1)
Pr.TL2(N=1001, rho=0.4, phi=1)
Pr.TL2(N=1001, rho=0.2, phi=1)

Pr.TL2(N=101, rho=0.99, phi=1)
Pr.TL2(N=101, rho=0.8, phi=1)
Pr.TL2(N=101, rho=0.4, phi=1)
Pr.TL2(N=101, rho=0.2, phi=1)

Pr.TL2(N=11, rho=0.99, phi=1)
Pr.TL2(N=11, rho=0.8, phi=1)
Pr.TL2(N=11, rho=0.4, phi=1)
Pr.TL2(N=11, rho=0.2, phi=1)

Pr.TL2(N=100000, phi=0.99, rho=0.99)
Pr.TL2(N=100000, phi=0.8,  rho=0.99)
Pr.TL2(N=100000, phi=0.4,  rho=0.99)
Pr.TL2(N=100000, phi=0.2,  rho=0.99)

# Explore more 
rho.vals <- seq(0, 0.99, length.out=50)
phi.vals <- seq(0, 0.99, length.out=50)
  #seq(1000, 1e+10, length.out=50)

m <- mesh(rho.vals, phi.vals)
Pr.vals <- with(m, Pr.TL2(N=10, rho=x, phi=y))
par(mfrow=c(1,3))
persp3D(x=rho.vals, y=phi.vals, z=Pr.vals,
        zlim=c(0,1),
        theta = 220,
        phi = 30,
        colkey=F,
        border='black',
        col=ramp.col (col = c("red", "darkcyan", "cyan", "lightcyan", "greenyellow"), n = 100, alpha = 1),
        ticktype="detailed",
        zlab='Pr( T(A,B) | L(A,B) )',
        xlab='rho',
        ylab='phi',
        main='N = 10')

Pr.vals <- with(m, Pr.TL2(N=1000, rho=x, phi=y))

persp3D(x=rho.vals, y=phi.vals, z=Pr.vals,
        zlim=c(0,1),
        theta = 220,
        phi = 30,
        colkey=F,
        border='black',
        col=ramp.col (col = c("red", "darkcyan", "cyan", "lightcyan", "greenyellow"), n = 100, alpha = 1),
        ticktype="detailed",
        zlab='Pr( T(A,B) | L(A,B) )',
        xlab='rho',
        ylab='phi',
        main='N = 1000')

Pr.vals <- with(m, Pr.TL2(N=100000, rho=x, phi=y))

persp3D(x=rho.vals, y=phi.vals, z=Pr.vals,
        zlim=c(0,1),
        theta = 220,
        phi = 30,
        colkey=F,
        border='black',
        col=ramp.col (col = c("red", "darkcyan", "cyan", "lightcyan", "greenyellow"), n = 100, alpha = 1),
        ticktype="detailed",
        zlab='Pr( T(A,B) | L(A,B) )',
        xlab='rho',
        ylab='phi',
        main='N = 100000')

## Contour plots
par(mfrow=c(1,3))

m <- matrix(NA, nrow=length(rho.vals), ncol=length(phi.vals))
for (i in seq_along(rho.vals)) {
  for(j in seq_along(phi.vals)) {
    
    m[i,j] <- Pr.TL2(N=100, rho=rho.vals[i], phi=phi.vals[j])
    
  }
}

contour(x=rho.vals, y=phi.vals, z=m, main='N = 100', xlab='rho', ylab='phi')

m <- matrix(NA, nrow=length(rho.vals), ncol=length(phi.vals))
for (i in seq_along(rho.vals)) {
  for(j in seq_along(phi.vals)) {
    
    m[i,j] <- Pr.TL2(N=10000, rho=rho.vals[i], phi=phi.vals[j])
    
  }
}

contour(x=rho.vals, y=phi.vals, z=m, main='N = 10000', xlab='rho', ylab='phi')

m <- matrix(NA, nrow=length(rho.vals), ncol=length(phi.vals))
for (i in seq_along(rho.vals)) {
  for(j in seq_along(phi.vals)) {
    
    m[i,j] <- Pr.TL2(N=100000, rho=rho.vals[i], phi=phi.vals[j])
    
  }
}

contour(x=rho.vals, y=phi.vals, z=m, main='N = 100000', xlab='rho', ylab='phi')
