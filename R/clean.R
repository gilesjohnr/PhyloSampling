
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




