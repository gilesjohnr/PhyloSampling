
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

Pr.TL2(N=1000, phi=0.99, rho=0.99)
Pr.TL2(N=1000, phi=0.8,  rho=0.99)
Pr.TL2(N=1000, phi=0.4,  rho=0.99)
Pr.TL2(N=1000, phi=0.2,  rho=0.99)

# Explore more 
rho.vals <- seq(0, 0.99, length.out=100)
phi.vals <- seq(0, 0.99, length.out=100)
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


# Plotly library
library(plotly)


plot_ly(x=rho.vals, y=phi.vals, z=m, type = "surface") %>%
  layout(
    title = "N=1000",
    scene = list(
      xaxis = list(title = "rho"),
      yaxis = list(title = "phi"),
      zaxis = list(title = "Pr")
    ))
  
### ggplot heatmap + contour
g <- expand.grid(rho.vals, phi.vals)
names(g) <- c('rho', 'phi')
g$Pr <- Pr.TL2(N=1000, rho=g$rho, phi=g$phi)

p <- ggplot(g, aes(x=phi, y=rho, z=Pr))
  
p + stat_contour()

p + geom_raster(aes(fill = Pr), hjust=0.5, vjust=0.5, interpolate=FALSE) + 
  scale_fill_gradient(low="blue", high="green", name='Pr(A=A*|L(A,B))') +
  geom_contour(aes(z = Pr), color='black', size=0.5) +
  geom_text_contour(aes(z = Pr), nudge_y=0.0175) +
  theme_bw() +
  ylab(expression(rho)) + 
  xlab(expression(phi^'+')) + 
  theme(axis.text=element_text(size=12),
        axis.title.x=element_text(size=15),
        axis.title.y=element_text(size=15),
        axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())

#########################################################################
# Checking Expected # true pairs

tp <- function(N, rho=1) (rho/2) * N
tp(2, 1)

tp2 <- function(N, rho=1) {
  ( ((N*(N - 1)) /2) ) * ( rho * (((rho*N) - 1) / (N - 1)) )
}
tp2(N=2, rho=1)

tp(10, 1)
tp2(10, 1)

N <- 100
M <- 0:100
k <- 5

rho <- M/N
tau <- k/N

#x <- (tau)# * ((1-tau)

#plot(x*5, type='l', ylim=c(0,1))
1 - (tau * (1 - tau)^(M-k))*5

choose(M, 1) * (tau* (1 - tau)^(M-1))

(tau) * (1 - (rho*tau))
dbinom(1, M, tau)
plot(1 - pbinom(0.5, M, tau), 
     type='l', ylim=c(0,1))

lines((M/N)^2, type='l', ylim=c(0,1))

#################################################################################

f1 <- function(chi, M=100) (1 - chi^(M-1))*(1-chi)
f2 <- function(chi, M=100) ((1 - chi)^(M-1))*(1-chi)
f3 <- function(chi, M=100) ((1 - (1 - chi))^(M-1))*(1-chi)
f4 <- function(chi, M=100) (chi^(M-1))*(1-chi)

plot(f1(seq(0,1,length.out=100)), lwd=2, type='l', xlab='chi', xaxt='n', ylab='Value of term')
lines(f2(seq(0,1,length.out=100)), lwd=2, col='blue')
lines(f3(seq(0,1,length.out=100)), lwd=2, col='red')
lines(f4(seq(0,1,length.out=100)), lwd=2, col='purple')
axis(1,seq(0, 100, 20), seq(0, 1, 0.2))

p1 <- function(chi, M=1000, rho=1) rho/(rho + ((1 - chi^(M-1))*(1-chi))*(1-rho))
p2 <- function(chi, M=1000, rho=1) rho/(rho + (((1 - chi)^(M-1))*(1-chi))*(1-rho))
p3 <- function(chi, M=1000, rho=1) rho/(rho + ( ( ((1 - (1 - chi))^(M-1)) * (1-chi) ) * (1-rho) ))
p4 <- function(chi, M=1000, rho=1) rho/(rho + ((chi^(M-1))*(1-chi))*(1-rho))

plot(p1(seq(0,1,length.out=100)), lwd=2, type='l', ylim=c(0,1), xlab='chi', xaxt='n', ylab='Value of term')
lines(p2(seq(0,1,length.out=100)), lwd=2, col='blue')
lines(p3(seq(0,1,length.out=100)), lwd=2, col='red')
lines(p4(seq(0,1,length.out=100)), lwd=2, col='purple')
axis(1,seq(0, 100, 20), seq(0, 1, 0.2))
abline(h=1, col='green', lty=2, lwd=2)

plot(p1(seq(0,1,length.out=100), rho=0.5), lwd=2, type='l', ylim=c(0,1), xlab='chi', xaxt='n', ylab='Pr(y|z)')
lines(p2(seq(0,1,length.out=100), rho=0.5), lwd=2, col='blue')
lines(p3(seq(0,1,length.out=100), rho=0.5), lwd=2, col='red')
lines(p4(seq(0,1,length.out=100), rho=0.5), lwd=2, col='purple')
axis(1,seq(0, 100, 20), seq(0, 1, 0.2))
abline(h=0.5, col='green', lty=2, lwd=2)

plot(p1(seq(0,1,length.out=100), rho=0.25), lwd=2, type='l', ylim=c(0,1), xlab='chi', xaxt='n', ylab='Pr(y|z)')
lines(p2(seq(0,1,length.out=100), rho=0.25), lwd=2, col='blue')
lines(p3(seq(0,1,length.out=100), rho=0.25), lwd=2, col='red')
lines(p4(seq(0,1,length.out=100), rho=0.25), lwd=2, col='purple')
axis(1,seq(0, 100, 20), seq(0, 1, 0.2))
abline(h=0.25, col='green', lty=2, lwd=2)

p1(1, rho=1)
p2(1, rho=1)
p3(1, rho=1)
p4(1, rho=1)

p1(1, rho=0.5)
p2(1, rho=0.5)
p3(1, rho=0.5)
p4(1, rho=0.5)

#################################################################################
# Probability of y given z using Poisson distributed k

Pr.yz <- function(N,         # final outbreak size
                  rho,       # sampling proportion 
                  eta,       # in sample sensitivity
                  chi,        # specificity of linkage criteria
                  lambda      # Average k per case i
){
  
  cons1 <- N > 0
  cons2 <- rho >= 0 & rho <= 1
  cons3 <- eta >= 0 & eta <= 1
  cons4 <- chi >= 0 & chi <= 1
  M <- N * rho
  
  if (all(cons1, cons2, cons3, cons4)) {
    
    num <- 1 - exp(-(rho*lambda*eta))
    den <- (1 - chi^(N-1)) * (exp((1/chi) * ((1-chi)*rho*lambda - rho*lambda*eta)))
    out <- num/den
  } else {
    out <- NA
  }
  return(out)
}

Pr.yz(N=100, rho=0.05, eta=0.9, chi=0.9, lambda=1)


