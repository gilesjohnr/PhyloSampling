### phylosampling plots ###

library(ggplot2)
library(reshape2)
library(gridExtra)
library(lattice)
library(colorspace)

###############################################################################

# A1
# Make four plots showing relationship between M, rho, false discovery rate, and expected number of links
# assumptions: single linkage, single transmission, perfect sensitivity
# assumptions: fixed values of specificity

chi <- 0.99 # specificity of the linkage criteria
eta <- 1 # perfect sensitivity

p1 <- plt.eq(chi = chi,
             eta = eta,
             rho = seq(0, 0.99, length.out = 100),
             M = c(10,50,100,200),
             x = "rho",
             eq = EO.A1,
             lbls = c("rho","expected # of links"))

p2 <- plt.eq(chi = chi,
             eta = eta,
             rho = c(0.25,0.5,0.75,0.9,0.95),
             M = seq(1,200),
             x = "M",
             eq = EO.A1,
             lbls = c("M","expected # of links"))

p3 <- plt.eq(chi = chi,
             eta = eta,
             rho = seq(0, 0.99, length.out = 100),
             M = c(10,50,100,200),
             x = "rho",
             eq = Pr.A1,
             lbls = c("rho","false discovery rate"),
             inverse = TRUE)

p4 <- plt.eq(chi = chi,
             eta = eta,
             rho = c(0.25,0.5,0.75,0.9,0.95),
             M = seq(1,200),
             x = "M",
             eq = Pr.A1,
             lbls = c("M","false discovery rate"),
             inverse = TRUE)

grid.title <- paste("specificity = ",chi," , sensitivity = ", eta, sep = "")
grid.arrange(p1,p2,p3,p4,nrow=2, top=grid.title)


###############################################################################

# A2
# Make four plots showing relationship between M, rho, false discovery rate, and expected number of links
# assumptions: single linkage, single transmission
# assumptions: fixed values of specificity and sensitivity

chi <- 0.85 # specificity of the linkage criteria
eta <- 0.9 # sensitivity of the linkage criteria

p1 <- plt.eq(chi = chi,
             eta = eta,
             rho = seq(0, 0.99, length.out = 100),
             M = c(10,50,100,200),
             x = "rho",
             eq = EO.A2,
             lbls = c("rho","expected # of links"))

p2 <- plt.eq(chi = chi,
             eta = eta,
             rho = c(0.25,0.5,0.75,0.9,0.95),
             M = seq(1,200),
             x = "M",
             eq = EO.A2,
             lbls = c("M","expected # of links"))

p3 <- plt.eq(chi = chi,
             eta = eta,
             rho = seq(0, 0.99, length.out = 100),
             M = c(10,50,100,200),
             x = "rho",
             eq = Pr.A2,
             lbls = c("rho","false discovery rate"),
             inverse = TRUE)

p4 <- plt.eq(chi = chi,
             eta = eta,
             rho = c(0.25,0.5,0.75,0.9,0.95),
             M = seq(1,200),
             x = "M",
             eq = Pr.A2,
             lbls = c("M","false discovery rate"),
             inverse = TRUE)

grid.title <- paste("specificity = ",chi," , sensitivity = ", eta, sep = "")
grid.arrange(p1,p2,p3,p4,nrow=2, top=grid.title)


###############################################################################

# A3
# Make four plots showing relationship between M, rho, false discovery rate, and expected number of links
# assumptions: single linkage
# assumptions: fixed values of specificity and sensitivity

chi <- 0.85 # specificity of the linkage criteria
eta <- 0.9 # sensitivity of the linkage criteria
R <- 1 # effective reproductive number

p1 <- plt.eq(chi = chi,
             eta = eta,
             R = R,
             rho = seq(0, 0.99, length.out = 100),
             M = c(10,50,100,200),
             x = "rho",
             eq = EO.A3,
             lbls = c("rho","expected # of links"))

p2 <- plt.eq(chi = chi,
             eta = eta,
             R = R,
             rho = c(0.25,0.5,0.75,0.9,0.95),
             M = seq(1,200),
             x = "M",
             eq = EO.A3,
             lbls = c("M","expected # of links"))

p3 <- plt.eq(chi = chi,
             eta = eta,
             R = R,
             rho = seq(0, 0.99, length.out = 100),
             M = c(10,50,100,200),
             x = "rho",
             eq = Pr.A3,
             lbls = c("rho","false discovery rate"),
             inverse = TRUE)

p4 <- plt.eq(chi = chi,
             eta = eta,
             R = R,
             rho = c(0.25,0.5,0.75,0.9,0.95),
             M = seq(1,200),
             x = "M",
             eq = Pr.A3,
             lbls = c("M","false discovery rate"),
             inverse = TRUE)

grid.title <- paste("specificity = ",chi," , sensitivity = ", eta, " , R = ", R, sep = "")
grid.arrange(p1,p2,p3,p4,nrow=2, top=grid.title)


###############################################################################

# Heatmap
# Make heatmap showing effect of changing sensitivity and specificity
chi <- seq(0.9,1,length.out = 100)
eta <- seq(0,1,length.out = 100)

plt.heatmap(chi,eta,rho=0.8,M=200,eq=Pr.A3,R=3)