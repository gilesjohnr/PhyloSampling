source('R/load.R')
source('R/func.R')

# Check first generation function of Pr(A=A*|L(a,B))
Pr.TL(x=0.99, rho=0.99, N=1001)
Pr.TL(x=0.99, rho=0.9, N=11)

# Explore
rho.vals <- seq(0, 0.99, length.out=50)
N.vals <- seq(1000, 1e+05, length.out=50)

m <- mesh(rho.vals, N.vals)
Pr.vals <- with(m, Pr.TL(x=0.99, rho=x, N=y))

persp3D(x=rho.vals, y=N.vals, z=Pr.vals,
        theta = 210,
        phi = 30,
        colkey=F,
        border='black',
        col=ramp.col (col = c("red", "darkcyan", "cyan", "lightcyan", "greenyellow"), 
                      n = 100, alpha = 1),
        ticktype="detailed",
        zlab='Pr( T(A,B) | L(A,B) )',
        xlab='rho',
        ylab='N',
        main='x = 0.99')

#############################################################
# Check second generation function of Pr(A=A*|L(a,B)) that uses sensitivity
Pr.TL2(N=1001, rho=0.99, phi=0.99)
Pr.TL2(N=11, rho=0.9, phi=0.99)

# Contour
g <- expand.grid(rho.vals, phi.vals)
names(g) <- c('rho', 'phi')
g$Pr <- Pr.TL2(N=1000, rho=g$rho, phi=g$phi)

p <- ggplot(g, aes(x=phi, y=rho, z=Pr))

p + stat_contour()

# Heatmap and contour
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