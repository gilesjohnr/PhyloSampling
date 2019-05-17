### phylosampling plotting functions ###

###############################################################################

# A function producing a plot of three variables
# using the phylosampling function specified
# and fixed sensititvity and specficity

plt.eq <- function(chi,                 # number: specificity of the linkage criteria
                   eta,                 # number: sensitivity of the linkage criteria
                   R=1,                 # [optional] number: effective reproductive number
                   rho,                 # vector: values of rho to evaluate
                   M,                   # vector: values of M to evaluate
                   x="rho",             # string: which variable to put on the x axis
                   eq,                  # string: phylosampling function to evaluate
                   lbls=c("",""),       # labels for plot as: c("xlab","ylab")
                   inverse=FALSE,       # [optional] TRUE to plot 1 minus result of equation
                   legend=TRUE          # [optional] TRUE to show legend to the right of the plot
){
  
  # set up the dataframe to be used in plotting
  
  if (x == "rho"){
    
    g <- expand.grid(rho)
    names(g) <- c('x')
    
    for (i in seq(1,length(M))){
      cname <- paste("M=",M[i],sep="") # set name for column to be added
      if (inverse == FALSE){
        g <- cbind(g, eq(chi, g$x, M[i], eta, R))
      }
      else {
        g <- cbind(g, 1-eq(chi, g$x, M[i], eta, R))
      }
      colnames(g)[length(colnames(g))] <- cname
    }
    
  }
  
  else if (x == "M"){
    
    g <- expand.grid(M)
    names(g) <- c('x')
    
    for (i in seq(1,length(rho))){
      cname <- paste("rho=",rho[i],sep="") # set name for column to be added
      if (inverse == FALSE){
        g <- cbind(g, eq(chi, rho[i], g$x, eta, R))
      }
      else {
        g <- cbind(g, 1-eq(chi, rho[i], g$x, eta, R))
      }
      colnames(g)[length(colnames(g))] <- cname
    }
    
  }
  
  else {
    return("Error: x axis variable must be either rho or M")
  }

  # set up the plot
  
  melted.g <- melt(g, id = 'x')
  ggplot(melted.g, aes(x = x, y = value, colour = variable)) +
    geom_line(show.legend = legend) +
    xlab(lbls[1]) +
    ylab(lbls[2])
}

###############################################################################

# A function producing a heatmap of the false discovery rate
# for different values of sensititvity and specficity

plt.heatmap <- function(chi,            # vector: specificity of the linkage criteria
                        eta,            # vector: sensitivity of the linkage criteria
                        R=0,            # number: effective reproductive number
                        rho,            # number: sampling proportion
                        M,              # number: sample size
                        eq              # string: phylosampling function to evaluate
){
  
  g <- expand.grid(chi,eta)
  names(g) <- c('chi','eta')
  
  g <- cbind(g, 1-eq(chi = g$chi, eta = g$eta, rho = rho, M = M, R = R))
  colnames(g)[length(colnames(g))] <- "FDR"
  
  levelplot(FDR ~ chi*eta, data = g,
            col.regions = sequential_hcl(100)[length(sequential_hcl(100)):1])
}

