### phylosampling plotting functions ###

###############################################################################

# A function producing a plot of three variables
# using the phylosampling function specified
# and fixed sensititvity and specficity

plt.eq <- function(chi=0.95,            # number: specificity of the linkage criteria
                   eta=0.95,            # number: sensitivity of the linkage criteria
                   rho,                 # vector: values of rho to evaluate
                   M,                   # vector: values of M to evaluate
                   x="rho",             # string: which variable to put on the x axis
                   eq,                  # string: phylosampling function to evaluate
                   lbls=c("",""),       # labels for plot as: c("xlab","ylab")
                   inverse=FALSE        # TRUE to plot 1 minus result of equation
){
  
  # set up the dataframe to be used in plotting
  
  if (x == "rho"){
    
    g <- expand.grid(rho)
    names(g) <- c('x')
    
    for (i in seq(1,length(M))){
      cname <- paste("M=",M[i],sep="") # set name for column to be added
      if (inverse == FALSE){
        g <- cbind(g, eq(chi, g$x, M[i]))
      }
      else {
        g <- cbind(g, 1-eq(chi, g$x, M[i]))
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
        g <- cbind(g, eq(chi, rho[i], g$x))
      }
      else {
        g <- cbind(g, 1-eq(chi, rho[i], g$x))
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
    geom_line() +
    xlab(lbls[1]) +
    ylab(lbls[2])
}