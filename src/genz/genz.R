# ================================================================
# Genz test functions. Input must be of type Matrix
# ================================================================

library(matrixStats)
library(MASS)
library(msm)

cont <- function(xx)
{
  ##########################################################################
  #
  # CONTINUOUS INTEGRAND FAMILY
  #
  # Authors: Sonja Surjanovic, Simon Fraser University
  #          Derek Bingham, Simon Fraser University
  # Questions/Comments: Please email Derek Bingham at dbingham@stat.sfu.ca.
  #
  # Copyright 2013. Derek Bingham, Simon Fraser University.
  #
  # THERE IS NO WARRANTY, EXPRESS OR IMPLIED. WE DO NOT ASSUME ANY LIABILITY
  # FOR THE USE OF THIS SOFTWARE.  If software is modified to produce
  # derivative works, such modified software should be clearly marked.
  # Additionally, this program is free software; you can redistribute it 
  # and/or modify it under the terms of the GNU General Public License as 
  # published by the Free Software Foundation; version 2.0 of the License. 
  # Accordingly, this program is distributed in the hope that it will be 
  # useful, but WITHOUT ANY WARRANTY; without even the implied warranty 
  # of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU 
  # General Public License for more details.
  #
  # For function details and reference information, see:
  # http://www.sfu.ca/~ssurjano/
  #
  ##########################################################################
  #
  # INPUTS:
  #
  # xx = c(x1, x2, ..., xd)
  # u = c(u1, u2, ..., ud) (optional), with default value
  #     c(0.5, 0.5, ..., 0.5)
  # a = 150/dim^3, where dim = number of column of xx
  #
  ##########################################################################
  
  if (is.matrix(xx) == FALSE) { 
    xx <- matrix(xx, nrow = 1)    
  }

  dim <- ncol(xx)
  a <- rep(150/dim^3, dim)  
  u <- rep(0.5, dim)

  sum <- abs(xx-u) %*% a
  y <- exp(-sum)
  
  return(y)
}


copeak <- function(xx)
{
  ##########################################################################
  #
  # CORNER PEAK INTEGRAND FAMILY
  #
  # Authors: Sonja Surjanovic, Simon Fraser University
  #          Derek Bingham, Simon Fraser University
  # Questions/Comments: Please email Derek Bingham at dbingham@stat.sfu.ca.
  #
  # Copyright 2013. Derek Bingham, Simon Fraser University.
  #
  # THERE IS NO WARRANTY, EXPRESS OR IMPLIED. WE DO NOT ASSUME ANY LIABILITY
  # FOR THE USE OF THIS SOFTWARE.  If software is modified to produce
  # derivative works, such modified software should be clearly marked.
  # Additionally, this program is free software; you can redistribute it 
  # and/or modify it under the terms of the GNU General Public License as 
  # published by the Free Software Foundation; version 2.0 of the License. 
  # Accordingly, this program is distributed in the hope that it will be 
  # useful, but WITHOUT ANY WARRANTY; without even the implied warranty 
  # of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU 
  # General Public License for more details.
  #
  # For function details and reference information, see:
  # http://www.sfu.ca/~ssurjano/
  #
  ##########################################################################
  #
  # INPUTS:
  #
  # xx = c(x1, x2, ..., xd)
  # u = c(u1, u2, ..., ud) (optional), with default value
  #     c(0.5, 0.5, ..., 0.5)
  # a = 600/dim^3, where dim = number of column of xx
  #
  #########################################################################
  
  if (is.matrix(xx) == FALSE) { 
    xx <- matrix(xx, nrow = 1) 
  }

  dim <- ncol(xx)
  u <- rep(0.5, dim)

  a <- matrix(rep(600/dim^3, dim), ncol = 1)
  
  sum <- xx %*% a
  y <- (1 + sum)^(-dim - 1)
  
  return(y)
}


disc <- function(xx)
{
  ##########################################################################
  #
  # DISCONTINUOUS INTEGRAND FAMILY
  #
  # Authors: Sonja Surjanovic, Simon Fraser University
  #          Derek Bingham, Simon Fraser University
  # Questions/Comments: Please email Derek Bingham at dbingham@stat.sfu.ca.
  #
  # Copyright 2013. Derek Bingham, Simon Fraser University.
  #
  # THERE IS NO WARRANTY, EXPRESS OR IMPLIED. WE DO NOT ASSUME ANY LIABILITY
  # FOR THE USE OF THIS SOFTWARE.  If software is modified to produce
  # derivative works, such modified software should be clearly marked.
  # Additionally, this program is free software; you can redistribute it 
  # and/or modify it under the terms of the GNU General Public License as 
  # published by the Free Software Foundation; version 2.0 of the License. 
  # Accordingly, this program is distributed in the hope that it will be 
  # useful, but WITHOUT ANY WARRANTY; without even the implied warranty 
  # of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU 
  # General Public License for more details.
  #
  # For function details and reference information, see:
  # http://www.sfu.ca/~ssurjano/
  #
  ##########################################################################
  #
  # INPUTS:
  #
  # xx = c(x1, x2, ..., xd)
  # u = c(u1, u2, ..., ud) (optional), with default value
  #     c(0.5, 0.5, ..., 0.5)
  # a = 100/dim^3, where dim = number of column of xx
  #
  #########################################################################
  
  if (is.matrix(xx) == FALSE) { 
    xx <- matrix(xx, nrow = 1)  
  }
  
  # Function only defined for dimension >= 2
  dim <- ncol(xx)
  if (dim == 1) {
    u <- rep(0.5, dim)
    a <- rep(10/dim^3, dim)
    
    x1 <- xx[ ,1]
    u1 <- u[1]
    
    xx[which(x1 > u1), ] <- 0
    
    sum <- xx * a
    y <- exp(sum)
    y[which (y == 1)] <- 0
    
    return(y)
  }
  
  u <- rep(0.5, dim)
  a <- rep(10/dim^3, dim)
  
  x1 <- xx[ ,1]
  x2 <- xx[ ,2]
  u1 <- u[1]
  u2 <- u[2]

  xx [which ( x1 > u1 | x2 > u2), ] <- 0

  sum <- xx %*% a
  y <- exp(sum)
  y[which (y == 1)] <- 0
  
  return(y)
  
}


gaussian <- function(xx)
{
  ##########################################################################
  #
  # GAUSSIAN PEAK INTEGRAND FAMILY
  #
  # Authors: Sonja Surjanovic, Simon Fraser University
  #          Derek Bingham, Simon Fraser University
  # Questions/Comments: Please email Derek Bingham at dbingham@stat.sfu.ca.
  #
  # Copyright 2013. Derek Bingham, Simon Fraser University.
  #
  # THERE IS NO WARRANTY, EXPRESS OR IMPLIED. WE DO NOT ASSUME ANY LIABILITY
  # FOR THE USE OF THIS SOFTWARE.  If software is modified to produce
  # derivative works, such modified software should be clearly marked.
  # Additionally, this program is free software; you can redistribute it 
  # and/or modify it under the terms of the GNU General Public License as 
  # published by the Free Software Foundation; version 2.0 of the License. 
  # Accordingly, this program is distributed in the hope that it will be 
  # useful, but WITHOUT ANY WARRANTY; without even the implied warranty 
  # of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU 
  # General Public License for more details.
  #
  # For function details and reference information, see:
  # http://www.sfu.ca/~ssurjano/
  #
  ##########################################################################
  #
  # INPUTS:
  #
  # xx = c(x1, x2, ..., xd)
  # u  = c(u1, u2, ..., ud) (optional), with default value
  #      c(0.5, 0.5, ..., 0.5)
  # a  = 100/dim^2, where dim = number of column of xx
  #
  ##########################################################################
  
  if (is.matrix(xx) == FALSE) { 
    xx <- matrix(xx, nrow = 1) 
  }
  
  dim <- ncol(xx)
  u <- rep(0.5, dim)
  a <- rep(100/dim^2, dim) 
  
  sum <- (xx - u)^2 %*% a^2
  y <- exp(-sum)
  
  return(y)
  
}

oscil <- function(xx)
{
  ##########################################################################
  #
  # OSCILLATORY INTEGRAND FAMILY
  #
  # Authors: Sonja Surjanovic, Simon Fraser University
  #          Derek Bingham, Simon Fraser University
  # Questions/Comments: Please email Derek Bingham at dbingham@stat.sfu.ca.
  #
  # Copyright 2013. Derek Bingham, Simon Fraser University.
  #
  # THERE IS NO WARRANTY, EXPRESS OR IMPLIED. WE DO NOT ASSUME ANY LIABILITY
  # FOR THE USE OF THIS SOFTWARE.  If software is modified to produce
  # derivative works, such modified software should be clearly marked.
  # Additionally, this program is free software; you can redistribute it 
  # and/or modify it under the terms of the GNU General Public License as 
  # published by the Free Software Foundation; version 2.0 of the License. 
  # Accordingly, this program is distributed in the hope that it will be 
  # useful, but WITHOUT ANY WARRANTY; without even the implied warranty 
  # of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU 
  # General Public License for more details.
  #
  # For function details and reference information, see:
  # http://www.sfu.ca/~ssurjano/
  #
  ##########################################################################
  #
  # INPUTS:
  #
  # xx = c(x1, x2, ..., xd)
  # u  = c(u1, u2, ..., ud) (optional), with default value
  #      c(0.5, 0.5, ..., 0.5)
  # a  = 110/dim^(5/2), where dim = number of column of xx
  #
  ##########################################################################
  
  if (is.matrix(xx) == FALSE) { 
    xx <- matrix(xx, nrow = 1) 
  }
  
  dim <- ncol(xx)
  u <- rep(0.5, 1, dim)
  a <- rep(110/dim^(5/2), dim)
  
  sum <- xx %*% a
  y <- cos(2 * pi * u[1] + sum)
  
  return(y)
  
}



prpeak <- function(xx)
{
  ##########################################################################
  #
  # PRODUC PEAK INTEGRAND FAMILY
  #
  # Authors: Sonja Surjanovic, Simon Fraser University
  #          Derek Bingham, Simon Fraser University
  # Questions/Comments: Please email Derek Bingham at dbingham@stat.sfu.ca.
  #
  # Copyright 2013. Derek Bingham, Simon Fraser University.
  #
  # THERE IS NO WARRANTY, EXPRESS OR IMPLIED. WE DO NOT ASSUME ANY LIABILITY
  # FOR THE USE OF THIS SOFTWARE.  If software is modified to produce
  # derivative works, such modified software should be clearly marked.
  # Additionally, this program is free software; you can redistribute it 
  # and/or modify it under the terms of the GNU General Public License as 
  # published by the Free Software Foundation; version 2.0 of the License. 
  # Accordingly, this program is distributed in the hope that it will be 
  # useful, but WITHOUT ANY WARRANTY; without even the implied warranty 
  # of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU 
  # General Public License for more details.
  #
  # For function details and reference information, see:
  # http://www.sfu.ca/~ssurjano/
  #
  ##########################################################################
  #
  # INPUTS:
  #
  # xx = c(x1, x2, ..., xd)
  # u  = c(u1, u2, ..., ud) (optional), with default value
  #      c(0.5, 0.5, ..., 0.5)
  # a  = 600/dim^3, where dim = number of column of xx
  #
  ##########################################################################

  if (is.matrix(xx) == FALSE) { 
    xx <- matrix(xx, nrow = 1)  
  }
  
  dim <- ncol(xx)
  u <- rep(0.5, dim)
  a <- rep(600/dim^3, dim)
  
  sum <- a^(-2) + (xx - u)^2
  y <- rowProds(1/sum)
  
  return(y)
}


step <- function(xx, regularizer=1e-5, jumps)
{
  ##########################################################################
  #
  # STEP FUNCTION
  #
  # f(x) = 0,             if x_1 < 1/(jumps + 1);
  #        j/(jumps + 1), if j/(jumps + 1) < x_1 <= (j + 1)/(jumps + 1),
  #                       for j = 1, ..., jumps
  # where x \in \R^d
  #
  ##########################################################################
  #
  # INPUTS:
  #
  # xx = c(x1, x2, ..., xd)
  # jumps = number of jumps
  #
  ##########################################################################
  if (is.matrix(xx) == FALSE) { 
    xx <- matrix(xx, nrow = 1)  
  }
  
  dim <- ncol(xx)
  
  y <- matrix(NA, nrow = dim(xx)[1])

  for (j in 2:jumps){
    y[xx[, 1] > (j - 1) / (jumps + 1) & xx[, 1] <= j / (jumps + 1), ] <- j / (jumps + 1)
  }
  y[xx[, 1] <= 1/(jumps + 1), ] <- regularizer
  y[xx[, 1] > jumps/(jumps + 1), ] <- 1 - regularizer

  return(y)
}


mix <- function(xx){
  
  ##########################################################################
  #
  # MIXTURE GENZ FUNCTION
  #
  # yi = copeak(xi) for xi <= 0.5, cont(xi) otherwise
  # y = sum(yi)
  #
  ##########################################################################
  
  
  
  if (is.matrix(xx) == FALSE) { 
    xx <- matrix(xx, nrow = 1)  
  }
  
  dim <- ncol(xx)
  
  yval <- function(x, dim){
    
    indices_0 <- which(x < 0.5)
    indices_1 <- which(x >= 0.5)
    
    y0 <- sum(apply(as.matrix(x[indices_0]), 1, copeak))
    y1 <- sum(apply(as.matrix(x[indices_1]), 1, cont))
    
    y <- y0+y1
    
    return(y)
  }
  # Compute Genz function values
  y <- apply(xx, 1, yval)
  
  return(as.matrix(y))
}


gaussian_weighted <- function(xx)
{
  ##########################################################################
  #
  # WEIGHTED GAUSSIAN PEAK INTEGRAND FAMILY
  #
  # Authors: Sonja Surjanovic, Simon Fraser University
  #          Derek Bingham, Simon Fraser University
  # Questions/Comments: Please email Derek Bingham at dbingham@stat.sfu.ca.
  #
  # Copyright 2013. Derek Bingham, Simon Fraser University.
  #
  # THERE IS NO WARRANTY, EXPRESS OR IMPLIED. WE DO NOT ASSUME ANY LIABILITY
  # FOR THE USE OF THIS SOFTWARE.  If software is modified to produce
  # derivative works, such modified software should be clearly marked.
  # Additionally, this program is free software; you can redistribute it 
  # and/or modify it under the terms of the GNU General Public License as 
  # published by the Free Software Foundation; version 2.0 of the License. 
  # Accordingly, this program is distributed in the hope that it will be 
  # useful, but WITHOUT ANY WARRANTY; without even the implied warranty 
  # of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU 
  # General Public License for more details.
  #
  # For function details and reference information, see:
  # http://www.sfu.ca/~ssurjano/
  #
  ##########################################################################
  #
  # INPUTS:
  #
  # xx = c(x1, x2, ..., xd)
  # u  = c(u1, u2, ..., ud) (optional), with default value
  #      c(0.5, 0.5, ..., 0.5)
  # a  = 100/dim^2, where dim = number of column of xx
  #
  ##########################################################################
  
  return(gaussian(xx) / dtnorm(xx, mean = 0.5, lower = 0, upper = 1))
  
}


prpeak_weighted <- function(xx)
{
  ##########################################################################
  #
  # WEIGHTED PRODUC PEAK INTEGRAND FAMILY
  #
  # Authors: Sonja Surjanovic, Simon Fraser University
  #          Derek Bingham, Simon Fraser University
  # Questions/Comments: Please email Derek Bingham at dbingham@stat.sfu.ca.
  #
  # Copyright 2013. Derek Bingham, Simon Fraser University.
  #
  # THERE IS NO WARRANTY, EXPRESS OR IMPLIED. WE DO NOT ASSUME ANY LIABILITY
  # FOR THE USE OF THIS SOFTWARE.  If software is modified to produce
  # derivative works, such modified software should be clearly marked.
  # Additionally, this program is free software; you can redistribute it 
  # and/or modify it under the terms of the GNU General Public License as 
  # published by the Free Software Foundation; version 2.0 of the License. 
  # Accordingly, this program is distributed in the hope that it will be 
  # useful, but WITHOUT ANY WARRANTY; without even the implied warranty 
  # of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU 
  # General Public License for more details.
  #
  # For function details and reference information, see:
  # http://www.sfu.ca/~ssurjano/
  #
  ##########################################################################
  #
  # INPUTS:
  #
  # xx = c(x1, x2, ..., xd)
  # u  = c(u1, u2, ..., ud) (optional), with default value
  #      c(0.5, 0.5, ..., 0.5)
  # a  = 600/dim^3, where dim = number of column of xx
  #
  ##########################################################################

  return(prpeak(xx) / dtnorm(xx, mean = 0.5, lower = 0, upper = 1))
}


step_weighted <- function(xx, regularizer=1e-5, jumps)
{
  ##########################################################################
  #
  # WEIGHTED STEP FUNCTION
  #
  # f(x) = 0,             if x_1 < 1/(jumps + 1);
  #        j/(jumps + 1), if j/(jumps + 1) < x_1 <= (j + 1)/(jumps + 1),
  #                       for j = 1, ..., jumps
  # where x \in \R^d
  #
  ##########################################################################
  #
  # INPUTS:
  #
  # xx = c(x1, x2, ..., xd)
  # jumps = number of jumps
  #
  ##########################################################################

  return(step(xx, regularizer, jumps) / dtnorm(xx, mean = 0.5, lower = 0, upper = 1))
}


additive_gaussian <- function(xx, a=NA)
{
  ##########################################################################
  #
  # ADDITIVE GAUSSIAN FUNCTION
  #
  # f(x) = \sum_{j=1}^{dim} gaussian(xx_j, a_j)
  # 
  # where xx_j is the jth column of xx, and a_j is the jth entry of a.
  #
  ##########################################################################
  #
  # INPUTS:
  #
  # xx = c(x1, x2, ..., xd)
  # a  = Parameter a in integrand. Can be a scalar or a vector of dimension dim. 
  #      If not given, a = 100/dim^2.
  #
  ##########################################################################

    if (is.matrix(xx) == FALSE) { 
        xx <- matrix(xx, nrow = 1) 
    }

    dim <- ncol(xx)
    u <- 0.5
    if ( any(is.na(a)) ){
        a <- 100/dim^2
    }

    if (length(a) == 1){
        a <- rep(a, dim)
    }

    sum <- 0
    for (j in 1:dim){
        sum <- sum + exp(-a[j]^2 * (matrix(xx[, j], ncol = 1) - u)^2)
    }

    return(sum)
}