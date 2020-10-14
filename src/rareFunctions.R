indicator_greater <- function(xx, threshold=3, regularizer=1e-5)
{
  ##########################################################################
  #
  # INCIDATOR GREATER THAN FUNCTION
  #
  # f(x) = I(x > threshold)
  # 
  # Strictly 1D
  #
  ##########################################################################
  #
  # INPUTS:
  #
  # xx = c(x1)
  # threshold = Parameter threshold in integrand. Scalar.
  #      If not given, threshold=3
  #
  ##########################################################################
  
  return(as.numeric(xx > threshold) + regularizer)
}

indicator_square_greater <- function(xx, threshold=3, regularizer=1e-5)
{
  ##########################################################################
  #
  # INCIDATOR X^2 GREATER THAN FUNCTION
  #
  # f(x) = I(x^2 > threshold)
  # 
  # Strictly 1D
  #
  ##########################################################################
  #
  # INPUTS:
  #
  # xx = c(x1)
  # threshold = Parameter threshold in integrand. Scalar.
  #      If not given, threshold=3
  #
  ##########################################################################
  
  return(as.numeric(xx^2 > threshold) + regularizer)
}

portfolio_loss <- function(xx, regularizer=1e-5) {
  ##########################################################################
  #
  # PORTFOLIO LOSS FUNCTION
  #
  # f(x) = c1*I(X1>x1) + c2*I(X2>x2) with X1, X2 ~ exp(1)
  # 
  # Strictly 2D
  # Source: https://pat-laub.github.io/rare-events/
  #
  ##########################################################################
  #
  # INPUTS:
  #
  # xx = c(x1)
  # threshold = Parameter threshold in integrand. Scalar.
  #      If not given, threshold=3
  #
  ##########################################################################
  
  # c1 <- exp(0.2)
  # c2 <- exp(0.4)
  
  # x1 <- exp(0.2)
  # x2 <- exp(0.4)


  if (is.matrix(xx) == FALSE) { 
    xx <- matrix(xx, nrow = 1)    
  }

  dim <- ncol(xx)
  a <- matrix(rep(exp(0.2 * (0:(dim - 1))), each = nrow(xx)), ncol = dim)
  c <- matrix(rep(exp(0.2 * (0:(dim - 1))), each = nrow(xx)), ncol = dim)

  return(rowSums(c * (xx > a)) + regularizer)
}
