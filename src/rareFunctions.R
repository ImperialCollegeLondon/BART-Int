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
