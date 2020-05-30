library(reticulate)


optimise_gp_r <- function(trainX, trainY, kernel, epochs)
#'
#'
#'kernel == "rbf" or "matern"
#'
#'
{
  use_virtualenv("r-BOBART")
  source_python("python/gp_tune.py")
  lengthscale <- optimise_gp(trainX, trainY, kernel, epochs)
  return (lengthscale)
}

install_python_env <- function()
{
  # create a new environment 
  
  virtualenv_create("r-BOBART")
  virtualenv_install("r-BOBART", "gpytorch")
  virtualenv_install("r-BOBART", "torch")
  # import SciPy (it will be automatically discovered in "r-reticulate")
  use_virtualenv("r-BOBART")
}
