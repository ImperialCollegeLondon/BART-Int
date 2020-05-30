# We compute the exact integrals for all the genz functions for dimensions
# 1,2,3,5,10,20 and then store them as CSV format.
# We use the notation from https://www.sfu.ca/~ssurjano/disc.html,
# thus alpha corresponds to the a's, and beta corresponds to the b's.
# We fix u = [0.5,....]
# and choose a according to the rescaling scheme (see paper).
# The range of integration is from 0 to 1.

source("src/genz/analyticalIntegrals.R")
source("src/genz/genz.R")
library("cubature")

numGenz <- 6
dimensions <- c(1, 2, 3, 5, 10, 20)
integrals <- matrix(rep(NA, length(dimensions) * numGenz), nrow = numGenz)
integralsCubature <- matrix(rep(NA, length(dimensions) * numGenz), nrow = numGenz)

for (k in 1:length(dimensions)){

    # Compute integrals
    dim <- dimensions[k]
    integrals[1, k] <- contIntegral(dim)
    integrals[2, k] <- copeakIntegral(dim)
    integrals[3, k] <- discIntegral(dim)
    integrals[4, k] <- gaussianIntegral(dim)
    integrals[5, k] <- oscillatoryIntegral(dim)
    integrals[6, k] <- productPeakIntegral(dim)
    # Discontinuous integrand is only defined for dim >= 2
    if (dim > 1){
        integrals[3, k] <- discIntegral(dim)
    }

    cat("dim:", dim, "completed. \n")

}

# Write results to file
write.table(integrals, file = "results/genz/integrals.csv", sep=",", row.names=FALSE, col.names=FALSE)

