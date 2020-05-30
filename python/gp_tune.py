import os
os.environ["KMP_DUPLICATE_LIB_OK"]="TRUE"
import gpytorch 
from torch import Tensor, linspace, sin, randn
from torch.optim import Adam
import math
import os
os.environ["KMP_DUPLICATE_LIB_OK"] = "True"
# We will use the simplest form of GP model, exact inference
class ExactGPModel(gpytorch.models.ExactGP):
    def __init__(self, train_x, train_y, likelihood, kernel="rbf"):
        super(ExactGPModel, self).__init__(train_x, train_y, likelihood)
        self.mean_module = gpytorch.means.ZeroMean()
        if kernel == "rbf":
            self.covar_module = gpytorch.kernels.RBFKernel()
        elif kernel == "matern32":
            self.covar_module = gpytorch.kernels.MaternKernel(nu=1.5)
        else:
            raise NotImplementedError()

    def forward(self, x):
        mean_x = self.mean_module(x)
        covar_x = self.covar_module(x)
        return gpytorch.distributions.MultivariateNormal(mean_x, covar_x)


def optimise_gp(train_x, train_y, kernel, epochs):
    # initialize likelihood and model
    train_x = Tensor(train_x)
    train_y = Tensor(train_y)
    train_y = train_y.reshape(train_y.shape[0])
    # print(train_x.shape, train_y.shape)
    likelihood = gpytorch.likelihoods.GaussianLikelihood()
    model = ExactGPModel(train_x, train_y, likelihood, kernel)
    model.train()
    likelihood.train()
    # Use the adam optimizer
    optimizer = Adam([
        {'params': model.parameters()},  # Includes GaussianLikelihood parameters
    ], lr=0.1)

    # "Loss" for GPs - the marginal log likelihood
    mll = gpytorch.mlls.ExactMarginalLogLikelihood(likelihood, model)

    for i in range(int(epochs)):
        # Zero gradients from previous iteration
        optimizer.zero_grad()
        # Output from model
        output = model(train_x)
        # Calc loss and backprop gradients
        loss = -mll(output, train_y)
        loss.backward()
        if (i + 1) % 50 == 0:
            print('Iter %d/%d - Loss: %.3f   lengthscale: %.3f   noise: %.3f' % (
                i + 1, epochs, loss.item(),
                model.covar_module.lengthscale.item(),
                model.likelihood.noise.item()
            ))
        optimizer.step()
    return model.covar_module.lengthscale.item()