import numpy as np
import nestle
from scipy.stats import beta
#from scipy.stats import poisson
from scipy.stats import norm
import corner
import matplotlib.pyplot as plt
import pandas as pd
import math

data = pd.read_csv("Q:/Technical/Python/abschildnew.csv")

df = data.loc[0]

carriage = df['carriage']
disease = df['disease']
n_swab = df['n.swab']
N = df['N']
time_int = df['time.int']

# Define a likelihood function
def loglike(theta):
    inv, carrprev = theta
    
    lambd = inv*carrprev*N*time_int
    #L = poisson.pmf(k = disease, mu = lambd)
    #LL = np.log(L)
    
    L = norm.pdf(x = disease, loc = lambd, scale = math.sqrt(lambd))
    LL = np.log(L)
    
    return LL

# Define a function mapping the unit cube to the prior space.
def prior_transform(theta):
    inv_prime, carrprev_prime = theta
    
    invmin = 0.
    invmax = 0.5
    inv = inv_prime*(invmax-invmin) + invmin
    
    carrprev_alpha = carriage + 1
    carrprev_beta = n_swab - carriage + 1
    carrprev = beta.ppf(carrprev_prime, carrprev_alpha, carrprev_beta)
    
    return np.array([inv, carrprev])

# Run nested sampling.
result = nestle.sample(loglike, prior_transform, ndim = 2, npoints = 500, method = 'classic')

#result.logz     # log evidence
#result.logzerr  # numerical (sampling) error on logz
#result.samples  # array of sample parameters
#result.weights  # array of weights associated with each sample

# weighted average and covariance:
p, cov = nestle.mean_and_cov(result.samples, result.weights)

print("inv = {0:5.3f} +/- {1:5.3f}".format(p[0], np.sqrt(cov[0, 0])))
print("carrprev = {0:5.3f} +/- {1:5.3f}".format(p[1], np.sqrt(cov[1, 1])))

theta_true = [df.invasiveness,df['carr.prev']] # dummy var for invasiveness
fig = corner.corner(result.samples, weights=result.weights, labels=['inv', 'carrprev'],
                    range=[0.99999, 0.99999], truths = theta_true, bins=30)
plt.show()
print(result.summary())