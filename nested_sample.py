import numpy as np
import nestle
from scipy.stats import beta
from scipy.special import gammaln
import corner
#import matplotlib.pyplot as plt
import pandas as pd
import time
import os

os.chdir("Q:/Technical/Python/")

# individual invasiveness evidence 

data = pd.read_csv("Q:/Technical/Python/abschildnew.csv")
logz_df = pd.DataFrame(data = np.nan, index = range(0,len(data)), columns = ["logz", "logzerr"])

for i in range(len(data)):
    df = data.loc[i]
    carriage = df['carriage']
    disease = df['disease']
    n_swab = df['n.swab']
    N = df['N']
    time_int = df['time.int']
    DS = df['DS']
    activ_pts = 500
    
    # Define a likelihood function
    def loglike(theta):
        inv, carrprev = theta
        lambd = inv*carrprev*N*time_int
        LL = disease*np.log(lambd) - lambd - gammaln(disease+1)
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
    result = nestle.sample(loglike, prior_transform, ndim = 2, npoints = activ_pts)
    infogainnestle = result.h
    logZerrnestle = np.sqrt(infogainnestle/activ_pts)

    logz_df.loc[i, ['logz']] = result.logz     # log evidence
    logz_df.loc[i, ['logzerr']] = logZerrnestle
    
    #result.logzerr  # numerical (sampling) error on logz
    #result.samples  # array of sample parameters
    #result.weights  # array of weights associated with each sample

    # weighted average and covariance:
    #p, cov = nestle.mean_and_cov(result.samples, result.weights)
    
    #print()
    #print("inv = {0:5.3f} +/- {1:5.3f}".format(p[0], np.sqrt(cov[0, 0])))
    #print("carrprev = {0:5.3f} +/- {1:5.3f}".format(p[1], np.sqrt(cov[1, 1])))
    
    #result.summary()
    
    #theta_true = [df.invasiveness,df['carr.prev']] # dummy var for invasiveness
    #fig = corner.corner(result.samples, weights=result.weights, labels=['Inv', 'Oxford'],
    #                    range=[0.99999, 0.99999], truths = theta_true, bins=30, label_kwargs=dict(fontsize=20))
    #plt.show()
    #print(result.summary())
    
    # Theoretical answer...
    #evidence_theory = (n_swab + 1)/(carriage * N * time_int * 0.5)
    #print("Log evidence from theory: ", np.log(evidence_theory))

#data['logz'] = logz
newdata = pd.concat([data, logz_df], axis=1)
newdata.to_csv('specific_inv.csv')

# consolidated invasiveness evidence

child_prePCV_DS =  ["Alabama.pre.PCV", "Bogota.pre.PCV", "Caracas.pre.PCV", "Czech.pre.PCV", "Goroka.pre.PCV", "Oxford.pre.PCV", "Sweden.pre.PCV", "E.W.pre.PCV", "Atlanta.pre.PCV", "Finland.pre.PCV",
              "Ontario.pre.PCV", "Morocco.pre.PCV"]

childrendat = data[data.DS.isin(child_prePCV_DS)].copy()
childrendat.columns = ['agegrp', 'DS', 'Serogroup', 'Serotype', 'carriage', 'disease',
       'n_swab', 'N', 'time_int', 'carr_prev', 'carr_prev_low',
       'carr_prev_high', 'lambda', 'lambda_low', 'lambda_high', 'invasiveness',
       'invasiveness_low', 'invasiveness_high', 'logz']

consol_logz = np.empty(len(childrendat)) #pd.DataFrame(data = np.nan, index = range(0,len(childrendat)), columns = ["consol_logz"])
consol_logz[:] = np.nan

childrendat['consol_logz'] = consol_logz
#adultdat = data[data.DS.isin(adult_prePCV_DS)]

unique_sero_child = childrendat.Serotype.unique()
#unique_sero_adult = adultdat.Serotype.unique()

evidence_df = pd.DataFrame(data = np.nan, index = range(len(unique_sero_child)), columns = ["Serotype", "evidence_indiv", "evidence_consol"])
evidence_df['Serotype'] = unique_sero_child

newsero = []
numDS = []
for i in unique_sero_child:
    sero_spec = childrendat.loc[childrendat['Serotype'] == i]
    n = len(sero_spec)
    if (n < 9 and n > 1):
        newsero.append(i)
        numDS.append(n)
        
for i in newsero: #unique_sero_child:
    sero_spec = childrendat.loc[childrendat['Serotype'] == i].reset_index(drop = True)
    start = time.time()
    n = len(sero_spec)
    if n > 9:
        activ_pts = 200
    else:
        activ_pts = 500
    
    def loglike(theta):
        LL = np.empty(n)
        inv = theta[0]
        carrprev = theta[1:n+1]
        
        for j in range(n):
            lambd = inv*carrprev[j]*sero_spec.N[j]*sero_spec.time_int[j]
            LL[j] = sero_spec.disease[j]*np.log(lambd) - lambd - gammaln(sero_spec.disease[j]+1)
            
        return sum(LL)
    
    def prior_transform(theta):
        inv_prime = theta[0]
        carrprev_prime = theta[1:n+1]
        
        invmin = 0.
        invmax = 0.5
        inv = inv_prime*(invmax-invmin) + invmin
        
        carrprev_alpha = np.empty(n)
        carrprev_beta = np.empty(n)
        carrprev = np.empty(n)
        
        for j in range(n):
            carrprev_alpha[j] = sero_spec.carriage[j] + 1
            carrprev_beta[j] = sero_spec.n_swab[j] - sero_spec.carriage[j] + 1
            carrprev[j] = beta.ppf(carrprev_prime[j], carrprev_alpha[j], carrprev_beta[j])
    
        return np.append(inv, carrprev)
        
    result_cons = nestle.sample(loglike, prior_transform, ndim = n+1, npoints = activ_pts)
    end = time.time()
    result_cons.logz
    
    logz = result_cons.logz  # log evidence
    #result.logzerr  # numerical (sampling) error on logz
    #result.samples  # array of sample parameters
    #result.weights  # array of weights associated with each sample

    # weighted average and covariance:
    p, cov = nestle.mean_and_cov(result_cons.samples, result_cons.weights)
    
    #print()
    #print("inv = {0:5.3f} +/- {1:5.3f}".format(p[0], np.sqrt(cov[0, 0])))
    #print("carrprev = {0:5.3f} +/- {1:5.3f}".format(p[1], np.sqrt(cov[1, 1])))
    
    #result.summary()
    
    theta_true = np.ndarray.tolist(np.insert(sero_spec.carr_prev.to_numpy(), 0, 0.000253134))
    fig = corner.corner(result_cons.samples, weights=result_cons.weights, labels=['Inv', 'Bogota', 'Oxford', 'Eng & Wales'], 
                        range=[0.99999, 0.99999, 0.99999, 0.99999], truths = theta_true, bins=30, label_kwargs=dict(fontsize=20))
    #theta_true = [df.invasiveness,df['carr.prev']] # dummy var for invasiveness
    #fig = corner.corner(result.samples, weights=result.weights, labels=['inv', 'carrprev'],
    #                    range=[0.99999, 0.99999], truths = theta_true, bins=30)
    #plt.show()
    #print(result.summary())
    
    #childrendat.consol_logz[childrendat['Serotype'] == i] = result_cons.logz
    evidence_df.loc[evidence_df['Serotype'] == i,['evidence_indiv']] = sum(childrendat.logz[childrendat['Serotype'] == i])
    evidence_df.loc[evidence_df['Serotype'] == i,['evidence_consol']] = result_cons.logz
    
    print(end-start)

evidence_df['bayesfac'] = evidence_df['evidence_indiv'] - evidence_df['evidence_consol']

