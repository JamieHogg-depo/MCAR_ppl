## Simple ICAR model 

# Load packages
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from plotnine import *
import geopandas as gpd
import libpysal as lps
import arviz as az
import pymc as pm
from pymc.distributions import continuous, distribution
import scipy.stats as st

## Shapefile and weight matrix ## ---------------------------------------------

aus_map = gpd.read_file("C:/Users/n9401849/OneDrive - Queensland University of Technology/DataLab/Requested file load/2016/2016_SA2_Shape/SA2_2016_AUST.shp")
aus2 = aus_map.loc[aus_map.STATE_NAME == "Australian Capital Territory"]
aus2 = aus2.reset_index(drop=True)
act_map = aus2.iloc[0:131,:]
qW = lps.weights.Queen.from_dataframe(act_map)
W = qW.full()[0]

## Add observations from the MCAR ## ------------------------------------------

# Fix the value of rho
rho = 0.9
M = W.shape[1]
K = 2

# get rowSums
wplus = W.sum(axis=1) #  axis = 1 is rows

# get L matrix
L = np.diag(1-rho + rho*wplus)

# Create identity matrix
I = np.identity(M)
D = np.diag(wplus)
C = I - D + W

# get eigen values
lambda_C, Q_c = np.linalg.eig(C)
ldet_C = np.log(1-rho*lambda_C)

# Correlation of two RVs
Cor_R = np.array([[1.0, 0.9],
                 [0.9, 1.0]])

# Convert correlation matrix to covariance matrix
Sum_R = np.outer(np.sqrt(np.diag(Cor_R)), np.sqrt(np.diag(Cor_R))) * Cor_R
# get lower cholesky
Sum_R_chol = np.linalg.cholesky(Sum_R)
log_Sum_R_chol_d = np.log(np.diag(Sum_R_chol))
# convert to precison
Omega_R = np.linalg.inv(Sum_R)

# Get the Omega precision matrix
#Omega_S = L - rho * W
Omega_S = I - rho * C
Omega_A = np.kron(Omega_S, Omega_R)
# Omega is a 262 x 262 prevision matrix
Sum_A = np.linalg.inv(Omega_A)

# Draw input vector for the 131 ACT SA2s
y_vec = np.random.multivariate_normal(np.zeros(Sum_A.shape[1]), cov = Sum_A)
y_mat = y_vec.reshape((M,K)) #byrow by default

# Get log probabililty
# pymc
mvn = pm.MvNormal.dist(mu=np.zeros(262), cov=Sum_A)
pm.logp(mvn, y_vec).eval()

# scipy
st.multivariate_normal.logpdf(x = y_vec, mean=np.zeros(262), cov=Sum_A)

# Fast implementation
A_S = np.matmul(Omega_S, y_mat)
A_R = np.matmul(Omega_R, np.transpose(y_mat))
-(M*K/2)*np.log(2*np.pi) + 0.5*(-2*M*sum(log_Sum_R_chol_d) + K * sum(ldet_C)) - 0.5 * np.matmul(A_R, A_S).trace()

## Plots ## -------------------------------------------------------------------

# map the simulated data
act_map['y1'] = y_mat[:,0]
act_map['y2'] = y_mat[:,1]

# map modelled estimates
fig, axs = plt.subplots(1, 2, figsize=(10, 4))
act_map.plot(column="y1", cmap="OrRd", ax=axs[0]); axs[0].set_title("y1")
act_map.plot(column="y2", cmap="OrRd", ax=axs[1]); axs[1].set_title("y2")
# Adjust spacing between subplots
plt.subplots_adjust(wspace=0.3); plt.show()

## Fit MCAR model in Pymc ## --------------------------------------------------

# Setup model
with pm.Model() as model0:

    # prior on rho
    rho = pm.Uniform("rho", 0, 1)
    sd_dist = pm.Exponential.dist(0.5, shape=(2,))

    # get back standard deviations and rho:
    chol, corr, stds = pm.LKJCholeskyCov("chol", n=2, eta=2.0, sd_dist=sd_dist)
    # LKJ prior is a distribution on the correlation matrix which combined with priors
    # on standard deviations induces a prior on the covaiance matrix
    # it uses the Cholesky decomposition sigma = L L^T where L is lower triangular
    Sum_R = pm.Deterministic("Sum_R", chol.dot(chol.T))
    Omega_R = pm.math.matrix_inverse(Sum_R)

    # construct precision matrix
    Sum_S = L - rho * W
    Omega_A = pm.math.kronecker(Omega_S, Omega_R)

    # Likelihood
    Y_obs = pm.MvNormal('y_vec', mu=np.zeros(262), tau=Omega_A, observed=y_vec)

    # run NUTS
    trace0 = pm.sample(draws = 1000, tune = 1000)

az.summary(trace0, round_to=2)















# Add observations
act_map.y = np.random.normal(0,1,131)
act_map.plot(column="y", cmap = "OrRd"); plt.show()

# test priors for rho
x = np.linspace(0,1,1000)
y = st.beta.pdf(x, 6, 1.25)
plt.plot(x, y); plt.show()

## Model 0: Fit the pCAR model ## ---------------------------------------------

with pm.Model() as model0:

    # Priors for unknown model parameters
    alpha_scalar = pm.Normal("alpha", mu=0, sigma=10)
    alpha = np.repeat(alpha_scalar, 131)

    # Spatial priors
    tau = pm.Gamma("tau", alpha=2.0, beta=2.0)
    rho = pm.Beta("rho", 6, 1.25)
    theta = pm.CAR("theta", mu=alpha, W=act_W, tau=tau, alpha = 0.99)

    # Likelihood (sampling distribution) of observations
    Y_obs = pm.Normal("y", mu=theta, sigma=1, observed=act_map.y)

    # Sample
    trace0 = pm.sample(draws = 1000, tune = 1000)
    rfvi = pm.fit(method="fullrank_advi")
    trace_rf = pm.sample_approx(rfvi, draws = 4000)
    approx = pm.Empirical(trace0)
    # use https://www.pymc.io/projects/examples/en/latest/variational_inference/empirical-approx-overview.html

# Summarize 
az.summary(trace0, round_to=2)
az.summary(trace_rf, round_to=2)

# Trace plots
az.plot_trace(trace0, ['tau']); plt.show()
az.plot_trace(approx_s, ['tau']); plt.show()

# get modelled estimates
dic0  = trace0.to_dict()['posterior']
y_mean = dic0['theta'].mean(axis=(0,1))
act_map.loc[:,"sm"] = y_mean

# map modelled estimates
fig, axs = plt.subplots(1, 2, figsize=(10, 4))
act_map.plot(column="y", cmap="OrRd", ax=axs[0]); axs[0].set_title("Observed data")
act_map.plot(column="sm", cmap="OrRd", ax=axs[1]); axs[1].set_title("Smoothed")
# Adjust spacing between subplots
plt.subplots_adjust(wspace=0.3); plt.show()
