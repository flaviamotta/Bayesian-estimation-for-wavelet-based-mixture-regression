# Dynamic Mixture Weights Estimation in Gaussian Mixture Models
This repository contains the code and implementation of the methods discussed in the paper titled "Dynamic Mixture Weights Estimation in Gaussian Mixture Models". The goal of this work is to relax the traditional assumption that the mixture weights in Gaussian mixture models (GMM) remain constant over time or space, which limits their applicability to more complex, real-world datasets. Our approach allows the mixture weights to vary dynamically, which is particularly useful in scenarios like time-series data, space-dependent data, or medical data, where the mixture weights change over time or spatial dimensions.

# Key Contributions

- **Dynamic Mixture Weights:** We propose a Bayesian framework for estimating dynamic mixture weights in a two-component Gaussian mixture model (GMM).
- **Bayesian Inference:** The method uses Bayesian inference, specifically Gibbs sampling, to estimate both the dynamic mixture weights and the parameters of the Gaussian components.
- **Wavelet Regularization:** The method leverages wavelet-based regularization to capture the time- or space-dependence of the mixture weights.
- **Applications to GBM Data:** The approach is demonstrated on simulated data as well as real tumor data from Glioblastoma Multiforme (GBM) studies, successfully identifying chromosomal abnormalities associated with GBM presence.
