---
title: "`cvCovEst`: Cross-validated covariance matrix estimator selection and evaluation in `R`"
tags:
  - R
  - covariance matrix
  - cross-validation
  - high-dimensional statistics
  - loss-based estimation
  - multivariate analysis
authors:
  - name: Philippe Boileau
    orcid: 0000-0002-4850-2507
    affiliation: "1, 2"
  - name: Nima S. Hejazi
    orcid: 0000-0002-7127-2789
    affiliation: "1, 2"
  - name: Brian Collica
    orcid: 0000-0003-1127-2557
    affiliation: 3
  - name: Mark J. van der Laan
    orcid: 0000-0003-1432-5511
    affiliation: "2, 3, 4"
  - name: Sandrine Dudoit
    orcid: 0000-0002-6069-8629
    affiliation: "2, 3, 4"
affiliations:
  - name: Graduate Group in Biostatistics, University of California, Berkeley
    index: 1
  - name: Center for Computational Biology, University of California, Berkeley
    index: 2
  - name: Department of Statistics, University of California, Berkeley
    index: 3
  - name: Division of Biostatistics, School of Public Health, University of California, Berkeley
    index: 4
date: "\\today"
bibliography: paper.bib
---

# Summary

Covariance matrices play fundamental roles in myriad statistical procedures.
When the observations in a dataset far outnumber the features, asymptotic
theory and empirical evidence have demonstrated the sample covariance matrix to
be the optimal estimator of this parameter. This assertion does not hold when
the number of observations is commensurate with or smaller than the number of
features.  Consequently, statisticians have derived many novel covariance
matrix estimators for the high-dimensional regime, often relying on additional
assumptions about the parameter's structural characteristics (e.g., sparsity).
While these estimators have greatly improved the ability to estimate covariance
matrices in high-dimensional settings, objectively selecting the best estimator
from among the many possible candidates remains a largely unaddressed
challenge. The `cvCovEst` package addresses this methodological gap through its
implementation of a cross-validated framework for covariance matrix estimator
selection. This data-adaptive procedure's selections are asymptotically optimal
under minimal assumptions -- in fact, they are equivalent to the selections
that would be made if given full knowledge of the true
data-generating processes (i.e., an oracle selector) [@laan_dudoit:2003].

# Statement of Need

When the number of observations in a dataset far exceeds the number of
features, the estimator of choice for the covariance matrix is the sample
covariance matrix. It is efficient under minimal regularity assumptions on the
data-generating distribution. In high-dimensional regimes, however, its
performance is unsatisfactory: the sample covariance matrix is highly variable,
and produces estimates with diverging condition numbers and over-dispersed
eigenvalues [@johnstone2001distribution]. Analyses employing this demonstrably
poor estimator may be negatively impacted.

As high-dimensional data have become widespread, researchers have derived many
novel covariance matrix estimators to remediate the sample covariance matrix's
shortcomings. These estimators come in many flavors, though most are
constructed by regularizing the sample covariance matrix. Comprehensive
reviews of these estimators are provided by @fan2016 and @pourahmadi2013.

This variety brings with it many challenges. Identifying an "optimal" estimator
from among a collection of candidates can prove a daunting task, one whose
objectivity is often compromised by the data analyst's decisions. Though
data-driven approaches for selecting an optimal estimator from among estimators
belonging to certain (limited) classes have been derived, the question of
selecting from a diverse collection of candidate procedures remains unaddressed.

# `cvCovEst` Framework

The solution provided by `cvCovEst` is a general, cross-validation-based,
estimator-agnostic framework for covariance matrix estimator selection. The
asymptotic optimality of selections are guaranteed under a few non-restrictive
assumptions by extending the seminal work of @laan_dudoit:2003, @dudoit2005, and
@vaart2006 on data-adaptive estimator selection to high-dimensional covariance
matrix estimation [@boileau2021]. Here, optimality is defined as choosing an
estimator with an equivalent risk difference to that which would have been
selected were the underlying data-generating distribution _completely known_.

The `cvCovEst` software package implements this framework for the `R` language
and environment for statistical computing [@R]. Included is a collection of
covariance matrix estimators spanning the work of many researchers (Table 1).
They may be employed independently of the cross-validation procedure.
`cvCovEst` also provides a variety of plotting and summary functions. These
diagnostic tools allow users to gauge the algorithm's performance, diagnose
issues that might arise during estimation procedures, and build intuition about
the many estimators' behaviors. Additionally, users have options to increase
the cross-validation method's computational efficiency via parallel
computation. Parallelization relies on the suite of `future` packages
[@future] by way of the `origami` package [@origami].

Table 1: Covariance matrix estimators implemented as of [version
0.3.5](https://github.com/PhilBoileau/cvCovEst).

|Estimator | Implementation | Description |
|----------|----------|-------------|
| Sample covariance matrix | `sampleCovEst()` | The sample covariance matrix. |
| Hard thresholding [@Bickel2008_thresh] | `thresholdingEst()` | Applies a hard thresholding operator to the entries of the sample covariance matrix. |
| SCAD thresholding [@rothman2009;@fan2001] | `scadEst()` | Applies the SCAD thresholding operator to the entries of the sample covariance matrix.|
| Adaptive LASSO [@rothman2009] | `adaptiveLassoEst()` | Applies the adaptive LASSO thresholding operator to the entries of the sample covariance matrix. |
| Banding [@bickel2008_banding] | `bandingEst()` | Replaces the sample covariance matrix's off-diagonal bands by zeros. |
| Tapering [@cai2010] | `taperingEst()` | Tapers the sample covariance matrix's off-diagonal bands, eventually replacing them by zeros. |
| Optimal Linear Shrinkage [@Ledoit2004] | `linearShrinkLWEst()` | Asymptotically optimal shrinkage of the sample covariance matrix towards the identity. |
| Linear Shrinkage [@Ledoit2004] | `linearShrinkEst()` | Shrinkage of the sample covariance matrix towards the identity, but the shrinkage is controlled by a hyperparameter. |
| Dense Linear Shrinkage [@shafer2005] | `denseLinearShrinkEst()` | Asymptotically optimal shrinkage of the sample covariance matrix towards a dense matrix whose diagonal elements are the mean of the sample covariance matrix's diagonal and whose off-diagonal elements are the mean of the sample covariance matrix's off-diagonal elements. |
| Nonlinear Shrinkage [@Ledoit2020] | `nlShrinkLWEst()` | Analytical estimator for the nonlinear shrinkage of the sample covariance matrix. |
| POET [@fan2013] | `poetEst()` | An estimator based on latent variable estimation and thresholding. |
| Robust POET [@fan2018] | `robustPoetEst()` | A robust (and more computationally taxing) take on the POET estimator. |


# Examples

We briefly showcase `cvCovEst`'s functionality through a toy example and an
application to single-cell transcriptomic data.

## Toy Dataset Example

Multivariate normal data are simulated using a covariance matrix with a Toeplitz
structure and then fed to the `cvCovEst` function. A summary of the
cross-validated estimation procedure is provided via the `plot` method.

```
library(MASS)
library(cvCovEst)
set.seed(1584)

# function to generate a toeplitz matrix
toep_sim <- function(p, rho, alpha) {

    times <- seq_len(p)
    H <- abs(outer(times, times, "-")) + diag(p)
    H <- H^-(1 + alpha) * rho
    covmat <- H + diag(p) * (1 - rho)

    sign_mat <- sapply(
      times,
      function(i) {
        sapply(
          times,
          function(j) {
            (-1)^(abs(i - j))
          }
        )
      }
    )
    return(covmat * sign_mat)
}

# generate a 100 x 100 covariance matrix
sim_covmat <- toep_sim(p = 100, rho = 0.6, alpha = 0.3)

# sample 75 observations from multivariate normal mean = 0, var = sim_covmat
sim_dat <-  MASS::mvrnorm(n = 75, mu = rep(0, 100), Sigma = sim_covmat)

# run CV-selector
cv_cov_est_sim <- cvCovEst(
  dat = sim_dat,
  estimators = c(
    linearShrinkEst, thresholdingEst, bandingEst, adaptiveLassoEst,
    sampleCovEst, taperingEst
  ),
  estimator_params = list(
    linearShrinkEst = list(alpha = seq(0.25, 0.75, 0.05)),
    thresholdingEst = list(gamma = seq(0.25, 0.75, 0.05)),
    bandingEst = list(k = seq(2L, 10L, 2L)),
    adaptiveLassoEst = list(lambda = c(0.1, 0.25, 0.5, 0.75, 1), n = seq(1, 5)),
    taperingEst = list(k = seq(2L, 10L, 2L))
  ),
  cv_scheme = "v_fold",
  v_folds = 5
)

# plot a summary of the results
plot(cv_cov_est_sim, data_in = sim_dat)
```

![A summary of the `cvCovEst` procedure's results. In the top left corner, the
selected estmators risk is plotted against its considered hyperparameters. In
the rop right, the eigenvalues of the selected estimator's estimate are
displayed. The bottom left plot presents the estimated covariance matrix.
Entries are colored based on their absolute values. Finally, the table in the
bottom right summarizes the performance of the best estimators from each class.
](summary_plot.png){ width=95% }


## Single Cell Transcriptomic Data

Single-cell transcriptome sequencing (scRNA-seq) measures the gene expression
profiles of individual cells within a given population, permitting the
identification of rare cell types and the study of developmental trajectories.
The datasets resulting from these experiments are typically high-dimensional:
expression data for hundreds or thousands of cells are collected for tens of
thousands of genes. A critical step in most analytic workflows is therefore
that of dimension reduction. In addition to facilitating visualization, this
reduction is thought to have a denoising effect. That is, the effects of
uninteresting biological variation are typically mitigated in these
lower-dimensional embeddings.

A popular method for the dimensionality reduction of scRNA-seq is uniform
manifold approximation and projection (UMAP) [@mcinnes2018], capable of
capturing non-linear relationships between features, applied to the dataset's
leading principal components. Since these principal components (PCs) are
derived from the sample covariance matrix, however, they are likely to be poor
estimates of the true PCs when the number of genes exceeds the number of cells.
Instead, the `cvCovEst` estimate could be used to compute the initial
dimensionality reduction.

Indeed, we find that the two-dimensional UMAP embedding resulting from the
`cvCovEst`-based approach improves upon that of the standard PCA-based approach
when applied to a dataset of 285 mouse visual coretex's cells' 1,000 most
variable genes [@tasic2016]. Fewer rare cells are misclustered, engendering a
47\% improvement in average silhouette width. For further discussion, see
@boileau2021.

![A comparison of UMAP embeddings using the 20 leading PCs from traditional PCA
and from `cvCovEst`-based PCA as initializations.](allen-umap.png){ width=95% }

# Availability

A stable release of the `cvCovEst` package is freely-available via the
[Comprehensive `R` Archive Network](https://CRAN.R-project.org/package=cvCovEst).
Its development version can be found on
[GitHub](https://github.com/PhilBoileau/cvCovEst). Documentation and examples
are contained in each version's manual pages, vignette, and `pkgdown` website
at https://philboileau.github.io/cvCovEst.

# Acknowledgments

Philippe Boileau's contribution to this work was supported by the Fonds de
recherche du Québec - Nature et technologies (B1X) and by the National
Institute of Environmental Health Sciences [P42ES004705] Superfund Research
Program at UC Berkeley.

We thank Jamarcus Liu for his contributions to the software package.

# References
