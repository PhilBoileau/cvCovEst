---
title: `cvCovEst`: Cross-validated covariance matrix estimation in `R`
tags:
  - R
  - covariance
  - cross-validation
authors:
  - name: Philippe Boileau
    orcid: 0000-0002-4850-2507
    affiliation: "1, 2"
  - name: Brian Collica
    orcid: 0000-0003-1127-2557
    affiliation: 3
  - name: Nima S. Hejazi
    orcid: 0000-0002-7127-2789
    affiliation: "1, 2"
  - name: Mark J. van der Laan
    orcid: 0000-0003-1432-5511
    affiliation: "4, 3, 2"
  - name: Sandrine Dudoit
    orcid: 0000-0002-6069-8629
    affiliation: "3, 4, 2"
affiliations:
  - name: Graduate Group in Biostatistics, University of California, Berkeley
    index: 1
  - name: Center for Computational Biology, University of California, Berkeley
    index: 2
  - name: Department of Statistics, University of California, Berkeley
    index: 3
  - name: Division of Biostatistics, School of Public Health, University of California, Berkeley
    index: 4
date: 03 February 2021
bibliography: paper.bib
---

# Summary

# Statement of Need

When the number of observations in a dataset far exceeds the number of
features, the estimator of choice for the covariance matrix is the sample
covariance matrix. It is efficient under minimal regularity assumptions on the
data-generating distribution. In high-dimensional regimes, however, its
performance is unsatisfactory: The sample covariance matrix is highly variable,
and produces estimates with diverging condition numbers and over-dispersed
eigenvalues. These issues amplify estimation error in resulting estimates. 

As high-dimensional data have become widespread, researchers have derived many
novel covariance matrix estimators to remediate the sample covariance matrix's
deficiencies. These estimators come in many flavours, though most are
constructed by regularizing the sample covariance matrix.

This variety brings with it many challenges. Identifying an "optimal" estimator
from among a collection of candidates can prove a daunting task, one whose
objectivity is often compromised by the analyst's decisions. Though data-driven
approaches for selecting an optimal estimator from among estimators belonging to
certain (limited) classes have been derived, the question of selecting from a
diverse collection of candidates remains unaddressed.

# `cvCovEst` Framework

Our response is a general, cross-validation-based framework for covariance
matrix estimator selection capable of accomplishing just that. The asymptotic
optimality of selections are guaranteed by extending the seminal work of
@laan_dudoit:2003, @dudoit2005, and @vaart2006 on data-adaptive estimator
selection to high-dimensional covariance matrix estimation [@boileau2021].

The `cvCovEst` software package implements this framework for the `R` language
and environment for statistical computing [@R]. Included is an accumulation of
covariance matrix estimators spanning the work of many researchers (Table 1)
They may be employed independently of the cross-validation procedure.
`cvCovEst` also provides a slew of plotting and summary functions. These
diagnostic tools allow users to gauge the algorithm's performance, diagnose
issues that might arise during estimation procedures, and build intuition about
the many estimators' behaviours. Given that the package was designed with
high-dimensional datasets in mind, users have options to increase the
cross-validation method's computational efficiency via parallel computation.
Parallelization relies on the suite of `future` packages [@future] by way of
the `origami` package [@origami].

Table 1: Covariance matrix estimators implemented as of version 0.3.4.

|Estimator | Implementation | Description |
|----------|----------|-------------|
| Sample covariance matrix | `sampleCovEst()` | The sample covariance matrix. |
| Hard thresholding [@Bickel2008_thresh] | `thresholdingEst()` | Applies a hard thresholding operator to the entries of the sample covariance matrix |
| SCAD thresholding [@rothman2009;@fan2001] | `scadEst()` | Applies the SCAD thresholding operator to the entries of the sample covariance matrix.|
| Adaptive LASSO [@rothman2009] | `adaptiveLassoEst()` | Applies the adaptive LASSO thresholding operator to the entries of the sample covariance matrix |
| Banding [@bickel2008_banding] | `bandingEst()` | Replaces the sample covariance matrix's off-diagonal bands by zeros. |
| Tapering [@cai2010] | `taperingEst()` | Tapers the sample covariance matrix's off-diagonal bands, eventually replacing them by zeros. |
| Optimal Linear Shrinkage [@Ledoit2004] | `linearShrinkLWEst()` | Asymptotically optimal shrinkage of the sample covariance matrix towards the identity. |
| Linear Shrinkage [@Ledoit2004] | `linearShrinkEst()` | Shrinkage of the sample covariance matrix towards the identity, but the shrinkage is controlled by a hyperparameter. |
| Dense Linear Shrinkage [@shafer2005] | `denseLinearShrinkEst()` | Asymptotically optimal shrinkage of the sample covariance matrix towards a dense matrix whose diagonal elements are the mean of the sample covariance matrix's diagonal, and whose off-diagonal elements are the mean of the sample covariance matrix's off-diagonal elements. |
| Nonlinear Shrinkage [@Ledoit2020] | `nlShrinkLWEst()` | Analytical estimator for the nonlinear shrinkage of the sample covariance matrix. |
| POET [@fan2013] | `poetEst()` | An estimator based on latent variable estimation and thresholding. |
| Robust POET [@fan2018] | `robustPoetEst()` | A robust (and more computationally taxing) take on the POET estimator. |


# Toy Dataset Example

We showcase `cvCovEst`'s functionality through a toy example.


# Availability

A stable release of the `cvCovEst` package is freely available via the
[Comprehensive R Archive Network](https://CRAN.R-project.org/package=cvCovEst);
its development version can be found on
[GitHub](https://github.com/PhilBoileau/cvCovEst). Documentation and examples
are contained in each version's manual pages and vignettes.

# Acknowledgments

Philippe Boileau's contribution to this work was supported by the Fonds de
recherche du Québec - Nature et technologies (B1X) by the National Institute of
Environmental Health Sciences [P42ES004705] Superfund Research Program at UC
Berkeley.

We thank Jamarcus Liu for his contributions to the software package.

# References
