---
title: "`cvCovEst`: Cross-validated covariance matrix estimation in `R`"
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
  - name: Division of Epidemiology and Biostatistics, School of Public Health, University of California, Berkeley
    index: 4
date: 03 February 2021
bibliography: paper.bib
---

# Summary

When the number of observations in a dataset far exceeds the number of
features, the estimator of choice for the covariance matrix is the sample
covariance matrix. It is an efficient estimator under minimal regularity
assumptions on the data-generating distribution. In high-dimensional regimes,
however, it leaves much to be desired: The sample covariance matrix is either
singular, numerically unstable, or both, thereby amplifying estimation error.

As high-dimensional data have become widespread, researchers have derived many
novel covariance matrix estimators to remediate the sample covariance matrix's
deficiencies. These estimators come in many flavours, though most are
constructed by regularizing the sample covariance matrix, or through the
estimation of latent factors. A comprehensive review is provided by @fan2016.

This variety brings with it many challenges. Identifying an "optimal" estimator
from among a collection of candidates can prove a daunting task, one whose
objectivity is often compromised by the analyst's decisions. Though data-driven
approaches for selecting an optimal estimator from among estimators belonging to
certain (limited) classes have been derived, the question of selecting an
estimator from among a diverse collection of candidates remains unaddressed.

Our response is a general, cross-validation-based framework for covariance
matrix estimator selection capable of accomplishing just that. The asymptotic
optimality of selections are guaranteed based upon extensions of the seminal
work of @laan_dudoit:2003, @dudoit2005, and @vaart2006 on data-adaptive
estimator selection to high-dimensional covariance matrix estimation
[@boileau2021]. The interested reader is invited to review theoretical
underpinnings of the methodology as described in @boileau2021.

The `cvCovEst` software package implements this framework for the `R` language
and environment for statistical computing [@R]. Given that the package was
specifically designed with high-dimensional datasets in mind, we provide users
with options to increase the method's computationally efficiency. In
particular, we make use of the suite of `future` software packages [@future] by
way of the `origami` software package [@origami] to facilitate parallel
computation. Finally, `cvCovEst` provides a slew of plotting and summary
functions. These tools allow users to gauge the algorithm's performance,
diagnose issues, and gain intuition about the behaviour of the many implemented
estimators.


# Acknowledgments

Philippe Boileau's contribution to this work was supported by the Fonds de
recherche du Québec - Nature et technologies (B1X) by the National Institute of
Environmental Health Sciences [P42ES004705] Superfund Research Program at UC
Berkeley.

We thank Jamarcus Liu for his contributions to the software package.

# References

