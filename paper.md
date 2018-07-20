---
title: 'logKDE: log-transformed kernel density estimation'
authors:
- affiliation: 1
  name: Andrew T. Jones
- affiliation: 2
  name: Hien D. Nguyen
  orcid: 0000-0002-9958-432X
- affiliation: 1
  name: Geoffrey J. McLachlan
date: "20 July 2018"
bibliography: paper.bib
tags:
- data visualization
- exploratory data analysis
- non-parametric
- positive data
- probability density function
- R
affiliations:
- index: 1
  name: School of Mathematics and Physics, University of Queensland, St. Lucia 4072,
    Queensland Australia
- index: 2
  name: Department of Mathematics and Statistics, La Trobe University, Bundoora 3086,
    Victoria Australia
---

# Summary

Exploratory data analysis, as proposed in @tukey1977exploratory, is an important paradigm for conducting meaningful and useful applied statistics. According to Chapters 5, 19, and 20 of @tukey1977exploratory, the visualization of data, and their density and distributional characterizations may provide practitioners with great insight into the necessary processes that are required in order to analyse the various features of the data.

Since its introduction in the pioneering works of @rosenblatt1956remarks and @parzen1962estimation, the method of kernel density estimation (KDE) has become among one of the most popular methods for estimation, interpolation, and visualization of probability density functions (PDFs), due to its effectiveness and simplicity. The manuscripts of @silverman2018density and @wand1994kernel provide thorough treatments on the topic of KDE.

In the base `stats` package of the R programming environment [@R2016], univariate KDE can be conducted via the `density` function. The KDE method implemented in the `density` function is the standard KDE technique for estimation of univariate PDFs, assuming that the data are real-valued. The `density` function is often taught in many introductory R classes and presented in many popular textbooks. See for example Section 5.6 of @venables1228bd, Section 6.3 of @keen2010graphics, Section 2.1 of @maindonald2010data, Section 2.3 @verzani2014using.

# References
