# AHR

Methods for estimation of multivariate average hazard ratios as
defined by Kalbfleisch and Prentice. The underlying survival functions
of the event of interest in each group can be estimated using either
the (weighted) Kaplan-Meier estimator or the Aalen-Johansen estimator
for the transition probabilities in Markov multi-state
models. Right-censored and left-truncated data is supported. Moreover,
the difference in restricted mean survival can be estimated.
    
## Installation

Get the released version from CRAN:

```R
install.packages("AHR")
```

Or from github:

```R
# install.packages("devtools")
devtools::install_github("mbrueckner/AHR")