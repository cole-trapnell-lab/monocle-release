## ----include=FALSE-------------------------------------------------------
library(knitr)
opts_chunk$set(
#concordance=TRUE
)

## ----package_loads, include=FALSE, eval=TRUE-----------------------------
library(Biobase)
library(knitr)
library(reshape2)
library(ggplot2)

knitr::opts_chunk$set(autodep=TRUE, cache=FALSE, warning=FALSE, dev='png', dpi=600)
set.seed(0)

