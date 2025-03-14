#############################
Instructions for Reproducing Application 5.3

###################################

The repository contains two distinct folders:

 - B-splines (Contains the codes for run the model of application 5.3 with B-spline basis.)
 - Bernstein (Contains the codes for run the model of application 5.3 with Bernstein basis.)

#############################################
#### 1. Files
#############################################

The files in both folders follow the same logic. First, we have the .R files:

   - `Fitting.R` (Fits the model and reproduces the information and figures from Application 5.3)

1) **Models folder** includes:
    - `model.R` (Contains the function that fit the model)
    - `funcs_aux.R` (Some auxiliary functions for model fitting)

2) **cpp folder** includes:
    - `arg_min.cpp` (The function defined in 3.36, page 48.)
    - `gradiente.cpp` (The function defined in 3.37, page 49.)
    - `funcs.cpp` (Functions for fitting model written in C++)
    - `cpo.cpp` (Functions to evaluate predictive measures written in C++)

##################################
#### 2. Installing Required Packages
###################################
- The following packages are essential for running the models, such as **Rcpp**, and **RcppArmadillo**.
Installation may vary depending on the operating system.

install.packages("Rcpp")
install.packages("RcppArmadillo")

***
Note: Before developing with Rcpp, you need to install a C++ compiler.
The book "Rcpp for everyone" helps you to install in your OS (https://teuder.github.io/rcpp4everyone_en/020_install.html)
***

- Also is important install other packages. We recommend to run the following sequence of codes:

install.packages("dplyr")
install.packages("tramME")
install.packages("lme4")
install.packages("splines2")
install.packages("Matrix")
install.packages("parallel")
install.packages("MASS")
install.packages("mgcv")
install.packages("scam")
install.packages("invgamma")
install.packages("pracma")
install.packages("fastmatrix")
install.packages("DescTools")
install.packages("plyr")
install.packages("LaplacesDemon")
install.packages("expm")
install.packages("fBasics")
install.packages("coda")
install.packages("interp")
install.packages("survival")
install.packages("tidyverse")
install.packages("ggplot2")

###############################
#### 3. Fitting the Model
###############################

To fit the model from 5.3 you just run the `Fitting.R` script. 
All necessary functions will be loaded from the script.
