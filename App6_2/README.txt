#############################
Instructions for Reproducing Application 6.2

#############################################
#### 1. Files
#############################################

The files in both folders follow the same logic. First, we have the files:

   - `fitting.R` (Fits Model 1 and reproduces the information and figures of Model 1 from Application 6.2)

   - `fitting_model2.R` (Fits Model 2 and reproduces the information and figures of Model 1 from Application 6.2)

   - `fitting_model3.R` (Fits Model 3 and reproduces the information and figures of Model 1 from Application 6.2)

    - `Allmodels.R` (Contains the functions that fit Models 1, 2, and 3)

    - `func_aux.R` (Some auxiliary functions for model fitting)

    - `funcs1.cpp` (Functions for fitting the model 1 written in C++)

    - `funcs2_3.cpp` (Functions for fitting the model 2 and 3 written in C++)

    - `arg_min.cpp` (The function defined in 3.36, page 48.)
   
    - `gradiente.cpp` (The function defined in 3.37, page 49.)
   
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
install.packages("coda")
install.packages("survival")
install.packages("tidyverse")
install.packages("ggplot2")

###############################
#### 3. Fitting the Model
###############################

To fit the models from 6.2 you just run the `fitting.R ` script to fit the model 1, and `fitting_model2.R` or `fitting_model3.R` to fit the model 2 or 3.
All necessary functions will be loaded from the script.
