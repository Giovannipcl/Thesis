#############################
Instructions for Reproducing Application 5.2 

###################################

The repository contains two distinct folders:

  - macOS (If your machine runs macOS, we suggest using this folder. If you are using a Linux distribution, you can try it, but some Linux versions may require additional adjustments. In that case, we recommend using the Other folder.)
  - Other (For Windows and Linux distributions.)

#############################################
#### 1. Files
#############################################

The files in both folders follow the same logic. First, we have the .R files:

   - `Fittingmodel1.R` (Fits Model 1 and reproduces the information and figures of Model 1 from Application 5.1)
   - `Fittingmodel2.R` (Fits Model 2 and reproduces the information and figures of Model 2 from Application 5.1)
   - `Fittingmodel3.R` (Fits Model 3 and reproduces the information and figures of Model 3 from Application 5.1)
   - `Fittingmodel4.R` (Fits Model 4 and reproduces the information and figures of Model 4 from Application 5.1)


1) **models folder** includes:

   - `All_models.R` (Contains the functions that fit Models 1, 2, 3, and 4)
   - `funcs_aux.R` (Some auxiliary functions for model fitting)


2) **cpp folder** includes:
   - `model1.cpp` (Functions for fitting Model 1 written in C++)
   - `model2.cpp` (Functions for fitting Model 2 written in C++)
   - `model3.cpp` (Functions for fitting Model 3 written in C++)
   - `model4.cpp` (Functions for fitting Model 4 written in C++)



##################################
#### 2. Installing Required Packages
###################################
- The following packages are essential for running the models, such as **Rcpp**, **RcppArmadillo**, and **RcppParallel**. 
Installation may vary depending on the operating system.

install.packages("Rcpp")
install.packages("RcppArmadillo")
install.packages("RcppParallel")

***
Note: Before developing with Rcpp, you may need to install a C++ compiler.
The book "Rcpp for everyone" helps you to install in your OS (https://teuder.github.io/rcpp4everyone_en/020_install.html)
***

- Also is important install other packages. We recommend to run the following sequence of codes:

install.packages("dplyr")
install.packages("splines2")
install.packages("optimParallel")    #Just for MacOS
install.packages("parallel")         #Just for MacOS
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
install.packages("RcppParallel")
install.packages("tidyverse")
install.packages("ggplot2")



###############################
#### 3. Fitting the Model
###############################


To fit any model from Section 5.1, simply choose one of the scripts: `model1.R`, `model2.R`, `model3.R`, or `model4.R`. 
All necessary functions will be loaded from the script.


