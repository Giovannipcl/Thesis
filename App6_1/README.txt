#############################
Instructions for Reproducing Application 6.1

###################################

The repository contains two distinct folders:

 - macOS (If your machine runs macOS, we suggest using this folder. If you are using a Linux distribution, you can try it, but some Linux versions may require additional adjustments. In that case, we recommend using the Other folder.)
  - Other (For Windows and Linux distributions.)

#############################################
#### 1. Files
#############################################

The files in both folders follow the same logic. First, we have the .R files:

   - `Fitting.R` (Fits the model and reproduces the information and figures from Application 6.1)
   - `Model.R` (Contains the functions that fit the Model)

the following .cpp file

   - `funcs.cpp` (Functions for fitting the model written in C++)

And the following .csv

   - `DataFinal.csv` (Final dataset for the application)


##################################
#### 2. Installing Required Packages
###################################
- The following packages are essential for running the models, such as **Rcpp**, **RcppArmadillo**, and **RcppParallel**. 
Installation may vary depending on the operating system.

install.packages("Rcpp")
install.packages("RcppArmadillo")
install.packages("RcppParallel")

***
Note: Before developing with Rcpp, you need to install a C++ compiler.
The book "Rcpp for everyone" helps you to install in your OS (https://teuder.github.io/rcpp4everyone_en/020_install.html)
***

- Also is important install other packages. We recommend to run the following sequence of codes:

install.packages("dplyr")
install.packages("splines2")
install.packages("Matrix")
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


To fit the model from 6.2 you just run the `Fitting.R` script. 
All necessary functions will be loaded from the script.


