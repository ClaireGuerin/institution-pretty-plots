#library(tidyverse)
library(readr)
library(fs)
library(stringr)
library(glue)

my_dir = "/home/claire/Desktop/technofirstbatch/technology_alphaResources0.1betaTech0.1p0.1atech2.0btech1.0q0.1gamma0.01rb2.0"
setwd(my_dir)
fit_pars = read_csv("/home/claire/Desktop/technofirstbatch/technology_alphaResources0.1betaTech0.1p0.1atech2.0btech1.0q0.1gamma0.01rb2.0/fitness_parameters.txt", col_names = c("name","value"))
output_files = path_file(dir_ls(my_dir,glob = "*out*"))

read_csv(output_files[5], col_names = c("name","value"))
