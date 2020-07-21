#library(tidyverse)
library(readr)
library(fs)
library(stringr)
library(glue)
library(dplyr)

my_dir = "/home/claire/Desktop/technofirstbatch/technology_alphaResources0.1betaTech0.1p0.1atech2.0btech1.0q0.1gamma0.01rb2.0"
setwd(my_dir)
fit_pars = read_csv("/home/claire/Desktop/technofirstbatch/technology_alphaResources0.1betaTech0.1p0.1atech2.0btech1.0q0.1gamma0.01rb2.0/fitness_parameters.txt", col_names = c("name","value"))
output_files = path_file(dir_ls(my_dir,glob = "*out*"))

variable_file = output_files[5]
variable_name = str_match(variable_file,"(?<=_).*(?=\\.)")
variable_tibble = read_csv(variable_file, col_names = c(paste(variable_name,"_mean",sep=""),paste(variable_name,"_variance",sep="")))


phenotypes = read_csv("out_phenotypes.txt",col_names = FALSE)
rename(phenotypes, c(phenotype1_mean,phenotype1_average,phenotype2_mean,phenotype2_average,phenotype3_mean,phenotype3_average,phenotype4_mean,phenotype4_average) = c(X1,X2,X3,X4,X5,X6,X7,X8))
phenotypes_renamed = phenotypes %>% rename(phenotype1_mean = X1, phenotype2_mean = X2, phenotype3_mean = X3, phenotype4_mean = X4, phenotype1_variance = X5, phenotype2_variance = X6, phenotype3_variance = X7, phenotype4_variance = X8)