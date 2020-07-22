library(tidyverse)
library(fs)
# library(readr)
# library(stringr)
# library(glue)
# library(dplyr)
# library(purrr)
# library(ggplot2)


fit_pars = read_csv("fitness_parameters.txt", col_names = c("name","value"))


# Read output files except phenotypes ====
read_non_phenotype_output <- function(path2dir){
  setwd(path2dir)
  # list output files of interest
  all_output_files <- path_file(dir_ls(my_dir, glob = "*out*"))
  sub_output_files <- all_output_files[all_output_files != "out_phenotypes.txt"]
  
  # WARNING: below is unfinished business
  variable_file <- output_files[5]
  variable_name <- str_match(variable_file, "(?<=_).*(?=\\.)")
  variable_tibble <- read_csv(variable_file, col_names = c(paste(variable_name,"_mean", sep = ""), paste(variable_name, "_variance", sep = "")))
  
}

output_files <- path_file(dir_ls(my_dir, glob = "*out*"))

# Read phenotype data ====
read_phenotype_data <- function(path2file, method = 3){
  phenotypes <- read_csv(path2file, col_names = FALSE)
  phenotypes_range <- 1:(dim(phenotypes)[2]/2)
  phenotypes_as_strings <- lapply(phenotypes_range, toString)
  phenotypes_names <- "phenotype" %>% 
    paste(phenotypes_as_strings, sep = "")
  phenotypes_mean_colnames <- phenotypes_names %>% 
    paste("_mean", sep="")
  phenotypes_variance_colnames <- phenotypes_names %>% 
    paste("_variance", sep="")
  column_new_names <- c(phenotypes_mean_colnames, phenotypes_variance_colnames)
  
  if (method == 1) {
    ## My method with dplyr rename
    ## WARNING: Fails
    column_old_names <- colnames(phenotypes)
    command_renaming <- column_new_names %>% paste("=") %>% paste(column_old_names)
    phenotypes_renamed <- phenotypes %>% rename(!!!rlang:parse_exprs(command_renaming))
  } else if (method == 2) {
    ## Raph's method with base R
    ## WARNING: might lead to errors down the line
    colnames(phenotypes) <- column_new_names
    phenotypes_renamed <- phenotypes
  } else if (method == 3) {
    ## Raph's method with dplyr rename
    #1: create a vector of the current column names
    new_names <- colnames(phenotypes) 
    #2: name the vector elements with their soon-to-be new names
    names(new_names) <- column_new_names 
    #3: evaluate the named vector to assign column names to phenotypes
    phenotypes_renamed <- phenotypes %>% 
      rename(!!new_names) 
  }
  
  generations_range <- 0:(dim(phenotypes)[1]-1)
  final_df <- bind_cols(generations_range, phenotypes_renamed) %>% rename(generation = ...1)
  
  return(final_df)
  
}

# PLOT ====
plot_phenotype <- function(phen_dat){
  ggplot(data = phen_dat) +
  geom_ribbon(mapping = aes(x = generation, ymin = phenotype1_mean - phenotype1_variance, ymax = phenotype1_mean + phenotype1_variance),
              colour = "blue",
              stat = "smooth") +
  geom_point(mapping = aes(x = generation, y = phenotype1_mean), 
            size = 0.1
            ) +
  theme_classic()
}

plot_phenotype_smooth <- function(phendat){
  ggplot(data = phendat, aes(x = generation, y = phenotype1_mean)) +
  geom_smooth() +
  theme_classic()
}

# Example ====
my_dir = "/home/claire/Desktop/technofirstbatch/technology_alphaResources0.1betaTech0.1p0.1atech2.0btech1.0q0.1gamma0.01rb2.0"
setwd(my_dir)
pheno_data <- read_phenotype_data("out_phenotypes.txt")
plot_phenotype_smooth(pheno_data)
plot_phenotype(pheno_data)
