library(tidyverse)
library(fs)
library(patchwork)
# library(readr)
# library(stringr)
# library(glue)
# library(dplyr)
# library(purrr)
# library(ggplot2)

# Read output files except phenotypes ====
make_tibble_from_file <- function(filename){
  var_name <- str_match(filename, "(?<=_).*(?=\\.)")
  var_tibble <- read_csv(filename, col_names = c(paste(var_name,"_mean", sep = ""), paste(var_name, "_variance", sep = "")))
}

read_non_phenotype_data <- function(path2dir){
  setwd(path2dir)
  # list output files of interest
  all_output_files <- path_file(dir_ls(my_dir, glob = "*out*"))
  sub_output_files <- all_output_files[all_output_files != "out_phenotypes.txt"]
  
  # Make full tibble
  my_tibble <- make_tibble_from_file(sub_output_files[1])
  indexed_tibble <- 1:dim(my_tibble)[1] %>% 
    bind_cols(my_tibble) %>%
    rename(generation = ...1)
  
  for(variable_file in sub_output_files[-1]){
    indexed_tibble <- indexed_tibble %>% 
      bind_cols(make_tibble_from_file(variable_file))
  }
  
  return(indexed_tibble)
}


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

# Plots ====
plot_non_phen <- function(dat, var_name){
  var_average <- paste(var_name, "_mean", sep = "")
  var_variance <- paste(var_name, "_variance", sep = "")
  p <- ggplot(dat) +
    geom_ribbon(mapping = aes(x = generation,
                              ymin = !!rlang::parse_expr(var_average) - !!rlang::parse_expr(var_variance),
                              ymax = !!rlang::parse_expr(var_average) + !!rlang::parse_expr(var_variance)),
                colour = "gray") +
    geom_point(mapping = aes(x = generation, y = !!rlang::parse_expr(var_average)),
               colour = "seagreen3") +
    xlab("Time (number of generations)") +
    ylab(str_to_sentence(var_name)) +
    theme_classic()
  
  return(p)
}

plot_phenotype <- function(phen_dat, phen_index){
  pheno_average <- paste("phenotype", phen_index, "_mean", sep = "")
  pheno_variance <- paste("phenotype", phen_index, "_variance", sep = "")
  
  p <- ggplot(data = phen_dat) +
    geom_ribbon(mapping = aes(x = generation, ymin = !!rlang::parse_expr(pheno_average) - !!rlang::parse_expr(pheno_variance), ymax = !!rlang::parse_expr(pheno_average) + !!rlang::parse_expr(pheno_variance)),
              colour = "gray"
              ) +
    geom_point(mapping = aes(x = generation, y = !!rlang::parse_expr(pheno_average)), 
            size = 0.1,
            colour = "indianred3"
            ) +
    xlab("Time (number of generations)") +
    ylab(paste("Phenotype", phen_index, sep = " ")) +
  theme_classic()
  
  return(p)
}

plot_phenotype_smooth <- function(phen_dat, phen_index){
  phen_average <- paste("phenotype", phen_index, "_mean", sep = "")
  
  p <- ggplot(data = phen_dat, 
              aes(x = generation, y = !!rlang::parse_expr(phen_average))) +
  geom_smooth(colour = "coral1") +
  theme_classic()
  
  return(p)
}

# Read data from multiple simulations ==== 
# ... and store some kind of summary in a tibble
 
# 1: loop over folders
extract_mean <- function(file_path){
  variable_file <- path_file(file_path)
  variable_name <- str_match(variable_file,"(?<=_).*(?=\\.)")
  mean_column <- paste(variable_name, "_mean", sep = "")
  variance_column <- paste(variable_name, "_variance", sep = "")
  variable_tibble <- read_csv(variable_file, col_names = c(mean_column, variance_column))
  summarize(variable_tibble, mean = mean(!!rlang::parse_expr("demography_mean")))
}
extract_meta_data <- function(dir_path){
  directory_list <- dir_ls(dir_path) %>% file_info()
  file_path <- directory_list$path
  file_type <- directory_list$type
  
  
}
directory_list <- dir_ls("/home/claire/Desktop/test-folder") %>% file_info()
# 2: extract important information:
## - mean value from generation x (user-defined) to end
## - fitness parameter values
fit_pars = read_csv("fitness_parameters.txt", col_names = c("name","value"))
# 3: contour plot of mean values for each variable, according to different parameter values (user-defined)

# Example: single simulation ====
my_dir = "/home/claire/Desktop/technofirstbatch/technology_alphaResources0.1betaTech0.1p0.1atech2.0btech1.0q0.1gamma0.01rb2.0"
# setwd(my_dir)
non_pheno_data <- read_non_phenotype_data(my_dir)
pheno_data <- read_phenotype_data("out_phenotypes.txt")

plot_phenotype_smooth(pheno_data, 1)
plot_phenotype(pheno_data, 1)

## Plot phenotypes
1:4 %>% 
  map(plot_phenotype, phen_dat = pheno_data) %>%
  wrap_plots()

1:4 %>% 
  map(plot_phenotype_smooth, phen_dat = pheno_data) %>%
  wrap_plots()

## Plot other variables
c("consensus", "resources", "demography", "technology") %>% 
  map(plot_non_phen, dat = non_pheno_data) %>%
  wrap_plots()

# Example: multiple simulations ====
