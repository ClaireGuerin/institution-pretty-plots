library(tidyverse)
library(fs)
library(patchwork)
library(GGally)
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

reorganize_phenotype_summary <- function(phen_dat){
  summary <- phen_data %>%
    pivot_longer(cols = colnames(.)[-1], names_to = "var_name") %>% 
    # organize by column names (except generation number)
    mutate(
      var_type = str_remove(var_name, "^.*_"),
      var_name = str_remove(var_name, "_.*$")
    ) %>% 
    # create new column "var_type" that indicates mean or variance, and remove the _mean or _variance indication in "name" column
    mutate(value = ifelse(var_type == "variance", sqrt(value), value)) %>% 
    # apply square root to value if var_type is variance
    mutate(var_type = if_else(var_type == "variance", "std", var_type)) %>%
    group_by(var_name, var_type) %>% 
    # group by both variable name and type, thus creating 8 groups (if there are 4 phenotypes accounted for)
    summarize(mean = mean(value)) %>% 
    # for each group, calculate the mean value
    ungroup() %>%
    pivot_wider(names_from = "var_type", values_from = "mean") 
  # re-organize data frame so as to have one row by phenotype, and average value for phen mean and variance in 2 separate columns
  
  return(summary)
}

### example
setwd("/home/claire/Desktop/technofirstbatch/technology_alphaResources0.1betaTech0.1p0.1atech2.0btech1.0q0.1gamma0.01rb2.0")
phen_data <- read_phenotype_data("out_phenotypes.txt", 3)
phen_data %>% head
reorganize_phenotype_summary(phen_data)
###

# 1: loop over folders
extract_global_mean_and_variance <- function(file_path, is_phen = FALSE){
  variable_name <- str_match(file_path,"(?<=/out_).*(?=\\.)")
  
  # mean # name
  # function(x) mean(x) # lambda simple
  # function(x) { mean(sqrt(x)) } # lambda complexe
  # ~ mean(.x) # formula
  
  if (is_phen) {
    phen_data <- read_phenotype_data(file_path, 3)
    data_summary <- reorganize_phenotype_summary(phen_data)
    
  } else {
    mean_column <- paste(variable_name, "_mean", sep = "")
    variance_column <- paste(variable_name, "_variance", sep = "")
    variable_tibble <- read_csv(file_path, col_names = c(mean_column, variance_column))
    m <- summarize(variable_tibble, mean = mean(!!rlang::parse_expr(mean_column)))
    v <- summarize(variable_tibble, std = mean(sqrt(!!rlang::parse_expr(variance_column))))
    data_summary <- bind_cols(var_name = variable_name[,1], m, v)
  }
  
  return(data_summary)
}

extract_global_mean_and_variance("out_phenotypes.txt", TRUE)
test_demography <- extract_global_mean_and_variance("out_demography.txt")
extract_global_mean_and_variance(output_files_list[1])


extract_meta_data <- function(dir_path = "."){
  files_list <- dir_ls(dir_path, type = "file")
  output_files_list <- files_list[str_detect(files_list,"(?<=/)out_(?!phenotypes)")]
  
  # Collect summaries of non-phenotypic data
  global_summary <- output_files_list %>%
    map_dfr(extract_global_mean_and_variance) %>%
    bind_rows(extract_global_mean_and_variance(files_list[str_detect(files_list,"(?<=/)out_phenotypes.txt")], TRUE))
  
  # Add phenotype summary
  #phenotype_output_file <- "out_phenotypes.txt"
  #phenotype_tibble <- extract_global_mean_and_variance(phenotype_output_file)
  
  # Join
  
  return(global_summary)
  
}

extract_meta_data(".") 
setwd("~")
my_dir = "/home/claire/Desktop/technofirstbatch/technology_alphaResources0.1betaTech0.1p0.1atech2.0btech1.0q0.1gamma0.01rb2.0"
extract_meta_data(my_dir)

# 2: extract important information:
## - mean value from generation x (user-defined) to end
## - fitness parameter values

gather_pars_and_output_in_single_row <- function(path_to_sim){
  output_data <- extract_meta_data(path_to_sim) %>%
    pivot_longer(cols = c(mean, std), names_to = "val_type") %>%
    mutate(var_name = paste(var_name, "_", val_type, sep = "")) %>%
    select(-val_type) %>% 
    pivot_wider(names_from = var_name, values_from = value)
  
  input_data <- read_csv(paste(path_to_sim,"fitness_parameters.txt", sep = "/"), col_names = c("par_name","par_val")) %>% 
    pivot_wider(names_from = par_name, values_from = par_val)
  
  return(bind_cols(input_data, output_data))
}

gather_pars_and_output_in_single_row(my_dir)

setwd("~")
my_dir = "/home/claire/Desktop/technofirstbatch/technology_alphaResources0.1betaTech0.1p0.1atech2.0btech1.0q0.1gamma0.01rb2.0"
gather_pars_and_output_in_single_row(my_dir)
extract_meta_data(my_dir)

# list all sim folders
collect_all_simulations_data <- function(path_to_all_dirs){
  dir_ls(path_to_all_dirs, type = "directory") %>%
    map_dfr(gather_pars_and_output_in_single_row)
}

#==== Contour plots ====
# contour plot of mean values for each variable, according to different parameter values (user-defined)

data_set <- collect_all_simulations_data("/home/claire/Desktop/technofirstbatch")

fitness_parameters <- colnames(data_set)[!str_detect(colnames(data_set),"_")]
par_comb <- tibble(par = fitness_parameters) %>% 
  expand(x = par, y = par) %>% 
  filter(x < y)

par_order <- par_comb %>% count(x)
slice(par_comb, par_order$x[which.max(par_order$n)])
first_row <- 
  par_comb[which(par_comb$x == par_order$x[which.max(par_order$n)]),]

named_vector_aes <- first_row[1,] %>% 
  bind_cols(z = "resources_mean")  %>%
  pivot_longer(cols = c(x,y,z)) %>%
  deframe()


ggplot(data_set, aes(x = !!rlang::parse_expr(named_vector_aes["x"]), y = !!rlang::parse_expr(named_vector_aes["y"]), z = !!rlang::parse_expr(named_vector_aes["z"]))) + 
  geom_contour_filled() +
  theme_classic()

ggplot(data_set, aes(x = atech, y = alphaResources, z = resources_mean)) + 
  geom_contour_filled() +
  theme_classic()

npars <- as.numeric(tally(par_comb))
1:floor(sqrt(npars)) %>% rev

for(i in 1:floor(sqrt(npars))) {
  print(c(i, as.numeric(npars) / i))
}

data_set %>%
  select(-matches("_")) %>%
  ggally_denstrip(mapping = aes(x = q, y = p, z = resources_mean))
  
ggpairs()
ggplot(data_set, aes(x = q, y = p, z = resources_mean)) + 
  geom_density_2d_filled(alpha = 0.5) +
  theme_classic()

ggplot(data_set, aes(x = q, y = p, z = resources_mean)) + 
  geom_contour_filled() +
  theme_classic()

#==== example from GGally ====
# Small function to display plots only if it's interactive
p_ <- GGally::print_if_interactive

data(tips, package = "reshape")
p_(ggally_facetdensity(tips, mapping = ggplot2::aes(x = total_bill, y = sex)))
p_(ggally_facetdensity(tips,mapping = ggplot2::aes_string(y = "total_bill", x = "sex", color = "sex")))
###

# Desired format:
# Each row is one simulation, each column a parameter value or mean / variance variable

file_path <- "out_demography.txt"
variable_file <- path_file(file_path)
variable_name <- str_match(variable_file,"(?<=_).*(?=\\.)")

# mean # name
# function(x) mean(x) # lambda simple
# function(x) { mean(sqrt(x)) } # lambda complexe
# ~ mean(.x) # formula

mean_column <- paste(variable_name, "_mean", sep = "")
variance_column <- paste(variable_name, "_variance", sep = "")
variable_tibble <- read_csv(variable_file, col_names = c(mean_column, variance_column))
m <- summarize(variable_tibble, mean = mean(!!rlang::parse_expr(mean_column)))
v <- summarize(variable_tibble, std = mean(sqrt(!!rlang::parse_expr(variance_column))))
data_summary <- bind_cols(var_name = variable_name[,1], m, v)



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

extract_mean("out_demography.txt")
extract_meta_data()