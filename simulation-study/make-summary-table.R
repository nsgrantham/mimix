library(argparse)
library(tidyverse)
library(xtable)

parser <- ArgumentParser()
parser$add_argument("input", help = "Directory in which to read output files from aggregate-results.jl")
parser$add_argument("output", help = "Directory in which to write output files")
args <- parser$parse_args()

stopifnot(dir.exists(args$output))

# Table 1: Local test and estimation results

local_estimates_results_path <- file.path(args$input, "local-estimates-results.tsv")
stopifnot(file.exists(local_estimates_results_path))

local_results <- read_tsv(local_estimates_results_path) %>%
  left_join(global_results, on=c("model", "setting", "rep"))

df_local <- local_results %>% 
  group_by(setting, model, dense, block_var, error_var) %>%
  summarize(root_mean_squared_error = 100 * mean((value - mean)^2), 
            credible_interval_coverage_95 = 100 * mean((value > `2.5%`) & (value < `97.5%`))) %>%
  ungroup() %>%
  select(dense, error_var, block_var, model, root_mean_squared_error, credible_interval_coverage_95) %>%
  arrange(dense, error_var, block_var)

med_block_var <- df_local %>% 
  filter(block_var == 1) %>% 
  ungroup() %>%
  select(-block_var)

high_block_var <- df_local %>% 
  filter(block_var == 4) %>% 
  ungroup() %>%
  select(-dense, -error_var, -block_var)

local_estimates_table_path <- file.path(args$output, "local-estimates-table.txt")
table <- xtable(cbind(med_block_var, high_block_var), digits = 1)
print(table, file = local_estimates_table_path, include.rownames = FALSE)
