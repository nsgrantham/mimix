library(argparse)
library(xtable)
library(tidyverse)

parser <- ArgumentParser()
parser$add_argument("input", help = "Directory in which to read output files from aggregate-results.jl")
parser$add_argument("output", help = "Directory in which to write output files")
args <- parser$parse_args()

stopifnot(dir.exists(args$output))

# Figure 2: Global hypothesis test results

global_test_results_path <- file.path(args$input, "global-test-results.tsv")
stopifnot(file.exists(global_test_results_path))

global_results <- read_tsv(global_test_results_path)

df_global <- global_results %>%
  group_by(model, dense, block_var, error_var) %>%
  summarize(prob_reject = mean(as.logical(reject_global_null)))

block_var_description <- c(
  "1" = "Medium block variance",
  "4" = "High block variance"
)

error_var_description <- c(
  "1" = "Medium error variance",
  "4" = "High error variance",
  "9" = "Very high error variance"
)

p <- ggplot(df_global, aes(round(100*dense), prob_reject, color = model)) + 
  geom_line(size = 0.7) +
  geom_hline(yintercept = 0.05, size = 0.2, linetype = 'longdash') +
  geom_point(size = 1) +
  facet_grid(block_var ~ error_var, 
             labeller = labeller(block_var = block_var_description,
                                 error_var = error_var_description)) +
  scale_color_manual(values = rev(c('#e2bddb','#c17bd5','#9932cc'))) +
  labs(x = "Dense (%)", y = "Estimated probability of rejection", color = NULL) +
  theme(legend.position = "top", text = element_text(size = 14)) +
  theme_minimal()

global_test_figure_path <- file.path(args$output, "global-test-figure.png")
ggsave(global_test_figure_path, plot = p, width = 7.5, height = 5)


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
