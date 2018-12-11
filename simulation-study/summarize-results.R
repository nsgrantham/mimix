library(argparse)
library(xtable)
library(tidyverse)

parser <- ArgumentParser()
parser$add_argument("input", help = "Directory in which to read output files from aggregate-results.jl")
parser$add_argument("output", help = "Directory in which to write output files")
args <- parser$parse_args()

stopifnot(dir.exists(args$output))

recode_models <- list(
  "permanova" = "PERMANOVA",
  "mimix-20-factors-G-priors" = "MIMIX Gaussian",
  "mimix-0-factors-DL-priors" = "MIMIX w/o Factors",
  "mimix-20-factors-DL-priors" = "MIMIX"
)

# Figure 2: Global hypothesis test results

global_test_results_path <- file.path(args$input, "global-test-results.tsv")
stopifnot(file.exists(global_test_results_path))

global_results <- read_tsv(global_test_results_path)

df_global <- global_results %>%
  mutate(model = fct_relevel(recode(model, !!! recode_models), rev(unlist(recode_models)))) %>%
  group_by(model, dense, block_var, error_var, form) %>%
  summarize(prob_reject = mean(as.logical(reject_global_null))) %>%
  ungroup()

block_var_description <- c(
  "1" = "Medium block variance",
  "4" = "High block variance"
)

error_var_description <- c(
  "1" = "Medium error variance",
  "4" = "High error variance",
  "9" = "Very high error variance"
)

forms <- df_global %>%
  select(form) %>%
  pull() %>%
  unique()

for (beta_form in forms) {
  df_form <- df_global %>%
    filter(form == beta_form)
  p <- ggplot(df_form, aes(round(100*dense), prob_reject, color = model)) + 
    geom_line(size = 0.7, alpha=0.7) +
    #geom_hline(yintercept = 0.05, size = 0.2, linetype = 'longdash') +
    geom_point(size = 1) +
    facet_grid(block_var ~ error_var, 
               labeller = labeller(block_var = block_var_description,
                                   error_var = error_var_description)) +
    scale_color_manual(values = rev(c('#e2bddb','#c17bd5','#9932cc', '#4b0082'))) +
    labs(x = "Dense (%)", y = "Proportion of global null rejections in 50 trials", color = NULL) +
    theme_minimal(base_size = 14) +
    theme(legend.position = "top")

  global_test_figure_path <- file.path(args$output, paste0("global-test-", beta_form, ".png"))
  ggsave(global_test_figure_path, plot = p, width = 7.5, height = 5)
}


# Table 1: Local test and estimation results

local_estimates_results_path <- file.path(args$input, "local-estimates-results.tsv")
stopifnot(file.exists(local_estimates_results_path))

recode_models[["permanova"]] <- NULL

local_results <- read_tsv(local_estimates_results_path) %>%
  left_join(global_results, on=c("model", "setting", "rep")) %>%
  mutate(model = fct_relevel(recode(model, !!! recode_models), rev(unlist(recode_models))))

df_local <- local_results %>% 
  mutate(real_pos = (value != 0), 
         pred_pos = (`2.5%` > 0) | (`97.5%` < 0)) %>%
  group_by(setting, model, dense, block_var, error_var, form) %>%
  summarize(root_mean_squared_error = 100 * mean((value - mean)^2), 
            credible_interval_coverage_95 = 100 * mean((value >= `2.5%`) & (value <= `97.5%`)),
            tp = sum(real_pos & pred_pos),
            fp = sum(!real_pos & pred_pos),
            fn = sum(real_pos & !pred_pos),
            tn = sum(!real_pos & !pred_pos),
            tpr = 100 * tp / (tp + fn),
            tnr = 100 * tn / (tn + fp)) %>%
  ungroup() 

for (beta_form in forms) {
  df_form <- df_local  %>%
    filter(form == beta_form) %>%
    select(dense, error_var, block_var, model, root_mean_squared_error, credible_interval_coverage_95, tpr, tnr) %>%
    arrange(dense, error_var, block_var)

  med_block_var <- df_form %>% 
    filter(block_var == 1) %>% 
    ungroup() %>%
    select(-block_var)

  high_block_var <- df_form %>% 
    filter(block_var == 4) %>% 
    ungroup() %>%
    select(-model, -dense, -error_var, -block_var)

  local_estimates_table_path <- file.path(args$output, paste0("local-estimates-table-", beta_form, ".txt"))
  table <- xtable(cbind(med_block_var, high_block_var), digits = 1)
  print(table, file = local_estimates_table_path, include.rownames = FALSE)
}

