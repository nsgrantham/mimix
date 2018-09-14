library(argparse)
library(tidyverse)

parser <- ArgumentParser()
parser$add_argument("input", help = "Directory in which to read output files from aggregate-results.jl")
parser$add_argument("output", help = "Directory in which to write output files")
args <- parser$parse_args()

stopifnot(dir.exists(args$output))

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
