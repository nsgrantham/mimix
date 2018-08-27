library(argparse)
library(tidyverse)
library(xtable)

global_results <- read_tsv("~/Desktop/sim-mimix/global-test-results.tsv")

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

p <- ggplot(df_global, aes(x=round(100*dense), y=prob_reject, color=model))
p <- p + geom_line(size=0.7)
p <- p + geom_hline(yintercept=0.05, size=0.2, linetype='longdash')
p <- p + geom_point(size=1)
p <- p + facet_grid(block_var ~ error_var, labeller = labeller(block_var = block_var_description, error_var = error_var_description))
p <- p + scale_color_manual(values=rev(c('#e2bddb','#c17bd5','#9932cc')))
p <- p + labs(x="Dense (%)", y="Estimated probability of rejection", color=NULL)
p <- p + theme(legend.position = "top", text = element_text(size = 14))
p <- p + theme_minimal()
p
ggsave(file.path("figures", "global-test.png"), width=7.5, height=5)

# Table 1: Local test and estimation results

local_results <- read_tsv("~/Desktop/sim-mimix/local-estimates-results.tsv") %>%
  left_join(global_results, on=c("model", "setting", "rep"))

df_local <- local_results %>% 
  group_by(setting, model, dense, block_var, error_var) %>%
  summarize(RMSE = 100 * mean((value - mean)^2), 
            C95 = 100 * mean((value > `2.5%`) & (value < `97.5%`))) %>%
  ungroup() %>%
  select(dense, error_var, block_var, model, RMSE, C95) %>%
  arrange(dense, error_var, block_var)

med_block_var <- df_local %>% 
  filter(block_var == 1) %>% 
  ungroup() %>%
  select(dense, error_var, model, RMSE, C95)

high_block_var <- df_local %>% 
  filter(block_var == 4) %>% 
  ungroup() %>%
  select(dense, error_var, model, RMSE, C95)

print(xtable(cbind(med_block_var, high_block_var), digits=1), include.rownames=FALSE)  

