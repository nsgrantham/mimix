library(ggplot2)
library(dplyr)
library(readr)
library(reshape2)
library(forcats)
library(stringr)
theme_set(theme_minimal())

# Figure 1 was created with Monodraw (https://monodraw.helftone.com)

## NutNet data analysis

Y <- read.csv(file.path("nutnet-analysis", "data", "Y.csv"), header=FALSE)
X <- read.csv(file.path("nutnet-analysis", "data", "X.csv"), header=FALSE)
Z <- read.csv(file.path("nutnet-analysis", "data", "Z.csv"), header=FALSE)
K <- ncol(Y)

if (file.exists(file.path("nutnet-analysis", "data", "tax.csv"))) {
  tax <- read.csv(file.path("data", "tax.csv"), header=FALSE, stringsAsFactors = FALSE)
  tax$V2[tax$V2 == "None"] <- "; ; ; ; ; ; Fungi sp"
  tax_table <- as.data.frame(do.call(rbind, strsplit(tax$V2, "; ")), stringsAsFactors = FALSE)
  otu_names <- make.unique(sub(" sp", "", gsub("_", " ", sub("s__", "", tax_table$V7))), sep=" ")
} else {
  otu_names <- sapply(1:K, toString)
}

# Posterior predictive checks

plot_post_pred_check <- function(obs, preds, group, xlab="", ylab="") {
  post_pred_bounds <- t(apply(preds, 1, function(x) quantile(x, c(0.025, 0.975))))
  obs_in_bounds <- (post_pred_bounds[, 1] < obs) & (obs < post_pred_bounds[, 2])
  df <- cbind.data.frame(
    x = 1:table(group)[1],
    y = unname(obs),
    signif = ordered(unname(obs_in_bounds), levels=c(TRUE, FALSE)),
    group = group
  )
  limits <- aes(
    ymin = post_pred_bounds[, 1], 
    ymax = post_pred_bounds[, 2], 
    alpha=signif
  )
  ggplot(df, aes(x=x, y=y)) + 
    geom_errorbar(limits, width=0.1, size=0.4) + 
    scale_alpha_manual(values=c(0.3, 1)) +
    geom_point(size = 0.7, stroke = 0, shape = 16) +
    facet_grid(. ~ group) +
    labs(x = xlab, y = ylab) +
    scale_x_continuous(breaks = scales::pretty_breaks(n = 10)) +
    scale_y_continuous(breaks = scales::pretty_breaks(n = 5)) +
    guides(alpha=FALSE)
}

obs_max_count <- read.delim(file.path("nutnet-analysis", "results", "mimix", "obs-max-count.tsv"), header=FALSE)
obs_prop_eq_zero <- read.delim(file.path("nutnet-analysis", "results", "mimix", "obs-prop-eq-zero.tsv"), header=FALSE)
obs_prop_leq_two <- read.delim(file.path("nutnet-analysis", "results", "mimix", "obs-prop-leq-two.tsv"), header=FALSE)
post_pred_max_count <- t(read.delim(file.path("nutnet-analysis", "results", "mimix", "post-pred-max-count.tsv"), header=FALSE))
post_pred_prop_eq_zero <- t(read.delim(file.path("nutnet-analysis", "results", "mimix", "post-pred-prop-eq-zero.tsv"), header=FALSE))
post_pred_prop_leq_two <- t(read.delim(file.path("nutnet-analysis", "results", "mimix", "post-pred-prop-leq-two.tsv"), header=FALSE))
#plot_post_pred_check(obs_max_count, post_pred_max_count)
#plot_post_pred_check(obs_prop_eq_zero, post_pred_prop_eq_zero, x="Sample", y="Proportion of OTUs in sample with no reads", title="Posterior predictive check of sparsity")
#ggsave("../figures/post-pred-prop-eq-zero.png", width=7.5, height=6)
#plot_post_pred_check(obs_prop_leq_one, post_pred_prop_leq_one)
obs_order <- order(obs_prop_eq_zero)
obs_prop_eq_zero <- obs_prop_eq_zero[obs_order, ]
post_pred_prop_eq_zero <- post_pred_prop_eq_zero[obs_order, ]
obs_order <- order(obs_prop_leq_two)
obs_prop_leq_two <- obs_prop_leq_two[obs_order, ]
post_pred_prop_leq_two <- post_pred_prop_leq_two[obs_order, ]
p <- plot_post_pred_check(c(obs_prop_eq_zero, obs_prop_leq_two), 
                     rbind(post_pred_prop_eq_zero, post_pred_prop_leq_two), 
                     group = c(rep("Exactly zero reads", length(obs_prop_eq_zero)), 
                               rep("Two or fewer reads", length(obs_prop_leq_two))),
                     xlab="Sample IDs", ylab="Proportion of OTUs in sample")#, 
ggsave(file.path("nutnet-analysis", "figures", "post-pred-sparsity.png"), plot = p, width=7.5, height=4)


## Global variable selection results

# MIMIX posterior probability of variable inclusion

omega <- read_tsv(file.path("nutnet-analysis", "results", "mimix", "omega.tsv"))
mean(rowSums(select(omega, ends_with("1]"))) > 0)  # herbivore exclusion
mean(rowSums(select(omega, ends_with("2]"))) > 0)  # nutrient amendment
mean(rowSums(select(omega, ends_with("3]"))) > 0)  # interaction

# MIMIX w/o Factors posterior probability of variable inclusion

omega <- read_tsv(file.path("nutnet-analysis", "results", "mimix-no-factors", "omega.tsv"))
mean(rowSums(select(omega, ends_with("1]"))) > 0)  # herbivore exclusion
mean(rowSums(select(omega, ends_with("2]"))) > 0)  # nutrient amendment
mean(rowSums(select(omega, ends_with("3]"))) > 0)  # interaction


## Get ordering of OTUs from hierarchical clustering of the factor correlation matrix

Lambda <- as.matrix(read.delim(file.path("nutnet-analysis", "results", "mimix", "Lambda-postmean.tsv"), header=FALSE))
LambdatLambda <- t(Lambda) %*% Lambda
R <- cov2cor(LambdatLambda)
C <- abs(R)
clusters <- hclust(as.dist(1-C))
ids <- rev(clusters$order)
R <- R[ids, ids]
LambdatLambda <- LambdatLambda[ids, ids]
otu_names <- otu_names[ids]

beta2 <- select(read_tsv(file.path("nutnet-analysis", "results", "mimix", "beta.tsv")), ends_with("2]"))
#beta2 <- read_tsv(file.path("results", "analyze", "mimix-beta2.csv"), col_names=FALSE)
beta2 <- beta2[, ids]      # order OTUs
alpha <- 0.05
quant_beta <- apply(beta2, 2, function(x) quantile(x, probs = c(alpha/2, 1-alpha/2)))
lower_beta <- quant_beta[1, ]
upper_beta <- quant_beta[2, ]
is_signif <- unname(!((lower_beta <= 0) & (0 <= upper_beta)))
mean_beta <- colMeans(beta2)

## Figure 4a: Posterior beta estimates for nutrient supplement

df_beta <- data.frame(
  otus = ordered(otu_names, levels=rev(otu_names)),
  beta = mean_beta,
  signif = factor(is_signif)
)

limits <- aes(ymax = upper_beta, ymin = lower_beta)
p <- ggplot(df_beta, aes(otus, beta, alpha = signif), size = 0.1) + 
  geom_errorbar(limits) +
  geom_hline(yintercept = 0, linetype = 2) +
  scale_alpha_manual(values = c(0.1, 1)) +
  guides(alpha = FALSE) +
  scale_x_discrete(breaks = NULL) + 
  scale_y_continuous(breaks = seq(-2, 4, by = 1), lim = c(-2.5, 4.5)) +
  labs(x = "OTUs", y = expression(paste("Posterior ", beta))) +
  theme(axis.ticks.x = element_blank()) +
  coord_flip()
ggsave(file.path("nutnet-analysis", "figures", "betarange.png"), width=3.75, height=9)


## Figure 4b: Posterior beta estimates for nutrient supplement (only significantly non-zero)

lower_beta_signif <- lower_beta[is_signif]
upper_beta_signif <- upper_beta[is_signif]
mean_beta_signif <- mean_beta[is_signif]

otu_names_signif <- otu_names[is_signif]

df_beta_signif <- data.frame(
  otus = ordered(otu_names_signif, levels = rev(otu_names_signif)),
  beta = mean_beta_signif
)

limits <- aes(ymax = upper_beta_signif, ymin = lower_beta_signif)
p <- ggplot(df_beta_signif, aes(otus, beta), size = 0.1) + 
  geom_point(size = 1) +
  geom_errorbar(limits, width = 0.1) +
  geom_hline(yintercept = 0, linetype = 2) +
  scale_y_continuous(breaks = seq(-2, 4, by = 1), lim = c(-2.5, 4.0)) +
  labs(x = "OTUs", y = expression(paste("Posterior ", beta))) +
  theme(axis.text.y = element_text(size=6)) +
  coord_flip()
ggsave(file.path("nutnet-analysis", "figures", "betarange-signif.png"), plot = p, width = 3.75, height = 9)


## Figure 5: Proportion of variance explained by site- and block-level effects

g_var <- read_tsv(file.path("nutnet-analysis", "results", "mimix", "g_var.tsv"))
theta_var <- read_tsv(file.path("nutnet-analysis", "results", "mimix", "theta_var.tsv"))
theta_var <- theta_var[, ids]

g_var_hat <- colMeans(g_var)
theta_var_hat <- colMeans(theta_var)
LLt_factor_var <- diag(LambdatLambda) * (1 + sum(g_var_hat))
eta <- LLt_factor_var / (LLt_factor_var + theta_var_hat)
mean(eta > 0.5)

df_var <- data.frame(
  x = theta_var_hat, 
  y = eta
)

df_var_signif <- data.frame(
  x = theta_var_hat[is_signif],
  y = eta[is_signif], 
  z = otu_names_signif
)

p <- ggplot(df_var) +
  geom_point(aes(x=x, y=y), size = 0.5, alpha = 0.2) +
  geom_point(data = df_var_signif, aes(x=x, y=y), size = 1) +
  geom_text(data = df_var_signif, aes(x=x, y=y, label=z), size = 2, 
            hjust = 0, vjust = 0, nudge_x = 0.2, nudge_y = 0.01, check_overlap = TRUE) +
  scale_x_continuous(lim = c(0, 35)) + 
  labs(x = expression(paste("Estimated error variance, ", hat(sigma)[k]^2)),
       y = expression(paste("Estimated prop. of variance due to site- & block-level effects, ", hat(eta)[k])))

ggsave(file.path("nutnet-analysis", "figures", "proportion-variance.png"), plot = p, width = 7.5, height = 5)


## Factor correlation matrix with significant OTUs marked on the diagonal
pal <- c('#ff8c00','#ffb663','#ffdca3','#ffffe0','#e2bddb','#c17bd5','#9932cc')
cuts <- c(-1, -0.7, -0.4, -0.1, 0.1, 0.4, 0.7, 1.0)

get_lower_tri <- function(cormat){
  cormat[upper.tri(cormat)] <- NA
  return(cormat)
}

R_melt <- melt(get_lower_tri(R), na.rm = TRUE)
R_melt$Var1 <- as.numeric(R_melt$Var1)
R_melt$Var2 <- as.numeric(R_melt$Var2)
R_melt$value2 <- cut(R_melt$value, breaks=cuts, include.lowest=TRUE)

idx_signif <- which(is_signif)
df_idx_signif <- data.frame(
  x = idx_signif,
  y = idx_signif,
  z = otu_names_signif
)

p <- ggplot(data = R_melt, aes(x = Var2, y = Var1, fill = value2)) + 
  geom_tile(size=0) +
  scale_fill_manual(values=pal, name="Correlation", drop=FALSE) +
  guides(fill=guide_legend(reverse=TRUE)) +
  scale_y_reverse() +
  geom_text(data = df_idx_signif, aes(x = x, y = y, label="—"), 
            size = 3, inherit.aes = FALSE, angle=45, nudge_x=12, nudge_y=18) +
  geom_text(data = df_idx_signif, aes(x = x, y = y, label=z), inherit.aes = FALSE, 
            size=2.5, hjust=0, vjust=0, nudge_x=30, nudge_y=35, check_overlap=TRUE) +
  labs(x = "OTUs", y = "OTUs") +
  theme(
    axis.text.y = element_blank(),
    axis.text.x = element_blank(),
    panel.grid = element_blank(),
    panel.background = element_blank(),
    axis.ticks = element_blank(),
    legend.position = c(0.9, 0.8)
  )
ggsave(file.path("nutnet-analysis", "figures", "correlation.png"), width = 10, height = 10)