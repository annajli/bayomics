#!/usr/bin/env Rscript
# =============================================================================
# 03_bn_leaf_queries.R
# -----------------------------------------------------------------------------
# HELD-OUT VALIDATION — Step 3: attach the held-out outcome as a leaf and query
#
# The flu response was held OUT of structure learning, so it is not a node in
# any fitted network. To produce the literal conditional probabilities the
# validation calls for — e.g. P(seroconversion | age = Older, CMV = Positive) —
# we take the L5 network's LEARNED structure as fixed, add the held-out outcome
# as a new leaf node whose parents are the discrete roots, fit ONLY that node's
# CPD on the 88 joined subjects, and query it. The structure the roots' effects
# were learned from (root -> immune-cell edges) is untouched; we are only asking
# the fitted network to report the conditional distribution of the outcome.
#
# Design notes
#   * Binary seroconversion attaches as a discrete node. In a conditional
#     linear Gaussian (CLG) network a discrete node may only have discrete
#     parents, so its parents are the three roots (age, sex, cmv). The fitted
#     CPD is then exactly P(seroconvert | age, sex, cmv) — an exact, smoothed
#     version of the cohort's stratified rate, which is precisely the quantity
#     the validation asks the network to produce.
#   * Continuous targets (igg_log2fc_mean, hai_peak_max) attach as CLG nodes
#     with the roots as parents; we report their conditional means.
#   * We query the two literature-motivated contrasts:
#       "suppressed" profile : age = Older,  cmv = Positive
#       "responsive" profile : age = Young,  cmv = Negative
#     and each single-root contrast, marginalizing the unqueried roots.
#
# Outputs (output/validation/):
#   bn_leaf_cpqueries.csv   conditional probabilities / means per profile
# =============================================================================

suppressMessages({ library(bnlearn); library(dplyr); library(readr) })
set.seed(46)

here_bn <- if (file.exists("config")) "." else ".."
out_dir <- file.path(here_bn, "output", "validation")
m <- readRDS(file.path(here_bn, "output", "models",
                       "L5_freq_immunophenotype_models.rds"))
d <- readRDS(file.path(out_dir, "L5_validation_joined.rds"))

roots <- c("age_group", "sex", "cmv")
d$age_group <- factor(d$age_group, levels = c("Young Adult", "Older Adult"))
d$sex       <- factor(d$sex,       levels = c("Female", "Male"))
d$cmv       <- factor(d$cmv,       levels = c("Negative", "Positive"))

# Learned DAG structure (tabu network) — the object we hold fixed.
dag <- m$bn_tabu
feat_nodes <- setdiff(nodes(dag), roots)

# The fitted-model feature columns (needed to fit the augmented CLG network).
model_feat <- d[, feat_nodes, drop = FALSE]

# ---- Helper: build augmented network with one outcome leaf, fit, and query ----
# outcome_name : column in d to attach
# is_binary    : TRUE -> discrete factor node (parents = roots only)
#                FALSE -> continuous CLG node (parents = roots)
query_leaf <- function(outcome_name, is_binary) {
  # Assemble the training frame: roots + all learned features + the outcome.
  train <- cbind(d[, roots, drop = FALSE], model_feat)
  if (is_binary) {
    train[[outcome_name]] <- factor(ifelse(d[[outcome_name]], "yes", "no"),
                                    levels = c("no", "yes"))
  } else {
    train[[outcome_name]] <- as.numeric(d[[outcome_name]])
  }
  train <- as.data.frame(train)

  # Augment the DAG: add the outcome node + arcs from each root.
  aug <- dag
  aug <- add.node(aug, outcome_name)
  for (r in roots) aug <- set.arc(aug, from = r, to = outcome_name)

  # Fit CLG parameters on the augmented structure (existing nodes re-fit
  # identically; only the new leaf's CPD is genuinely new).
  fit <- bn.fit(aug, train)

  # Query helper via likelihood weighting given root evidence (a named list).
  # cpquery evaluates `event` with non-standard evaluation inside the sampler,
  # so we splice the actual node name into the call with substitute().
  q <- function(ev) {
    if (is_binary) {
      call <- substitute(
        cpquery(fit, event = (OUT == "yes"), evidence = ev,
                method = "lw", n = 2e5),
        list(OUT = as.name(outcome_name)))
      eval(call)
    } else {
      sims <- cpdist(fit, nodes = outcome_name, evidence = ev,
                     method = "lw", n = 2e5)
      w <- attr(sims, "weights")
      weighted.mean(sims[[outcome_name]], w)
    }
  }

  profiles <- list(
    `Older+CMV+ (suppressed)` = list(age_group = "Older Adult", cmv = "Positive"),
    `Young+CMV- (responsive)` = list(age_group = "Young Adult", cmv = "Negative"),
    `Older (any CMV)`         = list(age_group = "Older Adult"),
    `Young (any CMV)`         = list(age_group = "Young Adult"),
    `CMV+ (any age)`          = list(cmv = "Positive"),
    `CMV- (any age)`          = list(cmv = "Negative"),
    `Male (any age/CMV)`      = list(sex = "Male"),
    `Female (any age/CMV)`    = list(sex = "Female")
  )
  # Build the evidence list expression for likelihood weighting.
  ev_of <- function(p) {
    ev <- as.list(rep(NA, 0))
    for (nm in names(p)) ev[[nm]] <- p[[nm]]
    ev
  }
  res <- lapply(names(profiles), function(pn) {
    val <- q(ev_of(profiles[[pn]]))
    data.frame(target = outcome_name,
               metric = if (is_binary) "P(seroconvert)" else "E[value]",
               profile = pn, value = val)
  })
  bind_rows(res)
}

cat("Attaching held-out outcomes as leaf nodes and querying the L5 network...\n")
out <- bind_rows(
  query_leaf("igg_seroconvert_any", is_binary = TRUE),
  query_leaf("igg_log2fc_mean",     is_binary = FALSE),
  query_leaf("hai_peak_max",        is_binary = FALSE)
)

out$value <- round(out$value, 4)
cat("\n=== BN leaf-node conditional queries ===\n")
print(as.data.frame(out), row.names = FALSE)
write_csv(out, file.path(out_dir, "bn_leaf_cpqueries.csv"))

# ---- Contrast table: suppressed vs responsive, and per-root ------------------
wide <- out %>% tidyr::pivot_wider(names_from = profile, values_from = value)
cat("\n=== Key contrasts (BN inference) ===\n")
contrasts <- out %>%
  tidyr::pivot_wider(names_from = profile, values_from = value) %>%
  mutate(
    `Older_vs_Young`   = `Older (any CMV)` - `Young (any CMV)`,
    `CMVpos_vs_neg`    = `CMV+ (any age)`  - `CMV- (any age)`,
    `Male_vs_Female`   = `Male (any age/CMV)` - `Female (any age/CMV)`,
    `Suppressed_vs_Responsive` =
      `Older+CMV+ (suppressed)` - `Young+CMV- (responsive)`) %>%
  select(target, metric, Older_vs_Young, CMVpos_vs_neg, Male_vs_Female,
         Suppressed_vs_Responsive)
print(as.data.frame(contrasts), row.names = FALSE, digits = 3)
write_csv(contrasts, file.path(out_dir, "bn_leaf_contrasts.csv"))

cat("\nDone. Wrote bn_leaf_cpqueries.csv and bn_leaf_contrasts.csv\n")
