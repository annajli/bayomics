library(bnlearn)
library(Rgraphviz)

base_dir <- here::here("bn_learning")
out_dir <- file.path(base_dir, "output", "networks", "visualized_networks")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

model_files <- list.files(file.path(base_dir, "output", "models"), pattern = "_models\\.rds$", full.names = TRUE)
cat(sprintf("Found %d model files\n\n", length(model_files)))

# Load each config to get ROOTS for highlighting
config_files <- list.files(file.path(base_dir, "config"), pattern = "\\.R$", full.names = TRUE)

for (f in sort(model_files)) {
  net_id <- sub("_models\\.rds$", "", basename(f))
  mod <- readRDS(f)

  # Find matching config to get ROOTS
  cfg_match <- config_files[grepl(net_id, basename(config_files), fixed = TRUE)]
  roots <- c("age_group", "sex", "cmv")
  if (length(cfg_match) == 1) {
    e <- new.env(parent = globalenv())
    tryCatch({
      source(cfg_match, local = e)
      if (exists("ROOTS", envir = e)) roots <- e$ROOTS
    }, error = function(err) NULL)
  }

  # Get nodes actually present in each network for highlight intersection
  plot_dag <- function(bn, title, out_file) {
    present_roots <- intersect(roots, nodes(bn))
    hl <- if (length(present_roots) > 0) list(nodes = present_roots, fill = "lightblue") else list()
    png(out_file, width = 1200, height = 1200)
    graphviz.plot(bn, main = title, shape = "ellipse", highlight = hl)
    dev.off()
  }

  # Tabu DAG
  if (!is.null(mod$bn_tabu)) {
    out_file <- file.path(out_dir, sprintf("%s_tabu.png", net_id))
    plot_dag(mod$bn_tabu, sprintf("%s — Tabu", net_id), out_file)
    cat(sprintf("OK  %s_tabu.png\n", net_id))
  }

  # Averaged (optimal threshold)
  if (!is.null(mod$avg_opt)) {
    thresh <- if (!is.null(mod$config$thresh_opt)) mod$config$thresh_opt
              else { t <- attr(mod$avg_opt, "threshold"); if (!is.null(t)) t else NA }
    out_file <- file.path(out_dir, sprintf("%s_avg_opt.png", net_id))
    plot_dag(mod$avg_opt, sprintf("%s — Averaged (t=%.3f)", net_id, thresh), out_file)
    cat(sprintf("OK  %s_avg_opt.png\n", net_id))
  }

  # Averaged (conservative threshold)
  if (!is.null(mod$avg_conservative)) {
    thresh <- mod$config$thresh_conservative
    out_file <- file.path(out_dir, sprintf("%s_avg_conservative.png", net_id))
    plot_dag(mod$avg_conservative, sprintf("%s — Averaged Conservative (t=%.3f)", net_id, thresh), out_file)
    cat(sprintf("OK  %s_avg_conservative.png\n", net_id))
  }

  cat("\n")
}

cat(sprintf("All plots saved to %s\n", out_dir))
