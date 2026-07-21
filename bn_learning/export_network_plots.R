library(bnlearn)

base_dir <- here::here("bn_learning")
out_dir <- file.path(base_dir, "output", "networks", "visualized_networks")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

model_files <- list.files(file.path(base_dir, "output", "models"), pattern = "_models\\.rds$", full.names = TRUE)
cat(sprintf("Found %d model files\n\n", length(model_files)))

for (f in sort(model_files)) {
  net_id <- sub("_models\\.rds$", "", basename(f))
  mod <- readRDS(f)

  # Tabu DAG
  if (!is.null(mod$bn_tabu)) {
    out_file <- file.path(out_dir, sprintf("%s_tabu.png", net_id))
    png(out_file, width = 1200, height = 1200)
    plot(mod$bn_tabu, main = sprintf("%s — Tabu", net_id))
    dev.off()
    cat(sprintf("OK  %s_tabu.png\n", net_id))
  }

  # Averaged (optimal threshold)
  if (!is.null(mod$avg_opt)) {
    thresh <- if (!is.null(mod$config$thresh_opt)) mod$config$thresh_opt
              else { t <- attr(mod$avg_opt, "threshold"); if (!is.null(t)) t else NA }
    out_file <- file.path(out_dir, sprintf("%s_avg_opt.png", net_id))
    png(out_file, width = 1200, height = 1200)
    plot(mod$avg_opt, main = sprintf("%s — Averaged (t=%.3f)", net_id, thresh))
    dev.off()
    cat(sprintf("OK  %s_avg_opt.png\n", net_id))
  }

  # Averaged (conservative threshold)
  if (!is.null(mod$avg_conservative)) {
    thresh <- mod$config$thresh_conservative
    out_file <- file.path(out_dir, sprintf("%s_avg_conservative.png", net_id))
    png(out_file, width = 1200, height = 1200)
    plot(mod$avg_conservative, main = sprintf("%s — Averaged Conservative (t=%.3f)", net_id, thresh))
    dev.off()
    cat(sprintf("OK  %s_avg_conservative.png\n", net_id))
  }

  cat("\n")
}

cat(sprintf("All plots saved to %s\n", out_dir))
