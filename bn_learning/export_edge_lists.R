library(bnlearn)
library(dplyr)

base_dir <- here::here("bn_learning")
out_dir <- file.path(base_dir, "output", "networks", "csv_networks")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

model_files <- list.files(file.path(base_dir, "output", "models"), pattern = "_models\\.rds$", full.names = TRUE)
cat(sprintf("Found %d model files\n\n", length(model_files)))

for (f in sort(model_files)) {
  net_id <- sub("_models\\.rds$", "", basename(f))
  mod <- readRDS(f)

  if (is.null(mod$boot_str) || is.null(mod$avg_conservative)) {
    cat(sprintf("SKIP  %s — missing boot_str or avg_conservative (old format)\n", net_id))
    next
  }

  final_edges <- data.frame(
    from = arcs(mod$avg_conservative)[, "from"],
    to   = arcs(mod$avg_conservative)[, "to"],
    stringsAsFactors = FALSE
  )

  literature_edges <- mod$boot_str %>%
    mutate(
      edge = paste(from, to, sep = " -> "),
      in_final_dag = edge %in% paste(final_edges$from, final_edges$to, sep = " -> ")
    ) %>%
    select(from, to, strength, direction, in_final_dag) %>%
    arrange(desc(in_final_dag), desc(strength))

  out_file <- file.path(out_dir, sprintf("%s_literature_edges.csv", net_id))
  write.csv(literature_edges, out_file, row.names = FALSE)

  n_dag <- sum(literature_edges$in_final_dag)
  cat(sprintf("OK    %-45s  %d rows, %d in_final_dag  -> %s\n",
              net_id, nrow(literature_edges), n_dag, out_file))
}
