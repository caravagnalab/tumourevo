#!/usr/bin/env Rscript
# Rescore subclonal clusters from a finished tumourevo results folder, without Nextflow.
#
#   Rscript /path/to/tumourevo/rescore_subclonal.R RESULTS_DIR [OUT_DIR] [ARGS] [PATIENTS]
#
#   RESULTS_DIR     the pipeline --outdir (contains assign_signature/, subclonal_deconvolution/, formatter/)
#   OUT_DIR         default RESULTS_DIR/rescored; outputs go to OUT_DIR/<dataset>/<patient>/
#   ARGS            scoring options, default "[tools:mobster,pyclone-vi,viber,sigprofiler]";
#                   add any option without the interpreter_ prefix, e.g. "[tools:mobster,pyclone-vi,sigprofiler, read_support:4]"
#   PATIENTS        optional comma-separated patient IDs
# Set VIBER_HEURISTIC=true to use the heuristic VIBER fit.
# The scoring code is read from modules/local/subclonal_interpretation in the repository holding this script.

a = commandArgs(trailingOnly = TRUE)
if (length(a) < 1) stop("Usage: Rscript rescore_subclonal.R RESULTS_DIR [OUT_DIR] [ARGS] [PATIENTS]")
script   = sub("^--file=", "", grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)[1])
repo     = dirname(normalizePath(script))
results  = normalizePath(a[1])
template = file.path(repo, "modules/local/subclonal_interpretation/templates/main_script.R")
out_dir  = if (length(a) >= 2 && nzchar(a[2])) a[2] else file.path(results, "rescored")
args     = if (length(a) >= 3 && nzchar(a[3])) a[3] else "[tools:mobster,pyclone-vi,viber,sigprofiler]"
patients = if (length(a) >= 4 && nzchar(a[4])) strsplit(a[4], ",")[[1]] else NULL
viber    = if (Sys.getenv("VIBER_HEURISTIC") == "true") "*_viber_best_heuristic_fit.rds" else "*_viber_best_fit.rds"
if (!file.exists(template)) stop("Template not found: ", template)

paths = function(pattern, dirs = FALSE) {
  hits = sort(Sys.glob(file.path(results, pattern)))
  paste(hits[dir.exists(hits) == dirs], collapse = " ")
}

for (pdir in Sys.glob(file.path(results, "assign_signature/prepare_cluster/*/*"))) {
  ds = basename(dirname(pdir)); pt = basename(pdir)
  if (!is.null(patients) && !pt %in% patients) next

  inputs = c(
    "$mutation_tables"     = paths(sprintf("assign_signature/prepare_cluster/%s/%s/*.txt", ds, pt)),
    "$results_sigprofiler" = paths(sprintf("assign_signature/assign_cluster/%s/%s/*", ds, pt), dirs = TRUE),
    "$fits"                = trimws(paste(
                               paths(sprintf("subclonal_deconvolution/pyclonevi/%s/%s/*_best_fit.txt", ds, pt)),
                               paths(sprintf("subclonal_deconvolution/viber/%s/%s/%s", ds, pt, viber)),
                               paths(sprintf("subclonal_deconvolution/mobster/%s/%s/*/*_mobster_best_fit.rds", ds, pt)))),
    "$joint_table"         = paths(sprintf("formatter/cnaqc2tsv/%s/%s/*_joint_table.tsv", ds, pt)))
  counts = sapply(strsplit(inputs, " +"), function(x) sum(nzchar(x)))
  message(sprintf("%s_%s: %s", ds, pt, paste(sub("^[$]", "", names(inputs)), counts, sep = "=", collapse = ", ")))
  if (counts[1] == 0 || counts[2] == 0) { message("  skipped: no cluster tables or SigProfiler results"); next }

  # Fill in the Nextflow template placeholders and undo its escapes
  src = paste(readLines(template), collapse = "\n")
  fill = c("$task.ext.prefix" = "null", "$task.ext.args" = args, "$meta.id" = paste0(ds, "_", pt), inputs)
  for (k in names(fill)) src = gsub(k, fill[[k]], src, fixed = TRUE)
  src = gsub("\\$", "$", src, fixed = TRUE)
  src = gsub("\\\\", "\\", src, fixed = TRUE)

  dest = file.path(out_dir, ds, pt)
  dir.create(dest, recursive = TRUE, showWarnings = FALSE)
  writeLines(src, file.path(dest, "run.R"))

  # One R process per patient, run inside its output folder; a failure does not stop the others
  old = setwd(dest)
  status = system2(file.path(R.home("bin"), "Rscript"), "run.R", stdout = "run.log", stderr = "run.log")
  setwd(old)
  message(if (status == 0) paste("  done:", dest) else paste("  FAILED, see", file.path(dest, "run.log")))
}
