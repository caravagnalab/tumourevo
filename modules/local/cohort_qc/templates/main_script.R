#!/usr/bin/env Rscript

pkgs <- c("tidyverse", "ComplexHeatmap", "circlize", "scales", "grid", "patchwork")
sapply(pkgs, require, character.only = TRUE)

parse_args <- function(x) {
  # Remove brackets
  x = gsub("\\\\[","",x)
  x = gsub("\\\\]","",x)

  # Split into key:value pairs
  args_list = unlist(strsplit(x, ", (?=[^)]*(?:\\\\(|\$))", perl=TRUE))

  # Convert to tibble and split key:value
  tibble(arg = args_list) %>%
    separate(arg, into = c("key", "value"), sep = ":", extra = "merge", fill = "right") %>%
    dplyr::mutate(across(everything(), ~str_trim(.))) %>%
    dplyr::filter(!is.na(key) & !is.na(value)) %>%
    tibble::deframe()
}

opt = list(
    prefix = ifelse('$task.ext.prefix' == 'null', '$meta.id', '$task.ext.prefix')
)

args_opt = parse_args('$task.ext.args')
for ( ao in names(args_opt)) opt[[ao]] = args_opt[[ao]]

prefix <- opt[["prefix"]]


### Helpers ###

safe_mean <- function(x) {
  if (is.null(x) || length(x) == 0 || all(is.na(x))) return(NA_real_)
  mean(x, na.rm = TRUE)
}

safe_median <- function(x) {
  if (is.null(x) || length(x) == 0 || all(is.na(x))) return(NA_real_)
  median(x, na.rm = TRUE)
}

# IDs parser
strip_suffix <- function(x, suffix) {
  ifelse(endsWith(x, suffix), substr(x, 1, nchar(x) - nchar(suffix)), x)
}

extract_id <- function(x, pattern, default = NA_character_) {
  out <- stringr::str_extract(x, pattern)
  ifelse(is.na(out), default, out)
}

get_ids <- function(file,
                    suffix,
                    patient_pattern = "U[0-9]+",
                    cohort_pattern = "^[^_]+") {
  fname <- basename(file)
  sample_id <- strip_suffix(fname, suffix)

  tibble::tibble(
    cohort_id  = extract_id(sample_id, cohort_pattern),
    patient_id = extract_id(sample_id, patient_pattern),
    sample_id  = sample_id
  )
}

# TIN classification
classify_tin <- function(x) {
  x <- 100 * x
  dplyr::case_when(
    is.na(x) ~ NA_character_,
    x < 1 ~ "No",
    x >= 1 & x < 7 ~ "Low",
    x >= 7 & x < 15 ~ "Contamination",
    TRUE ~ "High"
  )
}

### TINC summary extractor ###

get_tinc_summary <- function(fit, sample_id = NULL) {
  sid <- sample_id
  if (is.null(sid)) {
    sid <- if (!is.null(fit[["sample"]])) fit[["sample"]] else NA_character_
  }
  tin <- fit[["TIN"]]
  tin_pct <- if (!is.null(tin) && !is.na(tin)) 100 * tin else NA_real_
  tibble::tibble(
    sample_id = as.character(sid),
    TIN = tin,
    TIN_pct = tin_pct,
    TIN_class = classify_tin(tin)
  )
}

### CNAqc summary extractor ###

get_cnaqc_summary <- function(qc, sample_id = NULL) {
  sid <- sample_id
  if (is.null(sid)) {
    sid <- if (!is.null(qc[["sample"]])) qc[["sample"]] else NA_character_
  }
  
  # mutation QC
  mut_qc <- qc[["mutations"]][["QC_PASS"]]
  mut_tested <- sum(mut_qc %in% c(TRUE, FALSE), na.rm = TRUE)
  mut_pass <- sum(mut_qc == TRUE, na.rm = TRUE)
  mut_fail <- sum(mut_qc == FALSE, na.rm = TRUE)
  mut_na <- sum(is.na(mut_qc))
  mutation_pass_rate <- if (mut_tested > 0) mut_pass / mut_tested else NA_real_
  mutation_fail_rate <- if (mut_tested > 0) mut_fail / mut_tested else NA_real_
  mutation_na_fraction <- if (length(mut_qc) > 0) mut_na / length(mut_qc) else NA_real_
  
  # CNA QC
  cna_qc <- qc[["cna"]][["QC_PASS"]]
  cna_tested <- sum(cna_qc %in% c(TRUE, FALSE), na.rm = TRUE)
  cna_pass <- sum(cna_qc == TRUE, na.rm = TRUE)
  cna_fail <- sum(cna_qc == FALSE, na.rm = TRUE)
  cna_na <- sum(is.na(cna_qc))
  cna_pass_rate <- if (cna_tested > 0) cna_pass / cna_tested else NA_real_
  cna_fail_rate <- if (cna_tested > 0) cna_fail / cna_tested else NA_real_
  cna_na_fraction <- if (length(cna_qc) > 0) cna_na / length(cna_qc) else NA_real_
  
  # peaks analysis
  peak_score <- if (!is.null(qc[["peaks_analysis"]][["score"]])) qc[["peaks_analysis"]][["score"]] else NA_real_
  peak_qc <- if (!is.null(qc[["peaks_analysis"]][["QC"]])) qc[["peaks_analysis"]][["QC"]] else NA_character_
  
  matches <- qc[["peaks_analysis"]][["matches"]]
  if (!is.null(matches) && nrow(matches) > 0) {
    n_peaks <- nrow(matches)
    n_peaks_pass <- sum(matches[["QC"]] == "PASS", na.rm = TRUE)
    n_peaks_fail <- sum(matches[["QC"]] == "FAIL", na.rm = TRUE)
    peak_pass_rate <- n_peaks_pass / n_peaks
    mean_delta_vaf <- safe_mean(matches[["delta_vaf"]])
    median_counts_per_bin <- safe_median(matches[["counts_per_bin"]])
    
    dominant_karyotype <- matches[["karyotype"]][which.max(matches[["weight"]])][1]
    dominant_karyotype_weight <- max(matches[["weight"]], na.rm = TRUE)
  } else {
    n_peaks <- NA_integer_
    n_peaks_pass <- NA_integer_
    n_peaks_fail <- NA_integer_
    peak_pass_rate <- NA_real_
    mean_delta_vaf <- NA_real_
    median_counts_per_bin <- NA_real_
    dominant_karyotype <- NA_character_
    dominant_karyotype_weight <- NA_real_
  }
  
  # CN structure
  n_mutations <- qc[["n_mutations"]]
  n_cna <- qc[["n_cna"]]
  n_cna_clonal <- qc[["n_cna_clonal"]]
  n_cna_subclonal <- qc[["n_cna_subclonal"]] 
  purity <- qc[["purity"]]
  ploidy <- qc[["ploidy"]]
  most_prevalent_karyotype <- qc[["most_prevalent_karyotype"]]
   
  # Assign QC class
  qc_class <- dplyr::case_when(
    !is.na(mutation_pass_rate) && !is.na(cna_pass_rate) &&
      mutation_pass_rate >= 0.80 && cna_pass_rate >= 0.80 &&
      (is.na(peak_pass_rate) || peak_pass_rate >= 0.50) ~ "PASS",
    !is.na(mutation_pass_rate) && !is.na(cna_pass_rate) &&
      (mutation_pass_rate < 0.50 || cna_pass_rate < 0.50 ||
         (!is.na(peak_pass_rate) && peak_pass_rate < 0.50)) ~ "FAIL",
    TRUE ~ "WARN"
  )
  
  data.frame(
    sample_id = sid,
    n_mutations = n_mutations,
    purity = purity,
    ploidy = ploidy,
    n_cna = n_cna,
    n_cna_clonal = n_cna_clonal,
    n_cna_subclonal = n_cna_subclonal,
    most_prevalent_karyotype = most_prevalent_karyotype,
    mutation_pass_rate = mutation_pass_rate,
    mutation_fail_rate = mutation_fail_rate,
    mutation_na_fraction = mutation_na_fraction,
    cna_pass_rate = cna_pass_rate,
    cna_fail_rate = cna_fail_rate,
    cna_na_fraction = cna_na_fraction,
    peak_score = peak_score,
    peak_qc = peak_qc,
    n_peaks = n_peaks,
    n_peaks_pass = n_peaks_pass,
    n_peaks_fail = n_peaks_fail,
    peak_pass_rate = peak_pass_rate,
    mean_delta_vaf = mean_delta_vaf,
    median_counts_per_bin = median_counts_per_bin,
    dominant_karyotype = dominant_karyotype,
    dominant_karyotype_weight = dominant_karyotype_weight,
    qc_class = qc_class
  )
}


### RDS file loader ###

files <- list.files(".", full.names = TRUE)

cnaqc_rds_files <- files[endsWith(files, "_qc.rds")]
tinc_rds_files  <- files[endsWith(files, "_fit.rds")]

load_qc_objects <- function(files, suffix, extractor_fun) {
  purrr::map_dfr(files, function(f) {
    obj <- readRDS(f)
    ids <- get_ids(f, suffix = suffix)
    summary <- extractor_fun(obj, sample_id = ids[["sample_id"]])
    dplyr::left_join(ids, summary, by = "sample_id")
  })
}

### Load TINC and CNAqc objects and get cohort summary ###

cohort_qc_tinc <- load_qc_objects(
  files = tinc_rds_files,
  suffix = "_fit.rds",
  extractor_fun = get_tinc_summary
)

cohort_qc_cnaqc <- load_qc_objects(
  files = cnaqc_rds_files,
  suffix = "_qc.rds",
  extractor_fun = get_cnaqc_summary
)

# Merge into unified cohort summary 

cohort_qc_all <- cohort_qc_cnaqc %>%
  dplyr::full_join(cohort_qc_tinc, by = c("patient_id", "sample_id"))

cohort_qc <- cohort_qc_all %>%
  dplyr::mutate(
    TIN_class = factor(
      TIN_class,
      levels = c("No", "Low", "Contamination", "High")
    )
  )

saveRDS(cohort_qc, sprintf("%s.qc_summary.rds", prefix))


### Plot settings ###

qc_colors <- c(
  PASS = "#1b9e77",
  WARN = "#e6ab02",
  FAIL = "#F8766D"
)

theme_dash <- ggplot2::theme_bw(base_size = 11) +
  ggplot2::theme(
    plot.title = ggplot2::element_text(face = "bold", size = 11),
    axis.title = ggplot2::element_text(size = 10)
  )

# Summary strip

summary_df <- cohort_qc %>%
  dplyr::summarise(
    n_samples = dplyr::n(),
    n_pass = sum(qc_class == "PASS", na.rm = TRUE),
    n_warn = sum(qc_class == "WARN", na.rm = TRUE),
    n_fail = sum(qc_class == "FAIL", na.rm = TRUE),
    med_purity = median(purity, na.rm = TRUE),
    med_ploidy = median(ploidy, na.rm = TRUE),
    med_tin = median(TIN_pct, na.rm = TRUE)
  )

summary_text <- paste0(
  "Samples: ", summary_df["n_samples"],
  "    |    PASS: ", summary_df[["n_pass"]],
  "    |    WARN: ", summary_df[["n_warn"]],
  "    |    FAIL: ", summary_df[["n_fail"]],
  "\nMedian purity: ", round(summary_df[["med_purity"]], 2),
  "    |    Median ploidy: ", round(summary_df[["med_ploidy"]], 2),
  "    |    Median TIN: ", round(summary_df[["med_tin"]], 2), "%"
)

summary_plot <- ggplot2::ggplot() +
  ggplot2::theme_void() +
  ggplot2::annotate(
    "text",
    x = 0, y = 1,
    label = summary_text,
    hjust = 0, vjust = 1,
    size = 4.5,
    fontface = "plain"
  ) +
  ggplot2::xlim(0, 1) + ggplot2::ylim(0, 1)


### Main plots ###

p1 <- ggplot2::ggplot(cohort_qc, ggplot2::aes(qc_class, fill = qc_class)) +
  ggplot2::geom_bar() +
  ggplot2::scale_fill_manual(values = qc_colors) +
  ggplot2::geom_text(
    stat = "count",
    ggplot2::aes(label = after_stat(count)),
    vjust = -0.3,
    size = 4
  ) +
  theme_dash +
  ggplot2::labs(
    title = "QC classification",
    x = "QC class",
    y = "Number of samples"
  )

p2 <- ggplot2::ggplot(cohort_qc, ggplot2::aes(n_mutations)) +
  ggplot2::geom_histogram(bins = 30, fill = "#4682B4") +
  ggplot2::scale_x_log10(
    breaks = 10^(0:7),
    labels = scales::trans_format("log10", scales::math_format(10^.x))
  ) +
  theme_dash +
  ggplot2::labs(
    title = "Mutation burden",
    x = "Mutations (log scale)",
    y = "Number of samples"
  )

  p3 <- ggplot2::ggplot(cohort_qc, ggplot2::aes(purity, fill = qc_class)) +
  ggplot2::geom_histogram(
    bins = 30,
    position = "stack",
    alpha = 0.7,
    color = "grey20",
    linewidth = 0.2
  ) +
  ggplot2::scale_fill_manual(values = qc_colors, name = "QC class") +
  theme_dash +
  ggplot2::labs(
    title = "Purity",
    x = "Purity",
    y = "Number of samples"
  )

p4 <- ggplot2::ggplot(cohort_qc, ggplot2::aes(ploidy)) +
  ggplot2::geom_histogram(
    bins = 30,
    fill = "#6A6599FF",
    color = "grey20",
    linewidth = 0.2
  ) +
  theme_dash +
  ggplot2::labs(
    title = "Ploidy",
    x = "Ploidy",
    y = "Number of samples"
  )


p5 <- ggplot2::ggplot(cohort_qc, ggplot2::aes(mutation_pass_rate, cna_pass_rate)) +
  ggplot2::annotate("rect", xmin = 0, xmax = 0.5, ymin = 0, ymax = 1, fill = qc_colors["FAIL"], alpha = 0.08) +
  ggplot2::annotate("rect", xmin = 0, xmax = 1, ymin = 0, ymax = 0.5, fill = qc_colors["FAIL"], alpha = 0.08) +
  ggplot2::annotate("rect", xmin = 0.5, xmax = 0.8, ymin = 0.5, ymax = 0.8, fill = qc_colors["WARN"], alpha = 0.08) +
  ggplot2::annotate("rect", xmin = 0.8, xmax = 1, ymin = 0.8, ymax = 1, fill = qc_colors["PASS"], alpha = 0.1) +
  ggplot2::geom_point(ggplot2::aes(color = qc_class), alpha = 0.9, size = 1.6) +
  ggplot2::geom_vline(xintercept = c(0.5, 0.8), linetype = c("dotted", "dashed"), color = "grey30") +
  ggplot2::geom_hline(yintercept = c(0.5, 0.8), linetype = c("dotted", "dashed"), color = "grey30") +
  ggplot2::scale_color_manual(values = qc_colors, name = "QC class") +
  ggplot2::coord_equal(xlim = c(0, 1), ylim = c(0, 1)) +
  theme_dash +
  ggplot2::labs(
    title = "Mutation QC vs CNA QC",
    subtitle = "FAIL (<0.5), WARN (0.5–0.8), PASS (>0.8)",
    x = "Mutation QC pass rate",
    y = "CNA QC pass rate"
  )

p6 <- ggplot(cohort_qc, aes(peak_pass_rate)) +
  geom_histogram(bins = 30, fill = "#B07AA1") +
  geom_vline(xintercept = 0.49, linetype = "dashed", color = "grey40") +
  theme_dash +
  labs(title = "Peak pass rate", x = "Pass rate", y = "Number of samples")

 

p7 <- ggplot2::ggplot(cohort_qc, ggplot2::aes(peak_pass_rate)) +
  ggplot2::geom_histogram(bins = 30, fill = "#B07AA1") +
  ggplot2::geom_vline(xintercept = 0.49, linetype = "dashed", color = "grey40") +
  theme_dash +
  ggplot2::labs(title = "Peak pass rate", x = "Pass rate", y = "Number of samples")

  p8 <- ggplot2::ggplot(cohort_qc, ggplot2::aes(TIN_pct)) +
  ggplot2::geom_histogram(bins = 40, fill = "#D95F02", color = "white") +
  ggplot2::geom_vline(xintercept = c(1, 7, 15), linetype = "dashed", color = "grey40") +
  ggplot2::annotate("text", x = 0.5, y = Inf, label = "No contamination", vjust = 2, size = 3.5) +
  ggplot2::annotate("text", x = 4, y = Inf, label = "Low", vjust = 2, size = 3.5) +
  ggplot2::annotate("text", x = 11, y = Inf, label = "Contamination", vjust = 2, size = 3.5) +
  ggplot2::annotate("text", x = 20, y = Inf, label = "High", vjust = 2, size = 3.5) +
  theme_dash +
  ggplot2::labs(
    title = "Tumour-in-normal contamination",
    subtitle = "TIN classification: <1%, 1–7%, 7–15%, >15%",
    x = "TIN (%)",
    y = "Number of samples"
  )

### Wrap plots ###

first_grid  <- (p1 | p5 | p6)
second_grid <- (p2 | p3 | p4 | p7)

dashboard <- summary_plot / first_grid / second_grid +
  patchwork::plot_layout(heights = c(0.08, 0.45, 0.45)) +
  patchwork::plot_annotation(
    title = paste0(cohort_id, " cohort QC summary"),
    theme = ggplot2::theme(
      plot.title = ggplot2::element_text(size = 16, face = "bold", hjust = 0.5)
    )
  )

ggplot2::ggsave(
  filename = sprintf("%s.qc_report.pdf", prefix),
  plot = dashboard,
  width = 14,
  height = 16,
  dpi = 400
)

saveRDS(dashboard, sprintf("%s.qc_plot.rds", prefix))

# version export
f = file("versions.yml","w")
tidyverse_version = sessionInfo()\$otherPkgs\$tidyverse\$Version
ComplexHeatmap_version = sessionInfo()\$otherPkgs\$ComplexHeatmap\$Version
circlize_version = sessionInfo()\$otherPkgs\$circlize\$Version
scales_version = sessionInfo()\$otherPkgs\$scales\$Version
grid_version = sessionInfo()\$otherPkgs\$grid\$Version
patchwork_version = sessionInfo()\$otherPkgs\$patchwork\$Version
writeLines(paste0('"', "$task.process", '"', ":"), f)
writeLines(paste("    tidyverse:", tidyverse_version), f)
writeLines(paste("    complexheatmap:", ComplexHeatmap_version), f)
writeLines(paste("    circlize:", circlize_version), f)
writeLines(paste("    scales:", scales_version), f)
writeLines(paste("    grid:", grid_version), f)
writeLines(paste("    patchwork", patchwork_version), f)
close(f)
