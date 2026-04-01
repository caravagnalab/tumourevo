#!/usr/bin/env Rscript

pkgs <- c("tidyverse", "circlize", "scales", "grid", "patchwork")
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

get_id <- function(x, pattern, default = NA_character_) {
  out <- stringr::str_extract(x, pattern)
  ifelse(is.na(out), default, out)
}

get_ids <- function(file,
                    suffix,
                    patient_pattern = "U[0-9]+"
		    ){
  fname <- basename(file)
  sample_id <- strip_suffix(fname, suffix)

  tibble::tibble(
    patient_id = get_id(sample_id, patient_pattern),
    sample_id  = sample_id
  )
}

get_n_karyotype <- function(qc, sample_id = NULL) {
  sid <- sample_id
  if (is.null(sid)) {
    sid <- if (!is.null(qc[["sample"]])) qc[["sample"]] else NA_character_
  }

  nk <- qc[["n_karyotype"]]

  if (is.null(nk) || length(nk) == 0) {
    return(tibble::tibble(
      sample_id = character(),
      karyotype = character(),
      n_mutations = numeric()
    ))
  }

  tibble::tibble(
    sample_id = sid,
    karyotype = names(nk),
    n_mutations = as.numeric(nk)
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
  mut_na_fraction <- if (length(mut_qc) > 0) mut_na / length(mut_qc) else NA_real_
  
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
   
    dominant_karyotype <- matches[["karyotype"]][which.max(matches[["weight"]])][1]
    dominant_karyotype_weight <- max(matches[["weight"]], na.rm = TRUE)
  } else {
    n_peaks <- NA_integer_
    n_peaks_pass <- NA_integer_
    n_peaks_fail <- NA_integer_
    peak_pass_rate <- NA_real_
    dominant_karyotype <- NA_character_
    dominant_karyotype_weight <- NA_real_
  }
  
  # CN structure
  n_mutations <- qc[["n_mutations"]]
  n_cna <- qc[["n_cna"]]
  purity <- qc[["purity"]]
  ploidy <- qc[["ploidy"]]
  most_prevalent_karyotype <- qc[["most_prevalent_karyotype"]]
  
  # CNAqc-oriented QC classification
  qc_class <- dplyr::case_when(
    !is.na(cna_pass_rate) &&
      cna_pass_rate >= 0.80 &&
      (is.na(peak_pass_rate) || peak_pass_rate >= 0.50) ~ "PASS",
    
    !is.na(cna_pass_rate) &&
      (cna_pass_rate < 0.50 ||
         (!is.na(peak_pass_rate) && peak_pass_rate < 0.50)) ~ "FAIL",
    
    TRUE ~ "WARN"
  )
  
  tibble::tibble(
    sample_id = sid,
    n_mutations = n_mutations,
    purity = purity,
    ploidy = ploidy,
    n_cna = n_cna,
    most_prevalent_karyotype = most_prevalent_karyotype,
    mut_na_fraction = mut_na_fraction,
    cna_pass_rate = cna_pass_rate,
    cna_fail_rate = cna_fail_rate,
    cna_na_fraction = cna_na_fraction,
    peak_score = peak_score,
    peak_qc = peak_qc,
    n_peaks = n_peaks,
    n_peaks_pass = n_peaks_pass,
    n_peaks_fail = n_peaks_fail,
    peak_pass_rate = peak_pass_rate,
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

load_karyotype_objects <- function(files, suffix) {
  purrr::map_dfr(files, function(f) {
    obj <- readRDS(f)
    ids <- get_ids(f, suffix = suffix)

    get_n_karyotype(obj, sample_id = ids[["sample_id"]]) %>%
      dplyr::mutate(patient_id = ids[["patient_id"]]) %>%
      dplyr::select(patient_id, sample_id, karyotype, n_mutations)
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

cohort_karyotype <- load_karyotype_objects(
  files = cnaqc_rds_files,
  suffix = "_qc.rds"
)

# Simmarize across samples
all_samples <- unique(cohort_karyotype[["sample_id"]])
all_karyotypes <- unique(cohort_karyotype[["karyotype"]])

cohort_karyotype_complete <- tidyr::expand_grid(
  sample_id = all_samples,
  karyotype = all_karyotypes
) %>%
  dplyr::left_join(
    cohort_karyotype,
    by = c("sample_id", "karyotype")
  ) %>%
  dplyr::mutate(
    n_mutations = dplyr::coalesce(n_mutations, 0)
  )

karyotype_summary <- cohort_karyotype_complete %>%
  dplyr::group_by(karyotype) %>%
  dplyr::summarise(
    total_mutations = sum(n_mutations, na.rm = TRUE),
    mean_mutations = mean(n_mutations, na.rm = TRUE),
    median_mutations = median(n_mutations, na.rm = TRUE),
    n_samples = dplyr::n(),
    n_nonzero_samples = sum(n_mutations > 0, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  dplyr::arrange(dplyr::desc(total_mutations))

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
  summarise(
    n_samples = n(),
    mean_purity = mean(purity, na.rm = TRUE),
    mean_ploidy = mean(ploidy, na.rm = TRUE),
    mean_tin = mean(TIN_pct, na.rm = TRUE)
  )

# summary strip

summary_text <- paste0(
  "Samples: ", summary_df["n_samples"],
  "\nMean purity: ", round(summary_df["mean_purity"], 3),
  "    |    Mean ploidy: ", round(summary_df["mean_ploidy"], 3),
  "    |    Mean TIN: ", round(summary_df["mean_tin"], 2), "%"
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

# QC classification
p1 <- ggplot(cohort_qc, aes(qc_class, fill = qc_class)) +
  geom_bar() +
  scale_fill_manual(
    values = qc_colors
  ) +
  geom_text(
    stat = "count",
    aes(label = after_stat(count)),
    vjust = -0.3,
    size = 4
  ) +
  theme_dash +
  labs(
    title = "QC classification",
    x = "QC class",
    y = "Number of samples"
  )

# Mutation burden
p2 <- ggplot(cohort_qc, aes(n_mutations)) +
  geom_histogram(bins = 30,
                 fill = "#4682B4") +
  scale_x_log10(
    breaks = 10^(0:7),
    labels = trans_format("log10", math_format(10^.x))
  ) +
  theme_dash +
  labs(
    title = "Mutation burden",
    x = "Mutations (log scale)",
    y = "Number of samples"
  )

# Purity
p3 <- ggplot(cohort_qc, aes(purity, fill = qc_class)) +
  geom_histogram(
    bins = 30,
    position = "stack",
    alpha = 0.7,
    color = "grey20",
    linewidth = 0.2
  ) +
  scale_fill_manual(values = qc_colors,
                    name = "QC class") +
  theme_dash +
  labs(
    title = "Purity",
    x = "Purity",
    y = "Number of samples"
  )

# Ploidy
p4 <- ggplot(cohort_qc, aes(ploidy)) +
  geom_histogram(
    bins = 30,
    fill = "#6A6599FF",
    color = "grey20",
    linewidth = 0.2
  ) +
  theme_dash +
  labs(
    title = "Ploidy",
    x = "Ploidy",
    y = "Number of samples"
  )

# CNA segments per sample
p5 <- ggplot2::ggplot(cohort_qc, ggplot2::aes(n_cna, fill = qc_class)) +
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
    title = "CNA segments per sample",
    x = "Number of CNA segments",
    y = "Number of samples"
  )

# CNA QC vs Peak QC
p6 <- ggplot2::ggplot(cohort_qc, ggplot2::aes(cna_pass_rate, peak_pass_rate)) +
  ggplot2::annotate("rect", xmin = 0, xmax = 0.5, ymin = 0, ymax = 1,
                    fill = qc_colors["FAIL"], alpha = 0.08) +
  ggplot2::annotate("rect", xmin = 0.5, xmax = 0.8, ymin = 0.5, ymax = 1,
                    fill = qc_colors["WARN"], alpha = 0.08) +
  ggplot2::annotate("rect", xmin = 0.8, xmax = 1, ymin = 0.5, ymax = 1,
                    fill = qc_colors["PASS"], alpha = 0.10) +
  ggplot2::geom_point(ggplot2::aes(color = qc_class), alpha = 0.9, size = 1.8) +
  ggplot2::geom_vline(xintercept = c(0.5, 0.8), linetype = c("dotted", "dashed"), color = "grey30") +
  ggplot2::geom_hline(yintercept = 0.5, linetype = "dashed", color = "grey30") +
  ggplot2::scale_color_manual(values = qc_colors, name = "QC class") +
  ggplot2::coord_equal(xlim = c(0, 1), ylim = c(0, 1)) +
  theme_dash +
  ggplot2::labs(
    title = "CNA QC vs peak QC",
    x = "CNA QC pass rate",
    y = "Peak pass rate"
  )

# Karyotype stats
p7 <- ggplot2::ggplot(
  karyotype_summary,
  ggplot2::aes(
    x = reorder(karyotype, total_mutations),
    y = total_mutations
  )
) +
  ggplot2::geom_col(fill = "#CC6677") +
  ggplot2::coord_flip() +
  theme_dash +
  ggplot2::labs(
    title = "Mutations per karyotype",
    x = "Karyotype",
    y = "Total mutations across cohort"
  )

# TIN distribution
p8 <- ggplot(cohort_qc, aes(TIN_pct)) +
  geom_histogram(bins = 40, fill = "#D95F02", color = "white") +
  
  geom_vline(xintercept = c(1, 7, 15),
             linetype = "dashed", color = "grey40") +
  
  # labels
  annotate("text", x = 0.5, y = Inf, label = "No contamination",
           vjust = 2, size = 3.5) +
  annotate("text", x = 4, y = Inf, label = "Low",
           vjust = 2, size = 3.5) +
  annotate("text", x = 11, y = Inf, label = "Contamination",
           vjust = 2, size = 3.5) +
  annotate("text", x = 20, y = Inf, label = "High",
           vjust = 2, size = 3.5) +
  
  theme_dash +
  labs(
    title = "Tumour-in-normal contamination",
    subtitle = "TIN classification: <1%, 1–7%, 7–15%, >15%",
    x = "TIN (%)",
    y = "Number of samples"
  )

# Wrap plots

first_grid  <- (p1 | p6)
second_grid <- (p2 | p3 | p4)
third_grid  <- (p5 | p7 | p8)

dashboard <- summary_plot / first_grid / second_grid / third_grid +
  patchwork::plot_layout(heights = c(0.10, 0.28, 0.30, 0.32)) +
  patchwork::plot_annotation(
    title = paste0(prefix, " cohort CNAqc / TINC summary"),
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
