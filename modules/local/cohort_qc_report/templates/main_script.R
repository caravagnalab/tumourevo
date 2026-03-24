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
cohort_id <- prefix


### Helpers ###

safe_mean <- function(x) {
  if (is.null(x) || length(x) == 0 || all(is.na(x))) return(NA_real_)
  mean(x, na.rm = TRUE)
}

safe_median <- function(x) {
  if (is.null(x) || length(x) == 0 || all(is.na(x))) return(NA_real_)
  median(x, na.rm = TRUE)
}


get_ids <- function(file, suffix) {
  fname <- basename(file)

  sample_id <- if (endsWith(fname, suffix)) {
    substr(fname, 1, nchar(fname) - nchar(suffix))
  } else {
    fname
  }

  patient_id <- stringr::str_extract(sample_id, "U[0-9]+")

  tibble::tibble(
    patient_id = patient_id,
    sample_id = sample_id
  )
}

classify_tin <- function(x) {
  x <- 100 * x
  dplyr::case_when(
    is.na(x) ~ NA_character_,
    x < 1 ~ "No contamination",
    x >= 1 & x < 7 ~ "Low contamination",
    x >= 7 & x < 15 ~ "Contamination",
    TRUE ~ "High contamination"
  )
}

classify_tit <- function(x) {
  x <- 100 * x

  dplyr::case_when(
    is.na(x) ~ NA_character_,
    x < 15 ~ "Very low purity",
    x >= 15 & x < 45 ~ "Bad purity",
    x >= 45 & x < 65 ~ "Average purity",
    x >= 65 & x < 85 ~ "Good purity",
    x >= 85 ~ "Very high purity"
  )
}

### TINC summary extractor ###

get_tinc_summary <- function(fit, sample_id = NULL) {
  sid <- sample_id
  if (is.null(sid)) {
    sid <- if (!is.null(fit[["sample"]])) fit[["sample"]] else NA_character_
  }

  tin <- fit[["TIN"]]
  tit <- fit[["TIT"]]

  tin_pct <- if (!is.null(tin) && !is.na(tin)) 100 * tin else NA_real_
  tit_pct <- if (!is.null(tit) && !is.na(tit)) 100 * tit else NA_real_

  tibble::tibble(
    sample_id = as.character(sid),
    TIN = tin,
    TIT = tit,
    TIN_pct = tin_pct,
    TIT_pct = tit_pct,
    TIN_class = classify_tin(tin),
    TIT_class = classify_tit(tit)
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
  
  # FGA from basepairs_by_karyotype
  fga <- NA_real_
  
  bp_tbl <- qc[["basepairs_by_karyotype"]]
  if (!is.null(bp_tbl) && nrow(bp_tbl) > 0) {
    if (all(c("minor", "Major", "n", "karyotype") %in% names(bp_tbl))) {
      total_bp <- sum(bp_tbl[["n"]], na.rm = TRUE)
      total_cn <- bp_tbl[["minor"]] + bp_tbl[["Major"]]
      
      altered_bp <- sum(bp_tbl[["n"]][bp_tbl[["karyotype"]] != "1:1"], na.rm = TRUE)
      
      if (total_bp > 0) {
        fga <- altered_bp / total_bp
      }
    }
  }
  
  # combined QC class
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
    fga = fga,
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


### Generic file loader ###

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

### Load TINC and CNAqc tables ###

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

### Merge into unified cohort summary ###

cohort_qc_all <- cohort_qc_cnaqc %>%
  dplyr::full_join(cohort_qc_tinc, by = c("patient_id", "sample_id"))

cohort_qc_all[["sample_id"]] <- sub(".*_(CRC-.*)", "\\1", cohort_qc_all[["sample_id"]])

cohort_qc <- cohort_qc_all %>%
  dplyr::mutate(
    TIN_class = factor(
      TIN_class,
      levels = c("No contamination", "Low contamination", "Contamination", "High contamination")
    ),
    TIT_class = if ("TIT_class" %in% names(.)) TIT_class else classify_tit(TIT),
    TIT_class = factor(
      TIT_class,
      levels = c("Very low purity", "Bad purity", "Average purity", "Good purity", "Very high purity")
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

### Summary strip ###

summary_df <- cohort_qc %>%
  dplyr::summarise(
    n_samples = dplyr::n(),
    n_pass = sum(qc_class == "PASS", na.rm = TRUE),
    n_warn = sum(qc_class == "WARN", na.rm = TRUE),
    n_fail = sum(qc_class == "FAIL", na.rm = TRUE),
    med_purity = median(purity, na.rm = TRUE),
    med_ploidy = median(ploidy, na.rm = TRUE),
    med_fga = median(fga, na.rm = TRUE),
    med_tin = median(TIN_pct, na.rm = TRUE),
    med_tit = median(TIT_pct, na.rm = TRUE)
  )

summary_text <- paste0(
  "Samples: ", summary_df["n_samples"],
  "    |    PASS: ", summary_df[["n_pass"]],
  "    |    WARN: ", summary_df[["n_warn"]],
  "    |    FAIL: ", summary_df[["n_fail"]],
  "\nMedian purity: ", round(summary_df[["med_purity"]], 2),
  "    |    Median ploidy: ", round(summary_df[["med_ploidy"]], 2),
  "    |    Median FGA: ", round(summary_df[["med_fga"]], 2),
  "    |    Median TIN: ", round(summary_df[["med_tin"]], 2), "%",
  "    |    Median TIT: ", round(summary_df[["med_tit"]], 2), "%"
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

  p5 <- ggplot2::ggplot(cohort_qc, ggplot2::aes(fga)) +
  ggplot2::geom_histogram(bins = 30, fill = "#C49C94") +
  theme_dash +
  ggplot2::labs(title = "Fraction genome altered", x = "FGA", y = "Number of samples")

p6 <- ggplot2::ggplot(cohort_qc, ggplot2::aes(mutation_pass_rate, cna_pass_rate)) +
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

cohort_qc <- cohort_qc %>%
  dplyr::mutate(
    purity_cat = cut(
      TIT_pct,
      breaks = c(0, 15, 45, 65, 85, 100),
      labels = c("Very low", "Bad", "Average", "Good", "Very high"),
      include.lowest = TRUE
    )
  )

tit_counts <- cohort_qc %>%
  dplyr::mutate(
    TIT_class = dplyr::case_when(
      TIT_pct < 15 ~ "Very low",
      TIT_pct < 45 ~ "Bad",
      TIT_pct < 65 ~ "Average",
      TIT_pct < 85 ~ "Good",
      TRUE ~ "Very high"
    )
  ) %>%
  dplyr::count(TIT_class) %>%
  dplyr::mutate(perc = n / sum(n))

label_positions <- tibble::tibble(
  TIT_class = c("Very low", "Bad", "Average", "Good", "Very high"),
  x = c(7, 30, 55, 75, 92)
)

tit_counts <- dplyr::left_join(tit_counts, label_positions, by = "TIT_class")


p9 <- ggplot2::ggplot(cohort_qc, ggplot2::aes(TIT_pct, fill = qc_class)) +
  ggplot2::annotate("rect", xmin = 0, xmax = 15, ymin = 0, ymax = Inf, fill = "#d73027", alpha = 0.15) +
  ggplot2::annotate("rect", xmin = 15, xmax = 45, ymin = 0, ymax = Inf, fill = "#fc8d59", alpha = 0.15) +
  ggplot2::annotate("rect", xmin = 45, xmax = 65, ymin = 0, ymax = Inf, fill = "#fee08b", alpha = 0.15) +
  ggplot2::annotate("rect", xmin = 65, xmax = 85, ymin = 0, ymax = Inf, fill = "#91cf60", alpha = 0.15) +
  ggplot2::annotate("rect", xmin = 85, xmax = max(cohort_qc["TIT_pct"], na.rm = TRUE), ymin = 0, ymax = Inf, fill = "#1a9850", alpha = 0.15) +
  ggplot2::geom_histogram(
    bins = 40,
    position = "stack",
    color = "grey20",
    linewidth = 0.2,
    alpha = 0.7
  ) +
  ggplot2::scale_fill_manual(values = qc_colors, name = "QC class", breaks = c("FAIL", "WARN", "PASS")) +
  ggplot2::geom_vline(xintercept = c(15, 45, 65, 85), linetype = "dashed", color = "grey40") +
  ggplot2::annotate("text", x = 7,  y = Inf, label = "Very low", vjust = 2, size = 3.5) +
  ggplot2::annotate("text", x = 30, y = Inf, label = "Bad", vjust = 2, size = 3.5) +
  ggplot2::annotate("text", x = 55, y = Inf, label = "Average", vjust = 2, size = 3.5) +
  ggplot2::annotate("text", x = 75, y = Inf, label = "Good", vjust = 2, size = 3.5) +
  ggplot2::annotate("text", x = 92, y = Inf, label = "Very high", vjust = 2, size = 3.5) +
  ggplot2::geom_text(
    data = tit_counts,
    ggplot2::aes(x = x, y = Inf, label = paste0(n, "\n(", scales::percent(perc, accuracy = 1), ")")),
    inherit.aes = FALSE,
    vjust = 2.0,
    size = 4,
    fontface = "plain"
  ) +
  theme_dash +
  ggplot2::labs(
    title = "Tumour purity estimated by TINC",
    subtitle = "Purity classes: <15%, 15–45%, 45–65%, 65–85%, >85%",
    x = "TIT (%)",
    y = "Number of samples"
  )


p10 <- ggplot2::ggplot(cohort_qc, ggplot2::aes(purity, TIT, color = qc_class)) +
  ggplot2::geom_point(alpha = 0.7, size = 1.6) +
  ggplot2::geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey50") +
  ggplot2::scale_color_manual(values = qc_colors, drop = FALSE) +
  ggplot2::coord_equal(xlim = c(0, 1), ylim = c(0, 1)) +
  theme_dash +
  ggplot2::labs(
    title = "CNAqc purity vs TINC purity",
    x = "CNAqc purity",
    y = "TINC purity",
    color = "QC class"
  )


### QC heatmap ###

qc_heat <- cohort_qc %>%
  dplyr::select(
    sample_id,
    qc_class,
    mutation_pass_rate,
    cna_pass_rate,
    peak_pass_rate,
    mutation_na_fraction,
    cna_na_fraction,
    peak_score,
    mean_delta_vaf,
    purity,
    ploidy,
    fga,
    n_mutations,
    n_cna
  ) %>%
  dplyr::mutate(
    mut_na_good = 1 - mutation_na_fraction,
    cna_na_good = 1 - cna_na_fraction,
    peak_score_good = -peak_score,
    delta_vaf_good = -mean_delta_vaf,
    qc_score =
      0.35 * mutation_pass_rate +
      0.30 * cna_pass_rate +
      0.20 * peak_pass_rate +
      0.10 * mut_na_good +
      0.05 * cna_na_good
  ) %>%
  dplyr::select(
    sample_id,
    qc_class,
    qc_score,
    mutation_pass_rate,
    cna_pass_rate,
    peak_pass_rate,
    mut_na_good,
    cna_na_good,
    peak_score_good,
    delta_vaf_good,
    purity,
    ploidy,
    fga,
    n_mutations,
    n_cna
  ) %>%
  dplyr::arrange(factor(qc_class, levels = c("FAIL", "WARN", "PASS")), qc_score)

qc_heat_sub <- qc_heat %>%
  dplyr::filter(qc_class != "PASS" | qc_score < stats::quantile(qc_score, 0.15, na.rm = TRUE))

get_failure_driver <- function(df) {
  df %>%
    dplyr::mutate(
      failure_driver = dplyr::case_when(
        qc_class == "PASS" ~ "OK",
        mutation_pass_rate < 0.5 ~ "Low mutation QC",
        cna_pass_rate < 0.5 ~ "Low CNA QC",
        peak_pass_rate < 0.5 ~ "Low peak QC",
        mut_na_good < 0.7 ~ "High mutation NA",
        cna_na_good < 0.7 ~ "High CNA NA",
        delta_vaf_good < -0.05 ~ "Poor VAF concordance",
        purity < 0.3 ~ "Low purity",
        TRUE ~ "Mixed / unclear"
      )
    )
}

qc_heat_sub <- get_failure_driver(qc_heat_sub)

qc_heat_sub["qc_class"] <- factor(
  qc_heat_sub["qc_class"],
  levels = c("FAIL", "WARN", "PASS")
)

qc_heat_sub["failure_driver"] <- factor(
  qc_heat_sub["failure_driver"],
  levels = c(
    "Low mutation QC",
    "Low CNA QC",
    "Low peak QC",
    "High mutation NA",
    "High CNA NA",
    "Poor VAF concordance",
    "Low purity",
    "Mixed / unclear",
    "OK"
  )
)

# Matrix for heatmap
mat <- qc_heat_sub %>%
  dplyr::select(
    mutation_pass_rate,
    cna_pass_rate,
    peak_pass_rate,
    mut_na_good,
    cna_na_good,
    peak_score_good,
    delta_vaf_good,
    purity,
    ploidy,
    fga,
    n_mutations,
    n_cna
  ) %>%
  as.matrix()

rownames(mat) <- qc_heat_sub[["sample_id"]]

# z-score by column, then transpose
mat_scaled <- scale(mat)
mat_scaled <- t(mat_scaled)

# Pretty metric labels
pretty_names <- c(
  mutation_pass_rate = "Mutation QC",
  cna_pass_rate = "CNA QC",
  peak_pass_rate = "Peak QC",
  mut_na_good = "Callable mutations",
  cna_na_good = "Callable CNA segments",
  peak_score_good = "Peak fit",
  delta_vaf_good = "VAF concordance",
  purity = "Purity",
  ploidy = "Ploidy",
  fga = "FGA",
  n_mutations = "Mutations",
  n_cna = "CN segments"
)

rownames(mat_scaled) <- pretty_names[rownames(mat_scaled)]
colnames(mat_scaled) <- qc_heat_sub[["sample_id"]]

# Colors
driver_colors <- c(
  "Low mutation QC" = "#d73027",
  "Low CNA QC" = "#fc8d59",
  "Low peak QC" = "#fee08b",
  "High mutation NA" = "#4575b4",
  "High CNA NA" = "#74add1",
  "Poor VAF concordance" = "#984ea3",
  "Low purity" = "#1b9e77",
  "Mixed / unclear" = "grey70",
  "OK" = "#66bd63"
)

col_fun <- colorRamp2(
  c(-2, 0, 2),
  c("#1a476f", "white", "#90353b")
)

# Column annotations
col_ha <- ComplexHeatmap::HeatmapAnnotation(
  `QC class` = qc_heat_sub[["qc_class"]],
  `Failure driver` = qc_heat_sub[["failure_driver"]],
  col = list(
    `QC class` = qc_colors,
    `Failure driver` = driver_colors
  ),
  annotation_name_gp = gpar(fontsize = 10, fontface = "bold")
)

# Heatmap
ht <- Heatmap(
  mat_scaled,
  name = "Z-score",
  col = col_fun,
  top_annotation = col_ha,
  column_split = qc_heat_sub["qc_class"],
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  show_row_names = TRUE,
  row_names_gp = gpar(fontsize = 10),
  show_column_names = FALSE,
  row_title = NULL,
  column_title = "QC: problematic and low-score samples",
  column_title_gp = gpar(fontsize = 14, fontface = "bold"),
  heatmap_legend_param = list(
    title = "Scaled value",
    at = c(-2, -1, 0, 1, 2)
  )
)

ht_grob <- grid.grabExpr(
  draw(ht, heatmap_legend_side = "right", annotation_legend_side = "right")
)

ht_plot <- patchwork::wrap_elements(full = ht_grob)

### Wrap plots ###

first_grid  <- (p1 | p6 | p7)
second_grid <- (p2 | p3 | p4 | p5)
third_grid  <- (p8 | p9 | p10)

dashboard <- summary_plot / first_grid / second_grid / third_grid / ht_plot +
  patchwork::plot_layout(heights = c(0.08, 0.22, 0.20, 0.20, 0.34)) +
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
