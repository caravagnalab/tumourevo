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
  mutation_na_fraction <- if (length(mut_qc) > 0) mut_na / length(mut_qc) else NA_real_
  
  # CNA QC
  cna_tbl <- qc[["cna"]]
  cna_qc <- cna_tbl[["QC_PASS"]]
  
  cna_tested <- sum(cna_qc %in% c(TRUE, FALSE), na.rm = TRUE)
  cna_pass <- sum(cna_qc == TRUE, na.rm = TRUE)
  cna_fail <- sum(cna_qc == FALSE, na.rm = TRUE)
  cna_na <- sum(is.na(cna_qc))
  
  cna_pass_rate <- if (cna_tested > 0) cna_pass / cna_tested else NA_real_
  cna_fail_rate <- if (cna_tested > 0) cna_fail / cna_tested else NA_real_
  cna_na_fraction <- if (length(cna_qc) > 0) cna_na / length(cna_qc) else NA_real_
  
  # peaks analysis
  peak_score <- if (!is.null(qc[["peaks_analysis"]][["score"]])) {
    qc[["peaks_analysis"]][["score"]]
  } else {
    NA_real_
  }
  
  peak_qc <- if (!is.null(qc[["peaks_analysis"]][["QC"]])) {
    qc[["peaks_analysis"]][["QC"]]
  } else {
    NA_character_
  }
  
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
  n_cna_total <- qc[["n_cna"]]
  purity <- qc[["purity"]]
  ploidy <- qc[["ploidy"]]
  most_prev_karyotype <- qc[["most_prevalent_karyotype"]]
  n_karyotype <- qc[["n_karyotype"]]
  
  # mutation-based karyotype dominance
  if (!is.null(n_karyotype) &&
      length(n_karyotype) > 0 &&
      !is.null(most_prev_karyotype) &&
      !is.na(most_prev_karyotype) &&
      most_prev_karyotype %in% names(n_karyotype)) {
    
    n_mutations_prev_karyotype <- as.numeric(
      n_karyotype[[most_prev_karyotype]]
    )
    
    mutation_karyotype_fraction <- if (!is.na(n_mutations) && n_mutations > 0) {
      n_mutations_prev_karyotype / n_mutations
    } else {
      NA_real_
    }
    
  } else {
    n_mutations_prev_karyotype <- NA_real_
    mutation_karyotype_fraction <- NA_real_
  }
  
  # CNA segment / length dominance for the mutation-defined prevalent karyotype
  if (!is.null(cna_tbl) &&
      nrow(cna_tbl) > 0 &&
      all(c("Major", "minor") %in% names(cna_tbl)) &&
      !is.null(most_prev_karyotype) &&
      !is.na(most_prev_karyotype)) {
    
    cna_tbl <- cna_tbl %>%
      dplyr::mutate(
        segment_karyotype = paste0(.data[["Major"]], ":", .data[["minor"]])
      )
    
    n_cna_observed <-  n_cna_total
    
    n_cna_prev_karyotype <- sum(
      cna_tbl[["segment_karyotype"]] == most_prev_karyotype,
      na.rm = TRUE
    )
    
    cna_karyotype_frac <- n_cna_prev_karyotype / n_cna_observed
    
    if ("length" %in% names(cna_tbl)) {
      cna_total_length <- sum(cna_tbl[["length"]], na.rm = TRUE)
      
      cna_length_prev_karyotype <- sum(
        cna_tbl[["length"]][cna_tbl[["segment_karyotype"]] == most_prev_karyotype],
        na.rm = TRUE
      )
      
      cna_length_karyotype_frac <- if (cna_total_length > 0) {
        cna_length_prev_karyotype / cna_total_length
      } else {
        NA_real_
      }
    } else {
      cna_total_length <- NA_real_
      cna_length_prev_karyotype <- NA_real_
      cna_length_karyotype_frac <- NA_real_
    }
    
  } else {
    n_cna_observed <- NA_integer_
    n_cna_prev_karyotype <- NA_integer_
    cna_karyotype_frac <- NA_real_
    cna_total_length <- NA_real_
    cna_length_prev_karyotype <- NA_real_
    cna_length_karyotype_frac <- NA_real_
  }
  
  # QC classification
  qc_class <- dplyr::case_when(
    !is.na(cna_pass_rate) &&
      cna_pass_rate >= 0.80 &&
      (is.na(peak_pass_rate) || peak_pass_rate >= 0.50) ~ "PASS",
    
    !is.na(cna_pass_rate) &&
      (
        cna_pass_rate < 0.50 ||
          (!is.na(peak_pass_rate) && peak_pass_rate < 0.50)
      ) ~ "FAIL",
    
    TRUE ~ "WARN"
  )
  
  tibble::tibble(
    sample_id = sid,
    n_mutations = n_mutations,
    purity = purity,
    ploidy = ploidy,
    n_cna_total = n_cna_total,
    most_prev_karyotype = most_prev_karyotype,
    n_mutations_prev_karyotype = n_mutations_prev_karyotype,
    mutation_karyotype_fraction = mutation_karyotype_fraction,
    n_cna_prev_karyotype = n_cna_prev_karyotype,
    cna_karyotype_frac = cna_karyotype_frac,
    cna_length_prev_karyotype = cna_length_prev_karyotype,
    cna_total_length = cna_total_length,
    cna_length_karyotype_frac = cna_length_karyotype_frac,
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
  PASS = "#79AF97FF",
  WARN = "#DF8F44FF",
  FAIL = "#B24745FF"
)


theme_dash <- theme_bw(base_size = 12) +
  ggplot2::theme(
    plot.title = element_text(size = 14, face = "plain", color = "black"),
    plot.subtitle = element_text(size = 12, face = "plain", color = "black"),
    axis.title = element_text(size = 12, color = "black"),
    axis.text = element_text(size = 10),
    legend.title = element_text(size = 12, color = "black"),
    legend.text = element_text(size = 10, color = "black"),
    strip.text = element_text(size = 14, face = "plain")
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
qc_counts <- cohort_qc %>%
  count(qc_class) %>%
  dplyr::mutate(
    pct = 100 * n / sum(n),
    label = paste0(n, "\n", round(pct, 1), "%")
  )


p1 <- ggplot(qc_counts, aes(x = qc_class, y = n, fill = qc_class)) +
  geom_col(width = 0.7, color = "grey30") +
  geom_text(
    aes(label = label),
    vjust = -0.3,
    size = 4
  ) +
  scale_fill_manual(values = qc_colors,
		    name = "QC class", 
		    drop = FALSE) +
  theme_dash +
  labs(
    title = "QC classification",
    x = "QC class",
    y = "Number of samples"
  )  


# Mutation burden vs most prevalent karyotype

p2 <- cohort_qc %>%
  dplyr::filter(
    !is.na(most_prev_karyotype),
    !is.na(n_mutations),
    n_mutations > 0
  ) %>%
  ggplot(aes(
    x = reorder(most_prev_karyotype, n_mutations, FUN = median),
    y = n_mutations,
    color = qc_class
  )) +
  geom_boxplot(
    aes(group = most_prev_karyotype),
    outlier.shape = NA,
    color = "black",
    fill = "grey90",
    linewidth = 0.5
  ) +
  geom_jitter(
    width = 0.2,
    alpha = 0.5,
    size = 1.4
  ) +
  scale_y_log10(labels = scales::comma) +
  scale_color_manual(values = qc_colors,
		     name = "QC class",
		     drop = FALSE) +
  coord_flip() +
  theme_dash +
  labs(
    title = "Mutation burden by most prevalent karyotype",
    x = "Most prevalent karyotype",
    y = "Number of mutations per sample"
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
    fill = "#925E9FB2",
    color = "grey20",
    linewidth = 0.2
  ) +
  theme_dash +
  labs(
    title = "Ploidy",
    x = "Ploidy",
    y = "Number of samples"
  )

# CNA segment number vs most prevalent karyotype

p5 <- cohort_qc %>%
  filter(
    !is.na(most_prev_karyotype),
    !is.na(n_cna_total),
    n_cna_total > 0
  ) %>%
  ggplot(aes(
    x = reorder(most_prev_karyotype, n_cna_total, FUN = median),
    y = n_cna_total,
    color = qc_class
  )) +
  geom_boxplot(
    aes(group = most_prev_karyotype),
    outlier.shape = NA,
    color = "black",
    fill = "grey90",
    linewidth = 0.5
  ) +
  geom_jitter(
    width = 0.2,
    alpha = 0.5,
    size = 1.4
  ) +
  scale_y_log10(labels = scales::comma) +
  scale_color_manual(values = qc_colors,
		     name = "QC class",
		     drop = FALSE) +
  coord_flip() +
  theme_dash +
  labs(
    title = "CNA segments by most prevalent karyotype",
    x = "Most prevalent karyotype",
    y = "Number of CNA segments"
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

top_n_karyotypes <- 15

karyotype_freq <- cohort_qc %>%
  dplyr::filter(!is.na(most_prev_karyotype)) %>%
  count(most_prev_karyotype, name = "n_samples") %>%
  dplyr::mutate(
    pct_samples = 100 * n_samples / sum(n_samples)
  ) %>%
  dplyr::arrange(desc(pct_samples))

top_karyotype_freq <- karyotype_freq %>%
  slice_head(n = top_n_karyotypes)

p7 <- ggplot(top_karyotype_freq, aes(
  x = reorder(most_prev_karyotype, pct_samples),
  y = pct_samples
)) +
  geom_col(fill = "#4C78A8") +
  geom_text(
    aes(label = paste0(round(pct_samples, 1), "%")),
    hjust = -0.1,
    size = 3
  ) +
  coord_flip() +
  expand_limits(y = max(top_karyotype_freq[["pct_samples"]], na.rm = TRUE) * 1.15) +
  theme_dash +
  labs(
    title = "Most prevalent karyotypes",
    #subtitle = paste0("Top ", top_n_karyotypes, " karyotypes by percentage of samples"),
    x = "Most prevalent karyotype",
    y = "Samples (%)"
  )


# TIN distribution
p8 <- ggplot(cohort_qc, aes(TIN_pct)) +
  geom_histogram(bins = 40, fill = "#D95F02", color = "white") +
  
  geom_vline(xintercept = c(1, 7, 15),
             linetype = "dashed", color = "grey40") +
  
  # labels
  annotate("text", x = 0.5, y = Inf, label = "No",
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
second_grid <- (p7 | p2 | p5)
third_grid  <- (p4 | p3 | p8)

dashboard <- summary_plot / first_grid / second_grid / third_grid +
  patchwork::plot_layout(heights = c(0.10, 0.30, 0.28, 0.28)) +
  patchwork::plot_annotation(
    title = paste0(prefix, " cohort QC summary"),
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
