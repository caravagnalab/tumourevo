#!/usr/bin/env Rscript

library(tidyverse)
library(ggplot2)
library(RColorBrewer)


parse_args = function(x) {
  x = gsub("\\\\[","",x)
  x = gsub("\\\\]","",x)
  # giving errors when we have lists like c(xxx, xxx) since it will separate it
  # args_list = unlist(strsplit(x, ', ')[[1]])
  args_list = unlist(strsplit(x, ", (?=[^)]*(?:\\\\(|\$))", perl=TRUE))
  args_vals = lapply(args_list, function(x) {
    x_splt = strsplit(x, split=":")[[1]]
    c(x_splt[1],  paste(x_splt[2:length(x_splt)], collapse=":"))
  })

  # Ensure the option vectors are length 2 (key/ value) to catch empty ones
  args_vals = lapply(args_vals, function(z){ length(z) = 2; z})

  parsed_args = structure(lapply(args_vals, function(x) x[2]), names = lapply(args_vals, function(x) x[1]))
  parsed_args[! is.na(parsed_args)]
}

opt = list(
  prefix = ifelse('$task.ext.prefix' == 'null', '$meta.id', '$task.ext.prefix')
)
args_opt = parse_args('$task.ext.args')
for ( ao in names(args_opt)) opt[[ao]] = args_opt[[ao]]

# Scoring parameters, see "Genome interpreter: fit QC, data quality and revised
# expansion score" in docs/. Any of them can be overridden through task.ext.args.
defaults = list(
  tools = "",
  w_driver = 1,
  w_signature = 1,
  w_tail = 1,
  missing_value = 0.5,
  tier_thresholds = "0.9,0.55,0.2",
  clonal_min_ccf = 0.9,
  truncal_drivers = "APC,KRAS,TP53,SMAD4",
  read_support = 6,
  min_cluster_size = 100,
  min_cluster_frac = 0.02,
  qmode_bootstrap = 200,
  qmode_dip = 0.1,
  qmode_scale = "ccf",
  qfit_tv_low = 0.05,
  qfit_tv_high = 0.15,
  qfit_bin_width = 0.02,
  qfit_noise_correction = TRUE,
  qfit_nonconverged_cap = 0.5,
  tail_min_mutations = 1000,
  tail_min_clonal_ratio = 3,
  tail_shape_bin_width = 0.01,
  tail_shape_min_bins = 3,
  tail_shape_z = 3,
  tail_shape_min_frac = 0.02,
  driver_chance = TRUE,
  signature_null = "refit",
  signature_null_draws = 200,
  signature_min_mutations = 100,
  artefact_signatures = "SBS27,SBS43,SBS45-SBS60,SBS95",
  seed = 42
)
for (d in names(defaults)) {
  if (is.null(opt[[d]]) || is.na(opt[[d]]) || opt[[d]] %in% c("", "null", "NA")) opt[[d]] = defaults[[d]]
}

num = function(k) as.numeric(opt[[k]])
flag = function(k) toupper(as.character(opt[[k]])) %in% c("TRUE", "T", "1", "YES")
split_list = function(x) {
  x = trimws(unlist(strsplit(as.character(x), "[,;]")))
  x[x != ""]
}
split_paths = function(x) {
  x = unlist(strsplit(trimws(x), " +"))
  x[x != ""]
}

# "SBS45-SBS60" -> SBS45, SBS46, ..., SBS60
expand_signatures = function(x) {
  unlist(lapply(split_list(x), function(s) {
    if (!grepl("-", s, fixed=TRUE)) return(s)
    ends = strsplit(s, "-", fixed=TRUE)[[1]]
    prefix = gsub("[0-9]", "", ends[1])
    paste0(prefix, seq(as.integer(gsub("[^0-9]", "", ends[1])), as.integer(gsub("[^0-9]", "", ends[2]))))
  }))
}

set.seed(num("seed"))
tiers = sort(as.numeric(split_list(opt[["tier_thresholds"]])), decreasing=TRUE)
artefacts = expand_signatures(opt[["artefact_signatures"]])
truncal_genes = split_list(opt[["truncal_drivers"]])
missing_value = num("missing_value")
a_star = num("read_support")

mutation_tables = split_paths("$mutation_tables")
results_sigprofiler = split_paths("$results_sigprofiler")
fit_files = split_paths("$fits")
joint_files = split_paths("$joint_table")

patient_id = opt[["prefix"]]

mut_key = function(chr, pos, ref, alt) {
  paste(sub("^chr", "", as.character(chr)), as.integer(pos), ref, alt, sep=":")
}

col_or = function(df, name, default) {
  if (name %in% colnames(df)) df[[name]] else rep(default, nrow(df))
}

tier_of = function(x) {
  case_when(is.na(x) ~ NA_character_,
            x >= tiers[1] ~ "Tier 1",
            x >= tiers[2] ~ "Tier 2",
            x >= tiers[3] ~ "Tier 3",
            TRUE ~ "Tier 4")
}

# Model-free density helpers (Layer 1 and q_mode) ######

kde = function(x) {
  bw = tryCatch(bw.SJ(x), error=function(e) bw.nrd0(x))
  density(x, bw=bw, n=512)
}

# Local maxima of the KDE; a maximum counts as a mode only if the dip to its
# neighbour is at least `dip` of the lower of the two peaks.
find_modes = function(x, dip=num("qmode_dip")) {
  x = x[is.finite(x)]
  if (length(x) < 5 || length(unique(x)) < 3) return(tibble(x=numeric(), h=numeric()))
  d = kde(x)
  y = d[["y"]]
  gx = d[["x"]]
  n = length(y)
  idx = which(y[2:(n-1)] > y[1:(n-2)] & y[2:(n-1)] >= y[3:n]) + 1
  repeat {
    if (length(idx) < 2) break
    drop = NA
    for (i in seq_len(length(idx) - 1)) {
      low = min(y[idx[i]], y[idx[i+1]])
      valley = min(y[idx[i]:idx[i+1]])
      if (low - valley < dip * low) {
        drop = ifelse(y[idx[i]] < y[idx[i+1]], i, i + 1)
        break
      }
    }
    if (is.na(drop)) break
    idx = idx[-drop]
  }
  tibble(x=gx[idx], h=y[idx])
}

# Clonal 1:1 VAF peak: the highest-VAF mode holding at least a quarter of the
# tallest peak's density.
clonal_peak = function(vaf) {
  vaf = vaf[is.finite(vaf) & vaf > 0]
  if (length(vaf) < 50) return(NA_real_)
  m = find_modes(vaf)
  if (nrow(m) == 0) return(NA_real_)
  max(m[["x"]][m[["h"]] >= 0.25 * max(m[["h"]])])
}

# q_mode: bootstrap share of resamples in which a KDE mode, assigned to the
# nearest cluster centre, falls inside the cluster's interquartile range.
mode_support = function(df) {
  centres = df %>%
    group_by(cluster) %>%
    summarise(centre=median(x), q1=quantile(x, .25), q3=quantile(x, .75), n=n(), .groups="drop")
  own_peak = function(xs) {
    m = find_modes(xs)[["x"]]
    if (length(m) == 0) return(rep(FALSE, nrow(centres)))
    nearest = sapply(m, function(v) which.min(abs(centres[["centre"]] - v)))
    sapply(seq_len(nrow(centres)), function(j) {
      any(nearest == j & m >= centres[["q1"]][j] & m <= centres[["q3"]][j])
    })
  }
  hits = matrix(replicate(num("qmode_bootstrap"), own_peak(sample(df[["x"]], replace=TRUE))),
                nrow=nrow(centres))
  centres %>%
    mutate(q_mode=ifelse(n >= 3, rowMeans(hits), NA_real_)) %>%
    select(cluster, q_mode)
}

# q_fit: total-variation distance between the observed VAF histogram and the
# VAF histogram predicted by the fit at the observed depths.
fit_tv = function(DP, NV, p) {
  ok = is.finite(p) & is.finite(DP) & is.finite(NV) & DP > 0
  DP = DP[ok]
  NV = NV[ok]
  p = pmin(pmax(p[ok], 1e-9), 1 - 1e-9)
  if (length(DP) < 50) return(c(tv=NA_real_, tv_noise=NA_real_))
  w = num("qfit_bin_width")
  lowers = seq(0, 1 - w / 2, by=w)
  K = length(lowers)
  cdf = sapply(lowers[-1], function(b) mean(pbinom(ceiling(b * DP - 1e-9) - 1, DP, p)))
  expected = diff(c(0, cdf, 1))
  hist_of = function(nv) tabulate(findInterval(nv / DP, lowers), nbins=K) / length(DP)
  tv = 0.5 * sum(abs(hist_of(NV) - expected))
  tv_noise = mean(replicate(10, 0.5 * sum(abs(hist_of(rbinom(length(DP), DP, p)) - expected))))
  c(tv=tv, tv_noise=tv_noise)
}

fit_converged = function(fit) {
  if (!is.list(fit)) return(NA)
  for (k in c("converged", "convergence")) {
    if (!is.null(fit[[k]])) return(isTRUE(as.logical(fit[[k]][1])))
  }
  if (!is.null(fit[["status"]])) {
    status = paste(fit[["status"]], collapse=" ")
    return(grepl("converge", status, ignore.case=TRUE) && !grepl("not", status, ignore.case=TRUE))
  }
  NA
}

# Joint table: per-mutation, per-sample depth, VAF, karyotype, purity ######
# Purity and copy number are validated by CNAqc upstream and are not re-scored here.

joint = NULL
if (length(joint_files) > 0) {
  joint_raw = lapply(joint_files, function(f) readr::read_tsv(f, show_col_types=FALSE, guess_max=1e5)) %>%
    bind_rows()
  sample_col = if ("Indiv" %in% colnames(joint_raw)) "Indiv" else "sample_id"
  joint = tibble(
    mut_key = mut_key(joint_raw[["chr"]], joint_raw[["from"]], joint_raw[["ref"]], joint_raw[["alt"]]),
    sample = as.character(col_or(joint_raw, sample_col, "sample")),
    DP = as.numeric(joint_raw[["DP"]]),
    NV = as.numeric(joint_raw[["NV"]]),
    VAF = as.numeric(joint_raw[["VAF"]]),
    karyotype = as.character(joint_raw[["karyotype"]]),
    purity = as.numeric(joint_raw[["purity"]]),
    CCF_cnaqc = as.numeric(col_or(joint_raw, "CCF", NA)),
    multiplicity = as.numeric(col_or(joint_raw, "mutation_multiplicity", NA))
  ) %>%
    distinct(mut_key, sample, .keep_all=TRUE) %>%
    mutate(CN = as.numeric(sub(":.*", "", karyotype)) + as.numeric(sub(".*:", "", karyotype)),
           multiplicity = ifelse(is.na(multiplicity), 1, multiplicity),
           CCF = ifelse(is.na(CCF_cnaqc),
                        VAF * (purity * CN + 2 * (1 - purity)) / (purity * multiplicity),
                        CCF_cnaqc))
}

# Clonal 1:1 VAF peak per sample, used by the tail-estimability check
sample_peaks = tibble(sample=character(), purity=numeric(), v_clonal=numeric())
if (!is.null(joint)) {
  sample_peaks = joint %>%
    group_by(sample) %>%
    summarise(purity=median(purity, na.rm=TRUE),
              v_clonal=clonal_peak(VAF[karyotype == "1:1"]),
              .groups="drop")
}

# MOBSTER: tail labels, tail estimability and tail shape (T5) per sample ######

mobster_fit_files = fit_files[endsWith(basename(fit_files), "mobster_best_fit.rds")]
mobster_files = grep(mutation_tables, pattern="mobster", value=T)

tail_shape = function(fit, lab) {
  empty = tibble(vaf_from=numeric(), vaf_to=numeric(), n_tail=integer())
  tail_vaf = lab[["VAF"]][lab[["label"]] == "Tail"]
  shape = suppressWarnings(as.numeric(fit[["shape"]])[1])
  scale = suppressWarnings(as.numeric(fit[["scale"]])[1])
  if (length(tail_vaf) == 0 || !is.finite(shape) || !is.finite(scale)) return(empty)
  pi_tail = fit[["pi"]]
  pi_tail = if (!is.null(pi_tail) && "Tail" %in% names(pi_tail)) pi_tail[["Tail"]] else length(tail_vaf) / nrow(lab)
  w = num("tail_shape_bin_width")
  lo = seq(scale, 1 - w, by=w)
  if (length(lo) < 1) return(empty)
  hi = lo + w
  pareto_cdf = function(x) ifelse(x <= scale, 0, 1 - (scale / x)^shape)
  expected = nrow(lab) * pi_tail * (pareto_cdf(hi) - pareto_cdf(lo))
  observed = sapply(seq_along(lo), function(i) sum(tail_vaf >= lo[i] & tail_vaf < hi[i]))
  z = (observed - expected) / sqrt(pmax(expected, 1e-9))
  runs = rle(z > num("tail_shape_z"))
  run_end = cumsum(runs[["lengths"]])
  run_start = run_end - runs[["lengths"]] + 1
  keep = which(runs[["values"]] & runs[["lengths"]] >= num("tail_shape_min_bins"))
  if (length(keep) == 0) return(empty)
  tibble(vaf_from=lo[run_start[keep]], vaf_to=hi[run_end[keep]],
         n_tail=sapply(keep, function(k) sum(observed[run_start[k]:run_end[k]]))) %>%
    filter(n_tail >= num("tail_shape_min_frac") * length(tail_vaf))
}

mobster_labels = tibble(mut_key=character(), sample=character(), label=character(), VAF=numeric())
tail_samples = tibble(sample=character(), tail_estimable=logical(), tail_reason=character())
tail_modes = tibble(sample=character(), vaf_from=numeric(), vaf_to=numeric(), n_tail=integer())

if (length(mobster_fit_files) > 0) {
  for (f in mobster_fit_files) {
    fit = readRDS(f)
    d = fit[["data"]]
    s = if ("sample_id" %in% colnames(d)) as.character(d[["sample_id"]][1]) else sub("_mobster_best_fit.rds", "", basename(f), fixed=TRUE)
    lab = tibble(mut_key=mut_key(d[["chr"]], d[["from"]], d[["ref"]], d[["alt"]]),
                 sample=s, label=as.character(d[["cluster"]]), VAF=as.numeric(d[["VAF"]]))
    mobster_labels = bind_rows(mobster_labels, lab)

    v_min = min(lab[["VAF"]], na.rm=TRUE)
    v_cl = sample_peaks %>% filter(sample == s) %>% pull(v_clonal)
    if (length(v_cl) == 0 || is.na(v_cl[1])) {
      beta_means = as.numeric(fit[["a"]]) / (as.numeric(fit[["a"]]) + as.numeric(fit[["b"]]))
      v_cl = if (length(beta_means) > 0) max(beta_means, na.rm=TRUE) else NA_real_
    }
    v_cl = v_cl[1]
    has_tail = !isFALSE(fit[["fit.tail"]]) && any(lab[["label"]] == "Tail")
    reason = case_when(
      !has_tail ~ "no tail fitted",
      nrow(lab) < num("tail_min_mutations") ~ paste0("<", num("tail_min_mutations"), " mutations"),
      is.na(v_cl) || v_cl < num("tail_min_clonal_ratio") * v_min ~ "clonal peak too close to VAF floor",
      TRUE ~ "")
    tail_samples = bind_rows(tail_samples, tibble(sample=s, tail_estimable=(reason == ""), tail_reason=reason))
    if (has_tail) tail_modes = bind_rows(tail_modes, tail_shape(fit, lab) %>% mutate(sample=s))
  }
} else if (length(mobster_files) > 0) {
  # No fit objects: fall back to the per-sample signature tables, without the
  # estimability and tail-shape checks.
  for (f in mobster_files) {
    tb = readr::read_tsv(f, show_col_types=FALSE, col_types=readr::cols(chrom="c", Sample="c", driver_label="c"))
    lab = tibble(mut_key=mut_key(tb[["chrom"]], tb[["pos_start"]], tb[["ref"]], tb[["alt"]]),
                 sample=as.character(tb[["Project"]]), label=as.character(tb[["Sample"]]), VAF=NA_real_)
    mobster_labels = bind_rows(mobster_labels, lab)
    tail_samples = bind_rows(tail_samples,
                             tibble(sample=unique(lab[["sample"]]),
                                    tail_estimable=any(lab[["label"]] == "Tail"),
                                    tail_reason=ifelse(any(lab[["label"]] == "Tail"), "", "no tail fitted")))
  }
}
run_mobster = nrow(mobster_labels) > 0
w_tail = ifelse(run_mobster, num("w_tail"), 0)

# Signature comparison helpers ######

cos_dist = function(a, b) {
  d = sqrt(sum(a^2)) * sqrt(sum(b^2))
  if (d == 0) 1 else 1 - sum(a * b) / d
}

clean_exposure = function(h) {
  h = h[!names(h) %in% artefacts]
  if (sum(h) > 0) h / sum(h) else h
}

# Divergence D_C: cosine dissimilarity to the clonal exposures, share of
# differing signatures (exposure > 5%) and dissimilarity to the subclonal background.
divergence = function(e, e0, bg, thr=.05) {
  s1 = names(e)[e > thr]
  s0 = names(e0)[e0 > thr]
  u = union(s1, s0)
  n_rel = if (length(u) > 0) length(setdiff(u, intersect(s1, s0))) / length(u) else 0
  parts = c(cos_dist(e, e0), n_rel)
  if (!is.null(bg)) parts = c(parts, cos_dist(e, bg))
  mean(parts)
}

# Exposure refit with fixed signatures (KL-NMF multiplicative updates)
refit = function(v, W, iters=500, tol=1e-6) {
  h = rep(sum(v) / ncol(W), ncol(W))
  ws = colSums(W)
  for (i in seq_len(iters)) {
    wh = pmax(as.vector(W %*% h), 1e-12)
    h_new = h * as.vector(crossprod(W, v / wh)) / ws
    converged = max(abs(h_new - h)) < tol * max(sum(v), 1)
    h = h_new
    if (converged) break
  }
  setNames(h, colnames(W))
}

# Multinomial resample of n mutations from a spectrum
resample_spectrum = function(counts, n) {
  as.vector(rmultinom(1, max(round(n), 1), counts / sum(counts)))
}

read_signature_inputs = function(tool, sign_type) {
  roots = grep(results_sigprofiler, pattern=tool, value=T)
  if (length(roots) == 0) return(NULL)
  dirs = list.dirs(roots, recursive=T)
  act_dir = dirs[endsWith(dirs, paste0(sign_type, "/Assignment_Solution/Activities"))]
  if (length(act_dir) == 0) {
    act_dir = dirs %>%
      keep(function(x) grepl(pattern=sign_type, x=x)) %>%
      keep(function(x) endsWith(x, "/Activities"))
  }
  if (length(act_dir) == 0) return(NULL)
  act_dir = act_dir[1]

  activities = read.table(file=file.path(act_dir, "Assignment_Solution_Activities.txt"), header=T, check.names=F)
  colnames(activities)[1] = "cluster_tool"
  activities[["cluster_tool"]] = as.character(activities[["cluster_tool"]])

  matrix_ext = ifelse(sign_type == "SBS", ".SBS96.all", ".ID83.all")
  files = list.files(roots, recursive=T, full.names=T)
  matrix_file = files[endsWith(files, matrix_ext)]
  profile_file = file.path(dirname(act_dir), "Signatures", "Assignment_Solution_Signatures.txt")
  spectra = profiles = NULL
  if (length(matrix_file) > 0 && file.exists(profile_file)) {
    spectra = read.table(matrix_file[1], header=T, check.names=F, sep="", comment.char="")
    rownames(spectra) = spectra[[1]]
    spectra = as.matrix(spectra[, -1, drop=F])
    profiles = read.table(profile_file, header=T, check.names=F, sep="", comment.char="")
    rownames(profiles) = profiles[[1]]
    profiles = as.matrix(profiles[, -1, drop=F])
    shared = intersect(rownames(spectra), rownames(profiles))
    spectra = spectra[shared, , drop=F]
    profiles = profiles[shared, , drop=F]
  }
  list(activities=activities, spectra=spectra, profiles=profiles)
}

# phi_signature: share of null draws whose divergence is below the observed one.
# Each draw compares two independent resamples of the clonal spectrum, one of the
# cluster's size and one of the clonal cluster's size, so the null carries the
# sampling noise of both sides of the comparison.
signature_scores = function(sig, clonal, clusters) {
  act = sig[["activities"]]
  sig_cols = setdiff(colnames(act), "cluster_tool")
  activity_of = function(cl) {
    row = act[act[["cluster_tool"]] == cl, sig_cols, drop=F]
    if (nrow(row) == 0) return(setNames(rep(0, length(sig_cols)), sig_cols))
    unlist(row[1, ])
  }
  mode = opt[["signature_null"]]
  spectra = sig[["spectra"]]
  have_spectra = !is.null(spectra) && all(c(clonal, clusters) %in% colnames(spectra))
  if (!have_spectra) mode = "exposure"

  subclonal = setdiff(clusters, clonal)
  if (mode == "refit") {
    used = sig_cols[colSums(act[, sig_cols, drop=F]) > 0]
    used = intersect(used, colnames(sig[["profiles"]]))
    if (length(used) == 0) mode = "exposure"
    W = sig[["profiles"]][, used, drop=F]
    W = sweep(W, 2, pmax(colSums(W), 1e-12), "/")
    exposure_of = function(v) clean_exposure(refit(v, W))
  }

  e0 = bg = NULL
  if (mode == "refit") {
    e0 = exposure_of(spectra[, clonal])
    if (length(subclonal) > 0) bg = exposure_of(rowSums(spectra[, subclonal, drop=F]))
  } else if (mode == "exposure") {
    e0 = clean_exposure(activity_of(clonal))
    if (length(subclonal) > 0) bg = clean_exposure(Reduce(`+`, lapply(subclonal, activity_of)))
  }

  clonal_total = if (have_spectra) sum(spectra[, clonal]) else sum(activity_of(clonal))
  if (clonal_total == 0) return(tibble())

  lapply(subclonal, function(cl) {
    n_c = if (have_spectra) sum(spectra[, cl]) else sum(activity_of(cl))
    B = num("signature_null_draws")
    if (mode == "spectrum") {
      v0 = spectra[, clonal]
      d_obs = cos_dist(spectra[, cl], v0)
      d_null = replicate(B, cos_dist(resample_spectrum(v0, n_c), resample_spectrum(v0, clonal_total)))
    } else if (mode == "refit") {
      v0 = spectra[, clonal]
      d_obs = divergence(exposure_of(spectra[, cl]), e0, bg)
      d_null = replicate(B, divergence(exposure_of(resample_spectrum(v0, n_c)),
                                       exposure_of(resample_spectrum(v0, clonal_total)), bg))
    } else {
      d_obs = divergence(clean_exposure(activity_of(cl)), e0, bg)
      p0 = activity_of(clonal)
      d_null = replicate(B, {
        draw = setNames(resample_spectrum(p0, n_c), names(p0))
        ref = setNames(resample_spectrum(p0, clonal_total), names(p0))
        divergence(clean_exposure(draw), clean_exposure(ref), bg)
      })
    }
    tibble(cluster_tool=cl, n_signature=n_c, divergence=d_obs,
           phi_sig=mean(d_null < d_obs), signature_null=mode)
  }) %>% bind_rows()
}

# Scoring per tool ######

tool_list <- strsplit(opt[['tools']], ",")[[1]]
tool_list <- tool_list[tool_list %in% c("viber", "pyclone-vi")]
tool_list <- gsub("pyclone-vi", "pyclonevi", tool_list)

table_signatures <- tibble()
table_driver <- tibble()
cluster_sample_table <- tibble()
sample_quality <- tibble()

score_table = lapply(tool_list, function(tool) {
  mutations_file = grep(mutation_tables, pattern=tool, value=T)
  if (length(mutations_file) == 0) return(tibble())

  mutations_tool = readr::read_tsv(mutations_file, show_col_types=FALSE, guess_max=1e6,
                                   col_types=readr::cols(chrom="c", driver_label="c", Sample="c")) %>%
    mutate(patient_id=patient_id, tool=tool, chrom=as.character(chrom),
           mut_key=mut_key(chrom, pos_start, ref, alt),
           mutation_id=paste(chrom, pos_start, pos_end, ref, alt, sep=":")) %>%
    select(Project, patient_id, mutation_id, mut_key, everything(), -chrom, -pos_start,
           -pos_end, -ref, -alt, -Type, -ID, -Genome, -mut_type) %>%
    rename(cluster_tool=Sample) %>%
    mutate(cluster_tool=as.character(cluster_tool))

  table_driver <<- bind_rows(table_driver, mutations_tool %>% filter(is_driver) %>% select(-mut_key))

  clusters = mutations_tool %>%
    group_by(cluster_tool) %>%
    summarise(n=n_distinct(mut_key), is_clonal_tool=any(is_clonal %in% TRUE), .groups="drop")
  n_fit = n_distinct(mutations_tool[["mut_key"]])
  n_min = max(num("min_cluster_size"), num("min_cluster_frac") * n_fit)
  clusters = clusters %>% mutate(q_size=pmin(1, n / n_min))

  # Per-mutation, per-sample data for the mutations this tool clustered
  muts = NULL
  if (!is.null(joint)) {
    muts = mutations_tool %>%
      distinct(mut_key, cluster_tool) %>%
      inner_join(joint, by="mut_key", relationship="many-to-many")
    if (nrow(muts) == 0) muts = NULL
  }

  # Clonal anchor (T3): largest cluster with CCF >= clonal_min_ccf in every sample
  cluster_ccf = tibble(cluster_tool=character(), sample=character(), ccf=numeric())
  anchor = clusters %>% filter(is_clonal_tool) %>% slice_max(n, n=1, with_ties=FALSE) %>% pull(cluster_tool)
  anchor_note = "tool label"
  if (!is.null(muts)) {
    cluster_ccf = muts %>%
      group_by(cluster_tool, sample) %>%
      summarise(ccf=median(CCF, na.rm=TRUE), .groups="drop")
    candidates = cluster_ccf %>%
      group_by(cluster_tool) %>%
      summarise(min_ccf=min(ccf, na.rm=TRUE), .groups="drop") %>%
      filter(min_ccf >= num("clonal_min_ccf")) %>%
      left_join(clusters, by="cluster_tool") %>%
      arrange(desc(n))
    if (nrow(candidates) > 0) {
      if (length(anchor) == 0 || candidates[["cluster_tool"]][1] != anchor[1]) {
        anchor_note = paste0("moved from ", ifelse(length(anchor) == 0, "none", anchor[1]))
      } else {
        anchor_note = "confirmed"
      }
      anchor = candidates[["cluster_tool"]][1]
    } else {
      anchor_note = paste0("tool label (no cluster with CCF >= ", num("clonal_min_ccf"), ")")
    }
  }
  clusters = clusters %>% mutate(is_clonal=(cluster_tool %in% anchor))

  # Fitted VAF per mutation and sample, for q_fit
  fitted_vaf = NULL
  converged = NA
  if (!is.null(muts)) {
    if (tool == "viber") {
      viber_file = fit_files[grepl("viber", basename(fit_files))]
      if (length(viber_file) > 0) {
        fit = readRDS(viber_file[1])
        converged = fit_converged(fit)
        theta = fit[["theta_k"]]
        if (!is.null(theta)) {
          theta_long = as.data.frame(as.table(as.matrix(theta)), stringsAsFactors=FALSE)
          colnames(theta_long) = c("sample", "cluster_tool", "p")
          fitted_vaf = muts %>%
            left_join(theta_long, by=c("sample", "cluster_tool")) %>%
            select(mut_key, sample, p)
        }
      }
    } else if (tool == "pyclonevi") {
      pyclone_file = fit_files[grepl("best_fit.txt", basename(fit_files), fixed=TRUE)]
      if (length(pyclone_file) > 0) {
        prevalence = readr::read_tsv(pyclone_file[1], show_col_types=FALSE) %>%
          transmute(sample=as.character(sample_id), cluster_tool=paste0("C", cluster_id),
                    prevalence=cellular_prevalence) %>%
          distinct(sample, cluster_tool, .keep_all=TRUE)
        fitted_vaf = muts %>%
          left_join(prevalence, by=c("sample", "cluster_tool")) %>%
          mutate(p=prevalence * purity * multiplicity / (purity * CN + 2 * (1 - purity))) %>%
          select(mut_key, sample, p)
      }
    }
  }

  # Layer 1 and Layer 2, per sample and per cluster in each sample
  cs = tibble(cluster_tool=character(), sample=character())
  if (!is.null(muts)) {
    samples = unique(muts[["sample"]])
    sq = lapply(samples, function(s) {
      m = muts %>% filter(sample == s)
      eff_depth = median(m[["purity"]], na.rm=TRUE) * median(m[["DP"]][m[["DP"]] > 0], na.rm=TRUE)
      tvs = c(tv=NA_real_, tv_noise=NA_real_)
      if (!is.null(fitted_vaf)) {
        fm = m %>% left_join(fitted_vaf %>% filter(sample == s) %>% select(-sample), by="mut_key")
        tvs = fit_tv(fm[["DP"]], fm[["NV"]], fm[["p"]])
      }
      tv_used = if (flag("qfit_noise_correction")) max(0, tvs[["tv"]] - tvs[["tv_noise"]]) else tvs[["tv"]]
      q_fit = pmin(1, pmax(0, (num("qfit_tv_high") - tv_used) / (num("qfit_tv_high") - num("qfit_tv_low"))))
      if (isFALSE(converged)) q_fit = min(q_fit, num("qfit_nonconverged_cap"), na.rm=TRUE)
      tibble(tool=tool, sample=s, n_mutations=n_distinct(m[["mut_key"]]),
             effective_depth=eff_depth, ccf_star=2 * a_star / eff_depth,
             tv=tvs[["tv"]], tv_noise=tvs[["tv_noise"]], fit_converged=converged, q_fit=q_fit)
    }) %>% bind_rows() %>%
      left_join(sample_peaks, by="sample")
    sample_quality <<- bind_rows(sample_quality, sq %>% mutate(patient_id=patient_id))

    support = muts %>%
      left_join(cluster_ccf, by=c("cluster_tool", "sample")) %>%
      mutate(a=purity * DP * ccf * multiplicity / (purity * CN + 2 * (1 - purity))) %>%
      group_by(cluster_tool, sample) %>%
      summarise(ccf=first(ccf), reads=median(a, na.rm=TRUE), .groups="drop") %>%
      mutate(q_support=pmin(1, reads / a_star))

    modes = lapply(samples, function(s) {
      m = muts %>% filter(sample == s, VAF > 0)
      if (opt[["qmode_scale"]] == "vaf_1_1") {
        m = m %>% filter(karyotype == "1:1") %>% mutate(x=VAF)
      } else {
        m = m %>% mutate(x=pmin(CCF, 1.5))
      }
      m = m %>% filter(is.finite(x)) %>% transmute(cluster=cluster_tool, x)
      if (nrow(m) < 5) return(tibble())
      mode_support(m) %>% transmute(cluster_tool=cluster, sample=s, q_mode)
    }) %>% bind_rows()

    cs = support %>%
      left_join(modes, by=c("cluster_tool", "sample")) %>%
      left_join(sq %>% select(sample, q_fit), by="sample")
  }

  factor_names = c("q_size", "q_support", "q_mode", "q_fit")
  cs = clusters %>%
    select(cluster_tool, q_size) %>%
    left_join(cs, by="cluster_tool") %>%
    mutate(tool=tool)
  for (f in setdiff(c(factor_names, "ccf", "reads"), colnames(cs))) cs[[f]] = NA_real_
  cs = cs %>%
    rowwise() %>%
    mutate(R=min(c_across(all_of(factor_names)), na.rm=TRUE),
           limiting_factor={
             v = c_across(all_of(factor_names))
             if (all(is.na(v)) || min(v, na.rm=TRUE) >= 1) "-" else factor_names[which.min(v)]
           }) %>%
    ungroup() %>%
    mutate(limiting_factor=ifelse(limiting_factor == "q_support",
                                  sprintf("q_support (~%.1f reads)", reads), limiting_factor))
  cluster_sample_table <<- bind_rows(cluster_sample_table, cs %>% mutate(patient_id=patient_id))

  # A cluster is judged in the sample where it is best resolved
  reliability = cs %>%
    group_by(cluster_tool) %>%
    slice_max(R, n=1, with_ties=FALSE) %>%
    ungroup() %>%
    select(cluster_tool, best_sample=sample, ccf, all_of(factor_names), R, limiting_factor)

  # phi_tail: share of mutations never labelled Tail, in samples where the tail is
  # estimable and has no extra mode at the cluster's VAF (T5)
  tail_table = clusters %>% select(cluster_tool) %>% mutate(phi_tail=NA_real_, tail_note="")
  driver_tail = tibble(mut_key=character(), tail_driver=logical())
  if (run_mobster) {
    labelled = mutations_tool %>%
      distinct(mut_key, cluster_tool) %>%
      inner_join(mobster_labels, by="mut_key", relationship="many-to-many") %>%
      left_join(tail_samples, by="sample")
    cluster_vaf = labelled %>%
      group_by(cluster_tool, sample) %>%
      summarise(vaf=median(VAF, na.rm=TRUE), .groups="drop")
    shaped = cluster_vaf %>%
      inner_join(tail_modes, by="sample", relationship="many-to-many") %>%
      filter(vaf >= vaf_from, vaf < vaf_to) %>%
      distinct(cluster_tool, sample) %>%
      mutate(tail_mode=TRUE)
    valid = labelled %>%
      left_join(shaped, by=c("cluster_tool", "sample")) %>%
      mutate(valid=(tail_estimable %in% TRUE) & is.na(tail_mode))
    tail_table = clusters %>%
      select(cluster_tool) %>%
      left_join(valid %>%
                  group_by(cluster_tool) %>%
                  summarise(phi_tail={
                              ok = valid
                              if (!any(ok)) NA_real_ else {
                                never = tapply(label[ok] != "Tail", mut_key[ok], all)
                                mean(never)
                              }
                            },
                            tail_note=paste(unique(c(
                              tail_reason[!(tail_estimable %in% TRUE) & !is.na(tail_reason) & tail_reason != ""],
                              ifelse(any(!is.na(tail_mode)), "tail-shape mode at cluster VAF", NA)
                            )) %>% na.omit(), collapse="; "),
                            .groups="drop"),
                by="cluster_tool") %>%
      mutate(tail_note=ifelse(is.na(phi_tail) & (is.na(tail_note) | tail_note == ""),
                              "no MOBSTER mutations", tail_note))
    driver_tail = valid %>%
      filter(valid) %>%
      group_by(mut_key) %>%
      summarise(tail_driver=any(label == "Tail"), .groups="drop")
  }

  # phi_driver: base level by tail status, times the chance a driver count this
  # high is not expected from the sample's driver rate
  drivers = mutations_tool %>%
    filter(is_driver) %>%
    distinct(mut_key, cluster_tool, driver_label) %>%
    left_join(driver_tail, by="mut_key") %>%
    mutate(tail_driver=tail_driver %in% TRUE,
           gene=sub("_.*", "", driver_label))
  driver_rate = n_distinct(drivers[["mut_key"]]) / n_fit
  driver_table = clusters %>%
    select(cluster_tool, n, is_clonal) %>%
    left_join(drivers %>%
                group_by(cluster_tool) %>%
                summarise(k=n_distinct(mut_key),
                          all_tail=all(tail_driver),
                          genes=paste(sort(unique(gene)), collapse=","),
                          truncal=paste(sort(unique(gene[gene %in% truncal_genes])), collapse=","),
                          .groups="drop"),
              by="cluster_tool") %>%
    mutate(k=replace_na(k, 0L),
           b=case_when(k == 0 ~ 0, all_tail ~ 0.5, TRUE ~ 1),
           p_chance=ifelse(k > 0, pbinom(k - 1, n, driver_rate, lower.tail=FALSE), 1),
           phi_driver=if (flag("driver_chance")) b * (1 - p_chance) else b,
           truncal_flag=ifelse(!is_clonal & !is.na(truncal) & truncal != "",
                               paste0("truncal driver in subclonal cluster: ", truncal), ""))

  # phi_signature, SBS and ID
  sig_table = clusters %>% select(cluster_tool) %>% mutate(phi_sig=NA_real_, signature_null=NA_character_)
  sig_by_type = lapply(c("SBS", "ID"), function(sign_type) {
    sig = read_signature_inputs(tool, sign_type)
    if (is.null(sig)) return(tibble())

    table_signatures <<- table_signatures %>% bind_rows(
      sig[["activities"]] %>%
        pivot_longer(cols=-cluster_tool, names_to="Signature", values_to="Nmuts") %>%
        group_by(cluster_tool) %>%
        mutate(Ntot=sum(Nmuts), Exposure=Nmuts/Ntot, patient_id=patient_id) %>%
        ungroup() %>%
        left_join(clusters %>% select(cluster_tool, is_clonal), by="cluster_tool") %>%
        mutate(tool=tool, signature_type=sign_type)
    )
    if (length(anchor) == 0 || !anchor[1] %in% sig[["activities"]][["cluster_tool"]]) return(tibble())
    signature_scores(sig, anchor[1], intersect(clusters[["cluster_tool"]], sig[["activities"]][["cluster_tool"]])) %>%
      mutate(signature_type=sign_type)
  }) %>% bind_rows()

  if (nrow(sig_by_type) > 0) {
    sbs_ok = sig_by_type %>%
      filter(signature_type == "SBS", n_signature >= num("signature_min_mutations")) %>%
      pull(cluster_tool)
    sig_table = sig_by_type %>%
      filter(cluster_tool %in% sbs_ok, n_signature >= num("signature_min_mutations")) %>%
      group_by(cluster_tool) %>%
      summarise(phi_sig=max(phi_sig), signature_null=first(signature_null), .groups="drop") %>%
      right_join(clusters %>% select(cluster_tool), by="cluster_tool")
  }

  # Evidence Phi_C, reliability R_C and final score Phi*_C
  w = c(driver=num("w_driver"), signature=num("w_signature"), tail=w_tail)
  clusters %>%
    select(-q_size) %>%
    left_join(driver_table %>% select(cluster_tool, k, b, p_chance, phi_driver, genes, truncal_flag), by="cluster_tool") %>%
    left_join(tail_table, by="cluster_tool") %>%
    left_join(sig_table, by="cluster_tool") %>%
    left_join(reliability, by="cluster_tool") %>%
    mutate(phi_sig=ifelse(is_clonal, 1, phi_sig),
           score_driver=phi_driver,
           score_tail=ifelse(is.na(phi_tail), missing_value, phi_tail),
           score_sign=ifelse(is.na(phi_sig), missing_value, phi_sig),
           score_all=(w[["driver"]] * score_driver + w[["signature"]] * score_sign + w[["tail"]] * score_tail) /
             sum(w[w > 0]),
           R=ifelse(is_clonal, 1, R),
           limiting_factor=ifelse(is_clonal, "-", limiting_factor),
           score_final=R * score_all,
           tier=ifelse(is_clonal, "Tier 1", tier_of(score_final)),
           is_driver=k > 0,
           driver=as.numeric(k > 0),
           n_never_tail=phi_tail,
           weight_tail=w[["tail"]],
           clonal_anchor=ifelse(is_clonal, anchor_note, ""),
           tool=tool, patient_id=patient_id)
}) %>% bind_rows()

if (nrow(score_table) > 0) {
  score_table = score_table %>%
    rowwise() %>%
    mutate(flags=paste(c(truncal_flag, tail_note) %>% discard(~ is.na(.x) || .x == ""), collapse="; ")) %>%
    ungroup() %>%
    select(patient_id, cluster_tool, tool, is_clonal, is_clonal_tool, clonal_anchor, n, ccf,
           score_tail, score_sign, score_driver, score_all, R, limiting_factor, score_final, tier,
           best_sample, any_of(c("q_size", "q_support", "q_mode", "q_fit")),
           phi_tail, phi_sig, phi_driver, k, b, p_chance, genes, is_driver, driver,
           n_never_tail, weight_tail, signature_null, flags)
}

# Plots ######

get_signature_colors <- function(names) {
  n <- length(names)

  color_pool <- unique(c(
    brewer.pal(8, "Set1"),
    brewer.pal(12, "Set3")
  ))

  if (n > length(color_pool)) {
    color_pool <- colorRampPalette(color_pool)(n)
  }

  color_pool = setNames(color_pool[1:n], nm = names)
  return(color_pool)
}

get_cluster_colors <- function(cluster_names) {
  n <- length(cluster_names)

  color_pool <- unique(c(
    brewer.pal(8, "Dark2"),
    brewer.pal(12, "Paired")
  ))

  if (n > length(color_pool)) {
    color_pool <- colorRampPalette(color_pool)(n)
  }

  setNames(color_pool[1:n], nm = cluster_names)
}

cluster_names = unique(score_table[['cluster_tool']])
cluster_colors <- get_cluster_colors(cluster_names)

score_axes = c("score_driver"="Driver", "score_tail"="Tail", "score_sign"="Signature",
               "score_all"="Evidence", "R"="Reliability", "score_final"="Final")
if (w_tail == 0) score_axes = score_axes[names(score_axes) != "score_tail"]

pl_scores = ggplot() + theme_void()
pl_factors = ggplot() + theme_void()
if (nrow(score_table) > 0) {
  pl_scores = score_table %>%
    pivot_longer(cols=all_of(names(score_axes))) %>%
    mutate(name=factor(name, levels=names(score_axes))) %>%
    ggplot() +
    annotate("rect", xmin=-Inf, xmax=Inf, ymin=tiers[1], ymax=1, fill="palegreen4", alpha=0.2) +
    annotate("rect", xmin=-Inf, xmax=Inf, ymin=tiers[2], ymax=tiers[1], fill="goldenrod", alpha=0.2) +
    annotate("rect", xmin=-Inf, xmax=Inf, ymin=tiers[3], ymax=tiers[2], fill="salmon1", alpha=0.2) +
    annotate("rect", xmin=-Inf, xmax=Inf, ymin=0, ymax=tiers[3], fill="gainsboro", alpha=0.2) +
    geom_point(aes(x=name, y=value, col=cluster_tool, shape=is_driver), size=4) +
    geom_line(data=~filter(.x, !is.na(value)),
              aes(x=name, y=value, color=cluster_tool, group=cluster_tool), linewidth=.6) +
    geom_text(data=~filter(.x, is_clonal & name=="score_final"),
              aes(x=name, y=value, color=cluster_tool, group=cluster_tool, label="Clonal"), vjust=0.07, show.legend=F) +
    scale_color_manual("Cluster", values=cluster_colors) +
    scale_shape_manual('Contains Driver', values=c("FALSE"=4, "TRUE"=20)) +
    scale_x_discrete(labels=score_axes) +
    facet_grid(.~tool) +
    theme_bw() +
    theme(axis.title.x=element_blank())

  pl_factors = score_table %>%
    filter(!is_clonal) %>%
    pivot_longer(cols=any_of(c("q_size", "q_support", "q_mode", "q_fit"))) %>%
    mutate(name=factor(name, levels=c("q_size", "q_support", "q_mode", "q_fit"))) %>%
    ggplot(aes(x=name, y=cluster_tool, fill=value)) +
    geom_tile(color="white") +
    geom_text(aes(label=ifelse(is.na(value), "n/c", sprintf("%.2f", value))), size=3) +
    scale_fill_gradient("Factor", low="salmon1", high="palegreen4", limits=c(0, 1), na.value="gainsboro") +
    facet_grid(tool~., scales="free_y", space="free_y") +
    labs(y="Cluster", subtitle="Reliability factors in the best-resolved sample (R = minimum)") +
    theme_bw() +
    theme(axis.title.x=element_blank())
}

pl_signature = ggplot() + theme_void()
if (nrow(table_signatures) > 0) {
  names <- unique(table_signatures[["Signature"]])
  signature_colors <- get_signature_colors(names = names)

  pl_signature <- table_signatures %>%
    group_by(signature_type) %>%
    mutate(Nmut = sum(Nmuts)) %>%
    filter(Nmut > 50) %>%
    ungroup() %>%
    ggplot(aes(fill=Signature, y=Exposure, x=as.factor(cluster_tool))) +
    geom_bar(position="fill", stat="identity")+
    coord_flip()+
    scale_fill_manual(values=signature_colors)+
    theme_bw()+
    xlab("Cluster")+
    ylab("Exposures")+
    facet_grid(tool ~ signature_type, scales = 'free_y') +
    theme_bw() +
    theme(axis.title.x=element_blank())
}

if (length(unique(table_signatures[['tool']])) == 2){
  wd = 10
} else {
  wd = 5
}

if (length(unique(table_signatures[['signature_type']])) == 2){
  hg = 8
} else {
  hg = 4
}

pdf(paste0(opt[["prefix"]], "_scores_clusters.pdf"), width = wd, height = 4)
print(pl_scores)
print(pl_factors)
invisible(dev.off())
ggsave(pl_signature, filename=paste0(opt[["prefix"]], "_signature_clusters.pdf"), width = wd, height = hg, units = 'in')

saveRDS(object = list('score' = score_table,
                      'driver' = table_driver %>% distinct(),
                      'sample_quality' = sample_quality,
                      'cluster_sample' = cluster_sample_table,
                      'tail_samples' = tail_samples,
                      'tail_modes' = tail_modes,
                      'parameters' = opt),
        file = paste0(opt[["prefix"]], "_scores.rds"))
saveRDS(object = table_signatures, file = paste0(opt[["prefix"]], "_table_signature.rds"))
readr::write_tsv(score_table, paste0(opt[["prefix"]], "_scores.tsv"))
readr::write_tsv(sample_quality, paste0(opt[["prefix"]], "_sample_quality.tsv"))
