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


mutation_tables = strsplit("$mutation_tables", " ")[[1]]
results_sigprofiler = strsplit("$results_sigprofiler", " ")[[1]]

cosine_similarity = function(vec1, vec2) {
  sum(vec1 * vec2) / (sqrt(sum(vec1^2)) * sqrt(sum(vec2^2)))
}

compare_signatures = function(df1, df2) {
  df1_f = df1 %>% filter(Exposure>.05)
  df2_f = df2 %>% filter(Exposure>.05)

  # Check if signatures are identical
  sigs_match = setequal(df1_f[["Signature"]], df2_f[["Signature"]])
  ndiff_1 = length(setdiff(df2_f[["Signature"]], df1_f[["Signature"]]))
  ndiff_2 = length(setdiff(df1_f[["Signature"]], df2_f[["Signature"]]))
  ndiff = sum(ndiff_1 + ndiff_2)
  ntot = length(c(df1_f[["Signature"]], df2_f[["Signature"]]) %>% unique())

  comparison = full_join(df1 %>% select(Signature, Exposure, n_sig),
                         df2 %>% select(Signature, Exposure, n_sig),
                         by="Signature", suffix=c("_tmp", "_clonal")) %>%
    mutate(diff_exposure=Exposure_tmp-Exposure_clonal,
           diff_n_sig=n_sig_tmp-n_sig_clonal) %>%
    mutate(across(.cols=where(is.numeric),
                  .fns=function(i) replace_na(i, 0.)))

  comparison_data = full_join(df1, df2, by="Signature",
                              suffix=c("_tmp", "_clonal")) %>%
    select(Signature, Exposure_tmp, Exposure_clonal, n_sig_tmp, n_sig_clonal) %>%
    mutate(across(where(is.numeric), ~replace_na(., 0)))

  # Calculate similarities
  cos_sim_exposure = cosine_similarity(comparison_data[["Exposure_tmp"]],
                                       comparison_data[["Exposure_clonal"]])

  return(list("df"=comparison,
              "match"=sigs_match,
              "cs_exp"=cos_sim_exposure,
              "n_diff"=ndiff,
              "n_tot"=ntot))
}

# Signature interpretation ######
patient_id = opt[["prefix"]]

mobster_files = grep(mutation_tables, pattern="mobster", value=T)
if (length(mobster_files)>0){
  mutations_mobster = lapply(mobster_files, FUN = function(f){
    readr::read_tsv(f) %>% mutate(chrom = as.character(chrom))
    }) %>%
    bind_rows() %>%
    mutate(patient_id=patient_id, tool="mobster",
           mutation_id=paste(chrom, pos_start, pos_end, ref, alt, sep=":")) %>%
    select(Project, patient_id, mutation_id, everything(), -chrom, -pos_start,
           -pos_end, -ref, -alt, -Type, -ID, -Genome, -mut_type) %>%
    rename(cluster_mobster=Sample)

  if ( 'Tail' %in% unique((mutations_mobster[["cluster_mobster"]]))){
    w_tail = 1
  } else {
    w_tail = 0
  }
}



table_signatures <- tibble()
table_driver <- tibble()
tool_list <- strsplit(opt[['tools']], ",")[[1]]
tool_list <- tool_list[tool_list %in% c("viber", "pyclone-vi")]
tool_list <- gsub("pyclone-vi", "pyclonevi", tool_list)

score_table = lapply(tool_list, function(tool) {
  lapply(c("SBS", "ID"), function(sign_type) {
    mutations_file = grep(mutation_tables, pattern=tool, value=T)

    if (length(mutations_file) == 0) return(tibble())

    mutations_tool = readr::read_tsv(mutations_file) %>%
      mutate(patient_id=patient_id, tool=tool,
              mutation_id=paste(chrom, pos_start, pos_end, ref, alt, sep=":")) %>%
      select(Project, patient_id, mutation_id, everything(), -chrom, -pos_start,
              -pos_end, -ref, -alt, -Type, -ID, -Genome, -mut_type) %>%
      rename(cluster_tool=Sample)

    table_driver <<- bind_rows(table_driver, mutations_tool %>% filter(is_driver))

    if (length(mobster_files>0)){
      never_tail_muts = mutations_mobster %>%
        group_by(mutation_id) %>%
        summarise(never_tail=all(cluster_mobster != "Tail"))

      final_table_subclonal = mutations_tool %>%
        left_join(never_tail_muts) %>%
        group_by(cluster_tool) %>%
        reframe(n_tail = sum(!is.na(never_tail)),
                n=n(),
                n_never_tail=sum(never_tail, na.rm=T) / n_tail,
                is_driver=any(is_driver), is_clonal=all(is_clonal))
    } else {
      final_table_subclonal = mutations_tool %>%
        group_by(cluster_tool) %>%
        reframe(n_never_tail=NA,
                is_driver=any(is_driver),
                is_clonal=all(is_clonal))
    }

    # one dir per sample
    dir_path = list.dirs(grep(results_sigprofiler, pattern=tool, value=T), recursive=T) %>%
      keep(function(x) grepl(pattern=sign_type, x=x)) %>%
      keep(function(x) grepl(pattern="/Activities", x=x))

    if (length(dir_path) == 0) return(tibble())

    signatures = read.table(file=file.path(dir_path, "Assignment_Solution_Activities.txt"), header=T) %>%
      rename(cluster_tool=Samples) %>%
      pivot_longer(cols=starts_with(sign_type), names_to="Signature", values_to="Nmuts") %>%
      group_by(cluster_tool) %>%
      mutate(Ntot=sum(Nmuts), Exposure=Nmuts/Ntot, patient_id=patient_id) %>%
      ungroup() %>%
      left_join(mutations_tool %>% select(cluster_tool, is_clonal) %>% unique())

    table_signatures <<- table_signatures %>% bind_rows(
      signatures %>% mutate(tool=tool, signature_type=sign_type)
    )

    clonal_sigs = signatures %>% filter(is_clonal, Nmuts>0) %>% rename(n_sig=Nmuts)
    other_sigs = signatures %>% filter(!is_clonal, Nmuts>0) %>% rename(n_sig=Nmuts)

    background_sigs = other_sigs %>%
      group_by(Signature) %>%
      # n total of mutations per cluster
      summarise(n_sig=sum(n_sig)) %>%
      mutate(n=sum(n_sig)) %>%
      mutate(Exposure=n_sig/n)

    tmp_final_table = lapply(unique(other_sigs[["cluster_tool"]]), function(c) {
      tmp = other_sigs %>% filter(cluster_tool == c)
      result_table = compare_signatures(df1=tmp, df2=clonal_sigs)
      result = result_table[["df"]] %>% mutate(match=result_table[["match"]],
                                          cs_exp=result_table[["cs_exp"]],
                                          cs_nsig=result_table[["cs_n"]],
                                          n_diff=result_table[["n_diff"]],
                                          n_tot=result_table[["n_tot"]]) %>%
        mutate(cluster_tool=as.character(c))
      return(result)
    }) %>% bind_rows()

    tmp_final_table_background = lapply(unique(other_sigs[["cluster_tool"]]), function(c){
      tmp = other_sigs %>% filter(cluster_tool == c)
      result_table = compare_signatures(df1=tmp, df2=background_sigs)
      result = result_table[["df"]] %>% mutate(match=result_table[["match"]],
                                          cs_exp=result_table[["cs_exp"]],
                                          cs_nsig=result_table[["cs_n"]],
                                          n_diff=result_table[["n_diff"]],
                                          n_tot=result_table[["n_tot"]]) %>%
        mutate(cluster_tool=as.character(c))
      return(result)
    }) %>% bind_rows()

    if (nrow(tmp_final_table) > 0) {
      f_background = tmp_final_table_background %>%
        mutate(n_rel=n_diff/n_tot) %>%
        select(cluster_tool, match, cs_exp, n_rel) %>%
        distinct() %>%
        dplyr::rename(bg_match=match, bg_cs_exp=cs_exp, bg_n_rel=n_rel) %>%
        filter(!is.na(cluster_tool))

      final_table_signature = tmp_final_table %>%
        mutate(n_rel=n_diff/n_tot) %>%
        select(cluster_tool, match, cs_exp, n_rel) %>%
        distinct() %>%
        left_join(f_background)

      final_table_subclonal %>%
        left_join(final_table_signature) %>%
        mutate(driver=ifelse(is_driver == F, 0, 1),
               bg_cs= 1 - bg_cs_exp,
               cs_sign=1 - cs_exp,
               cs_sign=ifelse(is.na(cs_exp) & is_clonal == T, 1, cs_sign),
               bg_sign=ifelse(is.na(bg_cs) & is_clonal == T, 1, bg_cs),
               n_rel=ifelse(is.na(cs_exp) & is_clonal == T, 1, n_rel)) %>%
        mutate(tool=tool, signature_type=sign_type, patient_id=patient_id)
    }
  }) %>% bind_rows()
}) %>% bind_rows() %>% select(patient_id, cluster_tool, everything())


w_sig = 1
w_driver = 1
score_table = score_table %>%
  group_by(patient_id, cluster_tool, tool, is_clonal, is_driver, n_never_tail, driver) %>%
  reframe(cs_sign = ifelse(is_clonal == T, 1, min(cs_sign, na.rm = T)),
            n_rel = ifelse(is_clonal == T, 1, max(n_rel, na.rm = T)),
            bg_sign =  ifelse(is_clonal == T, 1, min(bg_sign, na.rm = T))) %>%
  rowwise() %>%
  mutate(
    score_driver = driver,
    score_sign = (cs_sign + n_rel + bg_sign) / 3,
    score_tail = n_never_tail,
    score_all = ifelse(is.na(n_never_tail),
                       (w_driver * score_driver + w_sig * score_sign) / (w_driver+w_sig),
                       (w_driver * score_driver + w_sig * score_sign + w_tail * score_tail) / (w_driver+w_sig+w_tail))) %>%
  mutate(weight_tail = w_tail)

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

if (w_tail == 0){
  pl_scores = score_table %>%
    pivot_longer(cols=c(score_driver, score_all, score_sign)) %>%
    mutate(name=factor(name, levels=c("score_driver", "score_sign", "score_all"))) %>%
    ggplot() +
    annotate("rect", xmin=-Inf, xmax=Inf, ymin=0.9, ymax=1, fill="palegreen4", alpha=0.2) +
    annotate("rect", xmin=-Inf, xmax=Inf, ymin=0.55, ymax=.9, fill="goldenrod", alpha=0.2) +
    annotate("rect", xmin=-Inf, xmax=Inf, ymin=0.2, ymax=0.55, fill="salmon1", alpha=0.2) +
    annotate("rect", xmin=-Inf, xmax=Inf, ymin=0.2, ymax=0, fill="gainsboro", alpha=0.2) +
    geom_point(aes(x=name, y=value, col=cluster_tool, shape = is_driver), size=4) +
    geom_line(data=~filter(.x, !is.na(value)),
              aes(x=name, y=value, color=cluster_tool, group=cluster_tool), linewidth=.6)  +
    geom_text(data=~filter(.x, is_clonal & name=="score_all"),
              aes(x=name, y=value, color=cluster_tool, group=cluster_tool, label="Clonal"), vjust=0.07, show.legend = F) +
    scale_color_manual("Cluster", values=cluster_colors) +
    scale_shape_manual('Contains Driver', values = c(4, 20)) +
    scale_x_discrete(labels=c("score_driver"="Driver",
                              #"score_tail"="Tail",
                              "score_sign"="Signature",
                              "score_all"="All")) +
    facet_grid(.~tool) +
    theme_bw() +
    theme(axis.title.x=element_blank())

} else {
  pl_scores = score_table %>%
    pivot_longer(cols=c(score_driver, score_all, score_tail, score_sign)) %>%
    mutate(name=factor(name, levels=c("score_driver", "score_tail", "score_sign", "score_all"))) %>%
    ggplot() +
    annotate("rect", xmin=-Inf, xmax=Inf, ymin=0.9, ymax=1, fill="palegreen4", alpha=0.2) +
    annotate("rect", xmin=-Inf, xmax=Inf, ymin=0.55, ymax=.9, fill="goldenrod", alpha=0.2) +
    annotate("rect", xmin=-Inf, xmax=Inf, ymin=0.2, ymax=0.55, fill="salmon1", alpha=0.2) +
    annotate("rect", xmin=-Inf, xmax=Inf, ymin=0.2, ymax=0, fill="gainsboro", alpha=0.2) +
    geom_point(aes(x=name, y=value, col=cluster_tool, shape = is_driver), size=4) +
    geom_line(data=~filter(.x, !is.na(value)),
              aes(x=name, y=value, color=cluster_tool, group=cluster_tool), linewidth=.6)  +
    geom_text(data=~filter(.x, is_clonal & name=="score_all"),
              aes(x=name, y=value, color=cluster_tool, group=cluster_tool, label="Clonal"), vjust=0.07, show.legend = F) +
    scale_color_manual("Cluster", values=cluster_colors) +
    scale_shape_manual('Contains Driver', values = c(4, 20)) +
    scale_x_discrete(labels=c("score_driver"="Driver",
                              "score_tail"="Tail",
                              "score_sign"="Signature",
                              "score_all"="All")) +
    facet_grid(.~tool) +
    theme_bw() +
    theme(axis.title.x=element_blank())

}

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

ggsave(pl_scores, filename=paste0(opt[["prefix"]], "_scores_clusters.pdf"), width = wd, height = 4, units = 'in')
ggsave(pl_signature, filename=paste0(opt[["prefix"]], "_signature_clusters.pdf"), width = wd, height = hg, units = 'in')
saveRDS(object = list('score' = score_table, 'driver' = table_driver %>% distinct()), file = paste0(opt[["prefix"]], "_scores.rds"))
saveRDS(object = table_signatures, file = paste0(opt[["prefix"]], "_table_signature.rds"))
