#!/usr/bin/env Rscript

parse_args = function(x) {
    x = gsub("\\\\[","",x)
    x = gsub("\\\\]","",x)
    # giving errors when we have lists like c(xxx, xxx) since it will separate it
    # args_list = unlist(strsplit(x, ', ')[[1]])
    args_list = unlist(strsplit(x, ", (?=[^)]*(?:\\\\(|\$))", perl=TRUE))
    # args_vals = lapply(args_list, function(x) strsplit(x, split=":")[[1]])
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


# Load packages
library(bascule)
library(dplyr)
library(ggplot2)

# Get input data
# tsv_join contains a set of dataframes with event counts
# and a column called "type" with "SBS"/"DBS"/etc
# rownames: samples; colnames: contexts + type
counts_tsv = strsplit("$tsv_join", " ")[[1]]
counts_tmp = lapply(counts_tsv, function(p_table) {
    read.delim(p_table, sep = "\\t", header=T) %>%
        mutate(across(everything(), as.character))
    })

counts = lapply(counts_tmp, function(x) x %>% dplyr::select(-type)) %>%
    setNames(counts_tmp %>% bind_rows() %>% pull(type) %>% unique())

x = fit(counts=counts, k_list=as.integer(opt[["k_list"]])
        # reference_cat=reference_cat
        )

x_refined = refine_denovo_signatures(x)
x_refined_cluster = fit_clustering(x_refined,
                                   cluster=as.integer(opt[["cluster"]]))
x_refined_cluster = merge_clusters(x_refined_cluster)

pl_signatures = plot_signatures(x_refined_cluster)
pl_exposures = plot_exposures(x_refined_cluster)

pl_all = patchwork::wrap_plots(pl_exposures, pl_signatures, ncol=2) +
    patchwork::plot_annotation(title = "$meta.id")
ggplot2::ggsave(plot=plt_all, filename=paste0(opt[["prefix"]], "_plot_all.pdf"),
                width=210, height=297, units="mm")
saveRDS(object=plt_all, file=paste0(opt[["prefix"]], "_plot_all.rds"))


# version export
f = file("versions.yml","w")
bascule_version = sessionInfo()\$otherPkgs\$bascule\$Version
dplyr_version = sessionInfo()\$otherPkgs\$dplyr\$Version
ggplot2_version = sessionInfo()\$otherPkgs\$ggplot2\$Version
patchwork_version = sessionInfo()\$otherPkgs\$patchwork\$Version
writeLines(paste0('"', "$task.process", '"', ":"), f)
writeLines(paste("    bascule:", bascule_version), f)
writeLines(paste("    dplyr:", dplyr_version), f)
writeLines(paste("    ggplot2:", ggplot2_version), f)
writeLines(paste("    patchwork:", patchwork_version), f)
close(f)
