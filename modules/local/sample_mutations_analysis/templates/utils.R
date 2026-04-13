# utils functions

compute_tmb = function(x, seq_length) {
  
  coding_muts = x %>% 
    separate(Consequence, into = "Consequence", sep = "&") %>% 
    filter(Consequence != 'synonymous_variant')
  
  coding_tmb = coding_muts %>% 
    group_by(sample) %>% 
    count() %>% 
    mutate(TMB = n/seq_length) 
  
  return(coding_tmb)
  
}

# utils --> colors for plots

mut_cols = c('SNV' = '#7FBC41', 'Indel' = '#DE77AE')
# 
# 
# color = grDevices::colors()[grep('gr(a|e)y', grDevices::colors(), invert = T)]
# 
# set.seed(1234)
# n=41
# pie(rep(1,n), col=sample(color, n))
# pie(rep(1,n), col=color[1:41])


# colors of mutations
consequences_colors = setNames(
  c('transcript_ablation', 
    'splice_acceptor_variant', 
    'splice_donor_variant', 
    'stop_gained',
    'frameshift_variant', 
    'stop_lost', 
    'start_lost',
    'transcript_amplification', 
    'feature_elongation', 
    'feature_truncation', 
    'inframe_insertion', 
    'inframe_deletion',
    'missense_variant', 
    'protein_altering_variant', 
    'splice_donor_5th_base_variant',
    'splice_region_variant', 
    'splice_donor_region_variant', 
    'splice_polypyrimidine_tract_variant', 
    'incomplete_terminal_codon_variant', 
    'start_retained_variant', 
    'stop_retained_variant', 
    'synonymous_variant', 
    'coding_sequence_variant', 
    'mature_miRNA_variant', 
    '5_prime_UTR_variant', 
    '3_prime_UTR_variant', 
    'non_coding_transcript_exon_variant', 
    'intron_variant', 
    'NMD_transcript_variant', 
    'non_coding_transcript_variant', 
    'coding_transcript_variant', 
    'upstream_gene_variant', 
    'downstream_gene_variant', 
    'TFBS_ablation', 
    'TFBS_amplification',
    'TF_binding_site_variant', 
    'regulatory_region_ablation', 
    'regulatory_region_amplification', 
    'regulatory_region_variant', 
    'intergenic_variant', 
    'sequence_variant'
    ), 
  object = c("wheat4",
             "darkseagreen2",
             "gold2",
             "mistyrose4",
             'coral',
             "bisque",
             "mediumpurple",
             "bisque3",
             "burlywood3",
             "blue3",
             "seashell4",
             "lightblue",
             "tan4",
             "orchid4",
             "darkorange",
             "indianred",
             "seashell2",
             "plum3",
             "thistle2",
             "skyblue4",
             "red1",
             "darkolivegreen3",
             "green",
             "tomato4",
             "turquoise4",
             "greenyellow",
             "cyan3",
             "slateblue3",
             "lightblue3",
             "tomato",
             "sandybrown",
             "blue",
             "violetred4",
             "yellowgreen",
             "lightskyblue1",
             "blue2",
             "salmon4",
             "darkseagreen1",
             "palegreen",
             "plum",
             "powderblue")
)
