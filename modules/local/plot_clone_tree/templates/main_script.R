#!/usr/bin/env Rscript

library(tidyverse)
library(ggplot2)
library(ggraph)
library(igraph)
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

plot_tree <- function(ctree, signature, score, sub_tool, s_type){
  if (s_type == 'SBS'){
    s_colors = sbs_colors
  } else {
    s_colors = id_colors
  }

  CCF <- ctree[['CCF']] %>%
    left_join(score %>% filter(tool == sub_tool))

  ctree[['CCF']] <- CCF

  tb_tree = ctree[["tb_adj_mat"]] %>%
    left_join(score %>% filter(tool == sub_tool))

  cex = 1
  clones_orderings = igraph::topo_sort(igraph::graph_from_adjacency_matrix(DataFrameToMatrix(ctree[["transfer"]][["clones"]])),
                                       mode = 'out')\$name

  nDrivers = length(clones_orderings) - 1 # avoid GL
  layout_ctree <- create_layout(tb_tree, layout = "tree")

  pie_data <- layout_ctree %>%
    select(x, y, cluster,nMuts) %>%
    distinct() %>%
    left_join(signature %>% filter(tool == sub_tool, signature_type == s_type), by = "cluster") %>%
    mutate(Exposure = replace_na(Exposure, 0)) %>%
    group_by(cluster) %>%
    mutate(
      frac = Exposure / sum(Exposure, na.rm = TRUE),
      ymax = cumsum(frac),
      ymin = lag(ymax, default = 0)
    ) %>%
    ungroup()

  pie_data <- pie_data %>%
    mutate(
      r_node = 0.1 * cex + 0.25 * cex * (nMuts / max(nMuts, na.rm = TRUE))
    )


  ctree_plot <- ggraph(layout_ctree) +
  geom_edge_link(
    arrow = arrow(length = unit(2 * cex, 'mm')),
    end_cap = circle(5 * cex, 'mm'),
    start_cap = circle(5 * cex, 'mm')
  ) +
  ggforce::geom_arc_bar(
    data = pie_data,
    aes(
      x0 = x,
      y0 = y,
      r0 = 0,
      r = 0.25 * cex,
      start = 2 * pi * ymin,
      end   = 2 * pi * ymax,
      fill = Signature
    ),
    color = "white",
    linewidth = 0.2,
    inherit.aes = FALSE
  ) +
  ggrepel::geom_label_repel(
    aes(
      x = x,
      y = y,
      label = driver,
      colour = cluster
    ),
    na.rm = TRUE,
    nudge_x = .4,
    nudge_y = .4,
    size = 2.5 * cex,
    show.legend = F
  ) +
  geom_node_text(
    aes(
      label = cluster,
      colour = cluster
    ),
    vjust = 0,
    fontface = "bold",
    hjust = -3,
    show.legend = F
  ) +
  scale_color_manual(values = cluster_colors) +
  ggnewscale::new_scale_colour() +
  ggforce::geom_circle(
    data = distinct(layout_ctree, x, y, cluster, tier),
    aes(
      x0 = x,
      y0 = y,
      r = 0.25 * cex,
      color = tier
    ),
    inherit.aes = FALSE,
    linewidth = 1,
    alpha = 0.6
  ) +
  scale_color_manual(values = tier_colors) +
  coord_cartesian(clip = "off") +

    theme_void(base_size = 8 * cex) +
    theme(
      legend.position = "right"
    ) +
  scale_fill_manual(values = s_colors) +
  guides(
    fill = guide_legend("Signature"),
    colour = guide_legend("Tier")
  )
  return(ctree_plot)
}


# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
# Functions that are used to manipulate trees within REVOLVER
# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
# Three possible representaitons for a tree:
# - as edges strings "A~B", "B~C" etc.
# - as dataframe with from/ to clumns
# - as an adjacency matrix
# We have marshalling functions to switch among these representations
# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

# A~B to a DataFrame (from/to)
edgesToDataFrame = function(edges)
{
  # We convert the list of edges to bnlearn's from/to format
  dfedges = data.frame(from = NULL, to = NULL, stringsAsFactors = F)

  if(length(edges) == 0) return(dfedges)

  for(j in 1:length(edges))
  {
    aux = strsplit(edges[j], '~')[[1]]
    dfedges = rbind(dfedges, data.frame(from = aux[1], to = aux[2], stringsAsFactors = FALSE))
  }
  return(dfedges)
}

# A~B to a Adjacency Matrix
edgesToMatrix = function(edges)
{
  df = edgesToDataFrame(edges)
  vars = unique(unlist(df))
  matrix = matrix(0, nrow = length(vars), ncol = length(vars))
  colnames(matrix) = vars
  rownames(matrix) = vars

  if(nrow(df) == 0) return(matrix)

  for(j in 1:nrow(df))
  {
    matrix[df[j, 'from'], df[j, 'to']] = 1
  }

  return(matrix)
}

MatrixToDataFrame = function(matr)
{
  dfedges = data.frame(stringsAsFactors = F)
  for(i in 1:nrow(matr)) {
    for(j in 1:ncol(matr)){
      if(matr[i,j] == 1)
        dfedges = rbind(dfedges, data.frame(
          from = rownames(matr)[i],
          to = colnames(matr)[j],
          stringsAsFactors = FALSE))

    }
  }
  return(dfedges)
}

# Adjacency Matrix to A~B
MatrixToEdges = function(matr){
  return(DataFrameToEdges(MatrixToDataFrame(matr)))
}

# DataFrame to A~B
DataFrameToEdges = function(edges)
{
  edg = NULL
  for(j in 1:nrow(edges))
    edg = c(edg, paste(edges[j, 'from'], edges[j, 'to'], sep = '~'))
  return(edg)
}

# DataFrame to Adjacency Matrix
DataFrameToMatrix = function(edges)
{
  return(edgesToMatrix(DataFrameToEdges(edges)))
}

# Tibble graph to data frame
TidyGraphToDataFrame = function(x)
{
  M = TidyGraphToMatrix(x)

  MatrixToDataFrame(M)
}

# Tibble graph to adjacency matrix
TidyGraphToMatrix = function(x)
{
  # Nodes
  nodes = x %>%
    activate(nodes) %>%
    as_tibble()

  # Get nodes names and size
  labels = nodes %>% pull(!!colnames(nodes)[1])
  N_nodes = length(labels)

  # Edges
  edges = x %>%
    activate(edges) %>%
    as_tibble() %>%
    mutate(
      from = labels[from],
      to = labels[to]
    )

  M = matrix(0, nrow = N_nodes, ncol = N_nodes)
  colnames(M) = rownames(M) = labels

  if(nrow(edges) > 0)
  {
    for(i in 1:nrow(edges))  M[edge[["from"]][i], edges[["to"]][i]] = 1
  }

  M
}


tier_colors <- c(
  "Tier 1" = "palegreen4",
  "Tier 2" = "goldenrod",
  "Tier 3" = alpha("indianred",0.6),
  "Tier 4" = alpha("grey",0.6),
  "Missing"  = alpha("gainsboro",0.6)
)


score = readRDS("$rds_score") %>%
  dplyr::rename(cluster = cluster_tool,
                is_clonal_tool = is_clonal) %>%
  mutate(tier=case_when(
    is_clonal_tool                         ~ 'Tier 1',
    score_all >= .9                        ~ 'Tier 1',
    score_all >= .55 & score_all <.9       ~ 'Tier 2',
    score_all >= .2  & score_all <.55      ~ 'Tier 3',
    score_all <  .2                        ~ 'Tier 4'
  )) %>%
  select(cluster,tier,score_all, tool)

signature = readRDS("$rds_signature") %>%
  dplyr::rename(cluster = cluster_tool)

signature_type = unique(signature[["signature_type"]])
signature_name = unique(signature[["Signature"]])

cluster_names = unique(score[['cluster']])

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

# Usage
cluster_colors <- get_cluster_colors(cluster_names)
sbs_names <- signature_name[grep("^SBS", signature_name)]
id_names <-  signature_name[grep("^ID", signature_name)]
if (length(sbs_names) > 0){
  sbs_colors <- get_signature_colors(names = sbs_names)
}

if (length(id_names) > 0){
  id_colors <- get_signature_colors(names = id_names)
}


if ("$rds_ctree_pyclone" != ''){
  ctree_pyclone = readRDS("$rds_ctree_pyclone")[[1]]
  for (sign_type in signature_type){
    plt = plot_tree(ctree = ctree_pyclone, signature = signature, score = score, sub_tool = 'pyclonevi', s_type = sign_type)
    ggsave(plt, filename=paste0(opt[["prefix"]], "_tree_pyclonevi_",sign_type,".pdf"), width = 5, height = 6, units = 'in')

  }
}

if ("$rds_ctree_viber" != ''){
  ctree_viber = readRDS("$rds_ctree_viber")[[1]]
  for (sign_type in signature_type){
    plt = plot_tree(ctree = ctree_viber, signature = signature, score = score, sub_tool = 'viber', s_type = sign_type)
    ggsave(plt, filename=paste0(opt[["prefix"]], "_tree_viber_",sign_type,".pdf"), width = 5, height = 6, units = 'in')
  }
}
