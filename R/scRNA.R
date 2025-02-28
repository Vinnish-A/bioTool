
# scRNA -------------------------------------------------------------------


## Annotation ----

#' Map clusters to cell types
#'
#' Assigns cluster numbers to cell type labels. Unassigned clusters are grouped into the remaining group.
#'
#' @param ... Expressions that map cluster numbers to cell types.
#' @param nCluster_ Integer. The total number of clusters.
#' @param remain2_ Integer. The target group for unassigned clusters. Defaults to the first empty group if NULL.
#'
#' @return A named vector where keys are cluster numbers, and values are corresponding cell type labels.
#'
#' @examples
#' cluster2major(cluster_A = c(0, 1), cluster_B = c(2, 3), cluster_C = NULL, nCluster_ = 5)
#'
#' @export
cluster2major = function(..., nCluster_, remain2_ = NULL) {

  exprs_ = exprs(...)
  exprs_ = lapply(exprs_, eval)

  clusters_ = 1:nCluster_ - 1
  existing_ = unlist(exprs_)

  assertthat::assert_that(all(existing_ %in% clusters_))

  null_ = setdiff(clusters_, unlist(exprs_))
  if (!all(clusters_ %in% existing_)) {

    if (is.null(remain2_)) {
      exprs_[[which(sapply(exprs_, is.null))]] = null_
    } else {
      exprs_[[remain2_]] = c(exprs_[[remain2_]], null_)
    }

  }

  cluster_ = as.character(unlist(exprs_))
  cell_ = factor(rep(names(exprs_), sapply(exprs_, length)), levels = names(exprs_))

  table_ = setNames(cell_, cluster_)

  return(table_)

}

#' Update cell type assignments
#'
#' Updates the mapping of clusters to cell types based on the given focus cell type.
#' Unspecified clusters are added to the target group.
#'
#' @param ... Expressions that map cluster numbers to target groups.
#' @param cl2ma_ Named vector. Current mapping of clusters to cell types.
#' @param attention_ Character. The target cell type to focus on.
#' @param remain2_ Character. The target group for unassigned clusters. Defaults to be NULL.
#'
#' @return An updated named vector of cell type assignments.
#'
#' @importFrom assertthat assert_that
#'
#' @examples
#' cl2ma = cluster2major(cluster_A = c(0, 1), cluster_B = c(2, 3), cluster_C = NULL, nCluster_ = 5)
#' major2minor(cluster_A1 = 0, cluster_A2 = 1, cl2ma_ = cl2ma, attention_ = "cluster_A")
#'
#' @export
major2minor = function(..., cl2ma_, attention_, remain2_ = NULL) {

  exprs_ = exprs(...)
  exprs_ = lapply(exprs_, eval)

  remain_ = as.numeric(names(cl2ma_[cl2ma_ == attention_]))
  existing_ = unlist(exprs_)

  assertthat::assert_that(all(existing_ %in% remain_))

  null_ = setdiff(remain_, unlist(exprs_))

  if (!all(sort(existing_) == remain_)) {

    if (is.null(remain2_)) {
      exprs_[[which(sapply(exprs_, is.null))]] = null_
    } else {
      exprs_[[remain2_]] = c(exprs_[[remain2_]], null_)
    }

  }

  cluster_ = as.character(unlist(exprs_))
  cell_ = factor(rep(names(exprs_), sapply(exprs_, length)), levels = names(exprs_))

  table_min_ = setNames(cell_, cluster_)

  levels_ = as.list(levels(cl2ma_))
  levels_[[which(sapply(levels_, \(level__) level__ == attention_))]] = levels(cell_)
  levels_ = unlist(levels_)

  table_ = factor(c(cl2ma_[cl2ma_ != attention_], table_min_), levels = levels_)

  return(table_)

}

## modify ----

#' Merge Seurat objects
#'
#' Combines multiple Seurat objects into a single Seurat object.
#'
#' @param ... One or more Seurat objects to merge.
#' @param lst_ List. A list of Seurat objects to merge. Defaults to NULL, which uses the objects passed through `...`.
#'
#' @return A merged Seurat object.
#'
#' @examples
#' mergeSeurat(obj1, obj2, obj3)
#'
#' @export
mergeSeurat = function(..., lst_ = NULL) {

  if (is.null(lst_)) lst_ = list(...)

  x_ = lst_[[1]]
  y_ = lst_[-1]

  merge(x_, y_, add.cell.ids = T)

}

#' get_expr
#'
#' @noRd
#'
#' @export
get_expr = function(seu_, fea_, rename2expr_ = F) {

  if ('cell' %in% colnames(seu_@meta.data)) {

    res_ =  inner_join(
      as_tibble(FetchData(seu_, fea_), rownames = 'cell'),
      as_tibble(seu_@meta.data)
    )

  } else {

    res_ =  inner_join(
      as_tibble(FetchData(seu_, fea_), rownames = 'cell'),
      as_tibble(seu_@meta.data, rownames = 'cell')
    )

  }

  if (length(fea_) == 1 & rename2expr_) {
    colnames(res_)[[2]] = 'expression'
  }

  return(res_)

}

#' get_umap
#'
#' @noRd
#'
#' @export
get_umap = function(seu_) {

  if ('cell' %in% colnames(seu_@meta.data)) {

    seu_@meta.data$cell = NULL

  }

  if ('umap_1' %in% colnames(seu_@meta.data)) {

    seu_@meta.data$umap_1 = NULL

  }

  if ('umap_2' %in% colnames(seu_@meta.data)) {

    seu_@meta.data$umap_2 = NULL

  }

  inner_join(
    as_tibble(seu_@reductions$umap@cell.embeddings, rownames = 'cell'),
    as_tibble(seu_@meta.data, rownames = 'cell')
  )

}

#' get_umap8expr
#'
#' @noRd
#'
#' @export
get_umap8expr = function(seu_, fea_, rename2expr_ = F) {

  if ('cell' %in% colnames(seu_@meta.data)) {

    seu_@meta.data$cell = NULL

  }

  if ('umap_1' %in% colnames(seu_@meta.data)) {

    seu_@meta.data$umap_1 = NULL

  }

  if ('umap_2' %in% colnames(seu_@meta.data)) {

    seu_@meta.data$umap_2 = NULL

  }

  res_ = list(
    as_tibble(seu_@reductions$umap@cell.embeddings, rownames = 'cell'),
    as_tibble(seu_@meta.data, rownames = 'cell'),
    as_tibble(FetchData(seu_, fea_), rownames = 'cell')
  ) |> reduce(inner_join, by = 'cell')

  if (length(fea_) == 1 & rename2expr_) {
    colnames(res_)[which(colnames(res_) == fea_)] = 'expression'
  }

  return(res_)

}

# unstable
drop_outlier_each = function(dp_, x_ = 'umap_1', y_ = 'umap_2', distance_ = 0.5, nNeighbor_ = 5) {

  checkReliance('dbscan')

  cols_ = c(x_, y_)

  vec_group_ = dp_ |>
    select(all_of(cols_)) |>
    dbscan::dbscan(eps = distance_, minPts = nNeighbor_) |>
    _[['cluster']]

  dp_ |>
    filter(vec_group_ == 1)

}

#' drop_outlier
#'
#' @noRd
#'
#' @export
drop_outlier = function(dp_, group_, x_ = 'umap_1', y_ = 'umap_2', distance_ = 0.5, nNeighbor_ = 5) {

  lst_dp_ = split(dp_, dp_[[group_]])

  lst_res_ = map(lst_dp_, drop_outlier_each, x_ = x_, y_ = y_, distance_ = distance_, nNeighbor_ = nNeighbor_)

  bind_rows(lst_res_)

}

## plot ----

plot_sc = new.env()

#' Plot UMAP visualization
#'
#' Creates a UMAP plot with points colored by the specified grouping variable.
#'
#' @param df_ Data frame. Contains UMAP coordinates and grouping information.
#' @param x_ Character. The column name for the UMAP x-axis.
#' @param y_ Character. The column name for the UMAP y-axis.
#' @param colorBy_ Character. The column name for the grouping variable used for coloring points.
#' @param label_ Logic. Whether to add label.
#' @param border_ Logic. Whether to add border.
#' @param arrowSize_ Numeric. Relative size of arrow.
#'
#' @return A ggplot object displaying the UMAP visualization.
#'
#' @importFrom grid arrow
#'
#' @examples
#' plot_umap(df_, x_ = "UMAP_1", y_ = "UMAP_2", colorBy_ = "CellType")
#'
#' @export
plot_umap = function(df_, colorBy_ = 'seurat_clusters', x_ = 'umap_1', y_ = 'umap_2', label_ = T, border_ = T, arrowSize_ = 0.5) {

  params_chr_ = c(x_, y_, colorBy_)

  x_ = sym(x_)
  y_ = sym(y_)
  colorBy_ = sym(colorBy_)

  axis_ = guide_axis_truncated(
    trunc_lower = unit(0, "npc"),
    trunc_upper = unit(8 * arrowSize_, "cm")
  )

  nColor_ = length(unique(df_[[params_chr_[[3]]]]))
  color_ = get_color(nColor_)

  if (border_) {

    p_ = df_ |>
      ggplot() +
      ggunchull::stat_unchull(aes(!!x_, !!y_, color = !!colorBy_, fill = !!colorBy_), alpha = 0.2, size = 1, lty = 2, qval = 0.8, delta = 1) +
      geom_point(aes(!!x_, !!y_, color = !!colorBy_), size = 1) +
      guides(fill = "none", x = axis_, y = axis_) +
      guides(color = guide_legend(override.aes = list(size = 5))) +
      labs(color = "", title = paste("UMAP by", colorBy_)) +
      scale_x_continuous(breaks = NULL) +
      scale_y_continuous(breaks = NULL)+
      scale_fill_manual(values = color_) +
      scale_color_manual(values = color_) +
      theme(
        plot.title = element_text(hjust = 0.5),
        legend.key = element_blank(),
        aspect.ratio = 1,
        legend.position = "bottom",
        panel.background = element_blank(),
        panel.grid = element_blank(),
        axis.title = element_text(hjust = 0.05, face = "italic"),
        axis.line.x = element_line(arrow = arrow(type = "open", length = unit(arrowSize_, "cm"))),
        axis.line.y = element_line(arrow = arrow(type = "open", length = unit(arrowSize_, "cm")))
      )

  } else {

    p_ = df_ |>
      ggplot() +
      geom_point(aes(!!x_, !!y_, color = !!colorBy_), size = 1) +
      guides(fill = "none", x = axis_, y = axis_) +
      guides(color = guide_legend(override.aes = list(size = 5))) +
      labs(color = "", title = paste("UMAP by", colorBy_)) +
      scale_x_continuous(breaks = NULL) +
      scale_y_continuous(breaks = NULL) +
      scale_fill_manual(values = color_) +
      scale_color_manual(values = color_) +
      theme(
        plot.title = element_text(hjust = 0.5),
        legend.key = element_blank(),
        aspect.ratio = 1,
        legend.position = "bottom",
        panel.background = element_blank(),
        panel.grid = element_blank(),
        axis.title = element_text(hjust = 0.05, face = "italic"),
        axis.line.x = element_line(arrow = arrow(type = "open", length = unit(arrowSize_, "cm"))),
        axis.line.y = element_line(arrow = arrow(type = "open", length = unit(arrowSize_, "cm")))
      )

  }

  if (label_) {

    dp_label_ = df_ |>
      group_by(!!colorBy_) |>
      summarise(x = mean(!!x_), y = mean(!!y_))

    p_ = p_ +
      geom_shadowtext(aes(x, y, label = !!colorBy_), data = dp_label_, color = 'black', bg.color = 'white', bg.r = 0.2)

  }


  return(p_)

}

plot_umap_numeric = function(df_, colorBy_ = 'seurat_clusters', x_ = 'umap_1', y_ = 'umap_2', arrowSize_ = 0.5, colors = c('grey90', '#f38684')) {

  params_chr_ = c(x_, y_, colorBy_)

  x_ = sym(x_)
  y_ = sym(y_)
  colorBy_ = sym(colorBy_)

  axis_ = guide_axis_truncated(
    trunc_lower = unit(0, "npc"),
    trunc_upper = unit(8 * arrowSize_, "cm")
  )

  p_ = ggplot() +
    scattermore::geom_scattermore(aes(!!x_, !!y_), data = df_[df_[[params_chr_[[3]]]] == 0, ], color = colors[[1]], size = 1, pixels = c(512, 512)) +
    geom_point(aes(!!x_, !!y_, color = !!colorBy_), data = df_[df_[[params_chr_[[3]]]] != 0, ], size = 1) +
    guides(fill = "none", x = axis_, y = axis_) +
    scale_color_gradientn(colors = colors) +
    labs(color = "", title = paste("UMAP by", colorBy_)) +
    scale_x_continuous(breaks = NULL) +
    scale_y_continuous(breaks = NULL) +
    theme(
      plot.title = element_text(hjust = 0.5),
      legend.key = element_blank(),
      aspect.ratio = 1,
      legend.position = "bottom",
      panel.background = element_blank(),
      panel.grid = element_blank(),
      axis.title = element_text(hjust = 0.05, face = "italic"),
      axis.line.x = element_line(arrow = arrow(type = "open", length = unit(arrowSize_, "cm"))),
      axis.line.y = element_line(arrow = arrow(type = "open", length = unit(arrowSize_, "cm")))
    )

  return(p_)

}

plot_umap_raster = function(df_, colorBy_ = 'seurat_clusters', x_ = 'umap_1', y_ = 'umap_2', label_ = T, border_ = T, arrowSize_ = 0.5, rasterDPI_ = c(512, 512)) {

  params_chr_ = c(x_, y_, colorBy_)

  x_ = sym(x_)
  y_ = sym(y_)
  colorBy_ = sym(colorBy_)

  axis_ = guide_axis_truncated(
    trunc_lower = unit(0, "npc"),
    trunc_upper = unit(8 * arrowSize_, "cm")
  )

  nColor_ = length(unique(df_[[params_chr_[[3]]]]))
  color_ = get_color(nColor_)

  if (border_) {

    p_ = df_ |>
      ggplot() +
      ggunchull::stat_unchull(aes(!!x_, !!y_, color = !!colorBy_, fill = !!colorBy_), alpha = 0.2, size = 1, lty = 2, qval = 0.8, delta = 1) +
      scattermore::geom_scattermore(aes(!!x_, !!y_, color = !!colorBy_), size = 1, pixels = rasterDPI_) +
      guides(fill = "none", x = axis_, y = axis_) +
      guides(color = guide_legend(override.aes = list(size = 5))) +
      labs(color = "", title = paste("UMAP by", colorBy_)) +
      scale_x_continuous(breaks = NULL) +
      scale_y_continuous(breaks = NULL) +
      scale_fill_manual(values = color_) +
      scale_color_manual(values = color_) +
      theme(
        plot.title = element_text(hjust = 0.5),
        legend.key = element_blank(),
        aspect.ratio = 1,
        legend.position = "bottom",
        panel.background = element_blank(),
        panel.grid = element_blank(),
        axis.title = element_text(hjust = 0.05, face = "italic"),
        axis.line.x = element_line(arrow = arrow(type = "open", length = unit(arrowSize_, "cm"))),
        axis.line.y = element_line(arrow = arrow(type = "open", length = unit(arrowSize_, "cm")))
      )

  } else {

    p_ = df_ |>
      ggplot() +
      scattermore::geom_scattermore(aes(!!x_, !!y_, color = !!colorBy_), size = 1, pixels = rasterDPI_) +
      guides(fill = "none", x = axis_, y = axis_) +
      guides(color = guide_legend(override.aes = list(size = 5))) +
      labs(color = "", title = paste("UMAP by", colorBy_)) +
      scale_x_continuous(breaks = NULL) +
      scale_y_continuous(breaks = NULL) +
      scale_fill_manual(values = color_) +
      scale_color_manual(values = color_) +
      theme(
        plot.title = element_text(hjust = 0.5),
        legend.key = element_blank(),
        aspect.ratio = 1,
        legend.position = "bottom",
        panel.background = element_blank(),
        panel.grid = element_blank(),
        axis.title = element_text(hjust = 0.05, face = "italic"),
        axis.line.x = element_line(arrow = arrow(type = "open", length = unit(arrowSize_, "cm"))),
        axis.line.y = element_line(arrow = arrow(type = "open", length = unit(arrowSize_, "cm")))
      )

  }

  if (label_) {

    dp_label_ = df_ |>
      group_by(!!colorBy_) |>
      summarise(x = mean(!!x_), y = mean(!!y_))

    p_ = p_ +
      geom_shadowtext(aes(x, y, label = !!colorBy_), data = dp_label_, color = 'black', bg.color = 'white', bg.r = 0.2)

  }

  return(p_)

}

#' Plot Single-Cell Differential Expression (DEG) Results with GO Enrichment
#'
#' This function creates a combined visualization of single-cell differential expression (DEG) results.
#' It includes a heatmap of log2 fold changes, a dot plot of cluster annotations, and a bar plot for Gene Ontology (GO) enrichment results.
#'
#' @param seu_ A Seurat object containing the single-cell RNA-seq data.
#' @param marker_all_ A data frame with all marker gene data containing columns `gene`, `cluster`, and `avg_log2FC`.
#' @param num_each_ An integer specifying the number of genes to display per cluster (default is 15).
#' @param ident_ A vector of cluster identities to use in the analysis (default is the cluster identities from the Seurat object).
#' @param breaks_ A numeric vector specifying the breaks for the color scale in the heatmap (default is c(-5, 0, 5)).
#' @param organism_ A character string specifying the organism for GO enrichment analysis. Options are `'hsa'` for human (default) or `'mmu'` for mouse.
#' @return A combined plot consisting of:
#'   - A heatmap showing the average log2 fold change for genes across clusters.
#'   - A dot plot showing the cluster annotations.
#'   - A bar plot showing the top 5 GO terms enriched in each cluster.
#'
#' @examples
#' plot_sc_deg(seurat_object, marker_data, num_each_ = 20, ident_ = c('Cluster1', 'Cluster2'), organism_ = 'hsa')
#'
#' @export
plot_sc_deg = function(seu_, marker_all_, num_each_ = 15, ident_ = Idents(seu_), breaks_ = c(-5, 0, 5), organism_ = 'hsa') {

  ident_ = levels(Idents(seu_))

  deg_sliced_ = marker_all_ |>
    distinct(gene, .keep_all = T) |>
    filter(avg_log2FC > 0) |>
    group_by(cluster) |>
    slice_max(avg_log2FC, n = 6*num_each_) |>
    slice_sample(n = num_each_) |>
    ungroup() |>
    arrange(cluster) |>
    pull(gene, cluster)

  dp_logfc_ = map(
    ident_, \(cluster_) {
      marker_each_ = FindMarkers(seu_, ident.1 = cluster_, only.pos = F, features = deg_sliced_)
      marker_each_$cluster = cluster_
      return(as_tibble(marker_each_, rownames = "gene"))
    }
  ) |> bind_rows() |>
    arrange(cluster) |>
    mutate(cluster = factor(cluster),
           gene = factor(gene, levels = deg_sliced_)) |>
    select(cluster, avg_log2FC, gene)

  patch_ = expand.grid(cluster = ident_, gene = deg_sliced_) |>
    as_tibble() |>
    mutate(V = paste(cluster, gene, sep = '_')) |>
    filter(!V %in% with(dp_logfc_, paste(cluster, gene, sep = '_'))) |>
    select(-V) |>
    mutate(avg_log2FC = runif(length(row_number()), -0.75, 0.25))

  dp_logfc_ = bind_rows(dp_logfc_, patch_)

  if (organism_ == 'hsa') {
    enrich_ = compareCluster(gene ~ cluster, data = marker_all_, fun = 'enrichGO', OrgDb = 'org.Hs.eg.db', keyType = 'SYMBOL')
  } else if (organism_ == 'mmu') {
    enrich_ = compareCluster(gene ~ cluster, data = marker_all_, fun = 'enrichGO', OrgDb = 'org.Mm.eg.db', keyType = 'SYMBOL')
  }

  dp_enrich_ = enrich_@compareClusterResult |>
    filter(cluster != "NA") |>
    group_by(cluster) |>
    slice_max(-p.adjust, n = 5, with_ties = F) |>
    ungroup() |>
    mutate(cluster = factor(cluster, levels = rev(ident_)))

  limits_ = min(abs(range(dp_logfc_$avg_log2FC)))
  limits_ = c(-limits_, limits_)

  p1_ = dp_logfc_ |>
    ggplot(aes(x = cluster,
               y = reorder(gene, -as.numeric(cluster)),
               fill = avg_log2FC,
               color = avg_log2FC)) +
    geom_tile() +
    xlab("") +
    ylab("") +
    labs(fill = "", color = "") +
    scale_color_gradient2(low = "#01665e", mid = "white", high = "#8c510a", breaks = breaks_, limits = limits_, oob = scales::squish) +
    scale_fill_gradient2(low = "#01665e", mid = "white", high = "#8c510a", breaks = breaks_, limits = limits_, oob = scales::squish) +
    scale_x_discrete(breaks = NULL) +
    theme_bw() +
    theme(legend.position = "top",
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_rect(color = "black", linewidth = 0.2, fill = NA),
          plot.margin = margin(0.2, 0.2, 0.2, 0.2, "cm"),
          legend.key.size = unit(0.5, "cm"),
          legend.key.width = unit(0.5, "cm"))

  ncolor_ = length(ident_)
  colors_ = get_color(ncolor_)

  p2_ = tibble(x = 0, y = ident_) |>
    ggplot(aes(x = y, y = 1, color = factor(y))) +
    geom_point(size = 6, show.legend = F) +
    scale_color_manual(values = colors_) +
    scale_y_continuous(expand = c(0,0)) +
    theme(legend.position = "none",
          panel.spacing = unit(0, "lines"),
          panel.background = element_blank(),
          panel.border = element_blank(),
          plot.background = element_blank(),
          plot.margin = margin(0, 0, 0, 0, "pt"),
          axis.text.x = element_text(angle = 30,
                                     size = 12,
                                     hjust = 1,
                                     vjust = 1,
                                     color = "black"),
          axis.title  = element_blank(),
          axis.ticks = element_blank(),
          panel.grid = element_blank(),
          axis.text.y = element_blank())

  limits_ = c(0, mean(sort(-log10(dp_enrich_$p.adjust), T)[1:6]))

  p3_ = dp_enrich_ |>
    ggplot(aes(x = reorder(Description, -log10(p.adjust)),
               y = -log10(p.adjust),
               fill = cluster)) +
    geom_bar(position = position_dodge(),
             stat = "identity",
             show.legend = F) +
    scale_x_discrete(position = "top") +
    scale_y_continuous(expand = c(0,0), limits = limits_, oob = scales::squish, labels = scales::number_format(accuracy = 1)) +
    scale_fill_manual(values = rev(colors_))+
    facet_wrap(~ cluster, ncol = 1, scales = "free_y") +
    theme_bw() +
    theme(strip.background = element_blank(),
          strip.text = element_blank(),
          axis.title.y = element_blank(),
          panel.spacing = unit(0, "lines"),
          axis.text = element_text(size = 12),
          panel.border = element_rect(colour = "grey50")) +
    ylab("-Log10(Padj)") +
    coord_flip()

  p1_ + p3_ + p2_ + plot_spacer() + plot_layout(widths = c(2, 1), height = c(12, 1))

}

### register ----

plot_sc = new.env()

register(plot_umap, plot_sc)
register(plot_sc_deg, plot_sc)
