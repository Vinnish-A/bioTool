
# plot --------------------------------------------------------------------

#' Register a Function Across Multiple Environments
#'
#' @export
register = function(fun_, ...) {

  name_full_ = as.character(substitute(fun_))
  name_last_ = strsplit(name_full_, '_')[[1]]
  name_last_ = name_last_[[length(name_last_)]]

  fun_ = get(name_full_, envir = parent.frame())
  envs_ = list(...)

  for (i_ in seq_along(envs_)) assign(name_last_, fun_, envir = envs_[[i_]])

}

#' theme_dropx
#'
#' @export
theme_dropx = function() {

  theme(
    axis.line.x = element_blank(),
    axis.text.x = element_blank(),
    axis.title.x = element_blank(),
    axis.ticks.x = element_blank()
  )

}

#' theme_dropy
#'
#' @export
theme_dropy = function() {

  theme(
    axis.line.y = element_blank(),
    axis.text.y = element_blank(),
    axis.title.y = element_blank(),
    axis.ticks.y = element_blank()
  )

}

## universal ----

#' plot_stacked
#'
#' @noRd
#'
#' @export
plot_stacked = function(df_, group_, subgroup_, n_ = 'n', thres_ = 0.1, attention_ = NULL, summarise_ = F) {

  if (summarise_) {

    dp_bar_ = tibble(group = df_[[group_]],
                     subgroup = df_[[subgroup_]]) |>
      group_by(group, subgroup) |>
      summarise(n = n())

  } else {

    dp_bar_ = tibble(group = df_[[group_]],
                     subgroup = df_[[subgroup_]],
                     n = df_[[n_]])

  }

  dp_bar_ = dp_bar_ |>
    arrange(group, subgroup) |>
    group_by(group) |>
    mutate(ratio = signif(n/sum(n), 3),
           label = ifelse(ratio > thres_, paste0(ratio * 100, "%"), NA),
           postion_top = rev(cumsum(rev(ratio))),
           postion_bottom = postion_top - ratio)

  if (!is.null(attention_)) {
    dp_bar_ = dp_bar_ |>
      mutate(label = ifelse(subgroup %in% attention_, label, NA))
  }

  nGroup_ = length(unique(dp_bar_[['subgroup']]))
  colors_ = get_color(nGroup_)

  dp_bar_ |>
    ggplot() +
    geom_bar(aes(group, n, fill = subgroup), color = "#f3f4f4", position = "fill", stat = "identity", size = 1) +
    geom_text(aes(group, (postion_top + postion_bottom)/2, label = label), size = 4, color = "white", fontface = 'bold') +
    scale_y_continuous(labels = paste0(100 * seq(0, 1, 0.25), "%")) +
    scale_fill_manual(values = colors_) +
    xlab(NULL) +
    ylab(NULL) +
    labs(fill = NULL) +
    theme_classic() +
    theme(
      legend.position = "top",
      axis.text.x = element_text(angle = 30, vjust = 0.5, hjust = 0.75)
    )

}

#' plot_rain
#'
#' @noRd
#'
#' @export
plot_rain = function(dp_, x_, group_) {

  nColor_ = length(unique(dp_[[group_]]))
  color_ = get_color(nColor_)

  test_ = kruskal.test(dp_[[x_]], dp_[[group_]])

  x_ = sym(x_)
  group_ = sym(group_)

  p_top_ = dp_ |>
    ggplot(aes(x = !!x_, color = !!group_, fill = !!group_)) +
    geom_density(alpha = 0.7) +
    scale_fill_manual(values = color_) +
    scale_color_manual(values = color_) +
    theme_classic() +
    labs(title = 'SLE Disease Activity Index', subtitle = paste0('Mann-Whitney p = ', signif(test_$p.value, 3))) +
    xlab(NULL) +
    ylab(NULL) +
    theme(legend.position = "none",
          legend.title = element_blank(),
          axis.text.x = element_text(size = 12,color = "black"),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          axis.line.y = element_blank(),
          panel.background = element_blank(),
          panel.grid.major = element_blank(),
          plot.title = element_text(hjust = 1),
          plot.subtitle = element_text(hjust = 1 ),
          panel.grid.minor = element_blank()) +
    geom_rug()

  p_bot_ = dp_ |>
    ggplot(aes(!!group_, !!x_, fill = !!group_)) +
    geom_boxplot(aes(col = !!group_)) +
    scale_fill_manual(values = color_) +
    scale_color_manual(values = color_) +
    xlab(NULL) +
    ylab(NULL) +
    theme_void() +
    theme(legend.position = "none",
          legend.title = element_blank(),
          axis.text.x = element_blank(), #
          axis.text.y = element_text(size = 11,color = "black"),
          panel.background = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank()) +
    coord_flip()

  dp2_ = ggplot_build(p_bot_)$data[[1]]
  p_bot_ = p_bot_ + geom_segment(data=dp2_, aes(x=xmin, xend=xmax, y=middle, yend=middle), color="white", inherit.aes = F)

  p_ = p_top_ |> insert_bottom(p_bot_, height = 0.4)
  p_

}

#' plot_violin
#'
#' @noRd
#'
#' @export
plot_violin = function(dp_, celltype_, expr_ = 'expression', dot_ = F) {

  nColor_ = length(unique(dp_[[celltype_]]))
  color_ = get_color(nColor_)

  params_chr_ = c(celltype_, expr_)

  celltype_ = sym(celltype_)
  expr_ = sym(expr_)

  stat = kruskal.test(dp_[[expr_]], dp_[[celltype_]])

  p_vln_ = dp_ |>
    ggplot(aes(!!celltype_, !!expr_, fill = !!celltype_)) +
    geom_jitter(alpha = 0.2, shape = 21, color = 'white', width = 0.2) +
    geom_violin(alpha = 0.4, scale = 'width') +
    scale_y_continuous(expand = c(0, 0)) +
    scale_fill_manual(values = color_) +
    labs(y = 'Norm Expression', subtitle = paste0('Kruskal-Wallis Pvalue: ', signif(stat$p.value, 3))) +
    theme_classic() +
    theme_dropx() +
    theme(legend.title = element_blank(),
          plot.subtitle = element_text(hjust = 1))

  if (dot_) {

    p_dot_ = tibble(x = sort(unique(dp_[[params_chr_[[1]]]]))) |>
      ggplot() +
      geom_point(aes(x, 1, color = x), size = 4) +
      scale_color_manual(values = color_) +
      theme_void() +
      theme(legend.position = 'none')

    p_vln_ = p_vln_ / p_dot_

  }

  return(p_vln_)

}

#' plot_bar
#'
#' @noRd
#'
#' @export
plot_bar = function(df, x, y, group, width = 0.6) {

  se = \(vec) sd(vec)/sqrt(length(vec))

  dp = df |>
    select(all_of(c(x, y, group))) |>
    set_names(c('x', 'y', 'group'))

  stat = dp |>
    group_by(x) |>
    wilcox_test(y ~ group) |>
    add_xy_position() |>
    mutate(y = max(y.position),
           label = p2sig(p, 'significance'))

  dp_bar = dp |>
    group_by(x, group) |>
    summarise(mean = mean(y), se = se(y))

  dp_bar |>
    ggplot() +
    geom_bar(aes(x, mean, fill = group), stat = 'identity', position = position_dodge(), width = width, color = 'black') +
    geom_errorbar(aes(x, mean, ymin = mean - se, ymax = mean + se, group = group), position = position_dodge(width), width = 0.25) +
    geom_text(aes(x, max(dp_bar$mean + dp_bar$se), label = label), data = stat, vjust = 0) +
    scale_fill_manual(values = get_color(2)) +
    labs(x = NULL, y = 'Mean Expr.', fill = NULL) +
    theme_classic() +
    theme(legend.position = 'top')

}

## volcano ----

#' part_non
#'
#' @keywords internal
part_non = function(data_, x_, y_, label_, color_ = 'black', alpha_ = 0.75) {

  data_point_ = data_

  x_ = sym(x_)
  y_ = sym(y_)

  list(
    geom_point(data = data_point_, aes(!!x_, !!y_), color = color_, alpha = alpha_)
  )

}

#' part_up
#'
#' @keywords internal
part_up = function(data_, x_, y_, label_, color_ = '#cc0303', max_ = 12, alpha_ = 0.75) {

  data_point_ = data_
  data_label_ = data_[!is.na(data_[[label_]]), ]

  nudge_x_ = max_ - data_label_[[x_]]
  x_ = sym(x_)
  y_ = sym(y_)
  label_ = sym(label_)

  list(
    geom_point(data = data_point_, aes(!!x_, !!y_), color = color_, alpha = alpha_),
    geom_text_repel(
      data = data_label_, aes(!!x_, !!y_, label = !!label_),
      nudge_x = nudge_x_, color = color_,
      max.overlaps = Inf, hjust = 1, vjust = 0,
      direction = 'y',
      segment.linetype = 2, segment.size = 0.1, segment.alpha = 0.5,
      min.segment.length = 0
    )
  )

}

#' part_down
#'
#' @keywords internal
part_down = function(data_, x_, y_, label_, color_ = '#026401', min_ = -12, alpha_ = 0.75) {

  data_point_ = data_
  data_label_ = data_[!is.na(data_[[label_]]), ]

  nudge_x_ = min_ - data_label_[[x_]]
  x_ = sym(x_)
  y_ = sym(y_)
  label_ = sym(label_)

  list(
    geom_point(data = data_point_, aes(!!x_, !!y_), color = color_, alpha = alpha_),
    geom_text_repel(
      data = data_label_, aes(!!x_, !!y_, label = !!label_),
      nudge_x = nudge_x_, color = color_,
      max.overlaps = Inf, hjust = 1, vjust = 0,
      direction = 'y',
      segment.linetype = 2, segment.size = 0.1, segment.alpha = 0.5,
      min.segment.length = 0
    )
  )

}


#' Create a Volcano Plot
#'
#' Combines up-, down-, and non-significant layers into a complete volcano plot.
#'
#' @param df_ A data frame containing `x_`, `y_`, `label_`, and a grouping column (e.g., `colorBy_`).
#' @param x_ Column name for the x-axis variable. Default is "logFC".
#' @param y_ Column name for the y-axis variable. Default is "log10padj".
#' @param label_ Column name for gene labels. Default is "label".
#' @param colorBy_ Column name used for splitting data into "up", "down", and "non" groups. Default is "direction".
#' @param colors_ Named vector of colors for "up", "down", and "non" categories.
#' @param alpha_ Point transparency. Default is 0.75.
#' @param range_ A numeric vector of length 2 specifying the x-axis range. Default is c(-12, 12).
#'
#' @return A ggplot object representing the volcano plot.
#'
#' @importFrom ggrepel geom_text_repel
#'
#' @examples
#' df = data.frame(
#'  logFC = c(-3, 1.5, 2, -4, 0),
#'  log10padj = c(1, 2, 3, 4, 0.5),
#'  label = c("Gene1", "Gene2", NA, "Gene4", NA),
#'  direction = c("down", "up", "up", "down", "non")
#' )
#' plot_volcano(df)
#'
#' @export
plot_volcano = function(
    df_, x_ = 'logFC', y_ = 'log10padj', label_ = 'label',
    colorBy_ = 'direction', colors_ = get_color(3, 'htmap'),
    alpha_ = 0.75, range_ = c(-12, 12)
) {

  lst_dp_ = split(df_, df_[[colorBy_]])

  # browser()
  ggplot() +
    part_non(lst_dp_$non, x_, y_, label_, colors_[[2]], alpha_) +
    part_down(lst_dp_$down, x_, y_, label_, colors_[[1]], range_[[1]], alpha_) +
    part_up(lst_dp_$up, x_, y_, label_, colors_[[3]], range_[[2]], alpha_) +
    xlab('log2(Fold Change)') +
    ylab('-log10(Adj Pvalue)') +
    theme_bw()

}

## enrich ----

### gsea ----

#' Plot Gene Set Enrichment Analysis (GSEA) Results
#'
#' Generates a GSEA plot for a specific term.
#'
#' @param gsea_ A GSEA result object from the `clusterProfiler` package.
#' @param term_ The term to visualize.
#' @param sizef_ A scaling factor for the line thickness. Default is 2.
#' @param color_ A vector of colors for the gradient scale.
#'
#' @return A ggplot object representing the GSEA plot.
#'
#' @export
plot_gsea = function(gsea_, term_, sizef_ = 2, color_ = get_color(2, 'ddlc')) {

  checkReliance('enrichplot')

  geneSetID_ = which(gsea_@result$ID == term_)

  ES_ = gsea_@result$enrichmentScore[geneSetID_]
  FDR_ = gsea_@result$pvalue[geneSetID_]
  tag_ = sprintf("ES=%.3f, False Discovery Rate=%.3f", ES_, FDR_)

  gsdata_ = enrichplot:::gsInfo(gsea_, geneSetID_) |> as_tibble()

  p1_ = gsdata_ |>
    ggplot(aes(x = x)) +
    geom_segment(aes(xend = x, y = 0, yend = -runningScore, color = x), linewidth = 0.1 * sizef_, data = subset(gsdata_, position == 1)) +
    theme_void() +
    scale_x_continuous(expand = c(0, 0)) +
    scale_y_continuous(expand = c(0, 0.1), breaks = NULL) +
    scale_color_gradient(low = color_[[1]], high = color_[[2]]) +
    xlab(term_) +
    ylab(NULL) +
    theme(
      legend.position = "none",
      plot.margin = margin(t = -0.1, b = 0, unit = "cm"),
      axis.title.x = element_text(),
      axis.ticks = element_blank(),
      axis.text = element_blank(),
      axis.line.x = element_blank()
    ) +
    annotation_custom(
      grob = textGrob(
        tag_,
        x = unit(1, "npc"), y = unit(0.05, "npc"),
        hjust = 1, vjust = 0,
        gp = gpar(col = "black", fontsize = 10, family = "italic", fontface = "italic")
      )
    )

  p1_

}

#' Plot Gene Set Enrichment Analysis (GSEA) Results
#'
#' Generates a GSEA plot for a specific term.
#'
#' @param gsea_ A GSEA result object from the `clusterProfiler` package.
#' @param term_ The term to visualize.
#'
#' @return A ggplot object representing the GSEA plot.
#'
#' @export
plot_gsea2 = function(gsea_, term_) {

  geneSetID_ = which(gsea_@result$ID == term_)

  ES_ = gsea_@result$enrichmentScore[geneSetID_]
  FDR_ = gsea_@result$qvalue[geneSetID_]
  tag_ = sprintf("ES=%.3f, False Discovery Rate=%.3f", ES_, FDR_)

  if (ES_ > 0) {

    p1_ = enrichplot:::gseaplot(gsea_, term_)[[2]] +
      xlab(term_) +
      annotation_custom(
        grob = grid::textGrob(
          tag_,
          x = unit(0.95, "npc"), y = unit(0.90, "npc"),
          hjust = 1, vjust = 0,
          gp = gpar(col = "black", fontsize = 12, family = "italic", fontface = "italic")
        )
      )

  } else {

    p1_ = enrichplot::gseaplot(gsea_, term_)[[2]] +
      xlab(term_) +
      annotation_custom(
        grob = grid::textGrob(
          tag_,
          x = unit(0.95, "npc"), y = unit(0.1, "npc"),
          hjust = 1, vjust = 0,
          gp = gpar(col = "black", fontsize = 12, family = "italic", fontface = "italic")
        )
      )

  }

  p1_

}

### pathway ----

#' Visualize KEGG Pathway with Differential Gene Expression Data
#'
#' This function overlays gene expression data onto a KEGG pathway graph.
#'
#' @param deg_ A data frame with gene expression results (e.g., log fold change and adjusted p-values).
#' @param id_ The KEGG pathway ID (e.g., "hsa00010").
#' @param logFC_ Column name for log fold change in `deg_`. Default is "log2FoldChange".
#' @param padj_ Column name for adjusted p-values in `deg_`. Default is "padj".
#' @param color_ A vector of colors for the gradient scale. Default is `get_color(3, 'htmap')`.
#'
#' @return A ggraph object visualizing the KEGG pathway.
#'
#' @export
plot_pathway = function(deg_, id_, logFC_ = 'logFC', padj_ = 'padj', color_ = get_color(3, 'htmap')) {


  graph_ = ggkegg::pathway(id_, 'tmp') |>
    ggraph::activate(nodes) |>
    mutate(converted_name=convert_id(substr(id_, 1, 3))) |>
    left_join(deg_, by = c('converted_name' = 'symbol'))

  ggraph::ggraph(graph_, layout="manual", x = x, y = y)+
    ggraph::geom_node_rect(aes(fill = !!sym(logFC_), filter = !is.na(!!sym(padj_)) & !!sym(padj_)<0.05))+
    # ggfx::with_outer_glow(geom_node_rect(aes(fill=log2FoldChange, filter=!is.na(padj) & padj<0.05)), colour="yellow", expand=2)+
    ggraph::overlay_raw_map(id_, transparent_colors = c("#cccccc","#FFFFFF","#BFBFFF","#BFFFBF"))+
    scale_fill_gradient2(low=color_[[1]], mid=color_[[2]], high=color_[[3]]) +
    theme_void()

}

### classic ----

#' Plot GO Terms
#'
#' This function creates a bar plot for the top GO terms based on their -log10(q-value).
#' It supports grouping the data by a specified feature, such as the ontology type.
#'
#' @param df_ A dataframe, result attr of a GO obj.
#' @param title_ A character string specifying the title of the plot.
#' @param nTerm_ An integer specifying the number of top terms to display (default is 5).
#' @param colorBy_ A character string indicating which feature to group and color the bars by. Defaults to 'ONTOLOGY'.
#' @return A `ggplot` object displaying the bar plot for the top GO terms, grouped and colored by the specified feature.
#'
#' @examples
#' plot_GO(res_GO@result, 'Test')
#'
#' @export
plot_GO = function(df_, title_, nTerm_ = 5, colorBy_ = 'ONTOLOGY') {

  facet_ = as.formula(paste('~', colorBy_))

  colorBy_ = sym(colorBy_)

  dp_ = df_ |>
    group_by(!!colorBy_) |>
    slice_min(qvalue, n = nTerm_, with_ties = F) |>
    mutate(Description = factor(Description, levels = rev(Description))) |>
    ungroup()

  dp_ |>
    ggplot() +
    geom_col(aes(-log10(qvalue), Description, fill = !!colorBy_)) +
    geom_text(aes(0, Description, label = capitalize(Description)), hjust = 0, nudge_x = 0.1, size = 4) +
    geom_vline(xintercept = -log10(0.05), linetype = 2, alpha = 0.5) +
    scale_x_continuous(expand = c(0, 0)) +
    scale_fill_manual(values = get_color(4, 'morandi')) +
    xlab('-Log10(FDR)') +
    ylab(NULL) +
    labs(title = paste(title_, '| GO Terms')) +
    facet_wrap(facet_, ncol = 1, scales = 'free_y', strip.position = 'left') +
    theme_classic() +
    theme(legend.position = 'none',
          axis.ticks.y = element_blank(),
          panel.grid.major = element_line(color = "#fafafa"),
          panel.grid.minor = element_line(color = "#fafafa"),
          axis.text.y = element_blank(),
          axis.line.y = element_blank(),
          plot.title = element_text(hjust = 0.5))

}

#' Plot KEGG Terms
#'
#' This function creates a bar plot for the top KEGG pathway terms
#' based on their -log10(q-value) for a given KEGG result object.
#'
#' @param df_ A dataframe, result attr of a KEGG obj.
#' @param title_ A character string specifying the title of the plot.
#' @param nTerm_ An integer specifying the number of top terms to display (default is 5).
#' @return A `ggplot` object displaying the bar plot for the top KEGG pathway terms.
#'
#' @examples
#' plot_KEGG(res_KEGG@result, 'Test')
#'
#' @export
plot_KEGG = function(df_, title_, nTerm_ = 5) {

  dp_ = df_ |>
    slice_min(qvalue, n = nTerm_, with_ties = F) |>
    mutate(ONTOLOGY = 'KEGG',
           Description = factor(Description, levels = rev(Description)))

  dp_ |>
    ggplot() +
    geom_col(aes(-log10(qvalue), Description, fill = ONTOLOGY)) +
    geom_text(aes(0, Description, label = capitalize(Description)), hjust = 0, nudge_x = 0.1, size = 4) +
    geom_vline(xintercept = -log10(0.05), linetype = 2, alpha = 0.5) +
    scale_x_continuous(expand = c(0, 0)) +
    scale_fill_manual(values = get_color(4, 'morandi')) +
    xlab('-Log10(FDR)') +
    ylab(NULL) +
    labs(title = paste(title_, '| KEGG Terms')) +
    theme_classic() +
    theme(legend.position = 'none',
          axis.ticks.y = element_blank(),
          panel.grid.major = element_line(color = "#fafafa"),
          panel.grid.minor = element_line(color = "#fafafa"),
          axis.text.y = element_blank(),
          axis.line.y = element_blank(),
          plot.title = element_text(hjust = 0.5))

}

#' Plot Gene Ontology (GO) and KEGG Pathway Terms Together
#'
#' This function combines the GO and KEGG results and creates a bar plot for the top terms based on their -log10(q-value).
#'
#' @param go_ A GO result object containing the GO data.
#' @param kegg_ A KEGG result object containing the KEGG pathway data.
#' @param title_ A character string specifying the title of the plot.
#' @param nTerm_ An integer specifying the number of top terms to display (default is 5).
#' @param avoid_ A character string to exclude specific terms from the plot. Default is an empty string, meaning no terms are excluded.
#' @return A `ggplot` object displaying the combined bar plot for the top GO and KEGG terms.
#'
#' @examples
#' plot_GOKEGG(res_GO, res_KEGG, 'Test')
#'
#' @export
plot_GOKEGG = function(go_, kegg_, title_, nTerm_ = 5, avoid_ = '') {

  dp_ = bind_rows(
    go_@result |>
      select(Description, qvalue, ONTOLOGY),
    kegg_@result |>
      mutate(ONTOLOGY = 'KEGG') |>
      select(Description, qvalue, ONTOLOGY)
  )

  if (avoid_ != '') {

    dp_ = dp_ |>
      filter(!map(avoid_, ~ str_detect(Description, .x)) |> pmap(`|`) |> unlist())

  }

  dp_ = dp_ |>
    group_by(ONTOLOGY) |>
    slice_min(qvalue, n = nTerm_, with_ties = F) |>
    mutate(ONTOLOGY = factor(ONTOLOGY, levels = c('BP', 'CC', 'MF', 'KEGG')),
           Description = factor(Description, levels = rev(Description))) |>
    ungroup()

  dp_ |>
    ggplot() +
    geom_col(aes(-log10(qvalue), Description, fill = ONTOLOGY)) +
    geom_vline(xintercept = -log10(0.05), linetype = 2, alpha = 0.5) +
    geom_text(aes(0, Description, label = capitalize(Description)), hjust = 0, nudge_x = 0.1, size = 4) +
    scale_x_continuous(expand = c(0, 0)) +
    scale_fill_manual(values = get_color(4, 'morandi')) +
    xlab('-Log10(FDR)') +
    ylab(NULL) +
    labs(title = title_) +
    facet_wrap(~ ONTOLOGY, ncol = 1, scales = 'free_y', strip.position = 'left') +
    theme_classic() +
    theme(legend.position = 'none',
          axis.ticks.y = element_blank(),
          panel.grid.major = element_line(color = "#fafafa"),
          panel.grid.minor = element_line(color = "#fafafa"),
          axis.text.y = element_blank(),
          axis.line.y = element_blank(),
          plot.title = element_text(hjust = 0.5))

}

### register ----

plot_enrich = new.env()

register(plot_gsea, plot_enrich)
register(plot_pathway, plot_enrich)
register(plot_GO, plot_enrich)
register(plot_KEGG, plot_enrich)
register(plot_GOKEGG, plot_enrich)

## heatmap ----

#' plot_heatmap_basic
#'
#' @noRd
#'
#' @export
plot_heatmap_basic = function(df, x, y, color, colors = get_color(3)) {

  dp_pre = tibble(
    x = df[[x]],
    y = df[[y]],
    color = df[[color]],
  )

  pic_tree = dp_pre |>
    pivot_wider(id_cols = y, names_from = x, values_from = color) |>
    column_to_rownames('y') |>
    dist() |>
    hclust() |>
    ggtree(color = 'white')

  orders = pic_tree$data |>
    filter(isTip) |>
    arrange(y) |>
    pull(label)

  pic_heatmap = dp_pre |>
    mutate(y = factor(y, levels = orders)) |>
    ggplot(aes(x, y, fill = color)) +
    geom_tile(color = 'white') +
    scale_fill_gradientn(
      colors = colors,
      values = scales::rescale(c(-2, 0, 2)),
      limits = c(-2, 2),
      oob = scales::squish
    ) +
    xlab(NULL) +
    ylab(NULL) +
    theme_classic() +
    theme(
      axis.line.x = element_blank(),
      axis.ticks.x = element_blank(),
      plot.title = element_text(hjust = 0.5)
    )

  return(pic_heatmap)

}

#' Plot Heatmap with Gene Set and Class Annotations
#'
#' Generates a heatmap for selected genes with hierarchical clustering and class labels.
#'
#' @param df_ A data frame containing expression data with `symbol` as row names.
#' @param table_x_ A named vector mapping sample names to x labels.
#' @param geneset_ A named list of gene sets to visualize.
#'
#' @return A ggplot object representing the heatmap.
#'
#' @import ggplot2
#' @import ggh4x
#'
#' @export
plot_heatmap = function(df_, table_x_, scale_ = T, color_ = get_color(3), range_ = NULL) {

  checkReliance('ggtree')

  if (scale_) {
    handleRow_ = \(vec__) as.numeric(scale(vec__))
  } else {
    handleRow_ = \(vec__) as.numeric(vec__)
  }

  dp_pre_ = df_ |>
    tidyT('sample', 'symbol') |>
    mutate(across(where(is.numeric), ~ as.numeric(scale(log2(.x + 1))))) |>
    tidyT('symbol', 'sample') |>
    drop_na()

  pic_tree_ = dp_pre_ |>
    column_to_rownames('symbol') |>
    dist() |>
    hclust() |>
    ggtree::ggtree(color = 'white')

  order_ = pic_tree_$data |>
    filter(isTip) |>
    arrange(y) |>
    pull(label)

  dp_ = dp_pre_|>
    pivot_longer(-symbol, names_to = 'sample', values_to = 'expr') |>
    mutate(sample = factor(sample, levels = names(table_x_)),
           class = table_x_[sample],
           symbol = factor(symbol, levels = order_))

  limits = range(range_) %||% dp_[[expr]]

  p_ = dp_ |>
    ggplot(aes(interaction(sample, class), symbol, fill = expr)) +
    geom_tile() +
    scale_fill_gradientn(
      colors = color_,
      values = scales::rescale(c(-2, 0, 2)),
      limits = c(-2, 2),
      oob = scales::squish
    ) +
    xlab(NULL) +
    ylab(NULL) +
    guides(x = "axis_nested") +
    theme_classic() +
    theme(
      axis.text.x = element_text(color = 'white'),
      axis.line.x = element_blank(),
      axis.ticks.x = element_blank(),
      ggh4x.axis.nesttext.x = element_text(color = "black"),
      axis.line.y = element_blank(),
      axis.ticks.y = element_blank(),
      plot.title = element_text(hjust = 0.5)
    )

  return(p_)

}

#' plot_heatmap_align
#'
#' @noRd
#'
#' @export
plot_heatmap_align = function(df, x_class, y_class, color = get_color(3), range = NULL) {

  checkReliance('aplot')

  colnames(df)[[1]] = 'symbol'

  # x
  dp_bar_x = tibble(
    sample = factor(names(x_class), levels = names(x_class)),
    group = x_class
  )

  color_x = c('#cc0303', '#026401')

  p_bar_x = dp_bar_x |>
    ggplot() +
    geom_tile(aes(sample, 1, fill = group)) +
    scale_fill_manual(values = color_x) +
    theme_void() +
    theme(axis.text.x = element_text(angle = 30, vjust = 0.85, hjust = 0.75))

  # y
  dp_bar_y = df |>
    dplyr::select(symbol) |>
    mutate(group = y_class[symbol]) |>
    distinct(symbol, .keep_all = T) |>
    arrange(group) |>
    mutate(symbol = factor(symbol, levels = symbol))

  color_y = get_color(length(unique(dp_bar_y[['group']])), 'macaron')

  p_bar_y = dp_bar_y |>
    ggplot() +
    geom_tile(aes(1, symbol, fill = group)) +
    scale_fill_manual(values = color_y) +
    theme_void() +
    theme(axis.text.y = element_text())

  # body
  dp_body = df |>
    pivot_longer(-symbol, names_to = 'sample', values_to = 'expr') |>
    drop_na() |>
    mutate(sample = factor(sample, levels = names(x_class)),
           class = x_class[sample],
           symbol = factor(symbol, levels = levels(dp_bar_y$symbol)))

  limits = range %||% range(dp_body[['expr']])

  p_body = dp_body |>
    ggplot(aes(sample, symbol, fill = expr)) +
    geom_tile() +
    scale_fill_gradientn(
      colors = color,
      values = scales::rescale(limits),
      limits = limits,
      oob = scales::squish
    ) +
    xlab(NULL) +
    ylab(NULL) +
    theme_classic() +
    theme(
      axis.text.x = element_blank(),
      axis.line.x = element_blank(),
      axis.ticks.x = element_blank(),
      axis.line.y = element_blank(),
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
      plot.title = element_text(hjust = 0.5)
    )

  p = p_body |> aplot::insert_left(p_bar_y, width = 0.1) |> aplot::insert_bottom(p_bar_x, height = 0.05)

  return(p)

}

#' Plot Complex Heatmap with Class Splits
#'
#' Uses the ComplexHeatmap package to generate a detailed heatmap with class splits.
#'
#' @param df_ A data frame containing expression data.
#' @param table_x_ A named vector of class annotations.
#' @param geneset_ A named list of gene sets.
#'
#' @return A Heatmap object from the ComplexHeatmap package.
#'
#' @export
plot_cheatmap = function(df_, table_x_) {

  checkReliance('ComplexHeatmap', 'circlize')

  dp_pre_ = df_ |>
    select(symbol, all_of(names(table_x_))) |>
    tTidy('symbol', 'sample') |>
    mutate(across(where(is.numeric), ~ as.numeric(scale(log2(.x + 1))))) |>
    tTidy('sample', 'symbol') |>
    drop_na() |>
    column_to_rownames('symbol') |>
    as.matrix()

  pic_tree_ = dp_pre_ |>
    dist() |>
    hclust() |>
    ggtree::ggtree(color = 'white')

  order_ = pic_tree_$data |>
    filter(isTip) |>
    arrange(y) |>
    pull(label)

  dp_ = dp_pre_[order_, ]

  col_fun = circlize::colorRamp2(c(-1, 0, 1), get_color(3))

  res_ = ComplexHeatmap::Heatmap(
    dp_,
    name = "expr",
    col = col_fun,
    cluster_columns = FALSE,
    cluster_rows = FALSE,
    show_column_names = FALSE,
    row_names_side = "left",
    column_split = table_x_[colnames(dp_pre_)]
  )

  return(res_)

}

## survival ----

plot_km = function(df, OS, OStime, fea, ...) {

  lst_param = list(...)

  if (is.numeric(df[[fea]])) {

    plot_km_numeric(df, OS, OStime, fea, ...)

  } else {

    plot_km_category(df, OS, OStime, fea)

  }

}

plot_km_category = function(df, OS, OStime, fea) {

  dp = data.frame(
    OS = df[[OS]],
    OStime = df[[OStime]],
    group = df[[fea]]
  ) |> drop_na()

  stat = dp |> with(transGI::cal_HR(group, OStime, OS)) |> signif(3)

  fitSurv = survival::survfit(Surv(event = OS, time = OStime) ~ group, data = dp)

  plotLst = survminer::ggsurvplot(fitSurv, data = dp)

  statistic = paste0(
    paste(names(table(dp$group)), collapse = ':'), ' = ', paste(table(dp$group), collapse = ':'), '\n',
    'Hazard Ratio = ', stat[[1]], '\nLog-Rank Pvalue = ', stat[[4]]
  )

  plotLst$plot +
    annotate('text', label = statistic, x = 0, y = 0, hjust = 0, vjust = 0.1) +
    theme(
      legend.position = 'top'
    )

}

plot_km_numeric = function(df, OS, OStime, fea, strategy = 'best') {

  match.arg(strategy, c('best', 'quantile'))

  dp = data.frame(
    OS = df[[OS]],
    OStime = df[[OStime]],
    fea = df[[fea]]
  ) |> drop_na()

  if (strategy == 'best') {

    cutPoint = surv_cutpoint(dp, time = 'OStime', event = 'OS', variables = 'fea')
    dp$group = factor(ifelse(dp$fea > cutPoint$cutpoint[[1]], "H", "L"), levels = c('L', 'H'))

  } else if (strategy == 'quantile') {

    dp$group = factor(ifelse(dp$fea > quantile(dp$fea, 0.5)[[1]], "H", "L"), levels = c('L', 'H'))

  }

  stat = dp |> with(transGI::cal_HR(group, OStime, OS)) |> signif(3)

  fitSurv = survival::survfit(Surv(event = OS, time = OStime) ~ group, data = dp)

  plotLst = survminer::ggsurvplot(fitSurv, data = dp)

  statistic = paste0(
    paste(names(table(dp$group)), collapse = ':'), ' = ', paste(table(dp$group), collapse = ':'), '\n',
    'Hazard Ratio = ', stat[[1]], '\nLog-Rank Pvalue = ', stat[[4]]
  )

  plotLst$plot +
    annotate('text', label = statistic, x = 0, y = 0, hjust = 0, vjust = 0.1) +
    theme(
      legend.position = 'top'
    )

}

## summary ----

pdf2png = function(from_, to_, patterns_ = NULL, prefix_ = 'none') {

  match.arg(prefix_, c('none', 'dirname'))

  path_plot_ = list.files(from_, pattern = '.pdf', full.names = T, recursive = T)

  if (!is.null(patterns_)) {
    for (pattern_ in patterns_) {
      path_plot_ = path_plot_[str_detect(path_plot_, pattern_)]
    }
  }

  if (length(path_plot_) == 0) stop()

  if (prefix_ == 'none') {

    walk2(
      path_plot_, path_plot_ |> tools::file_path_sans_ext() |> basename(),
      ~ pdftools::pdf_convert(.x, filenames = file.path(to_, paste0(.y, '.png')), format = 'png', dpi = 300)
    )

  } else if(prefix_ == 'dirname') {

    walk(
      path_plot_,
      \(rawPath__) {
        toConvert__ = R.utils::getRelativePath(normalizePath(rawPath__), from_)
        converted__ = file.path(to_, gsub('\\/', '_', paste0(tools::file_path_sans_ext(toConvert__), '.png')))

        pdftools::pdf_convert(rawPath__, filenames = converted__, format = 'png', dpi = 300)
      }
    )

  }

}
