
plot_loadings <- function(L_est, Pop, legendYN = TRUE, scales_option = "fixed"){
  # options for scales_option are "fixed" and "free_y"
  n <- nrow(L_est)
  k <- ncol(L_est)
  Idx <- rep(c(1:n), k)
  Loading <- c(L_est)
  Factor <- paste0('k=',c(sapply(c(1:k), function(x, n){rep(x, n)}, n = n)))
  tib <- data.frame(Idx, Loading, Factor, Pop)
  plt <- ggplot(tib, aes(x = Idx, y = Loading, col = Pop)) +
    geom_point(show.legend = legendYN) +
    geom_hline(yintercept = 0, linetype = "dashed") +
    facet_grid(cols = vars(Factor), scales = scales_option) +
    theme(axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          axis.ticks.x = element_blank(),
          axis.ticks.y = element_blank(),
          axis.text.x = element_blank(),
          axis.text.y = element_blank(),
          panel.spacing = unit(1, "lines"))
  plot(plt)
}

plot_heatmap <- function(L, title = "", colors_range = c("gray96", "red"), brks = NULL){
  ### define the color map
  cols <- colorRampPalette(colors_range)(49)
  if (is.null(brks) == TRUE){
    brks <- seq(min(L), max(L), length=50)
  }
  
  plt <- pheatmap(L, show_rownames = FALSE, show_colnames = FALSE, cluster_rows = FALSE, cluster_cols = FALSE, color = cols, breaks = brks, main = title)
  return(plt)
}