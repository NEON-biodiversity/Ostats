### Code to test Ostats_multivariate_plot

# V1. Points only, no volumes ---------------------------------------------



# Panel widths and heights
units <- 'cm'
panel_height <- grid::unit(3, units = units)
panel_width <- grid::unit(3, units = units)

# Plot only the points.

sp <- iris$Species
traits <- iris[, 1:4]
plots <- factor(rep(1:3, nrow(traits)/3))

trait_combs <- combn(names(traits), 2) # All combinations of traits

# Create triangular layout
layout_mat <- matrix(as.numeric(NA), ncol(traits) - 1, ncol(traits) - 1)
layout_mat[lower.tri(layout_mat, diag = TRUE)] <- 1:ncol(trait_combs)
layout_mat <- layout_mat[nrow(layout_mat):1, ncol(layout_mat):1]
layout_mat[nrow(layout_mat), 1] <- ncol(trait_combs) + 1 # Location of legend

# If no color values are provided, produce default colors.
if (is.null(colorvalues)) {
  colorvalues <- sample(viridis::viridis(length(unique(sp)), alpha = 1))
}

sp_names <- rev(sort(unique(sp)))
color_scale <- ggplot2::scale_color_manual(values = setNames(colorvalues, sp_names))

plot_list <- list()

# Generate common legend for all plots by writing species names in different colors.

legend_panel <- ggplot2::ggplot(data.frame(x = 1, y = seq_along(sp_names), sp = sp_names),
                                ggplot2::aes(x = x, y = y, label = sp, color = sp)) +
  ggplot2::geom_text(hjust = 0) +
  color_scale +
  ggplot2::scale_y_continuous(expand = c(0.5, 0.5)) +
  ggplot2::theme_void() +
  ggplot2::theme(legend.position = 'none')

for (p in unique(plots)) {

  sp_plot <- sp[plots == p]
  traits_plot <- traits[plots == p, ]

  # Find ranges of each trait so that panels will have a common range.
  traits_range <- apply(traits_plot, 2, range)

  plot_theme <- ggplot2::theme_bw() +
    ggplot2::theme(legend.position = 'none')

  # Generate plots
  trait_pairs_plot_list <- apply(trait_combs, 2, function(traits_to_plot) {
    dat_plot <- data.frame(sp = sp_plot, x = traits_plot[, traits_to_plot[1]], y = traits_plot[, traits_to_plot[2]])
    x_range <- traits_range[, traits_to_plot[1]]
    y_range <- traits_range[, traits_to_plot[2]]
    ggplot2::ggplot(dat_plot, ggplot2::aes(x = x, y = y, group = sp, color = sp)) +
      ggplot2::geom_point() +
      color_scale +
      plot_theme +
      ggplot2::labs(x = traits_to_plot[1], y = traits_to_plot[2])
  })

  # Remove axis text and titles from plots not along the edge.
  for (i in layout_mat[upper.tri(layout_mat)]) {
    trait_pairs_plot_list[[i]] <- trait_pairs_plot_list[[i]] + ggplot2::theme(axis.text.x = ggplot2::element_blank(),
                                                                              axis.text.y = ggplot2::element_blank(),
                                                                              axis.title.x = ggplot2::element_blank(),
                                                                              axis.title.y = ggplot2::element_blank())
  }

  # Add legend panel to plot list
  trait_pairs_plot_list <- c(trait_pairs_plot_list, list(legend_panel))

  # Arrange plots
  # Sequentially bind the rows together, then the columns.
  trait_pairs_plots_rows <- list()

  # Create dummy grob to fill in the blank spaces, using a zeroGrob
  dummy_grob <- ggplot2::ggplotGrob(trait_pairs_plot_list[[1]])
  dummy_grob$grobs <- replicate(length(dummy_grob$grobs), ggplot2::zeroGrob(), simplify = FALSE)


  for (i in 1:nrow(layout_mat)) {
    trait_pairs_plots_rows[[i]] <- do.call(gridExtra::gtable_cbind, lapply(layout_mat[i, ], function(n) if (is.na(n)) dummy_grob else ggplot2::ggplotGrob(trait_pairs_plot_list[[n]])))
  }

  trait_pairs_plots_arranged <- do.call(gridExtra::gtable_rbind, trait_pairs_plots_rows)

  plot_list[[length(plot_list) + 1]] <- trait_pairs_plots_arranged

}


# V2. Include volumes -----------------------------------------------------

# Calculate hypervolumes
library(hypervolume)

# Hypervolume for each species/site combination.
sp_plot_combs <- expand.grid(sp = unique(sp), plot = unique(plots))
hv_list <- replicate(nrow(sp_plot_combs), NA, simplify = FALSE)

for (i in 1:nrow(sp_plot_combs)) {
  sp_plot_dat <- traits[sp == sp_plot_combs$sp[i] & plots == sp_plot_combs$plot[i], ]
  if (nrow(sp_plot_dat) > 2) {
    sp_plot_dat <- scale(sp_plot_dat)
    hv_list[[i]] <- hypervolume::hypervolume(sp_plot_dat, method = 'gaussian')
  }
}

# Code to draw contours somewhat modified from hypervolume::plot.HypervolumeList
get_contours <- function(hv) {
  hv_density <- nrow(hv@RandomPoints)/hv@Volume
  hv_dimensionality <- hv@Dimensionality
  radius_critical <- hv_density^(-1/hv_dimensionality)
  # Calculate kernel density estimate for each combinations of two variables.
  contour_list <- list()
  for (i in 1:ncol(trait_combs)) {
    kde <- MASS::kde2d(hv@RandomPoints[, trait_combs[1, i]], hv@RandomPoints[, trait_combs[2, i]], n = 50, h = radius_critical, lims = c(range(hv@RandomPoints[, trait_combs[1, i]]) * c(0.9, 1.1), range(hv@RandomPoints[, trait_combs[2, i]]) * c(0.9, 1.1)))
    contour_lines <- contourLines(kde, levels = 0.01)
    contour_line_dfs <- list()
    for (j in 1:length(contour_lines)) {
      contour_line_dfs[[j]] <- with(contour_lines[[j]], data.frame(polygon_id = j, x = x, y = y))
    }
    contour_list[[i]] <- data.frame(trait_x = trait_combs[1, i], trait_y = trait_combs[2, i], do.call(rbind, contour_line_dfs))
  }
  do.call(rbind, contour_list)
}

# c_df1 <- get_contours(hv_list[[1]])

# this works:
# df1<-data.frame(x=c_list[[1]][[1]]$x,y=c_list[[1]][[1]]$y)
# ggplot(df1,aes(x=x,y=y))+geom_polygon()
# ggplot(df1,aes(x=x,y=y))+geom_polygon(alpha=0.5)

###################

# New loop to create plot list, now with hypervolumes

plot_list <- list()

for (p in unique(plots)) {

  sp_plot <- sp[plots == p]
  traits_plot <- traits[plots == p, ]

  plot_theme <- ggplot2::theme_bw() +
    ggplot2::theme(legend.position = 'none')

  # Generate hypervolumes for each species at this plot.
  # Hypervolume for each species/site combination.
  sp_in_plot <- unique(sp_plot)
  hv_list <- replicate(length(sp_in_plot), NA, simplify = FALSE)

  for (i in 1:length(sp_in_plot)) {
    sp_plot_dat <- traits_plot[sp_plot == sp_in_plot[i], ]
    # Only generate hypervolume with at least 3 measurements
    # Generate UNSCALED hypervolume.
    if (nrow(sp_plot_dat) > 2) {
      hv_list[[i]] <- hypervolume::hypervolume(sp_plot_dat, method = 'gaussian')
    }
  }

  # Generate contours for all hypervolumes
  contours_list <- lapply(hv_list, function(hv) if (class(hv) == 'Hypervolume') get_contours(hv) else NA)
  # Join contours to data frame
  for (i in 1:length(contours_list)) contours_list[[i]][, 'sp'] = sp_in_plot[i]
  contours_df <- do.call(rbind, contours_list)

  # Generate plots, with contours
  trait_pairs_plot_list <- list()
  for (i in 1:ncol(trait_combs)) {

    dat_points <- data.frame(sp = sp_plot, x = traits_plot[, trait_combs[1, i]], y = traits_plot[, trait_combs[2, i]])
    dat_polygons <- contours_df[contours_df$trait_x == trait_combs[1, i] & contours_df$trait_y == trait_combs[2, i], ]

    trait_pairs_plot_list[[i]] <- ggplot2::ggplot(dat_points, ggplot2::aes(x = x, y = y, group = sp, color = sp)) +
      ggplot2::geom_polygon(data = dat_polygons, ggplot2::aes(group = interaction(sp, polygon_id)), fill = 'transparent') +
      ggplot2::geom_point() +
      ggplot2::scale_x_continuous(limits = range(contours_df$x[contours_df$trait_x == trait_combs[1, i]])) +
      ggplot2::scale_y_continuous(limits = range(contours_df$y[contours_df$trait_y == trait_combs[2, i]])) +
      color_scale +
      plot_theme +
      ggplot2::labs(x = trait_combs[1, i], y = trait_combs[2, i])
  }

  # Remove axis text and titles from plots not along the edge.
  for (i in layout_mat[upper.tri(layout_mat)]) {
    trait_pairs_plot_list[[i]] <- trait_pairs_plot_list[[i]] + ggplot2::theme(axis.text.x = ggplot2::element_blank(),
                                                                              axis.text.y = ggplot2::element_blank(),
                                                                              axis.title.x = ggplot2::element_blank(),
                                                                              axis.title.y = ggplot2::element_blank())
  }

  # Add legend panel to plot list
  trait_pairs_plot_list <- c(trait_pairs_plot_list, list(legend_panel))

  # Arrange plots
  # Sequentially bind the rows together, then the columns.
  trait_pairs_plots_rows <- list()

  # Create dummy grob to fill in the blank spaces, using a zeroGrob
  dummy_grob <- ggplot2::ggplotGrob(trait_pairs_plot_list[[1]])
  dummy_grob$grobs <- replicate(length(dummy_grob$grobs), ggplot2::zeroGrob(), simplify = FALSE)


  for (i in 1:nrow(layout_mat)) {
    trait_pairs_plots_rows[[i]] <- do.call(gridExtra::gtable_cbind, lapply(layout_mat[i, ], function(n) if (is.na(n)) dummy_grob else ggplot2::ggplotGrob(trait_pairs_plot_list[[n]])))
  }

  trait_pairs_plots_arranged <- do.call(gridExtra::gtable_rbind, trait_pairs_plots_rows)

  plot_list[[length(plot_list) + 1]] <- ggpubr::as_ggplot(trait_pairs_plots_arranged)

}

names(plot_list) <- unique(plots)

if (length(plot_list) == 1) {
  return(plot_list[[1]])
} else {
  return(plot_list)
}


# Test the actual function. -----------------------------------------------

library(Ostats)

sp <- iris$Species
traits <- iris[, 1:4]
plots <- factor(rep(1:3, nrow(traits)/3))

overlap_iris <- Ostats_multivariate(traits = as.matrix(traits), plots = plots, sp = sp, nperm = 1, random_seed = 1)

p_list <- Ostats_multivariate_plot(plots = plots, sp = sp, traits = traits, overlap_dat = overlap_iris, use_plots = NULL, colorvalues = c('red','blue','green'), plot_points = FALSE)


#map(1:3, ~ as.character(glue::glue('~/Documents/temp/plot{.}.png')))
walk(1:3, ~ ggsave(as.character(glue::glue('~/Documents/temp/plot{.}.png')), p_list[[.]]))


# for debugging.
overlap_dat = overlap_iris; use_plots = NULL; colorvalues = c('red','blue','green'); plot_points = FALSE
panel_height = 3; panel_width = 3; units = "cm"
hypervolume_args <- list()
get_contours = Ostats:::get_contours
