#' Convert O-stats Output to Long Form
#'
#' This function converts the output of Ostats to long form.
#'
#' @param o the output of a call to \code{Ostats}.
#'
#' @return a dataframe with nrows = the names of communities and ncol =
#' outputs from Ostats function.
#'
#' @export
Ostat2longform <- function(o) {
  result_names <- c('site','trait','ostat_norm','ostat_norm_localnull_lower','ostat_norm_localnull_upper','ostat_norm_localnull_ses','ostat_norm_localnull_seslower','ostat_norm_localnull_sesupper','ostat_unnorm','ostat_unnorm_localnull_lower','ostat_unnorm_localnull_upper','ostat_unnorm_localnull_ses','ostat_unnorm_localnull_seslower','ostat_unnorm_localnull_sesupper')

  res_list <- list()

  nsite <- nrow(o[[1]])
  ntrait <- ncol(o[[1]])

  for (i in 1:nsite) {
    for (j in 1:ntrait) {
      res_list[[length(res_list)+1]] <- c(o$overlaps_norm[i,j],
                                          o$overlaps_norm_ses$raw_lower[i,j],
                                          o$overlaps_norm_ses$raw_upper[i,j],
                                          o$overlaps_norm_ses$ses[i,j],
                                          o$overlaps_norm_ses$ses_lower[i,j],
                                          o$overlaps_norm_ses$ses_upper[i,j],
                                          o$overlaps_unnorm[i,j],
                                          o$overlaps_unnorm_ses$raw_lower[i,j],
                                          o$overlaps_unnorm_ses$raw_upper[i,j],
                                          o$overlaps_unnorm_ses$ses[i,j],
                                          o$overlaps_unnorm_ses$ses_lower[i,j],
                                          o$overlaps_unnorm_ses$ses_upper[i,j])
    }
  }
  res <- as.data.frame(do.call('rbind', res_list))
  res <- cbind(site = rep(dimnames(o[[1]])[[1]], each=ntrait), trait = rep(dimnames(o[[1]])[[2]], times=nsite), res)
  names(res) <- result_names
  res
}
