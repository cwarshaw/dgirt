#' Extract and name parameters
#'
#' `name ` is a wrapper for rstan::extract that attaches names to parameters
#' using the data originally passed to `dgirt`.
#' @param dgirt_output Return value of `dgirt`, a `stanfit` object.
#' @param dgirt_input Return value of `wrangle`.
#' @return What `rstan::extract` returns, but with dimension names.
#' @export
name <- function(dgirt_output, dgirt_input) {
  assertthat::assert_that(inherits(dgirt_output, "stanfit"))
  assertthat::assert_that(assertthat::not_empty(dgirt_input$vars))

  dgirt_extract <- rstan::extract(dgirt_output)

  vars <- dgirt_input$vars
  dim2_indexed_t <- c('theta_bar', 'xi', 'gamma', 'delta_gamma', 'delta_tbar', 'nu_geo', 'sd_theta', 'sd_theta_bar',
    'sd_total', 'theta_l2', 'var_theta_bar_l2')
  if (!as.logical(dgirt_input$constant_item)) dim2_indexed_t <- c(dim2_indexed_t, "kappa")

  for (i in dim2_indexed_t) {
    names(attributes(dgirt_extract[[i]])$dimnames)[2] <- 'time'
    assertthat::assert_that(identical(dim(dgirt_extract[[i]])[2], length(vars$use_t)))
    dimnames(dgirt_extract[[i]])[[2]] <- vars$use_t
  }

  names(attributes(dgirt_extract[['theta_bar']])$dimnames)[3] <- 'group'
  groups_concat = concat_groups(vars$covariate_groups, vars$groups, vars$geo_id, "groups")$groups
  assertthat::assert_that(identical(dim(dgirt_extract[['theta_bar']])[3], length(groups_concat)))
  dimnames(dgirt_extract[['theta_bar']])[[3]] <- groups_concat

  names(attributes(dgirt_extract[['gamma']])$dimnames)[3] <- 'param'
  assertthat::assert_that(identical(dim(dgirt_extract[['gamma']])[3], length(vars$hier_names)))
  dimnames(dgirt_extract[['gamma']])[[3]] <- vars$hier_names

  names(attributes(dgirt_extract[['kappa']])$dimnames)[3] <- 'item'
  assertthat::assert_that(identical(dim(dgirt_extract[['kappa']])[3], length(vars$gt_items)))
  dimnames(dgirt_extract[['kappa']])[[3]] <- vars$gt_items

  names(attributes(dgirt_extract[['sd_item']])$dimnames)[2] <- 'item'
  assertthat::assert_that(identical(dim(dgirt_extract[['sd_item']])[2], length(vars$gt_items)))
  dimnames(dgirt_extract[['sd_item']])[[2]] <- vars$gt_items

  dgirt_extract <- lapply(dgirt_extract, reshape2::melt)
  dgirt_extract <- lapply(dgirt_extract, dplyr::rename, iteration = iterations)

  dgirt_extract
}
