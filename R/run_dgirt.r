#' Estimate dynamic group-level IRT model
#'
#' @param dgirt_data Data prepared for use with \code{run_dgirt} by \code{format_data}.
#' @param n_iter See \code{iter} in \code{rstan::stan}.
#' @param n_chain See \code{chains} in \code{rstan::stan}.
#' @param n_warm See \code{warmup} in \code{rstan::stan}.
#' @param max_save Maximum iterations to save; only used in the default value of n_thin.
#' @param n_thin See \code{thin} in \code{rstan::stan}.
#' @param init_range See \code{init} in \code{rstan::stan}.
#' @param seed See \code{seed} in \code{rstan::stan}.
#' @param save_pars See \code{pars} in \code{rstan::rstan_options}
#' @param parallel See \code{rstan::rstan_options(auto_write = parallel)}.
#' @param method By default, `rstan::stan` estimates the model using MCMC
#'        sampling. Alternatively, `cmdstan optimize` or `cmdstan variational`
#'        can be used if `CmdStan` is available. Note that these methods
#'        are faster than MCMC sampling but return only point estimates.
#'        See \url{http://mc-stan.org/interfaces/cmdstan.html} for `CmdStan`
#'        installation instructions.
#' @export
run_dgirt <- function(dgirt_data, n_iter = 2000, n_chain = 2, max_save = 2000, n_warm = min(10000, 
  floor(n_iter * 3/4)), n_thin = ceiling((n_iter - n_warm)/(max_save/n_chain)), 
  init_range = 1, seed = 1, save_pars = c("theta_bar", "xi", "gamma", "delta_gamma", 
    "delta_tbar", "nu_geo", "nu_geo_prior", "kappa", "sd_item", "sd_theta", "sd_theta_bar", 
    "sd_gamma", "sd_innov_gamma", "sd_innov_delta", "sd_innov_logsd", "sd_total", 
    "theta_l2", "var_theta_bar_l2"), parallel = TRUE, method = c("rstan", "optimize", "variational")) {

  requireNamespace("rstan", quietly = TRUE)
  rstan::rstan_options(auto_write = parallel)
  if (parallel) {
    options(mc.cores = parallel::detectCores())
  }

  message("Started:", date())
  if (identical(method, "rstan") || identical(method, c("rstan", "optimize", "variational"))) {
    message("Running", n_iter, "iterations in each of ", n_chain, "chains. Thinning at an interval of", 
      n_thin, "with", n_warm, "adaptation iterations.")
    stan_out <- rstan::stan(model_code = stan_code, data = dgirt_data, iter = n_iter,
      chains = n_chain, warmup = n_warm, thin = n_thin, verbose = FALSE, pars = save_pars,
      seed = seed, init = "random", init_r = init_range)
  } else if (identical(method, "optimize") || identical(method, "variational")) {
    stan_out <- run_cmdstan(method, n_iter, init_range)
  } else {
    stop("Didn't recognize run_dgirt method")
  }
  message("Ended:", date())
  return(stan_out)
}

run_cmdstan = function(dgirt_data, method, n_iter, init_range) {
  dump_dgirt(dgirt_data)
  stan_call <- paste0(get_dgirt_path(), " ", method, " iter=", n_iter,
    " init='", init_range, "' data file=", get_dump_path())
  system(stan_call)
  unlink(get_dump_path())
  if (file.exists(get_output_path())) {
    stan_output <- read_cmdstan_output()
    return(stan_output)
  } else {
    warning("cmdstan didn't write an output file; check its output for errors.")
    return(NULL)
  }
}

read_cmdstan_output = function() {
    output_path  <- get_output_path()
    cmdstan_output <- readLines(output_path)
    cmdstan_config <- cmdstan_output[stringr::str_sub(cmdstan_output, 1, 1) == '#']
    message("Reading sampled values from disk. (This may take some time.)")
    cmdstan_values <- read.csv(output_path, skip = length(cmdstan_config))
    unlink(outout_path)
    return(list(config = cmdstan_config, values = cmdstan_values))
}

dump_dgirt <- function(dgirt_data) {
  stopifnot(is.list(dgirt_data))
  stopifnot(length(dgirt_data) > 0)
  rstan::stan_rdump(names(dgirt_data), get_dump_path(),
    envir = list2env(dgirt_data))
}

get_dgirt_path <- function() {
  system.file("dgirt", package = "dgirt", mustWork = TRUE)
}

get_output_path <- function() {
  paste0(system.file(package = "dgirt"), "/output.csv")
}

get_dump_path <- function() {
  paste0(system.file(package = "dgirt"), "/dgirt_data.Rdump")
}