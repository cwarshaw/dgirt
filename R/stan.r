## Last edited 4-Mar-2016
stan_code <- "
data {
  int<lower=1> G; ## number of covariate groups
  int<lower=1> Gl2; ## number of level-two demographic groups
  int<lower=1> Q; ## number of items/questions
  int<lower=1> T; ## number of years
  int<lower=1> N; ## number of observed cells
  int<lower=1> S; ## number of geographic units (e.g., states)
  int<lower=1> P; ## number of hierarchical parameters, including geographic
  int<lower=1> H; ## number of predictors for geographic unit effects
  int<lower=1> Hprior; ## number of predictors for geographic unit effects (t=1)
  int<lower=1> Tdiff; ## number of difficulty parameters per question
  int<lower=1> Tdisc; ## number of discrimination parameters per question
  int<lower=1> Tkappa; ## number of threshold parameters per question
  int<lower=0,upper=1> constant_diff; ## indicator for constant difficulties
  int<lower=0,upper=1> constant_disc; ## indicator for constant discriminations
  int<lower=0,upper=1> separate_t; ## indicator for no over-time smoothing
  real delta_tbar_prior_mean;
  real<lower=0> delta_tbar_prior_sd;
  real<lower=0> scale_sd_innov_delta;
  real<lower=0> scale_sd_innov_gamma;
  real<lower=0> scale_sd_innov_item;
  real<lower=0> scale_sd_innov_logsd;
  int n_vec[N]; ## long vector of trials
  int s_vec[N]; ## long vector of successes
  int NNl2[T, Q, Gl2]; ## trials
  int SSl2[T, Q, Gl2]; ## successes
  int<lower=0> MMM[T, Q, G]; ## missingness array
  matrix<lower=0, upper=1>[G, P] XX; ## indicator matrix for hierarchical vars.
  matrix<lower=0, upper=1>[Gl2, G] WT[T]; ## weight array
  matrix[P, H] ZZ[T]; ## data for geographic model
  matrix[P, Hprior] ZZ_prior[T]; ## data for geographic model (prior)
  matrix<lower=0, upper=1>[T, Q] l2_only;
}
transformed data {
}
parameters {
  vector[Q] diff_raw[Tdiff]; ## raw difficulty
  vector<lower=0>[Q] disc_raw[Tdisc]; ## raw discrimination
  vector[T] xi; ## common intercept
  vector[P] gamma_raw[T]; ## hierarchical parameters (raw)
  vector[T] delta_gamma; ## weight placed on gamma from prev. period
  vector[H] nu_geo[T]; ## weight on geographic predictors
  vector[Hprior] nu_geo_prior; ## weight on geographic predictors (t=1)
  vector[T] delta_tbar; ##
  vector[G] theta_bar_raw[T]; ## group mean ability (raw) #!#
  #!# vector[G] theta_bar[T]; ## group means
  vector<lower=0>[T] sd_theta_bar; ## residual sd of group ability means
  vector<lower=0>[T] sd_theta; ## sd of abilities (by period)
  real<lower=0> sd_gamma_geo; ## prior sd of geographic coefficients
  real<lower=0> sd_gamma_demo; ## prior sd of demographic coefficients
  real<lower=0> sd_innov_delta; ## innovation sd of nu_geo and delta_gamma
  real<lower=0> sd_innov_logsd; ## innovation sd of sd_theta
  real<lower=0> sd_innov_gamma; ## innovation sd of gamma and xi
  real<lower=0> sd_innov_item; ## innovation sd diff/disc (optional)
}
transformed parameters {
  vector[G] theta_bar[T]; ## group means (transformed) #!#
  vector[Q] diff[Tdiff]; ## adjusted difficulty
  vector[Q] kappa[Tkappa]; ## threshold
  vector<lower=0>[Q] disc[Tdisc]; ## normalized discrimination
  vector<lower=0>[Q] sd_item[Tdisc]; ## item standard deviation
  vector<lower=0>[Q] var_item[Tdisc]; ## item variance
  vector<lower=0>[T] var_theta; ## within-group variance of theta
  ## var. of theta_bar w/in each level-two group **NOT CONSTRAINED TO BE POSITIVE**
  vector[Gl2] var_theta_bar_l2[T];
  vector[P] gamma[T]; ## hierarchical parameters (adjusted)
  vector[G] mu_theta_bar[T]; ## linear predictor for group means
  vector[P] mu_gamma[T];
  vector[G] z[T, Q]; ## array of vectors of group deviates
  vector[Gl2] z_l2[T, Q]; ##
  real<lower=0,upper=1> prob[T, Q, G]; ## array of probabilities
  vector[Gl2] prob_l2[T, Q]; ## array of probabilities
  vector[Gl2] theta_l2[T]; ## second-level group abililities
  ## Constant item parameters
  if (constant_diff == 1 && constant_disc == 1) {
    ## identify location (mean = 0)
    diff[1] <- diff_raw[1] - mean(diff_raw[1]); 
    ## identify scale (product = 1)
    disc[1] <- disc_raw[1] * pow(exp(sum(log(disc_raw[1]))), (-inv(Q)));
    sd_item[1] <- inv(disc[1]);
    var_item[1] <- sd_item[1] .* sd_item[1];
    kappa[1] <- diff[1] ./ disc[1];
  }
  ## Evolving difficulty parameters
  if (constant_diff == 0 && constant_disc == 1) { 
    for (d in 1:Tdiff) {
      ## identify location (mean in first year = 0)
      diff[d] <- diff_raw[d] - mean(diff_raw[1]);
      ## identify scale (product = 1)
      disc[1] <- disc_raw[1] * pow(exp(sum(log(disc_raw[1]))), (-inv(Q)));
      sd_item[1] <- inv(disc[1]);
      var_item[1] <- sd_item[1] .* sd_item[1];
      kappa[d] <- diff[d] ./ disc[1]; 
    }
  }
  ## Evolving discrimination parameters
  if (constant_diff == 1 && constant_disc == 0) { 
    for (d in 1:Tdiff) {
      ## identify location (mean = 0)
      diff[1] <- diff_raw[1] - mean(diff_raw[1]); 
      ## identify scale (product in first year = 1)
      disc[d] <- disc_raw[d] * pow(exp(sum(log(disc_raw[1]))), (-inv(Q)));
      sd_item[d] <- inv(disc[d]);
      var_item[d] <- sd_item[d] .* sd_item[d];
      kappa[d] <- diff[1] ./ disc[d];
    }
  }
  var_theta <- sd_theta .* sd_theta;
  for (t in 1:T) { ## loop over years
    if (t == 1 || separate_t == 1) {
      mu_gamma[t] <- ZZ_prior[t] * nu_geo_prior;
      for (p in 1:P) {
        if (p <= S) { ## if geographic coefficient
          gamma[t][p] <- mu_gamma[t][p] + sd_gamma_geo * gamma_raw[t][p];
        }
        if (p > S) { ## if demographic coefficient
          gamma[t][p] <- mu_gamma[t][p] + sd_gamma_demo * gamma_raw[t][p];
        }
      }
      mu_theta_bar[t] <- xi[t] + XX * gamma[t];
      ##mu_theta_bar[t] <- XX * gamma[t];
    }
    if (t > 1 && separate_t == 0) {
      if (t == 2) {
        ## 2016-02-05: need to think more about nu_geo_prior; make it different
        ## for geographic and demographic parameters.
        ##
        ## In the second year, again use uniformative prior for gamma, rather
        ## than one centered on its lagged value, because gamma is likely to be
        ## very different in periods 1 and 2 because only in 2 is
        ## theta_bar[t - 1] used to inform theta_bar[t].
        mu_gamma[t] <- ZZ_prior[t] * nu_geo_prior;
        for (p in 1:P) {
          if (p <= S) { ## if geographic coefficient
            gamma[t][p] <- mu_gamma[t][p] + sd_gamma_geo * gamma_raw[t][p];
          }
          if (p > S) { ## if demographic coefficient
            gamma[t][p] <- mu_gamma[t][p] + sd_gamma_demo * gamma_raw[t][p];
          }
        }
      } else {
        ## 2016-02-05: maybe delta_gamma should differ for geographic and
        ## demographic parameters; could also do random walk DLM for demographic
        ## parameters
        for (p in 1:P) {
          if (p <= S) { ## if geographic coefficient
            mu_gamma[t][p] <- gamma[t - 1][p]*delta_gamma[t] + ZZ[t][p]*nu_geo[t];
          }
          if (p > S) { ## if demographic coefficient
            mu_gamma[t][p] <- gamma[t - 1][p]; ## random walk
          }
        }
        gamma[t] <- mu_gamma[t] + sd_innov_gamma * gamma_raw[t];
      }
      mu_theta_bar[t] <- xi[t] + XX * gamma[t] + theta_bar[t - 1] * delta_tbar[t];
      ##mu_theta_bar[t] <- theta_bar[t - 1] * delta_tbar[t] + XX * gamma[t];
    }
    ## Matt trick for group means
    theta_bar[t] <- mu_theta_bar[t] + sd_theta_bar[t] * theta_bar_raw[t]; #!#
    ## Weighted average of group means (weights must sum to 1)
    theta_l2[t] <- WT[t] * theta_bar[t]; ## Gl2x1 = Gl2xG * Gx1
    for (n in 1:Gl2) {
      matrix[G, G] WTdiag;
      for (g in 1:G) {
        for (h in 1:G) {
          if (g == h) {
            WTdiag[g, h] <- WT[t][n][g];
          }
          if (g != h) {
            WTdiag[g, h] <- 0;
          }
        }
      }
      ## (y - w'y)' W (y - w'y) = weighted variance
      var_theta_bar_l2[t][n] <- quad_form(WTdiag, theta_bar[t] - theta_l2[t, n]);
      }
      for (q in 1:Q) { ## loop over questions
        real sd_tq;
        real sd_l2_tq[Gl2];
        ## Constant discrimination
        if (constant_disc == 1) {
          sd_tq <- sqrt(var_theta[t] + var_item[1][q]);
        }
        ## Evolving discrimination
        if (constant_disc == 0) {
          sd_tq <- sqrt(var_theta[t] + var_item[t][q]);
        }
        for (n in 1:Gl2) {
          sd_l2_tq[n] <- sqrt(square(sd_tq) + var_theta_bar_l2[t, n]);
        }
        ## Constant item parameters
        if (constant_diff == 1 && constant_disc == 1) {
          z[t, q] <- (theta_bar[t] - kappa[1][q]) / sd_tq;
          for (n in 1:Gl2) {
            z_l2[t, q, n] <- (theta_l2[t, n] - kappa[1][q]) / sd_l2_tq[n];
            prob_l2[t, q, n] <- Phi_approx(z_l2[t, q, n]);
          }
        }
        ## Evolving difficulty or discrimination
        if (constant_diff == 0 || constant_disc == 0) {
          z[t, q] <- (theta_bar[t] - kappa[t][q]) / sd_tq;
          for (n in 1:Gl2) {
            z_l2[t, q, n] <- (theta_l2[t, n] - kappa[t][q]) / sd_l2_tq[n];
            prob_l2[t, q, n] <- Phi_approx(z_l2[t, q, n]);
          }
        }
        for (g in 1:G) { ## loop over groups
          prob[t, q, g] <- Phi_approx(z[t, q, g]); ## fast normal CDF
        }
      } ## end question loop
    } ## end year loop
  }
  model {
    ## TEMPORARY VARIABLES
    real prob_vec[N]; ## long vector of probabilities (empty cells omitted)
    int pos;
    pos <- 0;
    ## PRIORS
    diff_raw[1] ~ normal(0, 1); ## difficulty (first period, if evolving)
    disc_raw[1] ~ lognormal(0, 1); ## discrimination (first period, if evolving)
    sd_gamma_geo ~ cauchy(0, 2.5); ## sd of geographic effects
    sd_gamma_demo ~ cauchy(0, 2.5); ## sd of demographic effects
    sd_innov_delta ~ cauchy(0, scale_sd_innov_delta); ## innovation sd of nu_geo, delta_gamma
    sd_innov_gamma ~ cauchy(0, scale_sd_innov_gamma); ## innovation sd. of gamma and xi
    sd_innov_item ~ cauchy(0, scale_sd_innov_item); ## innovation sd. of disc/diff
    sd_innov_logsd ~ cauchy(0, scale_sd_innov_logsd); ## innovation sd of theta_sd
    for (t in 1:T) { ## loop over years
      gamma_raw[t] ~ normal(0, 1);
      theta_bar_raw[t] ~ normal(0, 1); ## Matt trick done above #!#
      #!# theta_bar[t] ~ normal(mu_theta_bar[t], sd_theta_bar[t]); ## group means
      if (t == 1) {
        ## Priors for first period
        sd_theta_bar[t] ~ cauchy(0, 2.5);
        sd_theta[t] ~ cauchy(0, 2.5);
        nu_geo[t] ~ normal(0, 10);
        nu_geo_prior ~ normal(0, 10);
        delta_gamma[t] ~ normal(0.5, 0.5); ## 68% of prior mass btwn 0 and 1
        delta_tbar[t] ~ normal(delta_tbar_prior_mean, delta_tbar_prior_sd);
        xi[t] ~ normal(0, 10); ## intercept
      }
      if (t > 1) {
        ## TRANSITION MODEL
        ## Item parameters (if not constant)
        if (constant_diff == 0) {
          diff_raw[t] ~ normal(diff_raw[t - 1], sd_innov_item);
        }
        if (constant_disc == 0) {
          disc_raw[t] ~ lognormal(disc_raw[t - 1], sd_innov_item);
        }
        ## predictors in geographic models (random walk)
        delta_gamma[t] ~ normal(delta_gamma[t - 1], sd_innov_delta);
        nu_geo[t] ~ normal(nu_geo[t - 1], sd_innov_delta);
        delta_tbar[t] ~ normal(delta_tbar[t - 1], sd_innov_delta);
        sd_theta_bar[t] ~ lognormal(log(sd_theta_bar[t - 1]), sd_innov_logsd);
        sd_theta[t] ~ lognormal(log(sd_theta[t - 1]), sd_innov_logsd);
        if (separate_t == 0 && t > 2) {
          xi[t] ~ normal(xi[t - 1], sd_innov_gamma);
        }
        if (separate_t == 1 || t == 2) { ## Estimate model anew each period
          xi[t] ~ normal(0, 10);
        }
      }
      for (q in 1:Q) { ## loop over questions
        if (l2_only[t, q] == 1) {
          ## Second-level mean
          SSl2[t, q] ~ binomial(NNl2[t, q], prob_l2[t, q]);
        }
        for (g in 1:G) { ## loop over groups
          if (MMM[t, q, g] == 0) { ## Use only if not missing
            pos <- pos + 1;
            prob_vec[pos] <- prob[t, q, g];
          }
        } ## end group loop
      } ## end question loop
    } ## end time loop
    ## Sampling model for group responses
    s_vec ~ binomial(n_vec, prob_vec);
  }
  generated quantities {
    vector<lower=0>[T] sd_total;
    for (t in 1:T) {
      sd_total[t] <- sqrt(variance(theta_bar[t]) + square(sd_theta[t]));
    }
  }
"
