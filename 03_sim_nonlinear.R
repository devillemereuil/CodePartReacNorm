conflicted::conflict_prefer("filter", "dplyr")
conflicted::conflict_prefer("unpack", "tidyr")
conflicted::conflicts_prefer(purrr::set_names)

here::i_am("Scripts R/03_sim_nonlinear.R")

# Libraries
library(tidyverse)
library(patchwork)
library(grid)
library(mvtnorm)
library(cubature)
library(nlme)
library(furrr)
library(here)

setwd(here())

options(mc.cores = parallel::detectCores() - 2)
plan(multicore)

# Some ggplot options
theme_set(theme_bw())
theme_update(text = element_text(family = "Linux Biolinum O", size = 20),
             plot.title = element_text(face = "bold", hjust = 0.5))

## ------------------ Setting parameters

# Setting up a seed for the simulation
set.seed(2021)

# Number of simulations
n_sim   <- 1000
# Number of individuals
n_ind   <- 100
# Number of replicates
n_rep   <- 10

# Environment used
env_template <- seq(-2, 2, length.out = n_rep)

# Sigmoidal reaction norm parameters
avg_sig <-
    c(a = 1,
      b = 3)
mat_sigm <- matrix(c(0.1, 0.02, 0.02, 0.05), ncol = 2)
sigmoid <- function(x, pars) {
    pars[ , "a"] / (1 + exp(-(pars[ , "b"] * x)))
}
sigmoid_cuba <- function(e, pars) {
    pars[1, ] / (1 + exp(-(pars[2, ] * e)))
}
d_sigmoid_cuba <- function(e, pars) {
    matrix(c(
        1 / (1 + exp(-(pars[2, ] * e))),
        pars[1, ] * exp(-(pars[2, ] * e)) * e / (1 + exp(-(pars[2, ] * e)))^2
    ),
           nrow = 2,
           byrow = TRUE)
}

# Gaussian-Gompertz reaction norm parameter
avg_gg <-
    c(cmax       = 1,
      xopt       = 0.9,
      rho        = 8,
      sigma_gaus = 0.4)
mat_gausgompz <- matrix(c(0.1,      0.01,   0,          0,
                          0.01,     0.05,   0,          0,
                          0,        0,      0,          0,
                          0,        0,      0,          0),
                        ncol = 4)
gausgompz <- function(x, pars) {
    pars[ , "cmax"] * exp(
        - exp(pars[ , "rho"] * (x - pars[ , "xopt"]) - 6) -       # Gompertz part
            pars[ , "sigma_gaus"] * (x - pars[ , "xopt"])^2       # Gaussian part
    )
}
gausgompz_cuba <- function(e, pars) {
    pars[1, ] * exp(
        - exp(pars[3, ] * (e - pars[2, ]) - 6) -        # Gompertz part
            pars[4, ] * (e - pars[2, ])^2               # Gaussian part
    )
}
d_gausgompz_cuba <- function(e, pars) {
    matrix(c(
            exp(-exp(-6 - pars[2, ] * pars[3, ] + e * pars[3, ]) - (pars[2, ] - e)^2 * pars[4, ]),
            pars[1, ] * exp(-exp(-6 - pars[2, ] * pars[3, ] + e * pars[3, ]) - (pars[2, ] - e)^2 * pars[4, ]) * (exp(-6 - pars[2, ] * pars[3, ] + e * pars[3, ]) * pars[3, ] + 2 * (-pars[2, ] + e) * pars[4, ])#,
#             pars[1, ] * exp(-6 - exp(-6 - pars[2, ] * pars[3, ] + e * pars[3, ]) + (-pars[2, ] + e) * pars[3, ] - (pars[2, ] - e)^2 * pars[4, ]) * (pars[2, ] - e),
#             -pars[1, ] * exp(-exp(-6 - pars[2, ] * pars[3, ] + e * pars[3, ]) - (pars[2, ] - e)^2 * pars[4, ]) * (pars[2, ] - e)^2
    ),
           nrow = 2,
           byrow = TRUE)
}
fixed_gg <- which(diag(mat_gausgompz) == 0)

# Residuals
sigma <- 0.01

## ------------------ Functions for the simulations

## Function to simulate the reaction norm, with random variations
# Args: - n_sim: Number of simulations
#       - n_gen: Number of genotypes
#       - n_ind: Number of individuals
#       - env: vector of environmental values
#       - avg_pars: average parameters to be used for reaction norm (must be named)
#       - mat_pars: VCV matrix of the parameters
#       - shape_fn: function giving the shape of the reaction norm
# Value: A tbl with reaction norm for each env value (N * n_env)
simulate_rn <- function(n_sim, n_ind, n_rep, avg_pars, mat_pars, shape_fn) {
    # Generating the parameters for each indivudal
    pars <- rmvnorm(n_sim * n_ind, avg_pars, mat_pars)
    pars <- pars[rep(1:nrow(pars), each = n_rep), ]

    # Generating the environment for each measure
#     env <- rnorm(n_sim * n_ind * n_rep, env_template, 0.1)
#     env <- runif(n_sim * n_ind * n_rep, -2, 2)
    env <- rep(env_template, n_sim * n_ind)

    # Generating the reaction norms, then the individual phenotypes and formatting
    out <-
        tibble(Simulation = as_factor(rep(1:n_sim, each = n_ind * n_rep)),
               Individual = as_factor(rep(rep(1:n_ind, each = n_rep), n_sim)),
               Env        = env,
               Phen       = shape_fn(env, pars) +
                            rnorm(n_sim * n_ind * n_rep, 0, sigma)) |>
        nest_by(Simulation, .key = "Data") |>
        ungroup()

    return(out)
}

## Function to run the non-linear model on the reaction norm data using nlme
# Args: - df: A dataset containing the reaction norm data
#       - shape: the function of the true reaction norm
# Value: The model fit for the data
fit_nlme <- function(df, shape) {
    if (shape == "Sigmoid") {
        nlme(Phen ~ I(a / (1 + exp(-(b * Env)))),
             fixed  = a + b ~ 1,
             random = a + b ~ 1|Individual,
             data   = df,
             start  = avg_sig)
    } else {
        nlme(Phen ~ I(cmax * exp(- exp(rho * (Env - xopt) - 6) -
                                 sigma_gaus * (Env - xopt)^2)),
             fixed  = cmax + xopt + rho + sigma_gaus ~ 1,
             random = cmax + xopt ~ 1|Individual,
             data   = df,
             start  = avg_gg[c("cmax", "xopt", "rho", "sigma_gaus")])
    }
}

## Compute the G matrix from variances and correlation
# Args: - vars: vector of variances (diagonal elements of G)
#       - corr: correlation coefficient (only one as G are all dim 2)
# Value: the G-matrix
vc_to_vcv <- function(vars, corr) {
   # Removing the residual variance and generating a diagonal matrix from variances
   G <- diag(vars[-3])

   # Adding the covariance
   G[1, 2] <- G[2, 1] <- corr * sqrt(G[1, 1] * G[2, 2])

   return(G)
}

## Compute the average conditional to the environment (E_g_e)
# Args: - shape: the function of the reaction norm fitted by the model
#       - theta: the average parameters estimated by the model
#       - G_theta: the genetic variance-covariance matrix estimated by the model
#       - e: the environmental value to condition to
#       - width: the width over which the integral must be computed (10 is a generally a good value)
# Value: The value for E_g_e (numeric)
compute_Eg_e <- function(shape, theta, G_theta, e, width = 10, fixed = NA) {
    # This function is modified from the QGglmm package

    # Handling when some terms are fixed
    if (!any(is.na(fixed))) {
        var       <- setdiff(1:length(theta), fixed)
        full_theta <- theta
        var_theta  <- theta[-fixed]
        if (nrow(G_theta) == length(theta)) {
            G_theta <- G_theta[-fixed, -fixed]
        }
    } else {
        full_theta <- theta
        var_theta  <- theta
    }

    # Setting the integral width according to vcov (lower mean-w, upper mean+w)
    w <- sqrt(diag(G_theta)) * width

    # Number of dimensions
    d <- length(w)

    # Computing the logdet of vcov
    logdet <- calc_logdet(G_theta)

    # Average
    avg <- cubature::hcubature(
        f  = function(x) {
            full_x      <- matrix(full_theta, nrow = length(full_theta), ncol = ncol(x))
            if (!any(is.na(fixed))) { full_x[var, ] <- x } else { full_x <- x }
            shape(e, full_x) * vec_mvnorm(x, var_theta, G_theta, logdet)
        },
        lowerLimit = var_theta - w,
        upperLimit = var_theta + w,
        fDim       = 1,
        tol        = 0.001,
        absError   = 0.0001,
        vectorInterface = TRUE
    )$integral

    return(avg)
}

## Compute the plastic variance (V_plas)
# Args: - shape: the name of the shape used to simulate/estimate the RN
#       - theta: the average parameters estimated by the model
#       - vars: the vars estimated from the model
#       - env: the environmental values over which the model has been estimated
#       - width: the width over which the integral must be computed (10 is a generally a good value)
#       - fixed: which part of the parameters are fixed (no genetic variation)
#       - average: should the average of variances be returned?
#                  If FALSE, return the variance for each environmental value
# Value: The value for V_plas (numeric)
compute_vplas <- function(shape, theta, G_theta, env, width = 10, fixed = NA) {
    # Defining the shape function
    if (shape == "Sigmoid") {
        shape_func <- sigmoid_cuba
    } else {
        shape_func <- gausgompz_cuba
    }

    map_dbl(env,
            \(e) compute_Eg_e(e = e, shape = shape_func, theta = theta, G_theta = G_theta,
                              width = width, fixed = fixed),
            .progress = TRUE) |>
        rep(n_ind * n_rep) |>
        var()
}

## Log-determinant of a VCV matrix
calc_logdet <- function(Sigma) {
    sum(log(eigen(Sigma, symmetric = TRUE, only.values = TRUE)$values))
}

## Vectorised multivariate Gaussian density function
# Shamelessly stolen from cubature vignette (credit to Balasubramanian Narasimhan)
# logdet = sum(log(eigen(sigma, symmetric = TRUE, only.values = TRUE)$values))
vec_mvnorm <- function(x, mean, Sigma, logdet = NULL) {
    # If logdet not provided, compute it
    if (is.null(logdet)) {
        logdet <- calc_logdet(Sigma)
    }
    # Compute Mahalanobis distance (corresponds to the exp. part of the density)
    distval <- stats::mahalanobis(t(x), center = mean, cov = Sigma)
    # Compute the vectorised MVN density
    out <- exp(matrix(-(nrow(x) * log(2 * pi) + logdet + distval)/2, ncol = ncol(x)))
    return(out)
}


## Compute the genetic variance conditionnally to the environment (V_g_e)
# Args: - shape: the function of the reaction norm fitted by the model
#       - theta: the average parameters estimated by the model
#       - G_theta: the genetic variance-covariance matrix estimated by the model
#       - e: the environmental value to condition to
#       - width: the width over which the integral must be computed (10 is a generally a good value)
# Value: The value for V_g_e (numeric)
compute_vg_e <- function(shape, theta, G_theta, e, width = 10, fixed = NA) {
    # This function is modified from the QGglmm package

    # Handling when some terms are fixed
    if (!any(is.na(fixed))) {
        var       <- setdiff(1:length(theta), fixed)
        full_theta <- theta
        var_theta  <- theta[-fixed]
        if (nrow(G_theta) == length(theta)) {
            G_theta <- G_theta[-fixed, -fixed]
        }
    } else {
        full_theta <- theta
        var_theta  <- theta
    }

    # Setting the integral width according to vcov (lower mean-w, upper mean+w)
    w <- sqrt(diag(G_theta)) * width

    # Number of dimensions
    d <- length(w)

    # Computing the logdet of vcov
    logdet <- calc_logdet(G_theta)

    # Average
    avg <- cubature::hcubature(
        f  = function(x) {
            full_x      <- matrix(full_theta, nrow = length(full_theta), ncol = ncol(x))
            if (!any(is.na(fixed))) { full_x[var, ] <- x } else { full_x <- x }
            shape(e, full_x) * vec_mvnorm(x, var_theta, G_theta, logdet)
        },
        lowerLimit = var_theta - w,
        upperLimit = var_theta + w,
        fDim       = 1,
        tol        = 0.001,
        absError   = 0.0001,
        vectorInterface = TRUE
    )$integral

    # Computing the integral
    cubature::hcubature(
        f  = function(x) {
            full_x      <- matrix(full_theta, nrow = length(full_theta), ncol = ncol(x))
            if (!any(is.na(fixed))) { full_x[var, ] <- x } else { full_x <- x }
            (shape(e, full_x) - avg)^2 * vec_mvnorm(x, var_theta, G_theta, logdet)
        },
        lowerLimit = var_theta - w,
        upperLimit = var_theta + w,
        fDim       = 1,
        tol        = 0.001,
        absError   = 0.0001,
        vectorInterface = TRUE
    )$integral
}

## Compute the genetic variance (V_gen)
# Args: - shape: the name of the shape used to simulate/estimate the RN
#       - theta: the average parameters estimated by the model
#       - vars: the vars estimated from the model
#       - env: the environmental values over which the model has been estimated
#       - width: the width over which the integral must be computed (10 is a generally a good value)
#       - fixed: which part of the parameters are fixed (no genetic variation)
#       - average: should the average of variances be returned?
#                  If FALSE, return the variance for each environmental value
# Value: The value for V_gen (numeric)
compute_vgen <- function(shape, theta, G_theta, env, width = 10, fixed = NA, average = TRUE) {
    # Defining the shape function
    if (shape == "Sigmoid") {
        shape_func <- sigmoid_cuba
    } else {
        shape_func <- gausgompz_cuba
    }

    out <-
        map_dbl(env,
        \(e) compute_vg_e(e = e, shape = shape_func, theta = theta, G_theta = G_theta,
                            width = width, fixed = fixed),
        .progress = TRUE)

    if (average) { out <- mean(out) }

    return(out)
}

## Compute the additive genetic variance conditionnally to the environment (V_A_e)
# Args: - shape: the function of the reaction norm fitted by the model
#       - theta: the average parameters estimated by the model
#       - G_theta: the genetic variance-covariance matrix estimated by the model
#       - e: the environmental value to condition to
#       - width: the width over which the integral must be computed (10 is a generally a good value)
# Value: The value for V_A_e (numeric)
compute_va_e <- function(shape, theta, G_theta, e, width = 10, fixed = NA) {
    # This function is modified from the QGglmm package

    # Handling when some terms are fixed
    if (!any(is.na(fixed))) {
        var       <- setdiff(1:length(theta), fixed)
        full_theta <- theta
        var_theta  <- theta[-fixed]
        if (nrow(G_theta) == length(theta)) {
            G_theta <- G_theta[-fixed, -fixed]
        }
    } else {
        full_theta <- theta
        var_theta  <- theta
    }

    # Setting the integral width according to vcov (lower mean-w, upper mean+w)
    w <- sqrt(diag(G_theta)) * width

    # Number of dimensions
    d <- length(w)

    # Computing the logdet of vcov
    logdet <- calc_logdet(G_theta)

#     # Average
#     k <- shape(e, matrix(full_theta, nrow = length(full_theta)))

    # Computing the integral for Psi
    Psi <- cubature::hcubature(
        f  = function(x) {
            full_x      <- matrix(full_theta, nrow = length(full_theta), ncol = ncol(x))
            if (!any(is.na(fixed))) { full_x[var, ] <- x } else { full_x <- x }
            shape(e, full_x) * matrix(rep(vec_mvnorm(x, var_theta, G_theta, logdet), d),
                                      nrow = d,
                                      byrow = TRUE)
        },
        lowerLimit = var_theta - w,
        upperLimit = var_theta + w,
        fDim       = d,
        tol        = 0.001,
        absError   = 0.0001,
        vectorInterface = TRUE
    )$integral

    # Now, computing V_A_e and the gamma-decomposition
    out <-
        tibble(V_A = as.numeric(t(Psi) %*% G_theta %*% Psi),
               Gamma_1  = as.numeric(Psi[1]^2 * G_theta[1, 1]),
               Gamma_2  = as.numeric(Psi[2]^2 * G_theta[2, 2]),
               Gamma_12 = as.numeric(2 * Psi[1] * Psi[2] * G_theta[1, 2]))

    return(out)
}

## Compute the additive genetic variance (V_A)
# Args: - shape: the name of the shape used to simulate/estimate the RN
#       - theta: the average parameters estimated by the model
#       - vars: the vars estimated from the model
#       - env: the environmental values over which the model has been estimated
#       - width: the width over which the integral must be computed (10 is a generally a good value)
#       - fixed: which part of the parameters are fixed (no genetic variation)
#       - average: should the average of variances be returned?
#                  If FALSE, return the variance for each environmental value
# Value: The value for V_A (numeric)
compute_va <- function(shape, theta, G_theta, env, width = 10, fixed = NA, average = TRUE) {
    # Defining the shape function
    if (shape == "Sigmoid") {
        shape_func <- d_sigmoid_cuba
    } else {
        shape_func <- d_gausgompz_cuba
    }

    out <-
        map_dfr(env,
            \(e) compute_va_e(e = e, shape = shape_func, theta = theta, G_theta = G_theta,
                              width = width, fixed = fixed),
            .progress = TRUE)

    if (average) {
        out <-
            out |>
            summarise(V_A = mean(V_A),
                      across(starts_with("Gamma"), mean)) |>
            mutate(across(starts_with("Gamma"), \(v) { v / V_A }))
    } else {
        out <-
            out |>
            mutate(across(starts_with("Gamma"), \(v) { v / V_A }))
    }

    return(out)
}

## ------------------ Simulating the data

# Generating the datasets
tbl_sim <-
    bind_rows(
        bind_cols(Shape = "Sigmoid",
                  simulate_rn(n_sim    = n_sim,
                              n_ind    = n_ind,
                              n_rep    = n_rep,
                              avg_pars = avg_sig,
                              mat_pars = mat_sigm,
                              shape_fn = sigmoid)),
        bind_cols(Shape = "GausGompz",
                  simulate_rn(n_sim    = n_sim,
                              n_ind    = n_ind,
                              n_rep    = n_rep,
                              avg_pars = avg_gg,
                              mat_pars = mat_gausgompz,
                              shape_fn = gausgompz)),
    )

## ------------------ Running the statistical models

## Running the models
tbl_sim[["Mod"]] <- future_map2(tbl_sim[["Data"]],
                                tbl_sim[["Shape"]],
                                \(df, sh) { fit_nlme(df, sh) },
                                .progress = TRUE)

## Getting the estimates
tbl_sim <-
    tbl_sim |>
    mutate(Param = list(fixed.effects(Mod[[1]])),
           Vars  = VarCorr(Mod[[1]])[ , "Variance"] |>
                   as.numeric() |>
                   set_names(case_when(
                       Shape[[1]] == "Sigmoid" ~ c("V_a", "V_b", "V_R"),
                       Shape[[1]] == "GausGompz" ~ c("V_cmax", "V_xopt", "V_R"),
                   )) |>
                   list(),
           Corr  = VarCorr(Mod[[1]])[2 , "Corr"] |>
                   as.numeric(),
           .by = c(Shape, Simulation))

# A data.frame with only the estimates
tbl_estim <-
    tbl_sim |>
    reframe(Coef = enframe(c(Param[[1]], Vars[[1]]),
                           name = "Param",
                           value = "Value"),
            .by = c(Shape, Simulation)) |>
    unpack(Coef) |>
    mutate(Value = as.numeric(Value))

## ------------------ Computing the variance components

# Computing the variance components
tbl_vars <-
    tbl_sim |>
    mutate(G        = future_map2(Vars, Corr, vc_to_vcv),
           Fixed    = if_else(Shape == "GausGompz", list(fixed_gg), list(NA)),
           V_plas   = future_pmap_dbl(list(shape     = Shape,
                                           theta      = Param,
                                           G_theta    = G,
                                           fixed     = Fixed),
                                      compute_vplas,
                                      .options = furrr_options(seed = TRUE),
                                      env = env_template),
           V_g_e    = future_pmap(list(shape     = Shape,
                                       theta      = Param,
                                       G_theta    = G,
                                       fixed     = Fixed),
                                  compute_vgen,
                                  .options = furrr_options(seed = TRUE),
                                  env = env_template,
                                  average = FALSE),
           V_gen    = map_dbl(V_g_e, mean),
           V_A_e    = future_pmap(list(shape     = Shape,
                                       theta      = Param,
                                       G_theta    = G,
                                       fixed     = Fixed),
                                  compute_va,
                                  .options = furrr_options(seed = TRUE),
                                  env = env_template,
                                  average = FALSE),
           V_A      = future_pmap(list(shape     = Shape,
                                       theta      = Param,
                                       G_theta    = G,
                                       fixed     = Fixed),
                                  compute_va,
                                  .options = furrr_options(seed = TRUE),
                                  env = env_template,
                                  average = TRUE),
           V_res    = map_dbl(Vars, 3),
           V_tot    = V_plas + V_gen + V_res,
           V_phen   = map_dbl(Data, \(df) { var(df[["Phen"]]) })) |>
    select(Shape, Simulation, G, starts_with("V_")) |>
    unnest(V_A)

# Computing the true expected values
tbl_true_vars <-
    tibble(Shape    = c("Sigmoid", "GausGompz"),
           V_plas   = c(compute_vplas("Sigmoid", avg_sig, mat_sigm, env_template),
                        compute_vplas("GausGompz", avg_gg, mat_gausgompz, env_template,
                                      fixed = fixed_gg)),
           V_g_e    = list(compute_vgen("Sigmoid", avg_sig, mat_sigm, env_template,
                                        average = FALSE),
                           compute_vgen("GausGompz", avg_gg, mat_gausgompz, env_template,
                                        fixed = fixed_gg, average = FALSE)),
           V_gen    = map_dbl(V_g_e, mean),
           V_A_e    = list(compute_va("Sigmoid", avg_sig, mat_sigm, env_template,
                                      average = FALSE)[["V_A"]],
                           compute_va("GausGompz", avg_gg, mat_gausgompz, env_template,
                                      fixed = fixed_gg, average = FALSE)[["V_A"]]),
           V_A      = list(compute_va("Sigmoid", avg_sig, mat_sigm, env_template,
                                      average = TRUE),
                           compute_va("GausGompz", avg_gg, mat_gausgompz, env_template,
                                      fixed = fixed_gg, average = TRUE)),
           V_res    = sigma^2,
           V_tot    = V_plas + V_gen + V_res) |>
    unnest(V_A)

# Correlation between V_tot and V_phen
tbl_vars |>
    summarise(Corr = cor(V_tot, V_phen),
              .by = Shape) |>
    deframe()

## ------------------ Graphical outputs

options(scipen = 999)

## Raw parameter estimations
# Checking bias/precision on raw parameters

# Getting the relative errors on parameters
tbl_estim_gg <-
    tbl_estim |>
    mutate(
        Rel_Bias = case_when(
            Param == "a"            ~ (Value - avg_sig["a"]),
            Param == "b"            ~ (Value - avg_sig["b"]) / avg_sig["b"],
            Param == "V_a"          ~ (Value - mat_sigm[1, 1]) / mat_sigm[1, 1],
            Param == "V_b"          ~ (Value - mat_sigm[2, 2]) / mat_sigm[2, 2],
            Param == "cmax"         ~ (Value - avg_gg["cmax"]) / avg_gg["cmax"],
            Param == "xopt"         ~ (Value - avg_gg["xopt"]) / avg_gg["xopt"],
            Param == "rho"          ~ (Value - avg_gg["rho"]) / avg_gg["rho"],
            Param == "sigma_gaus"   ~ (Value - avg_gg["sigma_gaus"]) / avg_gg["sigma_gaus"],
            Param == "V_cmax"       ~ (Value - mat_gausgompz[1, 1]) / mat_gausgompz[1, 1],
            Param == "V_xopt"       ~ (Value - mat_gausgompz[2, 2]) / mat_gausgompz[2, 2],
            Param == "V_sig"        ~ (Value - mat_gausgompz[4, 4]) / mat_gausgompz[4, 4],
            Param == "V_R"          ~ (Value - sigma^2) / (sigma^2),
        ),
        Param =
            factor(Param,
                   levels = c("a", "b", "xopt", "cmax", "rho", "sigma_gaus",
                              "V_a", "V_b", "V_cmax", "V_xopt", "V_sig", "V_R")) |>
                fct_recode("V[a]" = "V_a",
                           "V[b]" = "V_b",
                           "V[Res]" = "V_R",
                           "C[max]" = "cmax",
                           "x[opt]" = "xopt",
                           "sigma" = "sigma_gaus",
                           "V[C[max]]" = "V_cmax",
                           "V[x[opt]]" = "V_xopt",
                           "V[sigma]" = "V_sig"),
        Shape = recode(Shape, Sigmoid = "Sigmoid", GausGompz = "Performance Curve") |>
                as_factor()
    )

# Plot of the relative errors
p_params <-
    ggplot(tbl_estim_gg) +
    geom_violin(aes(x = Param, y = Rel_Bias), fill = "grey", scale = "width") +
    geom_hline(yintercept = 0) +
    stat_summary(aes(x = Param, y = Rel_Bias),
                 geom = "point",
                 colour = "red",
                 size = 2,
                 fun = "mean") +
    facet_grid(~ Shape, scale = "free_x", space = "free_x") +
    scale_x_discrete(labels = scales::parse_format()) +
    labs(x = "Variance component", y = "Relative Error") +
    theme(text         = element_text(family = "Linux Biolinum O"),
          axis.text.y  = element_text(size = 14),
          axis.text.x  = element_text(size = 14),
          axis.title   = element_text(size = 16),
          strip.text   = element_text(size = 20))

# Saving the plot
cairo_pdf("Figs/NonLin_Parameters.pdf", height = 5, width = 8)
p_params
dev.off()

## Variance components
# Formatting the data
tbl_var_plot <-
    tbl_vars |>
    select(-G, -V_phen, -ends_with("_e"), -starts_with("Gamma")) |>
    pivot_longer(contains("V"),
                 names_to   = "Variance",
                 values_to  = "Value") |>
    left_join(tbl_true_vars |>
              select(-ends_with("_e")) |>
              pivot_longer(contains("V"),
                           names_to  = "Variance",
                           values_to = "Truth")) |>
    mutate(Rel_Bias = (Value - Truth) / Truth,
           Shape    = recode(Shape,
                             GausGompz = "Performance Curve") |>
                      as_factor(),
           Variance = recode(Variance,
                             V_plas = "V[Plas]",
                             V_gen  = "V[Gen]",
                             V_A    = "V[A]",
                             V_res  = "V[Res]",
                             V_tot  = "V[Tot]") |>
                      as_factor())

# Now making the plot
p_var <-
    ggplot(tbl_var_plot) +
    geom_violin(aes(x = Variance, y = Rel_Bias), fill = "grey", scale = "width") +
    geom_hline(yintercept = 0) +
    stat_summary(aes(x = Variance, y = Rel_Bias),
                 geom = "point",
                 colour = "red",
                 size = 2,
                 fun = "mean") +
    facet_grid(~ Shape, scale = "free_x", space = "free_x") +
    scale_x_discrete(labels = scales::parse_format()) +
    labs(x = "Variance component", y = "Relative Error") +
    theme(text         = element_text(family = "Linux Biolinum O"),
          axis.text.y  = element_text(size = 14),
          axis.text.x  = element_text(size = 14),
          axis.title   = element_text(size = 16),
          strip.text   = element_text(size = 20))

# Saving the plot
cairo_pdf("Figs/NonLin_Variances.pdf", height = 5, width = 8)
p_var
dev.off()

# Statistical test for bias
tbl_var_plot |>
    summarise(P_val = wilcox.test(Rel_Bias)[["p.value"]],
              .by = c(Shape, Variance)) |>
    pivot_wider(values_from = "P_val", names_from = "Variance")
tbl_var_plot |>
    summarise(Mean_Rel_Bias = mean(Rel_Bias), .by = c(Shape, Variance))

## Examples of reaction norms
## Sigmoid
tbl_annot <-
    tbl_true_vars |>
    select(-ends_with("_e"), -starts_with("Gamma"), -V_res) |>
    pivot_longer(contains("V"),
                 names_to = "Variance",
                 values_to = "Value") |>
    mutate(Variance = recode(Variance,
                             V_plas = "V[Plas]",
                             V_gen  = "V[Gen]",
                             V_A    = "V[A]",
                             V_tot  = "V[Tot]"),
           Value = signif(Value, digits = 3) |>
                   map_chr(format, scientific = FALSE),
           Annot = str_c(Variance, " == ", Value),
           X     = 3,
           Y     = seq(0.17, 0.13, length.out = n()),
           .by = Shape)

tbl_annot_gamma <-
    tbl_true_vars |>
    select(Shape, starts_with("Gamma")) |>
    pivot_longer(contains("Gamma"),
                 names_to = "Gamma",
                 values_to = "Value") |>
    mutate(X     = 3,
           Y     = seq(0.115, 0.09, length.out = n()),
           Gamma = case_when(Shape == "Sigmoid"   & Gamma == "Gamma_1"  ~ "γ[L]",
                             Shape == "Sigmoid"   & Gamma == "Gamma_2"  ~ "γ[r]",
                             Shape == "Sigmoid"   & Gamma == "Gamma_12" ~ "γ[Lr]",
                             Shape == "GausGompz" & Gamma == "Gamma_1"  ~ "γ[C]",
                             Shape == "GausGompz" & Gamma == "Gamma_2"  ~ "γ[ε[0]]",
                             Shape == "GausGompz" & Gamma == "Gamma_12" ~ "γ[Cε[0]]"),
           Value = signif(Value, digits = 2) |>
                   map_chr(format, scientific = FALSE),
           Annot = str_c(Gamma, " == ", Value),
           .by = Shape)

# Function of average parameters
fe_sig <- fixef(tbl_sim[["Mod"]][[1]])
func_sig_1 <- function(x) {
    fe_sig["a"] / (1 + exp(- fe_sig["b"] * x))
}

# Average phenotype for each environment
avg_phen_sig <-
    map_dbl(env_template,
             compute_Eg_e,
             shape   = sigmoid_cuba,
             theta    = fe_sig,
             G_theta  = tbl_vars[["G"]][[1]],
             fixed   = NA)

# Now generating the plot
p_sig <-
    ggplot() +
    geom_line(data = tbl_sim[["Data"]][[1]],
              aes(x = Env, y = Phen, group = Individual),
              colour = "grey",
              alpha = 0.5,
              linewidth = 1) +
    geom_function(fun = func_sig_1,
                  linewidth = 1.2,
                  colour = "black") +
    geom_smooth(data = NULL,
                aes(x = env_template, y = avg_phen_sig),
                method = "gam",
                linewidth = 1.2,
                colour = "#ff557f") +
    annotate("rect", xmin = -2, xmax = -0.3, ymin = 1.15, ymax = 1.8,
         fill = "white", colour = "#FFCC00", linewidth = 2) +
    annotate("text", x = -1.6, y = 1.7, label = "italic(bar(L)) == 1",
             size = 6,
             family = "Linux Libertine O",
             parse = TRUE) +
    annotate("text", x = -0.9, y = 1.7, label = "italic(bar(r)) == 3",
             size = 6,
             family = "Linux Libertine O",
             parse = TRUE) +
    annotate("text", x = -1.15, y = 1.4, label = "G[paste(L,',',r)] == bgroup('(', atop('0.1   0.02', '0.02 0.05'),')')",
             size = 6,
             family = "Linux Libertine O",
             parse = TRUE) +
    ylab("Phenotype") + xlab("Environment") +
    theme(legend.position = "none")

# Also generating the equation
p_eq_sig <-
    ggplot() +
    geom_text(aes(x = 0, y = 0.5, label = "Sigmoid"),
              family = "Linux Biolinum O",
              fontface = "bold",
              size = 10) +
    geom_text(aes(x = 0, y = 0, label = "italic(f(ɛ) == frac(L, 1 + e^{-r * ɛ}))"),
              size = 10,
              family = "Linux Libertine O",
              parse = TRUE) +
    ylim(c(-0.2, 0.6)) +
    theme_void() +
    theme(plot.margin = margin(0,0,0,0, "pt"))

## Performance curve

# Function of the average parameter
fe_gg <- fixef(tbl_sim[["Mod"]][[n_sim + 1]])
func_gg_1 <- function(x) {
    fe_gg["cmax"] * exp(
        - exp(fe_gg["rho"] * (x - fe_gg["xopt"]) - 6) -       # Gompertz part
            fe_gg["sigma_gaus"] * (x - fe_gg["xopt"])^2       # Gaussian part
    )
}

# Average phenotype in each environment
avg_phen_gg <-
    map_dbl(env_template,
            compute_Eg_e,
            shape   = gausgompz_cuba,
            theta    = fe_gg,
            G_theta  = tbl_vars[["G"]][[n_sim + 1]],
            fixed   = fixed_gg)

# Now generating the plot
p_gg <-
    ggplot() +
    geom_line(data = tbl_sim[["Data"]][[n_sim + 1]],
              aes(x = Env, y = Phen, group = Individual),
              colour = "grey",
              alpha = 0.5,
              linewidth = 1) +
    geom_function(fun = func_gg_1,
                  linewidth = 1.2,
                  colour = "black") +
    geom_smooth(data = NULL,
                aes(x = env_template, y = avg_phen_gg),
                method = "gam",
                linewidth = 1.2,
                colour = "#ff557f") +
    annotate("rect", xmin = -2, xmax = -0.3, ymin = 1.08, ymax = 1.7,
             fill = "white", colour = "#FFCC00", linewidth = 2) +
    annotate("text", x = -1.15, y = 1.6, label = "paste(italic(bar(C)) == 1, '    ', italic(bar(ε[0])) == 0.9)",
             size = 6,
             family = "Linux Libertine O",
             parse = TRUE) +
    annotate("text", x = -1.15, y = 1.45, label = "paste(italic(bar(ρ)) == 8, '    ', italic(bar(σ)) == 0.4)",
             size = 6,
             family = "Linux Libertine O",
             parse = TRUE) +
    annotate("text", x = -1.15, y = 1.25, label = "G[paste(C,',',ε[0])] == bgroup('(', atop('0.1   0.01', '0.01 0.05'),')')",
             size = 6,
             family = "Linux Libertine O",
             parse = TRUE) +
    ylab("Phenotype") + xlab("Environment") +
    theme(legend.position = "none")

# Also generating the equation
p_eq_gg <-
    ggplot() +
    geom_text(aes(x = 0, y = 0.5, label = "Performance Curve"),
              family = "Linux Biolinum O",
              fontface = "bold",
              size = 10) +
    geom_text(aes(x = 0, y = 0, label = "italic(f(ɛ) == C * e^{-e^{ρ * (ɛ - ɛ[0]) - 6} - σ * (ɛ - ɛ[0])^2})"),
              family = "Linux Libertine O",
              size = 10,
              parse = TRUE) +
    ylim(c(-0.2, 0.6)) +
    theme_void() +
    theme(plot.margin = margin(0,0,0,0, "pt"))

## Genetic variances across environments
# Getting variances across environments
tbl_vars_e <-
    tbl_vars |>
    select(Shape, Simulation, V_g_e, V_A_e) |>
    unnest(starts_with("V")) |>
    select(-starts_with("Gamma")) |>
    rename(V_A_e = V_A) |>
    mutate(Env = rep(env_template, 2 * n_sim)) |>
    pivot_longer(starts_with("V"),
                 names_to = "Variance",
                 values_to = "Value") |>
    left_join(tbl_true_vars |>
              select(Shape, ends_with("_e")) |>
              unnest(starts_with("V")) |>
              mutate(Env = rep(env_template, 2)) |>
              pivot_longer(starts_with("V"),
                           names_to = "Variance",
                           values_to = "Truth")) |>
    mutate(Variance = factor(Variance, levels = c("V_g_e", "V_A_e")))

# Computing the theoretical genetic variances for more values
tbl_true_vgen <-
    tibble(Env  = seq(-2, 2, length.out = 200),
           V_g_e_sig = map_dbl(Env,
                               compute_vg_e,
                               shape       = sigmoid_cuba,
                               theta        = avg_sig,
                               G_theta      = mat_sigm,
                               .progress   = TRUE),
           V_g_e_gg = map_dbl(Env,
                              compute_vg_e,
                              shape       = gausgompz_cuba,
                              theta        = avg_gg,
                              G_theta      = mat_gausgompz,
                              fixed       = fixed_gg,
                              .progress   = TRUE),
           V_A_e_sig = map_dfr(Env,
                               compute_va_e,
                               shape       = d_sigmoid_cuba,
                               theta        = avg_sig,
                               G_theta      = mat_sigm,
                               .progress   = TRUE)[["V_A"]],
           V_A_e_gg = map_dfr(Env,
                              compute_va_e,
                              shape       = d_gausgompz_cuba,
                              theta        = avg_gg,
                              G_theta      = mat_gausgompz,
                              fixed       = fixed_gg,
                              .progress   = TRUE)[["V_A"]]) |>
    pivot_longer(starts_with("V"),
                 names_pattern = "(V_[gA]_e)_(.*)",
                 names_to = c("Variance", "Shape"),
                 values_to = c("Value")) |>
    mutate(Shape = recode(Shape, sig = "Sigmoid", gg = "GausGompz"),
           Env_numfac = 5.5 + Env * (9/4),
           Variance = factor(Variance, levels = c("V_g_e", "V_A_e")))

# Now generating the plot for Sigmoid
p_sig_var_e <-
    ggplot(tbl_vars_e |> filter(Shape == "Sigmoid")) +
    geom_violin(aes(x = as_factor(signif(Env, digits =2)), y = Value, fill = Variance),
                scale = "width") +
    geom_line(data = tbl_true_vgen |> filter(Shape == "Sigmoid"),
              mapping= aes(x = Env_numfac,
                           y = Value,
                           colour = Variance,),
              position = position_dodge(width = 1),
              linewidth = 1.2) +
    stat_summary(mapping= aes(x = as_factor(signif(Env, digits =2)), y = Value, group = Variance),
                 position = position_dodge(width = 1),
                 fun    = mean,
                 shape  = 16,
                 colour = "red") +
    annotate("rect", xmin = 1, xmax = 5, ymin = 0.08, ymax = 0.18,
             fill = "white", colour = "#FFCC00", linewidth = 2) +
    geom_text(data = tbl_annot |> filter(Shape == "Sigmoid"),
              mapping = aes(x = X, y = Y, label = Annot, colour = Variance),
              parse = TRUE,
              family = "Linux Biolinum O",
              size  = 6) +
    geom_text(data = tbl_annot_gamma |> filter(Shape == "Sigmoid"),
              mapping = aes(x = X, y = Y, label = Annot),
              parse = TRUE,
              family = "Linux Biolinum O",
              size  = 6) +
    ylim(c(0, 0.19)) +
    scale_fill_manual(values = c("#55aaff", "#aaff7f"), guide = "none") +
    scale_colour_manual(values = c("#00DD00", "#0000AA", "#00DD00", "#0000AA", "red", "black"),
                        guide = "none") +
    xlab("Environment") + ylab("(Additive) Genetic Variance")

# Now generating the plot for the Performance Curve
p_gg_var_e <-
    ggplot(tbl_vars_e |> filter(Shape == "GausGompz")) +
    geom_violin(aes(x = as_factor(signif(Env, digits =2)), y = Value, fill = Variance),
                scale = "width") +
    geom_line(data = tbl_true_vgen |> filter(Shape == "GausGompz"),
              mapping= aes(x = Env_numfac,
                           y = Value,
                           colour = Variance,),
              position = position_dodge(width = 1),
              linewidth = 1.2) +
    stat_summary(mapping= aes(x = as_factor(signif(Env, digits =2)), y = Value, group = Variance),
                 position = position_dodge(width = 1),
                 fun    = mean,
                 shape  = 16,
                 colour = "red") +
    annotate("rect", xmin = 1, xmax = 5, ymin = 0.08, ymax = 0.18,
             fill = "white", colour = "#FFCC00", linewidth = 2) +
    geom_text(data = tbl_annot |> filter(Shape == "GausGompz"),
              mapping = aes(x = X, y = Y, label = Annot, colour = Variance),
              parse = TRUE,
              family = "Linux Biolinum O",
              size  = 6) +
    geom_text(data = tbl_annot_gamma |> filter(Shape == "GausGompz"),
              mapping = aes(x = X, y = Y, label = Annot),
              parse = TRUE,
              family = "Linux Biolinum O",
              size  = 6) +
    ylim(c(0, 0.19)) +
    scale_fill_manual(values = c("#55aaff", "#aaff7f"), guide = "none") +
    scale_colour_manual(values = c("#00DD00", "#0000AA", "#00DD00", "#0000AA", "red", "black"),
                        guide = "none") +
    xlab("Environment") + ylab("(Additive) Genetic Variance")


## Full graph with reaction norms and variance estimations
# Relative errors of the variance partitions for Sigmoid
p_var_sig <-
    ggplot(tbl_var_plot |> filter(Shape == "Sigmoid", Variance != "V[Res]")) +
    geom_violin(aes(x = Variance, y = Rel_Bias, fill = Variance), scale = "width") +
    geom_hline(yintercept = 0) +
    stat_summary(aes(x = Variance, y = Rel_Bias),
                 geom = "point",
                 colour = "red",
                 size = 2,
                 fun = "mean") +
    scale_x_discrete(labels = scales::parse_format()) +
    scale_fill_manual(values = c("#ff557f", "#55aaff", "#aaff7f", "grey"),
                      guide  = "none") +
    labs(x = "Variance component", y = "Relative Error") +
    theme(text         = element_text(family = "Linux Biolinum O"),
          axis.text.y  = element_text(size = 14),
          axis.text.x  = element_text(size = 14),
          axis.title   = element_text(size = 16),
          strip.text   = element_text(size = 20))

# Relative errors of the variance partitions for the Performance Curve
p_var_gg <-
    ggplot(tbl_var_plot |> filter(Shape == "Performance Curve", Variance != "V[Res]")) +
    geom_violin(aes(x = Variance, y = Rel_Bias, fill = Variance), scale = "width") +
    geom_hline(yintercept = 0) +
    stat_summary(aes(x = Variance, y = Rel_Bias),
                 geom = "point",
                 colour = "red",
                 size = 2,
                 fun = "mean") +
    scale_x_discrete(labels = scales::parse_format()) +
    scale_fill_manual(values = c("#ff557f", "#55aaff", "#aaff7f", "grey"),
                      guide  = "none") +
    labs(x = "Variance component", y = "Relative Error") +
    theme(text         = element_text(family = "Linux Biolinum O"),
          axis.text.y  = element_text(size = 14),
          axis.text.x  = element_text(size = 14),
          axis.title   = element_text(size = 16),
          strip.text   = element_text(size = 20))

# Now we can generate the full graph
cairo_pdf("Figs/NonLin_Full.pdf", width = 12, height = 16)
(p_eq_sig / p_sig / p_sig_var_e / p_var_sig +
 theme(plot.margin = unit(c(0,30,0,0), "pt"))) +
    plot_layout(heights = c(1, 2, 2, 2))  |
(p_eq_gg / p_gg / p_gg_var_e / p_var_gg) +
    plot_layout(heights = c(1, 2, 2, 2))
grid.draw(linesGrob(x = unit(c(0.5, 0.5), "npc"), y = unit(c(0.02, 0.98), "npc")))
dev.off()
