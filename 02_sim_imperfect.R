conflicted::conflict_prefer("filter", "dplyr")
conflicted::conflict_prefer("unpack", "tidyr")

here::i_am("Scripts R/02_sim_imperfect.R")

# Libraries
library(tidyverse)
library(patchwork)
library(grid)
library(lme4)
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
set.seed(2023)

# Number of simulations
n_sim   <- 1000
# Number of "replicates" (individuals within an environment)
n_rep   <- 10

# Environment
n_env <- c(4, 10)
env <- map(n_env, ~ seq(-2, 2, length.out = .))
env <- map(env, ~ set_names(., str_glue("Env{1:length(.)}")))
names(env) <- as.character(n_env)

# Sigmoidal reaction norm parameters
a <- 1
b <- 4
sigmoid <- function(x) {
    a / (1 + exp(- b * x))
}

# Gaussian-Gompertz reaction norm parameter
cmax        <- 1
xopt        <- 1.2
rho         <- 8
sigma_gaus  <- 0.4
gausgompz <- function(x) {
    cmax * exp(
        - exp(rho * (x - xopt) - 6) -       # Gompertz part
            sigma_gaus * (x - xopt)^2       # Gaussian part
    )
}

# All shapes of reaction norm
shapes <- list(Sigmoid = sigmoid, GausGompz = gausgompz)

# Residuals
sigma <- 0.1

## ------------------ Functions for the simulations

## Function to simulate the reaction norm, with random variations
# Args: - n_sim: Number of simulations
#       - n_rep: Number of replicates (individuals) for each environment
#       - env: vector of environmental values
#       - shape_fn: function giving the shape of the reaction norm
# Value: A tbl with reaction norm for each env value (N * n_env)
simulate_rn <- function(n_sim, n_rep, env, shape_fn) {
    n <- n_sim * n_rep
    mvtnorm::rmvnorm(n, shape_fn(env), diag(length(env)) * sigma^2) |>
        as_tibble(.name_repair = ~ str_glue("Env{1:length(env)}")) |>
        mutate(Simulation = rep(1:n_sim, each = n_rep), .before = 1) |>
        pivot_longer(contains("Env"),
                     values_to = "Phen",
                     names_to = "Env_ref") |>
        mutate(Env      = env[Env_ref],
               Env_ref  = as_factor(Env_ref)) |>
        nest_by(Simulation, .key = "Data") |>
        ungroup()
}

## Function to run the ANOVA model on the reaction norm data
# Args: - df: A dataset containing the reaction norm data
# Value: The model fit for the data
fit_anova <- function(df) {
    # Fitting the model
    lm(Phen ~ 0 + Env_ref, data = df)
}

## Function to run the linear model on the reaction norm data
# Args: - df: A dataset containing the reaction norm data
#       - shape: the function of the true reaction norm
# Value: The model fit for the data
fit_lm <- function(df, shape) {
    if (shape == "Sigmoid") {
        lm(Phen ~ Env + I(Env^2), data = df)
    } else {
        lm(Phen ~ Env + I(Env^2), data = df)
    }
}

## Function to estimate the plastic variance from ANOVA model
# Args: - model: A lm fit yielded by fit_anova
# Values: the plastic variance
compute_var_plas <- function(model) {
    var(summary(model)[["coefficients"]][ , "Estimate"]) -
        sum(diag(vcov(model) %*% cov(model.matrix(model))))
}

## Function to estimate the slope variance from lm model
# Args: - model: A lm fit yielded by fit_lm
#       - se: Should the estimate be corrected by the standard error?
# Values: A tbl containing the slope (V_b), curvature (V_c) and
#         total variance from the model (V_fix)
compute_var_fix <- function(model, se = TRUE) {
    # Getting the coefficients
    out   <- summary(model)[["coefficients"]]
    X     <- model.matrix(model) |>
             as_tibble() |>
             distinct()

    # Computing the total variance from fixed effects
    v_fix <- t(out[ , "Estimate"]) %*% cov(X) %*% out[ , "Estimate"]
    if (se) {
        v_fix <- v_fix - t(out[ , "Std. Error"]) %*% cov(X) %*% out[ , "Std. Error"]
    }

    # Computing the variance arising from slope
    v_b  <- out["Env", "Estimate"]^2 * var(X[["Env"]])
    if (se) {
        v_b <- v_b - out["Env", "Std. Error"]^2 * var(X[["Env"]])
    }

    # Computing the curvature variance as simply the difference
    # between the total and slope variances
    v_c <- v_fix - v_b

    tibble(V_b = v_b, V_c = as.vector(v_c), V_fix = as.vector(v_fix))
}

## ------------------ Now running the simulations

## Generating the reaction norms datasets to be analysed
tbl_sim <-
    # Generating a crossed design between all parameters
    crossing(Env    = env,
             N_Rep  = n_rep,
             Shape  = c("Sigmoid", "GausGompz")) |>
    # Simulating the reaction norms
    mutate(Data     = pmap(list(N_Rep, Env, Shape),
                           ~ simulate_rn(n_sim, ..1, ..2, shapes[[..3]]))) |>
    unnest(Data)

## Running the models
# Using the function fit_anova, see above
tbl_sim <-
    tbl_sim |>
    mutate(ANOVA = future_map(Data, fit_anova, .progress = TRUE),
           LM    = future_map2(Data, Shape, fit_lm, .progress = TRUE))

## Estimating the plastic variance
tbl_sim <-
    tbl_sim |>
    mutate(V_LM     = map(LM, compute_var_fix) |> list_rbind(),
           V_plas   = map_dbl(ANOVA, compute_var_plas),
           V_res    = map_dbl(ANOVA, \(mod) { summary(mod)$sigma^2 }),
           V_tot    = V_plas + V_res) |>
    unpack(V_LM)

## ------------------ Graphics of the reaction norms and the model

## Sigmoids

# Select a simulated reaction norm
rn_sig <-
    tbl_sim |>
    filter(map_int(Env, length) == 10,
           N_Rep == 10,
           Shape == "Sigmoid") |>
    slice(1) |>
    flatten()

# Predicted values from the model
predict_lm_sig <-
    tibble(Env      = seq(-2.2, 2.2, 0.1)) %>%
    mutate(Predict  = predict(rn_sig[["LM"]], newdata = .),
           Shape    = sigmoid(Env))

# Predicted values from the ANOVA
predict_anova_sig <-
    tibble(Env_ref  = names(env[[2]]),
           Env      = env[[2]]) %>%
    mutate(Predict  = predict(rn_sig[["ANOVA"]], newdata = .))

# Graphics of the reaction norm
p_rn_sig <-
    ggplot() +
    geom_point(data     = rn_sig[["Data"]],
               mapping  = aes(x = Env, y = Phen),
               colour   = "grey70") +
    geom_line(data      = predict_lm_sig,
              mapping   = aes(x = Env, y = Shape),
              colour    = "#CC0000",
              size      = 1) +
    geom_line(data      = predict_lm_sig,
              mapping   = aes(x = Env, y = Predict),
              colour    = "black",
              size      = 1) +
    geom_point(data     = predict_anova_sig,
               mapping  = aes(x = Env, y = Predict),
               colour   = "#ff55ff",
               shape    = "square",
               size     = 3) +
    labs(x = "Environment", y = "Phenotype")

## Gaussian-Gompertz

# Select a simulated reaction norm
rn_ggz <-
    tbl_sim |>
    filter(map_int(Env, length) == 10,
           N_Rep == 10,
           Shape == "GausGompz") |>
    slice(1) |>
    flatten()

# Predicted values from the model
predict_lm_ggz <-
    tibble(Env      = seq(-2.2, 2.2, 0.1)) %>%
    mutate(Predict  = predict(rn_ggz[["LM"]], newdata = .),
           Shape    = gausgompz(Env))

# Predicted values from the ANOVA
predict_anova_ggz <-
    tibble(Env_ref  = names(env[[2]]),
           Env      = env[[2]]) %>%
    mutate(Predict  = predict(rn_ggz[["ANOVA"]], newdata = .))

# Graphics of the reaction norm
p_rn_ggz <-
    ggplot() +
    geom_point(data     = rn_ggz[["Data"]],
               mapping  = aes(x = Env, y = Phen),
               colour   = "grey70") +
    geom_line(data      = predict_lm_ggz,
              mapping   = aes(x = Env, y = Shape),
              colour    = "#CC0000",
              size      = 1) +
    geom_line(data      = predict_lm_ggz,
              mapping   = aes(x = Env, y = Predict),
              colour    = "black",
              size      = 1) +
    geom_point(data     = predict_anova_ggz,
               mapping  = aes(x = Env, y = Predict),
               colour   = "#ff55ff",
               shape    = "square",
               size     = 3) +
    labs(x = "Environment", y = "Phenotype")

## ------------------ Graphics of the simulations

## Formatting a dataset for the graph on simulations
tbl_sim_graph <-
    tbl_sim |>
    select(!Simulation:LM) |>
    mutate(Phi_1  = V_b / V_plas,
           Phi_2  = V_c / V_plas,
           M2    = (V_b + V_c) / V_plas,
           P2    = V_plas / V_tot) |>
    mutate(N_Env = map_int(Env, length), Env = NULL, .before = 1) |>
    select(-starts_with("V")) |>
    pivot_longer(matches("^[PM]"), names_to = "Component", values_to = "Variance") |>
    mutate(Colour = case_when(Component == "M2"     ~ "black",
                              Component == "P2"     ~ "#AA0000",
                              Component == "Phi_1"  ~ "#AA0000",
                              Component == "Phi_2"  ~ "#AA0000"))

## Another dataset containing the true simulated values
tbl_truth <-
    tbl_sim |>
    select(Env:Shape) |>
    distinct() |>
    mutate(N_Env    = names(Env),
           V_Env    = map_dbl(Env, var),
           V_Env2   = map_dbl(Env, ~ var(. ^ 2)),
           # NOTE This correction factor is needed to "translate" between
              # short and long versions of the reaction norm
           Corr_fac  = (length(N_Env) - 1) * N_Rep / (N_Rep * length(N_Env) - 1),
           V_plas   = map2_dbl(Shape, N_Env,
                           ~ var(shapes[[..1]](env[[as.character(..2)]]))),
           LM       = map2(Env, Shape,
                           ~ fit_lm(tibble(Env = ..1, Phen = shapes[[..2]](..1)), ..2)),
           V_LM     = map(LM, compute_var_fix, se = FALSE) |>
                      list_rbind(),
           V_res    = sigma^2,
           V_tot    = V_plas + V_res,
            Corr_fac = NULL) |>
    select(-V_Env, -V_Env2, -LM) |>
    unpack(V_LM) |>
    mutate(Phi_1 = V_b / V_plas,
           Phi_2 = if_else((V_c / V_plas) < 1e-5, 0, V_c / V_plas),
           M2    = (V_b + V_c) / V_plas,
           P2    = V_plas / V_tot) |>
    select(-starts_with("V")) |>
    pivot_longer(matches("^[PM]"), names_to = "Component", values_to = "Variance")

## Final dataset containing the pi from each component
tbl_rs <-
    tbl_truth |>
    mutate(N_Env = as.numeric(N_Env)) |>
    filter(N_Env == 10) |>
    mutate(across(where(is.numeric), ~ format(., digits = 2))) |>
    mutate(X = -1,
           Y = case_when(Component == "M2"     ~ 1.15,
                         Component == "P2"     ~ 1.05,
                         Component == "Phi_1"  ~ 0.95,
                         Component == "Phi_2"  ~ 0.84),
           Colour    = case_when(Component == "M2"     ~ "black",
                                 Component == "P2"     ~ "#AA0000",
                                 Component == "Phi_1"  ~ "#AA0000",
                                 Component == "Phi_2"  ~ "#AA0000"),
           Component = recode(Component,
                              Phi_1   = "φ[1]",
                              Phi_2   = "φ[2]",
                              M2      = "M[Plas]^2",
                              P2      = "P[RN]^2") |>
                       factor(levels = c("M[Plas]^2", "P[RN]^2", "φ[1]", "φ[2]")),
           Label = str_glue("{Component} == {Variance}"))

## A function to format the labels properly
format_labels <-
    . %>%
    mutate(Component = recode(Component,
                              Phi_1   = "φ[1]",
                              Phi_2   = "φ[2]",
                              M2      = "M[Plas]^2",
                              P2      = "P[RN]^2") |>
                       factor(levels = c("M[Plas]^2", "P[RN]^2", "φ[1]", "φ[2]")),
           N_Env = str_c("N[env] == ", N_Env))
tbl_truth <- format_labels(tbl_truth)
tbl_sim_graph <- format_labels(tbl_sim_graph)


## Preparing the yellow rectangle containing the pi values on the top graphs
yellow_rect <- annotate("rect", xmin = -1.69, xmax = -0.33, ymin = 0.77, ymax = 1.23,
                        fill = "white", colour = "#FFCC00", linewidth = 2)

## Reaction norm for the sigmoid
p_rn_sig_ann <-
    p_rn_sig +
    yellow_rect +
    geom_text(data = tbl_rs |> filter(Shape == "Sigmoid"),
              mapping = aes(x = X, y = Y, label = Label, colour = I(Colour)),
              parse = TRUE,
              family = "Linux Biolinum O",
              size  = 6) +
    ylim(c(-0.25, 1.25))

## Reaction norm for the Gompertz-Gaussian
p_rn_ggz_ann <-
    p_rn_ggz +
    yellow_rect +
    geom_text(data = tbl_rs |> filter(Shape == "GausGompz"),
              mapping = aes(x = X, y = Y, label = Label, colour = I(Colour)),
              parse = TRUE,
              family = "Linux Biolinum O",
              size  = 6) +
    ylim(c(-0.25, 1.25))

## Result of the simulations for the sigmoid
p_var_sig <-
    ggplot(tbl_sim_graph |> filter(Shape == "Sigmoid")) +
    geom_violin(aes(x = Component, y = Variance, fill = I(Colour)),
                alpha = 0.5,
                scale = "width") +
    stat_summary(aes(x = Component, y = Variance),
                 fun = mean,
                 shape = 16,
                 colour = "grey") +
    geom_point(data = tbl_truth |> filter(Shape == "Sigmoid"),
               mapping = aes(x = Component, y = Variance),
               size = 3,
               shape = 4,
               colour = "black") +
    facet_grid(N_Env ~ ., labeller = "label_parsed") +
    scale_x_discrete(labels = scales::parse_format())

## Results of the simulations for the Gompertz-Gaussian
p_var_ggz <-
    ggplot(tbl_sim_graph |> filter(Shape == "GausGompz")) +
    geom_violin(aes(x = Component, y = Variance, fill = I(Colour)),
                alpha = 0.5,
                scale = "width") +
    stat_summary(aes(x = Component, y = Variance),
                 fun = mean,
                 shape = 16,
                 colour = "grey") +
    geom_point(data = tbl_truth |> filter(Shape == "GausGompz"),
               mapping = aes(x = Component, y = Variance),
               size = 3,
               shape = 4,
               colour = "black") +
    facet_grid(N_Env ~ ., labeller = "label_parsed") +
    scale_x_discrete(labels = scales::parse_format())

## Equations
p_eq_sig <-
    ggplot() +
    geom_text(aes(x = 0, y = 0.5, label = "Sigmoid"),
              family = "Linux Biolinum O",
              fontface = "bold",
              size = 10) +
    geom_text(aes(x = -0.2, y = 0.1, label = "italic(f(ɛ) == frac(L, 1 + e^{-r * ɛ}))"),
              size = 10,
              family = "Linux Libertine O",
              parse = TRUE) +
    geom_text(aes(x = 0.3, y = 0.2, label = "italic(L) == 1"),
              size = 8,
              family = "Linux Libertine O",
              parse = TRUE) +
    geom_text(aes(x = 0.3, y = 0, label = "italic(r) == 4"),
              size = 8,
              family = "Linux Libertine O",
              parse = TRUE) +
    ylim(c(-0.2, 0.6)) + xlim(c(-0.5, 0.5)) +
    theme_void() +
    theme(plot.margin = margin(0,0,0,0, "pt"))

p_eq_gg <-
    ggplot() +
    geom_text(aes(x = 0, y = 0.5, label = "Thermal Perf. Curve"),
              family = "Linux Biolinum O",
              fontface = "bold",
              size = 10) +
    geom_text(aes(x = 0, y = 0.2, label = "italic(f(ɛ) == C * e^{-e^{ρ * (ɛ - ɛ[0]) - 6} - σ * (ɛ - ɛ[0])^2})"),
              family = "Linux Libertine O",
              size = 10,
              parse = TRUE) +
    geom_text(aes(x = 0, y = -0.1, label = "paste(italic(C) == 1, ', ', italic(ɛ[0]) == 1.2, ', ', italic(ρ) == 8, ', ', italic(σ) == 0.4)"),
              size = 8,
              family = "Linux Libertine O",
              parse = TRUE) +
    ylim(c(-0.2, 0.6)) + xlim(c(-0.5, 0.5)) +
    ylim(c(-0.2, 0.6)) +
    theme_void() +
    theme(plot.margin = margin(0,0,0,0, "pt"))

## Now, generating the whole graph
cairo_pdf("Figs/Unknown_shape.pdf", width = 12, height = 14)
(p_eq_sig / p_rn_sig_ann / p_var_sig +
    theme(plot.margin = unit(c(0,30,0,0), "pt"))) +
    plot_layout(heights = c(1, 3, 3)) |
(p_eq_gg / p_rn_ggz_ann / p_var_ggz) +
    plot_layout(heights = c(1, 3, 3))
grid.draw(linesGrob(x = unit(c(0.5, 0.5), "npc"), y = unit(c(0.02, 0.98), "npc")))
dev.off()
