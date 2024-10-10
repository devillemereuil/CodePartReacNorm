conflicted::conflicts_prefer(dplyr::filter)
conflicted::conflicts_prefer(tidyr::unpack)

here::i_am("Scripts R/04_sim_perfect_continuous.R")

# Libraries
library(tidyverse)
library(furrr)
library(here)
library(mvtnorm)
library(lme4)
library(patchwork)

setwd(here())

options(mc.cores = min(parallel::detectCores() - 2, 10))
plan(multicore)

# Some ggplot options
theme_set(theme_bw())
theme_update(text = element_text(family = "Linux Biolinum O", size = 20),
             plot.title = element_text(face = "bold", hjust = 0.5))

## ------------------ Setting parameters

# Setting up a seed for the simulation
set.seed(2021)

## Number of simulations
n_sim <- 1000
n_chunk <- n_sim / 10

## Parameters of the reaction norm (RN: z = a + b * e + c * e²)
beta <- c(1.5, 0.5, -0.5)

## Variance-covariance matrix for the RN parameters
sd_beta <- c(0.3, 0.4, 0.2)
mat_Sigma <- diag(sd_beta^2)
mat_Sigma[1, 2] <- mat_Sigma[2, 1] <- -0.2 * sd_beta[1] * sd_beta[2]
mat_Sigma[1, 3] <- mat_Sigma[3, 1] <- -0.2 * sd_beta[1] * sd_beta[3]
mat_Sigma[2, 3] <- mat_Sigma[3, 2] <-  0.1 * sd_beta[2] * sd_beta[3]

## Residual variance
res_sd <- 0.5

# ## Environment moments from uniform(-2, 2)
# mean_env <- 0
# mean_env2 <- 4/3
# mean_env3 <- 0
# mean_env4 <- 16/5
# var_env <- 4/3
# var_env2 <- 16/5 - (4/3)^2
# M_env <- matrix(c(1, mean_env, mean_env2,
#                   mean_env, mean_env2, mean_env3,
#                   mean_env2, mean_env3, mean_env4),
#                 nrow = 3, ncol = 3)
# env_VCV <- matrix(c(0, 0, 0, 0, var_env, 0, 0, 0, var_env2), nrow = 3, ncol = 3)

## Environment moments from norm(0, 1)
mean_env <- 0
mean_env2 <- 1
mean_env3 <- 0
mean_env4 <- 3
var_env <- 1
var_env2 <- 2
M_env <- matrix(c(1, mean_env, mean_env2,
                  mean_env, mean_env2, mean_env3,
                  mean_env2, mean_env3, mean_env4),
                nrow = 3, ncol = 3)
env_VCV <- matrix(c(0, 0, 0, 0, var_env, 0, 0, 0, var_env2), nrow = 3, ncol = 3)
avg_grad <- c(1, mean_env, mean_env2)

## True variances
true_V_plas <- as.vector(t(beta) %*% env_VCV %*% beta)
true_V_Add  <- sum(M_env * mat_Sigma)
true_V_A    <- as.vector(t(avg_grad) %*% mat_Sigma %*% avg_grad)
true_V_AxE  <- true_V_Add - true_V_A

## Setting up the simulation parameters
tbl_p <-
    tribble(~ Scenario, ~ N_Env, ~ N_Gen,
            1, 10, 200,
            2, 10, 50,
            3, 4, 200,
            4, 4, 50)

## ------------------ Functions for the simulations

## Function to simulate the reaction norm, with random variations
# Args: - n: Number of replicates
#       - env: vector of environmental values
#       - beta: average parameters of the reaction norm
#       - mat_vcv: variance-covariance matrix around beta
# Value: A tbl with reaction norm for each env value (N * n_env)
sim_rn <- function(n, n_env, beta, mat_vcv) {
    # Generating the random variation around the parameters
    params <- rmvnorm(n, beta, mat_vcv)

    # Generating the environmental variable
    env <- rnorm(n * n_env) |>
           scale() |>
           as.vector() |>
           matrix(ncol = n_env, nrow = n)

    # Simulating the reaction norm
    rn <- params[ , 1] + params[ , 2] * env + params[ , 3] * env^2

    # Outputting everything
    tibble(Geno     = rep(1:n, each = n_env),
           Env      = t(env) |> as.vector(),
           Reac_Norm = t(rn) |> as.vector())
}

## Function to run the reaction norm model on the data
# Args: - df: A dataset containing the phenotypic data
# Value: The model fit for the data
fit_lme_rn <- function(df) {
    # Formatting the data
    df[["Env_Sq"]] <- (df[["Env"]] - mean(df[["Env"]]))^2

    # Fitting the model
    lmer(Phen ~ Env + Env_Sq + (1 + Env + Env_Sq | Genotype),
         data = df)
}

## Function to extract the parameters of interest from the RN lme4 fit
# Args: - model: A lme fit yielded by fit_lme_rn
# Values: List of relevant parameters
extract_params_rn <- function(model) {
    # Variance-covariance matrix
    mat_vcv <- VarCorr(model)[["Genotype"]]
    attr(mat_vcv, "stddev") <- attr(mat_vcv, "correlation") <- NULL

    # Average (fixed) effects
    out   <- summary(model)[["coefficients"]]
    coefs <- out[ , "Estimate"]
    se    <- out[ , "Std. Error"]

    # Residual variance
    res_var <- attr(VarCorr(model), "sc")^2

    # Returning the list of parameters
    tibble(a    = coefs[1],
           b    = coefs[2],
           c    = coefs[3],
           a_se = se[1],
           b_se = se[2],
           c_se = se[3],
           Va   = mat_vcv[1, 1],
           Vb   = mat_vcv[2, 2],
           Vc   = mat_vcv[3, 3],
           Cab  = mat_vcv[1, 2],
           Cac  = mat_vcv[1, 3],
           Cbc  = mat_vcv[2, 3],
           Vr   = res_var)
}

## Function to compute the expected genotypic variance
# Args: - vcv: the variance-covariance matrix between the RN parameters
#       - mat_env: the "design matrix" of the reaction norm
# Values: The expected genotypic variance
compute_var_add <- function(vcv, mat_env) {
    # Sum of the element-wise multiplication of products with mat_env
    # with the elements of the variance-covariance matrix vcv
    # Equivalent (but faster) than averaging over the t(E[i, ]) %*% vcv %*% E[i, ] products
    sum(diag((1/nrow(mat_env)) * (t(mat_env) %*% mat_env) %*% vcv))
}

## Function to compute the expected V_A
# Args: - vcv: the variance-covariance matrix between the RN parameters
#       - mat_env: the "design matrix" of the reaction norm
# Values: The expected genotypic variance
compute_var_A <- function(vcv, mat_env) {
    # Sum of the element-wise multiplication of products with mat_env
    # with the elements of the variance-covariance matrix vcv
    # Equivalent (but faster) than averaging over the t(E[i, ]) %*% vcv %*% E[i, ] products
    as.numeric(t(colMeans(mat_env)) %*% vcv %*% colMeans(mat_env))
}

## Function to compute the expected V_AxE
# Args: - vcv: the variance-covariance matrix between the RN parameters
#       - vcv_env: the variance-covariance matrix of the environment
# Values: The expected genotypic variance
compute_var_AxE <- function(vcv, vcv_env) {
    # Sum of the element-wise multiplication of products with mat_env
    # with the elements of the variance-covariance matrix vcv
    # Equivalent (but faster) than averaging over the t(E[i, ]) %*% vcv %*% E[i, ] products
    sum(diag(vcv_env * vcv))
}

## Set colours for graphics depending on the type of variance
# Args: - vec_pars: a vector of parameter names
# Value: A vector of HTML colour names corresponding to the names in the input (character)
set_colours <- function(vec_pars) {
    case_when(
        vec_pars == "V[Plas]"      ~ "#AA0000",
        vec_pars == "P[RN]^2"      ~ "#AA0000",
        str_detect(vec_pars, "π")  ~ "#AA0000",
        vec_pars == "V[Gen]"       ~ "#AA5500",
        vec_pars == "V[G,ε]"       ~ "#AA5500",
        vec_pars == "H[RN]^2"      ~ "#AA5500",
        vec_pars == "V[Add]"       ~ "#68349C",
        vec_pars == "V[A,ε]"       ~ "#68349C",
        vec_pars == "h[RN]^2"      ~ "#68349C",
        str_detect(vec_pars, "γ")  ~ "#68349C",
        vec_pars == "V[A]"         ~ "#228B22",
        vec_pars == "h^2"          ~ "#228B22",
        vec_pars == "V[AxE]"       ~ "#2D489A",
        vec_pars == "h[I]^2"       ~ "#2D489A",
        str_detect(vec_pars, "ι")  ~ "#2D489A"
    )
}

## ------------------ Visualising the reaction norm variability

## Simulating a few reaction norms
# Number of reaction norms to simulate
N_vis <- 30
tbl_simvis <-
    sim_rn(N_vis, 10, beta, mat_Sigma)

## Plotting the reaction norms
p_rn <-
    ggplot(tbl_simvis) +
    geom_line(aes(x = Env, y = Reac_Norm, group = Geno),
              alpha = 0.5,
              linewidth = 1)
cairo_pdf(here("Figs/Sim_cont.pdf"))
plot(p_rn)
dev.off()

## ------------------ Now running the actual simulations

## Generating the reaction norms datasets to be analysed
tbl_sim <-
    # Generating the reaction norms for the total number of simulations
    tbl_p |>
    mutate(Reac_Norm = map2(N_Env, N_Gen,
                            \(n_env, n_gen) {
                                sim_rn(n_sim * n_gen, n_env, beta, mat_Sigma)
                            })) |>
    unnest(Reac_Norm) |>
    mutate(Simulation = rep(1:n_sim, each = unique(N_Gen) * unique(N_Env)),
           Genotype   = Geno %% unique(N_Gen),
           Genotype   = if_else(Genotype == 0, unique(N_Gen), Genotype) |>
                        str_pad(2, pad = "0"),
           Geno       = NULL,
           .before    = 1,
           .by        = Scenario) |>
    mutate(Env = as.vector(scale(Env)),
           .by = c(Simulation, Scenario)) |>
    mutate(Phen = Reac_Norm + rnorm(n(), mean = 0, sd = res_sd)) |>
    nest_by(Scenario, Simulation, .key = "Data") |>
    ungroup() |>
    arrange(Simulation)

## Running the models
# Using the function fit_lme_rn, see above
gc()
tbl_sim[["Model"]] <-
    future_map(tbl_sim[["Data"]],
               fit_lme_rn,
               .progress = TRUE,
               .options = furrr_options(seed = TRUE,
                                        chunk_size = n_chunk))
gc()

## Saving the design matrix
mat_env <-
    tbl_sim |>
    summarise(M = list(model.matrix(Model[[1]])),
              .by = Scenario)


## Extracting the parameters of interests
# Again see above for function extract_params_rn
tbl_sim <-
    tbl_sim |>
    mutate(Params = map(Model, extract_params_rn) |> list_rbind()) |>
    unpack(Params)

## Remove the (heavy) Model columns
tbl_sim <-
    tbl_sim |>
    select(-Model) |>
    arrange(Scenario)

## Saving the results
saveRDS(tbl_sim, file = here("rn_sims_cont.rds"))

## Loading previously saved results
tbl_sim <- readRDS(here("rn_sims_cont.rds"))

## ------------------ Graphically checking for bias

## Creating a dataset to study bias
tbl_bias <-
    tbl_sim |>
    select(-Data) |>
    # Computing the bias for each parameters
    mutate(a         = a - beta[1],
           b         = b - beta[2],
           c         = c - beta[3],
           Va        = Va - mat_Sigma[1, 1],
           Vb        = Vb - mat_Sigma[2, 2],
           Vc        = Vc - mat_Sigma[3, 3],
           Cab       = Cab - mat_Sigma[1, 2],
           Cac       = Cac - mat_Sigma[1, 3],
           Cbc       = Cbc - mat_Sigma[2, 3],
           Vr        = Vr - res_sd^2,
           .by = Scenario) |>
    select(-Simulation) |>
    pivot_longer(-Scenario, names_to = "Parameter", values_to = "Bias") |>
    mutate(Parameter = as_factor(Parameter))

## Generating the plot
p_bias <-
    ggplot(tbl_bias) +
    geom_violin(aes(x = Parameter, y = Bias),
                scale = "width",
                fill = "grey") +
    geom_hline(yintercept = 0)  +
    stat_summary(aes(x = Parameter, y = Bias),
                 geom = "point",
                 colour = "red",
                 size = 2,
                 fun = "mean") +
    facet_wrap(~ Scenario, ncol = 1)

cairo_pdf(here("Figs/Sim_Bias_cont.pdf"), width = 7, height = 9)
plot(p_bias)
dev.off()

## ------------------ Phenotypic variance partitioning

## Compute the variance partitions (including "realised" phenotypic variance V_P)
# V_Plas (for Reaction Norm) = V(E(z|e)) and V_Add (for genotype) = E(V(z|e))
# V_P is the realised phenotypic variance and V_tot is the sum of variance components
# (note that V_P ~ V_tot, but strict equality isn't expected)
tbl_part <-
    tbl_sim |>
    mutate(X_Env     = list(cbind(1,
                                  Data[[1]][["Env"]],
                                  (Data[[1]][["Env"]] - mean(Data[[1]][["Env"]]))^2)),
           Env_VCV   = list(cov(X_Env[[1]])),
           V_Plas    = as.vector(t(c(a, b, c)) %*% Env_VCV[[1]] %*% c(a, b, c)) -
                       as.vector(t(c(a_se, b_se, c_se)) %*% Env_VCV[[1]] %*% c(a_se, b_se, c_se)),
           V_Add     = compute_var_add(matrix(c(Va, Cab, Cac, Cab, Vb, Cbc, Cac, Cbc, Vc),
                                               ncol = 3, nrow = 3),
                                        X_Env[[1]]),
           V_A       = compute_var_A(matrix(c(Va, Cab, Cac, Cab, Vb, Cbc, Cac, Cbc, Vc),
                                              ncol = 3, nrow = 3),
                                       X_Env[[1]]),
           V_AxE     = compute_var_AxE(matrix(c(Va, Cab, Cac, Cab, Vb, Cbc, Cac, Cbc, Vc),
                                            ncol = 3, nrow = 3),
                                     Env_VCV[[1]]),
           Pi_b  = (b^2 * Env_VCV[[1]][2, 2] - b_se^2) / V_Plas,
           Pi_c  = (c^2 * Env_VCV[[1]][3, 3] - c_se^2) / V_Plas,
           Gamma_a   = Va / V_Add,
           Gamma_b   = Vb * mean(Data[[1]][["Env"]]^2) / V_Add,
           Gamma_c   = Vc * mean(Data[[1]][["Env"]]^4) / V_Add,
           Gamma_ac  = 2 * Cac * mean(Data[[1]][["Env"]]^2) / V_Add,
           Iota_b    = Vb * Env_VCV[[1]][2, 2] / V_AxE,
           Iota_c    = Vc * Env_VCV[[1]][3, 3] / V_AxE,
           V_Res     = Vr,
           V_Tot     = V_Plas + V_Add + V_Res,
           P2        = V_Plas / V_Tot,
           H2_T      = V_Add / V_Tot,
           H2_M      = V_A / V_Tot,
           H2_I      = V_AxE / V_Tot,
           V_Phen    = var(Data[[1]][["Phen"]]),
           .by = c(Scenario, Simulation)) |>
    select(Scenario,
           Simulation,
           P2,
           starts_with("H2"),
           starts_with("V_"),
           starts_with("Pi"),
           starts_with("Gamma"),
           starts_with("Iota"))

## ------------------ Plot the results of the simulations

## Formatting in a longer format
tbl_plot <-
    tbl_part |>
    select(-starts_with("V_Res")) |>
    pivot_longer(-c("Scenario","Simulation", "V_Phen"),
                 names_to = "Parameter", values_to = "Estimate") |>
    # Adding the "true values" computed from the true parameter values
    mutate(
        True_Value = case_when(
            Parameter == "Pi_b" ~ (beta[2]^2 * var_env) / true_V_plas,
            Parameter == "Pi_c" ~ (beta[3]^2 * var_env2) / true_V_plas,
            Parameter == "Gamma_a"  ~ mat_Sigma[1, 1] / true_V_Add,
            Parameter == "Gamma_b"  ~ mat_Sigma[2, 2] * mean_env2 / true_V_Add,
            Parameter == "Gamma_c"  ~ mat_Sigma[3, 3] * mean_env4 / true_V_Add,
            Parameter == "Gamma_ac" ~ 2 * mat_Sigma[1, 3] * mean_env2 / true_V_Add,
            Parameter == "Iota_b" ~ mat_Sigma[2, 2] * var_env / true_V_AxE,
            Parameter == "Iota_c" ~ mat_Sigma[3, 3] * var_env2 / true_V_AxE,
            Parameter == "V_Plas" ~ true_V_plas,
            Parameter == "P2"     ~ true_V_plas / (true_V_plas + true_V_Add + res_sd^2),
            Parameter == "V_Add"  ~ true_V_Add,
            Parameter == "V_A"    ~ true_V_A,
            Parameter == "V_AxE"  ~ true_V_AxE,
            Parameter == "H2_T"   ~ true_V_Add / (true_V_plas + true_V_Add + res_sd^2),
            Parameter == "H2_M"   ~ true_V_A / (true_V_plas + true_V_Add + res_sd^2),
            Parameter == "H2_I"   ~ true_V_AxE / (true_V_plas + true_V_Add + res_sd^2),
            Parameter == "V_Tot"  ~ true_V_plas + true_V_Add + res_sd^2,
            TRUE                            ~ NA_real_),
        Bias = Estimate - True_Value,
        Relative_Bias = Bias / True_Value,
        Parameter = factor(Parameter,
                        levels = c("V_Plas",
                                      "P2",
                                      "Pi_b",
                                      "Pi_c",
                                      "V_Add",
                                      "H2_T",
                                      "Gamma_a",
                                      "Gamma_b",
                                      "Gamma_c",
                                      "Gamma_ac",
                                      "V_A",
                                      "H2_M",
                                      "V_AxE",
                                      "H2_I",
                                      "Iota_b",
                                      "Iota_c",
                                      "V_Tot")),
        Parameter = recode(Parameter,
                           V_Add     = "V[Add]",
                           V_A       = "V[A]",
                           V_AxE     = "V[AxE]",
                           H2_T      = "h[RN]^2",
                           H2_M      = "h^2",
                           H2_I      = "h[I]^2",
                           Pi_b      = "π[Sl]",
                           Pi_c      = "π[Cv]",
                           Gamma_a   = "γ[a]",
                           Gamma_b   = "γ[b]",
                           Gamma_c   = "γ[c]",
                           Gamma_ac  = "γ[ac]",
                           Iota_b    = "ι[b]",
                           Iota_c    = "ι[c]",
                           V_Plas    = "V[Plas]",
                           P2        = "P[RN]^2",
                           V_Tot     = "V[Tot]"),
        Colour = case_when(
            Parameter == "V[Plas]"      ~ "#AA0000",
            Parameter == "P[RN]^2"      ~ "#AA0000",
            str_detect(Parameter, "π")  ~ "#AA0000",
            Parameter == "V[Add]"       ~ "#68349C",
            Parameter == "h[RN]^2"      ~ "#68349C",
            str_detect(Parameter, "γ")  ~ "#68349C",
            Parameter == "V[A]"         ~ "#228B22",
            Parameter == "h^2"          ~ "#228B22",
            Parameter == "V[AxE]"       ~ "#2D489A",
            Parameter == "h[I]^2"       ~ "#2D489A",
            str_detect(Parameter, "ι")  ~ "#2D489A",
            Parameter == "N[e]"         ~ "#2D489A"
        ),
        .by = Scenario) |>
    left_join(tbl_p) |>
    mutate(Scenario_Name = str_glue("paste(N[Env]=={N_Env},', ',N[Gen]=={N_Gen})"))

tbl_true <-
    tbl_plot |>
    select(Scenario, Scenario_Name, Parameter, True_Value, Colour) |>
    filter(Scenario == 1) |>
    distinct() |>
    arrange(Parameter) |>
    filter(!str_detect(Parameter, "^V")) |>
    mutate(Label     = str_c(Parameter, " == ", round(True_Value, digits = 2)),
           Y         = seq(12, -12, length.out = n()),
           X         = 0) #|>
#     mutate(Y = case_when(
#         Scenario == 1 ~ Y,
#         Scenario == 2 ~ Y - 0.25,
#         Scenario == 3 ~ Y + 0.25,
#         Scenario == 4 ~ Y + 0.5,
#     ))

## Generating the graph
p_part <-
    ggplot(tbl_plot |> filter(!str_detect(Parameter, "^V"))) +
    geom_violin(aes(x = Parameter, y = Bias, fill = I(Colour)),
                alpha = 0.5,
                scale = "width") +
    geom_hline(yintercept = 0) +
    stat_summary(aes(x = Parameter, y = Bias),
                 geom = "point",
                 colour = "#444444",
                 size = 2,
                 fun = "mean") +
    facet_grid(Scenario_Name ~ ., scale = "free", space = "free_x", labeller = label_parsed) +
    scale_x_discrete(labels = scales::parse_format()) +
    labs(x = "Variance component", y = "Error") +
    theme(text         = element_text(family = "Linux Biolinum O"),
          axis.text.y  = element_text(size = 14),
          axis.text.x  = element_text(size = 14),
          axis.title   = element_text(size = 16),
          strip.text.x = element_text(size = 20),
          strip.text.y = element_text(size = 12))

p_box <-
    ggplot(tbl_true |>
           filter(!str_detect(Parameter, "^V")) |>
           mutate(Y = 0.5 * Y)) +
    annotate("rect", xmin = -4, xmax = 4, ymin = -7, ymax = 7,
             fill = "white", colour = "#FFCC00", linewidth = 2) +
    geom_text(aes(x = X, y = Y, label = Label, colour = I(Colour)),
              size = 5,
              family = "Linux Biolinum O",
              parse = TRUE) +
    ylim(c(-12, 12)) +
    facet_wrap(~ Scenario_Name, ncol = 1) +
    theme_void() +
    theme(strip.text = element_blank())


## Saving the plot
cairo_pdf(here("Figs/Sim_varpart_bias_cont.pdf"), height = 8, width = 12)
plot((p_part | p_box) + plot_layout(widths = c(4,1)))
dev.off()

## Saving a subset of the plot for the hybrid figure in the main text
tbl_plot_sub <-
    tbl_plot |>
    filter(Scenario == 2) |>
    mutate(Source = "paste('Curve Parameter (continuous)')")
tbl_true_sub <-
    tbl_true |>
    mutate(X         = if_else(str_detect(Parameter, "^[πγι]"),
                               0.15,
                               -0.15)) |>
    mutate(Y         = seq(9, -9, length.out = n()),
           .by       = c(X, Scenario))

p_part_sub <-
    ggplot(tbl_plot_sub |> filter(!str_detect(Parameter, "^V"))) +
    geom_violin(aes(x = Parameter, y = Bias, fill = I(Colour)),
                alpha = 0.5,
                scale = "width") +
    geom_hline(yintercept = 0) +
    stat_summary(aes(x = Parameter, y = Bias),
                 geom = "point",
                 colour = "#444444",
                 size = 2,
                 fun = "mean") +
    facet_grid(Scenario_Name ~ Source, scale = "free", space = "free_x", labeller = label_parsed) +
    scale_x_discrete(labels = scales::parse_format()) +
    labs(x = "Variance component", y = "Error") +
    theme(text         = element_text(family = "Linux Biolinum O"),
          axis.text.y  = element_text(size = 16),
          axis.text.x  = element_text(size = 16),
          axis.title   = element_text(size = 16),
          strip.text.x = element_text(size = 20),
          strip.text.y = element_text(size = 12))

p_box_sub <-
    ggplot(tbl_true_sub |> filter(!str_detect(Parameter, "^V"))) +
    annotate("rect", xmin = -0.31, xmax = 0.31, ymin = -11, ymax = 11,
             fill = "white", colour = "#FFCC00", linewidth = 2) +
    geom_text(aes(x = X, y = Y, label = Label, colour = I(Colour)),
              size = 5,
              family = "Linux Biolinum O",
              parse = TRUE) +
#     ylim(c(-6, 6)) +
    facet_wrap(~ Scenario_Name, ncol = 1) +
    theme_void() +
    theme(strip.text = element_blank())
saveRDS(p_part_sub, file = here("ggplot_part_perfect_cont_sub.rds"))
saveRDS(p_box_sub, file = here("ggplot_box_perfect_cont_sub.rds"))

## Testing for bias
safe_test <- function(vec) {
    test <- safely(wilcox.test, otherwise = NA)(vec)[["result"]]

    if (!all(is.na(test))) test <- test[["p.value"]]

    return(test)
}

test_bias <-
    tbl_plot |>
    select(Scenario, Simulation, Parameter, Bias) |>
    pivot_wider(names_from = Parameter, values_from = Bias) |>
    select(-Simulation) |>
    summarise(P = across(starts_with("V"), safe_test),
              Mean = across(starts_with("V"), mean),
              MSE  = across(starts_with("V"), \(vec) mean(vec^2)),
              .by = c(Scenario))
test_bias

## Looking at the correlation between V_tot and V_phen
tbl_plot |>
    filter(str_detect(Parameter, "Tot")) |>
    select(Scenario, Simulation, V_Phen, Estimate) |>
    rename(V_Tot = Estimate) |>
    summarise(Cor = num(cor(V_Phen, V_Tot), digits = 4),
              Mean_Diff = mean(V_Tot - V_Phen),
              .by = c(Scenario))

## ------------------ Plot along the environment

env_seq <- seq(-2.5, 2.5, length.out = 500)
X <- cbind(1, env_seq, env_seq^2)

tbl_along <-
    tbl_sim |>
    mutate(G = pmap(list(Va, Cab, Cac, Cab, Vb, Cbc, Cac, Cbc, Vc),
                    \(...) matrix(c(...), ncol = 3, nrow = 3))) |>
    select(Scenario, Simulation, G) |>
    mutate(Env = map(G, \(g) { env_seq }),
           V_A = map(G, \(g) {diag(X %*% g %*% t(X))}),
           Gamma_a  = map2(G, V_A, \(g, va) { (g[1, 1]) / va }),
           Gamma_b  = map2(G, V_A, \(g, va) { (g[2, 2] * env_seq^2) / va }),
           Gamma_c  = map2(G, V_A, \(g, va) { (g[3, 3] * env_seq^4) / va }),
           Gamma_ab = map2(G, V_A, \(g, va) { (2 * g[1, 2] * env_seq) / va }),
           Gamma_ac = map2(G, V_A, \(g, va) { (2 * g[1, 3] * env_seq^2) / va }),
           Gamma_bc = map2(G, V_A, \(g, va) { (2 * g[2, 3] * env_seq^3) / va })) |>
    select(-G) |>
    unnest(c(!Simulation, !Scenario)) |>
    pivot_longer(matches("^[VG]"), names_to = "Parameter", values_to = "Value") |>
    summarise(Mean = mean(Value),
          Low  = quantile(Value, probs = 0.025),
          Up   = quantile(Value, probs = 0.975),
          .by = c(Parameter, Scenario, Env))

tbl_along_plot <-
    tbl_along |>
    mutate(Parameter    = as_factor(Parameter) |>
                          recode(V_A      = "V[A,ε]",
                                 Gamma_a  = "γ[a]",
                                 Gamma_b  = "γ[b]",
                                 Gamma_c  = "γ[c]",
                                 Gamma_ac = "γ[ab]",
                                 Gamma_ab = "γ[ac]",
                                 Gamma_bc = "γ[bc]"),
           Category     = if_else(str_detect(Parameter, "γ"),
                                  "γ-decomposition",
                                  "Additive genetic variance",),
           Colour       = set_colours(Parameter)) |>
    filter(Scenario == 1,
           !str_detect(Parameter, "[abc]{2}"))

lty <- c("solid", "solid", "dotted", "twodash")
names(lty) <- levels(tbl_along_plot[["Parameter"]])[1:4]

p_along <-
    ggplot(tbl_along_plot) +
    geom_ribbon(aes(x = Env, ymin = Low, ymax = Up, fill = I(Colour), group = Parameter),
                alpha = 0.3) +
    geom_line(aes(x = Env, y = Mean, colour = I(Colour), linetype = Parameter),
              linewidth = 1.2) +
    facet_wrap(~ Category, scales = "free") +
    scale_linetype_manual(values = lty,
                          breaks = names(lty)[-1],
                          labels = scales::parse_format(),
                          name = NULL) +
    xlab("Environment") + ylab("Value")  +
    theme(legend.position = c(0.93, 0.85))

cairo_pdf(here("Figs/Along_env_cont.pdf"), height = 6, width = 10)
plot(p_along)
dev.off()
