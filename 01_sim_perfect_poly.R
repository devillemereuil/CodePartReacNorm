conflicted::conflict_prefer("pack",     "tidyr")
conflicted::conflict_prefer("unpack",   "tidyr")
conflicted::conflict_prefer("expand",   "tidyr")
conflicted::conflict_prefer("extract",  "tidyr")
conflicted::conflict_prefer("filter",   "dplyr")
conflicted::conflict_prefer("lag",      "dplyr")
conflicted::conflict_prefer("chol2inv", "Matrix")

library(tidyverse)
library(mvtnorm)
library(lme4)
library(furrr)
library(patchwork)
library(here)

n_cores <- min(parallel::detectCores() - 2, 30)
options(mc.cores = n_cores)
plan(multicore)

# Some ggplot options
theme_set(theme_bw())

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

## Setting up the simulation parameters
tbl_p <-
    tribble(~ Scenario, ~ N_Env, ~ N_Gen, ~ N_Ind,
            1, 10, 20, 20,
            2,  4, 20, 20,
            3,  4,  5,  5) |>
    mutate(Env = map(N_Env, \(n) {seq(-2, 2, length.out = n) |>
                                  set_names(str_glue("Env{1:n}"))}),
           V_e = map_dbl(Env, var),
           V_e_sq = map_dbl(Env, \(vec) {var(vec^2)}),
           Vec_rn = map(Env, \(vec) {beta[1] + vec * beta[2] + (vec^2) * beta[3]}))

## ------------------ Functions for the simulations

## Function to simulate the reaction norm, with random variations
# Args: - n: Number of replicates
#       - env: vector of environmental values
#       - beta: average parameters of the reaction norm
#       - mat_vcv: variance-covariance matrix around beta
# Value: A tbl with reaction norm for each env value (N * n_env)
sim_rn <- function(n, env, beta, mat_vcv) {
    # Generating the random variation around the parameters
    params <- rmvnorm(n, beta, mat_vcv)

    # Simulating the reaction norm
    params %*% t(cbind(1, env, env^2))
}

## Function to run the reaction norm model on the data
# Args: - df: A dataset containing the phenotypic data
# Value: The model fit for the data
fit_lme_rn <- function(df) {
    # Formatting the data
    df[["Env_Sq"]] <- df[["Env"]]^2

    # Fitting the model
    lmer(Phen ~ Env + Env_Sq + (1 + Env + Env_Sq | Genotype),
         data = df)
}

## Function to run the character-state model on the data
# Args: - df: A dataset containing the phenotypic data
# Value: The model fit for the data
fit_lme_cs <- function(df) {
    # Formatting the data
    df[["Env_fac"]] <- as_factor(round(df[["Env"]], digits = 1))

    # Fitting the model
    lmer(Phen ~ 0 + Env_fac + (0 + Env_fac | Genotype),
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

## Function to extract the parameters of interest from the CS lme4 fit
# Args: - model: A lme fit yielded by fit_lme_cs
# Values: List of relevant parameters
extract_params_cs <- function(model) {
    # Genotypic VCV
    G <- VarCorr(model)[["Genotype"]]
    attr(G, "stddev") <- NULL
    attr(G, "correlation") <- NULL

    # Average (fixed) effects
    var_p <- var(predict(model, re.form = ~ 0)) -
             sum(diag(vcov(model) %*% cov(model.matrix(model))))

    # Residual variance
    res_var <- attr(VarCorr(model), "sc")^2

    # Eigen values of G
    eig <- eigen(G)$values

    # Returning the list of parameters
    tibble(V_Plas_CS   = var_p,
           V_Add_CS    = mean(diag(G)),
           V_A_CS      = (1/nrow(G)^2) * sum(G),
           G           = list(G),
           Ne_CS       = sum(eig) / max(eig),
           V_Res_CS    = res_var)
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
    as.numeric(t(colMeans(mat_env)) %*% mat_Sigma %*% colMeans(mat_env))
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

## ------------------ Visualising the reaction norm variability

## Simulating a few reaction norms
# Number of reaction norms to simulate
N_vis <- 30
tbl_simvis <-
    tibble(Rep = 1:N_vis,
           # Using the sim_rn() function to simulate the reaction norms
           Reac_Norm = sim_rn(N_vis, tbl_p[["Env"]][[1]], beta, mat_Sigma) |>
                       as_tibble(.name_repair = ~ str_glue("Env{1:tbl_p[['N_Env']][[1]]}"))) |>
    unpack(Reac_Norm) |>
    pivot_longer(contains("Env"),
                 values_to = "Reac_Norm",
                 names_to  = "Env_ref") |>
    # Adding the environmental variable
    mutate(Env = tbl_p[["Env"]][[1]][Env_ref])

## Plotting the reaction norms
p_rn <-
    ggplot(tbl_simvis) +
     geom_line(aes(x = Env, y = Reac_Norm, group = Rep),
               alpha = 0.5,
               size = 1)
cairo_pdf(here("Figs/Sim_RN.pdf"))
plot(p_rn)
dev.off()

## ------------------ Now running the actual simulations

## Generating the reaction norms datasets to be analysed
tbl_sim <-
    # Generating the reaction norms for the total number of simulations
    tbl_p |>
    select(Scenario, N_Env, N_Gen, N_Ind) |>
    reframe(Simulation  = rep(1:n_sim, each = unique(N_Gen)),
            Genotype    = str_pad(rep(1:unique(N_Gen), n_sim), 2, pad = "0"),
            Reac_Norm   = sim_rn(unique(N_Gen) * n_sim,
                                       tbl_p[["Env"]][[unique(Scenario)]],
                                       beta, mat_Sigma) |>
                                as.data.frame(),
            .by = c(Scenario, N_Env, N_Gen, N_Ind)) |>
    reframe(Simulation = rep(Simulation, unique(N_Ind)),
            Genotype = rep(Genotype, unique(N_Ind)),
            N_Env = unique(N_Env),
            Phen = Reac_Norm[rep(1:n(), unique(N_Ind)), ] +
                   rnorm(n      = nrow(Reac_Norm) * ncol(Reac_Norm),
                         mean   = 0,
                         sd     = res_sd),
            .by = Scenario) |>
    unpack(Phen) |>
    pivot_longer(matches("Env[1-9]+"),
                 values_to = "Phen",
                 names_to = "Env") |>
    mutate(Env = tbl_p[["Env"]][[unique(Scenario)]][Env],
           .by = Scenario) |>
    drop_na(Env, Phen) |>
    nest_by(Scenario, Simulation, .key = "Data") |>
    ungroup() |>
    arrange(Simulation)

## Running the models
# Using the function fit_lme_rn, see above
gc()
tbl_sim[["Model_CS"]] <-
    future_map(tbl_sim[["Data"]],
               fit_lme_cs,
               .progress = TRUE,
               .options = furrr_options(seed = TRUE,
                                        chunk_size = n_chunk))
tbl_sim[["Model_RN"]] <-
    future_map(tbl_sim[["Data"]],
               fit_lme_rn,
               .progress = TRUE,
               .options = furrr_options(seed = TRUE,
                                        chunk_size = n_chunk))
gc()

## Saving the design matrix
mat_env <-
    tbl_sim |>
    summarise(M = list(model.matrix(Model_RN[[1]])),
              .by = Scenario)

## Extracting the parameters of interests
# Again see above for function extract_params_rn
tbl_sim <-
    tbl_sim |>
    mutate(Params_RN = map(Model_RN, extract_params_rn) |> list_rbind()) |>
    unpack(Params_RN) |>
    mutate(Params_CS = map(Model_CS, extract_params_cs) |> list_rbind()) |>
    unpack(Params_CS)

## Remove the (heavy) Models columns
tbl_sim <-
    tbl_sim |>
    select(-Model_RN, -Model_CS) |>
    arrange(Scenario)

## Saving the results
saveRDS(tbl_sim, file = here("rn_sims.rds"))
saveRDS(mat_env, file = here("mat_env.rds"))

## Loading previously saved results
tbl_sim <- readRDS(here("rn_sims.rds"))
mat_env <- readRDS(here("mat_env.rds"))

## Compute the variance-covariance matrix of the environmental values
mat_env <-
    mat_env |>
    left_join(tbl_p |> select(Scenario, N_Env, N_Gen)) |>
    mutate(VCV = map(M, \(mat) { ((nrow(mat) - 1) / nrow(mat)) *  cov(mat) }),
           G   = map2(M, N_Env,
                      \(mat, n) { mat[1:n, ] %*% mat_Sigma %*% t(mat[1:n, ]) }),
           N_e = map_dbl(G, \(mat) { eig <- eigen(mat)$values; sum(eig) / max(eig) }))

## ------------------ Graphically checking for bias

## Creating a dataset to study bias
tbl_bias <-
    tbl_sim |>
    select(-Data, -G, -Ne_CS) |>
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
           V_Plas_CS = V_Plas_CS - as.numeric(t(beta) %*%
                                              mat_env[["VCV"]][[unique(Scenario)]] %*%
                                              beta),
           V_AxE_CS   = (V_Add_CS - V_A_CS) -
                       compute_var_AxE(mat_Sigma,
                                      mat_env[["VCV"]][[unique(Scenario)]]),
           V_Add_CS  = V_Add_CS - compute_var_add(mat_Sigma,
                                                  mat_env[["M"]][[unique(Scenario)]]),
           V_A_CS    = V_A_CS - compute_var_A(mat_Sigma,
                                              mat_env[["M"]][[unique(Scenario)]]),
           V_Res_CS  = V_Res_CS - res_sd^2,
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
    facet_wrap(~ Scenario, ncol = 1) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))

cairo_pdf(here("Figs/Sim_Bias.pdf"), width = 9, height = 9)
plot(p_bias)
dev.off()

## ------------------ Phenotypic variance partitioning


## Compute the variance partitions (including "realised" phenotypic variance V_P)
# V_Plas (for Reaction Norm) = V(E(z|e)) and V_Add (for genotype) = E(V(z|e))
# V_P is the realised phenotypic variance and V_tot is the sum of variance components
# (note that V_P ~ V_tot, but strict equality isn't expected)
tbl_part <-
    tbl_sim |>
    mutate(V_Plas_RN    = as.vector(t(c(a, b, c)) %*%
                                    mat_env[["VCV"]][[unique(Scenario)]] %*%
                                    c(a, b, c)) -
                          as.vector(t(c(a_se, b_se, c_se)) %*%
                                    mat_env[["VCV"]][[unique(Scenario)]] %*%
                                    c(a_se, b_se, c_se)),
           V_Add_RN     = compute_var_add(matrix(c(Va, Cab, Cac, Cab, Vb, Cbc, Cac, Cbc, Vc),
                                                  ncol = 3, nrow = 3),
                                           mat_env[["M"]][[unique(Scenario)]]),
           V_A_RN   = Va +
                      2 * Cac * mat_env[["VCV"]][[unique(Scenario)]][2, 2] +
                      Vc * mat_env[["VCV"]][[unique(Scenario)]][2, 2]^2,
           V_AxE_RN  = Vb * mat_env[["VCV"]][[unique(Scenario)]][2, 2] +
                      Vc * mat_env[["VCV"]][[unique(Scenario)]][3, 3],
           Pi_RN_b  = ((b^2 * mat_env[["VCV"]][[unique(Scenario)]][2, 2] - b_se^2) / V_Plas_RN),
           Pi_RN_c  = ((c^2 * mat_env[["VCV"]][[unique(Scenario)]][3, 3] - c_se^2) / V_Plas_RN),
           Gamma_RN_a   = Va / V_Add_RN,
           Gamma_RN_b   = Vb * mean(tbl_p[["Env"]][[unique(Scenario)]]^2) / V_Add_RN,
           Gamma_RN_c   = Vc * mean(tbl_p[["Env"]][[unique(Scenario)]]^4) / V_Add_RN,
           Gamma_RN_ac  = 2 * Cac * mean(tbl_p[["Env"]][[unique(Scenario)]]^2) / V_Add_RN,
           Iota_RN_b    = Vb * mat_env[["VCV"]][[unique(Scenario)]][2, 2] / V_AxE_RN,
           Iota_RN_c    = Vc * mat_env[["VCV"]][[unique(Scenario)]][3, 3] / V_AxE_RN,
           V_Res_RN     = Vr,
           V_Tot_RN     = V_Plas_RN + V_Add_RN + V_Res_RN,
           V_Tot_CS     = V_Plas_CS + V_Add_CS + V_Res_CS,
           V_AxE_CS      = V_Add_CS - V_A_CS,
           V_Phen       = var(Data[[1]][["Phen"]]),
           P2_RN        = V_Plas_RN / V_Tot_RN,
           H2_T_RN      = V_Add_RN / V_Tot_RN,
           H2_M_RN      = V_A_RN / V_Tot_RN,
           H2_I_RN      = V_AxE_RN / V_Tot_RN,
           P2_CS        = V_Plas_CS / V_Tot_CS,
           H2_T_CS      = V_Add_CS / V_Tot_CS,
           H2_M_CS      = V_A_CS / V_Tot_CS,
           H2_I_CS      = V_AxE_CS / V_Tot_CS,
           .by = c(Scenario, Simulation)) |>
    # Some values are aberrant due to (rare) numerical issues with Scenario 3 (low power)
    # The issue comes from an aberrant estimate for the s.e. of b and/or c
    filter(between(Pi_RN_b, 0, 1), between(Pi_RN_c, 0, 1)) |>
    select(Scenario,
           Simulation,
           contains("V_"),
           starts_with("Pi"),
           starts_with("Gamma"),
           starts_with("Iota"),
           starts_with("P2"),
           starts_with("H2"),
           Ne_CS)

# Accordingly subsetting tbl_sim
tbl_sim <-
    inner_join(tbl_sim, tbl_part |> select(Scenario, Simulation))

## ------------------ Efficient number of dimensions

tbl_sim |>
    select(Scenario, Simulation, Ne_CS) |>
    summarise(Mean = mean(Ne_CS),
              Low  = quantile(Ne_CS, probs = 0.025),
              Up   = quantile(Ne_CS, probs = 0.975),
              .by = Scenario)

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
            Parameter == "Pi_RN_b" ~ (beta[2]^2 * mat_env[["VCV"]][[unique(Scenario)]][2, 2]) /
                                     (t(beta) %*% mat_env[["VCV"]][[unique(Scenario)]] %*% beta),
            Parameter == "Pi_RN_c" ~ (beta[3]^2 * mat_env[["VCV"]][[unique(Scenario)]][3, 3]) /
                                     (t(beta) %*% mat_env[["VCV"]][[unique(Scenario)]] %*% beta),
            Parameter == "Gamma_RN_a"  ~ mat_Sigma[1, 1] /
                                         compute_var_add(mat_Sigma, mat_env[["M"]][[unique(Scenario)]]),
            Parameter == "Gamma_RN_b"  ~ mat_Sigma[2, 2] * mean(tbl_p[["Env"]][[unique(Scenario)]]^2) /
                                         compute_var_add(mat_Sigma, mat_env[["M"]][[unique(Scenario)]]),
            Parameter == "Gamma_RN_c"  ~ mat_Sigma[3, 3] * mean(tbl_p[["Env"]][[unique(Scenario)]]^4) /
                                         compute_var_add(mat_Sigma, mat_env[["M"]][[unique(Scenario)]]),
            Parameter == "Gamma_RN_ac" ~ 2 * mat_Sigma[1, 3] * mean(tbl_p[["Env"]][[unique(Scenario)]]^2) /
                                         compute_var_add(mat_Sigma, mat_env[["M"]][[unique(Scenario)]]),
            Parameter == "Iota_RN_b" ~ mat_Sigma[2, 2] * mat_env[["VCV"]][[unique(Scenario)]][2, 2] /
                                         compute_var_AxE(mat_Sigma, mat_env[["VCV"]][[unique(Scenario)]]),
            Parameter == "Iota_RN_c" ~ mat_Sigma[3, 3] * mat_env[["VCV"]][[unique(Scenario)]][3, 3] /
                                         compute_var_AxE(mat_Sigma, mat_env[["VCV"]][[unique(Scenario)]]),
            Parameter == "Ne_CS" ~ mat_env[["N_e"]][[unique(Scenario)]],
            str_detect(Parameter, "V_Plas") ~ t(beta) %*% mat_env[["VCV"]][[unique(Scenario)]] %*% beta,
            str_detect(Parameter, "V_Add")  ~ compute_var_add(mat_Sigma, mat_env[["M"]][[unique(Scenario)]]),
            str_detect(Parameter, "V_AxE")  ~ compute_var_AxE(mat_Sigma, mat_env[["VCV"]][[unique(Scenario)]]),
            str_detect(Parameter, "V_A")    ~ compute_var_A(mat_Sigma, mat_env[["M"]][[unique(Scenario)]]),
#             str_detect(Parameter, "V_Tot")  ~ V_Phen, # NOTE Not the real Vtot explanation
            str_detect(Parameter, "V_Tot")  ~ t(beta) %*% mat_env[["VCV"]][[unique(Scenario)]] %*% beta +
                                              compute_var_add(mat_Sigma, mat_env[["M"]][[unique(Scenario)]]) +
                                              res_sd^2,
            str_detect(Parameter, "P2") ~ t(beta) %*% mat_env[["VCV"]][[unique(Scenario)]] %*% beta /
                                     (t(beta) %*% mat_env[["VCV"]][[unique(Scenario)]] %*% beta +
                                              compute_var_add(mat_Sigma, mat_env[["M"]][[unique(Scenario)]]) +
                                              res_sd^2),
            str_detect(Parameter, "H2_T") ~ compute_var_add(mat_Sigma, mat_env[["M"]][[unique(Scenario)]]) /
                                     (t(beta) %*% mat_env[["VCV"]][[unique(Scenario)]] %*% beta +
                                              compute_var_add(mat_Sigma, mat_env[["M"]][[unique(Scenario)]]) +
                                              res_sd^2),
            str_detect(Parameter, "H2_M") ~ compute_var_A(mat_Sigma, mat_env[["M"]][[unique(Scenario)]]) /
                                     (t(beta) %*% mat_env[["VCV"]][[unique(Scenario)]] %*% beta +
                                              compute_var_add(mat_Sigma, mat_env[["M"]][[unique(Scenario)]]) +
                                              res_sd^2),
            str_detect(Parameter, "H2_I") ~ compute_var_AxE(mat_Sigma, mat_env[["VCV"]][[unique(Scenario)]]) /
                                     (t(beta) %*% mat_env[["VCV"]][[unique(Scenario)]] %*% beta +
                                              compute_var_add(mat_Sigma, mat_env[["M"]][[unique(Scenario)]]) +
                                              res_sd^2),
            TRUE                            ~ NA_real_),
        Source = case_when(
            str_detect(Parameter, "CS") ~ "Character-State",
            str_detect(Parameter, "RN") ~ "Curve Parameter",
            TRUE                        ~ NA_character_),
        Parameter = str_remove(Parameter, "(_RN|_CS)"),
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
                                      "Ne",
                                      "V_Tot")),
        Source = factor(Source,
                        levels = c("Curve Parameter", "Character-State")),
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
                           V_Tot     = "V[Tot]",
                           Ne       = "N[e]"),
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
    mutate(Scenario_Name = recode(Scenario,
                                  `1` = "paste(N[Env]==10,', ',N[Gen]==20,', ',N[Rep]==20)",
                                  `2` = "paste(N[Env]==4,', ',N[Gen]==20,', ',N[Rep]==20)",
                                  `3` = "paste(N[Env]==4,', ',N[Gen]==5,', ',N[Rep]==5)"),
           Source_Fmt = fct_relabel(Source, \(chr) { str_c("paste('", chr, "')") })) |>
    filter(Parameter != "N[e]")

tbl_plot_with_V <- tbl_plot

tbl_plot <-
    tbl_plot |>
    filter(!str_detect(Parameter, "^V\\["))

tbl_true <-
    tbl_plot |>
    select(Scenario, Scenario_Name, Parameter, True_Value, Colour) |>
    distinct() |>
    arrange(Parameter) |>
    mutate(Label     = str_c(Parameter, " == ", round(True_Value, digits = 2)),
           X         = if_else(str_detect(Parameter, "^[πγι]"),
                               0.15,
                               -0.15)) |>
    mutate(Y         = seq(9, -9, length.out = n()),
           .by       = c(X, Scenario))

## Generating the graph
p_part <-
    ggplot(tbl_plot |> filter(Parameter != "N[e]")) +
    geom_violin(aes(x = Parameter, y = Bias, fill = I(Colour)),
                alpha = 0.5,
                scale = "width") +
    geom_hline(yintercept = 0) +
    stat_summary(aes(x = Parameter, y = Bias),
                 geom = "point",
                 colour = "#444444",
                 size = 2,
                 fun = "mean") +
    facet_grid(Scenario_Name ~ Source_Fmt, scale = "free", space = "free_x", labeller = label_parsed) +
    scale_x_discrete(labels = scales::parse_format()) +
    labs(x = "Variance component", y = "Error") +
    theme(text         = element_text(family = "Linux Biolinum O"),
          axis.text.y  = element_text(size = 14),
          axis.text.x  = element_text(size = 14),
          axis.title   = element_text(size = 16),
          strip.text.x = element_text(size = 20),
          strip.text.y = element_text(size = 12))

p_box <-
    ggplot(tbl_true) +
    annotate("rect", xmin = -0.3, xmax = 0.3, ymin = -11, ymax = 11,
             fill = "white", colour = "#FFCC00", linewidth = 2) +
    geom_text(aes(x = X, y = Y, label = Label, colour = I(Colour)),
              size = 5,
              family = "Linux Biolinum O",
              parse = TRUE) +
#     ylim(c(-8, 8)) +
    facet_wrap(~ Scenario_Name, ncol = 1) +
    theme_void() +
    theme(strip.text = element_blank())


## Saving the plot
cairo_pdf(here("Figs/Sim_varpart_bias.pdf"), height = 8, width = 14)
plot((p_part | p_box) + plot_layout(widths = c(4,1)))
dev.off()


## Saving a sub-graph for the hybrid graph in the main text
tbl_plot_sub <-
    tbl_plot |>
    filter(Scenario %in% c(1, 3)) |>
    mutate(Source_Fmt = if_else(Source_Fmt == "paste('Curve Parameter')",
                                "paste('Curve Parameter (discrete)')",
                                Source_Fmt) |>
                        fct_rev()) |>
    filter(Parameter != "N[e]")
tbl_true_sub <-
    tbl_true |>
    filter(Scenario %in% c(1, 3))

p_part_sub <-
    ggplot(tbl_plot_sub) +
    geom_violin(aes(x = Parameter, y = Bias, fill = I(Colour)),
                alpha = 0.5,
                scale = "width") +
    geom_hline(yintercept = 0) +
    stat_summary(aes(x = Parameter, y = Bias),
                 geom = "point",
                 colour = "#444444",
                 size = 2,
                 fun = "mean") +
    facet_grid(Scenario_Name ~ Source_Fmt, scale = "free", space = "free_x", labeller = label_parsed) +
    scale_x_discrete(labels = scales::parse_format()) +
    labs(x = "Variance component", y = "Error") +
    theme(text         = element_text(family = "Linux Biolinum O"),
          axis.text.y  = element_text(size = 16),
          axis.text.x  = element_text(size = 16),
          axis.title   = element_text(size = 16),
          strip.text.x = element_text(size = 20),
          strip.text.y = element_text(size = 12))
p_box_sub <-
    ggplot(tbl_true_sub) +
    annotate("rect", xmin = -0.31, xmax = 0.31, ymin = -11, ymax = 11,
             fill = "white", colour = "#FFCC00", linewidth = 2) +
    geom_text(aes(x = X, y = Y, label = Label, colour = I(Colour)),
              size = 5,
              family = "Linux Biolinum O",
              parse = TRUE) +
    #     ylim(c(-8, 8)) +
    facet_wrap(~ Scenario_Name, ncol = 1) +
    theme_void() +
    theme(strip.text = element_blank())
saveRDS(p_part_sub, file = here("ggplot_part_perfect_disc_sub.rds"))
saveRDS(p_box_sub, file = here("ggplot_box_perfect_disc_sub.rds"))


## Testing for bias
safe_test <- function(vec) {
    test <- safely(wilcox.test, otherwise = NA)(vec)[["result"]]

    if (!all(is.na(test))) test <- test[["p.value"]]

    return(test)
}

test_bias <-
    tbl_plot |>
    select(Scenario, Simulation, Parameter, Bias, Source) |>
    pivot_wider(names_from = Parameter, values_from = Bias) |>
    select(-Simulation) |>
    summarise(Mean = across(matches("^[PHN]"), mean),
              P = across(matches("^[PHN]"), safe_test),
              MSE  = across(matches("^[PHN]"), \(vec) mean(vec^2)),
              .by = c(Scenario, Source))
test_bias[["P"]]
test_bias[["Mean"]]


## Looking at how well the approaches are correlated
tbl_plot_with_V |>
    filter(!(str_detect(Parameter, "π|γ"))) |>
    select(Scenario, Simulation, Parameter, Estimate, Source) |>
    pivot_wider(names_from = Source, values_from = Estimate) |>
    select(-Simulation) |>
    summarise(Correlation = num(cor(`Character-State`, `Curve Parameter`), digits = 4),
              .by = c(Scenario, Parameter))

## Looking at the correlation between V_tot and V_phen
tbl_plot_with_V |>
    filter(str_detect(Parameter, "Tot")) |>
    select(Scenario, Simulation, V_Phen, Source, Estimate) |>
    rename(V_Tot = Estimate) |>
    summarise(Cor = num(cor(V_Phen, V_Tot), digits = 4),
              .by = c(Scenario, Source))
