conflicted::conflict_prefer("pack",     "tidyr")
conflicted::conflict_prefer("unpack",   "tidyr")
conflicted::conflict_prefer("expand",   "tidyr")
conflicted::conflict_prefer("extract",  "tidyr")
conflicted::conflict_prefer("filter",   "dplyr")
conflicted::conflict_prefer("lag",      "dplyr")
conflicted::conflict_prefer("chol2inv", "Matrix")

#########################################################################################################
##                                                                                                     ##
##                          Code to simulate and analyse a polynomial reaction norm                    ##
##                                      Corresponds to Figure 3                                        ##
##                                    Pierre de Villemereuil (2023)                                    ##
##                                                                                                     ##
#########################################################################################################

here::i_am("01_sim_perfect_poly.R")

library(tidyverse)
library(mvtnorm)
library(lme4)
library(furrr)
library(patchwork)
library(here)

options(mc.cores = parallel::detectCores() - 2)
plan(multicore)

# Some ggplot options
theme_set(theme_bw())

## ------------------ Setting parameters

# Setting up a seed for the simulation
set.seed(2021)

## Number of simulations
n_sim <- 1000

## Number of genotypes and individuals
n_gen <- 20
n_ind <- 20

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

## Environment
n_env <- 10
env <- seq(-2, 2, length.out = n_env)
# env <- c(sort(rexp(n_env)))
names(env) <- str_glue("Env{1:n_env}")

## Compute the environmental variance and squared variance
V_e     <- var(env)
V_e_sq  <- var(env^2)

## Compute average reaction norm
vec_rn <- beta[1] + env * beta[2] + (env^2) * beta[3]

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
    out   <- summary(model)[["coefficients"]]
    var_p <- var(predict(model, re.form = ~ 0)) - mean(out[, "Std. Error"]^2)

    # Residual variance
    res_var <- attr(VarCorr(model), "sc")^2

    # Returning the list of parameters
    tibble(V_Plas_CS   = var_p,
           V_Gen_CS    = mean(diag(G)),
           V_Res_CS    = res_var)
}

## Function to compute the expected genotypic variance
# Args: - vcv: the variance-covariance matrix between the RN parameters
#       - mat_env: the "design matrix" of the reaction norm
# Values: The expected genotypic variance
compute_var_geno <- function(vcv, mat_env) {
    # Sum of the element-wise multiplication of products with mat_env
    # with the elements of the variance-covariance matrix vcv
    # Equivalent (but faster) than averaging over the t(E[i, ]) %*% vcv %*% E[i, ] products
    sum(diag((1/nrow(mat_env)) * (t(mat_env) %*% mat_env) %*% vcv))
}

## ------------------ Visualising the reaction norm variability

## Simulating a few reaction norms
# Number of reaction norms to simulate
N_vis <- 30
tbl_simvis <-
    tibble(Rep = 1:N_vis,
           # Using the sim_rn() function to simulate the reaction norms
           Reac_Norm = sim_rn(N_vis, env, beta, mat_Sigma) |>
                       as_tibble(.name_repair = ~ str_glue("Env{1:n_env}"))) |>
    unpack(Reac_Norm) |>
    pivot_longer(contains("Env"),
                 values_to = "Reac_Norm",
                 names_to  = "Env_ref") |>
    # Adding the environmental variable
    mutate(Env = env[Env_ref])

## Plotting the reaction norms
p_rn <-
    ggplot(tbl_simvis) +
     geom_line(aes(x = Env, y = Reac_Norm, group = Rep),
               alpha = 0.5,
               size = 1)
cairo_pdf(here("Sim_RN.pdf"))
plot(p_rn)
dev.off()

## ------------------ Now running the actual simulations

## Generating the reaction norms datasets to be analysed
tbl_sim <-
    # Generating the n_gen reaction norms for the n_sim simulations
    tibble(Simulation   = rep(1:n_sim, each = n_gen),
           Genotype     = str_pad(rep(1:n_gen, n_sim), 2, pad = "0"),
           Reac_Norm    = sim_rn(n_gen * n_sim, env, beta, mat_Sigma)) |>
    # Now simulating n_ind individual phenotypes as noise around the reaction norms
    group_by(Simulation, Genotype) |>
    # Complex code, but this is just adding white Gaussian noise around each point of the RN
    summarise(Phen = matrix(rep(Reac_Norm, n_ind),
                            nrow  = n_ind,
                            ncol  = n_env,
                            byrow = TRUE) +
                     rmvnorm(n      = n_ind,
                             mean   = rep(0, n_env),
                             sigma  = (res_sd^2) * diag(n_env))) |>
    mutate(Phen = as_tibble(Phen, .name_repair = ~ str_glue("Env{1:n_env}"))) |>
    unpack(Phen) |>
    ungroup() |>
    pivot_longer(contains("Env"),
                 values_to = "Phen",
                 names_to = "Env") |>
    mutate(Env = env[Env]) |>
    nest_by(Simulation, .key = "Data") |>
    ungroup()

## Running the models
# Using the function fit_lme_rn, see above
gc()
tbl_sim[["Model_CS"]] <-
    future_map(tbl_sim[["Data"]],
               fit_lme_cs,
               .progress = TRUE,
               .options = furrr_options(seed = TRUE))
tbl_sim[["Model_RN"]] <-
    future_map(tbl_sim[["Data"]],
               fit_lme_rn,
               .progress = TRUE,
               .options = furrr_options(seed = TRUE))
gc()

## Saving the design matrix
mat_env <- model.matrix(tbl_sim[["Model_RN"]][[1]])

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
    select(-Model_RN, -Model_CS)

## Saving the results
saveRDS(tbl_sim, file = here("rn_sims.rds"))
saveRDS(mat_env, file = here("mat_env.rds"))

## Loading previously saved results
tbl_sim <- readRDS(here("rn_sims.rds"))
mat_env <- readRDS(here("mat_env.rds"))

## Compute the variance-covariance matrix of the environmental values
vcv_env <- cov(mat_env)

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
           V_Plas_CS = V_Plas_CS - as.numeric(t(beta) %*% vcv_env %*% beta),
           V_Gen_CS  = V_Gen_CS - compute_var_geno(mat_Sigma, mat_env),
           V_Res_CS  = V_Res_CS - res_sd^2) |>
    select(-Simulation) |>
    pivot_longer(everything(), names_to = "Parameter", values_to = "Bias") |>
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
                 fun = "mean")

cairo_pdf(here("Sim_Bias.pdf"), width = 7, height = 6)
plot(p_bias)
dev.off()

## ------------------ Phenotypic variance partitioning

## Compute the variance partitions (including "realised" phenotypic variance V_P)
# V_Plas (for Reaction Norm) = V(E(z|e)) and V_Gen (for genotype) = E(V(z|e))
# V_P is the realised phenotypic variance and V_tot is the sum of variance components
# (note that V_P ~ V_tot, but strict equality isn't expected)
tbl_part <-
    tbl_sim |>
    rowwise() |>
    mutate(V_Plas_RN    = as.vector(t(c(a, b, c)) %*% vcv_env %*% c(a, b, c)) -
                          as.vector(t(c(a_se, b_se, c_se)) %*% vcv_env %*% c(a_se, b_se, c_se)),
           V_Gen_RN     = compute_var_geno(matrix(c(Va, Cab, Cac, Cab, Vb, Cbc, Cac, Cbc, Vc),
                                                  ncol = 3, nrow = 3),
                                           mat_env),
           Pi_RN_b  = ((b^2 * vcv_env[2, 2] - b_se^2) / V_Plas_RN),
           Pi_RN_c  = ((c^2 * vcv_env[3, 3] - c_se^2) / V_Plas_RN),
           Gamma_RN_a   = Va / V_Gen_RN,
           Gamma_RN_b   = Vb * mean(env^2) / V_Gen_RN,
           Gamma_RN_c   = Vc * mean(env^4) / V_Gen_RN,
           Gamma_RN_ac  = 2 * Cac * mean(env^2) / V_Gen_RN,
           V_Res_RN     = Vr,
           V_Tot_RN     = V_Plas_RN + V_Gen_RN + V_Res_RN,
           V_Tot_CS     = V_Plas_CS + V_Gen_CS + V_Res_CS,
           V_Phen       = var(Data[["Phen"]])) |>
    ungroup() |>
    select(Simulation, contains("V_"), starts_with("Pi"), starts_with("Gamma"))

## ------------------ Plot the results of the simulations

## Formatting in a longer format
tbl_plot <-
    tbl_part |>
    select(-starts_with("V_Res")) |>
    pivot_longer(-c("Simulation", "V_Phen"), names_to = "Parameter", values_to = "Estimate") |>
    # Adding the "true values" computed from the true parameter values
    mutate(
        True_Value = case_when(
            Parameter == "Pi_RN_b" ~ (beta[2]^2 * vcv_env[2, 2]) / (t(beta) %*% vcv_env %*% beta),
            Parameter == "Pi_RN_c" ~ (beta[3]^2 * vcv_env[3, 3]) / (t(beta) %*% vcv_env %*% beta),
            Parameter == "Gamma_RN_a"  ~ mat_Sigma[1, 1] / compute_var_geno(mat_Sigma, mat_env),
            Parameter == "Gamma_RN_b"  ~ mat_Sigma[2, 2] * mean(env^2) / compute_var_geno(mat_Sigma, mat_env),
            Parameter == "Gamma_RN_c"  ~ mat_Sigma[3, 3] * mean(env^4) / compute_var_geno(mat_Sigma, mat_env),
            Parameter == "Gamma_RN_ac" ~ 2 * mat_Sigma[1, 3] * mean(env^2) / compute_var_geno(mat_Sigma, mat_env),
            str_detect(Parameter, "V_Plas") ~ t(beta) %*% vcv_env %*% beta,
            str_detect(Parameter, "V_Gen")  ~ compute_var_geno(mat_Sigma, mat_env),
#             str_detect(Parameter, "V_Tot")  ~ V_Phen, # NOTE Not the real Vtot explanation
            str_detect(Parameter, "V_Tot")  ~ t(beta) %*% vcv_env %*% beta +
                                              compute_var_geno(mat_Sigma, mat_env) +
                                              res_sd^2,
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
                                      "Pi_b",
                                      "Pi_c",
                                      "V_Gen",
                                      "Gamma_a",
                                      "Gamma_b",
                                      "Gamma_c",
                                      "Gamma_ac",
                                      "V_Tot")),
        Source = factor(Source,
                        levels = c("Curve Parameter", "Character-State")),
        Parameter = recode(Parameter,
                           V_Gen     = "V[Gen]",
                           Pi_b      = "π[b]",
                           Pi_c      = "π[c]",
                           Gamma_a   = "γ[a]",
                           Gamma_b   = "γ[b]",
                           Gamma_c   = "γ[c]",
                           Gamma_ac  = "γ[ac]",
                           V_Plas    = "V[Plas]",
                           V_Tot     = "V[Tot]")
    )

## Generating the graph
p_part <-
    ggplot(tbl_plot) +
    geom_violin(aes(x = Parameter, y = Relative_Bias), fill = "grey", scale = "width") +
    geom_hline(yintercept = 0) +
#     geom_point(aes(x = Parameter, y = True_Value), colour = "blue", shape = "square", size = 2) +
    stat_summary(aes(x = Parameter, y = Relative_Bias),
                 geom = "point",
                 colour = "red",
                 size = 2,
                 fun = "mean") +
    facet_grid(~ Source, scale = "free_x", space = "free_x") +
    scale_x_discrete(labels = scales::parse_format()) +
    labs(x = "Variance component", y = "Relative Error") +
    theme(text         = element_text(family = "Linux Biolinum O"),
          axis.text.y  = element_text(size = 14),
          axis.text.x  = element_text(size = 14),
          axis.title   = element_text(size = 16),
          strip.text   = element_text(size = 20))

tbl_true <-
    tbl_plot |>
    select(Parameter, True_Value) |>
    distinct() |>
    arrange(Parameter) |>
    mutate(Label     = str_c(Parameter, " == ", round(True_Value, digits = 2)),
           X         = 0,
           Y         = seq(2, -2, length.out = n()))

p_box <-
    ggplot(tbl_true) +
    annotate("rect", xmin = -0.3, xmax = 0.3, ymin = -2.5, ymax = 2.5,
             fill = "white", colour = "#FFCC00", linewidth = 2) +
    geom_text(aes(x = X, y = Y, label = Label),
              size = 6,
              family = "Linux Biolinum O",
              parse = TRUE) +
    ylim(c(-5, 5)) +
    theme_void()

## Saving the plot
cairo_pdf(here("Sim_varpart_bias.pdf"), height = 6, width = 10)
plot(p_part | p_box) + plot_layout(widths = c(4,1))
dev.off()

## Testing for bias
safe_rank_test <- function(vec) {
    test <- safely(wilcox.test, otherwise = NA)(vec)[["result"]]

    if (!all(is.na(test))) test <- test[["p.value"]]

    return(test)
}

test_bias <-
    tbl_plot |>
    select(Simulation, Parameter, Relative_Bias, Source) |>
    pivot_wider(names_from = Parameter, values_from = Relative_Bias) |>
    select(-Simulation) |>
    group_by(Source) |>
    summarise(across(contains("hat(V)"), safe_rank_test))
test_bias

## Looking at how well the approaches are correlated
tbl_plot |>
    filter(!(Parameter %in% c("hat(V)[b]", "hat(V)[c]"))) |>
    select(Simulation, Parameter, Estimate, Source) |>
    pivot_wider(names_from = Source, values_from = Estimate) |>
    select(-Simulation) |>
    group_by(Parameter) |>
    summarise(Correlation = num(cor(`Character-State`, `Curve Parameter`), digits = 4))

## Looking at the correlation between V_tot and V_phen
tbl_plot |>
    filter(str_detect(Parameter, "Tot")) |>
    select(Simulation, V_Phen, Source, Estimate) |>
    rename(V_Tot = Estimate) |>
    group_by(Source) |>
    summarise(Cor = num(cor(V_Phen, V_Tot), digits = 4))
