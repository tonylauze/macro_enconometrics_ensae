##############################################
### REPLICATION CODE FOR MACROECONOMETRICS ###
##############################################

# Packages ----------------------------------------------------------------

install.packages(c("BVAR", "dplyr", "tsibble", "tidyverse"))
install.packages("MinnesotaPrior")
library(BVAR)
library(dplyr)
library(tsibble)
library(tidyverse)
library(MinnesotaPrior)  # For Minnesota priors

# Import clean data -------------------------------------------------------

data <- readRDS("~/work/macro_enconometrics_ensae/1_clean_data/data_clean.rds")

# export data in csv
# write.csv(data, "data_clean.csv", row.names = FALSE)

# convert in time series
data_ts <- ts(data %>% select(-date), start = c(2002, 1), frequency = 4)

# name of variables used in the model
vars <- c("hicp_ecb", "gdp_growth_log", "gscpi", "avg_shadow_rate", "log_deflated_price", "deficit")

# convert data in matrix
data_matrix <- as.matrix(data[, c(
  "avg_shadow_rate", 
  "log_deflated_price", 
  "deficit", 
  "hicp_ecb", 
  "gscpi", 
  "gdp_growth_log"
)])

# Define priors and parameters --------------------------------------------

lags <- 4

# Define the standard Minnesota prior
mn_priors <- bv_mn(
  lambda = bv_lambda(mode = 0.2, sd = 0.4, min = 0.0001, max = 5),  # Overall tightness
  alpha = bv_alpha(mode = 2, sd = 0.25, min = 1, max = 3),         # Lag decay
  psi = bv_psi(scale = 0.004, shape = 0.004)                      # Cross-variable variance
)

# Combine the Minnesota prior with other settings into a prior object
priors <- bv_priors(
  hyper = "auto",  # Treat Minnesota parameters as hyperparameters
  mn = mn_priors
)

# Estimate the model ------------------------------------------------------

# The BVAR package uses Giannone, Lenza and Primiceri (2015) method by default.

# However, in our article, it's the Lenza and Primiceri (2022) model that is used
# Lenza and Primiceri (2022) model takes into account for COVID additionnal volatility

bvar_model <- bvar(
  data = data_matrix,
  lags = lags,
  n_draw = 10000L,         # Total draws
  n_burn = 5000L,          # Burn-in period
  n_thin = 1L,             # Thinning factor
  priors = priors,         # Priors
  verbose = TRUE
)

# Plot residuals
plot(residuals(bvar_model, type = "mean"))

# Plot de posterior density of lambda
lambda_density <- density(bvar_model, vars = "lambda")
plot(lambda_density, main = "Densité a posteriori de lambda")

# Define sign restrictions matrix -----------------------------------------

# definition of the same sign restriction matrix as in the article
sign_restriction_matrix <- matrix(c(
  1, 1, 1, -1, -1, -1,    # GDP growth
  1, 1, 1, 1, 1, 1,       # Inflation
  NA, NA, NA, 0, 1, NA,    # GSCPI
  NA, -1, 1, NA, NA, NA,   # Shadow rate
  NA, NA, NA, -1, 0, 1,   # Real oil price
  1, -1, -1, NA, NA, NA    # Primary deficit
), nrow = 6, byrow = TRUE)

colnames(sign_restriction_matrix) <- c("Fiscal policy", "Monetary policy", 
                                       "Non-policy demand", "Cost push", 
                                       "Supply chain", "Oil supply")

rownames(sign_restriction_matrix) <- c("GDP growth", "Inflation", "GSCPI", 
                                       "Shadow rate", "Real oil price", 
                                       "Primary deficit")

print(sign_restriction_matrix)


# Impulse response functions ----------------------------------------------

# define the settings
irf_settings <- bv_irf(
  horizon = 4,            # Horizon of 12 quarters
  identification = TRUE,   # Enable identification
  sign_restr = sign_restriction_matrix,
  sign_lim=3000)

# Calculate IRFs (NOT WORKING)
irf_results <- irf(bvar_model, irf_settings)

# Plot IRFs
plot(irf(bvar_model), area = TRUE,
     vars_impulse = c("avg_shadow_rate", "deficit"), vars_response = c(1:2, 6))


# Forecast ----------------------------------------------------------------

predict(bvar_model) <- predict(bvar_model, horizon = 16, conf_bands = c(0.05, 0.16))

plot(predict(bvar_model), area = TRUE, t_back = 32,
     vars = c("avg_shadow_rate", "deficit", "gscpi"))


# Historical decomposition ------------------------------------------------

# Historical decomposition
hist_decomp_results <- hist_decomp(bvar_model)

