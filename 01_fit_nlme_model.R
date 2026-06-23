# Piecewise-exponential NLME model of anti-Spike antibody kinetics, fitted by SAEM.
# Run once per vaccination-group dataset.
#
# sample1 : empirical antibody data, wide format, one row per participant.
#           Required columns:
#             BI                              breakthrough-infection date (blank/NA if none)
#             1stVaccine, 2ndVaccine, 3rdVaccine   dose dates
#             Sab1Val..Sab9Val                antibody measurements (BAU/mL)
#             Sab1Date..Sab9Date              measurement dates

library(nlmixr2)
library(rxode2)
library(dplyr)
library(tidyr)
library(readr)

infile  <- "sample1.csv"     # <- input data (see header)
out_fit <- "fit_model.rds"   # <- output: fitted model object
out_obs <- "obs_long.csv"    # <- output: long-format observations used for the fit

## ---- load and reshape wide -> long ---------------------------------
dat <- read_csv(infile, show_col_types = FALSE)
dat$id <- seq_len(nrow(dat))

sab_val_long <- dat %>%
  select(id, BI, `1stVaccine`, `2ndVaccine`, `3rdVaccine`, matches("^Sab\\d+Val$")) %>%
  pivot_longer(matches("^Sab\\d+Val$"), names_to = "idx",
               names_pattern = "^Sab(\\d+)Val$", values_to = "DV")

sab_date_long <- dat %>%
  select(id, matches("^Sab\\d+Date$")) %>%
  pivot_longer(matches("^Sab\\d+Date$"), names_to = "idx",
               names_pattern = "^Sab(\\d+)Date$", values_to = "meas_date")

dat_long <- sab_val_long %>%
  left_join(sab_date_long, by = c("id", "idx")) %>%
  mutate(
    meas_date = as.Date(meas_date),
    TIME = as.numeric(meas_date - as.Date(`1stVaccine`)),       # days since 1st dose
    t2   = as.numeric(as.Date(`2ndVaccine`) - as.Date(`1stVaccine`)),
    t3   = as.numeric(as.Date(`3rdVaccine`) - as.Date(`1stVaccine`)),
    tau_raw = as.numeric(as.Date(BI) - as.Date(`1stVaccine`)),
    hasBI = ifelse(is.na(tau_raw), 0, 1),
    tau   = ifelse(is.na(tau_raw), 1e9, tau_raw),               # 1e9 = no BI
    DV    = as.numeric(DV)
  ) %>%
  filter(!is.na(DV), !is.na(TIME), !is.na(t2), !is.na(t3), TIME >= 0) %>%
  select(id, TIME, DV, t2, t3, tau, hasBI)

## ---- NLME model: piecewise-exponential boosting/waning between events
antibody_model <- function() {
  ini({
    log_A0 <- log(5)
    log_k1 <- log(0.05)
    log_k2 <- log(0.10)
    log_k3 <- log(0.01)
    log_k4 <- log(0.50)
    log_k5 <- log(0.01)
    log_k6 <- log(0.50)
    log_k7 <- log(0.005)
    log_phi2 <- log(21)
    log_phi3 <- log(21)
    log_phit <- log(10)

    eta_A0 ~ 2.0
    eta_k1 ~ 1.0
    eta_k2 ~ 0.5
    eta_k3 ~ 1.0
    eta_k4 ~ 0.5
    eta_k5 ~ 1.0
    eta_k6 ~ 1.5
    eta_k7 ~ 1.5
    eta_phi2 ~ 0.5
    eta_phi3 ~ 0.5
    eta_phit ~ 0.5

    sigma <- 0.2
  })
  model({
    Ab(0) = exp(log_A0 + eta_A0)
    k1 <- exp(log_k1 + eta_k1)
    k2 <- exp(log_k2 + eta_k2)
    k3 <- exp(log_k3 + eta_k3)
    k4 <- exp(log_k4 + eta_k4)
    k5 <- exp(log_k5 + eta_k5)
    k6 <- exp(log_k6 + eta_k6)
    k7 <- exp(log_k7 + eta_k7)
    phi2 <- exp(log_phi2 + eta_phi2)
    phi3 <- exp(log_phi3 + eta_phi3)
    phit <- exp(log_phit + eta_phit)

    kRate <-
      ifelse(t <= t2, k1,
             ifelse(t <= t2 + phi2, k2,
                    ifelse(t <= t3, -k3,
                           ifelse(t <= t3 + phi3, k4,
                                  ifelse((hasBI == 0) | (t - tau <= 0), -k5,
                                         ifelse((t - tau) <= phit, k6, -k7))))))
    d/dt(Ab) = kRate * Ab
    DV = Ab
    DV ~ prop(sigma)
  })
}

## ---- fit (SAEM) and save -------------------------------------------
fit <- nlmixr2(antibody_model, dat_long, est = "saem",
               control = saemControl(nBurn = 500, nEm = 4500, print = 50))

saveRDS(fit, out_fit)
write.csv(dat_long, out_obs, row.names = FALSE)
