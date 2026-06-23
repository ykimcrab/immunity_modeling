# Simulation-recovery of individual-level random effects (identifiability check).
# The fitted population model is treated as the data-generating truth: sparse
# observations are simulated at each participant's real sampling times, then each
# subject's random effects are re-estimated by MAP with the population FIXED.
# No SAEM refit. Run once per vaccination-group dataset.
#
# sample1 : empirical antibody data, wide format (same file used as input to 01)
# sample2 : fitted model object (.rds) produced by 01_fit_nlme_model.R

suppressMessages({
  library(nlmixr2); library(rxode2)
  library(dplyr); library(tidyr); library(readr)
})

data_file <- "sample1.csv"               # <- empirical antibody data (wide)
fit_file  <- "sample2.rds"               # <- fitted model object (from 01)
out_raw   <- "simrecovery_raw.csv"       # <- output: per-subject true vs estimated
out_sum   <- "simrecovery_perparam.csv"  # <- output: per-parameter recovery summary

ERR_MODEL <- "prop"   # generating + MAP error model (proportional)
N_REP     <- 10       # simulated cohorts (total fits = N_REP * n_subjects)
MIN_OBS   <- 3
FLOOR     <- 1e-3
set.seed(20260618)

PAR <- c("A0","k1","k2","k3","k4","k5","k6","k7","phi2","phi3","phit")
P   <- length(PAR)
IDX_BI    <- seq_len(P)
IDX_NONBI <- which(!PAR %in% c("k6","k7","phit"))   # non-BI subjects lack post-BI params
LNAMES <- paste0("l", PAR)

## closed-form piecewise-exponential trajectory (exact; no ODE integrator)
solve_ab <- function(lp, cov, times){
  A0  <- exp(lp[["lA0"]])
  k1  <- exp(lp[["lk1"]]); k2 <- exp(lp[["lk2"]]); k3 <- exp(lp[["lk3"]])
  k4  <- exp(lp[["lk4"]]); k5 <- exp(lp[["lk5"]])
  k6  <- exp(lp[["lk6"]]); k7 <- exp(lp[["lk7"]])
  phi2<- exp(lp[["lphi2"]]); phi3 <- exp(lp[["lphi3"]]); phit <- exp(lp[["lphit"]])
  t2  <- cov[["t2"]]; t3 <- cov[["t3"]]; tau <- cov[["tau"]]; bi <- cov[["hasBI"]]
  rate_at <- function(t){
    if (t <= t2)              return( k1)
    if (t <= t2 + phi2)       return( k2)
    if (t <= t3)              return(-k3)
    if (t <= t3 + phi3)       return( k4)
    if (bi == 0 || t <= tau)  return(-k5)
    if (t <= tau + phit)      return( k6)
    -k7
  }
  bps <- c(t2, t2 + phi2, t3, t3 + phi3)
  if (bi == 1) bps <- c(bps, tau, tau + phit)
  bps <- sort(unique(bps))
  out <- numeric(length(times))
  for (i in seq_along(times)){
    th <- times[i]
    if (!is.finite(th)) { out[i] <- NA_real_; next }
    if (th <= 0)        { out[i] <- A0;       next }
    grid <- sort(unique(c(0, bps[bps > 0 & bps < th], th)))
    A <- A0
    for (j in seq_len(length(grid) - 1)){
      a <- grid[j]; b <- grid[j + 1]
      A <- A * exp(rate_at((a + b) / 2) * (b - a))
    }
    out[i] <- A
  }
  out
}

## population truth from the fit: fixed effects, diagonal Omega, sigma
get_pop <- function(fit){
  fe <- tryCatch(fixef(fit), error = function(e) NULL)
  if (is.null(fe)) { pf <- as.data.frame(fit$parFixedDf)
                     fe <- setNames(pf$Estimate, rownames(pf)) }
  pick <- function(nm){
    for (cand in c(paste0("log_", nm), nm, paste0("l", nm)))
      if (cand %in% names(fe)) return(as.numeric(fe[[cand]]))
    NA_real_ }
  lpop <- setNames(vapply(PAR, pick, numeric(1)), LNAMES)
  if (any(is.na(lpop)))
    stop("Cannot map population params; names(fixef): ", paste(names(fe), collapse=", "))
  om <- fit$omega
  od <- setNames(diag(om), gsub("^eta_|^log_", "", rownames(om)))
  omega_diag <- as.numeric(od[PAR]); names(omega_diag) <- PAR
  omega_diag[is.na(omega_diag) | omega_diag <= 0] <- 1e-6
  pf <- as.data.frame(fit$parFixedDf)
  sig_row <- grep("sigma|prop|add", rownames(pf), ignore.case = TRUE)
  sig <- if (length(sig_row)) as.numeric(pf$Estimate[sig_row[1]]) else 0.2
  list(lpop = lpop, omega_diag = omega_diag, sig = sig)
}

## MAP / empirical-Bayes estimation of one subject's random effects
fit_eta <- function(idx_sub, cov, tobs, yobs, pop){
  lpop <- pop$lpop; od <- pop$omega_diag; sig <- pop$sig
  obj <- function(es){
    ef <- rep(0, P); ef[idx_sub] <- es
    lp <- lpop + ef
    f  <- solve_ab(lp, cov, tobs)
    if (any(!is.finite(f))) return(1e12)
    f  <- pmax(f, 1e-8)
    if (ERR_MODEL == "lnorm"){
      r  <- (log(yobs) - log(f)) / sig
      dd <- sum(0.5 * r^2) + length(yobs) * log(sig)
    } else {
      sdv <- sig * f
      r   <- (yobs - f) / sdv
      dd  <- sum(0.5 * r^2 + log(sdv))
    }
    pp <- sum(0.5 * es^2 / od[idx_sub])
    dd + pp
  }
  op <- tryCatch(optim(rep(0, length(idx_sub)), obj,
                       method = "BFGS", control = list(maxit = 200)),
                 error = function(e) NULL)
  ef <- rep(0, P)
  if (!is.null(op)) ef[idx_sub] <- op$par
  ef
}

## templates (event + sampling times) from the wide data
build_long <- function(file){
  dat <- read_csv(file, show_col_types = FALSE)
  dat$id <- seq_len(nrow(dat))
  val_long <- dat %>%
    select(id, BI, `1stVaccine`, `2ndVaccine`, `3rdVaccine`, matches("^Sab\\d+Val$")) %>%
    pivot_longer(matches("^Sab\\d+Val$"), names_to="idx",
                 names_pattern="^Sab(\\d+)Val$", values_to="DV")
  date_long <- dat %>%
    select(id, matches("^Sab\\d+Date$")) %>%
    pivot_longer(matches("^Sab\\d+Date$"), names_to="idx",
                 names_pattern="^Sab(\\d+)Date$", values_to="meas_date")
  val_long %>% left_join(date_long, by=c("id","idx")) %>%
    mutate(meas_date=as.Date(meas_date),
           TIME=as.numeric(meas_date-as.Date(`1stVaccine`)),
           t2=as.numeric(as.Date(`2ndVaccine`)-as.Date(`1stVaccine`)),
           t3=as.numeric(as.Date(`3rdVaccine`)-as.Date(`1stVaccine`)),
           tau_raw=as.numeric(as.Date(BI)-as.Date(`1stVaccine`)),
           hasBI=ifelse(is.na(tau_raw),0,1),
           tau=ifelse(is.na(tau_raw),1e9,tau_raw),
           DV=as.numeric(DV)) %>%
    filter(!is.na(DV), DV>0, !is.na(TIME), !is.na(t2), !is.na(t3), TIME>=0) %>%
    select(id, TIME, t2, t3, tau, hasBI) %>% arrange(id, TIME)
}

## ---- run -----------------------------------------------------------
d_all  <- build_long(data_file)
fit    <- readRDS(fit_file)
pop    <- get_pop(fit)
ids    <- split(d_all, d_all$id)
sd_eta <- sqrt(pop$omega_diag)

rec <- list(); nfit <- 0L
for (rep in seq_len(N_REP)){
  for (ii in seq_along(ids)){
    d <- ids[[ii]]; n <- nrow(d); if (n < MIN_OBS) next
    cov <- c(t2=d$t2[1], t3=d$t3[1], tau=d$tau[1], hasBI=d$hasBI[1])
    idx_sub <- if (d$hasBI[1]==1) IDX_BI else IDX_NONBI
    times <- d$TIME
    eta_true <- rep(0, P)
    eta_true[idx_sub] <- rnorm(length(idx_sub), 0, sd_eta[idx_sub])
    lp_true <- pop$lpop + eta_true
    f <- solve_ab(lp_true, cov, times)
    y <- pmax(f * (1 + pop$sig * rnorm(length(f))), FLOOR)
    eta_hat <- fit_eta(idx_sub, cov, times, y, pop)
    nfit <- nfit + 1L
    rec[[length(rec)+1]] <- data.frame(
      rep = rep, id = d$id[1], hasBI = d$hasBI[1], n_obs = n,
      param = PAR[idx_sub],
      eta_true = eta_true[idx_sub], eta_hat = eta_hat[idx_sub],
      theta_true = exp(lp_true[idx_sub]),
      theta_est  = exp((pop$lpop + eta_hat)[idx_sub]))
  }
  message(sprintf("rep %d/%d | %d fits", rep, N_REP, nfit))
}
rec <- bind_rows(rec)
write.csv(rec, out_raw, row.names = FALSE)

perparam <- rec %>% group_by(param) %>% summarise(
  n           = n(),
  recov_cor   = cor(eta_true, eta_hat),              # 1 = perfect, ~0 = prior-dominated
  recov_slope = tryCatch(coef(lm(eta_hat ~ eta_true))[2], error=function(e) NA_real_),
  bias_eta    = mean(eta_hat - eta_true),
  sd_true_eta = sd(eta_true),
  sd_est_eta  = sd(eta_hat),
  sd_ratio    = sd(eta_hat) / sd(eta_true),          # <1 => variance shrunk away
  rmse_eta    = sqrt(mean((eta_hat - eta_true)^2)),
  cor_theta   = cor(theta_true, theta_est),
  .groups = "drop") %>% arrange(match(param, PAR))
write.csv(perparam, out_sum, row.names = FALSE)
print(perparam, digits = 3)
