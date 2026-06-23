# Sensitivity of post-infection estimates to uncertainty in the BI date (t_b).
# Population (fixed effects, Omega, sigma) is held FIXED at the original fit; for
# each BI subject the recorded BI date is shifted by delta days, the subject's
# random effects are re-estimated by MAP, and the post-BI quantities are recomputed.
# No SAEM refit. Run once per vaccination-group dataset.
#
# sample1 : empirical antibody data, wide format (same file used as input to 01)
# sample2 : fitted model object (.rds) produced by 01_fit_nlme_model.R

suppressMessages({ library(dplyr); library(tidyr); library(readr) })

data_file <- "sample1.csv"                   # <- empirical antibody data (wide)
fit_file  <- "sample2.rds"                   # <- fitted model object (from 01)
out_ind   <- "tb_perturbation_byindividual.csv"
out_sum   <- "tb_perturbation_summary.csv"

ERR_MODEL <- "prop"
DELTAS    <- c(-5, -3, 0, 3, 5)              # days to shift the BI date
MIN_OBS   <- 3

PAR    <- c("A0","k1","k2","k3","k4","k5","k6","k7","phi2","phi3","phit")
P      <- length(PAR)
LNAMES <- paste0("l", PAR)
IDX_BI <- seq_len(P)

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

## population quantities from the fit: fixed effects, diagonal Omega, sigma
get_pop <- function(fit){
  fe <- tryCatch(nlmixr2::fixef(fit), error = function(e) NULL)
  if (is.null(fe)) { pf <- as.data.frame(fit$parFixedDf); fe <- setNames(pf$Estimate, rownames(pf)) }
  pick <- function(nm){
    for (cand in c(paste0("log_", nm), nm, paste0("l", nm)))
      if (cand %in% names(fe)) return(as.numeric(fe[[cand]]))
    NA_real_ }
  lpop <- setNames(vapply(PAR, pick, numeric(1)), LNAMES)
  if (any(is.na(lpop))) stop("Could not map population log-parameters: ", paste(names(fe), collapse=", "))
  om <- fit$omega
  od <- setNames(diag(om), gsub("^eta_|^log_", "", rownames(om)))
  omega_diag <- as.numeric(od[PAR]); names(omega_diag) <- PAR
  omega_diag[is.na(omega_diag) | omega_diag <= 0] <- 1e-6
  pf <- as.data.frame(fit$parFixedDf)
  sig_row <- grep("sigma|prop|add", rownames(pf), ignore.case = TRUE)
  sig <- if (length(sig_row)) as.numeric(pf$Estimate[sig_row[1]]) else 0.2
  list(lpop = lpop, omega_diag = omega_diag, sig = sig)
}

## MAP estimate of one subject's random effects (population fixed)
fit_eta <- function(idx_sub, cov, tobs, yobs, pop){
  lpop <- pop$lpop; od <- pop$omega_diag; sig <- pop$sig
  obj <- function(es){
    ef <- rep(0, P); ef[idx_sub] <- es
    f  <- solve_ab(lpop + ef, cov, tobs)
    if (any(!is.finite(f))) return(1e12)
    f  <- pmax(f, 1e-8)
    if (ERR_MODEL == "lnorm"){
      r <- (log(yobs) - log(f)) / sig; dd <- sum(0.5 * r^2) + length(yobs) * log(sig)
    } else {
      sdv <- sig * f; r <- (yobs - f) / sdv; dd <- sum(0.5 * r^2 + log(sdv))
    }
    dd + sum(0.5 * es^2 / od[idx_sub])
  }
  op <- tryCatch(optim(rep(0, length(idx_sub)), obj, method = "BFGS",
                       control = list(maxit = 200)), error = function(e) NULL)
  ef <- rep(0, P); if (!is.null(op)) ef[idx_sub] <- op$par
  ef
}

## templates + observations from the wide data
build_long <- function(file){
  dat <- read_csv(file, show_col_types = FALSE)
  dat$id <- seq_len(nrow(dat))
  val_long <- dat %>%
    select(id, BI, `1stVaccine`, `2ndVaccine`, `3rdVaccine`, matches("^Sab\\d+Val$")) %>%
    pivot_longer(matches("^Sab\\d+Val$"), names_to="idx", names_pattern="^Sab(\\d+)Val$", values_to="DV")
  date_long <- dat %>%
    select(id, matches("^Sab\\d+Date$")) %>%
    pivot_longer(matches("^Sab\\d+Date$"), names_to="idx", names_pattern="^Sab(\\d+)Date$", values_to="meas_date")
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
    select(id, TIME, DV, t2, t3, tau, hasBI) %>% arrange(id, TIME)
}

## per-subject estimate under a shifted BI date
estimate_one <- function(d, delta, pop){
  cov <- c(t2 = d$t2[1], t3 = d$t3[1], tau = d$tau[1] + delta, hasBI = 1)
  if (cov[["tau"]] <= cov[["t3"]] + 1) return(NULL)   # BI must stay after the booster window
  ef <- fit_eta(IDX_BI, cov, d$TIME, d$DV, pop)
  lp <- pop$lpop + ef
  k6 <- exp(lp[["lk6"]]); k7 <- exp(lp[["lk7"]]); taub <- exp(lp[["lphit"]])
  peak_lvl <- solve_ab(lp, cov, cov[["tau"]] + taub)
  pre_lvl  <- solve_ab(lp, cov, cov[["tau"]])
  pred     <- solve_ab(lp, cov, d$TIME)
  ok <- is.finite(pred) & pred > 0 & is.finite(d$DV) & d$DV > 0
  mape <- if (any(ok)) 100*mean(abs(log10(d$DV[ok]) - log10(pred[ok])) / abs(log10(d$DV[ok]))) else NA_real_
  data.frame(id = d$id[1], delta = delta, k6 = k6, k7 = k7, tau_b_days = taub,
             peak_log10 = log10(pmax(peak_lvl, 1e-8)),
             peak_increase_log10 = log10(pmax(peak_lvl,1e-8)) - log10(pmax(pre_lvl,1e-8)),
             insample_MAPE_pct = mape)
}

## ---- run -----------------------------------------------------------
d_all  <- build_long(data_file)
pop    <- get_pop(readRDS(fit_file))
ids    <- split(d_all, d_all$id)
bi_ids <- Filter(function(d) d$hasBI[1] == 1 && nrow(d) >= MIN_OBS, ids)
message(sprintf("%d BI subjects x %d shifts", length(bi_ids), length(DELTAS)))

rows <- list()
for (ii in seq_along(bi_ids)){
  d <- bi_ids[[ii]]
  for (dl in DELTAS){
    r <- estimate_one(d, dl, pop)
    if (!is.null(r)) rows[[length(rows)+1]] <- r
  }
}
ind <- bind_rows(rows)
write.csv(ind, out_ind, row.names = FALSE)

base <- ind %>% filter(delta == 0) %>%
  select(id, k6_0 = k6, taub_0 = tau_b_days, peak_0 = peak_log10, mape_0 = insample_MAPE_pct)
summ <- ind %>% left_join(base, by = "id") %>%
  mutate(d_taub_days = tau_b_days - taub_0, pct_k6 = 100*(k6 - k6_0)/k6_0,
         d_peak_log10 = peak_log10 - peak_0, d_mape_pct = insample_MAPE_pct - mape_0) %>%
  group_by(delta) %>%
  summarise(n = n(),
            med_tau_b_days   = median(tau_b_days, na.rm = TRUE),
            med_d_tau_b_days = median(d_taub_days, na.rm = TRUE),
            med_k6           = median(k6, na.rm = TRUE),
            med_pct_k6       = median(pct_k6, na.rm = TRUE),
            med_peak_log10   = median(peak_log10, na.rm = TRUE),
            med_d_peak_log10 = median(d_peak_log10, na.rm = TRUE),
            mean_MAPE_pct    = mean(insample_MAPE_pct, na.rm = TRUE),
            .groups = "drop")
write.csv(summ, out_sum, row.names = FALSE)
print(summ)
