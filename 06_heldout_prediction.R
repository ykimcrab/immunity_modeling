# Internal validation by held-out prediction (leave-one-time-point-out, LOTO).
# The population model is held FIXED; for each individual one observation is hidden,
# the subject's random effects are re-estimated by MAP from the remaining points,
# and the hidden point is predicted (with a population-only baseline for comparison).
# No SAEM refit. Run once per vaccination-group dataset.
#
# sample1 : empirical antibody data, wide format (same file used as input to 01)
# sample2 : fitted model object (.rds) from 01_fit_nlme_model.R

suppressMessages({ library(nlmixr2); library(rxode2); library(dplyr); library(tidyr); library(readr) })

data_file <- "sample1.csv"
fit_file  <- "sample2.rds"
out_pred  <- "loto_predictions.csv"     # <- per-point held-out predictions
out_byind <- "loto_byindividual.csv"    # <- per-individual held-out MAPE
out_sum   <- "loto_summary.csv"         # <- headline held-out error
out_phase <- "loto_phase_coverage.csv"  # <- error by peak-bracketing

ERR_MODEL <- "prop"
MIN_OBS   <- 3

PAR <- c("A0","k1","k2","k3","k4","k5","k6","k7","phi2","phi3","phit")
P   <- length(PAR); LNAMES <- paste0("l", PAR)
IDX_BI    <- seq_len(P)
IDX_NONBI <- which(!PAR %in% c("k6","k7","phit"))

## closed-form piecewise-exponential trajectory (exact; no ODE integrator)
solve_ab <- function(lp, cov, times){
  A0 <- exp(lp[["lA0"]])
  k1<-exp(lp[["lk1"]]);k2<-exp(lp[["lk2"]]);k3<-exp(lp[["lk3"]]);k4<-exp(lp[["lk4"]])
  k5<-exp(lp[["lk5"]]);k6<-exp(lp[["lk6"]]);k7<-exp(lp[["lk7"]])
  phi2<-exp(lp[["lphi2"]]);phi3<-exp(lp[["lphi3"]]);phit<-exp(lp[["lphit"]])
  t2<-cov[["t2"]];t3<-cov[["t3"]];tau<-cov[["tau"]];bi<-cov[["hasBI"]]
  rate_at <- function(t){
    if (t <= t2) return(k1); if (t <= t2+phi2) return(k2); if (t <= t3) return(-k3)
    if (t <= t3+phi3) return(k4); if (bi==0 || t<=tau) return(-k5)
    if (t <= tau+phit) return(k6); -k7 }
  bps <- c(t2, t2+phi2, t3, t3+phi3); if (bi==1) bps <- c(bps, tau, tau+phit)
  bps <- sort(unique(bps)); out <- numeric(length(times))
  for (i in seq_along(times)){
    th <- times[i]
    if (!is.finite(th)){ out[i] <- NA_real_; next }
    if (th <= 0){ out[i] <- A0; next }
    grid <- sort(unique(c(0, bps[bps>0 & bps<th], th))); A <- A0
    for (j in seq_len(length(grid)-1)){ a<-grid[j]; b<-grid[j+1]; A <- A*exp(rate_at((a+b)/2)*(b-a)) }
    out[i] <- A }
  out
}

get_pop <- function(fit){
  fe <- tryCatch(fixef(fit), error=function(e) NULL)
  if (is.null(fe)){ pf <- as.data.frame(fit$parFixedDf); fe <- setNames(pf$Estimate, rownames(pf)) }
  pick <- function(nm){ for (cand in c(paste0("log_",nm),nm,paste0("l",nm)))
    if (cand %in% names(fe)) return(as.numeric(fe[[cand]])); NA_real_ }
  lpop <- setNames(vapply(PAR, pick, numeric(1)), LNAMES)
  if (any(is.na(lpop))) stop("Could not map population params: ", paste(names(fe), collapse=", "))
  om <- fit$omega
  od <- setNames(diag(om), gsub("^eta_|^log_","", rownames(om)))
  omega_diag <- as.numeric(od[PAR]); names(omega_diag) <- PAR
  omega_diag[is.na(omega_diag) | omega_diag <= 0] <- 1e-6
  pf <- as.data.frame(fit$parFixedDf); sr <- grep("sigma|prop|add", rownames(pf), ignore.case=TRUE)
  sig <- if (length(sr)) as.numeric(pf$Estimate[sr[1]]) else 0.2
  list(lpop=lpop, omega_diag=omega_diag, sig=sig)
}

fit_eta <- function(idx_sub, cov, tobs, yobs, pop){
  lpop<-pop$lpop; od<-pop$omega_diag; sig<-pop$sig
  obj <- function(es){
    ef<-rep(0,P); ef[idx_sub]<-es; f<-solve_ab(lpop+ef, cov, tobs)
    if (any(!is.finite(f))) return(1e12); f<-pmax(f,1e-8)
    if (ERR_MODEL=="lnorm"){ r<-(log(yobs)-log(f))/sig; dd<-sum(0.5*r^2)+length(yobs)*log(sig)
    } else { sdv<-sig*f; r<-(yobs-f)/sdv; dd<-sum(0.5*r^2+log(sdv)) }
    dd + sum(0.5*es^2/od[idx_sub]) }
  op <- tryCatch(optim(rep(0,length(idx_sub)), obj, method="BFGS", control=list(maxit=200)), error=function(e) NULL)
  ef <- rep(0,P); if (!is.null(op)) ef[idx_sub] <- op$par; ef
}

## kinetic phase of a time point (for peak-bracketing stratification)
seg_of <- function(x, cov){
  t2<-cov[["t2"]];t3<-cov[["t3"]];tau<-cov[["tau"]];bi<-cov[["hasBI"]]
  if (x <= t2) return(1L); if (x <= t3) return(2L)
  if (bi == 0) return(3L); if (x <= tau) return(3L); 4L
}

build_long <- function(file){
  dat <- read_csv(file, show_col_types=FALSE); dat$id <- seq_len(nrow(dat))
  v <- dat %>% select(id, BI, `1stVaccine`, `2ndVaccine`, `3rdVaccine`, matches("^Sab\\d+Val$")) %>%
    pivot_longer(matches("^Sab\\d+Val$"), names_to="idx", names_pattern="^Sab(\\d+)Val$", values_to="DV")
  dte <- dat %>% select(id, matches("^Sab\\d+Date$")) %>%
    pivot_longer(matches("^Sab\\d+Date$"), names_to="idx", names_pattern="^Sab(\\d+)Date$", values_to="meas_date")
  v %>% left_join(dte, by=c("id","idx")) %>%
    mutate(meas_date=as.Date(meas_date),
           TIME=as.numeric(meas_date-as.Date(`1stVaccine`)),
           t2=as.numeric(as.Date(`2ndVaccine`)-as.Date(`1stVaccine`)),
           t3=as.numeric(as.Date(`3rdVaccine`)-as.Date(`1stVaccine`)),
           tau_raw=as.numeric(as.Date(BI)-as.Date(`1stVaccine`)),
           hasBI=ifelse(is.na(tau_raw),0,1), tau=ifelse(is.na(tau_raw),1e9,tau_raw),
           DV=as.numeric(DV)) %>%
    filter(!is.na(DV), DV>0, !is.na(TIME), !is.na(t2), !is.na(t3), TIME>=0) %>%
    select(id, TIME, DV, t2, t3, tau, hasBI) %>% arrange(id, TIME)
}

predict_holdout <- function(d, ret_idx, ho_idx, pop){
  cov <- c(t2=d$t2[1], t3=d$t3[1], tau=d$tau[1], hasBI=d$hasBI[1])
  idx_sub <- if (d$hasBI[1]==1) IDX_BI else IDX_NONBI
  ef <- fit_eta(idx_sub, cov, d$TIME[ret_idx], d$DV[ret_idx], pop)
  lp <- pop$lpop + ef; th <- d$TIME[ho_idx]; rt <- d$TIME[ret_idx]
  data.frame(id=d$id[1], hasBI=d$hasBI[1], n_obs=nrow(d), n_retained=length(ret_idx),
             holdout_time=th, holdout_DV=d$DV[ho_idx],
             pred_ind=solve_ab(lp, cov, th), pred_pop=solve_ab(pop$lpop, cov, th),
             seg=vapply(th, seg_of, integer(1), cov=cov),
             bracketed_seg=vapply(th, function(x){
               s<-seg_of(x,cov); sr<-vapply(rt, seg_of, integer(1), cov=cov)
               any(rt<x & sr==s) & any(rt>x & sr==s) }, logical(1)))
}

add_errors <- function(df){
  df %>% mutate(.lo=log10(pmax(holdout_DV,1e-8)), .li=log10(pmax(pred_ind,1e-8)),
                .lp=log10(pmax(pred_pop,1e-8)),
                log10_err=abs(.lo-.li), mape_log10=abs(.lo-.li)/abs(.lo),
                mape_log10_pop=abs(.lo-.lp)/abs(.lo)) %>%
    select(-.lo, -.li, -.lp)
}

## ---- run: LOTO -----------------------------------------------------
d_all <- build_long(data_file)
pop   <- get_pop(readRDS(fit_file))
ids   <- split(d_all, d_all$id)

loto <- list()
for (ii in seq_along(ids)){
  d <- ids[[ii]]; n <- nrow(d); if (n < MIN_OBS) next
  for (j in seq_len(n))
    loto[[length(loto)+1]] <- predict_holdout(d, setdiff(seq_len(n), j), j, pop)
}
loto <- add_errors(bind_rows(loto))
write.csv(loto, out_pred, row.names = FALSE)

loto_byind <- loto %>% group_by(id, hasBI) %>%
  summarise(heldout_MAPE_pct = 100*mean(mape_log10, na.rm=TRUE), n_points = n(), .groups="drop")
write.csv(loto_byind, out_byind, row.names = FALSE)

loto_summary <- cbind(
  loto %>% summarise(n_points=n(), n_individuals=n_distinct(id),
                     heldout_MAPE_mean_pct=100*mean(mape_log10, na.rm=TRUE),
                     pop_baseline_MAPE_mean_pct=100*mean(mape_log10_pop, na.rm=TRUE)),
  loto_byind %>% summarise(indMAPE_mean_pct=mean(heldout_MAPE_pct),
                           indMAPE_median_pct=median(heldout_MAPE_pct),
                           pct_ind_MAPE_lt10=100*mean(heldout_MAPE_pct<10),
                           pct_ind_MAPE_lt15=100*mean(heldout_MAPE_pct<15)))
write.csv(loto_summary, out_sum, row.names = FALSE)

phase_strat <- loto %>% group_by(bracketed_seg, seg) %>%
  summarise(n=n(), median_abs_log10err=median(log10_err, na.rm=TRUE),
            mean_MAPE_pct=100*mean(mape_log10, na.rm=TRUE), .groups="drop")
write.csv(phase_strat, out_phase, row.names = FALSE)
print(loto_summary)
