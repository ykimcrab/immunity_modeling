# tau_b fixed refit: compare the full model (tau_b has a subject-level random
# effect) with a reduced model in which tau_b is FIXED at the full-model
# population value (random effect removed). Reports OFV/AIC/BIC and in-sample fit.
# One SAEM refit (the reduced model). Run once per vaccination-group dataset.
#
# sample1 : empirical antibody data, wide format (same file used as input to 01)
# sample2 : full fitted model object (.rds) from 01_fit_nlme_model.R
#
# IMPORTANT: tau_b is fixed via fix(), which requires a literal value in the ini
# block. On load, this script prints the full-model log_phit estimate; set the
# fix() value in the reduced model below to that value for THIS group, then run.

library(nlmixr2); library(rxode2); library(dplyr); library(tidyr); library(readr)

data_file <- "sample1.csv"
fit_file  <- "sample2.rds"
out_cmp   <- "taub_fixed_modelcompare.csv"
out_red   <- "taub_fixed_reduced_fit.rds"

## ---- data loader (same transform as 01) ----------------------------
load_long <- function(file){
  dat <- read_csv(file, show_col_types = FALSE); dat$id <- seq_len(nrow(dat))
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
    filter(!is.na(DV), !is.na(TIME), !is.na(t2), !is.na(t3), TIME>=0) %>%
    select(id, TIME, DV, t2, t3, tau, hasBI)
}

## read the full-model log_phit estimate (the value to fix tau_b at)
get_logphit_est <- function(fit){
  pf <- as.data.frame(fit$parFixedDf); i <- which(rownames(pf) == "log_phit")
  if (length(i) != 1) return(NA_real_)
  as.numeric(pf[i, intersect(c("Estimate","Est"), colnames(pf))[1]])
}

fit_full <- readRDS(fit_file)
message(sprintf("Full-model log_phit = %.8f  (tau_b = %.3f days)",
                get_logphit_est(fit_full), exp(get_logphit_est(fit_full))))
## ^^^ set LOG_PHIT_FIX below to this printed log_phit value, then continue.

## ---- REDUCED model: tau_b fixed, eta_phit removed ------------------
## EDIT the fix() value to the full-model log_phit printed above (this group).
antibody_model_taub_fixed <- function() {
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
    log_phit <- fix(2.21141441731798)   # <- REPLACE with the full-model log_phit (this group)

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
    # eta_phit removed -> tau_b common to all subjects

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
    phit <- exp(log_phit)
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

## ---- helpers: model-fit criteria + in-sample metrics ---------------
sc <- function(x) if (is.null(x) || length(x) == 0) NA_real_ else suppressWarnings(as.numeric(x[[1]])[1])

# force a finite objective (SAEM fits may need setOfv with FOCEi on this stiff model)
force_ofv <- function(fit){
  fin <- function(f){ od <- tryCatch(as.data.frame(f$objDf), error=function(e) NULL)
    if (is.null(od) || nrow(od) < 1) return(FALSE)
    j <- grep("objf|objective", colnames(od), ignore.case=TRUE)[1]
    !is.na(j) && is.finite(sc(od[1, j])) }
  if (fin(fit)) return(fit)
  for (m in c("focei","foce","laplace")){
    f2 <- tryCatch(suppressWarnings(setOfv(fit, m)), error=function(e) NULL)
    if (!is.null(f2) && fin(f2)) return(f2)
  }
  fit
}

get_ic <- function(fit, label){
  od <- tryCatch(as.data.frame(fit$objDf), error=function(e) NULL)
  OFV <- AICv <- BICv <- LL <- NA_real_
  if (!is.null(od) && nrow(od) >= 1){
    cn <- colnames(od)
    pick <- function(pat){ i <- grep(pat, cn, ignore.case=TRUE)[1]; if (is.na(i)) NA_real_ else sc(od[nrow(od), i]) }
    OFV <- pick("objf|objective"); AICv <- pick("^AIC$"); BICv <- pick("^BIC$"); LL <- pick("ikelihood|logLik")
  }
  if (is.na(AICv)) AICv <- sc(tryCatch(AIC(fit), error=function(e) NA))
  if (is.na(BICv)) BICv <- sc(tryCatch(BIC(fit), error=function(e) NA))
  ll <- tryCatch(logLik(fit), error=function(e) NULL)
  if (is.na(LL)) LL <- sc(tryCatch(as.numeric(ll), error=function(e) NA))
  data.frame(model=label, n_par=sc(tryCatch(attr(ll,"df"), error=function(e) NA)),
             logLik=LL, OFV=OFV, AIC=AICv, BIC=BICv, stringsAsFactors=FALSE)
}

# in-sample fit from stored DV/IPRED (log10-normalized MAPE; always finite)
insample_metrics <- function(fit, label){
  nm <- tryCatch(names(fit), error=function(e) NULL)
  getcol <- function(w){ i <- match(tolower(w), tolower(nm)); if (is.na(i)) NULL else fit[[nm[i]]] }
  DV <- suppressWarnings(as.numeric(getcol("DV"))); IP <- suppressWarnings(as.numeric(getcol("IPRED")))
  ID <- getcol("ID"); ID <- if (is.null(ID)) rep("all", length(DV)) else as.character(ID)
  ok <- is.finite(DV) & DV>0 & is.finite(IP) & IP>0; DV<-DV[ok]; IP<-IP[ok]; ID<-ID[ok]
  per_id <- tapply(abs(log10(DV)-log10(IP))/abs(log10(DV)), ID, function(z) mean(z, na.rm=TRUE)*100)
  data.frame(model=label, n_obs=length(DV),
             insample_MAPE_pct_mean=mean(per_id, na.rm=TRUE),
             insample_MAPE_pct_median=median(per_id, na.rm=TRUE))
}

## ---- run: refit reduced, compare to full ---------------------------
dat_long <- load_long(data_file)
fit_red  <- nlmixr2(antibody_model_taub_fixed, dat_long, est="saem",
                    control = saemControl(nBurn=500, nEm=4500, print=50))
saveRDS(fit_red, out_red)

fit_full <- force_ofv(fit_full); fit_red <- force_ofv(fit_red)
ic  <- bind_rows(get_ic(fit_full,"full_taub_random"), get_ic(fit_red,"reduced_taub_fixed"))
ism <- bind_rows(insample_metrics(fit_full,"full_taub_random"), insample_metrics(fit_red,"reduced_taub_fixed"))
cmp <- left_join(ic, ism, by="model")
cmp$delta_BIC_vs_full <- c(NA, cmp$BIC[2]-cmp$BIC[1])
cmp$delta_OFV_vs_full <- c(NA, cmp$OFV[2]-cmp$OFV[1])
write.csv(cmp, out_cmp, row.names = FALSE)
print(cmp)
