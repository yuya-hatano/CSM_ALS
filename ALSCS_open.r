library(survival)
library(survminer)
library(tidyverse)
library(ggplot2)
library(data.table)
library(coxphf)

#als_prepared.csv:CSV file of the source data
#als_prepared20260719.csv:CSV file2 of the source data

#########################
#Mantel–Byar test
#########################
whas100 <- fread("als_prepared.csv")
whas100 <- whas100[whas100$upper==1,]
whas100[, id := .I]
whas100$stat <- whas100$event
surgery <- data.frame(id=whas100$id, surgery=whas100$t_surg)
surgery <- na.omit(surgery)
surgery$surgery_st <- 1
nons <- data.frame(id=whas100$id)
nons$surgery_st <- 0
nons$surgery <- 0
surgery <- rbind(surgery, nons)
whas100[, c("t_surg") := NULL]
data1 <- tmerge(whas100, whas100, id=id, tstop=t_onset, death  = event(t_onset, stat))
data2 <- tmerge(data1, surgery, id=id, cov=tdc(surgery, surgery_st))
fit_mb <- coxph(
  Surv(tstart, tstop, death) ~ cov,
  data = data2
)
summary(fit_mb)

########################
#Time-dependent Cox proportional hazards model with Firth penalization
########################
whas100 <- fread("als_prepared.csv")
whas100 <- whas100[whas100$upper==1,]
whas100[, id := .I]
whas100$stat <- whas100$event
surgery <- data.frame(id=whas100$id, surgery=whas100$t_surg)
surgery <- na.omit(surgery)
surgery$surgery_st <- 1
nons <- data.frame(id=whas100$id)
nons$surgery_st <- 0
nons$surgery <- 0
surgery <- rbind(surgery, nons)
whas100[, c("t_surg") := NULL]
data1 <- tmerge(whas100, whas100, id=id, tstop=t_onset, death  = event(t_onset, stat))
data2 <- tmerge(data1, surgery, id=id, cov=tdc(surgery, surgery_st))
fit_firth <- coxphf(
  Surv(tstart, tstop, death) ~ cov + female + age,
  data = data2,
  firth = TRUE,
  pl = TRUE
)
summary(fit_firth)
#Regression coefficient
coef(fit_firth)["cov"]
#SE
sqrt(fit_firth$var[1,1])


############################
#Time-dependent Cox proportional hazards model (model-based variance)
############################
whas100 <- fread("als_prepared.csv")
whas100 <- whas100[whas100$upper==1,]
whas100[, id := .I]
whas100$stat <- whas100$event
surgery <- data.frame(id=whas100$id, surgery=whas100$t_surg)
surgery <- na.omit(surgery)
surgery$surgery_st <- 1
nons <- data.frame(id=whas100$id)
nons$surgery_st <- 0
nons$surgery <- 0
surgery <- rbind(surgery, nons)
whas100[, c("t_surg") := NULL]
data1 <- tmerge(whas100, whas100, id=id, tstop=t_onset, death  = event(t_onset, stat))
data2 <- tmerge(data1, surgery, id=id, cov=tdc(surgery, surgery_st))
fit <- coxph(Surv(tstart, tstop, death) ~ cov + female + age, data = data2)

summary(fit)
#Proportional hazards assumption
cox.zph(fit)
#Regression coefficient
coef(fit)["cov"]
#SE
sqrt(fit$var[1,1])

############################
#Time-dependent Cox proportional hazards model with cluster-robust variance
############################
whas100 <- fread("als_prepared.csv")
whas100 <- whas100[whas100$upper==1,]
whas100[, id := .I]
whas100$stat <- whas100$event
surgery <- data.frame(id=whas100$id, surgery=whas100$t_surg)
surgery <- na.omit(surgery)
surgery$surgery_st <- 1
nons <- data.frame(id=whas100$id)
nons$surgery_st <- 0
nons$surgery <- 0
surgery <- rbind(surgery, nons)
whas100[, c("t_surg") := NULL]
data1 <- tmerge(whas100, whas100, id=id, tstop=t_onset, death  = event(t_onset, stat))
data2 <- tmerge(data1, surgery, id=id, cov=tdc(surgery, surgery_st))
fit <- coxph(Surv(tstart, tstop, death) ~ cov + female + age, data = data2,
                        cluster = id,
                        robust = TRUE)

summary(fit)
#Proportional hazards assumption
cox.zph(fit)
#Regression coefficient
coef(fit)["cov"]
#SE
sqrt(fit$var[1,1])
##################
#Swimmer plot
##################
data <- fread("als_prepared20260719.csv")
data <- data[data$upper==1,]
data <- data[data$surg==1,]
data <- data[order(data$t_surg), ]
data$patient <- 1:nrow(data)

data2 <- data
data3 <- data
data$event <- "operation"
data2$event <- "death"
data2$event[data2$TPPV==1] <- "TPPV"
data3$event <- "ALS diagnosis"
data$point <- data$t_surg
data2$point <- data2$t_onset
data3$point <- data3$delay
data_all <- rbind(data, data2)
data_all <- rbind(data_all, data3)

survp <- ggplot(data = data, aes(x = 0, xend = t_onset, y = patient, yend = patient)) +
  geom_segment(size = 1, color = "black") +
  geom_point(data= data_all, aes(x = point, y = patient, shape = event, color = event), size = 3) +
  scale_x_continuous(breaks = seq(0, 70, by = 10), limits = c(0, 70)) + scale_y_reverse(breaks = 1:10)+
  labs(x = "months from onset",
       y = "surgical patient ID") +
  theme_minimal() +
  theme(panel.grid.minor = element_blank(),
        axis.text.y = element_text(size = 10))
survp

pdf("swimmer_plot.pdf")
print(survp, newpage = FALSE)
dev.off()


##############
#Simon–Makuch curve
##############
whas100 <- fread("als_prepared.csv")
whas100 <- whas100[whas100$upper==1,]
whas100[, id := .I]
whas100$stat <- whas100$event
surgery <- data.frame(id=whas100$id, surgery=whas100$t_surg)
surgery <- na.omit(surgery)
surgery$surgery_st <- 1
nons <- data.frame(id=whas100$id)
nons$surgery_st <- 0
nons$surgery <- 0
surgery <- rbind(surgery, nons)
whas100[, c("t_surg") := NULL]
data1 <- tmerge(whas100, whas100, id=id, tstop=t_onset, death  = event(t_onset, stat))
data2 <- tmerge(data1, surgery, id=id, cov=tdc(surgery, surgery_st))
data2$cov_plot <- factor(
  data2$cov,
  levels = c(0, 1),
  labels = c("No surgery / Preoperative", "Postoperative")
)

fit_sm <- survfit(
  Surv(tstart, tstop, death) ~ cov_plot,
  data = data2,
  id = id,
  conf.type = "log"
)

plot(
  fit_sm,
  lty = c(1, 2),
  col = c("blue", "red"),
  lwd = 2,
  xlab = "Time from diagnosis (months)",
  ylab = "Probability of survival free from death or TPPV",
  mark.time = TRUE,
  conf.int = FALSE,
  xaxs = "i",
  yaxs = "i"
)

legend(
  "topright",
  legend = c("No surgery / Preoperative", "Postoperative"),
  lty = c(1, 2),
  col = c("blue", "red"),
  lwd = 2,
  bty = "n"
)


##############
#Kaplan Meyer
##############
whas100 <- fread("ALSC20251010_josijoho.csv")
whas100 <- whas100[whas100$josi1==1,]


KM <- whas100
KM$cov_plot <- factor(
  KM$als0c1,
  levels = c(0, 1),
  labels = c("No surgery group", "Surgery group")
)

fit_sm <- survfit(
  Surv(survival, stat) ~ cov_plot,
  data = KM
)

plot(
  fit_sm,
  lty = c(1, 2),
  lwd = 2,
  xlab = "Time from diagnosis (months)",
  ylab = "Probability of survival free from death or permanent ventilation",
  mark.time = TRUE,
  conf.int = FALSE,
  xaxs = "i",
  yaxs = "i"
)

legend(
  "topright",
  legend = c("No surgery group", "Surgery group"),
  lty = c(1, 2),
  lwd = 2,
  bty = "n"
)

##############
#forest_plot
##############
#p=0.00704、HR=2.01449,(95%CI = 1.2105 - 3.352)）
forest_data <- data.frame(
  model = c(
    "Unadjusted Cox model (Model based variance)",
    "Age- and sex-adjusted Cox model (Model based variance)",
    "Age- and sex-adjusted Cox model (cluster-robust SEs)",
    "Adjusted Firth Cox model",
    "IPCW-adjusted Cox model"
  ),
  HR = c(2.01, 1.85, 1.85, 1.90, 2.01),
  lower = c(1.00, 0.91, 1.11, 0.91, 1.21),
  upper = c(4.06, 3.73, 3.07, 3.66, 3.35)
)

forest_data$model <- factor(
  forest_data$model,
  levels = rev(forest_data$model)
)

ggplot(
  forest_data,
  aes(x = HR, y = model)
) +
  geom_vline(
    xintercept = 1,
    linetype = 2
  ) +
  geom_errorbarh(
    aes(xmin = lower, xmax = upper),
    height = 0.15
  ) +
  geom_point(size = 3) +
  scale_x_log10(
    breaks = c(0.5, 1, 2, 4, 8)
  ) +
  labs(
    x = "Hazard ratio (log scale)",
    y = NULL
  ) +
  theme_classic(base_size = 13)

###################################
#Time-dependent Cox proportional hazards model with
# inverse probability of censoring weighting (IPCW)
###################################
library(ipw)
whas100 <- fread("als_prepared.csv")
whas100 <- whas100[whas100$upper==1,]
whas100[, id := .I]
whas100$stat <- whas100$event
surgery <- data.frame(id=whas100$id, surgery=whas100$t_surg)
surgery <- na.omit(surgery)
surgery$surgery_st <- 1
nons <- data.frame(id=whas100$id)
nons$surgery_st <- 0
nons$surgery <- 0
surgery <- rbind(surgery, nons)
whas100[, c("t_surg") := NULL]
data1 <- tmerge(whas100, whas100, id=id, tstop=t_onset, death  = event(t_onset, stat))
data2 <- tmerge(data1, surgery, id=id, cov=tdc(surgery, surgery_st))
setDT(data2)
setorder(data2, id, tstart, tstop)

data2[, final_interval := seq_len(.N) == .N, by = id]

data2[, censor_event :=
        as.integer(final_interval & death == 0L)]
ipcw_object <- ipwtm(
  exposure    = censor_event,
  family      = "survival",
  numerator   = NULL,                    # Unstabilized weights
  denominator = ~ cov + female + age,
  id          = id,
  tstart      = tstart,
  timevar     = tstop,
  type        = "cens",
  trunc       = NULL,                    # No weight truncation
  data        = data2
)
data2[, w_ipcw := ipcw_object$ipw.weights]
fit <- coxph(Surv(tstart, tstop, death) ~ cov + female + age, data = data2,
                        weights = w_ipcw,
                        cluster = id,
                        robust = TRUE,
  ties    = "efron")
# No weight truncation
summary(fit)

coef(fit)["cov"]
sqrt(fit$var[1,1])
# Basic distribution of IPCW weights
summary(data2$w_ipcw)


##99th percentile truncation
q99 <- quantile(data2$w_ipcw, 0.99, na.rm = TRUE)

data2[, w_ipcw_trunc99 :=
        pmin(w_ipcw, q99)]

fit_ipcw_trunc99 <- coxph(
  Surv(tstart, tstop, death) ~ cov + female + age,
  data = data2,
  weights = w_ipcw_trunc99,
  cluster = id,
  robust = TRUE,
  ties = "efron"
)

summary(fit_ipcw_trunc99)
coef(fit_ipcw_trunc99)["cov"]
sqrt(fit_ipcw_trunc99$var[1,1])
