library(readxl)
library(car)
library(MASS)
library(lme4)
library(ggplot2)
library(multcomp)
setwd('C:/Users/jbai5/OneDrive/First year project/Experiment_Onion/')
DATA <- read_excel('Onion_fit_pertTrials_LOSO12_RMSE_v_2018-8-9-15-33.xlsx')
modelList <- c("expansion", "expansion2")
DATA <- DATA[DATA$Model %in% modelList, ]
DATA <- DATA[DATA$w %in% c(1, 0.6), ]


qqp(DATA$RMSE_v, "norm")

# lnorm means lognormal
qqp(DATA$RMSE_v, "lnorm")

# qqp requires estimates of the parameters of the negative binomial, Poisson
# and gamma distributions. You can generate estimates using the fitdistr
# function. Save the output and extract the estimates of each parameter as I
# have shown below.
nbinom <- fitdistr(DATA$RMSE_v, "Negative Binomial")
qqp(DATA$RMSE_v, "nbinom", size = nbinom$estimate[[1]], mu = nbinom$estimate[[2]])

poisson <- fitdistr(DATA$RMSE_v, "Poisson")
qqp(DATA$RMSE_v, "pois", poisson$estimate)

gamma <- fitdistr(DATA$RMSE_v, "gamma")
qqp(DATA$RMSE_v, "gamma", shape = gamma$estimate[[1]], rate = gamma$estimate[[2]])

lmm0 <- lmer(RMSE_v ~ Model * w * dv + (1 | subject), 
             data = DATA, REML = FALSE)
lmm1 <- lmer(RMSE_v ~ Model + w + dv + (1 | subject), 
             data = DATA, REML = FALSE)
lmm2 <- lmer(RMSE_v ~ Model * w + dv + (1 | subject),
             data = DATA, REML = FALSE)
lmm3 <- lmer(RMSE_v ~ Model * dv + w + (1 | subject),
             data = DATA, REML = FALSE) # Model:dv is sig
lmm4 <- lmer(RMSE_v ~ Model + w * dv + (1 | subject),
             data = DATA, REML = FALSE) # w:dv is sig

lmm5 <- lmer(RMSE_v ~ Model + w + dv + Model:dv + dv:w +
               (1 | subject), data = DATA, REML = FALSE)

lmm6 <- lmer(RMSE_v ~ Model + w + dv + Model:dv + dv:w + Model:w:dv +
               (1 | subject), data = DATA, REML = FALSE)

lmm7 <- lmer(RMSE_v ~ Model + dv + Model:dv + dv:w +
               (1 | subject), data = DATA, REML = FALSE)

lmm8 <- lmer(RMSE_v ~ dv + Model:dv + dv:w +
               (1 | subject), data = DATA, REML = FALSE)

lmm9 <- lmer(RMSE_v ~ Model + Model:dv + dv:w +
               (1 | subject), data = DATA, REML = FALSE)

lmm10 <- lmer(RMSE_v ~ Model:dv + dv:w +
                (1 | subject), data = DATA, REML = FALSE)



summary(lmm9)


Anova(lmm0, type = 3)


anova(lmm1, lmm2)
anova(lmm1, lmm3) # BIC lmm3 < lmm1
anova(lmm1, lmm4) # BIC lmm4 < lmm1
anova(lmm3, lmm5) # BIC lmm5 < lmm3
anova(lmm4, lmm5) # BIC lmm5 < lmm4
anova(lmm5, lmm6) # BIC lmm5 < lmm6
anova(lmm7, lmm5) # BIC lmm7 = lmm5
anova(lmm8, lmm7) # BIC lmm7 < lmm8
anova(lmm9, lmm7) # BIC lmm9 = lmm7
anova(lmm9, lmm10) # BIC lmm9 < lmm10

data <- DATA[c("Model", "w", "dv", "RMSE_v")]
data$w <- as.character(data$w)
data$dv <- as.character(data$dv)
# Main effect
ggplot(data) + aes(Model, RMSE_v) + geom_boxplot()
# Interactions
ggplot(data) + aes(Model, RMSE_v) + geom_boxplot() + facet_wrap(~w)

summary(glht(lmm0, linfct = mcp(Model = "Tukey")), test = adjusted("holm"))

