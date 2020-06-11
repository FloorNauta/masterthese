#### Master thesis Floor Nauta
#code for analysis of hypotheses;
#plot per hypothesis

#loading packages
library(dplyr)
library(modelr)
library(parSim)
library(tidyr)
library(bootnet)
library(qgraph)
library(igraph)


load('results_total.rdata')
load('results_releveled.rdata')


# add noise to specificity & MAE for calculation Rsquared ----
##SPECIFICITY
data_spec <- results_releveled
set.seed(123)
data_spec$specificity <- data_spec$specificity + rnorm(50000, 0, 0.2)
#values between 0 and 1
data_spec$specificity[data_spec$specificity > 1] <- 1
data_spec$specificity[data_spec$specificity < 0] <- 0


##MAE
data_MAE <- results_releveled
set.seed(123)
data_MAE$MAE <- data_MAE$MAE + rnorm(50000, 0, 0.4)
data_MAE$MAE[data_MAE$MAE > 1] <- 1
data_MAE$MAE[data_MAE$MAE < 0] <- 0




# long version results for plots ----
long_results <- results_total %>%
  pivot_longer(c('sensitivity', 'specificity', 'correlation', 'MAE'), 
             names_to = "outcome_measure", values_to = "outcome_value")
long_results$outcome_value[long_results$outcome_value < 0] <- 0

long_results <- transform(long_results, 
                             outcome_measure = factor(outcome_measure, levels = c('sensitivity', 'specificity', 'correlation', 'MAE'), 
                                                      labels = c("Sensitivity", "Specificity", "Correlation", "MAE"))
                             ,
                             conn_level = factor(conn_level, levels = c(1, 2, 3, 4),
                                                 labels = c("Lowest", "Low", "Medium", "High")),
                             network_type = factor(network_type, levels = c("random", "scalefree", "smallworld"),
                                                   labels = c("Random", "Scale free", "Small world")))


#ANALYSES hypothesis 1  ----
#hyp 1 = main effect sample size (n)

##SENSITIVITY
sen_hyp1 <- glm(sensitivity ~ n, family = binomial(link = "logit"), data = results_releveled)
summary(sen_hyp1)

#McFadden’s Pseudo-Rsquared   #https://stats.idre.ucla.edu/other/mult-pkg/faq/general/faq-what-are-pseudo-r-squareds/
nullmod_sen <- glm(sensitivity ~ 1, family = binomial(link = "logit"), data = results_releveled)
sen_pR2_hyp1 <- 1-logLik(sen_hyp1)/logLik(nullmod_sen) #https://thestatsgeek.com/2014/02/08/r-squared-in-logistic-regression/
sen_pR2_hyp1


#odds ratios
sen_odds_ratios_hyp1 <- exp(cbind(OR = coef(sen_hyp1), confint(sen_hyp1))) #https://stats.idre.ucla.edu/r/dae/logit-regression/
sen_odds_ratios_hyp1


##SPECIFICITY
spec_hyp1 <- glm(specificity ~ n, family = binomial(link = "logit"), data = results_releveled)
summary(spec_hyp1)

#McFadden’s Pseudo-Rsquared
#use data + noise
TESTspec_hyp1 <- glm(specificity ~ n, family = binomial(link = "logit"), data = data_spec)
#create null model/interaction model
nullmod_spec <- glm(specificity ~ 1, family = binomial(link = "logit"), data = data_spec)
spec_pR2_hyp1 <- 1-logLik(TESTspec_hyp1)/logLik(nullmod_spec)
spec_pR2_hyp1


#odds ratios
spec_odds_ratios_hyp1 <- exp(cbind(OR = coef(spec_hyp1), confint(spec_hyp1)))
spec_odds_ratios_hyp1

##CORRELATION
cor_hyp1 <- glm(abs(correlation) ~ n, family = binomial(link = "logit"), data = results_releveled)
summary(cor_hyp1)

#McFadden’s Pseudo-Rsquared
nullmod_cor <- glm(abs(correlation) ~ 1, family = binomial(link = "logit"), data = results_releveled)
cor_pR2_hyp1 <- 1-logLik(cor_hyp1)/logLik(nullmod_cor)
cor_pR2_hyp1



#odds ratios
cor_odds_ratios_hyp1 <- exp(cbind(OR = coef(cor_hyp1), confint(cor_hyp1)))
cor_odds_ratios_hyp1

##MAE
MAE_hyp1 <- glm(MAE ~ n, family = binomial(link = "logit"), data = results_releveled)
summary(MAE_hyp1)

#McFadden’s Pseudo-Rsquared
#use data + noise
TESTMAE_hyp1 <- glm(MAE ~ n, family = binomial(link = "logit"), data = data_MAE)
#create null/model
nullmod_MAE <- glm(MAE ~ 1, family = binomial(link = "logit"), data = data_MAE)
MAE_pR2_hyp1 <- 1-logLik(TESTMAE_hyp1)/logLik(nullmod_MAE)
MAE_pR2_hyp1 



#odds ratios
MAE_odds_ratios_hyp1 <- exp(cbind(OR = coef(MAE_hyp1), confint(MAE_hyp1)))
MAE_odds_ratios_hyp1




#average outcome 
average_outcome_hyp1 <- results_total %>%
  group_by(n) %>%
  summarise(average_sen = mean(sensitivity),
            average_spec = mean(specificity),
            average_cor = mean(correlation),
            average_MAE = mean(MAE))
average_outcome_hyp1



#plot hyp1----
# ggplot(results_total, aes(x = n, y = sensitivity)) +
#   geom_boxplot() +
#   ylim(0,1)
# 
# ggplot(results_total, aes(x = n, y = specificity)) +
#   geom_boxplot() +
#   ylim(0,1)
# 
# ggplot(results_total, aes(x = n, y = abs(correlation))) +
#   geom_boxplot() +
#   ylim(0,1) +
#   labs(y = "correlation")
# 
# ggplot(results_total, aes(x = n, y = MAE)) +
#   geom_boxplot() +
#   ylim(0, 0.1)


ggplot(long_results, aes(x = as.factor(n), y = outcome_value)) +
  geom_boxplot(fill = '#A4A4A4') +
  ylim(0, 1) +
  labs(x = "Sample size", y = "") +
  theme_bw() +
  facet_grid(. ~ outcome_measure)


#ANALYSES hypothesis 2  ----
#hyp 2 = interaction network type x sample size (n)


##SENSITIVITY
sen_hyp2 <- glm(sensitivity ~ n + network_type + n:network_type, 
                family = binomial(link = "logit"), data = results_releveled)
summary(sen_hyp2)

#McFadden’s Pseudo-Rsquared
sen_pR2_hyp2 <- 1-logLik(sen_hyp2)/logLik(nullmod_sen)
sen_pR2_hyp2 
#compare with model without interaction
sen_hyp2b <- glm(sensitivity ~ n + network_type, 
                 family = binomial(link = "logit"), data = results_releveled)
sen_pR2_hyp2b <- 1-logLik(sen_hyp2b)/logLik(nullmod_sen)
sen_pR2_hyp2b 


#odds ratios
sen_odds_ratios_hyp2 <- exp(cbind(OR = coef(sen_hyp2), confint(sen_hyp2)))
sen_odds_ratios_hyp2


##SPECIFICITY
spec_hyp2 <- glm(specificity ~ n + network_type + n:network_type, 
                 family = binomial(link = "logit"), data = results_releveled)
summary(spec_hyp2)

#McFadden’s Pseudo-Rsquared
#use data + noise
TESTspec_hyp2 <- glm(specificity ~ n + network_type + n:network_type, 
                      family = binomial(link = "logit"), data = data_spec)
spec_pR2_hyp2 <- 1-logLik(TESTspec_hyp2)/logLik(nullmod_spec)
spec_pR2_hyp2
#compare with model without interaction
spec_hyp2b <- glm(specificity ~ n + network_type, 
                       family = binomial(link = "logit"), data = data_spec)
spec_pR2_hyp2b <- 1-logLik(spec_hyp2b)/logLik(nullmod_spec) 
spec_pR2_hyp2b 


#odds ratios
spec_odds_ratios_hyp2 <- exp(cbind(OR = coef(spec_hyp2), confint(spec_hyp2)))
spec_odds_ratios_hyp2


##CORRELATION
cor_hyp2 <- glm(abs(correlation) ~ n + network_type + n:network_type, 
                family = binomial(link = "logit"), data = results_releveled)
summary(cor_hyp2)

#McFadden’s Pseudo-Rsquared
cor_pR2_hyp2 <- 1-logLik(cor_hyp2)/logLik(nullmod_cor)
cor_pR2_hyp2
#compare with model without interaction
cor_hyp2b <- glm(abs(correlation) ~ n + network_type, 
                 family = binomial(link = "logit"), data = results_releveled)
cor_pR2_hyp2b <- 1-logLik(cor_hyp2b)/logLik(nullmod_cor)
cor_pR2_hyp2b 


#odds ratios
cor_odds_ratios_hyp2 <- exp(cbind(OR = coef(cor_hyp2), confint(cor_hyp2)))
cor_odds_ratios_hyp2



##MAE
MAE_hyp2 <- glm(MAE ~ n + network_type + n:network_type, 
                family = binomial(link = "logit"), data = results_releveled)
summary(MAE_hyp2)

#McFadden’s Pseudo-Rsquared
#use data + noise
TESTMAE_hyp2 <- glm(MAE ~ n + network_type + n:network_type, 
                     family = binomial(link = "logit"), data = data_MAE)

MAE_pR2_hyp2 <- 1-logLik(TESTMAE_hyp2)/logLik(nullmod_MAE)
MAE_pR2_hyp2
#compare with model without interaction
MAE_hyp2b <- glm(MAE ~ n + network_type, 
                      family = binomial(link = "logit"), data = data_MAE)
MAE_pR2_hyp2b <- 1-logLik(MAE_hyp2b)/logLik(nullmod_MAE)
MAE_pR2_hyp2b 



#odds ratios
MAE_odds_ratios_hyp2 <- exp(cbind(OR = coef(MAE_hyp2), confint(MAE_hyp2)))
MAE_odds_ratios_hyp2




#average outcome 
average_outcome_hyp2 <- results_total %>%
  group_by(n, network_type) %>%
  summarise(average_sen = mean(sensitivity),
            average_spec = mean(specificity),
            average_cor = mean(correlation),
            average_MAE = mean(MAE))
average_outcome_hyp2



#plot hyp2----
#oud
# ggplot(results_total, aes(x = n, y = sensitivity, fill = network_type)) +
#   geom_boxplot() +
#   ylim(0,1)
# 
# ggplot(results_total, aes(x = n, y = specificity, fill = network_type)) +
#   geom_boxplot() +
#   ylim(0,1)
# 
# ggplot(results_total, aes(x = n, y = abs(correlation), fill = network_type)) +
#   geom_boxplot() +
#   ylim(0,1) +
#   labs(y = "correlation")
# 
# ggplot(results_total, aes(x = n, y = MAE, fill = network_type)) +
#   geom_boxplot() +
#   ylim(0, 0.1)

ggplot(long_results, aes(x = as.factor(n), y = outcome_value, fill = network_type)) +
  geom_boxplot() +
  ylim(0, 1) +
  labs(x = "Sample size", y = "",  fill = "Network type") +
  theme_bw() +
  facet_grid(. ~ outcome_measure)


#ANALYSES hypothesis 3  ----
#hyp 3 = interaction connectivity level x sample size (n)

##SENSITIVITY
sen_hyp3 <- glm(sensitivity ~ n + conn_level + n:conn_level, 
                family = binomial(link = "logit"), data = results_releveled)
summary(sen_hyp3)

#McFadden’s Pseudo-Rsquared
sen_pR2_hyp3 <- 1-logLik(sen_hyp3)/logLik(nullmod_sen)
sen_pR2_hyp3 
#compare with model without interaction
sen_hyp3b <- glm(sensitivity ~ n + conn_level, 
                 family = binomial(link = "logit"), data = results_releveled)
sen_pR2_hyp3b <- 1-logLik(sen_hyp2b)/logLik(nullmod_sen)
sen_pR2_hyp3b 


#odds ratios
sen_odds_ratios_hyp3 <- exp(cbind(OR = coef(sen_hyp3), confint(sen_hyp3)))
sen_odds_ratios_hyp3


##SPECIFICITY
spec_hyp3 <- glm(specificity ~ n + conn_level + n:conn_level, 
                 family = binomial(link = "logit"), data = results_releveled)
summary(spec_hyp3)

#McFadden’s Pseudo-Rsquared
#use data + noise
TESTspec_hyp3 <- glm(specificity ~ n + conn_level + n:conn_level, 
                      family = binomial(link = "logit"), data = data_spec)

spec_pR2_hyp3 <- 1-logLik(TESTspec_hyp3)/logLik(nullmod_spec)
spec_pR2_hyp3
#compare with model without interaction
spec_hyp3b <- glm(specificity ~ n + conn_level, 
                       family = binomial(link = "logit"), data = data_spec)
spec_pR2_hyp3b <- 1-logLik(spec_hyp3b)/logLik(nullmod_spec)
spec_pR2_hyp3b 


#odds ratios
spec_odds_ratios_hyp3 <- exp(cbind(OR = coef(spec_hyp3), confint(spec_hyp3)))
spec_odds_ratios_hyp3


##CORRELATION
cor_hyp3 <- glm(abs(correlation) ~ n + conn_level + n:conn_level, 
                family = binomial(link = "logit"), data = results_releveled)
summary(cor_hyp3)

#McFadden’s Pseudo-Rsquared
cor_pR2_hyp3 <- 1-logLik(cor_hyp3)/logLik(nullmod_cor)
cor_pR2_hyp3
#compare with model without interaction
cor_hyp3b <- glm(abs(correlation) ~ n + conn_level, 
                 family = binomial(link = "logit"), data = results_releveled)
cor_pR2_hyp3b <- 1-logLik(cor_hyp3b)/logLik(nullmod_cor)
cor_pR2_hyp3b


#odds ratios
cor_odds_ratios_hyp3 <- exp(cbind(OR = coef(cor_hyp3), confint(cor_hyp3)))
cor_odds_ratios_hyp3



##MAE
MAE_hyp3 <- glm(MAE ~ n + conn_level + n:conn_level, 
                family = binomial(link = "logit"), data = results_releveled)
summary(MAE_hyp3)

#McFadden’s Pseudo-Rsquared
#use data + noise
TESTMAE_hyp3 <- glm(MAE ~ n + conn_level + n:conn_level, 
                     family = binomial(link = "logit"), data = data_MAE)

MAE_pR2_hyp3 <- 1-logLik(TESTMAE_hyp3)/logLik(nullmod_MAE)
MAE_pR2_hyp3
#compare with model without interaction
MAE_hyp3b <- glm(MAE ~ n + conn_level, 
                      family = binomial(link = "logit"), data = data_MAE)
MAE_pR2_hyp3b <- 1-logLik(MAE_hyp3b)/logLik(nullmod_MAE)
MAE_pR2_hyp3b 


#odds ratios
MAE_odds_ratios_hyp3 <- exp(cbind(OR = coef(MAE_hyp3), confint(MAE_hyp3)))
MAE_odds_ratios_hyp3





#average outcome 
outcome_hyp3 <- results_total %>%
  group_by(n, conn_level) %>%
  summarise(average_sen = mean(sensitivity),
            sd_sen = sd(sensitivity),
            average_spec = mean(specificity),
            sd_spec = sd(specificity),
            average_cor = mean(correlation),
            sd_cor = sd(correlation),
            average_MAE = mean(MAE),
            sd_MAE = sd(MAE))
outcome_hyp3


#plot hyp3 ----
#oud
# ggplot(results_total, aes(x = n, y = sensitivity, fill = conn_level)) +
#   geom_boxplot() +
#   ylim(0,1)
# 
# ggplot(results_total, aes(x = n, y = specificity, fill = conn_level)) +
#   geom_boxplot() +
#   ylim(0,1)
# 
# ggplot(results_total, aes(x = n, y = abs(correlation), fill = conn_level)) +
#   geom_boxplot() +
#   ylim(0,1) +
#   labs(y = "correlation")
# 
# ggplot(results_total, aes(x = n, y = MAE, fill = conn_level)) +
#   geom_boxplot() +
#   ylim(0, 0.1)


ggplot(long_results, aes(x = as.factor(n), y = outcome_value, fill = conn_level)) +
  geom_boxplot() +
  ylim(0, 1) +
  labs(x = "Sample size", y = "",  fill = "Connectivity\nlevel") +
  theme_bw() +
  facet_grid(. ~ outcome_measure)


#nieuw
average_outcome_hyp3 <- outcome_hyp3 %>% 
  pivot_longer(c('average_sen', 'average_spec', 'average_cor', 'average_MAE'), 
               names_to = "outcome_measure", values_to = "outcome_value")

sd_hyp3 <- outcome_hyp3 %>%
  pivot_longer(c('sd_sen', 'sd_spec', 'sd_cor', 'sd_MAE'), 
               names_to = "sd_measure", values_to = "sd_value")



average_sd_hyp3 <- cbind(as.data.frame(average_outcome_hyp3), sd = sd_hyp3$sd_value)
average_sd_hyp3 <- mutate(average_sd_hyp3, ymax = outcome_value + sd,
                          ymin = outcome_value - sd)

#values between 0 and 1
average_sd_hyp3$ymax[average_sd_hyp3$ymax > 1] <- 1
average_sd_hyp3$ymin[average_sd_hyp3$ymin < 0] <- 0


average_sd_hyp3 <- transform(average_sd_hyp3, 
                             outcome_measure = factor(outcome_measure, levels = c('average_sen', 'average_spec', 'average_cor', 'average_MAE'), 
                                                      labels = c("Sensitivity", "Specificity", "Correlation", "MAE")),
                             conn_level = factor(conn_level, levels = c(1, 2, 3, 4),
                                                 labels = c("Lowest", "Low", "Medium", "High")))

ggplot(average_sd_hyp3, aes(x = as.factor(n), y = outcome_value, colour = conn_level, group = conn_level)) +
  # geom_col(position = "dodge", colour = 1) +
  geom_ribbon(aes(ymin = ymin, ymax = ymax, fill = conn_level), alpha = 0.15) +  
  geom_point() +
  geom_line() +
  ylim(0, 1) +
  labs(x = "Sample size", y = "", colour = "Connectivity\nlevel", fill = "Connectivity\nlevel") +
  #  scale_colour_manual(values = c(3, 4, 5, 6)) +
  theme_bw() +
  facet_grid(. ~ outcome_measure)


#ANALYSES hypothesis 4  ----
#hyp 4 = interaction network type x connectivity level x sample size (n)
#      = full model

##SENSITIVITY
sen_hyp4 <- glm(sensitivity ~ n + conn_level + network_type +
                  n:network_type + n:conn_level + n:network_type:conn_level, 
                family = "binomial", data = results_releveled)
summary(sen_hyp4)

#McFadden’s Pseudo-Rsquared
sen_pR2_hyp4 <- 1-logLik(sen_hyp4)/logLik(nullmod_sen)
sen_pR2_hyp4
#compare with model without 3-way interaction
sen_hyp4b <- glm(sensitivity ~ n + conn_level + network_type +
                   n:network_type + n:conn_level, 
                 family = binomial(link = "logit"), data = results_releveled)
sen_pR2_hyp4b <- 1-logLik(sen_hyp4b)/logLik(nullmod_sen)
sen_pR2_hyp4b 


#odds ratios
sen_odds_ratios_hyp4 <- exp(cbind(OR = coef(sen_hyp4), confint(sen_hyp4)))
sen_odds_ratios_hyp4




##SPECIFICITY
spec_hyp4 <- glm(specificity ~ n + conn_level + network_type +
                   n:network_type + n:conn_level + n:network_type:conn_level, 
                 family = "binomial", data = results_releveled)
summary(spec_hyp4)

#McFadden’s Pseudo-Rsquared
#use data + noise
TESTspec_hyp4 <- glm(specificity ~ n + conn_level + network_type +
                        n:network_type + n:conn_level + n:network_type:conn_level, 
                      family = "binomial", data = data_spec)

spec_pR2_hyp4 <- 1-logLik(TESTspec_hyp4)/logLik(nullmod_spec)
spec_pR2_hyp4
#compare with model without 3-way interaction
spec_hyp4b <- glm(specificity ~ n + conn_level + network_type +
                         n:network_type + n:conn_level, 
                       family = binomial(link = "logit"), data = data_spec)
spec_pR2_hyp4b <- 1-logLik(spec_hyp4b)/logLik(nullmod_spec)
spec_pR2_hyp4b 


#odds ratios
spec_odds_ratios_hyp4 <- exp(cbind(OR = coef(spec_hyp4), confint(spec_hyp4)))
spec_odds_ratios_hyp4


##CORRELATION
cor_hyp4 <- glm(abs(correlation) ~ n + conn_level + network_type +
                  n:network_type + n:conn_level + n:network_type:conn_level, 
                family = "binomial", data = results_releveled)
summary(cor_hyp4)

#McFadden’s Pseudo-Rsquared
cor_pR2_hyp4 <- 1-logLik(cor_hyp4)/logLik(nullmod_cor)
cor_pR2_hyp4
#compare with model without 3-way interaction
cor_hyp4b <- glm(abs(correlation) ~ n + conn_level + network_type +
                   n:network_type + n:conn_level, 
                 family = binomial(link = "logit"), data = results_releveled)
cor_pR2_hyp4b <- 1-logLik(cor_hyp4b)/logLik(nullmod_cor)
cor_pR2_hyp4b 


#odds ratios
cor_odds_ratios_hyp4 <- exp(cbind(OR = coef(cor_hyp4), confint(cor_hyp4)))
cor_odds_ratios_hyp4



##MAE
MAE_hyp4 <- glm(MAE ~ n + conn_level + network_type +
                  n:network_type + n:conn_level + n:network_type:conn_level, 
                family = "binomial", data = results_releveled)
summary(MAE_hyp4)

#McFadden’s Pseudo-Rsquared
#use data + noise
TESTMAE_hyp4 <- glm(MAE ~ n + conn_level + network_type +
                       n:network_type + n:conn_level + n:network_type:conn_level, 
                     family = "binomial", data = data_MAE)

MAE_pR2_hyp4 <- 1-logLik(TESTMAE_hyp4)/logLik(nullmod_MAE) 
MAE_pR2_hyp4
#compare with model without 3-way interaction
MAE_hyp4b <- glm(MAE ~ n + conn_level + network_type +
                        n:network_type + n:conn_level, 
                      family = binomial(link = "logit"), data = data_MAE)
MAE_pR2_hyp4b <- 1-logLik(MAE_hyp4b)/logLik(nullmod_MAE)
MAE_pR2_hyp4b 


#odds ratios
MAE_odds_ratios_hyp4 <- exp(cbind(OR = coef(MAE_hyp4), confint(MAE_hyp4)))
MAE_odds_ratios_hyp4

#average outcome 
outcome_hyp4 <- results_total %>%
  group_by(n, network_type, conn_level) %>%
  summarise(average_sen = mean(sensitivity),
            sd_sen = sd(sensitivity),
            average_spec = mean(specificity),
            sd_spec = sd(specificity),
            average_cor = mean(correlation),
            sd_cor = sd(correlation),
            average_MAE = mean(MAE),
            sd_MAE = sd(MAE))
outcome_hyp4



min_sen06 <- filter(outcome_hyp4, average_sen > 0.6)
min_cor08 <- filter(outcome_hyp4, average_cor > 0.8)




#plot hyp4----
average_outcome_full <- outcome_hyp4 %>% 
  pivot_longer(c('average_sen', 'average_spec', 'average_cor', 'average_MAE'), 
               names_to = "outcome_measure", values_to = "outcome_value")

sd_full <- outcome_hyp4 %>%
  pivot_longer(c('sd_sen', 'sd_spec', 'sd_cor', 'sd_MAE'), 
               names_to = "sd_measure", values_to = "sd_value")



average_sd_full <- cbind(as.data.frame(average_outcome_full), sd = sd_full$sd_value)
average_sd_full <- mutate(average_sd_full, ymax = outcome_value + sd,
                          ymin = outcome_value - sd)

#values between 0 and 1
average_sd_full$ymax[average_sd_full$ymax > 1] <- 1
average_sd_full$ymin[average_sd_full$ymin < 0] <- 0


average_sd_full <- transform(average_sd_full, 
                             outcome_measure = factor(outcome_measure, levels = c('average_sen', 'average_spec', 'average_cor', 'average_MAE'), 
                                                      labels = c("Sensitivity", "Specificity", "Correlation", "MAE"))
                             ,
                             conn_level = factor(conn_level, levels = c(1, 2, 3, 4),
                                                 labels = c("Lowest", "Low", "Medium", "High")),
                             network_type = factor(network_type, levels = c("random", "scalefree", "smallworld"),
                                                   labels = c("Random", "Scale free", "Small world")))



ggplot(average_sd_full, aes(x = as.factor(n), y = outcome_value, colour = conn_level, group = conn_level)) +
  # geom_col(position = "dodge", colour = 1) +
  geom_ribbon(aes(ymin = ymin, ymax = ymax, fill = conn_level), alpha = 0.15) +  
  geom_point() +
  geom_line() +
  ylim(0, 1) +
  labs(x = "Sample size", y = "", colour = "Connectivity\nlevel", fill = "Connectivity\nlevel") +
  #  scale_colour_manual(values = c(3, 4, 5, 6)) +
  theme_bw() +
  facet_grid(network_type ~ outcome_measure)





#plot % empty networks  ----

empty_networks <- mutate(results_total, empty = ifelse(sensitivity == 0 & specificity == 1, 1, 0)) %>%
  group_by(n, conn_level, network_type) %>%
  summarise(perc_empty = sum(empty)/1000) %>%
  transform(conn_level = factor(conn_level, levels = c(1, 2, 3, 4),
                                labels = c("Lowest", "Low", "Medium", "High")),
            network_type = factor(network_type, levels = c("random", "scalefree", "smallworld"),
                                  labels = c("Random", "Scale free", "Small world")))

ggplot(empty_networks, aes(x = as.factor(n), y = perc_empty, colour = conn_level, group = conn_level)) +
  geom_point() +
  geom_line() +
  # ylim(0, 1) +
  labs(x = "Sample size", y = "", colour = "Connectivity\nlevel", fill = "Connectivity\nlevel") +
  scale_y_continuous(labels = scales::percent) +
  theme_bw() +
  facet_grid(. ~ network_type)

