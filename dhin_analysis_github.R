#################
# Analysis of mammography uptake and breast cancer
# Citation: Goldstein ND, Zhang Y, Burstyn I, Vasireddy K, Garland D, Enright M, Siegel SD. Towards improved methodology to measure the uptake of mammography in populations using insurance claims: a case study from Delaware, USA. Manuscript in preparation.
# 5/15/25 -- Neal Goldstein
#################


### FUNCTIONS ###

library(sf) #spatial functions
library(spdep) #spatial functions
library(spatialreg) #spatial regression


### READ DATA ###

#DHIN and census data (using Yuchen's scripts)
load("Final_data_before_simulation.RData")

#breast cancer cases (uses 2020 boundaries)
cancer = read.csv("DE_tract_aggr_by_SEER_stage.csv", as.is=T, stringsAsFactors=F, na.strings="")


### DATA PREPROCESSING ###

#two tracts in 2022 have populations receiving public insurance greater than population, so cap those based on minimum of simular tracts
final_data_before_simulation_2022$shape2[final_data_before_simulation_2022$GEOID=="10003010106"] = 10
final_data_before_simulation_2022$shape2[final_data_before_simulation_2022$GEOID=="10005051202"] = 10

#reshape breast cancer cases to wide
cancer_wide=reshape(cancer, idvar="GEOID", timevar="Merged_Summary_Stage", direction="wide")

#create a total breast cancer count
cancer_wide$count.Total = rowSums(cancer_wide[, c("count.Distant","count.Local","count.Regional","count.Unknown")], na.rm=T)

#identify late cancer cases (subtract unknown counts if available)
cancer_wide$count.NonLocal = cancer_wide$count.Total - cancer_wide$count.Local - ifelse(is.na(cancer_wide$count.Unknown), 0, cancer_wide$count.Unknown)

#join to mammography data
final_data_before_simulation_2022 = merge(x=final_data_before_simulation_2022, y=cancer_wide[, c("GEOID","count.Local","count.NonLocal","count.Total")], by="GEOID", all.x=T, duplicateGeoms=T)

#define neighbors using Queen's contiguity (spdep package)
nb_list = poly2nb(final_data_before_simulation_2022)

#define weighting scheme
wt_list = nb2listw(nb_list) #using row standardization (default)


### PERFORM QBA ###

#number of simulations
nsim = 10000

#dataframes to store results
simulation_results_2019 = data.frame("Simulation"=rep(1:nsim, each=length(unique(final_data_before_simulation_2019$GEOID))), "GEOID"=rep(unique(final_data_before_simulation_2019$GEOID), nsim), "sim_gamma"=NA, "sim_eta"=NA, "sim_m_eq1"=NA, "sim_m_eq2"=NA, "sim_m_avg"=NA, stringsAsFactors=F)
simulation_results_2022 = data.frame("Simulation"=rep(1:nsim, each=length(unique(final_data_before_simulation_2022$GEOID))), "GEOID"=rep(unique(final_data_before_simulation_2022$GEOID), nsim), "sim_gamma"=NA, "sim_eta"=NA, "sim_m_eq1"=NA, "sim_m_eq2"=NA, "sim_m_avg"=NA, stringsAsFactors=F)
simulation_results_arlm_2022 = data.frame("Simulation"=seq(1:nsim), "Total"=NA, "Local"=NA, "NonLocal"=NA, stringsAsFactors=F)

for (i in 1:nsim) {
  
  cat("\n\n************** ","Simulation: ",i," **************\n",sep="")
  
  ## EQUATION 1 ##
  #M_i = N_i * P_i / gamma_i
  #N_i = CensusFemaleOver40
  #O_i = DHIN_mammo_members
  #P_i = beta(O_i + 1, N_i - O_i + 1)
  #T_i = DHIN_total_members
  #gamma_i = unif(T_i/N_i, 1)
  
  #perform stochastic draws, 2019
  sim_p_2019 = rbeta(nrow(final_data_before_simulation_2019), (final_data_before_simulation_2019$DHIN_mammo_members + 1), (final_data_before_simulation_2019$CensusFemaleOver40 - final_data_before_simulation_2019$DHIN_mammo_members + 1))
  sim_gamma_2019 = runif(nrow(final_data_before_simulation_2019), (final_data_before_simulation_2019$DHIN_total_members / final_data_before_simulation_2019$CensusFemaleOver40), 1)
  
  #apply constraints to sensitivity (must be >= observed prevalence and <1)
  sim_gamma_2019 = ifelse(sim_gamma_2019<(final_data_before_simulation_2019$DHIN_mammo_members/final_data_before_simulation_2019$CensusFemaleOver40), (final_data_before_simulation_2019$DHIN_mammo_members/final_data_before_simulation_2019$CensusFemaleOver40), sim_gamma_2019)
  sim_gamma_2019 = ifelse(sim_gamma_2019>1, 1, sim_gamma_2019)
  
  #calculate eq 1
  sim_m_eq1_2019 = final_data_before_simulation_2019$CensusFemaleOver40 * sim_p_2019 / sim_gamma_2019
  
  #perform stochastic draws, 2022
  sim_p_2022 = rbeta(nrow(final_data_before_simulation_2022), final_data_before_simulation_2022$DHIN_mammo_members, (final_data_before_simulation_2022$CensusFemaleOver40 - final_data_before_simulation_2022$DHIN_mammo_members))
  sim_gamma_2022 = runif(nrow(final_data_before_simulation_2022), (final_data_before_simulation_2022$DHIN_total_members / final_data_before_simulation_2022$CensusFemaleOver40), 1)
  
  #apply constraints to sensitivity (must be >= observed prevalence and <1)
  sim_gamma_2022 = ifelse(sim_gamma_2022<(final_data_before_simulation_2022$DHIN_mammo_members/final_data_before_simulation_2022$CensusFemaleOver40), (final_data_before_simulation_2022$DHIN_mammo_members/final_data_before_simulation_2022$CensusFemaleOver40), sim_gamma_2022)
  sim_gamma_2022 = ifelse(sim_gamma_2022>1, 1, sim_gamma_2022)
  
  #calculate eq 1
  sim_m_eq1_2022 = final_data_before_simulation_2022$CensusFemaleOver40 * sim_p_2022 / sim_gamma_2022
  
  ## EQUATION 2 ##
  #M_i = N_i * P_i / eta_i
  #eta_i = beta_i + ((gamma_i_star) * (1 - beta_i))
  #beta_i = beta(shape1, shape2)
  #gamma_i_star = unif(0, 0.6)
  
  #perform stochastic draws, 2019
  sim_beta_2019 = rbeta(nrow(final_data_before_simulation_2019), final_data_before_simulation_2019$shape1, final_data_before_simulation_2019$shape2)
  sim_gamma_star_2019 = runif(nrow(final_data_before_simulation_2019), 0, 0.6)
  sim_eta_2019 = sim_beta_2019 + (sim_gamma_star_2019 * (1-sim_beta_2019))
  
  #apply constraints to sensitivity (must be >= observed prevalence and <1)
  sim_eta_2019 = ifelse(sim_eta_2019<(final_data_before_simulation_2019$DHIN_mammo_members/final_data_before_simulation_2019$CensusFemaleOver40), (final_data_before_simulation_2019$DHIN_mammo_members/final_data_before_simulation_2019$CensusFemaleOver40), sim_eta_2019)
  sim_eta_2019 = ifelse(sim_eta_2019>1, 1, sim_eta_2019)
  
  #calculate eq 2
  sim_m_eq2_2019 = final_data_before_simulation_2019$CensusFemaleOver40 * sim_p_2019 / sim_eta_2019
  
  #perform stochastic draws, 2022
  sim_beta_2022 = rbeta(nrow(final_data_before_simulation_2022), final_data_before_simulation_2022$shape1, final_data_before_simulation_2022$shape2)
  sim_gamma_star_2022 = runif(nrow(final_data_before_simulation_2022), 0, 0.6)
  sim_eta_2022 = sim_beta_2022 + (sim_gamma_star_2022 * (1-sim_beta_2022))
  
  #apply constraints to sensitivity (must be >= observed prevalence and <1)
  sim_eta_2022 = ifelse(sim_eta_2022<(final_data_before_simulation_2022$DHIN_mammo_members/final_data_before_simulation_2022$CensusFemaleOver40), (final_data_before_simulation_2022$DHIN_mammo_members/final_data_before_simulation_2022$CensusFemaleOver40), sim_eta_2022)
  sim_eta_2022 = ifelse(sim_eta_2022>1, 1, sim_eta_2022)
  
  #calculate eq 2
  sim_m_eq2_2022 = final_data_before_simulation_2022$CensusFemaleOver40 * sim_p_2022 / sim_eta_2022
  
  #average between two methods
  sim_m_avg_2019 = (sim_m_eq1_2019 + sim_m_eq2_2019) / 2
  sim_m_avg_2022 = (sim_m_eq1_2022 + sim_m_eq2_2022) / 2
  
  ## SPATIAL REG, 2022 only ##
  arlm_m_total = errorsarlm(final_data_before_simulation_2022$count.Total ~ sim_m_avg_2022, weights=final_data_before_simulation_2022$CensusFemaleOver40, listw=wt_list)
  arlm_m_local = errorsarlm(final_data_before_simulation_2022$count.Local ~ sim_m_avg_2022, weights=final_data_before_simulation_2022$CensusFemaleOver40, listw=wt_list)
  arlm_m_nonlocal = errorsarlm(final_data_before_simulation_2022$count.NonLocal ~ sim_m_avg_2022, weights=final_data_before_simulation_2022$CensusFemaleOver40, listw=wt_list)
  
  #sample point estimate
  sim_arlm_total = rnorm(1, arlm_m_total$coefficients[2], arlm_m_total$rest.se[2])
  sim_arlm_local = rnorm(1, arlm_m_local$coefficients[2], arlm_m_local$rest.se[2])
  sim_arlm_nonlocal = rnorm(1, arlm_m_nonlocal$coefficients[2], arlm_m_nonlocal$rest.se[2])
  
  ## SAVE RESULTS ##
  simulation_results_2019$sim_gamma[simulation_results_2019$Simulation==i] = sim_gamma_2019
  simulation_results_2019$sim_eta[simulation_results_2019$Simulation==i] = sim_eta_2019
  simulation_results_2019$sim_m_eq1[simulation_results_2019$Simulation==i] = sim_m_eq1_2019
  simulation_results_2019$sim_m_eq2[simulation_results_2019$Simulation==i] = sim_m_eq2_2019
  simulation_results_2019$sim_m_avg[simulation_results_2019$Simulation==i] = sim_m_avg_2019
  
  simulation_results_2022$sim_gamma[simulation_results_2022$Simulation==i] = sim_gamma_2022
  simulation_results_2022$sim_eta[simulation_results_2022$Simulation==i] = sim_eta_2022
  simulation_results_2022$sim_m_eq1[simulation_results_2022$Simulation==i] = sim_m_eq1_2022
  simulation_results_2022$sim_m_eq2[simulation_results_2022$Simulation==i] = sim_m_eq2_2022
  simulation_results_2022$sim_m_avg[simulation_results_2022$Simulation==i] = sim_m_avg_2022
  
  simulation_results_arlm_2022$Total[i] = sim_arlm_total
  simulation_results_arlm_2022$Local[i] = sim_arlm_local
  simulation_results_arlm_2022$NonLocal[i] = sim_arlm_nonlocal
  
  #clean up
  rm(sim_arlm_local, sim_arlm_nonlocal, sim_arlm_total, sim_beta_2019, sim_beta_2022, sim_eta_2019, sim_eta_2022, sim_gamma_2019, sim_gamma_2022, sim_gamma_star_2019, sim_gamma_star_2022, sim_m_avg_2019, sim_m_avg_2022, sim_m_eq1_2019, sim_m_eq1_2022, sim_m_eq2_2019, sim_m_eq2_2022, sim_p_2019, sim_p_2022, arlm_m_local, arlm_m_nonlocal, arlm_m_total)
  gc()
}

#save results
save.image("Mammo_simulation_results.RData")


### ANALYSES FOR PAPER ###

#read data
load("Mammo_simulation_results.RData")

#functions
library(ggplot2)
library(sf) #spatial functions
library(spatialreg) #spatial regression


### MAMMOGRAPHY ANALYSIS ###

## pre-adjustment ##
#total mammograms, raw and per 1k women
sum(final_data_before_simulation_2019$DHIN_mammo_members)
sum(final_data_before_simulation_2019$DHIN_mammo_members)/sum(final_data_before_simulation_2019$CensusFemaleOver40)*1000

sum(final_data_before_simulation_2022$DHIN_mammo_members)
sum(final_data_before_simulation_2022$DHIN_mammo_members)/sum(final_data_before_simulation_2022$CensusFemaleOver40)*1000

#by tract, per 1k women
min(final_data_before_simulation_2019$Mammo_per_1000_census)
max(final_data_before_simulation_2019$Mammo_per_1000_census)

min(final_data_before_simulation_2022$Mammo_per_1000_census)
max(final_data_before_simulation_2022$Mammo_per_1000_census)

## post-adjustment ##
quantile(by(simulation_results_2019$sim_m_avg, simulation_results_2019$Simulation, sum)/sum(final_data_before_simulation_2019$CensusFemaleOver40)*1000, probs=c(0.025,0.5,0.975))
quantile(by(simulation_results_2022$sim_m_avg, simulation_results_2022$Simulation, sum)/sum(final_data_before_simulation_2022$CensusFemaleOver40)*1000, probs=c(0.025,0.5,0.975))

#method-specific results
quantile(by(simulation_results_2019$sim_m_eq1, simulation_results_2019$Simulation, sum)/sum(final_data_before_simulation_2019$CensusFemaleOver40)*1000, probs=c(0.025,0.5,0.975))
quantile(by(simulation_results_2019$sim_m_eq2, simulation_results_2019$Simulation, sum)/sum(final_data_before_simulation_2019$CensusFemaleOver40)*1000, probs=c(0.025,0.5,0.975))
quantile(by(simulation_results_2022$sim_m_eq1, simulation_results_2022$Simulation, sum)/sum(final_data_before_simulation_2022$CensusFemaleOver40)*1000, probs=c(0.025,0.5,0.975))
quantile(by(simulation_results_2022$sim_m_eq2, simulation_results_2022$Simulation, sum)/sum(final_data_before_simulation_2022$CensusFemaleOver40)*1000, probs=c(0.025,0.5,0.975))

#sensitivity
median(by(simulation_results_2019$sim_gamma, simulation_results_2019$GEOID, mean))
min(by(simulation_results_2019$sim_gamma, simulation_results_2019$GEOID, mean))
max(by(simulation_results_2019$sim_gamma, simulation_results_2019$GEOID, mean))
median(by(simulation_results_2019$sim_eta, simulation_results_2019$GEOID, mean))
min(by(simulation_results_2019$sim_eta, simulation_results_2019$GEOID, mean))
max(by(simulation_results_2019$sim_eta, simulation_results_2019$GEOID, mean))
cor(by(simulation_results_2019$sim_gamma, simulation_results_2019$GEOID, mean), by(simulation_results_2019$sim_eta, simulation_results_2019$GEOID, mean), method="spearman")

plot(density(simulation_results_2019$sim_gamma), type="l", xlim=c(0,1), xlab="Sensitivity", main="")
par(new=T)
plot(density(simulation_results_2019$sim_eta), type="l", lty=2, xlim=c(0,1), ylim=c(0,10), yaxt="n", xaxt="n", ann=F)
legend("topleft", legend=c(expression(gamma),expression(eta)), lty=c(1,2))

median(by(simulation_results_2022$sim_gamma, simulation_results_2022$GEOID, mean))
min(by(simulation_results_2022$sim_gamma, simulation_results_2022$GEOID, mean))
max(by(simulation_results_2022$sim_gamma, simulation_results_2022$GEOID, mean))
median(by(simulation_results_2022$sim_eta, simulation_results_2022$GEOID, mean))
min(by(simulation_results_2022$sim_eta, simulation_results_2022$GEOID, mean))
max(by(simulation_results_2022$sim_eta, simulation_results_2022$GEOID, mean))
cor(by(simulation_results_2022$sim_gamma, simulation_results_2022$GEOID, mean), by(simulation_results_2022$sim_eta, simulation_results_2022$GEOID, mean), method="spearman")

plot(density(simulation_results_2022$sim_gamma), type="l", xlim=c(0,1), xlab="Sensitivity", main="")
par(new=T)
plot(density(simulation_results_2022$sim_eta), type="l", lty=2, xlim=c(0,1), ylim=c(0,10), yaxt="n", xaxt="n", ann=F)
legend("topleft", legend=c(expression(gamma),expression(eta)), lty=c(1,2))


### MAPS ###

#unadjusted, 2019
ggplot()+
  geom_sf(data = final_data_before_simulation_2019, aes( fill=Mammo_per_1000_census)) + 
  #geom_sf(test_map,mapping=aes(colour = 'red')) +
  #scale_colour_gradientn(colours = c("#BD0026", "#F03B20", "#FD8D3C", "#FEB24C", "#FED976", "#FFFFB2"))+
  scale_fill_viridis_c(direction=-1)+   # color shows from lightest to darkest, if 1(defualt), it's from darkest to lightest
  labs(fill='') +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

#unadjusted, 2022
ggplot()+
  geom_sf(data = final_data_before_simulation_2022, aes( fill=Mammo_per_1000_census)) + 
  #geom_sf(test_map,mapping=aes(colour = 'red')) +
  #scale_colour_gradientn(colours = c("#BD0026", "#F03B20", "#FD8D3C", "#FEB24C", "#FED976", "#FFFFB2"))+
  scale_fill_viridis_c(direction=-1)+   # color shows from lightest to darkest, if 1(defualt), it's from darkest to lightest
  labs(fill='') +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

#adjusted, 2019
adjusted_map_2019 = final_data_before_simulation_2019
adjusted_map_2019$m_avg = as.numeric(by(simulation_results_2019$sim_m_avg, simulation_results_2019$GEOID, median))
ggplot()+
  geom_sf(data = adjusted_map_2019, aes( fill=m_avg)) + 
  #geom_sf(test_map,mapping=aes(colour = 'red')) +
  #scale_colour_gradientn(colours = c("#BD0026", "#F03B20", "#FD8D3C", "#FEB24C", "#FED976", "#FFFFB2"))+
  scale_fill_viridis_c(direction=-1)+   # color shows from lightest to darkest, if 1(defualt), it's from darkest to lightest
  labs(fill='') +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

#adjusted, 2022
adjusted_map_2022 = final_data_before_simulation_2022
adjusted_map_2022$m_avg = as.numeric(by(simulation_results_2022$sim_m_avg, simulation_results_2022$GEOID, median))
ggplot()+
  geom_sf(data = adjusted_map_2022, aes( fill=m_avg)) + 
  #geom_sf(test_map,mapping=aes(colour = 'red')) +
  #scale_colour_gradientn(colours = c("#BD0026", "#F03B20", "#FD8D3C", "#FEB24C", "#FED976", "#FFFFB2"))+
  scale_fill_viridis_c(direction=-1)+   # color shows from lightest to darkest, if 1(defualt), it's from darkest to lightest
  labs(fill='') +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

#breast cancer
final_data_before_simulation_2022$cancer_rate = final_data_before_simulation_2022$count.Total/final_data_before_simulation_2022$CensusFemaleOver40*1000
ggplot()+
  geom_sf(data = final_data_before_simulation_2022, aes( fill=cancer_rate)) + 
  #geom_sf(test_map,mapping=aes(colour = 'red')) +
  #scale_colour_gradientn(colours = c("#BD0026", "#F03B20", "#FD8D3C", "#FEB24C", "#FED976", "#FFFFB2"))+
  scale_fill_viridis_c(direction=-1)+   # color shows from lightest to darkest, if 1(defualt), it's from darkest to lightest
  labs(fill='') +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


### BLANT ALTMAN PLOT ###

#2019
#see https://www.statology.org/bland-altman-plot-r/
ba_plot_data = data.frame("Eq1"=by(simulation_results_2019$sim_m_eq1, simulation_results_2019$GEOID, median), "Eq2"=by(simulation_results_2019$sim_m_eq2, simulation_results_2019$GEOID, median), "Avg"=by(simulation_results_2019$sim_m_avg, simulation_results_2019$GEOID, median))
ba_plot_data$Diff = ba_plot_data$Eq2 - ba_plot_data$Eq1

ggplot(ba_plot_data, aes(x = Avg, y = Diff)) +
  geom_point(size=2) +
  geom_hline(yintercept = mean(ba_plot_data$Diff)) +
  geom_hline(yintercept = (mean(ba_plot_data$Diff) - 1.96*sd(ba_plot_data$Diff)), linetype="dashed") +
  geom_hline(yintercept = (mean(ba_plot_data$Diff) + 1.96*sd(ba_plot_data$Diff)), linetype="dashed") +
  ylab("Difference between sensitivity specifications") +
  xlab("Mammograms by averaging two specifications (post-adjustment)")

#2022
ba_plot_data = data.frame("Eq1"=by(simulation_results_2022$sim_m_eq1, simulation_results_2022$GEOID, median), "Eq2"=by(simulation_results_2022$sim_m_eq2, simulation_results_2022$GEOID, median), "Avg"=by(simulation_results_2022$sim_m_avg, simulation_results_2022$GEOID, median))
ba_plot_data$Diff = ba_plot_data$Eq2 - ba_plot_data$Eq1

ggplot(ba_plot_data, aes(x = Avg, y = Diff)) +
  geom_point(size=2) +
  geom_hline(yintercept = mean(ba_plot_data$Diff)) +
  geom_hline(yintercept = (mean(ba_plot_data$Diff) - 1.96*sd(ba_plot_data$Diff)), linetype="dashed") +
  geom_hline(yintercept = (mean(ba_plot_data$Diff) + 1.96*sd(ba_plot_data$Diff)), linetype="dashed") +
  ylab("Difference between sensitivity specifications") +
  xlab("Mammograms by averaging two specifications (post-adjustment)")


### BREAST CANCER ANALYSIS ###

sum(final_data_before_simulation_2022$count.Total)
sum(final_data_before_simulation_2022$count.Total) - sum(final_data_before_simulation_2022$count.NonLocal)
sum(final_data_before_simulation_2022$count.NonLocal)

#Moran's I for rate data, use empirical Bayes index modification
#checking for both dependent and independent variables
EBImoran.mc(n=final_data_before_simulation_2022$count.Total, x=final_data_before_simulation_2022$CensusFemaleOver40, listw=wt_list, nsim=10000)
EBImoran.mc(n=final_data_before_simulation_2022$count.Local, x=final_data_before_simulation_2022$CensusFemaleOver40, listw=wt_list, nsim=10000)
EBImoran.mc(n=final_data_before_simulation_2022$count.NonLocal, x=final_data_before_simulation_2022$CensusFemaleOver40, listw=wt_list, nsim=10000)
EBImoran.mc(n=final_data_before_simulation_2022$DHIN_mammo_members, x=final_data_before_simulation_2022$CensusFemaleOver40, listw=wt_list, nsim=10000)

#OLS model
lm_naive = lm(count.Total ~ DHIN_mammo_members, weights=CensusFemaleOver40, data=final_data_before_simulation_2022)
lm_naive = lm(count.Local ~ DHIN_mammo_members, weights=CensusFemaleOver40, data=final_data_before_simulation_2022)
lm_naive = lm(count.NonLocal ~ DHIN_mammo_members, weights=CensusFemaleOver40, data=final_data_before_simulation_2022)

summary(lm_naive)

#test the residuals for spatial autocorrelation
moran.test(lm_naive$residuals, wt_list, randomisation=F, alternative="two.sided") 
lm.morantest(lm_naive, wt_list) #for linear models

## pre-adjustment ##
#spatial error model
arlm_naive = errorsarlm(count.Total ~ DHIN_mammo_members, weights=CensusFemaleOver40, data=final_data_before_simulation_2022, wt_list)
summary(arlm_naive)$coefficients*1000
confint(arlm_naive)*1000

#test for improved fit over OLS model
anova(arlm_naive, lm_naive)

#spatial error model, local cases
arlm_naive = errorsarlm(count.Local ~ DHIN_mammo_members, weights=CensusFemaleOver40, data=final_data_before_simulation_2022, wt_list)
summary(arlm_naive)$coefficients*1000
confint(arlm_naive)*1000

#spatial error model, non local cases
arlm_naive = errorsarlm(count.NonLocal ~ DHIN_mammo_members, weights=CensusFemaleOver40, data=final_data_before_simulation_2022, wt_list)
summary(arlm_naive)$coefficients*1000
confint(arlm_naive)*1000

## post-adjustment ##
quantile(simulation_results_arlm_2022$Total*1000, probs=c(0.025,0.5,0.975))
quantile(simulation_results_arlm_2022$Local*1000, probs=c(0.025,0.5,0.975))
quantile(simulation_results_arlm_2022$NonLocal*1000, probs=c(0.025,0.5,0.975))

#percent change
(77-126)/126 #total
(55-90)/90 #local
(22-37)/37 #nonlocal
