Are rhizobia under selection to cheat?
================
Megan Frederickson
2019-04-17

Fitness conflict or fitness alignment?
--------------------------------------

This repository re-analyzes the data in:

Gano-Cohen KA, Wendlandt CE, Stokes PJ, Blanton MA, Quides KW, Zomorrodian A, Adinata ES, Sachs JL (2019) Interspecific conflict and the evolution of ineffective rhizobia. Ecology Letters. <https://doi.org/10.1111/ele.13247>

The authors make the case that in their legume-rhizobium study system, the rhizobia are selected to cheat. In other words, they report a negative correlation between legume and rhizobium fitnesses.

Here, I re-analyze their data using a standard genetic selection analysis approach, in which relativized fitness is regressed against standardized family-mean trait values.

First we need to load some nifty R packages.

``` r
#Load useful packages
library(tidyverse) #includes ggplot2, dplyr, readr, stringr
library(cowplot)
library(car)
library(glmm)
library(nlme)
library(lme4)
library(emmeans)
```

I downloaded the data from Dryad on April 16, 2019. The citation for the data package is:

Gano-Cohen KA, Wendlandt CE, Stokes PJ, Blanton MA, Quides KW, Zomorrodian A, Adinata ES, Sachs JL (2019) Data from: Interspecific conflict and the evolution of ineffective rhizobia. Dryad Digital Repository. <https://doi.org/10.5061/dryad.cr65269>

``` r
data <- read_csv("Table_S4.csv")
```

Calculate genotype means
========================

The data set has two measures of plant fitness and four measures of rhizobium fitness. First, we need to calculate the appropriate strain means for each one. I did this using the emmeans package, which calculates estimated marginal means from a linear mixed model. I have tried to stick as close as possible to the analysis presented in the paper, which states that "all effects were coded as fixed except block, which was treated as a random effect" and that variables were log-transformed as needed to improve normality. I also calculated strain means only on data from sympatric hosts, again following the paper's general modelling approach.

``` r
#Exclude control, uninoculated plants

data <- subset(data, Strain != "control")

#Exclude plants that did not form nodules

data <- subset(data, `Total nodules` > 0)

#Make sure factors are factors and numbers are numbers

data$Strain <- as.factor(data$Strain)
data$`Host Line` <- as.factor(data$`Host Line`)
data$Population <- as.factor(data$Population)
data$`Total nodules` <- as.numeric(data$`Total nodules`)
data$`Mean individual  nodule biomass (mg)` <- as.numeric(data$`Mean individual  nodule biomass (mg)`)
data$`Shoots mass (g)` <- as.numeric(data$`Shoots mass (g)`)
data$`Roots mass (g)` <- as.numeric(data$`Roots mass (g)`)
data$`Relative Growth` <- as.numeric(data$`Relative Growth`)
data$Block <- as.factor(data$Block)

#Log transform as needed, as per paper

data$log_nodules <- log(data$`Total nodules`)
data$log_mean_nod_biomass <- log(data$`Mean individual  nodule biomass (mg)`)

#Combine shoot and root mass into one measure of plant biomass
data$plant_biomass <- data$`Shoots mass (g)` + data$`Roots mass (g)`

#Divide into sympatric and universal hosts, as per paper

data_sym <- subset(data, `Host Line` != "A. heermannii" & `Host Line` != "UnH: Cla12.04" & `Host Line` != "UnL: Anz13.04")
data_uni <- subset(data, `Host Line` == "A. heermannii" | `Host Line` == "UnH: Cla12.04" | `Host Line` == "UnL: Anz13.04")

#Model total nodules
lmm1 <- lmer(log_nodules~Population + Strain + `Host Line` + (1|Block), data=data_sym) 
#plot(lmm1)
#summary(lmm1)
Anova(lmm1, type=3)
```

    ## Analysis of Deviance Table (Type III Wald chisquare tests)
    ## 
    ## Response: log_nodules
    ##               Chisq Df Pr(>Chisq)    
    ## (Intercept) 244.202  1  < 2.2e-16 ***
    ## Population   15.514  5   0.008377 ** 
    ## Strain      117.317 20  8.931e-16 ***
    ## `Host Line`  18.413  6   0.005279 ** 
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
#Extract line means and save them to new dataframe

df.means <- as.data.frame(emmeans(lmm1, "Strain"))
colnames(df.means) <- c("Strain", "Population", "tot_nod_lsmean", "tot_nod_SE", "tot_nod_df", "tot_nod_lowCL", "tot_nod_highCL")

#Mean nodule mass
lmm2 <- lmer(log_mean_nod_biomass~Population + Strain + `Host Line` + (1|Block), data=data_sym)
#plot(lmm2)
#summary(lmm2)
Anova(lmm2, type=3)
```

    ## Analysis of Deviance Table (Type III Wald chisquare tests)
    ## 
    ## Response: log_mean_nod_biomass
    ##               Chisq Df Pr(>Chisq)    
    ## (Intercept) 188.837  1  < 2.2e-16 ***
    ## Population   31.473  5  7.554e-06 ***
    ## Strain       98.410 20  2.422e-12 ***
    ## `Host Line`  20.749  6   0.002035 ** 
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
#Extract line means and add them to the dataframe

tmp <- as.data.frame(emmeans(lmm2, "Strain"))
colnames(tmp) <- c("Strain", "Population", "nod_mass_lsmean", "nod_mass_SE", "nod_mass_df", "nod_mass_lowCL", "nod_mass_highCL")
df.means <- cbind(df.means, tmp[,3:7])

#Plant biomass

lmm3 <- lmer(plant_biomass~Population + Strain  + `Host Line` + (1|Block), data=data_sym)
#plot(lmm3)
#summary(lmm3)
Anova(lmm3, type=3)
```

    ## Analysis of Deviance Table (Type III Wald chisquare tests)
    ## 
    ## Response: plant_biomass
    ##               Chisq Df Pr(>Chisq)    
    ## (Intercept)  5.9654  1    0.01459 *  
    ## Population  31.4326  5  7.694e-06 ***
    ## Strain      71.8260 20  9.143e-08 ***
    ## `Host Line` 29.1729  6  5.642e-05 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
#Extract line means and add them to the dataframe

tmp <- as.data.frame(emmeans(lmm3, "Strain"))
colnames(tmp) <- c("Strain", "Population", "plant_biomass_lsmean", "plant_biomass_SE", "plant_biomass_df", "plant_biomass_lowCL", "plant_biomass_highCL")
df.means <- cbind(df.means, tmp[,3:7])

#Relative growth rate, as per paper
#The paper says they took the log of relative growth rate, but there are negative numbers, hmmm...

lmm4 <- lmer(log(`Relative Growth`)~Population + Strain  + `Host Line` + (1|Block), data=data_sym)
```

    ## Warning in log(`Relative Growth`): NaNs produced

    ## Warning in log(`Relative Growth`): NaNs produced

    ## Warning in log(`Relative Growth`): NaNs produced

``` r
#plot(lmm4)
#summary(lmm4)
Anova(lmm4, type=3)
```

    ## Analysis of Deviance Table (Type III Wald chisquare tests)
    ## 
    ## Response: log(`Relative Growth`)
    ##                Chisq Df Pr(>Chisq)    
    ## (Intercept) 289.3214  1  < 2.2e-16 ***
    ## Population   11.5919  5   0.040829 *  
    ## Strain       41.7485 20   0.002982 ** 
    ## `Host Line`   6.6951  6   0.349969    
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
#Extract line means and add them to the dataframe

tmp <- as.data.frame(emmeans(lmm4, "Strain"))
colnames(tmp) <- c("Strain", "Population", "RGR_lsmean", "RGR_SE", "RGR_df", "RGR_lowCL", "RGR_highCL")
df.means <- cbind(df.means, tmp[,3:7])

#Add CHR and SI values for each isolate; they are identical for all replicates of each isolate, so no need to model them

tmp <- data_sym %>% group_by(Strain) %>% summarize(CHR = mean(`CHR local abundance`), SI = mean(`SI local abunance`))
df.means <- merge(df.means, tmp, by="Strain")

#I can infer the biomass mass of plants in the controls for each strain just by simple math
#If RGR is inoculated biomass/control biomass than control biomass is just inoculated biomass/RGR
df.means$control_biomass <- df.means$plant_biomass_lsmean/df.means$RGR_lsmean

#Standardize trait and fitness values for each population separately
tmp <- df.means %>% group_by(Population) %>% summarize(
  pop_mean_tot_nod=mean(tot_nod_lsmean), 
  pop_sd_tot_nod=sd(tot_nod_lsmean),
  pop_mean_nod_mass=mean(nod_mass_lsmean),
  pop_sd_nod_mass=sd(nod_mass_lsmean),
  pop_mean_CHR=mean(CHR),
  pop_mean_SI=mean(SI),
  pop_mean_plant_biomass=mean(plant_biomass_lsmean),
  pop_sd_plant_biomass=sd(plant_biomass_lsmean),
  pop_mean_RGR=mean(RGR_lsmean),
  pop_sd_RGR=mean(RGR_lsmean)
  )
df.means <- merge(df.means, tmp, by="Population")

#Standardize traits by subtracting the population mean and dividing by the population standard deviation
df.means$plant_biomass_std <- (df.means$plant_biomass_lsmean - df.means$pop_mean_plant_biomass)/df.means$pop_sd_plant_biomass
df.means$RGR <- (df.means$RGR_lsmean - df.means$pop_mean_RGR)/df.means$pop_sd_RGR

#Relativize fitness by dividing by the mean fitness
df.means$tot_nod_std <- df.means$tot_nod_lsmean/df.means$pop_mean_tot_nod
df.means$nod_mass_std <- df.means$nod_mass_lsmean/df.means$pop_mean_nod_mass

#I did not relativize CHR and SI because I think they are already relative frequencies within each population
```

Calculate selection gradients
=============================

Next, I adopt standard genetic selection analyses to estimate selection gradients, S, by regressing relativized fitness against standardized trait values.

``` r
#Total nodules and plant biomass

lm1 <- lm(tot_nod_std~plant_biomass_std, data=df.means) #Model
summary(lm1)
```

    ## 
    ## Call:
    ## lm(formula = tot_nod_std ~ plant_biomass_std, data = df.means)
    ## 
    ## Residuals:
    ##      Min       1Q   Median       3Q      Max 
    ## -0.38710 -0.05423  0.01020  0.07300  0.17984 
    ## 
    ## Coefficients:
    ##                   Estimate Std. Error t value Pr(>|t|)    
    ## (Intercept)        1.00000    0.02202  45.403   <2e-16 ***
    ## plant_biomass_std  0.06043    0.02511   2.407   0.0242 *  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 0.1123 on 24 degrees of freedom
    ## Multiple R-squared:  0.1944, Adjusted R-squared:  0.1608 
    ## F-statistic: 5.791 on 1 and 24 DF,  p-value: 0.02416

``` r
S1 <- paste0("S = ", round(summary(lm1)$coefficients[2,1], 3))
pval1 <- round(summary(lm1)$coefficients[2,4],3)
pval1 <- ifelse(pval1 <= 0.05, paste0("p = ", pval1), "ns")
S1 <- paste0(S1, ", ", pval1)

p1 <- ggplot(data=df.means, aes(y=tot_nod_std, x=plant_biomass_std))+
      geom_point()+
      geom_smooth(method="lm", se=FALSE)+
      geom_text(aes(label=Strain),hjust=0, vjust=0, size=3)+
      ylab("Nodule number")+
      xlab("Plant biomass")+
      guides(color=FALSE)

#Nodule mass and plant biomass

lm2 <- lm(nod_mass_std~plant_biomass_std, data=df.means) #Model
summary(lm2)
```

    ## 
    ## Call:
    ## lm(formula = nod_mass_std ~ plant_biomass_std, data = df.means)
    ## 
    ## Residuals:
    ##      Min       1Q   Median       3Q      Max 
    ## -0.29520 -0.06817  0.04082  0.07238  0.39187 
    ## 
    ## Coefficients:
    ##                   Estimate Std. Error t value Pr(>|t|)    
    ## (Intercept)        1.00000    0.03281  30.483   <2e-16 ***
    ## plant_biomass_std  0.03226    0.03740   0.862    0.397    
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 0.1673 on 24 degrees of freedom
    ## Multiple R-squared:  0.03006,    Adjusted R-squared:  -0.01035 
    ## F-statistic: 0.7439 on 1 and 24 DF,  p-value: 0.397

``` r
S2 <- paste0("S = ", round(summary(lm2)$coefficients[2,1], 3))
pval2 <- round(summary(lm2)$coefficients[2,4],3)
pval2 <- ifelse(pval2 <= 0.05, paste0("p = ", pval2), "ns")
S2 <- paste0(S2, ", ", pval2)

p2 <- ggplot(data=df.means, aes(y=nod_mass_std, x=plant_biomass_std))+
      geom_point()+
      geom_smooth(method="lm", se=FALSE, linetype="dashed")+
      geom_text(aes(label=Strain),hjust=0, vjust=0, size=3)+
      ylab("Nodule mass")+
      xlab("Plant biomass")+
      guides(color=FALSE)

#CHR and plant biomass

lm3 <- lm(CHR~plant_biomass_std, data=df.means) #Model
summary(lm3)
```

    ## 
    ## Call:
    ## lm(formula = CHR ~ plant_biomass_std, data = df.means)
    ## 
    ## Residuals:
    ##      Min       1Q   Median       3Q      Max 
    ## -0.16589 -0.07438 -0.03063  0.00495  0.49706 
    ## 
    ## Coefficients:
    ##                   Estimate Std. Error t value Pr(>|t|)   
    ## (Intercept)        0.09249    0.02920   3.167  0.00416 **
    ## plant_biomass_std -0.05287    0.03329  -1.588  0.12539   
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 0.1489 on 24 degrees of freedom
    ## Multiple R-squared:  0.09508,    Adjusted R-squared:  0.05737 
    ## F-statistic: 2.522 on 1 and 24 DF,  p-value: 0.1254

``` r
S3 <- paste0("S = ", round(summary(lm3)$coefficients[2,1], 3))
pval3 <- round(summary(lm3)$coefficients[2,4],3)
pval3 <- ifelse(pval3 <= 0.05, paste0("p = ", pval3), "ns")
S3 <- paste0(S3, ", ", pval3)

p3 <- ggplot(data=df.means, aes(y=CHR, x=plant_biomass_std))+
      geom_point()+
      geom_smooth(method="lm", se=FALSE, linetype="dashed")+
      geom_text(aes(label=Strain),hjust=0, vjust=0, size=3)+
      ylab("CHR")+
      xlab("Plant biomass")+
      guides(color=FALSE)

#SI and plant biomass

lm4 <- lm(SI~plant_biomass_std, data=df.means) #Model
summary(lm4)
```

    ## 
    ## Call:
    ## lm(formula = SI ~ plant_biomass_std, data = df.means)
    ## 
    ## Residuals:
    ##      Min       1Q   Median       3Q      Max 
    ## -0.11700 -0.08748 -0.05674  0.11802  0.23694 
    ## 
    ## Coefficients:
    ##                     Estimate Std. Error t value Pr(>|t|)    
    ## (Intercept)        0.1313591  0.0220007   5.971 4.35e-06 ***
    ## plant_biomass_std -0.0007803  0.0248091  -0.031    0.975    
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 0.11 on 23 degrees of freedom
    ##   (1 observation deleted due to missingness)
    ## Multiple R-squared:  4.301e-05,  Adjusted R-squared:  -0.04343 
    ## F-statistic: 0.0009893 on 1 and 23 DF,  p-value: 0.9752

``` r
S4 <- paste0("S = ", round(summary(lm4)$coefficients[2,1], 3))
pval4 <- round(summary(lm4)$coefficients[2,4],3)
pval4 <- ifelse(pval4 <= 0.05, paste0("p = ", pval4), "ns")
S4 <- paste0(S4, ", ", pval4)

p4 <- ggplot(data=df.means, aes(y=SI, x=plant_biomass_std))+
      geom_point()+
      geom_smooth(method="lm", se=FALSE, linetype="dashed")+
      geom_text(aes(label=Strain),hjust=0, vjust=0, size=3)+
      ylab("SI")+
      xlab("Plant biomass")+
      guides(color=FALSE)

#Total nodules and RGR

lm5 <- lm(tot_nod_std~RGR, data=df.means) #Model
summary(lm5)
```

    ## 
    ## Call:
    ## lm(formula = tot_nod_std ~ RGR, data = df.means)
    ## 
    ## Residuals:
    ##      Min       1Q   Median       3Q      Max 
    ## -0.48293 -0.03392  0.00387  0.05128  0.20747 
    ## 
    ## Coefficients:
    ##             Estimate Std. Error t value Pr(>|t|)    
    ## (Intercept)  1.00000    0.02452  40.781   <2e-16 ***
    ## RGR         -0.06671    0.36065  -0.185    0.855    
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 0.125 on 24 degrees of freedom
    ## Multiple R-squared:  0.001423,   Adjusted R-squared:  -0.04018 
    ## F-statistic: 0.03421 on 1 and 24 DF,  p-value: 0.8548

``` r
S5 <- paste0("S = ", round(summary(lm5)$coefficients[2,1], 3))
pval5 <- round(summary(lm5)$coefficients[2,4],3)
pval5 <- ifelse(pval5 <= 0.05, paste0("p = ", pval5), "ns")
S5 <- paste0(S5, ", ", pval5)

p5 <- ggplot(data=df.means, aes(y=tot_nod_std, x=RGR))+
      geom_point()+
      geom_smooth(method="lm", se=FALSE, linetype="dashed")+
      geom_text(aes(label=Strain),hjust=0, vjust=0, size=3)+
      ylab("Nodule number")+
      xlab("RGR")+
      guides(color=FALSE)

#Nodule mass and RGR

lm6 <- lm(nod_mass_std~RGR, data=df.means) #Model
summary(lm6)
```

    ## 
    ## Call:
    ## lm(formula = nod_mass_std ~ RGR, data = df.means)
    ## 
    ## Residuals:
    ##      Min       1Q   Median       3Q      Max 
    ## -0.30834 -0.06279  0.01340  0.09198  0.35255 
    ## 
    ## Coefficients:
    ##             Estimate Std. Error t value Pr(>|t|)    
    ## (Intercept)  1.00000    0.03330  30.028   <2e-16 ***
    ## RGR          0.05256    0.48980   0.107    0.915    
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 0.1698 on 24 degrees of freedom
    ## Multiple R-squared:  0.0004795,  Adjusted R-squared:  -0.04117 
    ## F-statistic: 0.01151 on 1 and 24 DF,  p-value: 0.9154

``` r
S6 <- paste0("S = ", round(summary(lm6)$coefficients[2,1], 3))
pval6 <- round(summary(lm6)$coefficients[2,4],3)
pval6 <- ifelse(pval6 <= 0.05, paste0("p = ", pval6), "ns")
S6 <- paste0(S6, ", ", pval6)

p6 <- ggplot(data=df.means, aes(y=nod_mass_std, x=RGR))+
      geom_point()+
      geom_smooth(method="lm", se=FALSE, linetype="dashed")+
      geom_text(aes(label=Strain),hjust=0, vjust=0, size=3)+
      ylab("Nodule mass")+
      xlab("RGR")+
      guides(color=FALSE)

#CHR and RGR

lm7 <- lm(CHR~RGR, data=df.means) #Model
summary(lm7)
```

    ## 
    ## Call:
    ## lm(formula = CHR ~ RGR, data = df.means)
    ## 
    ## Residuals:
    ##      Min       1Q   Median       3Q      Max 
    ## -0.13020 -0.07043 -0.05992  0.00096  0.48845 
    ## 
    ## Coefficients:
    ##             Estimate Std. Error t value Pr(>|t|)   
    ## (Intercept)  0.09249    0.02916   3.171  0.00412 **
    ## RGR         -0.68991    0.42894  -1.608  0.12082   
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 0.1487 on 24 degrees of freedom
    ## Multiple R-squared:  0.0973, Adjusted R-squared:  0.05969 
    ## F-statistic: 2.587 on 1 and 24 DF,  p-value: 0.1208

``` r
S7 <- paste0("S = ", round(summary(lm7)$coefficients[2,1], 3))
pval7 <- round(summary(lm7)$coefficients[2,4],3)
pval7 <- ifelse(pval7 <= 0.05, paste0("p = ", pval7), "ns")
S7 <- paste0(S7, ", ", pval7)

p7 <- ggplot(data=df.means, aes(y=CHR, x=RGR))+
      geom_point()+
      geom_smooth(method="lm", se=FALSE, linetype="dashed")+
      geom_text(aes(label=Strain),hjust=0, vjust=0, size=3)+
      ylab("CHR")+
      xlab("RGR")

#SI and RGR

lm8 <- lm(SI~RGR, data=df.means) #Model
summary(lm8)
```

    ## 
    ## Call:
    ## lm(formula = SI ~ RGR, data = df.means)
    ## 
    ## Residuals:
    ##      Min       1Q   Median       3Q      Max 
    ## -0.12043 -0.09170 -0.04298  0.11322  0.24101 
    ## 
    ## Coefficients:
    ##             Estimate Std. Error t value Pr(>|t|)    
    ## (Intercept)  0.13136    0.02184   6.016 3.91e-06 ***
    ## RGR         -0.18190    0.31493  -0.578    0.569    
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 0.1092 on 23 degrees of freedom
    ##   (1 observation deleted due to missingness)
    ## Multiple R-squared:  0.0143, Adjusted R-squared:  -0.02856 
    ## F-statistic: 0.3336 on 1 and 23 DF,  p-value: 0.5691

``` r
S8 <- paste0("S = ", round(summary(lm8)$coefficients[2,1], 3))
pval8 <- round(summary(lm8)$coefficients[2,4],3)
pval8 <- ifelse(pval8 <= 0.05, paste0("p = ", pval8), "ns")
S8 <- paste0(S8, ", ", pval8)

p8 <- ggplot(data=df.means, aes(y=SI, x=RGR))+
      geom_point()+
      geom_smooth(method="lm", se=FALSE, linetype="dashed")+
      geom_text(aes(label=Strain),hjust=0, vjust=0, size=3)+
      ylab("SI")+
      xlab("RGR")

plot_grid(p1,p2,p3,p4,p5,p6,p7,p8, nrow=2, ncol=4, labels=c(S1, S2, S3, S4, S5, S6, S7, S8), scale = 0.9)
```

    ## Warning: Removed 1 rows containing non-finite values (stat_smooth).

    ## Warning: Removed 1 rows containing missing values (geom_point).

    ## Warning: Removed 1 rows containing missing values (geom_text).

    ## Warning: Removed 1 rows containing non-finite values (stat_smooth).

    ## Warning: Removed 1 rows containing missing values (geom_point).

    ## Warning: Removed 1 rows containing missing values (geom_text).

![](README_files/figure-markdown_github/Selection%20gradients-1.png)

Why do different fitness measures give different answers?
=========================================================

``` r
#Correlations among rhizobium fitness proxies
lm9 <- lm(SI~CHR, data=df.means) #Model
summary(lm9)
```

    ## 
    ## Call:
    ## lm(formula = SI ~ CHR, data = df.means)
    ## 
    ## Residuals:
    ##      Min       1Q   Median       3Q      Max 
    ## -0.13783 -0.08810 -0.04400  0.07736  0.23980 
    ## 
    ## Coefficients:
    ##             Estimate Std. Error t value Pr(>|t|)    
    ## (Intercept)  0.11498    0.02516   4.570 0.000136 ***
    ## CHR          0.17131    0.13977   1.226 0.232722    
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 0.1065 on 23 degrees of freedom
    ##   (1 observation deleted due to missingness)
    ## Multiple R-squared:  0.06131,    Adjusted R-squared:  0.0205 
    ## F-statistic: 1.502 on 1 and 23 DF,  p-value: 0.2327

``` r
p9 <- ggplot(data=df.means, aes(y=SI, x=CHR))+
      geom_point()+
      geom_smooth(method="lm", se=FALSE)+
      xlab("CHR")+
      ylab("SI")+
      geom_text(aes(label=Strain),hjust=0, vjust=0, size=3)

lm10 <- lm(SI~nod_mass_std, data=df.means) #Model
summary(lm10)
```

    ## 
    ## Call:
    ## lm(formula = SI ~ nod_mass_std, data = df.means)
    ## 
    ## Residuals:
    ##      Min       1Q   Median       3Q      Max 
    ## -0.11942 -0.09091 -0.05124  0.11540  0.24595 
    ## 
    ## Coefficients:
    ##              Estimate Std. Error t value Pr(>|t|)
    ## (Intercept)   0.19860    0.13334   1.489    0.150
    ## nod_mass_std -0.06718    0.13144  -0.511    0.614
    ## 
    ## Residual standard error: 0.1093 on 23 degrees of freedom
    ##   (1 observation deleted due to missingness)
    ## Multiple R-squared:  0.01123,    Adjusted R-squared:  -0.03176 
    ## F-statistic: 0.2612 on 1 and 23 DF,  p-value: 0.6142

``` r
p10 <- ggplot(data=df.means, aes(y=SI, x=nod_mass_std))+
      geom_point()+
      geom_smooth(method="lm", se=FALSE)+
      xlab("Nodule mass")+
      ylab("SI")+
      geom_text(aes(label=Strain),hjust=0, vjust=0, size=3)

lm11 <- lm(SI~tot_nod_std, data=df.means) #Model
summary(lm11)
```

    ## 
    ## Call:
    ## lm(formula = SI ~ tot_nod_std, data = df.means)
    ## 
    ## Residuals:
    ##      Min       1Q   Median       3Q      Max 
    ## -0.12527 -0.09481 -0.02060  0.10651  0.22264 
    ## 
    ## Coefficients:
    ##             Estimate Std. Error t value Pr(>|t|)
    ## (Intercept)  0.01055    0.18176   0.058    0.954
    ## tot_nod_std  0.12141    0.18132   0.670    0.510
    ## 
    ## Residual standard error: 0.1089 on 23 degrees of freedom
    ##   (1 observation deleted due to missingness)
    ## Multiple R-squared:  0.01912,    Adjusted R-squared:  -0.02353 
    ## F-statistic: 0.4483 on 1 and 23 DF,  p-value: 0.5098

``` r
p11 <- ggplot(data=df.means, aes(y=SI, x=tot_nod_std))+
      geom_point()+
      geom_smooth(method="lm", se=FALSE)+
      xlab("Nodule number")+
      ylab("SI")+
      geom_text(aes(label=Strain),hjust=0, vjust=0, size=3)

lm12 <- lm(CHR~nod_mass_std, data=df.means) #Model
summary(lm12)
```

    ## 
    ## Call:
    ## lm(formula = CHR ~ nod_mass_std, data = df.means)
    ## 
    ## Residuals:
    ##      Min       1Q   Median       3Q      Max 
    ## -0.15607 -0.08306 -0.04832  0.01887  0.44022 
    ## 
    ## Coefficients:
    ##              Estimate Std. Error t value Pr(>|t|)  
    ## (Intercept)    0.3584     0.1825   1.964   0.0612 .
    ## nod_mass_std  -0.2659     0.1801  -1.476   0.1529  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 0.1499 on 24 degrees of freedom
    ## Multiple R-squared:  0.08324,    Adjusted R-squared:  0.04504 
    ## F-statistic: 2.179 on 1 and 24 DF,  p-value: 0.1529

``` r
p12 <- ggplot(data=df.means, aes(y=CHR, x=nod_mass_std))+
      geom_point()+
      geom_smooth(method="lm", se=FALSE)+
      xlab("Nodule mass")+
      ylab("CHR")+
      geom_text(aes(label=Strain),hjust=0, vjust=0, size=3)

lm13 <- lm(CHR~tot_nod_std, data=df.means) #Model
summary(lm13)
```

    ## 
    ## Call:
    ## lm(formula = CHR ~ tot_nod_std, data = df.means)
    ## 
    ## Residuals:
    ##      Min       1Q   Median       3Q      Max 
    ## -0.08636 -0.07815 -0.06983  0.02307  0.52308 
    ## 
    ## Coefficients:
    ##             Estimate Std. Error t value Pr(>|t|)
    ## (Intercept)  0.08047    0.25717   0.313    0.757
    ## tot_nod_std  0.01202    0.25533   0.047    0.963
    ## 
    ## Residual standard error: 0.1565 on 24 degrees of freedom
    ## Multiple R-squared:  9.239e-05,  Adjusted R-squared:  -0.04157 
    ## F-statistic: 0.002218 on 1 and 24 DF,  p-value: 0.9628

``` r
p13 <- ggplot(data=df.means, aes(y=CHR, x=tot_nod_std))+
      geom_point()+
      geom_smooth(method="lm", se=FALSE)+
      xlab("Nodule number")+
      ylab("CHR")+
      geom_text(aes(label=Strain),hjust=0, vjust=0, size=3)

lm14 <- lm(tot_nod_std~nod_mass_std, data=df.means) #Model
summary(lm14)
```

    ## 
    ## Call:
    ## lm(formula = tot_nod_std ~ nod_mass_std, data = df.means)
    ## 
    ## Residuals:
    ##      Min       1Q   Median       3Q      Max 
    ## -0.35809 -0.03602  0.00038  0.05111  0.16656 
    ## 
    ## Coefficients:
    ##              Estimate Std. Error t value Pr(>|t|)    
    ## (Intercept)    0.5533     0.1212   4.567 0.000125 ***
    ## nod_mass_std   0.4467     0.1196   3.735 0.001026 ** 
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 0.0995 on 24 degrees of freedom
    ## Multiple R-squared:  0.3676, Adjusted R-squared:  0.3413 
    ## F-statistic: 13.95 on 1 and 24 DF,  p-value: 0.001026

``` r
p14 <- ggplot(data=df.means, aes(y=tot_nod_std, x=nod_mass_std))+
      geom_point()+
      geom_smooth(method="lm", se=FALSE)+
      xlab("Nodule mass")+
      ylab("Nodule number")+
      geom_text(aes(label=Strain),hjust=0, vjust=0, size=3)

#Corrlation between plant fitness proxies

lm15 <- lm(RGR~plant_biomass_std, data=df.means)
summary(lm15)
```

    ## 
    ## Call:
    ## lm(formula = RGR ~ plant_biomass_std, data = df.means)
    ## 
    ## Residuals:
    ##       Min        1Q    Median        3Q       Max 
    ## -0.173456 -0.019751  0.001007  0.013759  0.121792 
    ## 
    ## Coefficients:
    ##                    Estimate Std. Error t value Pr(>|t|)    
    ## (Intercept)       1.297e-17  1.070e-02   0.000  1.00000    
    ## plant_biomass_std 4.936e-02  1.220e-02   4.045  0.00047 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 0.05457 on 24 degrees of freedom
    ## Multiple R-squared:  0.4054, Adjusted R-squared:  0.3806 
    ## F-statistic: 16.36 on 1 and 24 DF,  p-value: 0.0004699

``` r
p15 <- ggplot(data=df.means, aes(y=RGR, x=plant_biomass_std))+
      geom_point()+
      geom_smooth(method="lm", se=FALSE)+
      geom_text(aes(label=Strain),hjust=0, vjust=0, size=3)+
      xlab("Plant biomass")+
      ylab("RGR")

plot_grid(p9, p10, p11, p12, p13, p14, p15)
```

    ## Warning: Removed 1 rows containing non-finite values (stat_smooth).

    ## Warning: Removed 1 rows containing missing values (geom_point).

    ## Warning: Removed 1 rows containing missing values (geom_text).

    ## Warning: Removed 1 rows containing non-finite values (stat_smooth).

    ## Warning: Removed 1 rows containing missing values (geom_point).

    ## Warning: Removed 1 rows containing missing values (geom_text).

    ## Warning: Removed 1 rows containing non-finite values (stat_smooth).

    ## Warning: Removed 1 rows containing missing values (geom_point).

    ## Warning: Removed 1 rows containing missing values (geom_text).

![](README_files/figure-markdown_github/Correlations%20among%20fitness%20proxies,%20-1.png)
