Much ado about fitness correlations
================
Megan Frederickson
2019-04-17

Fitness conflict or fitness alignment?
--------------------------------------

This repository re-analyzes the data in: Gano-Cohen KA, Wendlandt CE, Stokes PJ, Blanton MA, Quides KW, Zomorrodian A, Adinata ES, Sachs JL (2019) Interspecific conflict and the evolution of ineffective rhizobia. Ecology Letters. <https://doi.org/10.1111/ele.13247>

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

Calculate genotype (i.e., bacterial strain) means for each plant and rhizobium fitness measure
==============================================================================================

The paper has several measures of plant and rhizobium fitness. First, we need to calculate the appropriate strain means for each variable. I did this using the emmeans package, which calculates estimated marginal means from a linear mixed model. I have tried to stick as close as possible to the analysis presented in the paper, which states that "all effects were coded as fixed except block, which was treated as a random effect" and that variables were log-transformed as needed to improve normality. I also calculated strain means only on data from sympatric hosts, again following the paper's general modelling approach.

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
```

    ## fixed-effect model matrix is rank deficient so dropping 10 columns / coefficients

``` r
plot(lmm1)
```

![](README_files/figure-markdown_github/Genetic%20values-1.png)

``` r
summary(lmm1)
```

    ## Linear mixed model fit by REML ['lmerMod']
    ## Formula: log_nodules ~ Population + Strain + `Host Line` + (1 | Block)
    ##    Data: data_sym
    ## 
    ## REML criterion at convergence: 411.7
    ## 
    ## Scaled residuals: 
    ##     Min      1Q  Median      3Q     Max 
    ## -3.3368 -0.5193  0.1419  0.5232  2.5940 
    ## 
    ## Random effects:
    ##  Groups   Name        Variance Std.Dev.
    ##  Block    (Intercept) 0.01402  0.1184  
    ##  Residual             0.27504  0.5244  
    ## Number of obs: 248, groups:  Block, 5
    ## 
    ## Fixed effects:
    ##                     Estimate Std. Error t value
    ## (Intercept)          3.24562    0.20769  15.627
    ## PopulationBMR        0.55387    0.27334   2.026
    ## PopulationCLA        0.48821    0.27081   1.803
    ## PopulationGRI       -0.08345    0.28323  -0.295
    ## PopulationUCR        0.35194    0.27162   1.296
    ## PopulationYUC        0.81103    0.27334   2.967
    ## Strain134           -0.28279    0.23454  -1.206
    ## Strain135            0.17050    0.23454   0.727
    ## Strain136           -0.57888    0.23454  -2.468
    ## Strain137            0.02122    0.24723   0.086
    ## Strain138            0.43402    0.24723   1.756
    ## Strain139            0.00619    0.24723   0.025
    ## Strain142           -0.15110    0.23454  -0.644
    ## Strain143           -0.14475    0.23454  -0.617
    ## Strain144            0.28410    0.23454   1.211
    ## Strain145            0.48398    0.23454   2.064
    ## Strain147           -0.23910    0.23454  -1.019
    ## Strain149           -1.92890    0.25888  -7.451
    ## Strain150           -0.01976    0.23454  -0.084
    ## Strain151            0.27323    0.24122   1.133
    ## Strain153           -0.10837    0.23454  -0.462
    ## Strain154            0.12528    0.23454   0.534
    ## Strain155            0.22920    0.23454   0.977
    ## Strain157           -0.09869    0.24723  -0.399
    ## Strain158            0.25307    0.24723   1.024
    ## Strain159            0.23547    0.24723   0.952
    ## `Host Line`Anz11.01  0.69719    0.17654   3.949
    ## `Host Line`BMR01.03 -0.19616    0.16584  -1.183
    ## `Host Line`Cla01.04  0.03900    0.14834   0.263
    ## `Host Line`Gri01.01  0.15723    0.17654   0.891
    ## `Host Line`UCR02.07 -0.04094    0.15523  -0.264
    ## `Host Line`Yuc02.01 -0.12039    0.16584  -0.726

    ## 
    ## Correlation matrix not shown by default, as p = 32 > 12.
    ## Use print(x, correlation=TRUE)  or
    ##     vcov(x)        if you need it

    ## fit warnings:
    ## fixed-effect model matrix is rank deficient so dropping 10 columns / coefficients

``` r
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
```

    ## NOTE: A nesting structure was detected in the fitted model:
    ##     Strain %in% Population, Host Line %in% Population

``` r
colnames(df.means) <- c("Strain", "Population", "tot_nod_lsmean", "tot_nod_SE", "tot_nod_df", "tot_nod_lowCL", "tot_nod_highCL")

#Mean nodule mass
lmm2 <- lmer(log_mean_nod_biomass~Population + Strain + `Host Line` + (1|Block), data=data_sym)
```

    ## fixed-effect model matrix is rank deficient so dropping 10 columns / coefficients

``` r
plot(lmm2)
```

![](README_files/figure-markdown_github/Genetic%20values-2.png)

``` r
summary(lmm2)
```

    ## Linear mixed model fit by REML ['lmerMod']
    ## Formula: log_mean_nod_biomass ~ Population + Strain + `Host Line` + (1 |  
    ##     Block)
    ##    Data: data_sym
    ## 
    ## REML criterion at convergence: 330.9
    ## 
    ## Scaled residuals: 
    ##     Min      1Q  Median      3Q     Max 
    ## -4.6452 -0.5094  0.0064  0.4812  2.8175 
    ## 
    ## Random effects:
    ##  Groups   Name        Variance  Std.Dev.
    ##  Block    (Intercept) 0.0009515 0.03085 
    ##  Residual             0.1965576 0.44335 
    ## Number of obs: 245, groups:  Block, 5
    ## 
    ## Fixed effects:
    ##                      Estimate Std. Error t value
    ## (Intercept)         -2.335885   0.169984 -13.742
    ## PopulationBMR        0.858279   0.230811   3.719
    ## PopulationCLA        0.916141   0.228672   4.006
    ## PopulationGRI        0.797451   0.239435   3.331
    ## PopulationUCR        0.314670   0.230257   1.367
    ## PopulationYUC        1.055891   0.230811   4.575
    ## Strain134            0.430224   0.198271   2.170
    ## Strain135           -0.321738   0.198271  -1.623
    ## Strain136            0.381695   0.198271   1.925
    ## Strain137            0.078228   0.208996   0.374
    ## Strain138            0.129668   0.208996   0.620
    ## Strain139            0.163475   0.208996   0.782
    ## Strain142            0.170838   0.198271   0.862
    ## Strain143           -0.248303   0.198271  -1.252
    ## Strain144           -0.474541   0.198271  -2.393
    ## Strain145           -0.925369   0.198271  -4.667
    ## Strain147           -0.175081   0.203866  -0.859
    ## Strain149            0.490680   0.218743   2.243
    ## Strain150           -0.277968   0.211028  -1.317
    ## Strain151           -0.191761   0.203866  -0.941
    ## Strain153           -0.797237   0.198271  -4.021
    ## Strain154           -0.776597   0.198271  -3.917
    ## Strain155           -0.781678   0.198271  -3.942
    ## Strain157           -0.075734   0.208996  -0.362
    ## Strain158           -0.185174   0.208996  -0.886
    ## Strain159           -0.364380   0.208996  -1.743
    ## `Host Line`Anz11.01  0.461724   0.148836   3.102
    ## `Host Line`BMR01.03 -0.001436   0.140199  -0.010
    ## `Host Line`Cla01.04 -0.169875   0.125398  -1.355
    ## `Host Line`Gri01.01 -0.096271   0.148836  -0.647
    ## `Host Line`UCR02.07  0.404081   0.136643   2.957
    ## `Host Line`Yuc02.01 -0.049732   0.140199  -0.355

    ## 
    ## Correlation matrix not shown by default, as p = 32 > 12.
    ## Use print(x, correlation=TRUE)  or
    ##     vcov(x)        if you need it

    ## fit warnings:
    ## fixed-effect model matrix is rank deficient so dropping 10 columns / coefficients

``` r
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
```

    ## NOTE: A nesting structure was detected in the fitted model:
    ##     Strain %in% Population, Host Line %in% Population

``` r
colnames(tmp) <- c("Strain", "Population", "nod_mass_lsmean", "nod_mass_SE", "nod_mass_df", "nod_mass_lowCL", "nod_mass_highCL")
df.means <- cbind(df.means, tmp[,3:7])

#Plant biomass

lmm3 <- lmer(plant_biomass~Population + Strain  + `Host Line` + (1|Block), data=data_sym)
```

    ## fixed-effect model matrix is rank deficient so dropping 10 columns / coefficients

``` r
plot(lmm3)
```

![](README_files/figure-markdown_github/Genetic%20values-3.png)

``` r
summary(lmm3)
```

    ## Linear mixed model fit by REML ['lmerMod']
    ## Formula: plant_biomass ~ Population + Strain + `Host Line` + (1 | Block)
    ##    Data: data_sym
    ## 
    ## REML criterion at convergence: -641.1
    ## 
    ## Scaled residuals: 
    ##     Min      1Q  Median      3Q     Max 
    ## -3.5131 -0.5915 -0.0139  0.5272  2.9215 
    ## 
    ## Random effects:
    ##  Groups   Name        Variance  Std.Dev.
    ##  Block    (Intercept) 0.0002556 0.01599 
    ##  Residual             0.0020753 0.04556 
    ## Number of obs: 248, groups:  Block, 5
    ## 
    ## Fixed effects:
    ##                      Estimate Std. Error t value
    ## (Intercept)          0.046073   0.018864   2.442
    ## PopulationBMR        0.058607   0.023751   2.468
    ## PopulationCLA        0.088591   0.023532   3.765
    ## PopulationGRI        0.113654   0.024603   4.620
    ## PopulationUCR        0.036486   0.023599   1.546
    ## PopulationYUC        0.021975   0.023751   0.925
    ## Strain134            0.039550   0.020373   1.941
    ## Strain135            0.045760   0.020373   2.246
    ## Strain136            0.005550   0.020373   0.272
    ## Strain137           -0.025067   0.021475  -1.167
    ## Strain138            0.035667   0.021475   1.661
    ## Strain139           -0.021844   0.021475  -1.017
    ## Strain142           -0.001810   0.020373  -0.089
    ## Strain143           -0.018120   0.020373  -0.889
    ## Strain144            0.036380   0.020373   1.786
    ## Strain145           -0.043070   0.020373  -2.114
    ## Strain147           -0.005990   0.020373  -0.294
    ## Strain149           -0.028026   0.022491  -1.246
    ## Strain150           -0.009070   0.020373  -0.445
    ## Strain151            0.005730   0.020955   0.273
    ## Strain153            0.038720   0.020373   1.901
    ## Strain154            0.082290   0.020373   4.039
    ## Strain155           -0.027600   0.020373  -1.355
    ## Strain157            0.020511   0.021475   0.955
    ## Strain158           -0.004633   0.021475  -0.216
    ## Strain159            0.018644   0.021475   0.868
    ## `Host Line`Anz11.01  0.068752   0.015347   4.480
    ## `Host Line`BMR01.03 -0.022440   0.014406  -1.558
    ## `Host Line`Cla01.04 -0.008648   0.012885  -0.671
    ## `Host Line`Gri01.01  0.025074   0.015347   1.634
    ## `Host Line`UCR02.07  0.015502   0.013487   1.149
    ## `Host Line`Yuc02.01 -0.022395   0.014406  -1.555

    ## 
    ## Correlation matrix not shown by default, as p = 32 > 12.
    ## Use print(x, correlation=TRUE)  or
    ##     vcov(x)        if you need it

    ## fit warnings:
    ## fixed-effect model matrix is rank deficient so dropping 10 columns / coefficients

``` r
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
```

    ## NOTE: A nesting structure was detected in the fitted model:
    ##     Strain %in% Population, Host Line %in% Population

``` r
colnames(tmp) <- c("Strain", "Population", "plant_biomass_lsmean", "plant_biomass_SE", "plant_biomass_df", "plant_biomass_lowCL", "plant_biomass_highCL")
df.means <- cbind(df.means, tmp[,3:7])

#Relative growth rate, as per paper
#The paper says they took the log of relative growth rate, but there are negative numbers, hmmm...

hist(data_sym$`Relative Growth`) 
```

![](README_files/figure-markdown_github/Genetic%20values-4.png)

``` r
lmm4 <- lmer(log(`Relative Growth`)~Population + Strain  + `Host Line` + (1|Block), data=data_sym)
```

    ## Warning in log(`Relative Growth`): NaNs produced

    ## Warning in log(`Relative Growth`): NaNs produced

    ## Warning in log(`Relative Growth`): NaNs produced

    ## fixed-effect model matrix is rank deficient so dropping 10 columns / coefficients

``` r
plot(lmm4)
```

![](README_files/figure-markdown_github/Genetic%20values-5.png)

``` r
summary(lmm4)
```

    ## Linear mixed model fit by REML ['lmerMod']
    ## Formula: log(`Relative Growth`) ~ Population + Strain + `Host Line` +  
    ##     (1 | Block)
    ##    Data: data_sym
    ## 
    ## REML criterion at convergence: 633.2
    ## 
    ## Scaled residuals: 
    ##     Min      1Q  Median      3Q     Max 
    ## -5.7991 -0.4184  0.0753  0.5820  1.9702 
    ## 
    ## Random effects:
    ##  Groups   Name        Variance Std.Dev.
    ##  Block    (Intercept) 0.02286  0.1512  
    ##  Residual             0.84703  0.9203  
    ## Number of obs: 241, groups:  Block, 5
    ## 
    ## Fixed effects:
    ##                     Estimate Std. Error t value
    ## (Intercept)          6.62572    0.38953  17.009
    ## PopulationBMR        0.02117    0.50303   0.042
    ## PopulationCLA        0.19531    0.49881   0.392
    ## PopulationGRI       -0.21875    0.53934  -0.406
    ## PopulationUCR       -0.48111    0.50083  -0.961
    ## PopulationYUC       -1.18409    0.50511  -2.344
    ## Strain134            0.43257    0.41159   1.051
    ## Strain135            0.45630    0.41159   1.109
    ## Strain136            0.09052    0.41159   0.220
    ## Strain137           -0.35229    0.44788  -0.787
    ## Strain138            0.20094    0.44788   0.449
    ## Strain139           -0.25812    0.44788  -0.576
    ## Strain142            0.02030    0.41159   0.049
    ## Strain143           -0.19541    0.41159  -0.475
    ## Strain144            0.31345    0.41159   0.762
    ## Strain145           -0.39913    0.41159  -0.970
    ## Strain147           -0.44873    0.41159  -1.090
    ## Strain149            0.08148    0.45425   0.179
    ## Strain150            0.11685    0.42328   0.276
    ## Strain151            0.17119    0.42329   0.404
    ## Strain153            0.78806    0.42334   1.862
    ## Strain154            1.27938    0.41159   3.108
    ## Strain155           -0.96894    0.43845  -2.210
    ## Strain157            0.37614    0.46017   0.817
    ## Strain158           -0.18542    0.44792  -0.414
    ## Strain159            0.18862    0.44792   0.421
    ## `Host Line`Anz11.01 -0.29739    0.32262  -0.922
    ## `Host Line`BMR01.03  0.23868    0.29104   0.820
    ## `Host Line`Cla01.04  0.26527    0.26031   1.019
    ## `Host Line`Gri01.01  0.44888    0.31583   1.421
    ## `Host Line`UCR02.07 -0.08661    0.27589  -0.314
    ## `Host Line`Yuc02.01 -0.43158    0.30522  -1.414

    ## 
    ## Correlation matrix not shown by default, as p = 32 > 12.
    ## Use print(x, correlation=TRUE)  or
    ##     vcov(x)        if you need it

    ## fit warnings:
    ## fixed-effect model matrix is rank deficient so dropping 10 columns / coefficients

``` r
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
```

    ## NOTE: A nesting structure was detected in the fitted model:
    ##     Strain %in% Population, Host Line %in% Population

``` r
colnames(tmp) <- c("Strain", "Population", "RGR_lsmean", "RGR_SE", "RGR_df", "RGR_lowCL", "RGR_highCL")
df.means <- cbind(df.means, tmp[,3:7])

#Add CHR and SI values for each isolate; they are identical for each isolate, so no need to model them

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

#Standardize traits
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
p1
```

![](README_files/figure-markdown_github/Selection%20gradients-1.png)

``` r
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
p2
```

![](README_files/figure-markdown_github/Selection%20gradients-2.png)

``` r
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
p3
```

![](README_files/figure-markdown_github/Selection%20gradients-3.png)

``` r
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
p4
```

    ## Warning: Removed 1 rows containing non-finite values (stat_smooth).

    ## Warning: Removed 1 rows containing missing values (geom_point).

    ## Warning: Removed 1 rows containing missing values (geom_text).

![](README_files/figure-markdown_github/Selection%20gradients-4.png)

``` r
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
p5
```

![](README_files/figure-markdown_github/Selection%20gradients-5.png)

``` r
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
p6
```

![](README_files/figure-markdown_github/Selection%20gradients-6.png)

``` r
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
p7
```

![](README_files/figure-markdown_github/Selection%20gradients-7.png)

``` r
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
p8
```

    ## Warning: Removed 1 rows containing non-finite values (stat_smooth).

    ## Warning: Removed 1 rows containing missing values (geom_point).

    ## Warning: Removed 1 rows containing missing values (geom_text).

![](README_files/figure-markdown_github/Selection%20gradients-8.png)

``` r
plot_grid(p1,p2,p3,p4,p5,p6,p7,p8, nrow=2, ncol=4, labels=c(S1, S2, S3, S4, S5, S6, S7, S8), scale = 0.9)
```

    ## Warning: Removed 1 rows containing non-finite values (stat_smooth).

    ## Warning: Removed 1 rows containing missing values (geom_point).

    ## Warning: Removed 1 rows containing missing values (geom_text).

    ## Warning: Removed 1 rows containing non-finite values (stat_smooth).

    ## Warning: Removed 1 rows containing missing values (geom_point).

    ## Warning: Removed 1 rows containing missing values (geom_text).

![](README_files/figure-markdown_github/Selection%20gradients-9.png)

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
p9
```

    ## Warning: Removed 1 rows containing non-finite values (stat_smooth).

    ## Warning: Removed 1 rows containing missing values (geom_point).

    ## Warning: Removed 1 rows containing missing values (geom_text).

![](README_files/figure-markdown_github/Correlations%20among%20fitness%20proxies-1.png)

``` r
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
p10
```

    ## Warning: Removed 1 rows containing non-finite values (stat_smooth).

    ## Warning: Removed 1 rows containing missing values (geom_point).

    ## Warning: Removed 1 rows containing missing values (geom_text).

![](README_files/figure-markdown_github/Correlations%20among%20fitness%20proxies-2.png)

``` r
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
p11
```

    ## Warning: Removed 1 rows containing non-finite values (stat_smooth).

    ## Warning: Removed 1 rows containing missing values (geom_point).

    ## Warning: Removed 1 rows containing missing values (geom_text).

![](README_files/figure-markdown_github/Correlations%20among%20fitness%20proxies-3.png)

``` r
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
p12
```

![](README_files/figure-markdown_github/Correlations%20among%20fitness%20proxies-4.png)

``` r
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
p13
```

![](README_files/figure-markdown_github/Correlations%20among%20fitness%20proxies-5.png)

``` r
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
p14
```

![](README_files/figure-markdown_github/Correlations%20among%20fitness%20proxies-6.png)

``` r
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
p15
```

![](README_files/figure-markdown_github/Correlations%20among%20fitness%20proxies-7.png)

``` r
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

![](README_files/figure-markdown_github/Correlations%20among%20fitness%20proxies-8.png)
