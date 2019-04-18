Are rhizobia under selection to cheat?
================
Megan Frederickson
2019-04-18

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

#Make sure factors are factors and numbers are numbers

data$Strain <- as.factor(data$Strain)
data$`Host Line` <- as.factor(data$`Host Line`)
data$Population <- as.factor(data$Population)
data$`Total nodules` <- as.numeric(data$`Total nodules`)
data$`Mean individual  nodule biomass (mg)` <- as.numeric(data$`Mean individual  nodule biomass (mg)`)
data$`Shoots mass (g)` <- as.numeric(data$`Shoots mass (g)`)
```

    ## Warning: NAs introduced by coercion

``` r
data$`Roots mass (g)` <- as.numeric(data$`Roots mass (g)`)
```

    ## Warning: NAs introduced by coercion

``` r
data$`Relative Growth` <- as.numeric(data$`Relative Growth`)
```

    ## Warning: NAs introduced by coercion

``` r
data$Block <- as.factor(data$Block)
```

Calculate genotype means
========================

The dataset has two measures of plant fitness and four measures of rhizobium fitness. The two measures of plant fitness are total plant biomass (i.e., shoot plus root mass) and relative growth rate (i.e., inoculated plant biomass divided by control plant biomass, hereafter "RGR"). The four measures of rhizobium fitness are total number of nodules, mean nodule mass, and two genotypic frequencies, which the authors abbreviate CHR and SI. See original paper for details.

First, we need to calculate the appropriate strain means for each variable. I did this using the emmeans package, which calculates estimated marginal means from a linear mixed model. I have tried to stick as close as possible to the analysis presented in the paper, which states that "all effects were coded as fixed except block, which was treated as a random effect" and that variables were log-transformed as needed to improve normality. I also calculated strain means only on data from sympatric hosts, again following the approach in the paper.

``` r
#Clean up dataset a bit
data <- subset(data, Strain != "control") #Exclude control, uninoculated plants
data <- subset(data, `Total nodules` > 0) #Exclude inoculated plants that did not form nodules, as these are cases where the inoculation may have failed
data$log_nodules <- log(data$`Total nodules`) #Log transform as needed, as per paper
data$log_mean_nod_biomass <- log(data$`Mean individual  nodule biomass (mg)`) #Log transform as needed, as per paper
data$plant_biomass <- data$`Shoots mass (g)` + data$`Roots mass (g)` #Combine shoot and root mass into one measure of plant biomass

#Divide into sympatric and universal hosts, as per paper
data_sym <- subset(data, `Host Line` != "A. heermannii" & `Host Line` != "UnH: Cla12.04" & `Host Line` != "UnL: Anz13.04")
data_uni <- subset(data, `Host Line` == "A. heermannii" | `Host Line` == "UnH: Cla12.04" | `Host Line` == "UnL: Anz13.04")

#Fit models 
lmm1 <- lmer(log_nodules~Population + Strain + `Host Line` + (1|Block), data=data_sym) #Total nodule number model
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
df.means <- as.data.frame(emmeans(lmm1, "Strain")) #Extract line means and save them to new dataframe
colnames(df.means) <- c("Strain", "Population", "tot_nod_lsmean", "tot_nod_SE", "tot_nod_df", "tot_nod_lowCL", "tot_nod_highCL") #Fix column names

lmm2 <- lmer(log_mean_nod_biomass~Population + Strain + `Host Line` + (1|Block), data=data_sym) #Mean nodule mass model
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
tmp <- as.data.frame(emmeans(lmm2, "Strain")) #Extract line means and add them to the dataframe
colnames(tmp) <- c("Strain", "Population", "nod_mass_lsmean", "nod_mass_SE", "nod_mass_df", "nod_mass_lowCL", "nod_mass_highCL")
df.means <- cbind(df.means, tmp[,3:7])

lmm3 <- lmer(plant_biomass~Population + Strain  + `Host Line` + (1|Block), data=data_sym) #Plant biomass model
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
tmp <- as.data.frame(emmeans(lmm3, "Strain")) #Extract line means and add them to the dataframe
colnames(tmp) <- c("Strain", "Population", "plant_biomass_lsmean", "plant_biomass_SE", "plant_biomass_df", "plant_biomass_lowCL", "plant_biomass_highCL")
df.means <- cbind(df.means, tmp[,3:7])

lmm4 <- lmer(log(`Relative Growth`)~Population + Strain  + `Host Line` + (1|Block), data=data_sym) #Relative growth rate, log-transformed as per the paper, but note this generates NAs because there are negative values
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
tmp <- as.data.frame(emmeans(lmm4, "Strain")) #Extract line means and add them to the dataframe
colnames(tmp) <- c("Strain", "Population", "RGR_lsmean", "RGR_SE", "RGR_df", "RGR_lowCL", "RGR_highCL")
df.means <- cbind(df.means, tmp[,3:7])

#Add CHR and SI values for each isolate; they are identical for all replicates of each isolate, so no need to model them
tmp <- data_sym %>% group_by(Strain) %>% summarize(CHR = mean(`CHR local abundance`), SI = mean(`SI local abunance`))
df.means <- merge(df.means, tmp, by="Strain")

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
  pop_sd_RGR=sd(RGR_lsmean)
  )
df.means <- merge(df.means, tmp, by="Population") #Merge data frames

#Relativize fitness by dividing by the mean fitness
df.means$tot_nod_std <- df.means$tot_nod_lsmean/df.means$pop_mean_tot_nod
df.means$nod_mass_std <- df.means$nod_mass_lsmean/df.means$pop_mean_nod_mass
df.means$CHR_std <- df.means$CHR/df.means$pop_mean_CHR
df.means$SI_std <- df.means$SI/df.means$pop_mean_SI

#Standardize traits by subtracting the population mean and dividing by the population standard deviation
df.means$plant_biomass_std <- (df.means$plant_biomass_lsmean - df.means$pop_mean_plant_biomass)/df.means$pop_sd_plant_biomass
df.means$RGR_std <- (df.means$RGR_lsmean - df.means$pop_mean_RGR)/df.means$pop_sd_RGR
```

Calculate selection gradients
=============================

Next, I adopt standard genetic selection analyses to estimate selection gradients by regressing relativized fitness against standardized trait values.

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
      geom_smooth(method="lm", se=TRUE)+
      geom_text(aes(label=Strain),hjust=0, vjust=0, size=3, check_overlap=TRUE)+
      ylab("Nodule number")+
      xlab("Plant biomass")

ggsave("fitness_correlation.png", p1, dpi=600)

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

lm3 <- lm(CHR_std~plant_biomass_std, data=df.means) #Model
summary(lm3)
```

    ## 
    ## Call:
    ## lm(formula = CHR_std ~ plant_biomass_std, data = df.means)
    ## 
    ## Residuals:
    ##     Min      1Q  Median      3Q     Max 
    ## -1.2634 -0.6186 -0.4071  0.3447  2.1260 
    ## 
    ## Coefficients:
    ##                   Estimate Std. Error t value Pr(>|t|)    
    ## (Intercept)         1.0000     0.1991   5.021 3.94e-05 ***
    ## plant_biomass_std  -0.4478     0.2271  -1.972   0.0602 .  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 1.015 on 24 degrees of freedom
    ## Multiple R-squared:  0.1394, Adjusted R-squared:  0.1036 
    ## F-statistic: 3.889 on 1 and 24 DF,  p-value: 0.06023

``` r
S3 <- paste0("S = ", round(summary(lm3)$coefficients[2,1], 3))
pval3 <- round(summary(lm3)$coefficients[2,4],3)
pval3 <- ifelse(pval3 <= 0.05, paste0("p = ", pval3), "ns")
S3 <- paste0(S3, ", ", pval3)

p3 <- ggplot(data=df.means, aes(y=CHR_std, x=plant_biomass_std))+
      geom_point()+
      geom_smooth(method="lm", se=FALSE, linetype="dashed")+
      geom_text(aes(label=Strain),hjust=0, vjust=0, size=3)+
      ylab("CHR")+
      xlab("Plant biomass")+
      guides(color=FALSE)

#SI and plant biomass

lm4 <- lm(SI_std~plant_biomass_std, data=df.means) #Model
summary(lm4)
```

    ## 
    ## Call:
    ## lm(formula = SI_std ~ plant_biomass_std, data = df.means)
    ## 
    ## Residuals:
    ##      Min       1Q   Median       3Q      Max 
    ## -0.82811 -0.54796 -0.02637  0.37176  1.81230 
    ## 
    ## Coefficients:
    ##                   Estimate Std. Error t value Pr(>|t|)    
    ## (Intercept)        1.00000    0.14758   6.776  1.8e-06 ***
    ## plant_biomass_std  0.07864    0.16907   0.465    0.647    
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 0.6763 on 19 degrees of freedom
    ##   (5 observations deleted due to missingness)
    ## Multiple R-squared:  0.01126,    Adjusted R-squared:  -0.04078 
    ## F-statistic: 0.2163 on 1 and 19 DF,  p-value: 0.6471

``` r
S4 <- paste0("S = ", round(summary(lm4)$coefficients[2,1], 3))
pval4 <- round(summary(lm4)$coefficients[2,4],3)
pval4 <- ifelse(pval4 <= 0.05, paste0("p = ", pval4), "ns")
S4 <- paste0(S4, ", ", pval4)

p4 <- ggplot(data=df.means, aes(y=SI_std, x=plant_biomass_std))+
      geom_point()+
      geom_smooth(method="lm", se=FALSE, linetype="dashed")+
      geom_text(aes(label=Strain),hjust=0, vjust=0, size=3)+
      ylab("SI")+
      xlab("Plant biomass")

#Total nodules and RGR

lm5 <- lm(tot_nod_std~RGR_std, data=df.means) #Model
summary(lm5)
```

    ## 
    ## Call:
    ## lm(formula = tot_nod_std ~ RGR_std, data = df.means)
    ## 
    ## Residuals:
    ##      Min       1Q   Median       3Q      Max 
    ## -0.48398 -0.03390  0.00029  0.04925  0.20545 
    ## 
    ## Coefficients:
    ##               Estimate Std. Error t value Pr(>|t|)    
    ## (Intercept)  1.000e+00  2.454e-02  40.752   <2e-16 ***
    ## RGR_std     -5.287e-05  2.798e-02  -0.002    0.999    
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 0.1251 on 24 degrees of freedom
    ## Multiple R-squared:  1.488e-07,  Adjusted R-squared:  -0.04167 
    ## F-statistic: 3.571e-06 on 1 and 24 DF,  p-value: 0.9985

``` r
S5 <- paste0("S = ", round(summary(lm5)$coefficients[2,1], 3))
pval5 <- round(summary(lm5)$coefficients[2,4],3)
pval5 <- ifelse(pval5 <= 0.05, paste0("p = ", pval5), "ns")
S5 <- paste0(S5, ", ", pval5)

p5 <- ggplot(data=df.means, aes(y=tot_nod_std, x=RGR_std))+
      geom_point()+
      geom_smooth(method="lm", se=FALSE, linetype="dashed")+
      geom_text(aes(label=Strain),hjust=0, vjust=0, size=3)+
      ylab("Nodule number")+
      xlab("RGR")

#Nodule mass and RGR

lm6 <- lm(nod_mass_std~RGR_std, data=df.means) #Model
summary(lm6)
```

    ## 
    ## Call:
    ## lm(formula = nod_mass_std ~ RGR_std, data = df.means)
    ## 
    ## Residuals:
    ##      Min       1Q   Median       3Q      Max 
    ## -0.31423 -0.06893 -0.00134  0.09635  0.33829 
    ## 
    ## Coefficients:
    ##             Estimate Std. Error t value Pr(>|t|)    
    ## (Intercept)  1.00000    0.03324  30.081   <2e-16 ***
    ## RGR_std     -0.01171    0.03790  -0.309     0.76    
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 0.1695 on 24 degrees of freedom
    ## Multiple R-squared:  0.003959,   Adjusted R-squared:  -0.03754 
    ## F-statistic: 0.09538 on 1 and 24 DF,  p-value: 0.7601

``` r
S6 <- paste0("S = ", round(summary(lm6)$coefficients[2,1], 3))
pval6 <- round(summary(lm6)$coefficients[2,4],3)
pval6 <- ifelse(pval6 <= 0.05, paste0("p = ", pval6), "ns")
S6 <- paste0(S6, ", ", pval6)

p6 <- ggplot(data=df.means, aes(y=nod_mass_std, x=RGR_std))+
      geom_point()+
      geom_smooth(method="lm", se=FALSE, linetype="dashed")+
      geom_text(aes(label=Strain),hjust=0, vjust=0, size=3)+
      ylab("Nodule mass")+
      xlab("RGR")

#CHR and RGR

lm7 <- lm(CHR_std~RGR_std, data=df.means) #Model
summary(lm7)
```

    ## 
    ## Call:
    ## lm(formula = CHR_std ~ RGR_std, data = df.means)
    ## 
    ## Residuals:
    ##     Min      1Q  Median      3Q     Max 
    ## -1.2209 -0.6395 -0.3932  0.3639  2.3859 
    ## 
    ## Coefficients:
    ##             Estimate Std. Error t value Pr(>|t|)    
    ## (Intercept)   1.0000     0.2034   4.917 5.12e-05 ***
    ## RGR_std      -0.3843     0.2319  -1.657     0.11    
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 1.037 on 24 degrees of freedom
    ## Multiple R-squared:  0.1027, Adjusted R-squared:  0.06531 
    ## F-statistic: 2.747 on 1 and 24 DF,  p-value: 0.1105

``` r
S7 <- paste0("S = ", round(summary(lm7)$coefficients[2,1], 3))
pval7 <- round(summary(lm7)$coefficients[2,4],3)
pval7 <- ifelse(pval7 <= 0.05, paste0("p = ", pval7), "ns")
S7 <- paste0(S7, ", ", pval7)

p7 <- ggplot(data=df.means, aes(y=CHR_std, x=RGR_std))+
      geom_point()+
      geom_smooth(method="lm", se=FALSE, linetype="dashed")+
      geom_text(aes(label=Strain),hjust=0, vjust=0, size=3)+
      ylab("CHR")+
      xlab("RGR")

#SI and RGR

lm8 <- lm(SI_std~RGR_std, data=df.means) #Model
summary(lm8)
```

    ## 
    ## Call:
    ## lm(formula = SI_std ~ RGR_std, data = df.means)
    ## 
    ## Residuals:
    ##      Min       1Q   Median       3Q      Max 
    ## -0.79047 -0.57746 -0.01811  0.34575  1.81801 
    ## 
    ## Coefficients:
    ##             Estimate Std. Error t value Pr(>|t|)    
    ## (Intercept)  1.00000    0.14815   6.750 1.89e-06 ***
    ## RGR_std      0.04443    0.16972   0.262    0.796    
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 0.6789 on 19 degrees of freedom
    ##   (5 observations deleted due to missingness)
    ## Multiple R-squared:  0.003594,   Adjusted R-squared:  -0.04885 
    ## F-statistic: 0.06853 on 1 and 19 DF,  p-value: 0.7963

``` r
S8 <- paste0("S = ", round(summary(lm8)$coefficients[2,1], 3))
pval8 <- round(summary(lm8)$coefficients[2,4],3)
pval8 <- ifelse(pval8 <= 0.05, paste0("p = ", pval8), "ns")
S8 <- paste0(S8, ", ", pval8)

p8 <- ggplot(data=df.means, aes(y=SI_std, x=RGR_std))+
      geom_point()+
      geom_smooth(method="lm", se=FALSE, linetype="dashed")+
      geom_text(aes(label=Strain),hjust=0, vjust=0, size=3)+
      ylab("SI")+
      xlab("RGR")

plot_grid(p1,p2,p3,p4, labels=c(S1, S2, S3, S4), hjust = c(-0.8, -1, -1, -1), scale=0.9)
```

    ## Warning: Removed 5 rows containing non-finite values (stat_smooth).

    ## Warning: Removed 5 rows containing missing values (geom_point).

    ## Warning: Removed 5 rows containing missing values (geom_text).

![](README_files/figure-markdown_github/Selection%20gradients-1.png)

``` r
plot_grid(p5,p6,p7,p8, labels=c(S5, S6, S7, S8), hjust = c(-1.7, -1, -1, -1), scale = 0.9)
```

    ## Warning: Removed 5 rows containing non-finite values (stat_smooth).

    ## Warning: Removed 5 rows containing missing values (geom_point).

    ## Warning: Removed 5 rows containing missing values (geom_text).

![](README_files/figure-markdown_github/Selection%20gradients-2.png)

Why do different fitness measures give different answers?
=========================================================

``` r
#Correlations among rhizobium fitness proxies
lm9 <- lm(SI_std~CHR_std, data=df.means) #Model
summary(lm9)
```

    ## 
    ## Call:
    ## lm(formula = SI_std ~ CHR_std, data = df.means)
    ## 
    ## Residuals:
    ##      Min       1Q   Median       3Q      Max 
    ## -0.85915 -0.56948  0.01595  0.27896  1.85360 
    ## 
    ## Coefficients:
    ##             Estimate Std. Error t value Pr(>|t|)    
    ## (Intercept)  0.96279    0.20532   4.689  0.00016 ***
    ## CHR_std      0.03721    0.14216   0.262  0.79635    
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 0.6789 on 19 degrees of freedom
    ##   (5 observations deleted due to missingness)
    ## Multiple R-squared:  0.003592,   Adjusted R-squared:  -0.04885 
    ## F-statistic: 0.0685 on 1 and 19 DF,  p-value: 0.7963

``` r
p9 <- ggplot(data=df.means, aes(y=SI_std, x=CHR_std))+
      geom_point()+
      geom_smooth(method="lm", se=FALSE)+
      xlab("CHR")+
      ylab("SI")+
      geom_text(aes(label=Strain),hjust=0, vjust=0, size=3)
p9
```

    ## Warning: Removed 5 rows containing non-finite values (stat_smooth).

    ## Warning: Removed 5 rows containing missing values (geom_point).

    ## Warning: Removed 5 rows containing missing values (geom_text).

![](README_files/figure-markdown_github/Correlations%20among%20fitness%20proxies,%20-1.png)

``` r
lm10 <- lm(SI_std~nod_mass_std, data=df.means) #Model
summary(lm10)
```

    ## 
    ## Call:
    ## lm(formula = SI_std ~ nod_mass_std, data = df.means)
    ## 
    ## Residuals:
    ##      Min       1Q   Median       3Q      Max 
    ## -0.97665 -0.44368 -0.00975  0.19633  1.56611 
    ## 
    ## Coefficients:
    ##              Estimate Std. Error t value Pr(>|t|)   
    ## (Intercept)    2.5878     0.8256   3.134  0.00546 **
    ## nod_mass_std  -1.5878     0.8144  -1.950  0.06614 . 
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 0.6208 on 19 degrees of freedom
    ##   (5 observations deleted due to missingness)
    ## Multiple R-squared:  0.1667, Adjusted R-squared:  0.1228 
    ## F-statistic: 3.801 on 1 and 19 DF,  p-value: 0.06614

``` r
p10 <- ggplot(data=df.means, aes(y=SI_std, x=nod_mass_std))+
      geom_point()+
      geom_smooth(method="lm", se=FALSE)+
      xlab("Nodule mass")+
      ylab("SI")+
      geom_text(aes(label=Strain),hjust=0, vjust=0, size=3)
p10
```

    ## Warning: Removed 5 rows containing non-finite values (stat_smooth).

    ## Warning: Removed 5 rows containing missing values (geom_point).

    ## Warning: Removed 5 rows containing missing values (geom_text).

![](README_files/figure-markdown_github/Correlations%20among%20fitness%20proxies,%20-2.png)

``` r
lm11 <- lm(SI_std~tot_nod_std, data=df.means) #Model
summary(lm11)
```

    ## 
    ## Call:
    ## lm(formula = SI_std ~ tot_nod_std, data = df.means)
    ## 
    ## Residuals:
    ##      Min       1Q   Median       3Q      Max 
    ## -0.85654 -0.56652 -0.03333  0.31532  1.79756 
    ## 
    ## Coefficients:
    ##             Estimate Std. Error t value Pr(>|t|)
    ## (Intercept)    2.188      2.561   0.854    0.403
    ## tot_nod_std   -1.188      2.557  -0.465    0.647
    ## 
    ## Residual standard error: 0.6763 on 19 degrees of freedom
    ##   (5 observations deleted due to missingness)
    ## Multiple R-squared:  0.01124,    Adjusted R-squared:  -0.0408 
    ## F-statistic: 0.216 on 1 and 19 DF,  p-value: 0.6474

``` r
p11 <- ggplot(data=df.means, aes(y=SI_std, x=tot_nod_std))+
      geom_point()+
      geom_smooth(method="lm", se=FALSE)+
      xlab("Nodule number")+
      ylab("SI")+
      geom_text(aes(label=Strain),hjust=0, vjust=0, size=3)
p11
```

    ## Warning: Removed 5 rows containing non-finite values (stat_smooth).

    ## Warning: Removed 5 rows containing missing values (geom_point).

    ## Warning: Removed 5 rows containing missing values (geom_text).

![](README_files/figure-markdown_github/Correlations%20among%20fitness%20proxies,%20-3.png)

``` r
lm12 <- lm(CHR_std~nod_mass_std, data=df.means) #Model
summary(lm12)
```

    ## 
    ## Call:
    ## lm(formula = CHR_std ~ nod_mass_std, data = df.means)
    ## 
    ## Residuals:
    ##     Min      1Q  Median      3Q     Max 
    ## -0.9856 -0.7460 -0.4616  0.5164  2.5514 
    ## 
    ## Coefficients:
    ##              Estimate Std. Error t value Pr(>|t|)
    ## (Intercept)     2.017      1.316   1.532    0.139
    ## nod_mass_std   -1.017      1.299  -0.783    0.441
    ## 
    ## Residual standard error: 1.081 on 24 degrees of freedom
    ## Multiple R-squared:  0.0249, Adjusted R-squared:  -0.01573 
    ## F-statistic: 0.6128 on 1 and 24 DF,  p-value: 0.4414

``` r
p12 <- ggplot(data=df.means, aes(y=CHR_std, x=nod_mass_std))+
      geom_point()+
      geom_smooth(method="lm", se=FALSE)+
      xlab("Nodule mass")+
      ylab("CHR")+
      geom_text(aes(label=Strain),hjust=0, vjust=0, size=3)
p12
```

![](README_files/figure-markdown_github/Correlations%20among%20fitness%20proxies,%20-4.png)

``` r
lm13 <- lm(CHR_std~tot_nod_std, data=df.means) #Model
summary(lm13)
```

    ## 
    ## Call:
    ## lm(formula = CHR_std ~ tot_nod_std, data = df.means)
    ## 
    ## Residuals:
    ##     Min      1Q  Median      3Q     Max 
    ## -0.9292 -0.7286 -0.4693  0.5858  2.5298 
    ## 
    ## Coefficients:
    ##             Estimate Std. Error t value Pr(>|t|)
    ## (Intercept)   0.7937     1.7981   0.441    0.663
    ## tot_nod_std   0.2063     1.7853   0.116    0.909
    ## 
    ## Residual standard error: 1.094 on 24 degrees of freedom
    ## Multiple R-squared:  0.0005562,  Adjusted R-squared:  -0.04109 
    ## F-statistic: 0.01336 on 1 and 24 DF,  p-value: 0.909

``` r
p13 <- ggplot(data=df.means, aes(y=CHR_std, x=tot_nod_std))+
      geom_point()+
      geom_smooth(method="lm", se=FALSE)+
      xlab("Nodule number")+
      ylab("CHR")+
      geom_text(aes(label=Strain),hjust=0, vjust=0, size=3)
p13
```

![](README_files/figure-markdown_github/Correlations%20among%20fitness%20proxies,%20-5.png)

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


plot_grid(p9, p10, p11, p12, p13, p14)
```

    ## Warning: Removed 5 rows containing non-finite values (stat_smooth).

    ## Warning: Removed 5 rows containing missing values (geom_point).

    ## Warning: Removed 5 rows containing missing values (geom_text).

    ## Warning: Removed 5 rows containing non-finite values (stat_smooth).

    ## Warning: Removed 5 rows containing missing values (geom_point).

    ## Warning: Removed 5 rows containing missing values (geom_text).

    ## Warning: Removed 5 rows containing non-finite values (stat_smooth).

    ## Warning: Removed 5 rows containing missing values (geom_point).

    ## Warning: Removed 5 rows containing missing values (geom_text).

![](README_files/figure-markdown_github/Correlations%20among%20fitness%20proxies,%20-6.png)

``` r
#Corrlation between plant fitness proxies

lm15 <- lm(RGR_std~plant_biomass_std, data=df.means)
summary(lm15)
```

    ## 
    ## Call:
    ## lm(formula = RGR_std ~ plant_biomass_std, data = df.means)
    ## 
    ## Residuals:
    ##      Min       1Q   Median       3Q      Max 
    ## -1.82251 -0.27160  0.05819  0.22504  1.63749 
    ## 
    ## Coefficients:
    ##                    Estimate Std. Error t value Pr(>|t|)    
    ## (Intercept)       2.057e-16  1.124e-01   0.000        1    
    ## plant_biomass_std 7.783e-01  1.282e-01   6.073 2.85e-06 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 0.5732 on 24 degrees of freedom
    ## Multiple R-squared:  0.6058, Adjusted R-squared:  0.5894 
    ## F-statistic: 36.88 on 1 and 24 DF,  p-value: 2.849e-06

``` r
p15 <- ggplot(data=df.means, aes(y=RGR_std, x=plant_biomass_std))+
      geom_point()+
      geom_smooth(method="lm", se=FALSE)+
      geom_text(aes(label=Strain),hjust=0, vjust=0, size=3)+
      xlab("Plant biomass")+
      ylab("RGR")
p15
```

![](README_files/figure-markdown_github/Correlations%20among%20fitness%20proxies,%20-7.png)
