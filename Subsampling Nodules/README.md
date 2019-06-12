Sampling nodules versus plants?
================
Megan Frederickson
2019-06-11

This code expands on my re-analyzis of the data in: Gano-Cohen KA, Wendlandt CE, Stokes PJ, Blanton MA, Quides KW, Zomorrodian A, Adinata ES, Sachs JL (2019) Interspecific conflict and the evolution of ineffective rhizobia. Ecology Letters. <https://doi.org/10.1111/ele.13247>

I downloaded their data from Dryad on April 16, 2019. The citation for the data package is: Gano-Cohen KA, Wendlandt CE, Stokes PJ, Blanton MA, Quides KW, Zomorrodian A, Adinata ES, Sachs JL (2019) Data from: Interspecific conflict and the evolution of ineffective rhizobia. Dryad Digital Repository. <https://doi.org/10.5061/dryad.cr65269>

What about sub-sampling nodules, instead of plants?
---------------------------------------------------

In my technical comment, I down-sampled the data by choosing 2 or 4 plants at random from each population, iterating 100 times per population, and calculating mean genotype frequencies. What happens if, instead, I down-sampled by the smallest number of nodules sampled in any given population (39 for CHR frequencies, or 4 for SI frequencies)?

Here is what I get if I choose 39 nodules at random 1000 times from each population:

``` r
min <- 39 #Minimum number of nods sampled 
pop <- c("UCR", "CLA", "GRI", "BMR", "YUC", "ANZ") #Create vector of population names
df.nods <- data.frame(Population=character(), haplotype = character(), locus=character(), n=double(), tot_n=double(), freq=double(), stringsAsFactors=FALSE) #Initialize empty frame to store all the sub-sampled data in

#Two for loops that loop through 100 iterations for each population
for (i in 1:6) {
  tmp.pop <- pop[i]
for (j in 1:100) {
  tmp <- subset(table_S1, Population == tmp.pop)
  tmp.nods <- sample(unique(tmp$Full_Strain_Name), min, replace = FALSE)
  tmp.data <- subset(tmp, tmp$Full_Strain_Name %in% tmp.nods)
  tmp.data.long <- gather(tmp.data, locus, haplotype, glnII_Haplotype:CHR_haplotype, factor_key=TRUE)
  tmp.data.long <- subset(tmp.data.long, locus == "CHR_haplotype" | locus == "SI_haplotype")
  tmp.data.long <- subset(subset(tmp.data.long, haplotype != "n/an/a"), haplotype != "n/a_n/a")
  tmp2 <- tmp.data.long %>% group_by(Population, haplotype, locus) %>% summarize(n=n())
  n_SI <- length(subset(tmp.data.long, locus == "SI_haplotype")[[1]])
  n_CHR <- length(subset(tmp.data.long, locus == "CHR_haplotype")[[1]])
  tmp2$tot_n <- ifelse(tmp2$locus == "CHR_haplotype", n_CHR, n_SI)
  tmp2$freq <- tmp2$n/tmp2$tot_n
  df.nods<-rbind(df.nods, as.data.frame(tmp2))
  }
}

#Merge re-sampled and original data by haplotype and population
df$popCHR <- paste0(df$Population, df$CHR_Haplotype)
df$popSI <- paste0(df$Population, df$SI_haplotye)
df.nods$pophaplotype <- paste0(df.nods$Population, df.nods$haplotype)
df.nods.CHR <- merge(df.nods, df, by.x="pophaplotype", by.y="popCHR")
df.nods.SI <- merge(df.nods, df, by.x="pophaplotype", by.y="popSI")

#Calculate mean genotype frequencies per strain from 100 iterations
CHR.sum <- df.nods.CHR %>% group_by(Population.x, haplotype) %>% summarize(n=n(), mean_CHR_freq=mean(freq, na.rm=TRUE), sd=sd(freq, na.rm=TRUE), se=sd/sqrt(n))
CHR.sum$popCHR <- paste0(CHR.sum$Population.x, CHR.sum$haplotype)
CHR.sum <- merge(CHR.sum, df, by="popCHR")

SI.sum <- df.nods.SI %>% group_by(Population.x, haplotype) %>% summarize(n=n(), mean_SI_freq=mean(freq, na.rm=TRUE), sd=sd(freq, na.rm=TRUE), se=sd/sqrt(n))
SI.sum$popSI <- paste0(SI.sum$Population.x, SI.sum$haplotype)
SI.sum <- merge(SI.sum, df, by="popSI")

#Fit model for mean CHR genotype frequency from re-sampling and RGR
model1 <- lm(mean_CHR_freq~mean_log10_RGR, data=CHR.sum)
summary(model1) 
```

    ## 
    ## Call:
    ## lm(formula = mean_CHR_freq ~ mean_log10_RGR, data = CHR.sum)
    ## 
    ## Residuals:
    ##      Min       1Q   Median       3Q      Max 
    ## -0.14036 -0.07995 -0.04212  0.01062  0.41087 
    ## 
    ## Coefficients:
    ##                Estimate Std. Error t value Pr(>|t|)   
    ## (Intercept)      0.3360     0.1167   2.880  0.00824 **
    ## mean_log10_RGR  -0.2711     0.1298  -2.089  0.04748 * 
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 0.1402 on 24 degrees of freedom
    ## Multiple R-squared:  0.1539, Adjusted R-squared:  0.1186 
    ## F-statistic: 4.364 on 1 and 24 DF,  p-value: 0.04748

``` r
#Fit model for mean SI genotype frequency from re-sampling and RGR
model2 <- lm(mean_SI_freq~mean_log10_RGR, data=SI.sum)
summary(model2) 
```

    ## 
    ## Call:
    ## lm(formula = mean_SI_freq ~ mean_log10_RGR, data = SI.sum)
    ## 
    ## Residuals:
    ##      Min       1Q   Median       3Q      Max 
    ## -0.15050 -0.10551 -0.03290  0.03176  0.31547 
    ## 
    ## Coefficients:
    ##                Estimate Std. Error t value Pr(>|t|)  
    ## (Intercept)      0.3392     0.1246   2.723   0.0131 *
    ## mean_log10_RGR  -0.1946     0.1362  -1.429   0.1685  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 0.1454 on 20 degrees of freedom
    ## Multiple R-squared:  0.09261,    Adjusted R-squared:  0.04724 
    ## F-statistic: 2.041 on 1 and 20 DF,  p-value: 0.1685

Here is what I get if I choose 4 nodules at random 100 times from each population:

``` r
min <- 4 #Minimum number of nods sampled 
df.nods <- data.frame(Population=character(), haplotype = character(), locus=character(), n=double(), tot_n=double(), freq=double(), stringsAsFactors=FALSE) #Initialize empty frame to store all the sub-sampled data in

#Two for loops that loop through 100 iterations for each population
for (i in 1:6) {
  tmp.pop <- pop[i]
for (j in 1:100) {
  tmp <- subset(table_S1, Population == tmp.pop)
  tmp.nods <- sample(unique(tmp$Full_Strain_Name), min, replace = FALSE)
  tmp.data <- subset(tmp, tmp$Full_Strain_Name %in% tmp.nods)
  tmp.data.long <- gather(tmp.data, locus, haplotype, glnII_Haplotype:CHR_haplotype, factor_key=TRUE)
  tmp.data.long <- subset(tmp.data.long, locus == "CHR_haplotype" | locus == "SI_haplotype")
  tmp.data.long <- subset(subset(tmp.data.long, haplotype != "n/an/a"), haplotype != "n/a_n/a")
  tmp2 <- tmp.data.long %>% group_by(Population, haplotype, locus) %>% summarize(n=n())
  n_SI <- length(subset(tmp.data.long, locus == "SI_haplotype")[[1]])
  n_CHR <- length(subset(tmp.data.long, locus == "CHR_haplotype")[[1]])
  tmp2$tot_n <- ifelse(tmp2$locus == "CHR_haplotype", n_CHR, n_SI)
  tmp2$freq <- tmp2$n/tmp2$tot_n
  df.nods<-rbind(df.nods, as.data.frame(tmp2))
  }
}

#Merge re-sampled and original data by haplotype and population
df.nods$pophaplotype <- paste0(df.nods$Population, df.nods$haplotype)
df.nods.CHR <- merge(df.nods, df, by.x="pophaplotype", by.y="popCHR")
df.nods.SI <- merge(df.nods, df, by.x="pophaplotype", by.y="popSI")

#Calculate mean genotype frequencies per strain from 100 iterations
CHR.sum <- df.nods.CHR %>% group_by(Population.x, haplotype) %>% summarize(n=n(), mean_CHR_freq=mean(freq, na.rm=TRUE), sd=sd(freq, na.rm=TRUE), se=sd/sqrt(n))
CHR.sum$popCHR <- paste0(CHR.sum$Population.x, CHR.sum$haplotype)
CHR.sum <- merge(CHR.sum, df, by="popCHR")

SI.sum <- df.nods.SI %>% group_by(Population.x, haplotype) %>% summarize(n=n(), mean_SI_freq=mean(freq, na.rm=TRUE), sd=sd(freq, na.rm=TRUE), se=sd/sqrt(n))
SI.sum$popSI <- paste0(SI.sum$Population.x, SI.sum$haplotype)
SI.sum <- merge(SI.sum, df, by="popSI")

#Fit model for mean CHR genotype frequency from re-sampling and RGR
model3 <- lm(mean_CHR_freq~mean_log10_RGR, data=CHR.sum)
summary(model3) 
```

    ## 
    ## Call:
    ## lm(formula = mean_CHR_freq ~ mean_log10_RGR, data = CHR.sum)
    ## 
    ## Residuals:
    ##      Min       1Q   Median       3Q      Max 
    ## -0.08521 -0.04653 -0.01962  0.01696  0.27275 
    ## 
    ## Coefficients:
    ##                Estimate Std. Error t value Pr(>|t|)    
    ## (Intercept)     0.45023    0.07039   6.396  1.3e-06 ***
    ## mean_log10_RGR -0.18344    0.07830  -2.343   0.0278 *  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 0.08463 on 24 degrees of freedom
    ## Multiple R-squared:  0.1861, Adjusted R-squared:  0.1522 
    ## F-statistic: 5.488 on 1 and 24 DF,  p-value: 0.02777

``` r
#Fit model for mean SI genotype frequency from re-sampling and RGR
model4 <- lm(mean_SI_freq~mean_log10_RGR, data=SI.sum)
summary(model4) 
```

    ## 
    ## Call:
    ## lm(formula = mean_SI_freq ~ mean_log10_RGR, data = SI.sum)
    ## 
    ## Residuals:
    ##      Min       1Q   Median       3Q      Max 
    ## -0.26128 -0.14504 -0.08692  0.09016  0.58375 
    ## 
    ## Coefficients:
    ##                Estimate Std. Error t value Pr(>|t|)    
    ## (Intercept)      0.9036     0.2158   4.187 0.000454 ***
    ## mean_log10_RGR  -0.4781     0.2360  -2.025 0.056369 .  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 0.2519 on 20 degrees of freedom
    ## Multiple R-squared:  0.1702, Adjusted R-squared:  0.1287 
    ## F-statistic: 4.103 on 1 and 20 DF,  p-value: 0.05637
