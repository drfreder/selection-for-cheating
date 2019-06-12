Sampling nodules versus plants?
================
Megan Frederickson
2019-06-11

This code expands on my re-analyzis of the data in: Gano-Cohen KA, Wendlandt CE, Stokes PJ, Blanton MA, Quides KW, Zomorrodian A, Adinata ES, Sachs JL (2019) Interspecific conflict and the evolution of ineffective rhizobia. Ecology Letters. <https://doi.org/10.1111/ele.13247>

I downloaded their data from Dryad on April 16, 2019. The citation for the data package is: Gano-Cohen KA, Wendlandt CE, Stokes PJ, Blanton MA, Quides KW, Zomorrodian A, Adinata ES, Sachs JL (2019) Data from: Interspecific conflict and the evolution of ineffective rhizobia. Dryad Digital Repository. <https://doi.org/10.5061/dryad.cr65269>

What about sub-sampling nodules, instead of plants?
---------------------------------------------------

In my technical comment, I down-sampled the data by repeatedly choosing 2 or 4 plants at random from each population. What happens if, instead, I re-sample by the smallest number of nodules sampled in any population? I first re-sample to 39 nodules, as in the response to my technical comment:

``` r
min <- 39 #Minimum number of nods sampled 
pop <- c("UCR", "CLA", "GRI", "BMR", "YUC", "ANZ") #Create vector of population names
df.nods <- data.frame(Population=character(), haplotype = character(), locus=character(), n=double(), tot_n=double(), freq=double(), stringsAsFactors=FALSE) #Initialize empty frame to store all the sub-sampled data in

#Two for loops that loop through 10000 iterations for each population
for (i in 1:6) {
  tmp.pop <- pop[i]
for (j in 1:10000) {
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
    ## -0.14016 -0.07988 -0.04139  0.01049  0.41599 
    ## 
    ## Coefficients:
    ##                Estimate Std. Error t value Pr(>|t|)   
    ## (Intercept)      0.3354     0.1173   2.859  0.00865 **
    ## mean_log10_RGR  -0.2705     0.1305  -2.073  0.04909 * 
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 0.141 on 24 degrees of freedom
    ## Multiple R-squared:  0.1518, Adjusted R-squared:  0.1165 
    ## F-statistic: 4.297 on 1 and 24 DF,  p-value: 0.04909

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
    ## -0.14792 -0.10080 -0.03078  0.03620  0.28413 
    ## 
    ## Coefficients:
    ##                Estimate Std. Error t value Pr(>|t|)  
    ## (Intercept)      0.3387     0.1193   2.839   0.0101 *
    ## mean_log10_RGR  -0.1980     0.1305  -1.517   0.1449  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 0.1393 on 20 degrees of freedom
    ## Multiple R-squared:  0.1032, Adjusted R-squared:  0.05834 
    ## F-statistic: 2.301 on 1 and 20 DF,  p-value: 0.1449

Here is what I get if I choose 4 nodules at random 10000 times from each population:

``` r
min <- 4 #Minimum number of nods sampled 
df.nods <- data.frame(Population=character(), haplotype = character(), locus=character(), n=double(), tot_n=double(), freq=double(), stringsAsFactors=FALSE) #Initialize empty frame to store all the sub-sampled data in

#Two for loops that loop through 10000 iterations for each population
for (i in 1:6) {
  tmp.pop <- pop[i]
for (j in 1:10000) {
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
    ##       Min        1Q    Median        3Q       Max 
    ## -0.080383 -0.045908 -0.022372  0.007118  0.271173 
    ## 
    ## Coefficients:
    ##                Estimate Std. Error t value Pr(>|t|)    
    ## (Intercept)     0.43098    0.07173   6.009 3.34e-06 ***
    ## mean_log10_RGR -0.16043    0.07979  -2.011   0.0557 .  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 0.08623 on 24 degrees of freedom
    ## Multiple R-squared:  0.1442, Adjusted R-squared:  0.1085 
    ## F-statistic: 4.043 on 1 and 24 DF,  p-value: 0.05572

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
    ## -0.26639 -0.14653 -0.07187  0.10033  0.50931 
    ## 
    ## Coefficients:
    ##                Estimate Std. Error t value Pr(>|t|)    
    ## (Intercept)      0.8652     0.2036   4.249 0.000393 ***
    ## mean_log10_RGR  -0.4298     0.2227  -1.930 0.067965 .  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 0.2377 on 20 degrees of freedom
    ## Multiple R-squared:  0.1569, Adjusted R-squared:  0.1148 
    ## F-statistic: 3.723 on 1 and 20 DF,  p-value: 0.06797
