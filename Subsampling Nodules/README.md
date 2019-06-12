Sampling nodules versus plants?
================
Megan Frederickson
2019-06-12

This code expands on my re-analyzis of the data in: Gano-Cohen KA, Wendlandt CE, Stokes PJ, Blanton MA, Quides KW, Zomorrodian A, Adinata ES, Sachs JL (2019) Interspecific conflict and the evolution of ineffective rhizobia. Ecology Letters. <https://doi.org/10.1111/ele.13247>

I downloaded their data from Dryad on April 16, 2019. The citation for the data package is: Gano-Cohen KA, Wendlandt CE, Stokes PJ, Blanton MA, Quides KW, Zomorrodian A, Adinata ES, Sachs JL (2019) Data from: Interspecific conflict and the evolution of ineffective rhizobia. Dryad Digital Repository. <https://doi.org/10.5061/dryad.cr65269>

In my technical comment, I down-sampled the data by repeatedly choosing 2 or 4 plants at random from each population. What happens if, instead, I re-sample by the smallest number of nodules sampled in any population? I first re-sample to 39 nodules, as in the response to my technical comment:

``` r
min <- 39 #Minimum number of nods sampled 
pop <- c("UCR", "CLA", "GRI", "BMR", "YUC", "ANZ") #Create vector of population names
df.nods <- data.frame(Population=character(), haplotype = character(), locus=character(), n=double(), tot_n=double(), freq=double(), stringsAsFactors=FALSE) #Initialize empty frame to store all the sub-sampled data in

#Two for loops that loop through 1000 iterations for each population
for (i in 1:6) { #Loop through each of 6 populations
  tmp.pop <- pop[i]
for (j in 1:1000) { #For each population, sample nodules at random 10000 times
  tmp <- subset(table_S1, Population == tmp.pop) #Subset haplotype data to just the current population
  tmp.nods <- sample(unique(tmp$Full_Strain_Name), min, replace = FALSE) #Sample nodules at random
  tmp.data <- subset(tmp, tmp$Full_Strain_Name %in% tmp.nods) #Subset to just sampled nodules
  tmp.data.long <- gather(tmp.data, locus, haplotype, glnII_Haplotype:CHR_haplotype, factor_key=TRUE) #Make wide data long
  tmp.data.long <- subset(tmp.data.long, locus == "CHR_haplotype" | locus == "SI_haplotype") #Use on concatenated haplotypes
  tmp.data.long <- subset(subset(tmp.data.long, haplotype != "n/an/a"), haplotype != "n/a_n/a") #Remove missing data (NAs)
  tmp2 <- tmp.data.long %>% group_by(Population, haplotype, locus) %>% summarize(n=n()) #Count number of times each haplotype was chosen
  n_SI <- length(subset(tmp.data.long, locus == "SI_haplotype")[[1]]) #Count total number of samples (should always be 39 for CHR)
  n_CHR <- length(subset(tmp.data.long, locus == "CHR_haplotype")[[1]]) #Count total number of samples (often less than 39 for SI because of missing data)
  tmp2$tot_n <- ifelse(tmp2$locus == "CHR_haplotype", n_CHR, n_SI) #Fill in appropriate denominator (i.e., CHR or SI)
  tmp2$freq <- tmp2$n/tmp2$tot_n #Divide number of times each haplotype was sampled by total number of samples to obtain frequency
  df.nods<-rbind(df.nods, as.data.frame(tmp2)) #Add to dataframe
  }
}

#Merge re-sampled and original data by haplotype and population
df$popCHR <- paste0(df$Population, df$CHR_Haplotype)
df$popSI <- paste0(df$Population, df$SI_haplotye)
df.nods$pophaplotype <- paste0(df.nods$Population, df.nods$haplotype)
df.nods.CHR <- merge(df.nods, df, by.x="pophaplotype", by.y="popCHR")
df.nods.SI <- merge(df.nods, df, by.x="pophaplotype", by.y="popSI")

#Calculate mean genotype frequencies per strain from 1000 iterations
CHR.sum <- df.nods.CHR %>% group_by(Population.x, haplotype) %>% summarize(n=n(), mean_CHR_freq=mean(freq, na.rm=TRUE), sd=sd(freq, na.rm=TRUE), se=sd/sqrt(n))
CHR.sum$popCHR <- paste0(CHR.sum$Population.x, CHR.sum$haplotype)
CHR.sum <- merge(CHR.sum, df, by="popCHR")
SI.sum <- df.nods.SI %>% group_by(Population.x, haplotype) %>% summarize(n=n(), mean_SI_freq=mean(freq, na.rm=TRUE), sd=sd(freq, na.rm=TRUE), se=sd/sqrt(n))
SI.sum$popSI <- paste0(SI.sum$Population.x, SI.sum$haplotype)
SI.sum <- merge(SI.sum, df, by="popSI")
```

OK, now I model the relationship between the re-sampled rhizobium fitness data and the original symbiotic effectiveness values, for both the chromosome (CHR) and symbiosis island (SI) genotype frequencies. The output shows a significant relationship between re-sampled CHR genotype frequency and symbiotic effectiveness, but a non-significant relationship between re-sampled SI genotype frequency and symbiotic effectiveness.

``` r
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
    ## -0.14028 -0.08000 -0.04167  0.01038  0.41471 
    ## 
    ## Coefficients:
    ##                Estimate Std. Error t value Pr(>|t|)   
    ## (Intercept)      0.3355     0.1172   2.863  0.00857 **
    ## mean_log10_RGR  -0.2705     0.1304  -2.075  0.04887 * 
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 0.1409 on 24 degrees of freedom
    ## Multiple R-squared:  0.1521, Adjusted R-squared:  0.1168 
    ## F-statistic: 4.306 on 1 and 24 DF,  p-value: 0.04887

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
    ## -0.14846 -0.10107 -0.03100  0.03552  0.29173 
    ## 
    ## Coefficients:
    ##                Estimate Std. Error t value Pr(>|t|)  
    ## (Intercept)      0.3392     0.1203   2.820   0.0106 *
    ## mean_log10_RGR  -0.1978     0.1316  -1.504   0.1483  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 0.1404 on 20 degrees of freedom
    ## Multiple R-squared:  0.1016, Adjusted R-squared:  0.05666 
    ## F-statistic: 2.261 on 1 and 20 DF,  p-value: 0.1483

Here is what I get if instead I choose 4 nodules at random 5000 times from each population.

``` r
min <- 4 #Minimum number of nods sampled 
df.nods <- data.frame(Population=character(), haplotype = character(), locus=character(), n=double(), tot_n=double(), freq=double(), stringsAsFactors=FALSE) #Initialize empty frame to store all the sub-sampled data in

#Two for loops that loop through 1000 iterations for each population
for (i in 1:6) {
  tmp.pop <- pop[i]
for (j in 1:1000) {
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

#Calculate mean genotype frequencies per strain from 1000 iterations
CHR.sum <- df.nods.CHR %>% group_by(Population.x, haplotype) %>% summarize(n=n(), mean_CHR_freq=mean(freq, na.rm=TRUE), sd=sd(freq, na.rm=TRUE), se=sd/sqrt(n))
CHR.sum$popCHR <- paste0(CHR.sum$Population.x, CHR.sum$haplotype)
CHR.sum <- merge(CHR.sum, df, by="popCHR")
SI.sum <- df.nods.SI %>% group_by(Population.x, haplotype) %>% summarize(n=n(), mean_SI_freq=mean(freq, na.rm=TRUE), sd=sd(freq, na.rm=TRUE), se=sd/sqrt(n))
SI.sum$popSI <- paste0(SI.sum$Population.x, SI.sum$haplotype)
SI.sum <- merge(SI.sum, df, by="popSI")
```

The output shows non-significant relationships between rhizobium fitness and symbiotic effectiveness regardless of locus.

``` r
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
    ## -0.081165 -0.046378 -0.022959  0.006756  0.284777 
    ## 
    ## Coefficients:
    ##                Estimate Std. Error t value Pr(>|t|)    
    ## (Intercept)      0.4327     0.0739   5.855 4.87e-06 ***
    ## mean_log10_RGR  -0.1619     0.0822  -1.969   0.0606 .  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 0.08884 on 24 degrees of freedom
    ## Multiple R-squared:  0.1391, Adjusted R-squared:  0.1033 
    ## F-statistic: 3.878 on 1 and 24 DF,  p-value: 0.06056

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
    ## -0.27441 -0.14310 -0.06634  0.10033  0.47543 
    ## 
    ## Coefficients:
    ##                Estimate Std. Error t value Pr(>|t|)    
    ## (Intercept)      0.8805     0.2010   4.381 0.000289 ***
    ## mean_log10_RGR  -0.4457     0.2198  -2.027 0.056155 .  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 0.2346 on 20 degrees of freedom
    ## Multiple R-squared:  0.1705, Adjusted R-squared:  0.129 
    ## F-statistic:  4.11 on 1 and 20 DF,  p-value: 0.05616
