Are rhizobia under selection to cheat?
================
Megan Frederickson
2019-04-20

Is there fitness conflict between legumes and rhizobia?
-------------------------------------------------------

This repository re-analyzes the data in:

Gano-Cohen KA, Wendlandt CE, Stokes PJ, Blanton MA, Quides KW, Zomorrodian A, Adinata ES, Sachs JL (2019) Interspecific conflict and the evolution of ineffective rhizobia. Ecology Letters. <https://doi.org/10.1111/ele.13247>

The authors make the case that in their legume-rhizobium study system, the rhizobia are selected to cheat. In other words, they report a negative correlation between legume and rhizobium fitnesses.

Here, I re-analyze their data using a standard genetic selection analysis approach, in which relativized fitness is regressed against standardized family-mean trait values.

I downloaded the data from Dryad on April 16, 2019. The citation for the data package is:

Gano-Cohen KA, Wendlandt CE, Stokes PJ, Blanton MA, Quides KW, Zomorrodian A, Adinata ES, Sachs JL (2019) Data from: Interspecific conflict and the evolution of ineffective rhizobia. Dryad Digital Repository. <https://doi.org/10.5061/dryad.cr65269>

First we need to read in the data, which is in three different tables in the Dryad package, and wrangle it in to a usable form.

``` r
data <- read_csv("Table_S4.csv") #Read in single inoculation experiment data
data$Strain <- as.factor(data$Strain) #Factors are factors and numbers are numbers
data$`Host Line` <- as.factor(data$`Host Line`)
data$Population <- as.factor(data$Population)
data$`Total nodules` <- as.numeric(data$`Total nodules`)
data$`Mean individual  nodule biomass (mg)` <- as.numeric(data$`Mean individual  nodule biomass (mg)`)
data$`Shoots mass (g)` <- as.numeric(data$`Shoots mass (g)`)
data$`Roots mass (g)` <- as.numeric(data$`Roots mass (g)`)
data$`Relative Growth` <- as.numeric(data$`Relative Growth`)
data$Block <- as.factor(data$Block)

#Ignore strain 156, just for interest sake
#data <- subset(data, Strain != 156)

gf.data <- read_csv("Table_S1.csv") #Read in genotype frequency data
gf.data$SI <- paste0(gf.data$`nodZ haplotype`, gf.data$`nolL haplotype`) #Concatenate SI haplotypes, as per paper
gf.data$CHR <- paste0(gf.data$`glnII Haplotype`, gf.data$`recA Haplotype`) #Concatenate CHR haplotypes, as per paper
gf.data <- subset(gf.data, `nolL haplotype` != "X") #Remove Xs, which are failures to amplify/sequence, as per Dryad data description
gf.data$plant <- gsub('_.*', "", gf.data$`Full Strain name`) #Make a column for unique plant ids
gf.data$plant <- toupper(gsub('R.*', "", gf.data$plant))  #Make plant ids upper case
gf.data <- gf.data[,c(1, 10, 2:9)] #Reorder columns
gf.data.long <- gather(gf.data, locus, haplotype, `glnII Haplotype`:CHR, factor_key=TRUE) #Make wide data into long format
gf.data.long <- subset(gf.data.long, haplotype != "n/an/a") #Remove NAs
gf.data.long <- subset(gf.data.long, locus == "SI" | locus == "CHR") #Prune to just CHR and SI haplotypes

#Summarize frequencies per plant
sum.gf <- gf.data.long %>% group_by(plant, `Plant Collection Site`, Year, locus, haplotype) %>% summarize(n = n())
tmp <- gf.data %>% group_by(plant) %>% summarize(total_n=n())
sum.gf <- merge(sum.gf, tmp, by="plant")

#Split by locus
sum.CHR <- subset(sum.gf, locus == "CHR")
sum.SI <- subset(sum.gf, locus == "SI")

#Make lists of which haplotypes are in each site
chr.haplos <- sum.CHR %>% group_by(`Plant Collection Site`, haplotype) %>% summarize(n=n())
SI.haplos <- sum.SI %>% group_by(`Plant Collection Site`, haplotype) %>% summarize(n=n())

#Initialize two empty dataframes, one for CHR and one for SI, to put all the data in 
full.chr <- data.frame(plant = character(),
                 `Plant Collection Site` = character(), 
                 Year = character(),
                 locus = character(),
                 haplotype = character(),
                 n = character(),
                 total_n = character(),
                 stringsAsFactors=FALSE) 

full.SI <- data.frame(plant = character(),
                 `Plant Collection Site` = character(), 
                 Year = character(),
                 locus = character(),
                 haplotype = character(),
                 n = character(),
                 total_n = character(),
                 stringsAsFactors=FALSE) 

#Add zeros for strains present in a population but not sampled in a particular plant
for(i in 1:length(unique(sum.CHR$plant))){
  
  plant.id = as.character(unique(sum.CHR$plant)[i])
  tmp2 <- subset(sum.CHR, plant == plant.id)
  tmp3 <- chr.haplos[chr.haplos$`Plant Collection Site` == tmp2[1,2], ]
  tmp4 <- data.frame(matrix(ncol = 7, nrow = length(setdiff(tmp3$haplotype, tmp2$haplotype))))
  colnames(tmp4) <- colnames(tmp2)  
  tmp4$plant <- plant.id
  tmp4$`Plant Collection Site` <- tmp2[1,2]
  tmp4$Year <- tmp2[1,3]
  tmp4$locus <- "CHR"
  tmp4$haplotype <- setdiff(tmp3$haplotype, tmp2$haplotype)
  tmp4$n <- 0
  tmp4$total_n = tmp2[1,7]
  tmp5 <- rbind(tmp2, tmp4)
  full.chr <- rbind(full.chr, tmp5)
}

for(i in 1:length(unique(sum.SI$plant))){
  
  plant.id = as.character(unique(sum.SI$plant)[i])
  tmp2 <- subset(sum.SI, plant == plant.id)
  tmp3 <- SI.haplos[SI.haplos$`Plant Collection Site` == tmp2[1,2], ]
  tmp4 <- data.frame(matrix(ncol = 7, nrow = length(setdiff(tmp3$haplotype, tmp2$haplotype))))
  colnames(tmp4) <- colnames(tmp2)  
  tmp4$plant <- plant.id
  tmp4$`Plant Collection Site` <- tmp2[1,2]
  tmp4$Year <- tmp2[1,3]
  tmp4$locus <- "SI"
  tmp4$haplotype <- setdiff(tmp3$haplotype, tmp2$haplotype)
  tmp4$n <- 0
  tmp4$total_n = tmp2[1,7]
  tmp5 <- rbind(tmp2, tmp4)
  full.SI <- rbind(full.SI, tmp5)
}
  
#Calculate frequencies
full.chr$freq <- full.chr$n/full.chr$total_n  
full.SI$freq <- full.SI$n/full.SI$total_n  

#Summarize frequences by population, haplotype, and locus
chr.sum <- full.chr %>% group_by(`Plant Collection Site`, haplotype) %>% summarize(mean_CHR_freq = mean(freq), n_CHR_freq = n(), 
                                                                                   sd_CHR_freq = sd(freq), se_CHR=sd_CHR_freq/sqrt(n_CHR_freq))
chr.sum$upper_CHR_CI <- chr.sum$mean_CHR_freq + (1.96*chr.sum$se_CHR/sqrt(chr.sum$n_CHR_freq))
chr.sum$lower_CHR_CI <- chr.sum$mean_CHR_freq - (1.96*chr.sum$se_CHR/sqrt(chr.sum$n_CHR_freq))
chr.sum$CHRpop <- paste0(chr.sum$`Plant Collection Site`, chr.sum$haplotype)
SI.sum <- full.SI %>% group_by(`Plant Collection Site`, haplotype) %>% summarize(mean_SI_freq = mean(freq), n_SI_freq = n(), 
                                                                                   sd_SI_freq = sd(freq), se_SI=sd_SI_freq/sqrt(n_SI_freq))
SI.sum$upper_SI_CI <- SI.sum$mean_SI_freq + (1.96*SI.sum$se_SI/sqrt(SI.sum$n_SI_freq))
SI.sum$lower_SI_CI <- SI.sum$mean_SI_freq - (1.96*SI.sum$se_SI/sqrt(SI.sum$n_SI_freq))
SI.sum[77,2] <- "Z62L79"  #Fix a typo in the original table
SI.sum$SIpop <- paste0(SI.sum$`Plant Collection Site`, SI.sum$haplotype)

#Match strain names to genotypes
strain_names <- read_csv("Table_S2.csv")
strain_names <- strain_names[-1, c(-7, -9, -11, -13, -18)]

#Fix population names
strain_names$Population <- ifelse(strain_names$Population == "Bodega Marine Reserve", "BMR",
                                  ifelse(strain_names$Population == "Griffith Park", "GRI",
                                  ifelse(strain_names$Population == "Robert J. Bernard Biological Field Station", "CLA",
                                  ifelse(strain_names$Population == "University of California Riverside", "UCR", 
                                  ifelse(strain_names$Population == "Anza Borrego Desert State Park", "ANZ", "YUC")))))

#Remove underscores from haplotype names to merge between datasets
strain_names$`chromosome haplotype (CHR)` <- gsub("_", "", strain_names$`chromosome haplotype (CHR)`)
strain_names$`symbiosis island  haplotye (SI)` <- gsub("_", "", strain_names$`symbiosis island  haplotye (SI)`)
strain_names$SIpop <- paste0(strain_names$Population, strain_names$`symbiosis island  haplotye (SI)`)
strain_names$CHRpop <- paste0(strain_names$Population, strain_names$`chromosome haplotype (CHR)`)

#Merge datasets
full.geno <- merge(strain_names, SI.sum, by="SIpop", all.x = TRUE)
full.geno <- merge(full.geno, chr.sum, by = "CHRpop", all.x =TRUE)

#Calculate nodules and plants sampled per population
sum <- sum.gf %>% group_by(`Plant Collection Site`, locus) %>% summarize(total_nods_sampled=sum(n), total_plants_sampled=length(unique(plant)))
```

How many nodules and plants were sampled per site?
==================================================

The sampling is uneven, with few plants sampled in ANZ and YUC.

``` r
table <- sum
colnames(table) <- c("Plant Collection Site", "Locus", "Nodules sampled (no.)", "Plants sampled (no.)")
kable(table)
```

| Plant Collection Site | Locus |  Nodules sampled (no.)|  Plants sampled (no.)|
|:----------------------|:------|----------------------:|---------------------:|
| ANZ                   | SI    |                     43|                     4|
| ANZ                   | CHR   |                     44|                     4|
| BMR                   | SI    |                    107|                    16|
| BMR                   | CHR   |                    136|                    16|
| CLA                   | SI    |                     68|                    20|
| CLA                   | CHR   |                     68|                    20|
| GRI                   | SI    |                      4|                     4|
| GRI                   | CHR   |                     68|                    18|
| UCR                   | SI    |                     17|                     5|
| UCR                   | CHR   |                     88|                    31|
| YUC                   | SI    |                     14|                     2|
| YUC                   | CHR   |                     38|                     7|

Plot sampling of CHR and SI frequencies
=======================================

``` r
lm1 <- lm(as.numeric(`CHR genotype frequency`)~n_CHR_freq, data=full.geno)
```

    ## Warning in eval(predvars, data, env): NAs introduced by coercion

``` r
summary(lm1)
```

    ## 
    ## Call:
    ## lm(formula = as.numeric(`CHR genotype frequency`) ~ n_CHR_freq, 
    ##     data = full.geno)
    ## 
    ## Residuals:
    ##      Min       1Q   Median       3Q      Max 
    ## -0.15864 -0.06893 -0.02852  0.03546  0.45500 
    ## 
    ## Coefficients:
    ##              Estimate Std. Error t value Pr(>|t|)   
    ## (Intercept)  0.209254   0.058857   3.555  0.00161 **
    ## n_CHR_freq  -0.006979   0.003098  -2.253  0.03368 * 
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 0.1422 on 24 degrees of freedom
    ##   (4 observations deleted due to missingness)
    ## Multiple R-squared:  0.1746, Adjusted R-squared:  0.1402 
    ## F-statistic: 5.075 on 1 and 24 DF,  p-value: 0.03368

``` r
CHR <- ggplot(data=full.geno, aes(y=as.numeric(`CHR genotype frequency`), x=n_CHR_freq, color=Population))+
       geom_jitter(width=0.25)+
        geom_smooth(method="lm", se=FALSE, color=1)+
        xlab("Plants sampled (no.)")+
        ylab("CHR genotype frequency")+
        geom_text(aes(label=`Inoculation #`),hjust=0, vjust=0, size=2.5)
CHR
```

    ## Warning in FUN(X[[i]], ...): NAs introduced by coercion

    ## Warning in FUN(X[[i]], ...): NAs introduced by coercion

    ## Warning in FUN(X[[i]], ...): NAs introduced by coercion

    ## Warning in FUN(X[[i]], ...): NAs introduced by coercion

    ## Warning: Removed 4 rows containing non-finite values (stat_smooth).

    ## Warning: Removed 4 rows containing missing values (geom_point).

    ## Warning: Removed 4 rows containing missing values (geom_text).

![](README_files/figure-markdown_github/Data%20distributions-1.png)

``` r
lm2 <- lm(as.numeric(`SI genotype frequency`)~n_SI_freq, data=full.geno)
```

    ## Warning in eval(predvars, data, env): NAs introduced by coercion

``` r
summary(lm2)
```

    ## 
    ## Call:
    ## lm(formula = as.numeric(`SI genotype frequency`) ~ n_SI_freq, 
    ##     data = full.geno)
    ## 
    ## Residuals:
    ##      Min       1Q   Median       3Q      Max 
    ## -0.15576 -0.06492  0.01539  0.07204  0.19793 
    ## 
    ## Coefficients:
    ##              Estimate Std. Error t value Pr(>|t|)    
    ## (Intercept)  0.207921   0.032609   6.376 2.54e-06 ***
    ## n_SI_freq   -0.007491   0.002757  -2.717   0.0129 *  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 0.09656 on 21 degrees of freedom
    ##   (7 observations deleted due to missingness)
    ## Multiple R-squared:  0.2601, Adjusted R-squared:  0.2249 
    ## F-statistic: 7.383 on 1 and 21 DF,  p-value: 0.01291

``` r
SI <- ggplot(data=full.geno, aes(y=as.numeric(`SI genotype frequency`), x=as.numeric(n_SI_freq), color=Population))+
        geom_smooth(method="lm", se=FALSE, color=1)+
        geom_jitter(width = 0.25)+
        xlab("Plants sampled (no.)")+
        ylab("SI genotype frequency")+
        geom_text(aes(label=`Inoculation #`),hjust=0, vjust=0, size=2.5)

fig1 <- plot_grid(CHR, SI, nrow=2, labels="auto")
```

    ## Warning in FUN(X[[i]], ...): NAs introduced by coercion

    ## Warning in FUN(X[[i]], ...): NAs introduced by coercion

    ## Warning in FUN(X[[i]], ...): NAs introduced by coercion

    ## Warning in FUN(X[[i]], ...): NAs introduced by coercion

    ## Warning: Removed 4 rows containing non-finite values (stat_smooth).

    ## Warning: Removed 4 rows containing missing values (geom_point).

    ## Warning: Removed 4 rows containing missing values (geom_text).

    ## Warning in FUN(X[[i]], ...): NAs introduced by coercion

    ## Warning in FUN(X[[i]], ...): NAs introduced by coercion

    ## Warning in FUN(X[[i]], ...): NAs introduced by coercion

    ## Warning in FUN(X[[i]], ...): NAs introduced by coercion

    ## Warning: Removed 7 rows containing non-finite values (stat_smooth).

    ## Warning: Removed 7 rows containing missing values (geom_point).

    ## Warning: Removed 7 rows containing missing values (geom_text).

``` r
save_plot("Fig1.png", fig1, base_width=8, base_height=8)
```

Calculate genotype means
========================

This dataset has two measures of plant fitness and four measures of rhizobium fitness. The two measures of plant fitness are total plant biomass (i.e., shoot plus root mass) and relative growth rate (i.e., inoculated plant biomass divided by control plant biomass, hereafter "RGR"). The four measures of rhizobium fitness are total number of nodules, mean nodule mass, and the two genotypic frequencies, which the authors abbreviate CHR and SI. See original paper for details.

First, we need to calculate the appropriate strain means for each variable. I think the authors did this just by taking the average. I did this using the emmeans package, which calculates estimated marginal means from a linear mixed model. Otherwise, I have tried to stick as close as possible to the analysis presented in the paper, which states that "all effects were coded as fixed except block, which was treated as a random effect" and that variables were log-transformed as needed to improve normality. I also calculated strain means only on data from sympatric hosts, again following the approach in the paper.

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

#Set dataset to use
df <- data_sym

#Fit models 
lmm1 <- lmer(log_nodules~Population + Strain + `Host Line` + (1|Block), data=df) #Total nodule number model
plot(lmm1)
```

![](README_files/figure-markdown_github/Genotype%20means-1.png)

``` r
summary(lmm1)
```

    ## Linear mixed model fit by REML ['lmerMod']
    ## Formula: log_nodules ~ Population + Strain + `Host Line` + (1 | Block)
    ##    Data: df
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
df.means <- as.data.frame(emmeans(lmm1, "Strain")) #Extract line means and save them to new dataframe
colnames(df.means) <- c("Strain", "Population", "tot_nod_lsmean", "tot_nod_SE", "tot_nod_df", "tot_nod_lowCL", "tot_nod_highCL") #Fix column names

lmm2 <- lmer(log_mean_nod_biomass~Population + Strain + `Host Line` + (1|Block), data=df) # Mean nodule mass model
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

lmm3 <- lmer(plant_biomass~Population + Strain  + `Host Line` + (1|Block), data=df) #Plant biomass model
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

lmm4 <- lmer(log(`Relative Growth`)~Population + Strain  + `Host Line` + (1|Block), data=df) #Relative growth rate, log-transformed as per the paper, but note this generates NAs because there are negative values
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

df.means <- merge(df.means, full.geno, by.x="Strain", by.y = "Inoculation #") #Add CHR and SI values for each isolate
```

Relative fitness within populations and standardize trait values
================================================================

Ideally, to compare across studies and fitness measures, we should calculate selection gradients in the standard way, as we would for any continuous phenotype. Here, the phenotype of interest is mutualist quality, or the amount of fitness benefit that a rhizobium strain provides to its legume host. The authors provide two measures of this trait in their dataset: plant biomass and RGR. Generally, trait values are standardized by subtracting the mean and dividing by the population standard deviation, and I did the same in my re-analysis of the data. The authors also report four fitness proxies in their dataset, namely the total number of nodules, mean nodule mass, and the two genotype frequencies CHR and SI. Fitness is typically relativized by dividing by the population mean, allowing comparisons of the strength of selection across analyses. In all cases, trait values and fitnesses are genotype means, i.e., breeding values (although this is an odd term when applied to bacteria). This is how Porter and Simms (2014) analyzed their data, when they found evidence for fitness conflict in a different legume-rhizobium mutualism, and we could thus directly compare the results.

``` r
#Calculate means and SDs for each fitness/trait measure for each population
tmp <- df.means %>% group_by(Population.x) %>% summarize(
  pop_mean_CHR=mean(as.numeric(`CHR genotype frequency`), na.rm=TRUE),
  pop_mean_SI=mean(as.numeric(`SI genotype frequency`), na.rm=TRUE)
  )
```

    ## Warning in mean(as.numeric(`SI genotype frequency`), na.rm = TRUE): NAs
    ## introduced by coercion

``` r
df.means <- merge(df.means, tmp, by="Population.x") #Merge data frames

#Using global means and standard deviations for traits measured in common garden
#Relativize fitness by dividing by the mean fitness, the global mean for fitness measures from the common garden but the pop mean 
#for CHR and SI
df.means$tot_nod_std <- df.means$tot_nod_lsmean/mean(df.means$tot_nod_lsmean)
df.means$nod_mass_std <- df.means$nod_mass_lsmean/mean(df.means$nod_mass_lsmean)
df.means$CHR_std <- as.numeric(df.means$`CHR genotype frequency`)/df.means$pop_mean_CHR
df.means$SI_std <- as.numeric(df.means$`SI genotype frequency`)/df.means$pop_mean_SI
```

    ## Warning: NAs introduced by coercion

``` r
#Standardize traits by subtracting the mean and dividing by the standard deviation
df.means$plant_biomass_std <- (df.means$plant_biomass_lsmean - mean(df.means$plant_biomass_lsmean))/sd(df.means$plant_biomass_lsmean)
df.means$RGR_std <- (df.means$RGR_lsmean - mean(df.means$RGR_lsmean))/sd(df.means$RGR_lsmean)
```

Estimate selection gradients
============================

Next, I adopt standard genetic selection analyses to estimate selection gradients by regressing relativized fitness against standardized trait values. I estimated eight selection gradients, because there are four fitness measures (nodule number, nodule mass, CHR, and SI) and two traits (plant biomass and RGR) in the paper.

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
    ## -0.54068 -0.04724  0.00767  0.07151  0.17623 
    ## 
    ## Coefficients:
    ##                   Estimate Std. Error t value Pr(>|t|)    
    ## (Intercept)       1.000000   0.027888  35.857   <2e-16 ***
    ## plant_biomass_std 0.001843   0.028441   0.065    0.949    
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 0.1422 on 24 degrees of freedom
    ## Multiple R-squared:  0.000175,   Adjusted R-squared:  -0.04148 
    ## F-statistic: 0.004201 on 1 and 24 DF,  p-value: 0.9489

``` r
S1 <- paste0("S = ", round(summary(lm1)$coefficients[2,1], 3))
pval1 <- round(summary(lm1)$coefficients[2,4],3)
linetype <- ifelse(pval1 <= 0.05, "solid", "dashed")
pval1 <- ifelse(pval1 <= 0.05, paste0("p = ", pval1), "ns")
S1 <- paste0(S1, ", ", pval1)

p1 <- ggplot(data=df.means, aes(y=tot_nod_std, x=plant_biomass_std, color=Population.x))+
      geom_point()+
      geom_smooth(method="lm", se=FALSE, linetype=linetype, color=1)+
      geom_text(aes(label=Strain),hjust=0, vjust=0, size=3, check_overlap=TRUE)+
      ylab("Nodule number")+
      xlab("Plant biomass")

#Nodule mass and plant biomass

lm2 <- lm(nod_mass_std~plant_biomass_std, data=df.means) #Model
summary(lm2)
```

    ## 
    ## Call:
    ## lm(formula = nod_mass_std ~ plant_biomass_std, data = df.means)
    ## 
    ## Residuals:
    ##     Min      1Q  Median      3Q     Max 
    ## -0.4023 -0.1287  0.0164  0.1508  0.3692 
    ## 
    ## Coefficients:
    ##                   Estimate Std. Error t value Pr(>|t|)    
    ## (Intercept)        1.00000    0.04347  23.005   <2e-16 ***
    ## plant_biomass_std -0.06764    0.04433  -1.526     0.14    
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 0.2216 on 24 degrees of freedom
    ## Multiple R-squared:  0.08844,    Adjusted R-squared:  0.05046 
    ## F-statistic: 2.328 on 1 and 24 DF,  p-value: 0.1401

``` r
S2 <- paste0("S = ", round(summary(lm2)$coefficients[2,1], 3))
pval2 <- round(summary(lm2)$coefficients[2,4],3)
linetype <- ifelse(pval2 <= 0.05, "solid", "dashed")
pval2 <- ifelse(pval2 <= 0.05, paste0("p = ", pval2), "ns")
S2 <- paste0(S2, ", ", pval2)


p2 <- ggplot(data=df.means, aes(y=nod_mass_std, x=plant_biomass_std, color=Population.x))+
      geom_point()+
      geom_smooth(method="lm", se=FALSE, linetype=linetype, color=1)+
      geom_text(aes(label=Strain),hjust=0, vjust=0, size=3)+
      ylab("Nodule mass")+
      xlab("Plant biomass")

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
    ## -1.0135 -0.7299 -0.4674  0.5305  2.3455 
    ## 
    ## Coefficients:
    ##                   Estimate Std. Error t value Pr(>|t|)    
    ## (Intercept)         1.0000     0.2098   4.766 7.52e-05 ***
    ## plant_biomass_std  -0.2259     0.2140  -1.056    0.302    
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 1.07 on 24 degrees of freedom
    ## Multiple R-squared:  0.04437,    Adjusted R-squared:  0.004552 
    ## F-statistic: 1.114 on 1 and 24 DF,  p-value: 0.3017

``` r
S3 <- paste0("S = ", round(summary(lm3)$coefficients[2,1], 3))
pval3 <- round(summary(lm3)$coefficients[2,4],3)
linetype <- ifelse(pval3 <= 0.05, "solid", "dashed")
pval3 <- ifelse(pval3 <= 0.05, paste0("p = ", pval3), "ns")
S3 <- paste0(S3, ", ", pval3)


p3 <- ggplot(data=df.means, aes(y=CHR_std, x=plant_biomass_std, color=Population.x))+
      geom_point()+
      geom_smooth(method="lm", se=FALSE, linetype=linetype, color=1)+
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
    ## -0.78329 -0.58917 -0.02558  0.35574  1.81858 
    ## 
    ## Coefficients:
    ##                   Estimate Std. Error t value Pr(>|t|)    
    ## (Intercept)        0.99967    0.15099   6.621 9.37e-07 ***
    ## plant_biomass_std  0.01658    0.15173   0.109    0.914    
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 0.7548 on 23 degrees of freedom
    ##   (1 observation deleted due to missingness)
    ## Multiple R-squared:  0.0005188,  Adjusted R-squared:  -0.04294 
    ## F-statistic: 0.01194 on 1 and 23 DF,  p-value: 0.9139

``` r
S4 <- paste0("S = ", round(summary(lm4)$coefficients[2,1], 3))
pval4 <- round(summary(lm4)$coefficients[2,4],3)
linetype <- ifelse(pval4 <= 0.05, "solid", "dashed")
pval4 <- ifelse(pval4 <= 0.05, paste0("p = ", pval4), "ns")
S4 <- paste0(S4, ", ", pval4)

p4 <- ggplot(data=df.means, aes(y=SI_std, x=plant_biomass_std, color=Population.x))+
      geom_point()+
      geom_smooth(method="lm", se=FALSE, linetype=linetype, color=1)+
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
    ## -0.54833 -0.03214  0.01904  0.06870  0.17767 
    ## 
    ## Coefficients:
    ##             Estimate Std. Error t value Pr(>|t|)    
    ## (Intercept)  1.00000    0.02775  36.034   <2e-16 ***
    ## RGR_std     -0.01389    0.02830  -0.491    0.628    
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 0.1415 on 24 degrees of freedom
    ## Multiple R-squared:  0.009942,   Adjusted R-squared:  -0.03131 
    ## F-statistic: 0.241 on 1 and 24 DF,  p-value: 0.6279

``` r
S5 <- paste0("S = ", round(summary(lm5)$coefficients[2,1], 3))
pval5 <- round(summary(lm5)$coefficients[2,4],3)
linetype <- ifelse(pval5 <= 0.05, "solid", "dashed")
pval5 <- ifelse(pval5 <= 0.05, paste0("p = ", pval5), "ns")
S5 <- paste0(S5, ", ", pval5)

p5 <- ggplot(data=df.means, aes(y=tot_nod_std, x=RGR_std, color=Population.x))+
      geom_point()+
      geom_smooth(method="lm", se=FALSE, linetype=linetype, color=1)+
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
    ## -0.35871 -0.16105  0.01479  0.16774  0.40330 
    ## 
    ## Coefficients:
    ##             Estimate Std. Error t value Pr(>|t|)    
    ## (Intercept)  1.00000    0.04464  22.402   <2e-16 ***
    ## RGR_std     -0.04475    0.04552  -0.983    0.335    
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 0.2276 on 24 degrees of freedom
    ## Multiple R-squared:  0.03871,    Adjusted R-squared:  -0.001349 
    ## F-statistic: 0.9663 on 1 and 24 DF,  p-value: 0.3354

``` r
S6 <- paste0("S = ", round(summary(lm6)$coefficients[2,1], 3))
pval6 <- round(summary(lm6)$coefficients[2,4],3)
linetype <- ifelse(pval6 <= 0.05, "solid", "dashed")
pval6 <- ifelse(pval6 <= 0.05, paste0("p = ", pval6), "ns")
S6 <- paste0(S6, ", ", pval6)


p6 <- ggplot(data=df.means, aes(y=nod_mass_std, x=RGR_std, color=Population.x))+
      geom_point()+
      geom_smooth(method="lm", se=FALSE, linetype=linetype, color=1)+
      geom_text(aes(label=Strain),hjust=0, vjust=0, size=3)+
      ylab("Nodule mass")+
      xlab("RGR")

#CHR and RGR
lm7 <- lm(CHR_std~RGR_std, data=df.means)
summary(lm7)
```

    ## 
    ## Call:
    ## lm(formula = CHR_std ~ RGR_std, data = df.means)
    ## 
    ## Residuals:
    ##     Min      1Q  Median      3Q     Max 
    ## -1.0128 -0.6893 -0.5369  0.6664  2.4992 
    ## 
    ## Coefficients:
    ##             Estimate Std. Error t value Pr(>|t|)    
    ## (Intercept)   1.0000     0.2117   4.724 8.37e-05 ***
    ## RGR_std      -0.1768     0.2159  -0.819    0.421    
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 1.079 on 24 degrees of freedom
    ## Multiple R-squared:  0.02718,    Adjusted R-squared:  -0.01335 
    ## F-statistic: 0.6706 on 1 and 24 DF,  p-value: 0.4209

``` r
S7 <- paste0("S = ", round(summary(lm7)$coefficients[2,1], 3))
pval7 <- round(summary(lm7)$coefficients[2,4],3)
linetype <- ifelse(pval7 <= 0.05, "solid", "dashed")
pval7 <- ifelse(pval7 <= 0.05, paste0("p = ", pval7), "ns")
S7 <- paste0(S7, ", ", pval7)


p7 <- ggplot(data=df.means, aes(y=CHR_std, x=RGR_std, color=Population.x))+
      geom_point()+
      geom_smooth(method="lm", se=FALSE, linetype=linetype, color=1)+
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
    ## -0.77954 -0.60289 -0.00088  0.31093  1.83186 
    ## 
    ## Coefficients:
    ##              Estimate Std. Error t value Pr(>|t|)    
    ## (Intercept)  1.000137   0.151025   6.622 9.33e-07 ***
    ## RGR_std     -0.006645   0.151839  -0.044    0.965    
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 0.755 on 23 degrees of freedom
    ##   (1 observation deleted due to missingness)
    ## Multiple R-squared:  8.325e-05,  Adjusted R-squared:  -0.04339 
    ## F-statistic: 0.001915 on 1 and 23 DF,  p-value: 0.9655

``` r
S8 <- paste0("S = ", round(summary(lm8)$coefficients[2,1], 3))
pval8 <- round(summary(lm8)$coefficients[2,4],3)
linetype <- ifelse(pval8 <= 0.05, "solid", "dashed")
pval8 <- ifelse(pval8 <= 0.05, paste0("p = ", pval8), "ns")
S8 <- paste0(S8, ", ", pval8)


p8 <- ggplot(data=df.means, aes(y=SI_std, x=RGR_std, color=Population.x))+
      geom_point()+
      geom_smooth(method="lm", se=FALSE, linetype=linetype, color=1)+
      geom_text(aes(label=Strain),hjust=0, vjust=0, size=3)+
      ylab("SI")+
      xlab("RGR")
```

Selection on symbiont quality, measured as RGR
==============================================

After standardizing relative growth within each population and relativizing all four bacterial fitness measures, there is no significant selection on symbiont quality when it is measured as RGR, as in the paper. Dashed lines show non-significant regressions.

``` r
plot_grid(p5,p6,p7,p8, labels=c(S5, S6, S7, S8), hjust=c(-1,-1,-0.76,-0.76), scale = 0.9)
```

![](README_files/figure-markdown_github/Visualize%20selection%20gradients,%20using%20RGR-1.png)

Selection on symbiont quality, measured as plant biomass
========================================================

Performing the same analysis but measuring symbiont quality using plant biomass, in the more traditional fashion, shows all non-significant relationships too. Again, dashed lines show non-signficiant selection gradients.

``` r
plot_grid(p1,p2,p3,p4, labels=c(S1, S2, S3, S4), hjust = c(-1,-1,-0.76,-0.84), scale=0.9)
```

![](README_files/figure-markdown_github/Visualize%20selection%20gradients,%20using%20plant%20biomass-1.png)

Re-create original analyses
===========================

``` r
original <- df %>% group_by(Population, Strain) %>% summarize(log_mean_RGR = mean(log(as.numeric(`Relative Growth`))))
```

    ## Warning in log(as.numeric(`Relative Growth`)): NaNs produced

    ## Warning in log(as.numeric(`Relative Growth`)): NaNs produced

    ## Warning in log(as.numeric(`Relative Growth`)): NaNs produced

    ## Warning in log(as.numeric(`Relative Growth`)): NaNs produced

    ## Warning in log(as.numeric(`Relative Growth`)): NaNs produced

    ## Warning in log(as.numeric(`Relative Growth`)): NaNs produced

``` r
original <- merge(original, df.means, by="Strain")

lm.orig1 <- lm(log_mean_RGR~as.numeric(`CHR genotype frequency`), data=original)
summary(lm.orig1)
```

    ## 
    ## Call:
    ## lm(formula = log_mean_RGR ~ as.numeric(`CHR genotype frequency`), 
    ##     data = original)
    ## 
    ## Residuals:
    ##     Min      1Q  Median      3Q     Max 
    ## -1.0192 -0.3374  0.0497  0.3109  0.7403 
    ## 
    ## Coefficients:
    ##                                      Estimate Std. Error t value Pr(>|t|)
    ## (Intercept)                            6.6881     0.1242  53.841   <2e-16
    ## as.numeric(`CHR genotype frequency`)  -1.4293     0.6420  -2.226    0.039
    ##                                         
    ## (Intercept)                          ***
    ## as.numeric(`CHR genotype frequency`) *  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 0.4748 on 18 degrees of freedom
    ##   (6 observations deleted due to missingness)
    ## Multiple R-squared:  0.2159, Adjusted R-squared:  0.1724 
    ## F-statistic: 4.957 on 1 and 18 DF,  p-value: 0.03899

``` r
orig.fig3a <- ggplot(data=original, aes(x=log_mean_RGR, y=as.numeric(`CHR genotype frequency`)))+
              geom_point()
orig.fig3a
```

    ## Warning: Removed 6 rows containing missing values (geom_point).

![](README_files/figure-markdown_github/Original%20analyses-1.png)

``` r
lm.orig2 <- lm(log_mean_RGR~as.numeric(`SI genotype frequency`), data=original)
```

    ## Warning in eval(predvars, data, env): NAs introduced by coercion

``` r
summary(lm.orig2)
```

    ## 
    ## Call:
    ## lm(formula = log_mean_RGR ~ as.numeric(`SI genotype frequency`), 
    ##     data = original)
    ## 
    ## Residuals:
    ##      Min       1Q   Median       3Q      Max 
    ## -1.06441 -0.27457  0.04902  0.38369  0.66552 
    ## 
    ## Coefficients:
    ##                                     Estimate Std. Error t value Pr(>|t|)
    ## (Intercept)                            6.763      0.194  34.862   <2e-16
    ## as.numeric(`SI genotype frequency`)   -1.655      1.301  -1.272     0.22
    ##                                        
    ## (Intercept)                         ***
    ## as.numeric(`SI genotype frequency`)    
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 0.5166 on 17 degrees of freedom
    ##   (7 observations deleted due to missingness)
    ## Multiple R-squared:  0.08692,    Adjusted R-squared:  0.03321 
    ## F-statistic: 1.618 on 1 and 17 DF,  p-value: 0.2205

``` r
orig.fig3b <- ggplot(data=original, aes(x=log_mean_RGR, y=as.numeric(`SI genotype frequency`)))+
              geom_point()
orig.fig3b
```

    ## Warning in FUN(X[[i]], ...): NAs introduced by coercion

    ## Warning in FUN(X[[i]], ...): NAs introduced by coercion

    ## Warning: Removed 7 rows containing missing values (geom_point).

![](README_files/figure-markdown_github/Original%20analyses-2.png)

``` r
lm.orig1 <- lm(log_mean_RGR~as.numeric(`CHR genotype frequency`), data=subset(original, Population != "YUC"))
summary(lm.orig1)
```

    ## 
    ## Call:
    ## lm(formula = log_mean_RGR ~ as.numeric(`CHR genotype frequency`), 
    ##     data = subset(original, Population != "YUC"))
    ## 
    ## Residuals:
    ##      Min       1Q   Median       3Q      Max 
    ## -0.98094 -0.31051  0.07014  0.30382  0.63437 
    ## 
    ## Coefficients:
    ##                                      Estimate Std. Error t value Pr(>|t|)
    ## (Intercept)                            6.6361     0.1241   53.48   <2e-16
    ## as.numeric(`CHR genotype frequency`)  -0.2300     0.8851   -0.26    0.798
    ##                                         
    ## (Intercept)                          ***
    ## as.numeric(`CHR genotype frequency`)    
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 0.4555 on 16 degrees of freedom
    ##   (4 observations deleted due to missingness)
    ## Multiple R-squared:  0.004202,   Adjusted R-squared:  -0.05804 
    ## F-statistic: 0.06751 on 1 and 16 DF,  p-value: 0.7983

``` r
orig.fig3a <- ggplot(data=original, aes(x=log_mean_RGR, y=as.numeric(`CHR genotype frequency`)))+
              geom_point()
orig.fig3a
```

    ## Warning: Removed 6 rows containing missing values (geom_point).

![](README_files/figure-markdown_github/Original%20analyses-3.png)

``` r
lm.orig2 <- lm(log_mean_RGR~as.numeric(`SI genotype frequency`), data=subset(original, Population != "YUC"))
```

    ## Warning in eval(predvars, data, env): NAs introduced by coercion

``` r
summary(lm.orig2)
```

    ## 
    ## Call:
    ## lm(formula = log_mean_RGR ~ as.numeric(`SI genotype frequency`), 
    ##     data = subset(original, Population != "YUC"))
    ## 
    ## Residuals:
    ##     Min      1Q  Median      3Q     Max 
    ## -1.0014 -0.3476  0.1017  0.3093  0.6146 
    ## 
    ## Coefficients:
    ##                                     Estimate Std. Error t value Pr(>|t|)
    ## (Intercept)                          6.65756    0.17417  38.225  2.3e-16
    ## as.numeric(`SI genotype frequency`) -0.06913    1.32289  -0.052    0.959
    ##                                        
    ## (Intercept)                         ***
    ## as.numeric(`SI genotype frequency`)    
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 0.4508 on 15 degrees of freedom
    ##   (5 observations deleted due to missingness)
    ## Multiple R-squared:  0.000182,   Adjusted R-squared:  -0.06647 
    ## F-statistic: 0.002731 on 1 and 15 DF,  p-value: 0.959

``` r
orig.fig3b <- ggplot(data=original, aes(x=log_mean_RGR, y=as.numeric(`SI genotype frequency`)))+
              geom_point()
orig.fig3b
```

    ## Warning in FUN(X[[i]], ...): NAs introduced by coercion

    ## Warning in FUN(X[[i]], ...): NAs introduced by coercion

    ## Warning: Removed 7 rows containing missing values (geom_point).

![](README_files/figure-markdown_github/Original%20analyses-4.png)
