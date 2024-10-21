SES & the microbiome 3.1: Beta Diversity and Permanova; CORRECTION 2021
================
Ruth CE Bowyer
2024-10-21

# **0. About**

This files contains the correctly specified adjusted models for measures
of gut microbiome dissimilarity. The univariate models with just the key
variable of interest (income, education or IMD) have not been re-run in
this file as they were correctly specified

This file was initially ran in 2021 and was subsequently updated for
github publication(see ReadMe for more info)

# **1. Data and set up**

## **Libraries**

``` r
library(phyloseq)
library(data.table)
```

    ## Warning: package 'data.table' was built under R version 4.3.3

``` r
library(ggplot2)
library(ape)
library(vegan)
```

    ## Warning: package 'vegan' was built under R version 4.3.3

    ## Loading required package: permute

    ## Loading required package: lattice

    ## This is vegan 2.6-8

## **Phyloseq object creation**

``` r
#dd is the root directory
otus <- read.table(paste0(dd, "Variance_transformed_SEStrimmed5%OTU_table.txt"), sep = " ", row.names = 1, header = T) 

tree <- read.tree(paste0(dd,"microbiome_tables/phyloseq_tables_aug2017/reduced_tree.tre"))
tax <- read.csv(paste0(dd,"taxa_SES_trimmed5%.csv"), row.names = 1)
mapping <- read.csv(paste0(dd,"mapping_ses_correctpcdes041218.csv"))

mapping$IMD5f <- as.factor(mapping$IMD5f)
mapping$eduff <- as.factor(mapping$eduff)
mapping$Income4F <- as.factor(mapping$Income4F)
mapping$FIsqrt <- sqrt(mapping$FI)

OTU <- otu_table(otus, taxa_are_rows = T)
tax <- as.matrix(tax)
TAX <- tax_table(tax)

physeq <- phyloseq(OTU, TAX)
mapping <- sample_data(mapping)
row.names(mapping) <- mapping$SequencingSpecificName
physeq1 <- merge_phyloseq(physeq, mapping, tree)
```

``` r
metadata <- as(sample_data(physeq1), "data.frame")
metadata1 <- metadata[complete.cases(metadata$FIsqrt, metadata$Age, metadata$HEI, metadata$BMI),] 
map <- sample_data(physeq1)
bray <- vegdist(t(otus), "bray") 
bray <- as.dist(as(bray, "matrix"))
Ord.bray <- ordinate(physeq1, method="MDS", distance=bray)

braym <- as.matrix(bray)
braym <- as.data.frame(braym)
bray1 <- braym[which(row.names(braym)%in%metadata1$SequencingSpecificName),]
bray1 <- bray1[,which(names(bray1)%in%metadata1$SequencingSpecificName)]
bray1<- as.dist(as(bray1, "matrix"))

set.seed(1)

a1 <- adonis(bray1~IMD5f, data=metadata1, permutations = 5000)
```

    ## Warning: 'adonis' is deprecated.
    ## Use 'adonis2' instead.
    ## See help("Deprecated") and help("vegan-deprecated").

``` r
a1$aov.tab
```

    ## Permutation: free
    ## Number of permutations: 5000
    ## 
    ## Terms added sequentially (first to last)
    ## 
    ##             Df SumsOfSqs  MeanSqs F.Model      R2   Pr(>F)   
    ## IMD5f        4     0.174 0.043505  1.3729 0.00328 0.007798 **
    ## Residuals 1667    52.823 0.031687         0.99672            
    ## Total     1671    52.997                  1.00000            
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

# **2. Models**

## **IMD - Bray**

There was a misspecification of a model relating to the adjusted results
reported as page 5 in the original publication, referenced on page 5 of
the manuscript as ‘adjusted NPMANOVA for Bray–Curtis dissimilarity… and
weighted UniFrac’

The model was incorrectly specified as follows (not run):

``` r
adonis(bray1 ~ IMD5f + FI + Age + BMI + HEI + Library_sizelog10, data = metadata1,
```

However, as the variables are added sequentially in this model, the key
factor should have been included as the last variable (as in the
subsequent model ouputs below)

### **The correctly specified model - IMD/Bray**

``` r
imd.a <- adonis(bray1 ~ FI + Age + BMI + HEI + Library_sizelog10 + IMD5f, data = metadata1, permutations = 5000)

imd.a$aov.tab
```

    ## Permutation: free
    ## Number of permutations: 5000
    ## 
    ## Terms added sequentially (first to last)
    ## 
    ##                     Df SumsOfSqs MeanSqs F.Model      R2  Pr(>F)    
    ## FI                   1     0.278 0.27849  8.9942 0.00525 0.00020 ***
    ## Age                  1     0.215 0.21507  6.9461 0.00406 0.00020 ***
    ## BMI                  1     0.239 0.23861  7.7064 0.00450 0.00020 ***
    ## HEI                  1     0.144 0.14414  4.6552 0.00272 0.00020 ***
    ## Library_sizelog10    1     0.518 0.51801 16.7302 0.00977 0.00020 ***
    ## IMD5f                4     0.142 0.03562  1.1503 0.00269 0.09998 .  
    ## Residuals         1662    51.460 0.03096         0.97100            
    ## Total             1671    52.997                 1.00000            
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

If one wanted to understand the relative contribution of each factor to
the outcome (ie a marginal model) this is also offered within the vegan
package and can give interesting insight

``` r
imd.m <- adonis2(bray1 ~ FI + Age + BMI + HEI + Library_sizelog10 + IMD5f, by="margin",
        data = metadata1, permutations = 5000)
imd.m
```

    ## Permutation test for adonis under reduced model
    ## Marginal effects of terms
    ## Permutation: free
    ## Number of permutations: 5000
    ## 
    ## adonis2(formula = bray1 ~ FI + Age + BMI + HEI + Library_sizelog10 + IMD5f, data = metadata1, permutations = 5000, by = "margin")
    ##                     Df SumOfSqs      R2       F Pr(>F)    
    ## FI                   1    0.188 0.00356  6.0858 0.0002 ***
    ## Age                  1    0.221 0.00417  7.1433 0.0002 ***
    ## BMI                  1    0.199 0.00376  6.4384 0.0002 ***
    ## HEI                  1    0.141 0.00266  4.5480 0.0002 ***
    ## Library_sizelog10    1    0.515 0.00971 16.6204 0.0002 ***
    ## IMD5f                4    0.142 0.00269  1.1503 0.1050    
    ## Residual          1662   51.460 0.97100                   
    ## Total             1671   52.997 1.00000                   
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

Below are the correctly specified models only as they should have been
included in the paper.

## **IMD - W. Unifrac**

``` r
metadata1 <- metadata[complete.cases(metadata$FIsqrt, metadata$Age, metadata$HEI, metadata$BMI),] 
w.unifrac <- read.table(paste0(dd, "distance_matrixes/weigthed_unifrac_phyloseq_rooted_tree.txt"))

wum <- as.matrix(w.unifrac)
wum <- as.data.frame(wum)
wum1 <- wum[which(row.names(wum)%in%metadata1$SequencingSpecificName),] 
wum1 <- wum1[,which(names(wum1)%in%metadata1$SequencingSpecificName)] 
wum1<- as.dist(as(wum1, "matrix"))
```

### **The correctly specified model - IMD/W. UniFrac**

``` r
full.wum.IMD.adonis <- adonis(wum1 ~ FI + Age + BMI + HEI + Library_sizelog10 + IMD5f, 
                              data = metadata1, permutations = 5000) 
full.wum.IMD.adonis$aov.tab
```

    ## Permutation: free
    ## Number of permutations: 5000
    ## 
    ## Terms added sequentially (first to last)
    ## 
    ##                     Df SumsOfSqs  MeanSqs F.Model      R2 Pr(>F)    
    ## FI                   1    0.1168 0.116785 10.8532 0.00635 0.0002 ***
    ## Age                  1    0.0974 0.097393  9.0511 0.00530 0.0002 ***
    ## BMI                  1    0.0977 0.097746  9.0839 0.00532 0.0002 ***
    ## HEI                  1    0.0463 0.046320  4.3047 0.00252 0.0002 ***
    ## Library_sizelog10    1    0.0980 0.098039  9.1111 0.00533 0.0002 ***
    ## IMD5f                4    0.0481 0.012028  1.1178 0.00262 0.2104    
    ## Residuals         1662   17.8838 0.010760         0.97257           
    ## Total             1671   18.3882                  1.00000           
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

## **Education - Bray**

``` r
metadata2 <- metadata[complete.cases(metadata$FIsqrt, metadata$Age, metadata$HEI, metadata$eduff),]
bray2 <- braym[which(row.names(braym)%in%metadata2$SequencingSpecificName),] 
bray2 <- bray2[,which(names(bray2)%in%metadata2$SequencingSpecificName)] 
bray2<- as.dist(as(bray2, "matrix"))
```

### **The correctly specified model**

``` r
full.bray.eduff.adonis <- adonis(bray2 ~ FI + Age + BMI + HEI + Library_sizelog10 + eduff, 
                                 data = metadata2, permutations = 5000)

full.bray.eduff.adonis$aov.tab
```

    ## Permutation: free
    ## Number of permutations: 5000
    ## 
    ## Terms added sequentially (first to last)
    ## 
    ##                     Df SumsOfSqs MeanSqs F.Model      R2   Pr(>F)    
    ## FI                   1     0.278 0.27761  8.9865 0.00616 0.000200 ***
    ## Age                  1     0.163 0.16311  5.2802 0.00362 0.000200 ***
    ## BMI                  1     0.158 0.15840  5.1275 0.00352 0.000200 ***
    ## HEI                  1     0.123 0.12273  3.9729 0.00272 0.000200 ***
    ## Library_sizelog10    1     0.416 0.41623 13.4737 0.00924 0.000200 ***
    ## eduff                3     0.137 0.04575  1.4808 0.00305 0.006999 ** 
    ## Residuals         1417    43.774 0.03089         0.97169             
    ## Total             1425    45.049                 1.00000             
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

## **Education - W.unifrac**

``` r
wum2 <- wum[which(row.names(wum)%in%metadata2$SequencingSpecificName),] 
wum2 <- wum2[,which(names(wum2)%in%metadata2$SequencingSpecificName)] 
wum2<- as.dist(as(wum2, "matrix"))
```

### **The correctly specified model - Education/W. Unifrac**

``` r
full.wum.eduff.adonis <- adonis(wum2 ~ FI + Age + BMI + HEI + Library_sizelog10 + eduff, 
                                 data = metadata2, permutations = 5000)

full.wum.eduff.adonis$aov.tab
```

    ## Permutation: free
    ## Number of permutations: 5000
    ## 
    ## Terms added sequentially (first to last)
    ## 
    ##                     Df SumsOfSqs  MeanSqs F.Model      R2    Pr(>F)    
    ## FI                   1    0.1096 0.109594 10.2292 0.00703 0.0002000 ***
    ## Age                  1    0.0724 0.072445  6.7618 0.00465 0.0002000 ***
    ## BMI                  1    0.0654 0.065392  6.1035 0.00419 0.0002000 ***
    ## HEI                  1    0.0397 0.039667  3.7024 0.00254 0.0003999 ***
    ## Library_sizelog10    1    0.0775 0.077482  7.2320 0.00497 0.0002000 ***
    ## eduff                3    0.0445 0.014828  1.3840 0.00285 0.0425915 *  
    ## Residuals         1417   15.1814 0.010714         0.97376              
    ## Total             1425   15.5905                  1.00000              
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

## **Income - Bray**

``` r
metadata3 <- metadata[complete.cases(metadata$FIsqrt, metadata$Age, metadata$HEI, metadata$BMI, metadata$Income4F),]

bray3 <- braym[which(row.names(braym)%in%metadata3$SequencingSpecificName),] 
bray3 <- bray3[,which(names(bray3)%in%metadata3$SequencingSpecificName)] 
bray3<- as.dist(as(bray3, "matrix"))
```

### **The correctly specified model**

``` r
full.bray.inc.adonis <- adonis(bray3 ~ FI + Age + BMI + HEI + Library_sizelog10 + Income4F, 
                                 data = metadata3, permutations = 5000)

full.bray.inc.adonis$aov.tab
```

    ## Permutation: free
    ## Number of permutations: 5000
    ## 
    ## Terms added sequentially (first to last)
    ## 
    ##                    Df SumsOfSqs  MeanSqs F.Model      R2 Pr(>F)    
    ## FI                  1    0.1820 0.182014  5.8371 0.00716 0.0002 ***
    ## Age                 1    0.0803 0.080287  2.5748 0.00316 0.0002 ***
    ## BMI                 1    0.0953 0.095337  3.0574 0.00375 0.0002 ***
    ## HEI                 1    0.0704 0.070393  2.2575 0.00277 0.0020 ** 
    ## Library_sizelog10   1    0.2409 0.240898  7.7255 0.00948 0.0002 ***
    ## Income4F            3    0.1008 0.033600  1.0775 0.00397 0.2464    
    ## Residuals         790   24.6338 0.031182         0.96970           
    ## Total             798   25.4035                  1.00000           
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

## **Income - W.Unifrac**

``` r
wum3 <- wum[which(row.names(wum)%in%metadata3$SequencingSpecificName),] 
wum3 <- wum3[,which(names(wum3)%in%metadata3$SequencingSpecificName)] 
wum3<- as.dist(as(wum3, "matrix"))
```

### **The correctly specified model**

``` r
full.wum.inc.adonis <- adonis(wum3 ~ FI + Age + BMI + HEI + Library_sizelog10 + Income4F, 
                                 data = metadata3, permutations = 5000)

full.wum.inc.adonis$aov.tab
```

    ## Permutation: free
    ## Number of permutations: 5000
    ## 
    ## Terms added sequentially (first to last)
    ## 
    ##                    Df SumsOfSqs  MeanSqs F.Model      R2   Pr(>F)    
    ## FI                  1    0.0683 0.068274  6.3189 0.00778 0.000200 ***
    ## Age                 1    0.0237 0.023676  2.1913 0.00270 0.009798 ** 
    ## BMI                 1    0.0359 0.035868  3.3197 0.00409 0.000200 ***
    ## HEI                 1    0.0234 0.023427  2.1683 0.00267 0.007598 ** 
    ## Library_sizelog10   1    0.0531 0.053133  4.9177 0.00606 0.000200 ***
    ## Income4F            3    0.0333 0.011108  1.0281 0.00380 0.375725    
    ## Residuals         790    8.5356 0.010805         0.97291             
    ## Total             798    8.7733                  1.00000             
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

# **3. Packages Used**

``` r
citation("phyloseq")
```

    ## To cite phyloseq in publications, or otherwise credit, please use:
    ## 
    ##   phyloseq: An R package for reproducible interactive analysis and
    ##   graphics of microbiome census data. Paul J. McMurdie and Susan Holmes
    ##   (2013) PLoS ONE 8(4):e61217.
    ## 
    ## A BibTeX entry for LaTeX users is
    ## 
    ##   @Article{,
    ##     author = {Paul J. McMurdie and Susan Holmes},
    ##     journal = {PLoS ONE},
    ##     pages = {e61217},
    ##     title = {phyloseq: An R package for reproducible interactive analysis and graphics of microbiome census data},
    ##     volume = {8},
    ##     number = {4},
    ##     year = {2013},
    ##     url = {http://dx.plos.org/10.1371/journal.pone.0061217},
    ##   }

``` r
citation("vegan")
```

    ## To cite package 'vegan' in publications use:
    ## 
    ##   Oksanen J, Simpson G, Blanchet F, Kindt R, Legendre P, Minchin P,
    ##   O'Hara R, Solymos P, Stevens M, Szoecs E, Wagner H, Barbour M,
    ##   Bedward M, Bolker B, Borcard D, Carvalho G, Chirico M, De Caceres M,
    ##   Durand S, Evangelista H, FitzJohn R, Friendly M, Furneaux B, Hannigan
    ##   G, Hill M, Lahti L, McGlinn D, Ouellette M, Ribeiro Cunha E, Smith T,
    ##   Stier A, Ter Braak C, Weedon J (2024). _vegan: Community Ecology
    ##   Package_. R package version 2.6-8,
    ##   <https://CRAN.R-project.org/package=vegan>.
    ## 
    ## A BibTeX entry for LaTeX users is
    ## 
    ##   @Manual{,
    ##     title = {vegan: Community Ecology Package},
    ##     author = {Jari Oksanen and Gavin L. Simpson and F. Guillaume Blanchet and Roeland Kindt and Pierre Legendre and Peter R. Minchin and R.B. O'Hara and Peter Solymos and M. Henry H. Stevens and Eduard Szoecs and Helene Wagner and Matt Barbour and Michael Bedward and Ben Bolker and Daniel Borcard and Gustavo Carvalho and Michael Chirico and Miquel {De Caceres} and Sebastien Durand and Heloisa Beatriz Antoniazi Evangelista and Rich FitzJohn and Michael Friendly and Brendan Furneaux and Geoffrey Hannigan and Mark O. Hill and Leo Lahti and Dan McGlinn and Marie-Helene Ouellette and Eduardo {Ribeiro Cunha} and Tyler Smith and Adrian Stier and Cajo J.F. {Ter Braak} and James Weedon},
    ##     year = {2024},
    ##     note = {R package version 2.6-8},
    ##     url = {https://CRAN.R-project.org/package=vegan},
    ##   }

``` r
## For creation of Rmarkdown files
citation("rmarkdown")
```

    ## To cite package 'rmarkdown' in publications use:
    ## 
    ##   Allaire J, Xie Y, Dervieux C, McPherson J, Luraschi J, Ushey K,
    ##   Atkins A, Wickham H, Cheng J, Chang W, Iannone R (2024). _rmarkdown:
    ##   Dynamic Documents for R_. R package version 2.28,
    ##   <https://github.com/rstudio/rmarkdown>.
    ## 
    ##   Xie Y, Allaire J, Grolemund G (2018). _R Markdown: The Definitive
    ##   Guide_. Chapman and Hall/CRC, Boca Raton, Florida. ISBN
    ##   9781138359338, <https://bookdown.org/yihui/rmarkdown>.
    ## 
    ##   Xie Y, Dervieux C, Riederer E (2020). _R Markdown Cookbook_. Chapman
    ##   and Hall/CRC, Boca Raton, Florida. ISBN 9780367563837,
    ##   <https://bookdown.org/yihui/rmarkdown-cookbook>.
    ## 
    ## To see these entries in BibTeX format, use 'print(<citation>,
    ## bibtex=TRUE)', 'toBibtex(.)', or set
    ## 'options(citation.bibtex.max=999)'.

``` r
citation("kableExtra")
```

    ## To cite package 'kableExtra' in publications use:
    ## 
    ##   Zhu H (2024). _kableExtra: Construct Complex Table with 'kable' and
    ##   Pipe Syntax_. R package version 1.4.0,
    ##   <https://CRAN.R-project.org/package=kableExtra>.
    ## 
    ## A BibTeX entry for LaTeX users is
    ## 
    ##   @Manual{,
    ##     title = {kableExtra: Construct Complex Table with 'kable' and Pipe Syntax},
    ##     author = {Hao Zhu},
    ##     year = {2024},
    ##     note = {R package version 1.4.0},
    ##     url = {https://CRAN.R-project.org/package=kableExtra},
    ##   }

``` r
citation("knitr")
```

    ## To cite package 'knitr' in publications use:
    ## 
    ##   Xie Y (2024). _knitr: A General-Purpose Package for Dynamic Report
    ##   Generation in R_. R package version 1.48, <https://yihui.org/knitr/>.
    ## 
    ##   Yihui Xie (2015) Dynamic Documents with R and knitr. 2nd edition.
    ##   Chapman and Hall/CRC. ISBN 978-1498716963
    ## 
    ##   Yihui Xie (2014) knitr: A Comprehensive Tool for Reproducible
    ##   Research in R. In Victoria Stodden, Friedrich Leisch and Roger D.
    ##   Peng, editors, Implementing Reproducible Computational Research.
    ##   Chapman and Hall/CRC. ISBN 978-1466561595
    ## 
    ## To see these entries in BibTeX format, use 'print(<citation>,
    ## bibtex=TRUE)', 'toBibtex(.)', or set
    ## 'options(citation.bibtex.max=999)'.
