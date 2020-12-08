# Comparisons of sensitivity rankings for heavy metals  

This file contains part of the R code for the article:

**Malaj, E.**, Grote M, Schäfer RB, Brack W, von der Ohe PC. (2012). Physiological sensitivity of freshwater macroinvertebrates to heavy metals. Environmental Toxicology and Chemistry; 31: 1754-1764. https://doi.org/10.1002/etc.1868 

## Depends

To run this analysis two files are require:  
`Org_File.csv` original data downloaded from ECOTOX (https://cfpub.epa.gov/ecotox/);  
`Outlier.csv` file with lab conditions to be removed.  

R version 3.6.1

Packages: lme4, reshape, doBy

## Approach

The aim of this work was to derive a unique heavy metal ranking for a representative set of bivalent metals (Cd, Cu, Cr, Ni, Pb, Zn, and Hg) using acute laboratory assays of aquatic invertebrate species. This parameter served as a basis for developing indices to investigate metal pollution. The goal here was not to investigate if toxicity of metals differs (which it does!), but to investigate whether there are differences after scaling. The lack of differences means that we can aggregate metals into a combined index.


## Data cleaning

Data was retrieved from US EPA ECOTOX (16,827 observations x 97 features) in 2011 (https://cfpub.epa.gov/ecotox/). Note that more recent versions of ECOTOX dataset might yield different results. ECOTOX is continuously updated as new results are published in peer-reviewed journals. This is an example of how to format and use a dataset in R.

Data wrangling included reshaping, summarizing, making new variables, and combining different datasets using R. Duplicate entries were removed, features were scaled, and independent variable log10-transformed. Where no mean values were reported, the average of the maximum and minimum values was taken. Extreme test conditions were removed, and each species needed to have toxicity data available for at least three metals. In total, 75.7% of the original data were omitted. Subsequently, 4,103 entries remained for seven selected heavy metals, with a diverse taxonomic range of 114 species, based on 412 different sources over the last 50 years. 

## GLMM

Exploratory data analysis included Pearson correlation coefficient (**Figure 1**) to examine the strength of relationships variables, model diagnostics (qq-plots, histograms) to check for normality and heterogeneity, outlier detection by using outlier tests and estimating Cook’s distance, and general plotting of parameters in base R. A generalized mixed effect regression algorithm (‘glmer’ in the R package lme4) was used to check for differences between metals groups under the condition of unbalanced design.

Results from the glmer showed that the scaled ranking were comparable across different species, therefore a Smetal ranking was developed and used in other studies (e.g., https://doi.org/10.1016/j.envpol.2017.05.017).


![Pair_Plot](https://user-images.githubusercontent.com/54320408/101551227-55013a00-3976-11eb-9c12-cb7a55b03ddc.png)

**Figure 1**: Pairwise mean comparison of physiological sensitivity values for each taxa of the seven heavy metals.






