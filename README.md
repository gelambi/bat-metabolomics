# Untargeted metabolomics reveals fruit secondary metabolites alter bat nutrient absorption

This repository contains scripts, information, and figures related to the statistical analyses conducted for the research paper titled **"Untargeted metabolomics reveals fruit secondary metabolites alter bat nutrient absorption."** Using a mutualistic fruit bat (*Carollia perspicillata*), our research explores how four secondary metabolites (piperine, tannin acid, eugenol, and phytol) commonly found in plant tissues affect the foraging behavior and induce changes in the fecal metabolome. In this study, bats were captured and housed in flight cages. Nightly trials exposed them to varying concentrations of secondary metabolites. Objective 1 involved non-choice trials to measure food consumption, while Objective 2 evaluated the impact of metabolite consumption on the bat fecal metabolome. Fecal samples were collected, stored, and later analyzed to comprehend how secondary metabolites influence bat behavior and metabolism. All the analyses were performed in R v. 4.2.1.

## Scripts for objective 1. The effects of secondary metabolites on the foraging behavior of captive bats

### 1. `script1_objective1.R`

`script1_objective1.R` employs a Wilcoxon signed-rank test to assess the impact of four secondary metabolites on bat feeding preferences at three different concentrations. Raw consumption data is expressed as the amount of food (g) consumed in 30 minutes. Data is cleaned and outliers are addressed. The script calculates the ratio of consumption in the treatment group to the average consumption in the control group. Non-parametric Wilcoxon tests are conducted for each concentration, evaluating whether the observed differences in consumption ratios are statistically different from 1.

## Scripts for objective 2. The effects of secondary metabolites consumption on the bat fecal metabolome

### 1. `script2_objective2_metaMS.R`

`script2_objective2_metaMS.R` uses the `metaMS` package for chromatogram processing and metabolomic data analysis. The script initiates by configuring parameters for peak picking and grouping. Critical settings for matching peaks across samples are highlighted. The subsequent section involves reading and processing data files, utilizing the configured settings for peak identification. The resulting peak table and pseudospectra are saved for further analysis using the NIST library. Finally, the script calculates normalized instrument response by incorporating metadata and peak tables.

### 2. `script3_objective2_NMDS`

`script3_objective2_NMDS` explores the impact of four secondary metabolite consumption on composition of the bat fecal metabolome. For each concentration (0.1%, 2%, and 3%) a Non-Metric Multidimensional Scaling (NMDS) with subsequent visualization and statistical assessments were conducted. The NMDS results, stress plots, and ordination plots are generated for each concentration. Further, the script conducts PERMANOVA and BETADISP analyses to evaluate treatment effects. NMDS allows to visualize multivariate data patterns, PERMANOVA to assess treatment effects on sample groups, and BETADISP to evaluate multivariate homogeneity of group dispersion.

### 3. `script4_objective2_RF&GLMMs.R`

`script4_objective2_RF&GLMMs.R` focuses on exploring the impact of secondary metabolite consumption on the bat fecal metabolome through Random Forest (RF) and Boruta analyses. The script subsets the data based on different concentrations (0.1%, 2%, and 3%). The subsequent application of RF involves training models for each concentration subset to assess variable importance. The Boruta algorithm is employed for feature selection, identifying significant attributes among compounds.

Following the Boruta analysis, the script constructs Generalized Linear Mixed Models (GLMMs) to evaluate the significance of the identified compounds. In each model, the response variable aere the compound suggested in the Boruta analysis, and the predictor variable are the four different metabolites ingested by the bat at a given concentration. BatID was incorporated as a random effect, accounting for potential variations among individual bats.

Finally, the script uses the classyfireR package to obtain compound classification. The script extracts the InChIKeys from the dataset and creates a vector named InChI_Keys. Subsequently, it employs the purrr::map function to apply the get_classification function to each InChIKey, generating a list named Classification_List. This list holds the classification results for each compound.

## fecal_metabolome_CDF

Raw chromatograms in AIA Format (\*.CDF).

## Plots and Table folders

Contain different outputs from the analyses conducted.
