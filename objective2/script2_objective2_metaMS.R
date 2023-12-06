#################################################################################################
### Objective 2. The effects of secondary metabolites consumption on the bat fecal metabolome ###
#################################################################################################

rm(list=ls())

library(metaMS)
library(tidyverse)
library(dplyr)

# Installation of metaMS

# if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
# BiocManager::install("metaMS")
# if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
# BiocManager::install("metaMSdata")

################################
### Chromatogram processing ###
################################

## Tutorial: https://rdrr.io/bioc/metaMS/f/inst/doc/runGC.pdf

# Set parameters 
agilent.GC <-
  metaMSsettings("protocolName" = "agilent_single_quad.GC",          #// Name of the protocol, here is a Thermo Scientific Instrument with a tripple quadrupole detector...
                 "chrom" = "GC",               #// defines the kind of analisys LC vs GC...
                 PeakPicking = list(           #// Calls the function "PeakPicking" from the package "XCMS"...
                   method = "matchedFilter",   #// Defines what METHOD or algorithim will be used by the function "PeakPicking" to find peaks across samples...
                   #// 8 methods are available:
                   # "centWave" = for LCMS "Peak density and wavelet based feature detection for high resolution LC/MS data in centroid mode. This algorithm is most suitable for high resolution LC/{TOF,OrbiTrap,FTICR}-MS data in centroid mode."
                   # "centWaveWithPredictedIsotopeROIs"  = for LCMS "Peak density and wavelet based feature detection for high resolution LC/MS data in centroid mode with additional peak picking of isotope features on basis of isotope peak predictions"
                   # "massifquant" = for LCMS "Massifquant is a Kalman filter (KF) based feature detection for XC-MS data in centroid mode (currently in experimental stage)"
                   # "matchedFilter_orig" = old version of "matchedFilter"
                   # "MS1" = for tandem mass spec XML data
                   # "MSW" = for tandem mass spec XML data
                   # "addPredictedIsotopeFeatures" = for LCMS, Peak density and wavelet based feature detection aiming at isotope peaks for high resolution LC/MS data in centroid mode
                   # "matchedFilter"
                   #// binSize = numeric, specifying the width of the bins/slices in m/z dimension.
                   #// impute = Specify the method to be used for missing values("none" (no linear interpolation), "lin" (linear interpolation), "linbase" (linear interpolation within a certain bin-neighborhood) and "intlin". See imputeLinInterpol for more details.).
                   #// baseValue = The base value to which empty elements should be set. This is only considered for method = "linbase" from the impute setting above.
                   #// distance = Also for method "linbase": number of non-empty neighboring element that should be considered for linear interpolation.
                   #// fwhm = numeric, specify the full width at half maximum height of matched filtration gaussian model peak. Only used to calculate the actual sigma, see below.
                   #// sigma = fwhm/2.3548 // numeric, specifying the standard deviation (width) of the matched filtration model peak. Iether claculated with the formula provided or empirically found from the data...
                   #// max = numeric, representing the maximum number of peaks that are expected/will be identified per m/z slice. FIVE is the default and might be too small for the work we do (e.g. sesquiterpenoids).
                   #// snthresh = numeric, defining the signal to noise cutoff to be used in the chromatographic peak detection step. You can play with this number and a sunset of your data to find the ideal setting...
                   #// steps = numeric, defining the number of consecutive BINS to be MERGED into a single peak before filtration.
                   #// mzdiff = numeric, defining the minimum difference in m/z for peaks with overlapping retention times
                   #// index = logical, specifying whether indicies should be returned instead of values for m/z and retention times.
                   
                   step = 0.5,                 #// Same as above
                   steps = 2, 
                   mzdiff = .8, 
                   fwhm = 5,                 #// this is Full Width at Half Maxima. this means the width of the peaks at half their maximun hight/singal/intensity.abundance...
                   snthresh = 1,
                   max = 75),
                 CAMERA = list(perfwhm = 1))

###

metaSetting(agilent.GC, "DBconstruction") <- list(         #// This is are the setting to make the list of pseudospectra. Not terribly important I think...
  minintens = 0.0,
  rttol = .1,
  intensityMeasure = "maxo",
  DBthreshold = .80,
  minfeat = 5)
metaSetting(agilent.GC, "match2DB") <- list(                #// Setting to match to the data base of internal library. Not important unless you have the library...
  simthresh = 0.80,
  timeComparison = "rt",
  rtdiff = .5,
  RIdiff = 5,
  minfeat = 2)
metaSetting(agilent.GC, "matchIrrelevants") <- list(       #// Same as above...
  irrelevantClasses = c("Bleeding", "Plasticizers"),
  timeComparison = "rt",
  RIdiff = 2,
  rtdiff = .05,
  simthresh = 0.70)
metaSetting(agilent.GC, "betweenSamples") <- list(       #// THIS MATCHES THE PEAKS ACROSS THE SAMPLES AND THIS IS VERY IMPORTANT. THE CRITICAL SETTINGS HERE ARE THE RTDIFF AND SIMTHRESH.
  min.class.fraction = .01,                              #// The "rtdiff" is the maximum retention time difference between samples allowed to still be considered the same compound...
  min.class.size = 2,        # ****!!!!!!!!!****        #// The simthresh os the minimum similarity value in the fragmentation specta to still be considered the same compounds... I suggest to relax this a bit compare to the default...
  timeComparison = "rt",                                #// min.class.size = the minimum number of samples that need to have this compounds in order to be reported. E.G.: if =3 then, if a compound is only found in 2 samples it will be disregarded as a contamination or artifact.
  rtdiff = .2,        # ****!!!!!!!!!****               #// min.class.fraction = fraction of samples in which a compound is present before it is regarded as an unknown. my recomendation, control this factor with "min.class.size" and keep "min.class.fraction" as small as possible.
  RIdiff = 2,                                           #// timeComparison = choose how to compare retention time, via raw retention time "rt" or via standardized retention index "RI"...
  simthresh = .80)  

###

# check settings object...
agilent.GC
metaSetting(agilent.GC, "PeakPicking") 
metaSetting(agilent.GC, "betweenSamples") 

###

# Now, read data
# Files were converter to CDF in ChemStation

mzml_fecalmetabolome_dir <- "~/Desktop/bat_fecal_metabolome/fecal_metabolome_CDF"
mzml_fecalmetabolome_dir
mzml_fecalmetabolome_files <- list.files(mzml_fecalmetabolome_dir, full.names = TRUE, ignore.case = TRUE)
mzml_fecalmetabolome_files #check that all the names of all files are now in this list
result_all_fecalmetabolome <- runGC(files = mzml_fecalmetabolome_files, settings = agilent.GC, DB = NULL, findUnknowns = TRUE,  bpnworkers: 6)
all_fecalmetabolome_comp_table <- result_all_fecalmetabolome$PeakTable
str(all_fecalmetabolome_comp_table)
all_fecalmetabolome_comp_spectra <- result_all_fecalmetabolome$PseudoSpectra
all_fecalmetabolome_comp_spectra

write.csv(result_all_fecalmetabolome$PeakTable, file = "PeakTable_fecalmetabolome_1snthresh_simthresh80.csv")
write.msp(all_fecalmetabolome_comp_spectra, file = "fecalmetabolome_1snthresh_simthresh.80.msp") ### read this in the NIST search

######################################################
### Calculations of Normalized instrument response ###
######################################################

#After ID of compounds
metadata <- read.csv("metadata.csv")
head(metadata)
peak <- read.csv("peaktable_BS_dec2.csv")
head(peak)

data <- merge(metadata, peak, by = "sampleID") # merge the metadata with the peak table
head(data)

### Expressed instrument response as Normalized instrument response 
data_standard <- (data[ , 7:302])/data$Ribitol
data_final <- (data_standard)/data$weight

# Add the columns from the metadata
data_final$treatment <- data$treatment
data_final$batID <- data$batID
data_final$concentration <- data$concentration
data_final$data <- data$date
data_final <- na.omit(data_final)
data_final$concentration <- as.factor(data_final$concentration)
write.csv(data_final, file = "data_final_metaMS_dec2_normalizedinsresponse.csv")

####################################
####################################
####################################
####################################

