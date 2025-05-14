library(OlinkAnalyze)
library(dplyr)

npx_data1<-read.csv("//idnas11.d.uzh.ch/cdiazc$/HK postdoc/OLINKS/Proteomics of bariatric surgery/Yanghua inflammation/Final data long format/npx_long format_wt_covariates.csv", sep = ",")


npx_data1<-read.csv("//idnas11.d.uzh.ch/cdiazc$/HK postdoc/OLINKS/Proteomics of bariatric surgery/Yanghua metabolism/long format/Metabolism.csv", sep = ",")


npx_data1<-read.csv("//idnas11.d.uzh.ch/cdiazc$/HK postdoc/OLINKS/Proteomics of bariatric surgery/Yanghua cardiometabolic/cardiometabolic.csv", sep = ",")

npx_data1<-read.csv("//idnas11.d.uzh.ch/cdiazc$/HK postdoc/OLINKS/Proteomics of bariatric surgery/Yanghua CVD3/long format/CVD3.csv", sep = ",")

npx_data1<-read.csv("//idnas11.d.uzh.ch/cdiazc$/HK postdoc/OLINKS/Proteomics of bariatric surgery/Yangua CVD2/long format/CVD2.csv", sep = ",")




npx_data1$Time<- factor( npx_data1$Time, levels=c("Base", "1M", "3M","6M", "12M", "24M"))
npx_data1$Obesity<- factor( npx_data1$Obesity, levels=c("OW", "OB", "SO"))

anova_results <- olink_anova(df = npx_data1 , variable = c("Obesity","Time") , covariates = c( "Age","Gender","Type.of.surgery") )

#Filtering out significant and relevant results.
significant_assays <- anova_results %>%
  filter(Threshold == 'Significant' & term == "Obesity") %>%
  select(OlinkID) %>%
  distinct() %>%
  pull()



#Posthoc, all pairwise comparisons
anova_posthoc_results <- olink_anova_posthoc(npx_data1,
                                             variable=c("Time:Obesity"),
                                             covariates = c("Age","Type.of.surgery","Gender"),
                                             olinkid_list = significant_assays,
                                             
                                            effect = ("Time:Obesity"))







