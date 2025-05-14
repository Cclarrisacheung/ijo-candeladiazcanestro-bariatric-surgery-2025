library(lme4)
library(lmerTest)

library(OlinkAnalyze)

npx_data1<-read.csv("//idnas11.d.uzh.ch/cdiazc$/HK postdoc/OLINKS/Proteomics of bariatric surgery/Yanghua proteomics/npx_long format_wt_inflammation.csv", sep = ";")

npx_data1<-read.csv("//idnas11.d.uzh.ch/cdiazc$/HK postdoc/OLINKS/Proteomics of bariatric surgery/Yanghua metabolism/long format/Metabolism.csv", sep = ";")

npx_data1<-read.csv("//idnas11.d.uzh.ch/cdiazc$/HK postdoc/OLINKS/Proteomics of bariatric surgery/Yanghua cardiometabolic/long format.csv", sep = ";")

npx_data1<-read.csv("//idnas11.d.uzh.ch/cdiazc$/HK postdoc/OLINKS/Proteomics of bariatric surgery/Yanghua CVD3/long format/CVD3.csv", sep = ";")

npx_data1<-read.csv("//idnas11.d.uzh.ch/cdiazc$/HK postdoc/OLINKS/Proteomics of bariatric surgery/Yangua CVD2/long format/CVD2.csv", sep = ";")


npx_data1$Time<- factor( npx_data1$Time, levels=c("Baseline", "1 month", "3 months","6 months", "12 months", "24 months"))

library(dplyr)


# Results in model NPX~Time*Treatment+(1|Subject)+(1|Site) 

lmer_results <- olink_lmer(df = npx_data1, variable="Time",covariates=c("Gender","Age","Baseline.Weight","Baseline.BMI"), random = 'Subject', verbose = TRUE)

lmer_results <- olink_lmer(df = npx_data1, variable=c("Time","Type.of.surgery"),covariates=c("Gender","Age","Baseline.Weight","Baseline.BMI"), random = 'Subject', verbose = FALSE)

assay_list <- lmer_results %>%
  filter(Threshold == 'Significant' & term == 'Time') %>%
  select(OlinkID) %>%
  distinct() %>%
  pull()




results_lmer_posthoc <- olink_lmer_posthoc(df = npx_data1,
                                           olinkid_list = assay_list,
                                           model_formula = "NPX~Time+Gender+Age+Baseline.Weight+Baseline.BMI+(1|Subject)",
                                           effect_formula = "pairwise~Time",
                                           verbose = TRUE)

assay_list


write.csv(results_lmer_posthoc ,"C:/Users/USUARIO/Desktop/lmer CVD2.csv")
write.csv(lmer_results ,"C:/Users/USUARIO/Desktop/lmer inflammation.csv")

# Reorder following a precise order

library(ggplot2)
library(dplyr)



assay_list<- c('OID01130')


plots <- olink_lmer_plot(df = npx_data1,variable=c("Time","Type.of.surgery"),covariates=c("Gender","Age","Baseline.Weight","Baseline.BMI"),
                         random = 'Subject',
                         x_axis_variable = 'Time',
                         col_variable = 'Type.of.surgery',
                         verbose=TRUE,
                         olinkid_list = assay_list,
                         number_of_proteins_per_plot = 10)

plots

