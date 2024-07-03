library(survival)
library(survminer)

#Read OS and PFS individual patient data
os <- read.csv("IPD_both_drug_os.csv")
pfs <- read.csv("IPD_both_drug_pfs.csv")


#create survival objects, which are standard objects compatible with many of the time-to-event functions in R
s_os <- Surv(time = os$time, event = os$status, type = "right")
s_pfs <- Surv(time = pfs$time, event = pfs$status, type = "right")

#Kaplan-Meier curves separated by drug
#The survfit function takes the survival object". The ~ is used to define the variable by which the curves should be grouped
km_fit_os <- survfit(s_os ~ os$drug)
km_fit_pfs <- survfit(s_pfs ~ pfs$drug)

#Plot the KM curve - OS
ggsurvplot(fit = km_fit_os, 
           data = os,
           xlab = "Months",
           ylab = "Survival Probability",
           ggtheme = theme_minimal(),
           pval = TRUE)  # Add p-value for the log-rank test

#Plot the KM curve - PFS
ggsurvplot(fit = km_fit_pfs, 
           data = pfs,
           xlab = "Months",
           ylab = "Progression Free Survival Probability",
           ggtheme = theme_minimal(),
           pval = TRUE)  # Add p-value for the log-rank test