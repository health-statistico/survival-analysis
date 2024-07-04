library(survival)
library(survminer)

#Read OS and PFS individual patient data
#These datasets contain right-censoring."status"=1 denotes death/progression, "status"=0 denotes censoring.
os <- read.csv("IPD_both_drug_os.csv")
pfs <- read.csv("IPD_both_drug_pfs.csv")


#create survival objects, which are standard objects compatible with many of the time-to-event functions in R
s_os <- Surv(time = os$time, event = os$status, type = "right")
s_pfs <- Surv(time = pfs$time, event = pfs$status, type = "right")

#Generate Kaplan-Meier curves separated by drug
#The survfit function takes the survival object". The ~ is used to define the variable by which the curves should be grouped
km_fit_os <- survfit(s_os ~ os$drug)
km_fit_pfs <- survfit(s_pfs ~ pfs$drug)

#Plot the KM curve - OS
ggsurvplot(fit = km_fit_os, 
           data = os,
           xlab = "Months",
           ylab = "Overall Survival Probability",
           ggtheme = theme_minimal(),
           pval = TRUE)  # Add p-value for the log-rank test

#Plot the KM curve - PFS
ggsurvplot(fit = km_fit_pfs, 
           data = pfs,
           xlab = "Months",
           ylab = "Progression Free Survival Probability",
           ggtheme = theme_minimal(),
           pval = TRUE)  # Add p-value for the log-rank test

#Plot the log-log survival curves for OS to check if the two groups have proportional hazards. If both curves are approximately parallel then the proportional hazards assumtion holds
ggsurvplot(fit = km_fit_os, 
           data = os,
           fun = "cloglog",
           xlab = "Months",
           ylab = "Complementary Log-Log Survival",
           ggtheme = theme_minimal())  # Add p-value for the log-rank test

#Plot log-log PFS curves
ggsurvplot(fit = km_fit_pfs, 
           data = pfs,
           fun = "cloglog",
           xlab = "Months",
           ylab = "Complementary Log-Log PFS",
           ggtheme = theme_minimal())  # Add p-value for the log-rank test

#Fit the cox proportional hazards model to the data. Here we use the drug variable as the only covariate in the cox regression
cox_fit_os <- coxph(s_os ~ os$drug)
summary(cox_fit_os)
#Test the proportional hazards assumption using the Schoenfeld Residuals Test for proportional hazards on the cox regression results. Small p-values indicate a violation of the proportional hazards assumption.
cox.zph(cox_fit_os)
#You may also visually inspect the Schoenfeld residuals like so
ggcoxzph(cox.zph(cox_fit_os))

#now do the same for PFS
cox_fit_pfs <- coxph(s_pfs ~ pfs$drug)
summary(cox_fit_pfs)
cox.zph(cox_fit_pfs)
ggcoxzph(cox.zph(cox_fit_pfs))
