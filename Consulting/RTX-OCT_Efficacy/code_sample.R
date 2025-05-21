library(ggplot2)
library(ggcorrplot)
library(MASS)
library(dplyr)

# install.packages(c("survival", "survminer", "icenReg"))
# packages for time-to-event analysis
library(survival)
library(survminer)


# example uses IgG but analysis can also be applied to other lab values
# read data
dataMS = read.csv("client_data.csv", header = T)
# remove patients with no treatment label or no after treatment lab values
data_filtered = dataMS %>% filter(Length.of.Tx > 12, Treatment != "", !is.na(IgG.2))
data_filtered[data_filtered == ""] = NA # set empty cells to NA
for (j in 2:8) { # set labs with date but no lab value to NA
  month_col = paste0("Months.",j)
  lab_col = paste0("IgG.",j)
  for (i in 1:dim(data_filtered)[1]) {
    if (is.na(data_filtered[lab_col][i,1])) {
      data_filtered[month_col][i,1] = NA
    }
  }
}
data_filtered = data_filtered[-62,] # remove obs with conflicting timeline
# rename treatments
data_filtered$Treatment[data_filtered$Treatment=="O"] = "OCR"
data_filtered$Treatment[data_filtered$Treatment=="R"] = "RTX"

# Dataset columns are in the form
# Study.ID | Age | Sex | Income.Quantile | Treatment | Lab1.Date | Lab1.IgG | Lab2.Date | Lab2.IgG ...
# Here each row contains 1 patient's info and all lab values
# To do time-to-event analysis, we need to obtain the censored interval of each patient


## Reformat data for analysis and spaghetti plot 
id = c()
treatment = c()
months = c()
IgG = c()
IgG.Low = c() # IgG low indicator for each recorded lab
IgG.Low2 = c() # IgG low indicator for Cox model
left.interval = c() # left side of censored interval
right.interval = c() # right side of censored interval


# function to record and check if current patient-lab pair shows deficiency
deficiency_check = function(month_col, lab_col, i){
  if ((!is.na(data_filtered[month_col][i,1])) || (month_col == "Months.1")){
    id <<- c(id,data_filtered$Study.ID[i])
    treatment <<- c(treatment,data_filtered$Treatment[i])
    months <<- c(months,data_filtered[month_col][i,1])
    IgG <<- c(IgG,data_filtered[lab_col][i,1])
    IgG.Low <<- c(IgG.Low,data_filtered[lab_col][i,1]<7) # 7 is the IgG deficiency threshold
  }
}

# This code loops over all patients and their lab visits to record IgG deficiencies
for (i in 1:dim(data_filtered)[1]) { # loop over all patients
  IgG.Low_temp = 0
  left.interval_temp = 0
  right.interval_temp = 0
  for (j in 1:8) { # loop over all lab visits of a patient
    deficiency_check(paste0("Months.",j), paste0("IgG.",j), i)
    if (j>1){
      if (!is.na(data_filtered$IgG.1[i]) & (IgG.Low[length(IgG.Low)-1]==0) & 
          (IgG.Low[length(IgG.Low)]==1)){
        IgG.Low[length(IgG.Low)] = 3 # interval censored
        IgG.Low_temp = 3
        left.interval_temp = months[length(months)-1]
        right.interval_temp = months[length(months)]
      } else if ((j==2) & is.na(data_filtered$IgG.1[i]) & (IgG.Low[length(IgG.Low)]==1)) { 
        IgG.Low[length(IgG.Low)] = 3
        IgG.Low_temp = 3
        left.interval_temp = months[length(months)-1]
        right.interval_temp = months[length(months)]
      }
    }
  }
  IgG.Low2 = c(IgG.Low2,IgG.Low_temp)
  if (right.interval_temp == 0){ # no deficiency detected
    left.interval_temp = months[length(months)]
    right.interval_temp = Inf
  }
  left.interval = c(left.interval,left.interval_temp)
  right.interval = c(right.interval,right.interval_temp)
}


# data for spaghetti plot (each row is a patient-lab visit pair)
data_spaghetti = data.frame(id=id,treatment=treatment,months=months,IgG=IgG)
# spaghetti plot data view
#   id treatment months   IgG
# 1  8       RTX      0 10.00
# 2  8       RTX      7  9.80
# 3  8       RTX     12  9.20
# 4  8       RTX     25  8.56
# 5  9       RTX      0  9.40
# 6  9       RTX     47  7.84

# data reformatted for analysis (get rid of Lab value, Lab date columns and
# compute censored intervals)
data_reformatted = data.frame(age=data_filtered$Age, sex=data_filtered$Sex,
        income=data_filtered$Income.Quantile, treatment=data_filtered$Treatment,
        left.interval=left.interval, right.interval=right.interval, IgG.Low=IgG.Low2, 
        IgG.Low.2y=(right.interval<24)) # last col indicate deficiency within 2 years

# Reformatted dataset view
# IgG.Low: 0=right-censored, 3=interval-censored
#   age sex income treatment left.interval right.interval IgG.Low IgG.Low.2y
# 1  40   F      3       RTX            25            Inf       0      FALSE
# 2  51   M      3       RTX            59             62       3      FALSE
# 3  45   F      1       RTX            22            Inf       0      FALSE
# 4  48   M      2       RTX            30            Inf       0      FALSE
# 5  42   M      1       RTX            66            Inf       0      FALSE
# 6  40   F      2       RTX            23            Inf       0      FALSE









## Spaghetti plot
idx <- c()
for (i in data_filtered$Study.ID) { # remove patients with only 1 lab visit
  if (sum(!is.na(data_spaghetti$IgG[data_spaghetti$id == i])) > 1){
    idx <- c(idx,i)
  }
}
data_spaghetti = data_spaghetti[data_spaghetti$id %in% idx,]

# extract 15 patients from each treatment for spaghetti plot
set.seed(123)
id_OCR = sample(unique(data_spaghetti[data_spaghetti$treatment=="OCR",]$id),15)
id_RTX = sample(unique(data_spaghetti[data_spaghetti$treatment=="RTX",]$id),15)
idx = data_spaghetti$id %in% c(id_OCR,id_RTX)
p = ggplot(data = data_spaghetti[idx,], aes(x = months, y = IgG, group = id)) +
  geom_line(aes(color=treatment)) + geom_point(aes(color=treatment)) +
  geom_hline(yintercept=7) +
  xlab("Time (in months) since treatment began for each patient") +
  ylab("IgG level") + facet_grid(. ~ treatment) + theme_bw()
p  



## Logistic regression
# assume that age and sex are used along with treatment
res.logistic = glm(IgG.Low.2y ~ age + sex + treatment, 
                   data = data_reformatted, family = binomial)
summary(res.logistic)

# Call:
# glm(formula = IgG.Low.2y ~ age + sex + treatment, family = binomial, 
#       data = data_reformatted)
# 
# Coefficients:
#             Estimate Std. Error z value Pr(>|z|)    
# (Intercept) -5.41633    1.56838  -3.453 0.000553 ***
# age          0.04746    0.02957   1.605 0.108443    
# sexM        -0.74716    0.78950  -0.946 0.343960    
# treatmentR   0.09785    0.58301   0.168 0.866716    
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# (Dispersion parameter for binomial family taken to be 1)
# 
# Null deviance: 110.75  on 344  degrees of freedom
# Residual deviance: 107.21  on 341  degrees of freedom
# AIC: 115.21
# 
# Number of Fisher Scoring iterations: 6


## Time-to-event analysis
# package for interval censored data
library (icenReg)

# assume that age and sex are used along with treatment
res.cox = ic_sp(Surv(left.interval, right.interval, IgG.Low,type = "interval") ~ age + sex + treatment, 
                data = data_reformatted, bs_samples = 100)
summary(res.cox)

# Model:  Cox PH
# Dependency structure assumed: Independence
# Baseline:  semi-parametric 
# Call: ic_sp(formula = Surv(start, stop, IgG.Low, type = "interval") ~ 
#               age + sex + treatment, data = data_Cox, bs_samples = 100)
# 
#            Estimate Exp(Est) Std.Error z-value       p
# age         0.03954   1.0400   0.01707  2.3160 0.02058
# sexM       -0.58090   0.5594   0.40140 -1.4470 0.14790
# treatmentR -0.10690   0.8986   0.34610 -0.3089 0.75740
# 
# final llk =  -142.6641 
# Iterations =  15 
# Bootstrap Samples =  100

# coefficient interpretation: a positive sign means 
# that the hazard (risk of deficiency) is higher and vice versa