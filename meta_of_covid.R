############################################################
# Read in May 30 2020.
# 
# Author: 
# Yijia Weng
# 
# Environment:
# - use R 4.1.0 
# 
# Input:
# - Analysis_1.xlsx
# - Analysis_2.xlsx
# - Analysis_4.xlsx
# - risk_of_bias.xlsx
#
# Output:
# - characteristic_1.csv
# - characteristic_2.csv
# - characteristic_4.csv
# - forest1.png
# - forest3.png
# - funnel1.png
# - funnel3.png
# - Table of subgroup analysis 
############################################################

#### Initial setting ####
setwd("E:/me/Master Project/Single Collect/Papers")
# if (!require("devtools")) {
#   install.packages("devtools")
# }
# devtools::install_github("MathiasHarrer/dmetar")
library(dmetar) # only available on github, for subgroup analysis and rob summary
library(metafor) # meta analysis model
library(meta) # meta analysis (contain metagen)
library(readxl) # read xlsx file
library(tidyverse) # data manipulating
library(ggplot2) # for analysis 4

##### read #### 
dat1 <- read_excel("Analysis_1.xlsx", sheet=1) # only mean
dat2 <- read_excel("Analysis_2.xlsx", sheet=1) # only median
dat4 <- read_excel("Analysis_4.xlsx", sheet=1) # analysis based on meta-analysis
rob.10 <- read_excel("risk_of_bias.xlsx", sheet=2) # risk of bias assessment of 10 pts


#### 1. First Sets of packages #####

#### 1.1 Preprocessing ####


####  Analysis 1 ------------------------------------------------------------------------------
## Transform CI to SD (seperate CI -> drop CI if SD exists -> calculate SD in two ways and combine)
sum(is.na(dat1$`Sample Size`)) # no sample size value is missing
dat1.sep <- dat1 %>%
  separate(CI, into <- c('LB','UB'), sep = ',') %>%
  mutate(LB = as.numeric(gsub('[(]', '', LB)),
         UB = as.numeric(gsub('[)]', '', UB))) # seperate CI 

# delete CI if SD exists
dat1.1 <- dat1.sep %>% 
  mutate(UB = UB*is.na(SD),
         LB = LB*is.na(SD)) %>% 
  na_if(0) 

#### Sensitivity 2 ####
## exam if CI is symmetric around Mean
dat1.sym <- dat1.1 %>%
  mutate(diff1 = UB - Mean,
         diff2 = Mean-LB,
         sym.ind = (abs(diff1-diff2)>0.03) ) 
# mostly gamma and log normal is not symmetric
# only one Weibull

## gather only symmetric CI
sens2 <- dat1.sym %>% filter(sym.ind == F) %>% 
  mutate(sd = (UB-LB)*sqrt(`Sample Size`)/(2*qt(df=`Sample Size`, 0.975)), # t assumption
         SD = replace_na(SD,0),
         sd = replace_na(sd,0),
         SD = (SD + sd)/sqrt(`Sample Size`)) %>% # do not round any decimal
  select(-c(LB,UB,sd))
  
dim(sens2)

# Transform CI to SD, keep only mean and sd
dat1.cleaned <- dat1.1 %>%
  mutate(sd = (UB-LB)*sqrt(`Sample Size`)/(2*qt(df=`Sample Size`, 0.975)), # t assumption
         SD = replace_na(SD,0),
         sd = replace_na(sd,0),
         SD = (SD + sd)/sqrt(`Sample Size`)) %>% # do not round any decimal
  select(-c(LB,UB,sd))

# to deal with the asymmetric (naive one, choose the shorter part) Do not run!
# dat1.cleaned <- dat1.1 %>%
#   mutate(sd = pmin(UB-Mean,Mean-LB)/qt(df=`Sample Size`, 0.975),
#          sd = replace_na(sd,0),
#          SD = replace_na(SD,0)/sqrt(`Sample Size`),
#          SD = SD + sd) %>% # do not round any decimal
#   select(-c(LB,UB,sd))

# to deal with the asymmetric (simply delete them)
`%notin%` = function(x,y) { !(x %in% y) }
  
ana1.cleaned <- dat1.cleaned %>%
  filter(No. %notin% c(22,24,70,96)) # delete the asymmetric ci

ana1.pre <- ana1.cleaned %>%
  mutate(SD = round(SD,2)) %>% # save 2 decimal places only for presenting
  select(c(3:9,14))
# not that the round() function is not simply 4 down and 5 up, the IEC 60559 standard: ¡®go to the even digit¡¯ is used.


####  Analysis 2--------------------------------------------------------------------------------
# keep meadian and CI -> transform iqr and range -> combine
# use meand and CI -> Mean and SD


# Split CI, IQR and Range
dat2.split <- dat2 %>%
  separate(IQR, into <- c('lb_iqr','ub_iqr'), sep = ',') %>%
  mutate(lb_iqr = as.numeric(gsub('[(]', '', lb_iqr)),
         ub_iqr = as.numeric(gsub('[)]', '', ub_iqr))) %>%
  separate(Range, into <- c('lb_r','ub_r'), sep = ',') %>%
  mutate(lb_r = as.numeric(gsub('[(]', '', lb_r)),
         ub_r = as.numeric(gsub('[)]', '', ub_r))) %>%
  separate(CI2, into <- c('lb_c','ub_c'), sep = ',') %>%
  mutate(lb_c = as.numeric(gsub('[(]', '', lb_c)),
         ub_c = as.numeric(gsub('[)]', '', ub_c)))

# drop IQR and Range if CI exists
# then drop range if IQR exists
dat2.1 <- dat2.split %>% 
  mutate(lb_iqr = lb_iqr*is.na(lb_c),
         ub_iqr = ub_iqr*is.na(lb_c),
         lb_r = lb_r*is.na(lb_c),
         ub_r = ub_r*is.na(lb_c)
         ) %>% 
  mutate(lb_r = lb_r*is.na(lb_iqr),
         ub_r = ub_r*is.na(lb_iqr)) %>%
  na_if(0)
# Use IQR, range, CI -> Mean and sd
dat2.2 <- dat2.1 %>%
  mutate(`Sample Size` = as.numeric(`Sample Size`)) %>% # change sample size to numeric
  filter(!is.na(`Sample Size`)) %>% # drop 2 study that do not have the exact sample size
  mutate(mean_iqr = (Median + ub_iqr + lb_iqr)/3, # 1) iqr transform 
         sd_iqr = (ub_iqr - lb_iqr)/
           (2*qnorm((0.75*`Sample Size` - 0.125)/(`Sample Size` + 0.25))),
         sd_iqr = sd_iqr/sqrt(`Sample Size`),
         mean_r = (lb_r + 2*Median + ub_r)/4, # 2) range transform
         sd_r = (ub_r-lb_r)/
           (2*qnorm((`Sample Size` - 0.375)/(`Sample Size` + 0.25))), 
         sd_r = sd_r/sqrt(`Sample Size`),
         mean_c = Median * (!is.na(ub_c)), # 3) ci transform 
         sd_c = (ub_c-lb_c)/(2*qt(df=`Sample Size`, 0.975))) %>%
  replace(is.na(.),0) %>% # replace all NA with 0
  mutate(Mean = mean_iqr + mean_r + mean_c,
         SD = sd_iqr + sd_r +sd_c)

dat2.cleaned <- dat2.2 %>%
  select(c(1:8,10,11,21))

dat2.pre <- dat2.cleaned %>%
  mutate(Mean = round(Mean,2),
         SD = round(SD,2)) %>%
  select(c(3:11))



#### Analysis 3 ####
dat3.cleaned <- rbind(ana1.cleaned[,1:9],dat2.cleaned[1:9])

#### Analysis 4 ####
# The data from Banka was first transformed using SE=SD/sqrt(n) and then the CIs
# directly use dat4, normal assumption for formula, thus no sample size is needed

dat4.cleaned <- dat4 %>% separate(CI, into <- c('LB','UB'), sep = ',') %>%
  mutate(LB = as.numeric(gsub('[(]', '', LB)),
         UB = as.numeric(gsub('[)]', '', UB)),
         SD = (UB-LB)/(2*qnorm(0.975)))

ana4.pre <- dat4.cleaned %>%
  select(Author, Period, Region, Methodology, `Sample Size`, Mean, SD)

# the meta-analysis is changed to a summary graph: estimates (with CI) vs sample size
# extract number of studies
dat4.plot <- dat4.cleaned %>%
  mutate(`number of studies` = as.numeric(substr(`Sample Size`,1,2)))

# create a new data for jitter
dat4.jitter <- dat4.plot %>%
  mutate(duplicate_id = as.numeric(duplicated(dat4.plot$`number of studies`))) %>%
  mutate(`number of studies` = `number of studies` - 0.75* duplicate_id) %>%
  mutate(duplicate_id = as.factor(duplicate_id))

# create a new data for clusters with the same sample size
dat4.clustered.same <- dat4.plot %>%
  mutate(duplicate_id = as.numeric(duplicated(dat4.plot$`number of studies`))) %>%
  mutate(`number of studies` = `number of studies` - 0.75* duplicate_id) %>%
  mutate(duplicate_id = duplicate_id + duplicated(dat4.plot$`number of studies`,
                                                  fromLast = TRUE)) %>%
  mutate(duplicate_id = as.factor(duplicate_id))

# create a new data for clusters with similar sample size (threshold = 1)
dat4.clustered.same <- dat4.plot %>%
  mutate(duplicate_id = as.numeric(duplicated(dat4.plot$`number of studies`))) %>%
  mutate(`number of studies` = `number of studies` - 0.75* duplicate_id) %>%
  mutate(duplicate_id = duplicate_id + duplicated(dat4.plot$`number of studies`,
                                                  fromLast = TRUE)) %>%
  mutate(duplicate_id = duplicate_id + c(0,1,rep(0,4),1,rep(0,9))) %>% 
  mutate(duplicate_id = as.factor(duplicate_id))

# ggplot2 plot with jittered confidence intervals
ggplot(dat4.clustered.same, aes(x = `number of studies`, y = `Mean`,color = duplicate_id)) +
  geom_point() +
  geom_errorbar(aes(ymin = LB, ymax = UB)) +
  ylab('estimates of mean incubation time with 95% CI (days)') +
  theme_bw() +
  theme(legend.position = 'none') +
  scale_colour_manual(values = c('DodgerBlue','Darkorange')) 
  
#### Sensitivity 1 ####
# directly use dat1.cleaned 
# then combine dat1.cleaned with dat2.cleaned
sens1.1 <- dat1.cleaned # contain the asymetric CIs
sens1.2 <- rbind(dat1.cleaned[,1:9],dat2.cleaned[1:9])
head(dat1.cleaned)
head(dat2.cleaned)
dim(sens1.1)
dim(sens1.2)

#### Summary Statistics ####
# analysis 1

quantile(ana1.cleaned$Mean)

# function return the count and proportion
count.pro = function(var){
  out <- cbind(count = table(var),
              pro = table(var)/length(var))
  return(out)
    
}
count.pro(ana1.cleaned$Region_C)
count.pro(ana1.cleaned$Region_H)
count.pro(ana1.cleaned$Method2)

# analysis 2
sum.ana2 <-  dat2 %>%  
  mutate(`Sample Size` = as.numeric(`Sample Size`)) %>% # change sample size to numeric
  filter(!is.na(`Sample Size`))
dim(sum.ana2)

count.pro(sum.ana2$Region_C)
count.pro(sum.ana2$Region_H)
count.pro(sum.ana2$Method2)



# save the characteristic table

write_excel_csv(ana1.pre, 'characteristic_1.csv')
write_excel_csv(dat2.pre, 'characteristic_2.csv')
write_excel_csv(ana4.pre, 'characteristic_4.csv')
# write_excel_csv(sens1.1, 'sens1.1.csv')

#### 1.2 Summary Effects ####
## Method 1: use metagen(), tutorial on https://bookdown.org/MathiasHarrer/Doing_Meta_Analysis_in_R/pooling-es.html
set.seed(11)
sum1 <- metagen(TE = Mean,
                 seTE = SD,
                 studlab = Author, # study labels
                 data = ana1.cleaned,
                 sm = "Mean", # summary measure
                 comb.fixed = FALSE,
                 comb.random = TRUE, # random model
                 method.tau = "REML", # method of tau: restricted maximum liklihood
                 hakn = TRUE,
                 title = "Analysis 1")
# method.tau: estimators for tau^2: ML, REML, EB(empirical Bayes)
# https://www.metafor-project.org/doku.php/tips:convergence_problems_rma

sum1




## Analysis 2
sum2 <- metagen(TE = Mean,
                seTE = SD,
                studlab = Author, # study labels
                data = dat2.cleaned,
                sm = "Mean", # summary measure
                comb.fixed = FALSE,
                comb.random = TRUE, # random model
                method.tau = "REML",
                hakn = TRUE,
                title = "Analysis 2")

sum2

## Analysis 3
sum3 <- metagen(TE = Mean,
                seTE = SD,
                studlab = Author, # study labels
                data = dat3.cleaned,
                sm = "Mean", # summary measure
                comb.fixed = FALSE,
                comb.random = TRUE, # random model
                method.tau = "REML",
                hakn = TRUE,
                title = "Analysis 2")

sum3

## Analysis 4
sum4 <- metagen(TE = Mean,
                seTE = SD,
                studlab = Author, # study labels
                data = dat4.cleaned,
                sm = "Mean", # summary measure
                comb.fixed = FALSE,
                comb.random = TRUE, # random model
                method.tau = "REML",
                hakn = TRUE,
                title = "Analysis 2")

sum4

## sensitivity analysis 1

sum.sens1.1 <- metagen(TE = Mean,
                seTE = SD,
                studlab = Author, # study labels
                data = sens1.1,
                sm = "Mean", # summary measure
                comb.fixed = FALSE,
                comb.random = TRUE, # random model
                method.tau = "REML",
                hakn = TRUE,
                title = "Analysis 2")

sum.sens1.1


## sensitivity analysis 1.2
sum.sens1.2 <- metagen(TE = Mean,
                       seTE = SD,
                       studlab = Author, # study labels
                       data = sens1.2,
                       sm = "Mean", # summary measure
                       comb.fixed = FALSE,
                       comb.random = TRUE, # random model
                       method.tau = "REML",
                       hakn = TRUE,
                       title = "Analysis 2")

sum.sens1.2

## sensitivity analysis 2
sum.sens2 <- metagen(TE = Mean,
                       seTE = SD,
                       studlab = Author, # study labels
                       data = sens2,
                       sm = "Mean", # summary measure
                       comb.fixed = FALSE,
                       comb.random = TRUE, # random model
                       method.tau = "REML",
                       hakn = TRUE,
                       title = "Analysis 2")

sum.sens2
        
#### 1.3 Forest Plot ####

# size: analysis 1: 800*1250, analysis 2: 800*750, analysis 3: 800*1950

### Analysis 1

forest.meta(sum1,
            leftcols = c('SD'))

forest.meta(sum1, 
            xlim = c(-5,20),
            sortvar = Mean,
            rightcols = c("ci","w.random"),
            leftcols = c('Author','Mean','SD'),
            col.diamond = "lightblue",
            col.diamond.lines = "blue",
            col.predict = "black",
            prediction = TRUE,
            print.I2 = T,
            print.tau2 = FALSE)

forest.meta(sum1, 
            xlim = c(-2,18),
            sortvar = Mean,
            rightlabs = c("95% CI","weight"),
            leftlabs = c('Author','SD'),
           # col.diamond = "lightblue",
           # col.diamond.lines = "blue",
           #  col.predict = "black",
            prediction = TRUE,
            print.I2 = T,
            print.I2.ci = TRUE,
            print.tau2 = FALSE) 

# used
forest.meta(sum1, 
            xlim = c(-2,18),
            sortvar = Mean,
            layout = "JAMA",
            leftlabs = c('Author','Mean','CI'),
            text.random = 'synthetic estimate', # change the pooled result description
            text.predict = '95% PI',
            col.predict = "black",
            colgap.forest.left = unit(20,"mm"))
            

# Analysis 2
forest.meta(sum2, 
            xlim = c(-2,18),
            sortvar = Mean,
            layout = "JAMA",
            text.predict = "95% PI",
            col.predict = "black",
            colgap.forest.left = unit(20,"mm"))

forest.meta(sum2,sortvar = Mean,
            predict = TRUE, 
            print.tau2 = FALSE,
            leftlabs = c("Author", "g", "SE"))

# Analysis 3
forest.meta(sum3, 
            xlim = c(-2,18),
            sortvar = Mean,
            layout = "JAMA",
            text.random = 'synthetic estimate', # change the pooled result description
            text.predict = "95% PI",
            col.predict = "black",
            colgap.forest.left = unit(20,"mm"))

forest.meta(sum3, 
            sortvar = Mean,
            predict = TRUE, 
            print.tau2 = FALSE,
            leftlabs = c("Author", "g", "SE"))

#### 1.4 Subgroup Analysis ####

## type 1: test for region
# test for all regions (hardly useful)
subgroup.analysis.mixed.effects(x = sum1,
                                subgroups = ana1.cleaned$Region)

# test for China and Outside China
an1.sub.c <- subgroup.analysis.mixed.effects(x = sum1,
                                subgroups = ana1.cleaned$Region_C)

s1 <- an1.sub.c$within.subgroup.results[,c(1,2,4,5,8)]



# test for Wuhan and outside Hubei
an1.sub.h <- subgroup.analysis.mixed.effects(x = sum1,
                                subgroups = ana1.cleaned$Region_H)
s2 <- an1.sub.h$within.subgroup.results[,c(1,2,4,5,8)]

# random model between groups
update.meta(sum1,byvar = Region_C,
            tau.common = F)

update.meta(sum1,byvar = Region_C,
            tau.common = T)

update.meta(sum1,byvar = Region_C,
            comb.random = F,
            comb.fixed = T)
            
# test for seven continents: Asia, Europe, North America, South America, Africa and Antarctica



## type 2: test for study design category

# test for all methodology
subgroup.analysis.mixed.effects(x = sum1,
                                subgroups = ana1.cleaned$Methodology)
# test for: parametric, non-para, descriptive
an1.sub.m2 <- subgroup.analysis.mixed.effects(x = sum1,
                                subgroups = ana1.cleaned$Method2)
s3 <- an1.sub.m2$within.subgroup.results[,c(1,2,4,5,8)]


update.meta(sum1,byvar = Method2,
            comb.random = TRUE,
            comb.fixed = FALSE)

## type 3: test for different level of risk
ana.sub.risk <- subgroup.analysis.mixed.effects(x = sum1,
                                                subgroups = ana1.cleaned$ROB)

s4 <- ana.sub.risk$within.subgroup.results[,c(1,2,4,5,8)]
# compute the I^2 CI for specific studies
ana.sub.risk # no CI is shown since n is too small
sums4.lower <- metagen(TE = Mean,
                seTE = SD,
                studlab = Author, # study labels
                data = ana1.cleaned[ana1.cleaned$ROB=='Low',],
                sm = "Mean", # summary measure
                comb.fixed = TRUE,
                comb.random = FALSE, # random model
                method.tau = "REML",
                hakn = TRUE,
                title = "Analysis subgroup_low")

sums4.lower


# gather all 
subgroup1 <- rbind(s1,s2,s3,s4)
subgroup1 <- cbind(rownames(subgroup1),subgroup1)
subgroup1[2:6] <- apply(subgroup1[2:6],2,round,2)
write_excel_csv(subgroup1, 'subgroup_1.csv') # only for demonstration


#### metareg ####
# test if multiple factors affect the results simultaneously
library(metafor) # for multiple metareg

# none of the regression models have a small value for I^2 (residual heterogeneity / unaccounted variability)
model1 <- rma(yi = Mean,
              sei = SD,
              data = ana1.cleaned,
              method = "ML",
              mods = ~ Method2 + ROB + Region_C + Method2*ROB + Method2*Region_C + `Sample Size`,
              test = "knha")
model1

rma(yi = Mean,
    sei = SD,
    data = ana1.cleaned,
    method = "ML",
    mods = ~ Methodology + ROB + Region_C + `Sample Size`,
    test = "knha")

rma(yi = Mean,
    sei = SD,
    data = ana1.cleaned,
    method = "ML",
    mods = ~ Method2,
    test = "knha")


# how about using original ROB score?
rename.dat <- ana1.cleaned[,-3] %>%
  rename('Author' = `Name of Article`)
joint_data <- merge(rename.dat,rob.10,by = 'Author')

rma(yi = Mean,
    sei = SD,
    data = joint_data,
    method = "ML",
    mods = ~ Method2 + All + Region_C + All*Region_C,
    test = "knha") # still high

rma(yi = Mean,
    sei = SD,
    data = joint_data,
    method = "ML",
    mods = ~ All,
    test = "knha")



# auto model selection
multimodel.inference(TE = 'Mean',
                     seTE = 'SD',
                     data = joint_data,
                     predictors = c('All', 'Region_C', 'Method2'),
                     interaction = TRUE)

multimodel.inference(TE = 'Mean',
                     seTE = 'SD',
                     data = ana1.cleaned,
                     predictors = c('ROB', 'Region_H', 'Method2'),
                     interaction = TRUE)

metareg(sum1, ROB+Region_H+Method2)


#### 1.5 Publication Bias ####

# figure size: 600*450

funnel(sum1,xlab = "Hedges'g")
# funnel(sum1,xlab = "g",studlab = TRUE)
funnel(sum3,xlab = "Hedges' g")
# Egger's test
eggers.test(x = sum1)
eggers.test(x = sum2)
eggers.test(x = sum3)
# No substantial asymmetry in the Funnel plot, no publication bias.

# trim the results
trimfill(m.gen)
trimmed.1 <- trimfill(m.gen)
funnel(trimmed.1,xlab = "Hedges' g")
funnel(trimmed.1,xlab = "g",studlab = TRUE)
m.gen$TE.random;trimmed.1$TE.random # compare summary effects

# pcurve
pcurve(sum1)
pcurve(sum3)




#### 1.6 Risk of Bias Summary ####


str(rob.10)
table(rob.10$Risk_Level)
table(rob.10$Risk_Level)/dim(rob.10)[1]

rob.10.n <- rob.10 %>%
  as.data.frame() %>%
  select(-All,Risk_Level)

rob.10.n <- rob.10 %>%
  select(-All,Risk_Level) %>%
  mutate(across(everything(), as.character)) %>%
  as.data.frame()

str(rob.c)


rob.summary(rob.10.n, studies = rob.10.n$Author, table = F,
            name.low = '1',
            name.high = '0')

rob.summary(rob.10.n, studies = rob.10.n$Author, table = T,
            name.high = 0, name.low = 1)


# rob levels was used for subgroup analysis 
# analysis 1

# add a level label


# save the output
write.csv(rob.10,'ROB.csv')

#### self modified plot function #### (some problem with legend label for only high and low)

rob.self <- function (data, name.high = "High", name.unclear = "Unclear", 
                      name.low = "Low", studies, name.missing, table = FALSE) 
{
  if (class(data) != "data.frame") {
    stop("'data' must be of class 'data.frame'.")
  }
  if (missing(name.missing)) {
    colnames.rob = character()
    for (i in 1:ncol(data)) {
      vect = as.character(data[, i])
      for (j in 1:length(data[, i])) {
        if (vect[j] %in% c(name.high, name.unclear, name.low)) {
          colnames.rob[i] = TRUE
        }
        else {
          colnames.rob[i] = FALSE
          message(cat("Column '", colnames(data)[i], 
                      "' removed from plot because it did not contain the specified RoB ratings (only). \n", 
                      sep = ""))
          break
        }
      }
    }
    rob = data[, as.logical(colnames.rob)]
    for (i in 1:ncol(rob)) {
      rob[, i] = as.character(rob[, i])
      rob[rob[, i] == name.high, i] = "High"
      rob[rob[, i] == name.low, i] = "Low"
    }
    if (table == TRUE) {
      if (missing(studies)) {
        stop("'studies' has to be specified when 'table = TRUE'.")
      }
      if (length(as.vector(studies)) != nrow(data)) {
        stop("'studies' vector is not of equal length as the data.")
      }
      if (length(unique(studies)) != length(studies)) {
        stop("'studies' cannot contain duplicate study labels.")
      }
      robby = rob
      robby = data.frame(study = studies, condition = rep(colnames(robby), 
                                                          each = length(studies)), measurement = unlist(robby))
      rownames(robby) = NULL
      robby$condition = gsub("_", " ", robby$condition)
      robby$condition = gsub("\\.", " ", robby$condition)
      robby[robby$measurement == "Low", "measurement"] = "+"
      robby[robby$measurement == "Unclear", "measurement"] = "?"
      robby[robby$measurement == "High", "measurement"] = "-"
      robby$study = factor(robby$study, levels = unique(studies)[rev(order(unique(robby$study)))])
      rob.table = ggplot(data = robby, aes(y = study, x = condition)) + 
        geom_tile(color = "black", fill = "white", 
                  size = 0.8) + geom_point(aes(color = as.factor(measurement)), 
                                           size = 20) + geom_text(aes(label = measurement), 
                                                                  size = 8) + scale_x_discrete(position = "top") + 
        scale_color_manual(values = c(`?` = "#E2DF07", 
                                      `-` = "#BF0000", `+` = "#02C100")) + 
        theme_minimal() + coord_equal() + theme(axis.title.x = element_blank(), 
                                                axis.title.y = element_blank(), axis.ticks.y = element_blank(), 
                                                axis.text.y = element_text(size = 15, color = "black"), 
                                                axis.text.x = element_text(size = 13, color = "black", 
                                                                           angle = 90, hjust = 0), legend.position = "none", 
                                                panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
                                                panel.background = element_blank())
    }
    rob.long = data.frame(condition = rep(colnames(rob), 
                                          each = nrow(rob)), measurement = unlist(rob))
    rownames(rob.long) = NULL
    rob.long$condition = gsub("_", " ", rob.long$condition)
    rob.long$condition = gsub("-", " ", rob.long$condition)
    rob.long$condition = gsub("\\.", " ", rob.long$condition)
    rob.long$measurement = as.factor(rob.long$measurement)
    rob.long$measurement = factor(rob.long$measurement, levels(rob.long$measurement)[c(1, 
                                                                                       3, 2)])
    rob.plot = ggplot(data = rob.long) + geom_bar(mapping = aes(x = condition, 
                                                                fill = measurement), width = 0.7, position = "fill", 
                                                  color = "black") + coord_flip(ylim = c(0, 1)) + 
      guides(fill = guide_legend(reverse = TRUE)) + scale_fill_manual("Risk of Bias", 
                                                                      labels = c("    High risk of bias          ", 
                                                                                  "    Low risk of bias  "), 
                                                                      values = c(High = "#BF0000", 
                                                                                 Low = "#02C100")) + scale_y_continuous(labels = scales::percent) + 
      theme(axis.title.x = element_blank(), axis.title.y = element_blank(), 
            axis.ticks.y = element_blank(), axis.text.y = element_text(size = 18, 
                                                                       color = "black"), axis.line.x = element_line(colour = "black", 
                                                                                                                    size = 0.5, linetype = "solid"), legend.position = "bottom", 
            panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
            panel.background = element_blank(), legend.background = element_rect(linetype = "solid", 
                                                                                 colour = "black"), legend.title = element_blank(), 
            legend.key.size = unit(0.75, "cm"), legend.text = element_text(size = 14))
    plot(rob.plot)
    if (table == TRUE) {
      plot(rob.table)
    }
  }
  else {
    data = as.data.frame(data)
    colnames.rob = character()
    for (i in 1:ncol(data)) {
      vect = as.character(data[, i])
      for (j in 1:length(data[, i])) {
        if (vect[j] %in% c(name.high, name.unclear, name.low, 
                           name.missing)) {
          colnames.rob[i] = TRUE
        }
        else {
          colnames.rob[i] = FALSE
          message(cat("Column '", colnames(data)[i], 
                      "' removed from plot because it did not contain the specified RoB ratings (only). \n", 
                      sep = ""))
          break
        }
      }
    }
    rob = data[, as.logical(colnames.rob)]
    for (i in 1:ncol(rob)) {
      rob[, i] = as.character(rob[, i])
      rob[rob[, i] == name.high, i] = "High"
      rob[rob[, i] == name.low, i] = "Low"
      rob[rob[, i] == name.missing, i] = "Missing"
    }
    if (table == TRUE) {
      if (missing(studies)) {
        stop("'studies' has to be specified when 'table = TRUE'.")
      }
      if (length(as.vector(studies)) != nrow(data)) {
        stop("'studies' vector is not of equal length as the data.")
      }
      robby = rob
      robby = data.frame(study = as.factor(studies), condition = rep(colnames(robby), 
                                                                     each = length(studies)), measurement = unlist(robby))
      rownames(robby) = NULL
      robby$condition = gsub("_", " ", robby$condition)
      robby$condition = gsub("-", " ", robby$condition)
      robby$condition = gsub("\\.", " ", robby$condition)
      robby[robby$measurement == "Low", "measurement"] = "+"
      robby[robby$measurement == "High", "measurement"] = "-"
      robby[robby$measurement == "Missing", "measurement"] = " "
      robby$study = factor(robby$study, levels = unique(studies)[rev(order(unique(robby$study)))])
      rob.table = ggplot(data = robby, aes(y = study, x = condition)) + 
        geom_tile(color = "black", fill = "white", 
                  size = 0.8) + geom_point(aes(color = as.factor(measurement)), 
                                           size = 20) + geom_text(aes(label = measurement), 
                                                                  size = 8) + scale_x_discrete(position = "top") + 
        scale_color_manual(values = c( 
                                      `-` = "#BF0000", `+` = "#02C100", 
                                      ` ` = "white")) + theme_minimal() + 
        coord_equal() + theme(axis.title.x = element_blank(), 
                              axis.title.y = element_blank(), axis.ticks.y = element_blank(), 
                              axis.text.y = element_text(size = 15, color = "black"), 
                              axis.text.x = element_text(size = 13, color = "black", 
                                                         angle = 90, hjust = 0), legend.position = "none", 
                              panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
                              panel.background = element_blank())
    }
    rob.long = data.frame(condition = rep(colnames(rob), 
                                          each = nrow(rob)), measurement = unlist(rob))
    rownames(rob.long) = NULL
    rob.long$condition = gsub("_", " ", rob.long$condition)
    rob.long$condition = gsub("-", " ", rob.long$condition)
    rob.long$condition = gsub("\\.", " ", rob.long$condition)
    rob.long$measurement = as.factor(rob.long$measurement)
    rob.long$measurement = factor(rob.long$measurement, levels(rob.long$measurement)[c(3, 
                                                                                       1, 4, 2)])
    rob.plot = ggplot(data = rob.long) + geom_bar(mapping = aes(x = condition, 
                                                                fill = measurement), width = 0.7, position = "fill", 
                                                  color = "black") + coord_flip(ylim = c(0, 1)) + 
      guides(fill = guide_legend(reverse = TRUE)) + scale_fill_manual("Risk of Bias", 
                                                                      labels = c("  Missing information  ", "  High risk of bias   ", 
                                                                                 "  Unclear risk of bias  ", "  Low risk of bias  "), 
                                                                      values = c(Unclear = "#E2DF07", High = "#BF0000", 
                                                                                 Low = "#02C100", Missing = "white")) + 
      scale_y_continuous(labels = scales::percent) + theme(axis.title.x = element_blank(), 
                                                           axis.title.y = element_blank(), axis.ticks.y = element_blank(), 
                                                           axis.text.y = element_text(size = 18, color = "black"), 
                                                           axis.line.x = element_line(colour = "black", 
                                                                                      size = 0.5, linetype = "solid"), legend.position = "bottom", 
                                                           panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
                                                           panel.background = element_blank(), legend.background = element_rect(linetype = "solid", 
                                                                                                                                colour = "black"), legend.title = element_blank(), 
                                                           legend.key.size = unit(0.75, "cm"), legend.text = element_text(size = 14))
    plot(rob.plot)
    if (table == TRUE) {
      plot(rob.table)
    }
  }
}

rob.self(rob.10.n, studies = rob.10.n$Author, table = F,
         name.high = 0, name.low = 1)



table(rob.10$All)

#### 2. Another sets of Packages ####
# https://wviechtb.github.io/meta_analysis_books/borenstein2009.html#20)_Meta-Regression

# use escalc and rma instead, the aim is to remain the original ci.
# check help file for escalc to include certain data structure before meta-analysis
e1 <- escalc(measure = 'MN', yi = Mean, sei = SD, data = dat1.cleaned)
res1 <- rma(yi = Mean, sei = SD, data = e1) # same results as before

forest(sum2,
       leftcols=c("studlab", "effect", "ci"),
       leftlabs=c("Study", "Incubation time", "95% C.I."),
       rightcols=FALSE)
forest(res1)

# using both SD and CI automatically
m1 <- metagen(TE = Mean,
              seTE = SD,
              lower = LB,
              upper = UB,
              studlab = Author, # study labels
              data = dat1.1,
              method.ci = 't',
              sm = "MN", # summary measure: quantative value for single group
              comb.fixed = FALSE,
              comb.random = TRUE, # random model
              method.tau = "REML",
              hakn = TRUE,
              title = "Analysis 1(using CI)")


forest.meta(m1)


#### Above all: Flow diagram of Study Selection ####
# Original link: http://prisma-statement.org/prismastatement/flowdiagram.aspx
# APP and package: https://www.eshackathon.org/software/PRISMA2020.html
# package site: https://github.com/nealhaddaway/PRISMA2020
library(DiagrammeR)

devtools::install_github("nealhaddaway/PRISMA2020")
library(PRISMA2020)
data <- read.csv(file.choose(), stringsAsFactors=FALSE);
data <- read_PRISMAdata(data);
attach(data); 
plot <- PRISMA_flowdiagram(data,
                           fontsize = 12,
                           interactive = F,
                           previous = FALSE,
                           other = TRUE);
plot

# or use the online app
# https://estech.shinyapps.io/prisma_flowdiagram/

# style for 2009
# https://guides.lib.unc.edu/prisma





