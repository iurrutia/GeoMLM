#############################################################################################
#### CLEANING UP AND TIDYING DATA ###########################################################
#############################################################################################

cchs <- read.csv("CCHS.csv")
# setwd("~/src/JPG1400/Project Data")

summary(cchs)
summary(cchs$ucn_010)
str(cchs)
str(table(cchs$geodpmf))
table(cchs$incghh)
table(cchs$gendmhi)

# Note that we never want values of 6, 7, 8 or 9 to weight down an average,
# thus these are all replaced my NAs

cchs[  cchs== 6 |cchs == 7 |cchs == 8 |cchs == 9 ] <- NA
cchs[  cchs== 96 |cchs == 97 |cchs == 98 |cchs == 99 ] <- NA

# HCS_# responses 6,7,8,9 correspond to NA (not applicable, don't know, refusal, not stated, respectively)
# I am using the omit function (instead of rm) to handle NA's in this instance: https://stackoverflow.com/questions/41588315/the-difference-of-na-rm-and-na-omit-in-r
# Responses originally correspond to excellent (1), good (2), fair (3), and poor (4). I change them to reverse the direction of the relationship to:
# Final: excellent (3), good (2), fair (1), and poor (0) with abs(score - 4)
cchs$hcs_1 <-abs(cchs$hcs_1 -4)
cchs$hcs_2 <-abs(cchs$hcs_2 -4)
cchs$hcs_3 <-abs(cchs$hcs_3 -4)
cchs$hcs_4 <-abs(cchs$hcs_4 -4)


# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
# # # # # # # Variables relating to immigration/race/culture  # # # # # # # # # # # # # # # # # #

# sdcgcgt - Racial origin: 
# White - 1; Other - 2; other - not stated - 9
cchs$sdcgcgt[  cchs$sdcgcgt== 1 ] <- 0
cchs$sdcgcgt[  cchs$sdcgcgt== 2 ] <- 1
# FINAL: White - 0, Other - 1

# sdcgcb13 - Country of birth: 
# CANADA - 1; OTHER - 2; NOT STATED - 9
cchs$sdcgcb13[  cchs$sdcgcb13== 1 ] <- 0
cchs$sdcgcb13[  cchs$sdcgcb13== 2 ] <- 1
# FINAL: Canada - 0, Other - 1

# sdcfimm - Immigrant
# Yes - 1; No - 2; NOT STATED - 9
cchs$sdcfimm[  cchs$sdcfimm== 2 ] <- 0
# FINAL: Yes - 1, No - 0

# sdcglhm - Languages spoken at home
# Original variables for responses are condensed down to 'neither English nor French' or 'not'
# ENGLISH (WITH OR WITHOUT OTHER)- 1; FRENCH (WITH OR WITHOUT OTHER) - 2; ENGLISH & FRENCH (WITH OR W/O OTHER) - 3; NEITHER ENGLISH NOR FRENCH (OTHER)- 4; NOT STATED - 9
cchs$sdcglhm[  cchs$sdcglhm== 1 |cchs$sdcglhm == 2 |cchs$sdcglhm == 3] <- 0
cchs$sdcglhm[  cchs$sdcglhm== 4 ] <- 1
# FINAL: Neither English nor French (other) - 1; Not 'English and/or French' - 0; and NA

# sdcgres - Length of time in Canada since immigration
# 0 TO 9 YEARS - 1; 10 OR MORE YEARS - 2; NOT APPLICABLE - 6; NOT STATED - 9 
cchs$sdcgres[  cchs$sdcgres== 2 ] <- 0
# FINAL: 0 to 9 years - 1, 10 years or more - 0

# cmh_01k - Consulted mental health professional
# 0 TO 9 YEARS - 1; 10 OR MORE YEARS - 2; NOT APPLICABLE - 6; NOT STATED - 9 
cchs$cmh_01k[  cchs$cmh_01k== 2 ] <- 0
# FINAL: Yes - 1; No - 0



# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
# # # # # # # Variables relating to other general socioeconomic factors # # # # # # # # # # # # #
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 

# incghh - Total household income from all sources
# NO INCOME OR LESS THAN $20,000 - 1; $20,000 TO $39,999 - 2; $40,000 TO $59,999 - 3; $60,000 TO $79,999 - 4; $80,000 OR MORE - 5; NOT STATED - 9
# No changes

# dhh_sex- Sex
cchs$dhh_sex[  cchs$dhh_sex== 2 ] <- 0
# FINAL: Male - 1, Female - 0

# dhhgage - Age
# 1  12 TO 14 YEARS
# 2  15 TO 17 YEARS
# 3  18 TO 19 YEARS
# 4  20 TO 24 YEARS
# 5  25 TO 29 YEARS
# 6  30 TO 34 YEARS
# 7  35 TO 39 YEARS
# 8  40 TO 44 YEARS
# 9  45 TO 49 YEARS
# 10  50 TO 54 YEARS
# 11  55 TO 59 YEARS
# 12  60 TO 64 YEARS
# 13  65 TO 69 YEARS
# 14  70 TO 74 YEARS
# 15  75 TO 79 YEARS
# 16  80 YEARS OR MORE


# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
# # # # # # # Other: Variables relating to community, health, mental health # # # # # # # # # # #

# gen_10 - Sense of belonging to local community
# VERY STRONG - 1; SOMEWHAT STRONG - 2; SOMEWHAT WEAK - 3; VERY WEAK - 4; 
# DON'T KNOW - 7; REFUSAL - 8; NOT STATED - 9 (Condensed to 'don't know/refusal/not stated - NA)
# No changes

# gen_02b - Self-perceived mental health
# Note that the scale is 'opposite' to gen_02b, thus changed so they align.
# EXCELLENT - 1; VERY GOOD - 2; GOOD - 3; FAIR - 4; POOR - 5
# DON'T KNOW - 7; REFUSAL - 8; NOT STATED - 9 (Condensed to 'don't know/refusal/not stated - NA)
cchs$gen_02b[  abs(cchs$gen_02b - 5)] # This flips everything to the correct scale. E.g. abs(2-5)=3 (turns 2 to 3), so it turns 5 to 0, 4 to 1, 3 to 2, 2 to 3, and 1 to 5
# Final: poor - 0; fair - 1; good - 2; very good - 3; excellent - 4

# gendmhi - Perceived Mental Health
# POOR - 0; FAIR - 1; GOOD - 2; VERY GOOD - 3; EXCELLENT - 4
# NOT APPLICABLE - 6; DON'T KNOW - 7; REFUSAL - 8; NOT STATED - 9
# Final: poor - 0; fair - 1; good - 2; very good - 3; excellent - 4


# ucn_010 Unmet health care needs - self-perceived (During the past 12 months, was there ever a time when you felt that you needed health care but you didn't receive it?)
# YES - 1; NO - 2; DON'T KNOW - 7; REFUSAL - 8
cchs$ucn_010[  cchs$ucn_010== 2 ] <- 0
# FINAL: Yes - 1, No - 0

# ucn_020i Care not received - dr didn't think it was necessary
# YES - 1; NO - 2; NOT APPLICABLE - 6; DON'T KNOW - 7; REFUSAL - 8; NOT STATED - 9

# ucn_020d Care not received - felt would be inadequate
# YES - 1; NO - 2; NOT APPLICABLE - 6; DON'T KNOW - 7; REFUSAL - 8; NOT STATED - 9

# ucn_020h Care not received - decided not to seek care
# YES - 1; NO - 2; NOT APPLICABLE - 6; DON'T KNOW - 7; REFUSAL - 8; NOT STATED - 9


#############################################################################################
#############################################################################################
########### EXPLORING THE DATA ##############################################################
#############################################################################################
#############################################################################################


# install.packages(c("pgirmess", "rgdal", "spdep", "maptools", "classInt", "RColorBrewer"))

require(maptools)
require(RColorBrewer)
require(classInt)
require(spdep)
require(rgdal)
require(pgirmess)

# Loading and Plotting: STATSCAN Boundary files and simplified boundary file (made with QGIS)

HR_boundaries <- readOGR(dsn = "esriHR", layer = "lhrp000b06a_e_Oct2011") #directory and then name of shapefile w/o extension (.shp)
summary(HR_boundaries)
# plot(HR_boundaries)

# A simplified shapefile 
HR_bsim <- readOGR(dsn = "esriHR", layer = "simplified") #directory and then name of shapefile w/o extension (.shp)
# summary(HR_bsim)
str(HR_bsim)
# plot(HR_bsim)

# CONVERSION OF HEALTH REGION CODES 
# install.packages("varhandle")
library(varhandle)
?unfactor

# To update the HR codes from the boundary files to the format of the CCHS region codes, we need to 
# instert a 9 into the middle of each number. E.g. 1011 <- 10911, 1102 <-11902, etc.
# In order to do this, we could insert the 9 in as though the number was a string, or ... we can use math!

conversion <- function(old_code){
  new_code <- floor(old_code*(0.01))*1000 + 900 + (old_code - (floor(old_code*(0.01)))*100)
  return(new_code)
}

# Here's an example:
conversion(1011)

# We convert teh HR codes in HR_bsim to the CCHS code format
# summary(HR_bsim)
HR_bsim$HR_new <- HR_bsim$PR_HRUID
HR_bsim$HR_conv <- unfactor(HR_bsim$HR_new)
HR_bsim$HR_new <- sapply(HR_bsim$HR_conv, conversion)

# Now HR_bsim$_new contains the HR code values as per the CCHS file.
# There are a few health regions that are merged in the CCHS data:
HR_bsim$HR_new [ HR_bsim$HR_new == 10914 ] <- 10913
HR_bsim$HR_new [ HR_bsim$HR_new == 11902 |HR_bsim$HR_new == 11903] <- 11901
HR_bsim$HR_new [ HR_bsim$HR_new == 12911 |HR_bsim$HR_new == 12912] <- 12910
HR_bsim$HR_new [ HR_bsim$HR_new == 12934 |HR_bsim$HR_new == 12935] <- 12930
HR_bsim$HR_new [ HR_bsim$HR_new == 12946 |HR_bsim$HR_new == 12947] <- 12940
HR_bsim$HR_new [ HR_bsim$HR_new == 13905 ] <- 13904
HR_bsim$HR_new [ HR_bsim$HR_new == 13907 ] <- 13906
HR_bsim$HR_new [ HR_bsim$HR_new == 24910 ] <- 24909
HR_bsim$HR_new [ HR_bsim$HR_new == 35954 ] <- 35939
HR_bsim$HR_new [ HR_bsim$HR_new == 46945 ] <- 46915
HR_bsim$HR_new [ HR_bsim$HR_new == 46925 ] <- 46920
HR_bsim$HR_new [ HR_bsim$HR_new == 46930 |HR_bsim$HR_new == 46960] <- 46931
HR_bsim$HR_new [ HR_bsim$HR_new == 46985 ] <- 46970
HR_bsim$HR_new [ HR_bsim$HR_new == 47903 |HR_bsim$HR_new == 47902] <- 47901
HR_bsim$HR_new [ HR_bsim$HR_new == 47908 ] <- 47905
HR_bsim$HR_new [ HR_bsim$HR_new == 47910 ] <- 47907
HR_bsim$HR_new [ HR_bsim$HR_new == 47914 ] <- 47909
HR_bsim$HR_new [ HR_bsim$HR_new == 61901 |HR_bsim$HR_new == 62901] <- 60901

HR_bsim$HR_new

# Looks good. Now we save HR_bsim$HR_new back into the HR_bsim$PR_HRUID variable

HR_bsim$PR_HRUID <- HR_bsim$HR_new

# Unfortunately it seems data is only available in Alberta for the variables I wanted to use as dependent variables.


aggregate(cchs$hcs_1, list(cchs$geodpmf), mean, na.rm=TRUE)
aggregate(cchs$hcs_2, list(cchs$geodpmf), mean, na.rm=TRUE)
aggregate(cchs$hcs_3, list(cchs$geodpmf), mean, na.rm=TRUE)
aggregate(cchs$hcs_4, list(cchs$geodpmf), mean, na.rm=TRUE)

table(cchs$hcs_1,cchs$geodpmf)
table(cchs$hcs_2,cchs$geodpmf)
table(cchs$hcs_3,cchs$geodpmf)
table(cchs$hcs_4,cchs$geodpmf)

# Instead, I now examine the data set to see how other potential dependent variables are distributed geographically

HRgen_10 <- aggregate(cchs$gen_10, list(cchs$geodpmf), mean, na.rm=TRUE)
HRgen_02b <- aggregate(cchs$gen_02b, list(cchs$geodpmf), mean, na.rm=TRUE)
HRgen_07 <- aggregate(cchs$gen_07, list(cchs$geodpmf), mean, na.rm=TRUE)
HRgendmhi <- aggregate(cchs$gendmhi, list(cchs$geodpmf), mean, na.rm=TRUE)
HRcmh_01k <- aggregate(cchs$cmh_01k, list(cchs$geodpmf), mean, na.rm=TRUE)
HRcmhg01l <- aggregate(cchs$cmhg01l, list(cchs$geodpmf), mean, na.rm=TRUE)
HRucn_010 <- aggregate(cchs$ucn_010, list(cchs$geodpmf), mean, na.rm=TRUE)
HRdhh_sex <- aggregate(cchs$dhh_sex, list(cchs$geodpmf), mean, na.rm=TRUE)
HRdhhgage <- aggregate(cchs$dhhgage, list(cchs$geodpmf), mean, na.rm=TRUE)
HRincghh <- aggregate(cchs$incghh, list(cchs$geodpmf), mean, na.rm=TRUE)
HRsdcgcgt <- aggregate(cchs$sdcgcgt, list(cchs$geodpmf), mean, na.rm=TRUE)
HRsdcglhm <- aggregate(cchs$sdcglhm, list(cchs$geodpmf), mean, na.rm=TRUE)
HRsdcgcb13 <- aggregate(cchs$sdcgcb13, list(cchs$geodpmf), mean, na.rm=TRUE)
HRsdcfimm <- aggregate(cchs$sdcfimm, list(cchs$geodpmf), mean, na.rm=TRUE)
HRsdcgres <- aggregate(cchs$sdcgres, list(cchs$geodpmf), mean, na.rm=TRUE)


# We re-name the columns of the new frame created above so that they can be merged with the boundary files.

names(HRgen_10)[names(HRgen_10)=="Group.1"] <- "HR_new"
names(HRgen_10)[names(HRgen_10)=="x"] <- "gen_10"

names(HRgen_02b)[names(HRgen_02b)=="Group.1"] <- "HR_new"
names(HRgen_02b)[names(HRgen_02b)=="x"] <- "gen_02b"

names(HRgen_07)[names(HRgen_07)=="Group.1"] <- "HR_new"
names(HRgen_07)[names(HRgen_07)=="x"] <- "gen_07"

names(HRgendmhi)[names(HRgendmhi)=="Group.1"] <- "HR_new"
names(HRgendmhi)[names(HRgendmhi)=="x"] <- "gendmhi"

names(HRcmh_01k)[names(HRcmh_01k)=="Group.1"] <- "HR_new"
names(HRcmh_01k)[names(HRcmh_01k)=="x"] <- "cmh_01k"

names(HRcmhg01l)[names(HRcmhg01l)=="Group.1"] <- "HR_new"
names(HRcmhg01l)[names(HRcmhg01l)=="x"] <- "cmhg01l"

names(HRucn_010)[names(HRucn_010)=="Group.1"] <- "HR_new"
names(HRucn_010)[names(HRucn_010)=="x"] <- "ucn_010"

names(HRdhh_sex)[names(HRdhh_sex)=="Group.1"] <- "HR_new"
names(HRdhh_sex)[names(HRdhh_sex)=="x"] <- "dhh_sex"

names(HRdhhgage)[names(HRdhhgage)=="Group.1"] <- "HR_new"
names(HRdhhgage)[names(HRdhhgage)=="x"] <- "dhhgage"

names(HRincghh)[names(HRincghh)=="Group.1"] <- "HR_new"
names(HRincghh)[names(HRincghh)=="x"] <- "incghh"

names(HRsdcgcgt)[names(HRsdcgcgt)=="Group.1"] <- "HR_new"
names(HRsdcgcgt)[names(HRsdcgcgt)=="x"] <- "sdcgcgt"

names(HRsdcglhm)[names(HRsdcglhm)=="Group.1"] <- "HR_new"
names(HRsdcglhm)[names(HRsdcglhm)=="x"] <- "sdcglhm"

names(HRsdcgcb13)[names(HRsdcgcb13)=="Group.1"] <- "HR_new"
names(HRsdcgcb13)[names(HRsdcgcb13)=="x"] <- "sdcgcb13"

names(HRsdcfimm)[names(HRsdcfimm)=="Group.1"] <- "HR_new"
names(HRsdcfimm)[names(HRsdcfimm)=="x"] <- "sdcfimm"

names(HRsdcgres)[names(HRsdcgres)=="Group.1"] <- "HR_new"
names(HRsdcgres)[names(HRsdcgres)=="x"] <- "sdcgres"


# While we're at it, we sohuld check how the survey respondents for our potential dependent variables are
# distributed accross different provinces and health regions too.

table(cchs$gen_10, cchs$geodpmf) # All regions GOOD!
table(cchs$gen_02b, cchs$geodpmf) # All regions GOOD!
table(cchs$gen_07, cchs$geodpmf) # All regions GOOD!
table(cchs$gendmhi, cchs$geodpmf) # All regions GOOD!
table(cchs$cmh_01k, cchs$geodpmf) # Missing responses from 6 province regions MARGINAL
table(cchs$cmh_01k, cchs$geogprv) # Checking provinces for same variable as above
table(cchs$cmhg01l, cchs$geodpmf) # Missing responses from 6 province regions MARGINAL
table(cchs$cmhg01l, cchs$geogprv) # Checking provinces for same variable as above
table(cchs$ucn_010, cchs$geodpmf) # All regions GOOD!
table(cchs$gen_10, cchs$geodpmf) # gen_10 - Already considered above
table(cchs$sdcgcgt, cchs$geodpmf) # sdcgcgt GOOD!
table(cchs$sdcglhm, cchs$geodpmf) # sdcglhm GOOD!
table(cchs$sdcgcb13, cchs$geodpmf) # sdcgcb13 GOOD!
table(cchs$sdcfimm, cchs$geodpmf) # sdcfimm GOOD!
table(cchs$sdcgres, cchs$geodpmf) # sdcgres GOOD!
table(cchs$dhh_sex, cchs$geodpmf) # dhh_sex GOOD!
table(cchs$dhhgage, cchs$geodpmf) # dhhgage GOOD!
table(cchs$incghh, cchs$geodpmf) # incghh GOOD!

# We note that some variables are well represented in the survey in each health region
# These variables are all suitable (in that respect) as independent variables
# For the sake of exploring the variables we have, and their variation by HR, we create some plots
# Left outer join: merge(x = df1, y = df2, by = "common_variable", all.x = TRUE)

test <- merge(x = HR_bsim, y = HRgen_10, by = "HR_new", all.x = TRUE)
test <- merge(x = test, y = HRgen_02b, by = "HR_new", all.x = TRUE)
test <- merge(x = test, y = HRgen_07, by = "HR_new", all.x = TRUE)
test <- merge(x = test, y = HRgendmhi, by = "HR_new", all.x = TRUE)
test <- merge(x = test, y = HRgen_10, by = "HR_new", all.x = TRUE)
test <- merge(x = test, y = HRsdcgcgt, by = "HR_new", all.x = TRUE)
test <- merge(x = test, y = HRsdcglhm, by = "HR_new", all.x = TRUE)
test <- merge(x = test, y = HRsdcgcb13, by = "HR_new", all.x = TRUE)
test <- merge(x = test, y = HRsdcfimm, by = "HR_new", all.x = TRUE)
test <- merge(x = test, y = HRsdcgres, by = "HR_new", all.x = TRUE)
test <- merge(x = test, y = HRdhh_sex, by = "HR_new", all.x = TRUE)
test <- merge(x = test, y = HRdhhgage, by = "HR_new", all.x = TRUE)
test <- merge(x = test, y = HRincghh, by = "HR_new", all.x = TRUE)

nclr<-4
plotvar <- test$sdcgcgt
class <- classIntervals(plotvar, nclr, style = "jenks",dataPrecision = 2)
plotclr <- brewer.pal(nclr, "Blues")
colcode <- findColours(class, plotclr, digits = 3)
plot(columbus, col = colcode, pch = 19, axes = T,cex=1.5)
legend("topleft", legend = names(attr(colcode, "table")),fill = attr(colcode, "palette"), cex = 0.8)





#############################################################################################
#############################################################################################
########### MULTILEVEL MODEL ################################################################
#############################################################################################
#############################################################################################



#Import the library for multilevel logistic modeling
library(lme4)

# Checking correlation:
cor(cchs$sdcfimm, cchs$incghh, use="complete.obs")
cor(cchs$sdcfimm, cchs$gendmhi, use="complete.obs")


###############################
######## FIRST MODEL ##########
###############################

######## STEP 0: 

# Calculating the grand-mean centered variables.
# The new variables are "variable_gmc".
# In our case the primary variable of interest is wether or not if someone is an immigrant (sdcfimm) affects
# whether or not they have reported unmet mental health care needs, which is a level-1 variable
# but it is also reasonable to consider whether their response for gendmhi (perceived mental health) affects 
# their reported unmet mental health care needs ucn_010
grand_mean_gendmhi <- mean(cchs$gendmhi, na.rm=T)
cchs$gendmhi_gmc <- cchs$gendmhi - grand_mean_gendmhi

grand_mean_sdcfimm <- mean(cchs$sdcfimm, na.rm=T)
cchs$sdcfimm_gmc <- cchs$sdcfimm - grand_mean_sdcfimm

######## STEP 1: 

M0 <- glmer(ucn_010 ~ ( 1 | geodpmf), data=cchs, family=binomial)
summary(M0)
icc <- M0@theta[1]^2/ (M0@theta[1]^2 + (3.14159^2/3))
icc
# icc = 0.01698266
# And we note that the random intercept variance is very, very low: 0.05684
######## STEP 2: 

# We focus on the between-observation effect of the (level-1) variable, thus we use the grand-mean centered variables.


# Below are the commands to run the constrained intermediate model (CIM); 

CIM <- glmer(ucn_010 ~ gendmhi_gmc + sdcfimm_gmc + (1 | geodpmf), data = cchs, family = "binomial")
#summary(CIM)
paste("FYI: The deviance of the CIM is:", CIM@devcomp$cmp[[8]])

# Below are the commands to run the augmented intermediate model (AIM); 
# the model similar to the constrained intermediate model with the exception that it includes random slope terms).

# Each lower level variable is tested by itself:
AIM1 <- glmer(ucn_010 ~ gendmhi_gmc + sdcfimm_gmc  + (1 + gendmhi_gmc || geodpmf), data = cchs, family = "binomial")
#summary(AIM1)
paste("FYI: The deviance of the AIM1 is:", AIM@devcomp$cmp[[8]])

AIM2 <- glmer(ucn_010 ~ gendmhi_gmc + sdcfimm_gmc + (1 + sdcfimm_gmc || geodpmf), data = cchs, family = "binomial")
#summary(AIM2)
paste("FYI: The deviance of the AIM2 is:", AIM@devcomp$cmp[[8]])


#Below is the command to determine whether including the random slopes of centered variables improves the model. 
#The software performs a likelihood-ratio test LR X(1)^2,  comparing the deviance of the CIM with the deviance of the AIM.
anova(CIM, AIM1)
anova(CIM, AIM2)

# We note from anova() that AIC2 is better in comparicon to CIM, and statistically significant (we can also see that the deviance for AIM1 and AIM2 are slightly smaller than for CIM.)


#####################################################################
#STEP #3: Building the Final Model
#####################################################################

#Below is the command to run the final model (adding the cross-level interaction(s)). 
FM <- glmer(ucn_010 ~ gendmhi_gmc + sdcfimm_gmc + (1 + sdcfimm_gmc | geodpmf), data = cchs, family = "binomial")
sumFM <- summary(FM)

capture.output(sumFM, file = "sumFM.txt")

#Note:
# we have decided here to keep the random slope component of "sdcfimm_gmc"


#Below are the command to compare the coefficient estimates obtained in the final model,
#with (glmer) or without (glm) the use of multilevel modelling.
GLMER <- glmer(ucn_010 ~ gendmhi_gmc + sdcfimm_gmc + (1 + sdcfimm_gmc | geodpmf), data = cchs, family = "binomial")
GLM <-  glm(ucn_010 ~ gendmhi_gmc + sdcfimm_gmc, data = cchs, family = "binomial")
GLM2 <-  glm(ucn_010 ~ gendmhi_gmc , data = cchs, family = "binomial")


summary(GLMER)
summary(GLM)
summary(GLM2)


#####################################################################
# Interpreting the results
#####################################################################

OR <- exp(fixef(FM))
CI <- exp(confint(FM,parm="beta_")) 
OR
CI

OR_m1SD <- exp(fixef(FM_m1SD))
CI_m1SD <- exp(confint(FM_m1SD,parm="beta_")) 
OR_m1SD
CI_m1SD

OR_p1SD <- exp(fixef(FM_p1SD))
CI_p1SD <- exp(confint(FM_p1SD,parm="beta_")) 
OR_p1SD
CI_p1SD

OR.CI<-rbind(cbind(OR,CI), cbind(OR_m1SD,CI_m1SD)[3,], cbind(OR_p1SD,CI_p1SD)[3,])
rownames(OR.CI)<-c(rownames(cbind(OR,CI)), "var1", "var2")
OR.CI
OR.CI

