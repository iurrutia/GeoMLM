# README


---

[1. Overview](#overview)

[2. Data](#data)

[3. Model](#model)

[4. Conclusions and next steps](#concl)


---

## <a name="overview">Overview </a>

This project examines  differences in respondents'  assessments of healthcare quality and access, to assess if these varied by immigration history and ethnicity, using a multilevel model to examine differences between health regions. Specifically, I use data from the Canadian Community Health Survey (CCHS) to explore whether certain population groups are more/less adequately served by healthcare service options in Canada, using a multilevel model to account for possible differences in healthcare between health regions.

This script was used to create a model and findings which informed a report.



## <a name="data">Data </a>

The Canadian Community Health Survey (CCHS) can be accessed through the University of Toronto here: 
http://sda.chass.utoronto.ca.myaccess.library.utoronto.ca/sdaweb/html/cchs.htm (CHASS Microdata Analysis and Subsetting with SDA, Faculty of Arts & Sciences, University of Toronto)

This data set has 63522 observations, which are distributed among 97 simplified health regions. My final model considered following variables:

      geodpmf  :   Health Region  &  97 health regions  
      incghh   :   Total household income  & Five income buckets 
      sdcgcgt  :   Cultural/racial origin & Binary: white or not white
      sdcfimm  :   Immigrant & Binary: immigrant or non-immigrant
      gendmhi  :   Perceived mental health & Five options (`poor' to `excellent')   
      Dependent variable:
      ucn_010  : Unmet healthcare needs - self-perceived (Binary: y/n)


## <a name="model">Model </a>

The script follows these steps:
- Preliminary phase: Cluster- or grand-mean centering variables. 
- Step 1: Running an empty model and calculating the intraclass correlation coefficient (ICC) 
- Step 2: Running a constrained and an augmented intermediate model and performing a likelihood ratio test to determine whether considering the cluster-based variation of the effect of the lower- level variable improves the model fit
- Step 3: Running a final model and interpreting the odds ratio and confidence intervals to determine whether data support the hypothesis

##### Preliminary phase:
For ease of interpretation, this preliminary step involves centering variables (shifting them so the mean corresponds to zero). The choice of whether to center variables by the cluster-mean of the dependent variable, or by the overall mean of is based on whether we wish to test the between-observation or within-cluster effect, respectively. (Sommet and Morselli call variables *level 1* if they measure observation-level characteristics, *level 2* if they refer to first order cluster related characteristics, and so on.) In other words, the variable is centered on the level it is located. In our case we want to investigate wether or not  someone being an immigrant (sdcfimm) affects their reported unmet mental health needs. (We also consider their perceived mental health.) These are level-1 variables measuring observation-level characteristics, therefore we center the variables by the overall mean.
##### Step 1:
We note here that the intraclass correlation coefficient is extremely low (less than 2\%). (This is a sign that a multilevel model is not the best model to use for this project.) In this step, we also find  that the random intercept variance is very low (0.05684) meaning, again, that there is not much variation in the log(odds) of *ucn_010* based on which health region a respondent resides in.
##### Step 2: 
We create constrained and augmented intermediate models to determine wether or not we should keep the random slope variance in our final model. We have two level 1 variables we want to consider including in the model (*sdcfimm* and *gendmhi*), so we create an augmented model for each of these. We note that deviance of all three models is very high, and very similar. We also note that the both AIC models have better results than CIM, but only the augmented model with the (centered) random slope variance component for *sdcfimm* is significant.

##### Step 3:
Based on our results in Step 2, we add a term to account for the slope of the grand  mean centered *sdcfimm* variable to our final model. In this step, we would include any interaction terms between level 1 and 2 variables, which is not applicable to our data set.


Having centered the variables, we interpret the OR value for *sdcfimm_gmc* to mean that survey respondents who indicated they were immigrants were 1.07 times more likely to report unmet healthcare needs. We note that the value 1 is contained in the 95% confidence interval for *sdcfimm_gmc*, so we conclude that we cannot reject the null statement of our stated hypothesis. In other words, we do not conclude any significant effects of respondents' immigration history on their reported unmet healthcare needs.

## <a name="concl">Conclusions and next steps </a>

The model included terms relating to respondents' reported mental health and whether or not they were immigrants as fixed effects, a random slope term for *sdcfimm*, and included no level-2 variables (and thus no interaction terms with level 2 variables).

The example of the multilevel model outlined above is one of several multilevel models that I attempted to fit to the data. Other attempted models used different combinations of the variables, including breaking the data into clusters by *sdcfimm* and *sdcgcgt*. The only other potential multilevel model with a higher ICC value than the example above was ``glmer(cmh_01k ~ ( 1 | sdcfimm), data=cchs, family = "binomial")``, with an ICC of 2.1%. In other words, out of all the multilevel models I attempted to fit, at most 2.1% of the chances of a survey respondent indicating they had unmet healthcare needs was explained by breaking the data into clusters using the variables I proposed. These ICC values are very low, indicating that survey respondents regarding unmet healthcare needs within the same health region, or into clusters of immigrant and non-immigrant survey respondents, vary widely. Thus, we conclude that for this data set, using these variables, grouping survey responses into clusters contributes very little towards understanding the relationship between unmet healthcare needs and the independent variables chosen. Therefore a multilevel model is not appropriate for this data set. Indeed, running a simple logistic regression (GLM and GLM2 in the R script in this repo) shows that the multilinear model shows no big improvement over a basic logistic regression for this data set.

It is likely the case that the health regions are too large to constitute meaningfully different clusters for a multilevel analysis. A future endeavour could consider a similar project set-up using less homogenous clusters such as census tracks. Such an approach might find more meaningful results by capturing the variation between clusters based on more granular demographic characteristics into a model.