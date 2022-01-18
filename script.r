#installing relevant packages
#install.packages("tidyverse")
#install.packages("meta")
#dir.create("/Testing_6Jan")
# setwd("/Users/icmr7372/Desktop/LAB/Meta/Subgroup/Meta_Excel/CC/")
#setwd("/Users/icmr7372/Desktop/LAB/Meta/Subgroup/Meta_Excel/PC+CC/")
#load these installed packages
library(tidyverse)
library(meta)
library(ggplot2)
#library(xlsx)
# # generating data for all models
# # load excel containg data
library(readxl)
library(dplyr)
# allele model
# Selection For Cancer and Excel file


# setwd("../../../LAB/Endo-coding+noncoding-meta/8JUNE2021_Thesis_Meta/")

eggers.test = function(x) {
  
  # Validate
  x = x
  
  if (x$k < 10) {
    
    warning(paste("Your meta-analysis contains k =",
                  x$k, "studies. Egger's test may lack the statistical power to detect bias when the number of studies is small (i.e., k<10)."))
    
  }
  
  if (class(x)[1] %in% c("meta", "metabin", "metagen", "metacont", "metacor", "metainc", "metaprop")) {
    
    # Conduct metabias
    eggers = meta::metabias(x, k.min = 3, method = "linreg")
    
    # Get Intercept
    intercept = as.numeric(eggers$estimate[1])
    
    # Get SE
    se = as.numeric(eggers$estimate[2])
    
    # Calculate 95CI
    llci = intercept - qnorm(0.975) * se
    ulci = intercept + qnorm(0.975) * se
    
    # Get t
    t = as.numeric(eggers$statistic)
    
    # Get df
    df = as.numeric(eggers$parameters)
    
    # Get p
    p = as.numeric(eggers$p.value)
    
    # Make df
    returnlist = list(intercept = intercept,
                      llci = llci,
                      ulci = ulci,
                      t = t,
                      p = p,
                      meta.obj = x)
    
  } else {
    
    stop("x must be of type 'metabin', 'metagen', 'metacont', 'metainc' or 'metaprop'")
    
  }
  
  class(returnlist) = "eggers.test"
  
  return(returnlist)
  
}


Cancer ='_All_'
all_model = data.frame()
data_dom = data.frame()
data_res = data.frame()
data_het = data.frame()
data_hom = data.frame()

data_eggers_allele_all <- data.frame()
data_eggers_dominant_all <- data.frame()
data_eggers_recessive_all <- data.frame()
data_eggers_Heterozygous_all <- data.frame()
data_eggers_Homozygous_all <- data.frame()

#Add path of input file

RSID <- excel_sheets(path = "Cervical_All_Genotype.xlsx")

dir.create(Cancer) 

for( i in RSID) {
  arsid = gsub(".*_","",i)
  agene = gsub("_.*","",i)
  
  #Add path of input file
  print(i)
	
  names <- read_excel(path = "Cervical_All_Genotype.xlsx", sheet =i)
  length <- nrow(names)
  data_genotype <- names %>%
    transmute(
      PMID = (names[[3]]),
      Author = (names[[4]]),
      Year = (names[[5]]),
      Case_11 = (names[[8]]),
      Case_12 = (names[[9]]),
      Case_22 = (names[[10]]),
      Control_11 = (names[[11]]),
      Control_12 = (names[[12]]),
      Control_22 = (names[[13]])
      
    )
  data_all <- data_genotype %>% 
    mutate(Case_1 = Case_11+Case_11+Case_12) %>%
    mutate(Case_2 = Case_22+Case_22+Case_12) %>%
    mutate(Control_1 = Control_11+Control_11+Control_12) %>%
    mutate(Control_2 = Control_22+Control_22+Control_12)
  
  # calculating allele model
  data_allele <- data_all %>% transmute(
    Author_Name = Author,
    Year = Year,
    Case_Events_A = Case_2,
    Case_Total_A = Case_1+Case_2,
    Control_Events_A = Control_2,
    Control_Total_A = Control_1+Control_2
  )
  
  #calculating dominant model
  data_dominant <- data_all %>% transmute(
    Author_Name = Author,
    Year = Year,
    Case_Events_D = Case_22+Case_12,
    Case_Total_D = Case_11+Case_12+Case_22,
    Control_Events_D = Control_22+Control_12,
    Control_Total_D = Control_11+Control_12+Control_22
  )
  
  # calculating recessive model
  data_recessive <- data_all %>% transmute(
    Author_Name = Author,
    Year = Year,
    Case_Events_R = Case_22,
    Case_Total_R = Case_11+Case_12+Case_22,
    Control_Events_R = Control_22,
    Control_Total_R = Control_11+Control_12+Control_22
  )
  
  # calculating heterozygous model
  data_heterozygous <- data_all %>% transmute(
    Author_Name = Author,
    Year = Year,
    Case_Events_He = Case_12,
    Case_Total_He = Case_12+Case_11,
    Control_Events_He = Control_12,
    Control_Total_He = Control_12+Control_11
  )
  
  # calculating homozygous model
  data_homozygous <- data_all %>% transmute(
    Author_Name = Author,
    Year = Year,
    Case_Events_Ho = Case_22,
    Case_Total_Ho = Case_22+Case_11,
    Control_Events_Ho = Control_22,
    Control_Total_Ho = Control_22+Control_11
  )
  
  #files.list <- list.files(path = "Meta_input", pattern = "*.xlsx",full.names = TRUE)
  #files
  #all_data <- lapply(files,read_excel,sheet = "Allele") 
  #all_data
  
  #Author = paste(Author_Name,'(',Year,')',sep = '')
  
  # creating directory specific to RSID
  #dir.create(i)
  
  #Exporting Data For all 5 Models
  #y= paste(i,i,sep = '/')
  #y=paste(y,'xlsx',sep = '.')
  #write.xlsx(data_allele, file = y, sheetName = "Allele")
  #write.xlsx(data_dominant, file = y, sheetName = "Dominant", append = TRUE)
  #write.xlsx(data_recessive, file = y, sheetName = "Recessive", append = TRUE)
  #write.xlsx(data_heterozygous, file = y, sheetName = "Heterozygous", append = TRUE)
  #write.xlsx(data_homozygous, file = y, sheetName = "Homozygous", append = TRUE)
  
  # calculating OR with 95% CI for allele model
  Allele <- paste(i, "OR_A", sep=Cancer)
  
  
  OR_A <- metabin(Case_Events_A,Case_Total_A, Control_Events_A, Control_Total_A, title = "Allele Model",
                  data = data_allele,
                  studlab = paste(Author_Name,'(',Year,')',sep = ''),
                  comb.fixed = TRUE,
                  comb.random = TRUE,
                  method.tau = "SJ",
                  hakn = TRUE,
                  prediction = TRUE,
                  incr = 0.1,
                  sm = "OR")
  if (is.na(OR_A["I2"])){
    OR_A <- "This Value is NA"
    
    data_A <- data.frame(
      RSID = arsid,
      Gene = agene,
      Model = "NA",
      Study = "NA",
      OR_CI = "NA",
      Weight = "NA",
      Cochran.Q = "NA",
      I2 = "NA",
      T2 = "NA",
      Pooled.OR_CI = "NA",
      pvalue = "NA",
      Eggers = 'NA'
    )
  }else if (OR_A$I2 <= 0.50) {
    OR_A <- metabin(Case_Events_A,Case_Total_A, Control_Events_A, Control_Total_A, title = "Allele Model",
                    data = data_allele,
                    studlab = paste(Author_Name,'(',Year,')',sep = ''),
                    comb.fixed = TRUE,
                    comb.random = FALSE,
                    method.tau = "SJ",
                    hakn = TRUE,
                    prediction = TRUE,
                    incr = 0.1,
                    sm = "OR")
    
    A.et <- eggers.test(x = OR_A)
    
    data_A <- data.frame(
      RSID = arsid,
      Gene = agene,
      Model = "Fixed",
      Study = OR_A$studlab,
      OR_CI = paste(round(exp(OR_A$TE), 
                          digits = 2),'(',round(exp(OR_A$lower), digits = 2),'-',
                    round(exp(OR_A$upper), digits = 2),')'),
      Weight = paste(round(((OR_A$w.fixed/sum(OR_A$w.fixed))*100), digits = 2), '%') ,
      Cochran.Q = round(OR_A$pval.Q, digits = 2),
      I2 = round((OR_A$I2*100), digits = 2),
      T2 = round(OR_A$tau2, digits = 2),
      Pooled.OR_CI = paste(round(exp(OR_A$TE.fixed), digits = 2),'(',round(exp(OR_A$lower.fixed), digits = 2),'-',round(exp(OR_A$upper.fixed), digits = 2),')'),
      pvalue = round(OR_A$pval.fixed, digits = 2),
      Eggers = round(ifelse(length(A.et$p) == 0,0, A.et$p),digits = 2)
    )
    
  } else {
    OR_A <- metabin(Case_Events_A,Case_Total_A, Control_Events_A, Control_Total_A, title = "Allele Model",
                    data = data_allele,
                    studlab = paste(Author_Name,'(',Year,')',sep = ''),
                    comb.fixed = FALSE,
                    comb.random = TRUE,
                    method.tau = "SJ",
                    hakn = TRUE,
                    prediction = TRUE,
                    incr = 0.1,
                    sm = "OR")
    
    A.et <- eggers.test(x = OR_A)
    
    data_A <- data.frame(
      RSID = arsid,
      Gene = agene,
      Model = "Random",
      Study = OR_A$studlab,
      OR_CI = paste(round(exp(OR_A$TE), 
                          digits = 2),'(',round(exp(OR_A$lower), digits = 2),'-',
                    round(exp(OR_A$upper), digits = 2),')'),
      Weight = paste(round(((OR_A$w.random/sum(OR_A$w.random))*100), digits = 2), '%') ,
      Cochran.Q = round(OR_A$pval.Q, digits = 2),
      I2 = round((OR_A$I2*100), digits = 2),
      T2 = round(OR_A$tau2, digits = 2),
      Pooled.OR_CI = paste(round(exp(OR_A$TE.random), digits = 2),'(',round(exp(OR_A$lower.random), digits = 2),'-',round(exp(OR_A$upper.random), digits = 2),')'),
      pvalue = round(OR_A$pval.random, digits = 2),
      Eggers = round(ifelse(length(A.et$p) == 0,0, A.et$p),digits = 2)
    )
  }
  all_model <- rbind(all_model,data_A)
  s= paste(i,Allele,sep = '/')
  s=paste(s,'txt',sep = '.')
  #sink(s)
  print(OR_A)
  #sink()
  
  
  
  # calculating OR with 95% CI for dominant model
  
  Dominant <- paste(i, "OR_D", sep=Cancer)
  
  
  OR_D <- metabin(Case_Events_D,Case_Total_D, Control_Events_D, Control_Total_D, title = "Dominant Model",
                  data = data_dominant,
                  studlab = paste(Author_Name,'(',Year,')',sep = ''),
                  comb.fixed = TRUE,
                  comb.random = TRUE,
                  method.tau = "SJ",
                  hakn = TRUE,
                  prediction = TRUE,
                  incr = 0.1,
                  sm = "OR")
  if (is.na(OR_D["I2"])){
    OR_D <- "This Value is NA"
    data_D <- data.frame(
      RSID = arsid,
      Gene = agene,
      Model = "NA",
      Study = "NA",
      OR_CI = "NA",
      Weight = "NA",
      Cochran.Q = "NA",
      I2 = "NA",
      T2 = "NA",
      Pooled.OR_CI = "NA",
      pvalue = "NA",
      Eggers = 'NA'
    )
  }else if (OR_D$I2 <= 0.50) {
    OR_D <- metabin(Case_Events_D,Case_Total_D, Control_Events_D, Control_Total_D, title = "Dominant Model",
                    data = data_dominant,
                    studlab = paste(Author_Name,'(',Year,')',sep = ''),
                    comb.fixed = TRUE,
                    comb.random = FALSE,
                    method.tau = "SJ",
                    hakn = TRUE,
                    prediction = TRUE,
                    incr = 0.1,
                    sm = "OR")
    
    D.et <- eggers.test(x = OR_D)
    
    data_D <- data.frame(
      RSID = arsid,
      Gene = agene,
      Model = "Fixed",
      Study = OR_D$studlab,
      OR_CI = paste(round(exp(OR_D$TE), 
                          digits = 2),'(',round(exp(OR_D$lower), digits = 2),'-',
                    round(exp(OR_D$upper), digits = 2),')'),
      Weight = paste(round(((OR_D$w.fixed/sum(OR_D$w.fixed))*100), digits = 2), '%') ,
      Cochran.Q = round(OR_D$pval.Q, digits = 2),
      I2 = round((OR_D$I2*100), digits = 2),
      T2 = round(OR_D$tau2, digits = 2),
      Pooled.OR_CI = paste(round(exp(OR_D$TE.fixed), digits = 2),'(',round(exp(OR_D$lower.fixed), digits = 2),'-',round(exp(OR_D$upper.fixed), digits = 2),')'),
      pvalue = round(OR_D$pval.fixed, digits = 2),
      Eggers = round(ifelse(length(D.et$p) == 0,0, D.et$p),digits = 2)
    )
    
  } else {
    OR_D <- metabin(Case_Events_D,Case_Total_D, Control_Events_D, Control_Total_D, title = "Dominant Model",
                    data = data_dominant,
                    studlab = paste(Author_Name,'(',Year,')',sep = ''),
                    comb.fixed = FALSE,
                    comb.random = TRUE,
                    method.tau = "SJ",
                    hakn = TRUE,
                    prediction = TRUE,
                    incr = 0.1,
                    sm = "OR")
    
    D.et <- eggers.test(x = OR_D)
    
    data_D <- data.frame(
      RSID = arsid,
      Gene = agene,
      Model = "Random",
      Study = OR_D$studlab,
      OR_CI = paste(round(exp(OR_D$TE), 
                          digits = 2),'(',round(exp(OR_D$lower), digits = 2),'-',
                    round(exp(OR_D$upper), digits = 2),')'),
      Weight = paste(round(((OR_D$w.random/sum(OR_D$w.random))*100), digits = 2), '%') ,
      Cochran.Q = round(OR_D$pval.Q, digits = 2),
      I2 = round((OR_D$I2*100), digits = 2),
      T2 = round(OR_D$tau2, digits = 2),
      Pooled.OR_CI = paste(round(exp(OR_D$TE.random), digits = 2),'(',round(exp(OR_D$lower.random), digits = 2),'-',round(exp(OR_D$upper.random), digits = 2),')'),
      pvalue = round(OR_D$pval.random, digits = 2),
      Eggers = round(ifelse(length(D.et$p) == 0,0, D.et$p),digits = 2)
    )
    
  }
  data_dom <- rbind(data_dom,data_D)
  s2= paste(i,Dominant,sep = '/')
  s2=paste(s2,'txt',sep = '.')
  #sink(s2)
  print(OR_D)
  #sink()
  
  # calculating OR with 95% CI for Recessive model
  
  Recessive <- paste(i, "OR_R", sep=Cancer)
  
  OR_R <- metabin(Case_Events_R,Case_Total_R, Control_Events_R, Control_Total_R, title = "Recessive Model",
                  data = data_recessive,
                  studlab = paste(Author_Name,'(',Year,')',sep = ''),
                  comb.fixed = TRUE,
                  comb.random = TRUE,
                  method.tau = "SJ",
                  hakn = TRUE,
                  prediction = TRUE,
                  incr = 0.1,
                  sm = "OR")
  
  if (is.na(OR_R["I2"])){
    OR_R <- "This Value is NA"
    data_R <- data.frame(
      RSID = arsid,
      Gene = agene,
      Model = "NA",
      Study = "NA",
      OR_CI = "NA",
      Weight = "NA",
      Cochran.Q = "NA",
      I2 = "NA",
      T2 = "NA",
      Pooled.OR_CI = "NA",
      pvalue = "NA",
      Eggers = 'NA'
    )
    
  } else if (OR_R$I2 <= 0.50) {
    OR_R <- metabin(Case_Events_R,Case_Total_R, Control_Events_R, Control_Total_R, title = "Recessive Model",
                    data = data_recessive,
                    studlab = paste(Author_Name,'(',Year,')',sep = ''),
                    comb.fixed = TRUE,
                    comb.random = FALSE,
                    method.tau = "SJ",
                    hakn = TRUE,
                    prediction = TRUE,
                    incr = 0.1,
                    sm = "OR")
    
    R.et <- eggers.test(x = OR_R)
    
    data_R <- data.frame(
      RSID = arsid,
      Gene = agene,
      Model = "Fixed",
      Study = OR_R$studlab,
      OR_CI = paste(round(exp(OR_R$TE), 
                          digits = 2),'(',round(exp(OR_R$lower), digits = 2),'-',
                    round(exp(OR_R$upper), digits = 2),')'),
      Weight = paste(round(((OR_R$w.fixed/sum(OR_R$w.fixed))*100), digits = 2), '%') ,
      Cochran.Q = round(OR_R$pval.Q, digits = 2),
      I2 = round((OR_R$I2*100), digits = 2),
      T2 = round(OR_R$tau2, digits = 2),
      Pooled.OR_CI = paste(round(exp(OR_R$TE.fixed), digits = 2),'(',round(exp(OR_R$lower.fixed), digits = 2),'-',round(exp(OR_R$upper.fixed), digits = 2),')'),
      pvalue = round(OR_R$pval.fixed, digits = 2),
      Eggers = round(ifelse(length(R.et$p) == 0,0, R.et$p),digits = 2)
    )
    
  } else {
    OR_R <- metabin(Case_Events_R,Case_Total_R, Control_Events_R, Control_Total_R, title = "Recessive Model",
                    data = data_recessive,
                    studlab = paste(Author_Name,'(',Year,')',sep = ''),
                    comb.fixed = FALSE,
                    comb.random = TRUE,
                    method.tau = "SJ",
                    hakn = TRUE,
                    prediction = TRUE,
                    incr = 0.1,
                    sm = "OR")
    
    R.et <- eggers.test(x = OR_R)
    
    data_R <- data.frame(
      RSID = arsid,
      Gene = agene,
      Model = "Random",
      Study = OR_R$studlab,
      OR_CI = paste(round(exp(OR_R$TE), 
                          digits = 2),'(',round(exp(OR_R$lower), digits = 2),'-',
                    round(exp(OR_R$upper), digits = 2),')'),
      Weight = paste(round(((OR_R$w.random/sum(OR_R$w.random))*100), digits = 2), '%') ,
      Cochran.Q = round(OR_R$pval.Q, digits = 2),
      I2 = round((OR_R$I2*100), digits = 2),
      T2 = round(OR_R$tau2, digits = 2),
      Pooled.OR_CI = paste(round(exp(OR_R$TE.random), digits = 2),'(',round(exp(OR_R$lower.random), digits = 2),'-',round(exp(OR_R$upper.random), digits = 2),')'),
      pvalue = round(OR_R$pval.random, digits = 2),
      Eggers = round(ifelse(length(R.et$p) == 0,0, R.et$p),digits = 2)
    )
    
  }
  data_res <- rbind(data_res,data_R)
  s3= paste(i,Recessive,sep = '/')
  s3=paste(s3,'txt',sep = '.')
  #sink(s3)
  print(OR_R)
  #sink()
  
  # calculating OR with 95% CI for Heterozygous model
  
  Hetero <- paste(i, "OR_He", sep=Cancer)
  
  OR_He <- metabin(Case_Events_He,Case_Total_He, Control_Events_He, Control_Total_He, title = "Heterozygous Model",
                   data=data_heterozygous,
                   studlab = paste(Author_Name,'(',Year,')',sep = ''),
                   comb.fixed = TRUE,
                   comb.random = TRUE,
                   method.tau = "SJ",
                   hakn = TRUE,
                   prediction = TRUE,
                   incr = 0.1,
                   sm = "OR")
  if (is.na(OR_He["I2"])){
    OR_He <- "This Value is NA"
    data_He <- data.frame(
      RSID = arsid,
      Gene = agene,
      Model = "NA",
      Study = "NA",
      OR_CI = "NA",
      Weight = "NA",
      Cochran.Q = "NA",
      I2 = "NA",
      T2 = "NA",
      Pooled.OR_CI = "NA",
      pvalue = "NA",
      Eggers = 'NA'
    )
  }else if (OR_He$I2 <= 0.50) {
    OR_He <- metabin(Case_Events_He,Case_Total_He, Control_Events_He, Control_Total_He, title = "Heterozygous Model",
                     data=data_heterozygous,
                     studlab = paste(Author_Name,'(',Year,')',sep = ''),
                     comb.fixed = TRUE,
                     comb.random = FALSE,
                     method.tau = "SJ",
                     hakn = TRUE,
                     prediction = TRUE,
                     incr = 0.1,
                     sm = "OR")
    
    He.et <- eggers.test(x = OR_He)
    
    data_He <- data.frame(
      RSID = arsid,
      Gene = agene,
      Model = "Fixed",
      Study = OR_He$studlab,
      OR_CI = paste(round(exp(OR_He$TE), 
                          digits = 2),'(',round(exp(OR_He$lower), digits = 2),'-',
                    round(exp(OR_He$upper), digits = 2),')'),
      Weight = paste(round(((OR_He$w.fixed/sum(OR_He$w.fixed))*100), digits = 2), '%') ,
      Cochran.Q = round(OR_He$pval.Q, digits = 2),
      I2 = round((OR_He$I2*100), digits = 2),
      T2 = round(OR_He$tau2, digits = 2),
      Pooled.OR_CI = paste(round(exp(OR_He$TE.fixed), digits = 2),'(',round(exp(OR_He$lower.fixed), digits = 2),'-',round(exp(OR_He$upper.fixed), digits = 2),')'),
      pvalue = round(OR_He$pval.fixed, digits = 2),
      Eggers = round(ifelse(length(He.et$p) == 0,0, He.et$p),digits = 2)
    )
    
  } else {
    OR_He <- metabin(Case_Events_He,Case_Total_He, Control_Events_He, Control_Total_He, title = "Heterozygous Model",
                     data=data_heterozygous,
                     studlab = paste(Author_Name,'(',Year,')',sep = ''),
                     comb.fixed = FALSE,
                     comb.random = TRUE,
                     method.tau = "SJ",
                     hakn = TRUE,
                     prediction = TRUE,
                     incr = 0.1,
                     sm = "OR")
    
    He.et <- eggers.test(x = OR_He)
    
    data_He <- data.frame(
      RSID = arsid,
      Gene = agene,
      Model = "Random",
      Study = OR_He$studlab,
      OR_CI = paste(round(exp(OR_He$TE), 
                          digits = 2),'(',round(exp(OR_He$lower), digits = 2),'-',
                    round(exp(OR_He$upper), digits = 2),')'),
      Weight = paste(round(((OR_He$w.random/sum(OR_He$w.random))*100), digits = 2), '%') ,
      Cochran.Q = round(OR_He$pval.Q, digits = 2),
      I2 = round((OR_He$I2*100), digits = 2),
      T2 = round(OR_He$tau2, digits = 2),
      Pooled.OR_CI = paste(round(exp(OR_He$TE.random), digits = 2),'(',round(exp(OR_He$lower.random), digits = 2),'-',round(exp(OR_He$upper.random), digits = 2),')'),
      pvalue = round(OR_He$pval.random, digits = 2),
      Eggers = round(ifelse(length(He.et$p) == 0,0, He.et$p),digits = 2)
    )
    
  }
  data_het <- rbind(data_het,data_He)
  s4= paste(i,Hetero,sep = '/')
  s4=paste(s4,'txt',sep = '.')
  #sink(s4)
  print(OR_He)
  #sink()
  
  # calculating OR with 95% CI for Homozygous model
  
  Homo <- paste(i, "OR_Ho", sep=Cancer)
  
  OR_Ho <- metabin(Case_Events_Ho,Case_Total_Ho, Control_Events_Ho, Control_Total_Ho, title = "Homozygous Model",
                   data=data_homozygous,
                   studlab = paste(Author_Name,'(',Year,')',sep = ''),
                   comb.fixed = TRUE,
                   comb.random = TRUE,
                   method.tau = "SJ",
                   hakn = TRUE,
                   prediction = TRUE,
                   incr = 0.1,
                   sm = "OR")
  if (is.na(OR_Ho["I2"])){
    OR_Ho <- "This Value is NA"
    data_Ho <- data.frame(
      RSID = arsid,
      Gene = agene,
      Model = "NA",
      Study = "NA",
      OR_CI = "NA",
      Weight = "NA",
      Cochran.Q = "NA",
      I2 = "NA",
      T2 = "NA",
      Pooled.OR_CI = "NA",
      pvalue = "NA",
      Eggers = 'NA'
    )
  }else if (OR_Ho$I2 <= 0.50) {
    OR_Ho <- metabin(Case_Events_Ho,Case_Total_Ho, Control_Events_Ho, Control_Total_Ho, title = "Homozygous Model",
                     data=data_homozygous,
                     studlab = paste(Author_Name,'(',Year,')',sep = ''),
                     comb.fixed = TRUE,
                     comb.random = FALSE,
                     method.tau = "SJ",
                     hakn = TRUE,
                     prediction = TRUE,
                     incr = 0.1,
                     sm = "OR")
    
    Ho.et <- eggers.test(x = OR_Ho)
    
    data_Ho <- data.frame(
      RSID = arsid,
      Gene = agene,
      Model = "Fixed",
      Study = OR_Ho$studlab,
      OR_CI = paste(round(exp(OR_Ho$TE), 
                          digits = 2),'(',round(exp(OR_Ho$lower), digits = 2),'-',
                    round(exp(OR_Ho$upper), digits = 2),')'),
      Weight = paste(round(((OR_Ho$w.fixed/sum(OR_Ho$w.fixed))*100), digits = 2), '%') ,
      Cochran.Q = round(OR_Ho$pval.Q, digits = 2),
      I2 = round((OR_Ho$I2*100), digits = 2),
      T2 = round(OR_Ho$tau2, digits = 2),
      Pooled.OR_CI = paste(round(exp(OR_Ho$TE.fixed), digits = 2),'(',round(exp(OR_Ho$lower.fixed), digits = 2),'-',round(exp(OR_Ho$upper.fixed), digits = 2),')'),
      pvalue = round(OR_Ho$pval.fixed, digits = 2),
      Eggers = round(ifelse(length(Ho.et$p) == 0,0, Ho.et$p),digits = 2)
    )
    
  } else {
    OR_Ho <- metabin(Case_Events_Ho,Case_Total_Ho, Control_Events_Ho, Control_Total_Ho, title = "Homozygous Model",
                     data=data_homozygous,
                     studlab = paste(Author_Name,'(',Year,')',sep = ''),
                     comb.fixed = FALSE,
                     comb.random = TRUE,
                     method.tau = "SJ",
                     hakn = TRUE,
                     prediction = TRUE,
                     incr = 0.1,
                     sm = "OR")
    
    Ho.et <- eggers.test(x = OR_Ho)
    
    data_Ho <- data.frame(
      RSID = arsid,
      Gene = agene,
      Model = "Random",
      Study = OR_Ho$studlab,
      OR_CI = paste(round(exp(OR_Ho$TE), 
                          digits = 2),'(',round(exp(OR_Ho$lower), digits = 2),'-',
                    round(exp(OR_Ho$upper), digits = 2),')'),
      Weight = paste(round(((OR_Ho$w.random/sum(OR_Ho$w.random))*100), digits = 2), '%') ,
      Cochran.Q = round(OR_Ho$pval.Q, digits = 2),
      I2 = round((OR_Ho$I2*100), digits = 2),
      T2 = round(OR_Ho$tau2, digits = 2),
      Pooled.OR_CI = paste(round(exp(OR_Ho$TE.random), digits = 2),'(',round(exp(OR_Ho$lower.random), digits = 2),'-',round(exp(OR_Ho$upper.random), digits = 2),')'),
      pvalue = round(OR_Ho$pval.random, digits = 2),
      Eggers = round(ifelse(length(Ho.et$p) == 0,0, Ho.et$p),digits = 2)
    )
    
  }
  data_hom <- rbind(data_hom,data_Ho)
  s5= paste(i,Homo,sep = '/')
  s5=paste(s5,'txt',sep = '.')
  #sink(s5)
  print(OR_Ho)
  #sink()
  
  # generating forest plot for Allele model
  if (length > 12){
  h = length/10/1.2*500}
  else{
  h=500}
  
  if (is.na(OR_A["I2"])){
    OR_A <- "This Value is NA"
    figA <- paste(i, "forest_A", sep=Cancer)
    s5= paste(i,figA,sep = '/')
    s5=paste(s5,'txt',sep = '.')
    #sink(s5)
    print(OR_A)
    #sink()
  }else{
    Allele2 = paste(i,'_forest_A',sep = Cancer)
    figA= paste(Cancer,Allele2,sep = '/')
    figA=paste(figA,'jpg',sep = '.')
    
    jpeg(filename = figA,width = 1300, height = h, res = 120)
    fplot_A <- forest(OR_A, smlab = "Allele Model",lab.e = "Case", test.overall = TRUE,
                      col.diamond.fixed = "darkblue", col.diamond.random = "lightseagreen", print.stat = TRUE)
    dev.off()
  }
  # generating forest plot for Dominant model
  if (is.na(OR_D["I2"])){
    OR_D <- "This Value is NA"
    figD <- paste(i, "_forest_D", sep=Cancer)
    s5= paste(i,figD,sep = '/')
    s5=paste(s5,'txt',sep = '.')
    #sink(s5)
    print(OR_D)
    #sink()
  }else{
    Dominant2 = paste(i,'_forest_D',sep = Cancer)
    figD= paste(Cancer,Dominant2,sep = '/')
    figD=paste(figD,'jpg',sep = '.')
    
    jpeg(filename = figD,width = 1300, height = h, res = 120)
    fplot_D <- forest(OR_D, smlab = "Dominant Model",lab.e = "Case", test.overall = TRUE,
                      col.diamond.fixed = "darkblue", col.diamond.random = "lightseagreen", print.stat = TRUE)
    dev.off()
  }
  # generating forest plot for Recessive model
  if (is.na(OR_R["I2"])){
    OR_R <- "This Value is NA"
    figR <- paste(i, "_forest_R", sep=Cancer)
    s5= paste(i,figR,sep = '/')
    s5=paste(s5,'txt',sep = '.')
    #sink(s5)
    print(OR_R)
    #sink()
  }else{
    Recessive2 = paste(i,'_forest_R',sep = Cancer)
    figR= paste(Cancer,Recessive2,sep = '/')
    figR=paste(figR,'jpg',sep = '.')
    
    jpeg(filename = figR,width = 1300, height = h, res = 120)
    fplot_R <- forest(OR_R, smlab = "Recessive Model",lab.e = "Case", test.overall = TRUE,
                      col.diamond.fixed = "darkblue", col.diamond.random = "lightseagreen", print.stat = TRUE)
    dev.off()
  }
  # generating forest plot for Heterozygous model
  if (is.na(OR_He["I2"])){
    OR_He <- "This Value is NA"
    figHe <- paste(i, "forest_He", sep=Cancer)
    s5= paste(i,figHe,sep = '/')
    s5=paste(s5,'txt',sep = '.')
    #sink(s5)
    print(OR_He)
    #sink()
  }else{
    Hetero2 = paste(i,'_forest_He',sep = Cancer)
    figHe= paste(Cancer,Hetero2,sep = '/')
    figHe=paste(figHe,'jpg',sep = '.')
    
    jpeg(filename =figHe,width = 1300, height = h, res = 120)
    fplot_He <- forest(OR_He, smlab = "Heterozygous Model",lab.e = "Case", test.overall = TRUE,
                       col.diamond.fixed = "darkblue", col.diamond.random = "lightseagreen", print.stat = TRUE)
    dev.off()
  }
  # generating forest plot for Homozygous model
  if (is.na(OR_Ho["I2"])){
    OR_Ho <- "This Value is NA"
    figHo <- paste(i, "_forest_Ho", sep=Cancer)
    s5= paste(i,figHo,sep = '/')
    s5=paste(s5,'txt',sep = '.')
    #sink(s5)
    print(OR_Ho)
    #sink()
  }else{
    Homo2 = paste(i,'_forest_Ho',sep = Cancer)
    figHo= paste(Cancer,Homo2,sep = '/')
    figHo=paste(figHo,'jpg',sep = '.')
    
    jpeg(filename = figHo,width = 1300, height = h, res = 120)
    fplot_Ho <- forest(OR_Ho, smlab = "Homozygous Model",lab.e = "Case", test.overall = TRUE,
                       col.diamond.fixed = "darkblue", col.diamond.random = "lightseagreen", print.stat = TRUE)
    dev.off()
  }
  # Eggers test and funnel plots for all models
  
  # A.et <- eggers.test(x = OR_A)
  # A.et$p
  # D.et <- eggers.test(x = OR_D)
  # D.et$p
  # R.et <- eggers.test(x = OR_R)
  # R.et$p
  # He.et <- eggers.test(x = OR_He)
  # He.et$p
  # Ho.et <- eggers.test(x = OR_Ho)
  # Ho.et$p
  # 
  
  
  
  
  
  #data_eggers_allele <- data.frame()
  #Allele Model Funnel Plot:
  if (is.na(OR_A["I2"])){
    OR_A <- "This Value is NA"
    figA <- paste(i, "_funnel_A", sep=Cancer)
    s5= paste(i,figA,sep = '/')
    s5=paste(s5,'txt',sep = '.')
    #sink(s5)
    print(OR_A)
    #sink()
    
    data_eggers_allele <- data.frame(
      RSID = arsid,
      Gene = agene,
      Eggers = 'NA'
    )
    
  }else{
    A.et <- eggers.test(x = OR_A)
    #A.et$p
    Allele3 = paste(i,'_funnel_AA',sep = Cancer)
    funA= paste(Cancer,Allele3,sep = '/')
    funA=paste(funA,'jpg',sep = '.')
    
    jpeg(filename = funA,width = 450, height = 400, res = 120)
    funnel(x = OR_A, studlab = TRUE)
    title(A.et$p)
    dev.off()
    
    data_eggers_allele <- data.frame(
      RSID = arsid,
      Gene = agene,
      Eggers = ifelse(length(A.et$p) == 0,0, A.et$p)
    )
  }
  data_eggers_allele_all <- rbind(data_eggers_allele_all, data_eggers_allele)
  #Dominant Model Funnel Plot:
  if (is.na(OR_D["I2"])){
    OR_D <- "This Value is NA"
    figD <- paste(i, "_funnel_D", sep=Cancer)
    s5= paste(i,figD,sep = '/')
    s5=paste(s5,'txt',sep = '.')
    #sink(s5)
    print(OR_D)
    #sink()
    
    data_eggers_dominant <- data.frame(
      RSID = arsid,
      Gene = agene,
      Eggers = 'NA'
    )
  }else{
    D.et <- eggers.test(x = OR_D)
    D.et$p
    Dominant3 = paste(i,'_funnel_D',sep = Cancer)
    funD= paste(Cancer,Dominant3,sep = '/')
    funD=paste(funD,'jpg',sep = '.')
    
    
    jpeg(filename = funD,width = 450, height = 400, res = 120)
    funnel(x = OR_D, studlab = TRUE)
    title(D.et$p)
    dev.off() 
    
    data_eggers_dominant <- data.frame(
      RSID = arsid,
      Gene = agene,
      Eggers = ifelse(length(D.et$p) == 0,0, D.et$p)
    )
  }
  data_eggers_dominant_all <- rbind(data_eggers_dominant_all, data_eggers_dominant)
  
  #Recessive Model Funnel Plot:
  if (is.na(OR_R["I2"])){
    OR_R <- "This Value is NA"
    figR <- paste(i, "funnel_R", sep=Cancer)
    s5= paste(i,figR,sep = '/')
    s5=paste(s5,'txt',sep = '.')
    #sink(s5)
    print(OR_R)
    #sink()
    
    data_eggers_recessive <- data.frame(
      RSID = arsid,
      Gene = agene,
      Eggers = 'NA'
    )
    
  }else{
    R.et <- eggers.test(x = OR_R)
    R.et$p
    Recessive3 = paste(i,'_funnel_R',sep = Cancer)
    funR= paste(Cancer,Recessive3,sep = '/')
    funR=paste(funR,'jpg',sep = '.')
    
    
    jpeg(filename = funR,width = 450, height = 400, res = 120)
    funnel(x = OR_R, studlab = TRUE)
    title(R.et$p)
    dev.off()
    
    data_eggers_recessive <- data.frame(
      RSID = arsid,
      Gene = agene,
      Eggers = ifelse(length(R.et$p) == 0,0, R.et$p)
    )
  }
  data_eggers_recessive_all <- rbind(data_eggers_recessive_all, data_eggers_recessive)
  
  #Heterozygous Model Funnel Plot:
  if (is.na(OR_He["I2"])){
    OR_He <- "This Value is NA"
    figHe <- paste(i, "funnel_He", sep=Cancer)
    s5= paste(i,figHe,sep = '/')
    s5=paste(s5,'txt',sep = '.')
    #sink(s5)
    print(OR_He)
    #sink()
    
    data_eggers_Heterozygous <- data.frame(
      RSID = arsid,
      Gene = agene,
      Eggers = 'NA'
    )
    
  }else{
    He.et <- eggers.test(x = OR_He)
    He.et$p
    Hetero3 = paste(i,'_funnel_He',sep = Cancer)
    funHe= paste(Cancer,Hetero3,sep = '/')
    funHe=paste(funHe,'jpg',sep = '.')
    
    
    jpeg(filename = funHe,width = 450, height = 400, res = 120)
    funnel(x = OR_He, studlab = TRUE)
    title(He.et$p)
    dev.off()
    
    data_eggers_Heterozygous <- data.frame(
      RSID = arsid,
      Gene = agene,
      Eggers = ifelse(length(He.et$p) == 0,0, He.et$p)
    )
  }
  data_eggers_Heterozygous_all <- rbind(data_eggers_Heterozygous_all, data_eggers_Heterozygous)
  
  #Homozygous Model Funnel Plot:
  if (is.na(OR_Ho["I2"])){
    OR_Ho <- "This Value is NA"
    figHo <- paste(i, "funnel_Ho", sep=Cancer)
    s5= paste(i,figHo,sep = '/')
    s5=paste(s5,'txt',sep = '.')
    #sink(s5)
    print(OR_Ho)
    #sink()
    
    data_eggers_Homozygous <- data.frame(
      RSID = arsid,
      Gene = agene,
      Eggers = 'NA'
    )
  }else{
    Ho.et <- eggers.test(x = OR_Ho)
    Ho.et$p
    Homo3 = paste(i,'_funnel_Ho',sep = Cancer)
    funHo= paste(Cancer,Homo3,sep = '/')
    funHo=paste(funHo,'jpg',sep = '.')
    
    
    jpeg(filename = funHo,width = 450, height = 400, res = 120)
    funnel(x =OR_Ho, studlab = TRUE)
    title(Ho.et$p)
    dev.off()
    
    data_eggers_Homozygous <- data.frame(
      RSID = arsid,
      Gene = agene,
      Eggers = ifelse(length(Ho.et$p) == 0,0, Ho.et$p)
    )
  }
  
  data_eggers_Homozygous_all <- rbind(data_eggers_Homozygous_all, data_eggers_Homozygous)
  
}

#write.xlsx(all_model, file = "Model_Table.xlsx", sheetName = "Allele", row.names = FALSE)
#write.xlsx(data_dom, file = "Model_Table.xlsx", sheetName = "Dominant", append = TRUE, row.names = FALSE)
#write.xlsx(data_res, file = "Model_Table.xlsx", sheetName = "Recessive", append = TRUE, row.names = FALSE)
#write.xlsx(data_het, file = "Model_Table.xlsx", sheetName = "Heterozygous", append = TRUE, row.names = FALSE)
#write.xlsx(data_hom, file = "Model_Table.xlsx", sheetName = "Homozygous", append = TRUE, row.names = FALSE)

#y1=paste("Eggers_test_all",'xlsx',sep = '.')
#write.xlsx(data_eggers_allele_all, file = y1, sheetName = "Allele",row.names = FALSE)
#write.xlsx(data_eggers_dominant_all, file = y1, sheetName = "Dominant", append = TRUE,row.names = FALSE)
#write.xlsx(data_eggers_recessive_all, file = y1, sheetName = "Recessive", append = TRUE,row.names = FALSE)
#write.xlsx(data_eggers_Heterozygous_all, file = y1, sheetName = "Heterozygous", append = TRUE,row.names = FALSE)
#write.xlsx(data_eggers_Homozygous_all, file = y1, sheetName = "Homozygous", append = TRUE,row.names = FALSE)


