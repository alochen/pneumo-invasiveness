###########################################################################################################################################################
###########################################################################################################################################################
############################## Script for case-to-carrier likelihood model ################################################################################
###########################################################################################################################################################
###########################################################################################################################################################

setwd("Q:/Technical/R/Case-to-carrier")
library(plyr)
library(ggplot2)
library(dplyr)
library(tidyr)
library(Rmisc)
library(grid)
library(gridExtra)
library(gibbs.met)
library(RColorBrewer)
library(purrr)
library(readxl)
library(reshape)
library(DescTools)
library(gmp)
library(lemon)
library(ggpubr)
library(forcats)
library(stringr)

# upload output files
# abs.child <- read.csv("Q:/Technical/R/Case-to-carrier/abs.child-new.csv")
# abs.child$X <- NULL
# abs.adult <- read.csv("abs.adult-new.csv")
# abs.adult$X <- NULL
# consol.child <- read.csv("consol.child-new.csv")
# consol.child$X <- NULL
# consol.adult <- read.csv("consol.adult-new.csv")
# consol.adult$X <- NULL

infant <- read.csv("Q:/Technical/R/Case-to-carrier/pneumo_invasiveness_inf.csv")
adult <- read.csv("Q:/Technical/R/Case-to-carrier/pneumo_invasiveness_adu.csv")

# Weinberger priors
# inv_priors_df <- read.csv("Q:/Technical/R/Case-to-carrier/Weinberger_invpriors.csv")
# inv_priors_df <- inv_priors_df %>% mutate(stddev = sqrt(1/log.inv.prec.age1))

#### Rename column names of each dataset -------------------------------------------------------------------------------------------------------------------
colnames(infant)[2:24] <- paste("carriage_",colnames(infant)[2:24], sep = "")
colnames(infant)[26:48] <- paste("disease_",colnames(infant)[26:48], sep = "")
colnames(infant)[26:48] <- substr(colnames(infant)[26:48],1, nchar(colnames(infant)[26:48])-2)

colnames(adult)[2:8] <- paste("carriage_",colnames(adult)[2:8], sep = "")
colnames(adult)[10:16] <- paste("disease_",colnames(adult)[10:16], sep = "")
colnames(adult)[10:16] <- substr(colnames(adult)[10:16],1, nchar(colnames(adult)[10:16])-2)

Portugaldiseasedat <- read.csv("Q:/Technical/R/Case-to-carrier/portugal_invasivedis_aggregage.csv")
adult$disease_Portugal.pre.PCV <- NULL
adult <- full_join(adult, Portugaldiseasedat[,c(1,3)])
colnames(adult)[length(adult)] <- "disease_Portugal.pre.PCV"
adult$disease_Portugal.pre.PCV[is.na(adult$disease_Portugal.pre.PCV)] <- 0
adult$carriage_Portugal.pre.PCV[is.na(adult$carriage_Portugal.pre.PCV)] <- 0

infant<- full_join(infant, Portugaldiseasedat[,c(1,2)])
colnames(infant)[length(infant)] <- "disease_Portugal.pre.PCV"
infant$disease_Portugal.pre.PCV[is.na(infant$disease_Portugal.pre.PCV)] <- 0
infant <- full_join(infant, adult[,c('Serotype', 'carriage_Portugal.pre.PCV')])
infant$carriage_Portugal.pre.PCV[is.na(infant$carriage_Portugal.pre.PCV)] <- 0

NLprePCV <- read.csv('Q:/Technical/R/Case-to-carrier/NL_prePCV.csv')
colnames(NLprePCV)[c(3,4)] <- c('carriage_Netherlands.pre.PCV', 'disease_Netherlands.pre.PCV')
NLprePCV$X <- NULL
NLpostPCV7 <- read.csv('Q:/Technical/R/Case-to-carrier/NL_postPCV7.csv')
colnames(NLpostPCV7)[c(3,4)] <- c('carriage_Netherlands.post.PCV7', 'disease_Netherlands.post.PCV7')
NLpostPCV7$X <- NULL
NLpostPCV10 <- read.csv('Q:/Technical/R/Case-to-carrier/NL_postPCV10.csv')
colnames(NLpostPCV10)[c(3,4)] <- c('carriage_Netherlands.post.PCV10', 'disease_Netherlands.post.PCV10')
NLpostPCV10$X <- NULL
NLds <- full_join(NLprePCV, NLpostPCV7, by = 'Serotype') %>% full_join(NLpostPCV10, by = 'Serotype')
NLds[is.na(NLds)] <- 0
NLds <- NLds[!str_detect(NLds$Serotype, ','),]

# join NL and infant DS
infant <- full_join(infant, NLds, by = 'Serotype')
infant[is.na(infant)] <- 0

#adult[is.na(adult)] <- 0
#infant[is.na(infant)] <- 0

#### -------------------------------------------------------------------------------------------------------------------------------------------------------
#### -------------------------------------------------------------------------------------------------------------------------------------------------------
#### Function to get rid of second serotype column in raw data, separate/gather data by carriage or disease, add serogroup column, eliminate ---------------
#### serotypes with 0 cases of carriage and disease in each dataset ----------------------------------------------------------------------------------------
#### -------------------------------------------------------------------------------------------------------------------------------------------------------
#### -------------------------------------------------------------------------------------------------------------------------------------------------------

initial.bigdata.cleanup <- function(initial.data) {

  data.name <- deparse(substitute(initial.data))
  
  # initial data cleanup: remove second serotype column in raw data, split number of cases by carriers or disease for each dataset--------------------------
  initial.data <- initial.data[, -which(names(initial.data) %in% c("Serotype.1"))] #remove second serotype label column
  initial.data <- tidyr::gather(initial.data, "dataset", "num.cases", 2:ncol(initial.data)) # make each dataset column into a row
  initial.data <- initial.data %>% separate(dataset, c("dis.carriage", "DS"), sep = "\\_", remove = TRUE) %>%
    spread(dis.carriage, num.cases) # column for # disease and column for # carriage
  
  # add serogroup column -----------------------------------------------------------------------------------------------------------------------------------
  VT7_index <- which(initial.data$Serotype %in% c("4", "6B", "9V", "14", "18C", "19F", "23F"))
  VT10_index <- which(initial.data$Serotype %in% c("1", "5", "7F"))
  VT13_index <- which(initial.data$Serotype %in% c("3", "6A", "19A"))
  VT_index <- c(VT7_index, VT10_index, VT13_index)
  
  initial.data$Serogroup <- rep("NVT", nrow(initial.data))
  initial.data$Serogroup[VT7_index] <- "VT7"
  initial.data$Serogroup[VT10_index] <- "VT10"
  initial.data$Serogroup[VT13_index] <- "VT13"
  
  # eliminate serotypes in each dataset with 0 carriers and 0 cases (i.e. no data) -------------------------------------------------------------------------
  initial.data <- initial.data %>% dplyr::mutate(sum.dis.carr = disease + carriage)
  ind.zero <- which(initial.data$sum.dis.carr %in% 0)
  initial.data <- initial.data[-ind.zero, ]
  initial.data$sum.dis.carr <- NULL
  
  #### return function -------------------------------------------------------------------------------------------------------------------------------------
  return(initial.data)
}

new.infant <- initial.bigdata.cleanup(infant)
new.infant$agegrp <- rep("children", nrow(new.infant))
new.adult <- initial.bigdata.cleanup(adult)
new.adult$agegrp <- rep("adults", nrow(new.adult))
tot.ds <- dplyr::bind_rows(new.infant, new.adult)
tot.ds <- tot.ds[,c(6,2,5,1,3,4)] # rearrange columns into following order: age group, dataset, serogroup, serotype, carriage, disease
tot.ds$Serotype <- as.character(tot.ds$Serotype)
tot.ds <- tot.ds[which(tot.ds$Serotype != "NT"),] # remove NT
tot.ds <- tot.ds[complete.cases(tot.ds), ] # removes 19C NA in non-Portugal DS

#### Make population dataframe with study time interval, n.swab, population numbers ------------------------------------------------------------------------
# sources of population numbers available upon request #

DS <- c("Alabama.pre.PCV", "Alaska.pre.PCV", "Barcelona.post.PCV7", "Bogota.post.PCV7", "Bogota.pre.PCV", "Caracas.pre.PCV", "Czech.pre.PCV","France.post.PCV7",
        "Goroka.pre.PCV", "Kenya.pre.PCV", "Navajo.post.PCV7", "Oxford.pre.PCV", "Sweden.pre.PCV", "E.W.pre.PCV", "France.post.PCV13", "Massachusetts.post.PCV7", 
        "Atlanta.post.PCV7", "Atlanta.pre.PCV", "Finland.pre.PCV", "Greece.pre.PCV", "Ontario.pre.PCV", "Iceland.pre.PCV", "Morocco.pre.PCV", "Portugal.pre.PCV",
        "Stockholm.pre.PCV", 'Netherlands.pre.PCV', 'Netherlands.post.PCV7', 'Netherlands.post.PCV10')
time.int <- c(3.4, 5, 4, 1, 4, 1, 9, 2, # France wrong: 6
              6, 4/2, 6, 7, 8, 1, 2, 3, # France wrong: 6
              1, 1, 5, 5, 1, 9, 1, 2,
              1, 2, 4, 4)
n.swab <- c(827, NA, 209, 246, 197, 1004, 425, 1212, 
            2844, NA, 6541, NA, 550, 3752, 1212, 2969, # Ox originally 501 (wrong, now excluded)
            451, 231, 3290, 464, 1139, NA, 200, 1170,
            611, 321, 660, 659)
N.children <- c(19316, NA, 228000, 357200, 357200, 146125, 478177, 838866, #cz wrong: 784863, france wrong: 1592988,
                96207, NA, 65048, 37467, 397289, 3091000, 842076, 820000, #oxford wrong: 7489; goroka wrong: 3957; navajo wrong: 30558; france wrong: 1597519, sweden wrong: 376838; oxford wrong:3778636
                298831, 204680, 318083, NA, 580507, NA, 212566, 2071223,
                NA, 250924, 232251, 222671)
N.adults <- c(232373, NA, NA, NA, NA, NA, NA, NA, 
              NA, NA, 201553, NA, 7066020, 48702414, NA, NA, #navajo wrong: 108789
              NA, NA, NA, NA, NA, NA, NA, 8284894,
              2004152, NA, NA, NA)

population.data <- data.frame(DS, n.swab, N.children, N.adults, time.int)

#### Absolute rates data clean up --------------------------------------------------------------------------------------------------------------------------

## Data cleaning for children carriage vs children disease

child_only <- population.data %>% filter(!is.na(N.children)) %>% filter(!is.na(n.swab)) %>% dplyr::select(-4) # keep only popln data for which children have data on n.swabs and N.children
# write.csv(child_only, "children_DS_description.csv") # DESCRIPTION OF DATASETS
colnames(child_only)[3] <- "N" # change column name of imported DF so it's easier to append population values onto main dataframe

keepers <- paste(c(as.character(child_only$DS)), collapse = '|') # reformat dataset names that we are keeping 

child_df <- tot.ds %>% group_by(agegrp, DS) %>% nest() %>%
  filter(!grepl("adult", agegrp)) %>% # remove adult age group
  filter(grepl(keepers, DS)) %>% # keep only the keepers
  unnest(c(data))

child.abs.ds <- dplyr::left_join(child_df, child_only, by = "DS") %>% filter(Serotype != "TOTAL") # append dataframes together + remove total, working dataframe
# remove Finland and Sweden
child.abs.ds <- child.abs.ds %>% filter(!DS %in% c('Finland.pre.PCV', 'Sweden.pre.PCV'))

# Table: which serotypes are in each DS?
childdssero <- reshape2::dcast(child.abs.ds %>% select(DS, Serotype),formula = Serotype ~ DS)
childdssero <- cbind(childdssero$Serotype, data.frame(ifelse(is.na(childdssero[,2:ncol(childdssero)]),'','X')))
colnames(childdssero)[1] <- 'Serotype'
childdssero <- childdssero[,c(1, 2,4,7:10,13,15,19:21, 3,5,6,12,14,16,18, 11,17)]
write.csv(childdssero, 'children_whichseroineachDS.csv')

## Data cleaning for children carriage vs adults disease

adult_only <- population.data %>% filter(!is.na(N.adults)) %>% filter(!is.na(n.swab)) %>% dplyr::select(-3)
colnames(adult_only)[3] <- "N"

keepers.adult <- paste(c(as.character(adult_only$DS)), collapse = '|')

adult_df <- tot.ds %>% group_by(agegrp, DS) %>% nest() %>%
  filter(!grepl("children", agegrp)) %>%
  filter(grepl(keepers.adult, DS)) %>% 
  unnest(c(data))

adult.abs.ds <- dplyr::left_join(adult_df, adult_only, by = "DS") %>% filter(Serotype != "TOTAL") # append dataframes together + remove total, working dataframe
adult.abs.ds <- adult.abs.ds %>% filter(!DS %in% c('Finland.pre.PCV', 'Sweden.pre.PCV'))

# Table: which serotypes are in each DS?
adultdssero <- reshape2::dcast(adult.abs.ds %>% select(DS, Serotype),formula = Serotype ~ DS)
adultdssero <- cbind(adultdssero$Serotype, data.frame(ifelse(is.na(adultdssero[,2:ncol(adultdssero)]),'','X')))
colnames(adultdssero)[1] <- 'Serotype'
adultdssero <- adultdssero[,c(1,2,3,5,6,4)]
write.csv(adultdssero, 'adults_whichseroineachDS.csv')

#### -------------------------------------------------------------------------------------------------------------------------------------------------------
#### -------------------------------------------------------------------------------------------------------------------------------------------------------
#### Function to perform Markov chain sampling for carriage prevalence and IPD incidence distribution and then estimate invasiveness for each serotype------
#### + Function to clean invasiveness data and draw a heatmap of the invasiveness per serotype and dataset/study--------------------------------------------
#### -------------------------------------------------------------------------------------------------------------------------------------------------------
#### -------------------------------------------------------------------------------------------------------------------------------------------------------
#### -------------------------------------------------------------------------------------------------------------------------------------------------------

invasiveness <- function(row_val) {
  set.seed(123)
  
  agegroup <- row_val[["agegrp"]]
  DS <- row_val[["DS"]] # dataset/study location
  serogroup <- row_val[["Serogroup"]]
  serotype <- as.character(row_val[["Serotype"]])
  n.carriage <- as.numeric(row_val[["carriage"]]) # number of carriers
  n.disease <- as.numeric(row_val[["disease"]]) # number of disease cases
  N <- as.numeric(row_val[["N"]]) # population
  n.swab <- as.numeric(row_val[["n.swab"]]) # number of swabs
  time.int <- as.numeric(row_val[["time.int"]]) # time interval
  
  # carriage prevalence estimate
  carr_prev_samp <- seq(0,(n.carriage/n.swab)+0.05, length.out = 560000)
  carr_prev_LL <- dbeta(x = carr_prev_samp, shape1 = n.carriage + 1, shape2 = n.swab - n.carriage + 1, log = TRUE)
  carr_prev_MLE <- carr_prev_samp[which.max(carr_prev_LL)]
  
  carr_prev_maxLL <- max(carr_prev_LL)
  carr_prev_prob.dens <- exp(carr_prev_LL-carr_prev_maxLL)
  #plot(carr_prev_samp,carr_prev_prob.dens, type = "l", main = serotype, xlab = "carr_prev", ylab = "Probability Density")
  
  carr_prev_cp <- cumsum(carr_prev_prob.dens) / sum(carr_prev_prob.dens)
  carr_prev.low <- carr_prev_samp[which(carr_prev_cp >= 0.025)[1]]
  carr_prev.high <- carr_prev_samp[which(carr_prev_cp >= 0.975)[1]]
  
  # disease hazard estimate
  lambda_samp <- seq(0,(n.disease*10)+10, by = 0.01)#length.out = 560000)
  lambda_LL <- dgamma(x = lambda_samp, shape = n.disease + 1, rate = 1, log = TRUE)
  lambda_MLE <- lambda_samp[which.max(lambda_LL)]
  
  lambda_maxLL <- max(lambda_LL)
  lambda_prob.dens <- exp(lambda_LL-lambda_maxLL)
  #plot(lambda_samp,lambda_prob.dens, type = "l", main = serotype, xlab = "lambda", ylab = "Probability Density")

  lambda_cp <- cumsum(lambda_prob.dens) / sum(lambda_prob.dens)
  lambda.low <- lambda_samp[which(lambda_cp >= 0.025)[1]]
  lambda.high <- lambda_samp[which(lambda_cp >= 0.975)[1]]
  
  ######## inv and credible intervals in new bin way 30/08/2019 -----------------------------------------------------------------------------
  
  set.seed(123)
  num.dev <- 100000
  
  carr_prev_dev <- rbeta(n = num.dev, shape1 = n.carriage + 1, shape2 = n.swab - n.carriage + 1)
  
  # # prior distribution i.e. weinberger
  # if (serotype %in% inv_priors_df$st) 
  #   {invas_dev <- rlnorm(n = num.dev, meanlog = inv_priors_df$log.inv.age1[which(inv_priors_df$st == serotype)]-log(1000), 
  #                     sdlog = 3*inv_priors_df$stddev[which(inv_priors_df$st == serotype)])
  # } else {invas_dev <- runif(n = num.dev, min = 0.00001, max = 0.5)}
  
  # uniform prior distribution between 0 and 1
  if (serotype %in% c("1","46")) {
    inv.max = 1.4
  } else {
    inv.max = 1
  }

  invas_dev <- runif(n = num.dev, min = 0.0000001, max = inv.max)
  
  # # LOG PRIOR DISTRIBUTION
  # invas_dev <- runif(n = num.dev, min = -9,max = 1)
  # invas_dev <- 10^invas_dev
  
  # posterior distribution
  lambda <- carr_prev_dev*invas_dev*N*time.int
  posterior_samp <- dpois(n.disease, lambda = carr_prev_dev*invas_dev*N*time.int)
  bins <- seq(0, max(invas_dev), length.out = 5000)
  bin.width <- bins[2] - bins[1]
  inv.cut <- cut(invas_dev, bins)
  dist.sum <- tapply(posterior_samp, INDEX = inv.cut, FUN = sum)
  dist.mean <- dist.sum/num.dev
  dist.mean[is.na(dist.mean)] <- 0

  invas.MLE <- sum(bins[which.max(dist.mean)],bins[which.max(dist.mean)+1])/2
  
  ## likelihood; posterior distribution if prior is flat/uninformative ------------------------------------------------------------------------------------
  # 
  # invas_flat_dev <- runif(n = num.dev, min = 0.00001, max = 0.1)
  # lambda_flat <- carr_prev_dev*invas_flat_dev*N*time.int
  # bins.flat <- seq(0, max(invas_flat_dev), length.out = 1000)
  # bin.width.flat <- bins.flat[2] - bins.flat[1]
  # inv.cut.flat <- cut(invas_flat_dev, bins.flat)
  # posterior_samp_flat <- dpois(n.disease, lambda = lambda_flat)
  # dist.sum.flat <- tapply(posterior_samp_flat, INDEX = inv.cut.flat, FUN = sum)
  # dist.mean.flat <- dist.sum.flat/num.dev
  # dist.mean.flat[is.na(dist.mean.flat)] <- 0
  # 
  # # prior distribution
  # if (serotype %in% inv_priors_df$st)
  # {something <- dlnorm(invas_dev, meanlog = inv_priors_df$log.inv.age1[which(inv_priors_df$st == serotype)]-log(1000), #lognormal distribution
  #                     sdlog = 3*inv_priors_df$stddev[which(inv_priors_df$st == serotype)])
  # } else {something <- dunif(invas_dev, min = 0.00001, max = 0.2)}
  # 
  # plot(invas_dev, something*max(dist.mean)/max(something), col = "blue", main = paste(DS, agegroup, serotype, sep = " "),
  #      xlim = c(0, max(invas_dev, bins.flat[which.max(dist.mean.flat)]))) # lognormal distribution (prior)
  # lines(bins.flat[-length(bins.flat)], dist.mean.flat*max(dist.mean)/max(dist.mean.flat), type = "s") # likelihood (flat prior)
  # lines(bins[-length(bins)], dist.mean, type = "s", col = "red") # posterior distribution (w informative prior)
  
  #### bins & credible intervals of posterior inv ----
  bin.df <- data.frame(bins = bins[-length(bins)], distrib = dist.mean, DS = rep(DS, nrow(dist.mean)))
  invas_cp <- cumsum(bin.df$distrib/sum(bin.df$distrib))
  invas.low <- bin.df$bins[which(invas_cp >= 0.025)[1]] #sum(bin.df$bins[which(invas_cp >= 0.025)[1]],bin.df$bins[which(invas_cp >= 0.025)[2]])/2
  invas.high <- sum(bin.df$bins[which(invas_cp >= 0.975)[1]],bin.df$bins[which(invas_cp >= 0.975)[2]])/2
  
  # to plot vs weinberger distributions; and for model comparison
  #bin.df.new <- data.frame(bins = invas_dev, distrib = posterior_samp) # new bin.df
  #bin.df.new <- bin.df.new[order(bin.df.new$bins),]
  df.name <- paste(DS, "distrib", agegroup, serotype, sep = ".")
  assign(df.name, bin.df, envir = .GlobalEnv) #bin.df.new
   
  
  # ####### the original way -----------------------------------------------------------------------------------------------------------------
  # # invasiveness credible intervals
  # set.seed(123)
  # # deviates
  # carr_prev_dev <- rbeta(n = 560000, shape1 = n.carriage + 1, shape2 = n.swab - n.carriage + 1)
  # lambda_dev <- rgamma(n = 560000, shape = n.disease + 1, rate = 1)
  # 
  # carr_prev_dev_LL <- dbeta(carr_prev_dev, shape1 = n.carriage + 1, shape2 = n.swab - n.carriage + 1, log = TRUE)
  # lambda_dev_LL <- dgamma(lambda_dev, shape = n.disease + 1, rate = 1, log = TRUE)
  # 
  # invas_dev2 <- lambda_dev/(carr_prev_dev*N*time.int)
  # invas_dev_LL2 <- carr_prev_dev_LL+lambda_dev_LL
  # invas_maxLL2 <- max(invas_dev_LL2)
  # invas_prob.dens2 <- exp(invas_dev_LL2-invas_maxLL2)
  # new.invas_dev2 <- invas_dev2[order(invas_dev2)]
  # new.invas_prob.dens2 <- invas_prob.dens2[order(invas_dev2)]
  # plot(new.invas_dev2,new.invas_prob.dens2, xlab = "invas", ylab = "Probability Density", type = "l")
  # 
  # invas_cp <- cumsum(new.invas_prob.dens) / sum(new.invas_prob.dens)
  # invas.low <- new.invas_dev[which(invas_cp >= 0.025)[1]]
  # invas.high <- new.invas_dev[which(invas_cp >= 0.975)[1]]
  # 
  if (carr_prev_MLE == 0) {carr_prev.low <- 0}
  if (lambda_MLE == 0) {lambda.low <- 0}
  if (invas.MLE == 0) {invas.low <- 0}
  # #if (is.infinite(invas.MLE) == TRUE) {invas.low <- 0 # take away infinite invasiveness credible intervals
  # #invas.high <- 0}
  
  # return function ----------------------------------------------------------------------------------------------------------------------------------------
  
  new_row <- data.frame(row_val, carr.prev = carr_prev_MLE, carr.prev.low = carr_prev.low, carr.prev.high = carr_prev.high, 
                        lambda = lambda_MLE, lambda.low = lambda.low, lambda.high = lambda.high, 
                        invasiveness = invas.MLE, invasiveness.low = invas.low, invasiveness.high = invas.high)
  return(new_row)
  
}

# Functions that estimates one invasiveness per serotype, consolidating all datasets. This is done by estimating a serotype-specific carriage prevalence for each
# DS, using it as a prior for estimating serotype-specific invasiveness. Disease distribution is estimated in the same way but this time we make the substitution 
# disease rate = invasiveness * carriage prevalence with carriage prevalence being the prior.

overall_inv_2 <- function(df, invas) { # function that, given an invasiveness calculates carriage prev and lambda deviates 
                                       # for each dataset and then sums the mean log likelihood
  set.seed(123)
  n <- nrow(df)
  
  carr_prev_dev <- sapply(seq(1:n), function(s) {rbeta(n = 10000, shape1 = df[s,5] + 1, shape2 = df[s,7] - df[s,5] + 1)})
  lambda <- sweep(carr_prev_dev*invas, 2, df[,8]*df[,9], FUN = "*")
  L <- colMeans(sapply(seq(1:n), function(s) {dpois(df[s, 6], lambda[,s])}))
  LL <- sum(log(L))
  L.prime <- exp(LL)

  return(L.prime)
}

calc_invas <- function(df, invas) { # function that returns a single invasiveness with credible intervals for a serotype
                             # given all the datasets that contain it
  invas <- invas
        #seq(0.00000001, 0.001, # children pre and postPCV; adults
        #seq(0.00000001, 0.002, # children prePCV sero 14, 27;
                                # children postPCV sero 4, 8, 10B, 33F, 38, 24F, 18C, 18A, 3, 14
        #seq(0.00000001, 0.02, #0.01, # children prePCV sero 4, 7F, 18B, 13, 31, 12F; 
                               # adults sero 1, 4, 8, 20, 24F, 17F, 35A; 
                               # children postPCV sero 5, 7F, 12F, 25A, 12B
        #seq(0.00000001, 0.1, #children prePCV sero 1, 5
                             #children postPCV 33A, 33F, 37
        #seq(0.00000001, 0.04, # children postPCV sero 1
        #length.out = 2000)
  likelihood <- sapply(invas, function(s) overall_inv_2(df, s))
  #plot(invas, LL, main = df$Serotype[1])
  
  Serotype <- as.character(unique(df$Serotype))
  
  # Weinberger prior
  # if (Serotype %in% inv_priors_df$st) {
  # prior_invas <- dlnorm(invas, meanlog = inv_priors_df$log.inv.age1[which(inv_priors_df$st == Serotype)]-log(1000),
  #                              sdlog = inv_priors_df$stddev[which(inv_priors_df$st == Serotype)])
  # } else {
   prior_invas <- dunif(invas, min = 0.00001, max = 0.5)
  #}

  posteriormarginal <- likelihood*prior_invas
  maxposterior.index <- which.max(posteriormarginal) # (LL)
  invas.opt <- invas[maxposterior.index]
  
  cp <- cumsum(posteriormarginal / sum(posteriormarginal))
  invas.low <- invas[which(cp >= 0.025)[1]]
  invas.high <- invas[which(cp >= 0.975)[1]]
  
  inv_row <- data.frame(Serotype = df$Serotype[1], overall.inv = invas.opt, inv.low = invas.low, inv.high = invas.high)
  
  #plot(invas, likelihood, type = "l", main = df$Serotype[1])
  # plots to compare with weinberger's distributions
  likelihooddf <- data.frame(invas = invas, likelihood = likelihood, distribution = "Global invasiveness")
  weinbergerdf <- data.frame(invas = invas, likelihood = prior_invas*max(likelihood)/max(prior_invas), distribution = "Weinberger et al")
  consol.df <- bind_rows(likelihooddf, weinbergerdf)
  plot <- ggplot(consol.df) + geom_line(aes(x = invas, y = likelihood, colour = distribution)) +
  coord_cartesian(xlim = c(0,0.1)) +
   ggtitle(Serotype) + theme_classic() + labs(x = "Invasiveness", y = "Probability Density") +
  scale_y_continuous(labels = scales::scientific) + theme_classic() + theme(legend.position = "none")


  bin.df <- data.frame(bins = invas, distrib = likelihood)
  df.name <- paste(paste("sero", df$Serotype[1], sep = ""), "distrib", df$agegrp[1], sep = ".")
  assign(df.name, bin.df, envir = .GlobalEnv)
  
  return(inv_row) #plot) #
}

#### -------------------------------------------------------------------------------------------------------------------------------------------------------
#### CODE FOR PLOTS   --------------------------------------------------------------------------------------------------------------------------------------
#### -------------------------------------------------------------------------------------------------------------------------------------------------------
#cols <- c("VT7" = "#70AD47", "VT10"= "#4472C4", "VT13" = "#ED7D31", "NVT" = "black")
cols <- c("VT7" = "#35B779FF", "VT10"= "#31688EFF", "VT13" = "#E69F00", "NVT" = "#440154FF")

g_legend <- function(a.gplot){ # function for legend grob (same as get_legend from ggpubr package)
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)}

plot.inv.allages.allDS <- function(sero) {
  VT7 <- c("4","6B","9V", "14", "18C", "19F", "23F")
  VT10 <- c("1","5","7F")
  VT13 <- c("3", "6A", "19A")
  NVT <- as.character(unique(sero$Serotype[!sero$Serotype %in% c(VT7, VT10, VT13)]))

  sero$Serotype <- factor(sero$Serotype, levels = c(VT7, VT10, VT13, NVT))
  sero$Serogroup <- factor(sero$Serogroup, levels = c("VT7", "VT10", "VT13", "NVT"))
  sero$agegrp <- factor(sero$agegrp, levels = c("children","adults"), labels = c("Children", "Adults"))
  sero$vaccinationera <- factor(sero$vaccinationera, levels = c("pre.PCV", "post.PCV7", "post.PCV10/13"), labels = c("Pre-PCV", "Post-PCV7", "Post-PCV10/13"))
  sero$DS <- factor(sero$DS, levels = c("Alabama.pre.PCV", "Atlanta.pre.PCV", "Bogota.pre.PCV","Caracas.pre.PCV","Czech.pre.PCV","E.W.pre.PCV",#"Finland.pre.PCV", 
                                        "Goroka.pre.PCV","Morocco.pre.PCV", "Netherlands.pre.PCV", "Ontario.pre.PCV",#"Oxford.pre.PCV",
                                        "Portugal.pre.PCV","Stockholm.pre.PCV",#"Sweden.pre.PCV",  
                                        "Atlanta.post.PCV7", "Barcelona.post.PCV7", "Bogota.post.PCV7","France.post.PCV7","Massachusetts.post.PCV7","Navajo.post.PCV7",
                                        "Netherlands.post.PCV7",
                                        "France.post.PCV13", "Netherlands.post.PCV10" ),
                    labels = c("Alabama", "Atlanta", "Bogota","Caracas","Czech","E&W",#"Finland", 
                               "Goroka","Morocco", "Netherlands", "Ontario",#"Oxford",
                               "Portugal","Stockholm",#"Sweden",  
                               "Atlanta ", "Barcelona", "Bogota ","France","Massachusetts","Navajo", "Netherlands ",
                               "France ", "Netherlands  "))
  
  ggplot(sero, aes(x = Serotype, y = invasiveness, colour = Serogroup, # x = Serotype
                   group = vaccinationera)) + #interaction(Serogroup, agegrp))) + #
   geom_errorbar(aes(ymin = invasiveness.low, ymax = invasiveness.high), 
                 position = position_dodge(0.5),
                  width = 0.01) + 
    geom_point(aes(shape = vaccinationera), #agegrp), #
               size = 2, position = position_dodge(0.5)) +
    labs(x =  "", #"Dataset", #
         #y = "Invasiveness (case per carrier per yr)", #title = "Invasiveness by dataset", #"Study-specific invasiveness", #
         color = "Serotype category", 
         shape = "Vaccination era") + #"Age group") +
    ylab(expression(Invasiveness~(case~carrier^{-1}~year^{-1}))) +
    #guides(color = guide_legend(order = 1), shape = guide_legend(order = 2)) + 
    scale_color_manual(values= cols) +
    scale_shape_manual(values = c(16, 0, 17)) +
    scale_y_log10() +
    # comment out for x = Serotypes
    # scale_x_discrete(labels = ##parse(text = levels(sero$DS))) +
    #                    c(expression(Alabama[0]), expression(Atlanta[0]), expression(Bogota[0]),
    #                      expression(Caracas[0]),expression(Czech[0]),expression(E*"&"*W[0]),
    #                      #expression(Finland[0]),
    #                      expression(Goroka[0]),expression(Morocco[0]), expression(Netherlands[0]),
    #                      expression(Ontario[0]),expression(Portugal[0]),
    #                      expression(Stockholm[0]),
    #                      #expression(Sweden[0]),
    #                      expression(Atlanta[1] ), expression(Barcelona[1]), expression(Bogota[1] ),
    #                      expression(France[1]),expression(Massachusetts[1]),expression(Navajo[1]), expression(Netherlands[1]),
    #                      expression(France[2]), expression(Netherlands[2]))) +
    theme_bw() + theme_light() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5), # sero angle = 90, hjust = 1, vjust = 0.5; DS angle = 30, hjust = 1
          text = element_text(size = 17)) #poster text size  25
}
#
plot.carr.allages.allDS <- function(sero) { # plots carriage prevalence
  VT7 <- c("4","6B","9V", "14", "18C", "19F", "23F")
  VT10 <- c("1","5","7F")
  VT13 <- c("3", "6A", "19A")
  NVT <- as.character(unique(sero$Serotype[!sero$Serotype %in% c(VT7, VT10, VT13)]))
  sero$DS <- factor(sero$DS, levels = c("Alabama.pre.PCV", "Atlanta.pre.PCV", "Bogota.pre.PCV","Caracas.pre.PCV","Czech.pre.PCV","E.W.pre.PCV",#"Finland.pre.PCV", 
                                      "Goroka.pre.PCV","Morocco.pre.PCV", "Netherlands.pre.PCV", "Ontario.pre.PCV",#"Oxford.pre.PCV",
                                      "Portugal.pre.PCV","Stockholm.pre.PCV",#"Sweden.pre.PCV",  
                                      "Atlanta.post.PCV7", "Barcelona.post.PCV7", "Bogota.post.PCV7","France.post.PCV7","Massachusetts.post.PCV7","Navajo.post.PCV7",
                                      "Netherlands.post.PCV7",
                                      "France.post.PCV13", "Netherlands.post.PCV10" ),
                   labels = c("Alabama", "Atlanta", "Bogota","Caracas","Czech","E&W",#"Finland", 
                              "Goroka","Morocco", "Netherlands", "Ontario",#"Oxford",
                              "Portugal","Stockholm",#"Sweden",  
                              "Atlanta ", "Barcelona", "Bogota ","France","Massachusetts","Navajo", "Netherlands ",
                              "France ", "Netherlands  "))
  sero$Serotype <- factor(sero$Serotype, levels = c(VT7, VT10, VT13, NVT))
  sero$Serogroup <- factor(sero$Serogroup, levels = c("VT7", "VT10", "VT13", "NVT"))
  sero$vaccinationera <- factor(sero$vaccinationera, levels = c("pre.PCV", "post.PCV7", "post.PCV10/13"), labels = c("Pre-PCV", "Post-PCV7", "Post-PCV10/13"))

  ggplot(sero, aes(x = DS, y = carr.prev)) + #x = DS
    geom_errorbar(aes(ymin = carr.prev.low, ymax = carr.prev.high, group = Serogroup, color = Serogroup), # change grouping (group = Serogroup for DS, group = vaccinationera for Sero)
                  width = 0.01,
                  position = position_dodge(0.5), 
                  stat = 'identity') +
    geom_point(aes(color = Serogroup, shape = vaccinationera, group = Serogroup), #
               position = position_dodge(0.5), 
               stat = 'identity',
               size = 2) + 
    labs(x = "", y = "Carriage prevalence", 
         color = "Serotype category",
         shape = "Vaccination era") +
    scale_color_manual(values=cols)+#c("#70AD47", "#4472C4", "#ED7D31", "black")) + # PCV7 green, PCV10 blue, PCV13 orange, NVT black
    scale_shape_manual(values = c(16, 0, 17)) +
    
    # comment out for x = Serotypes
    scale_x_discrete(labels = ##parse(text = levels(sero$DS))) +
                       c(expression(Alabama[0]), expression(Atlanta[0]), expression(Bogota[0]),
                         expression(Caracas[0]),expression(Czech[0]),expression(E*"&"*W[0]),
                         #expression(Finland[0]),
                         expression(Goroka[0]),expression(Morocco[0]), expression(Netherlands[0]),
                         expression(Ontario[0]),expression(Portugal[0]),
                         expression(Stockholm[0]),
                         #expression(Sweden[0]),
                         expression(Atlanta[1] ), expression(Barcelona[1]), expression(Bogota[1] ),
                         expression(France[1]),expression(Massachusetts[1]),expression(Navajo[1]), expression(Netherlands[1]),
                         expression(France[2]), expression(Netherlands[2]))) +
    
    #scale_y_continuous(trans = 'log10', labels = scales::comma) +
    theme_bw() + theme_light() +
    theme(axis.text.x = element_text(angle = 30, hjust = 1), # DS: angle = 30, hjust = 1; sero: angle = 90, hjust = 1, vjust = 0.5
          text = element_text(size = 17))
}
#
plot.dis.allages.allDS <- function(sero) { # plots all datasets' IPD incidence for all serotypes in a graph
  VT7 <- c("4","6B","9V", "14", "18C", "19F", "23F")
  VT10 <- c("1","5","7F")
  VT13 <- c("3", "6A", "19A")
  NVT <- as.character(unique(sero$Serotype[!sero$Serotype %in% c(VT7, VT10, VT13)]))
  sero$Serotype <- factor(sero$Serotype, levels = c(VT7, VT10, VT13, NVT))
  
  # depending on whether loading Rdata or running for first time, these levels need to change (pre.PCV <--> pre-PCV)
  sero$DS <- factor(sero$DS, levels = c("Alabama.pre.PCV", "Atlanta.pre.PCV", "Bogota.pre.PCV","Caracas.pre.PCV","Czech.pre.PCV","E.W.pre.PCV",#"Finland.pre-PCV",
                                        "Goroka.pre.PCV","Morocco.pre.PCV", "Netherlands.pre.PCV", "Ontario.pre.PCV",#"Oxford.pre.PCV",
                                        "Portugal.pre.PCV","Stockholm.pre.PCV",#"Sweden.pre.PCV",
                                        "Atlanta.post.PCV7", "Barcelona.post.PCV7", "Bogota.post.PCV7","France.post.PCV7","Massachusetts.post.PCV7","Navajo.post.PCV7",
                                        "Netherlands.post.PCV7",
                                        "France.post.PCV13", "Netherlands.post.PCV10" ))#,
                    # labels = c("Alabama", "Atlanta", "Bogota","Caracas","Czech","E&W",#"Finland",
                    #            "Goroka","Morocco", "Netherlands", "Ontario","Portugal","Stockholm",#"Sweden",
                    #            "Atlanta ", "Barcelona", "Bogota ","France","Massachusetts","Navajo", "Netherlands ",
                    #            "France ", "Netherlands  "))
  
  sero$IPDinc <- (sero$lambda/(sero$N*sero$time.int))*100000

  sero$agegrp <- factor(sero$agegrp, levels = c("children","adults"), labels = c("Children", "Adults"))
  sero$Serogroup <- factor(sero$Serogroup, levels = c("VT7", "VT10", "VT13", "NVT"))
  ggplot(sero, aes(x = Serotype, #DS, #
                   y = IPDinc, colour = Serogroup, group = interaction(Serogroup, agegrp))) +
    geom_errorbar(aes(ymin = (sero$lambda.low/(sero$N*sero$time.int))*100000, 
                      ymax = (sero$lambda.high/(sero$N*sero$time.int))*100000, 
                      color = Serogroup), width = 0, 
                      position = position_dodge(0.5)) +
                      #position = position_jitterdodge(jitter.width = 0.7, dodge.width = 1, seed = 1)) +
    geom_point(aes(color = Serogroup, shape = agegrp),
               size = 2, #size = 5, 
               position = position_dodge(0.5)) +
               #position = position_jitterdodge(jitter.width = 0.7, dodge.width = 1, seed = 1)) +
    labs(x = "", 
         y = "IPD incidence per 100,000 ppl", 
         color = "Serotype category", shape = "Age group", size = element_blank()) +
    scale_color_manual(values= cols) + 
    scale_y_continuous(trans = 'log10') + #labels = scales::comma
    scale_shape_manual(values = c(16, 0)) +
    
    # comment out for x = Serotypes
    # scale_x_discrete(labels = ##parse(text = levels(sero$DS))) +
    #                    c(expression(Alabama[0]), expression(Atlanta[0]), expression(Bogota[0]),
    #                      expression(Caracas[0]),expression(Czech[0]),expression(E*"&"*W[0]),
    #                      #expression(Finland[0]),
    #                      expression(Goroka[0]),expression(Morocco[0]), expression(Netherlands[0]),
    #                      expression(Ontario[0]),expression(Portugal[0]),
    #                      expression(Stockholm[0]),
    #                      #expression(Sweden[0]),
    #                      expression(Atlanta[1] ), expression(Barcelona[1]), expression(Bogota[1] ),
    #                      expression(France[1]),expression(Massachusetts[1]),expression(Navajo[1]), expression(Netherlands[1]),
    #                      expression(France[2]), expression(Netherlands[2]))) +
    theme_bw() + theme_light() + guides(size = FALSE) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5), # sero: angle = 90, hjust = 1, vjust = 0.5; DS: angle = 30, hjust = 1
          text = element_text(size = 17))
}
#
plot.consol.inv <- function(sero) { # plots the overall/consolidated invasiveness of all serotypes in one graph
  #sero$agegrp <- factor(sero$agegrp, levels = c("children","adults"), labels = c("Children", "Adults"))
  sero$Serogroup <- factor(sero$Serogroup, levels = c("VT7", "VT10", "VT13","NVT"))
  
  ggplot(sero, aes(x = reorder(Serotype, -overall.invasiveness), y = overall.invasiveness))+ 
    geom_errorbar(aes(
      ymin = overall.invasiveness.low, ymax = overall.invasiveness.high, color = Serogroup),#, group = agegrp), 
                  width = 0, position = position_dodge(0.4)) + # width = 0.01 for paper; width = 0.025 for poster
    geom_point(aes(
      color = Serogroup),#, shape = agegrp, group = agegrp), 
               position = position_dodge(0.4), size = 3) +
    labs(x = "Serotype", color = "Serotype category",shape = "Age group") +
    ylab(expression(Invasiveness~(case~carrier^{-1}~year^{-1}))) +
    scale_color_manual(values= cols) +
    scale_y_continuous(trans = 'log10') +
    theme_classic() + theme_bw() + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5), text = element_text(size = 17)) # text size 30 for poster
} # pdf 5x12

# # setwd("C:/Users/kmcraemc/Box/ITS Team/Phase II - AD/Projects/Peri-threshold research paper/Risk score")
# tiff("posterfig1.tiff", height=480, width=2*480, units="px")
# plot.consol.inv(consol.tot)
# dev.off()

plot.DS.overall <- function(sero, consolidated.inv) { # plots one DS invasiveness with overall/consolidated invasiveness
  
  #### scatterplot w overall invasiveness vs DS-specific invasiveness
  new.sero <- sero %>% select(DS, Serogroup, Serotype, invasiveness, invasiveness.low, invasiveness.high)
  new.sero <- data.frame(lapply(new.sero, FUN = unlist))
  
  consolidated.inv.dat <- consolidated.inv %>% filter(Serotype %in% new.sero$Serotype)
  consolidated.inv.dat$DS <- rep("overall", nrow(consolidated.inv.dat))
  consolidated.inv.dat$agegrp <- NULL
  consolidated.inv.dat <- data.frame(lapply(consolidated.inv.dat, FUN = unlist))
  consolidated.inv.dat$X <- NULL
  colnames(consolidated.inv.dat) <- c("Serotype", "invasiveness", "invasiveness.low", "invasiveness.high", "Serogroup", "DS")
  
  new.df <- dplyr::left_join(new.sero, consolidated.inv.dat, by = c("Serotype", "Serogroup"))
  new.df$Serogroup <- factor(new.df$Serogroup, levels = c("VT7", "VT10", "VT13", "NVT"))
  DSnam <- sub("*.p.*","", unique(sero$DS))
  vaccera <- sub(pattern = "\\.", replacement = "-", sub("^.*\\.p","p", unique(sero$DS)))
  
  ggplot(data = new.df, aes(x = invasiveness.x, y = invasiveness.y, color = Serogroup)) + 
    geom_point() + # (size = 5.5) + if geom_text included 
    geom_abline(slope = 1, intercept = 0) +
    geom_errorbar(aes(ymin = invasiveness.low.y, ymax = invasiveness.high.y), width = 0.0001) +
    geom_errorbarh(aes(xmin = invasiveness.low.x, xmax = invasiveness.high.x), height = 0.001) +
    scale_y_continuous(trans = 'log10') +
    scale_x_continuous(trans = 'log10') +
    scale_color_manual(values= cols) +
    #geom_text(aes(label = Serotype), size = 3, color = "white", fontface = "bold") +
    labs(x = "Global invasiveness", y = paste(DSnam, vaccera, sep = " ")) +
    theme_bw() + theme_light()

  #### scatterplot w Serotype vs Invasiveness per 100,000 ppl ----
  #x.sero <- unlist(sero$Serotype)
  #y.inv <- unlist(sero$invasiveness)*100000
  #consolidated.inv.dat <- consolidated.inv %>% filter(Serotype %in% x.sero)
    
  #ggplot() +
  #  geom_point(data = sero, 
  #             aes(x = x.sero, y = y.inv), color = "red") + 
  #  geom_errorbar(data = sero, 
  #                aes(x = x.sero, y = y.inv, 
  #                    ymin = unlist(invasiveness.low)*100000, 
  #                    ymax = unlist(invasiveness.high)*100000), width = 0.01) + 
  #  geom_point(data = consolidated.inv.dat, 
  #             aes(x = unlist(Serotype), y = unlist(overall.invasiveness)*100000), color = "blue") + 
  #  geom_errorbar(data = consolidated.inv.dat, 
  #                aes(x = unlist(Serotype), y = overall.invasiveness*100000, 
  #                    ymin = unlist(overall.invasiveness.low)*100000, 
  #                    ymax = unlist(overall.invasiveness.high)*100000), width = 0.01) +
  #  labs(x = "Serotype", y = "Invasiveness per 100,000 ppl", title = paste(sero$DS, "comparison with overall invasiveness", sep = " ")) +
  #  #guides(color = guide_legend(override.aes = list(color = c("red", "blue")))) +
  #  #scale_color_manual(name = "Dataset") + 
  #  #scale_color_discrete(name = "Source", breaks = c("Dataset" = "red", "Overall" = "blue"))  +
  #  #ylim (0, 2000)+
  #  theme_bw() + theme_light()
}
#
#### -------------------------------------------------------------------------------------------------------------------------------------------------------
#### -------------------------------------------------------------------------------------------------------------------------------------------------------
#### Invasiveness Rates ------------------------------------------------------------------------------------------------------------------------------------
#### -------------------------------------------------------------------------------------------------------------------------------------------------------
#### -------------------------------------------------------------------------------------------------------------------------------------------------------
##### Sero and DS specific invasiveness --------------------------------------------------------------------------------------------------------------------
# Invasiveness for children only
abs.child <- as.data.frame(t(sapply(1:nrow(child.abs.ds), function(x) invasiveness(child.abs.ds[x,]))))
abs.child <- data.frame(map(abs.child, unlist))

# Invasiveness for adults relative to children's carriage
abs.adult <- as.data.frame(t(sapply(1:nrow(adult.abs.ds), function(x) invasiveness(adult.abs.ds[x,]))))
abs.adult <- data.frame(map(abs.adult, unlist))

# Invasiveness for children and adults in same DF
abs.tot <- bind_rows(abs.child, abs.adult)
abs.tot$agegrp <- as.factor(abs.tot$agegrp)
abs.tot$DS <- as.factor(abs.tot$DS)
abs.tot$Serotype <- as.factor(abs.tot$Serotype)
abs.tot$vaccinationera <- as.factor(sub("^.*\\.p","p", abs.tot$DS))
abs.tot$vaccinationera <- factor(abs.tot$vaccinationera, levels = c('pre.PCV', 'post.PCV7', 'post.PCV10', 'post.PCV13'))
levels(abs.tot$vaccinationera)[match("post.PCV10", levels(abs.tot$vaccinationera))] <- "post.PCV10/13"
levels(abs.tot$vaccinationera)[match("post.PCV13", levels(abs.tot$vaccinationera))] <- "post.PCV10/13"
abs.child <- abs.tot %>% filter(agegrp == "children")
abs.adult <- abs.tot %>% filter(agegrp == "adults")
write.csv(abs.child, "abs.child-new.csv")
write.csv(abs.adult, "abs.adult-new.csv")

# TABLE: sero ordering from largest inv to smallest (point estimate)
serorank <- arrange(abs.tot, DS, desc(invasiveness)) %>% select(Serotype, DS, invasiveness, agegrp)
sf <- serorank %>% filter(agegrp == 'children') %>% select(-invasiveness, -agegrp) %>% 
  pivot_wider(names_from = DS, values_from = Serotype)
sf_a <- serorank %>% filter(agegrp == 'adults') %>% select(-invasiveness, -agegrp) %>% 
  pivot_wider(names_from = DS, values_from = Serotype)
sf_child <- sapply(map(sf, unlist), '[', seq(max(sapply(map(sf,unlist), length))))
sf_child <- sf_child[,c(1,3,6:9, 12,14,18:20, 2,4,5, 11,13,15,17, 10,16)]
sf_adult <- sapply(map(sf_a, unlist), '[', seq(max(sapply(map(sf_a,unlist), length))))
sf_adult <- sf_adult[,c(1,2,4,5,3)]
write.csv(sf_child, 'sero_rankDS_child.csv')
write.csv(sf_adult, 'sero_rankDS_adult.csv')

# Carriage prevalence of children for all DS in one graph
plot.carr.allages.allDS(abs.tot) + theme(legend.position = 'bottom') # pdf 8 x 15 (by sero, leg bott) and pdf 6x12 (by DS)

# Disease incidence for all DS 
plot.dis.allages.allDS(abs.tot) + theme(legend.position = 'bottom') # pdf 8 x 15 (by sero, leg bott) and pdf 6 by 16 (by DS)

# overall incidence by dataset
childDS.incidence <- abs.child %>% group_by(DS) %>% nest()
childDS.incidencenum <- unlist(lapply(childDS.incidence$data, function(x) (sum(x$disease)/(x$N[1]*x$time.int[1]))*100000))
child_overallincidence <- data.frame(childDS.incidence$DS, childDS.incidencenum)
colnames(child_overallincidence) <- c('DS', 'incidence')
child_overallincidence$DS <- factor(child_overallincidence$DS, 
                                    levels = c("Alabama.pre.PCV", "Atlanta.pre.PCV", "Bogota.pre.PCV","Caracas.pre.PCV","Czech.pre.PCV","E.W.pre.PCV",#"Finland.pre.PCV", 
                                      "Goroka.pre.PCV","Morocco.pre.PCV", "Netherlands.pre.PCV", "Ontario.pre.PCV",#"Oxford.pre.PCV",
                                      "Portugal.pre.PCV","Stockholm.pre.PCV",#"Sweden.pre.PCV",  
                                      "Atlanta.post.PCV7", "Barcelona.post.PCV7", "Bogota.post.PCV7","France.post.PCV7","Massachusetts.post.PCV7","Navajo.post.PCV7",
                                      "Netherlands.post.PCV7",
                                      "France.post.PCV13", "Netherlands.post.PCV10" ),
                                    labels = c("Alabama pre-PCV", "Atlanta pre-PCV", "Bogota pre-PCV","Caracas pre-PCV","Czech pre-PCV","E&W pre-PCV",#"Finland", 
                                               "Goroka pre-PCV","Morocco pre-PCV", "Netherlands pre-PCV", "Ontario pre-PCV",#"Oxford pre-PCV",
                                               "Portugal pre-PCV","Stockholm pre-PCV",#"Sweden pre-PCV",  
                                               "Atlanta post-PCV7", "Barcelona post-PCV7", "Bogota post-PCV7","France post-PCV7","Massachusetts post-PCV7","Navajo post-PCV7", "Netherlands post-PCV7",
                                               "France post-PCV13", "Netherlands post-PCV10"))

ggplot(child_overallincidence) + 
  geom_point(aes(x = DS, y = incidence)) + 
  labs(x = "", y = 'Incidence per 100,000ppl') + 
  #scale_y_log10() +
  theme_bw() + theme_light() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

adultDS.incidence <- abs.adult %>% group_by(DS) %>% nest()
adultDS.incidencenum <- unlist(lapply(adultDS.incidence$data, function(x) (sum(x$disease)/(x$N[1]*x$time.int[1]))*100000))
adult_overallincidence <- data.frame(adultDS.incidence$DS, adultDS.incidencenum)
colnames(adult_overallincidence) <- c('DS', 'incidence')
adult_overallincidence$DS <- factor(adult_overallincidence$DS, 
                                    levels = c("Alabama.pre.PCV", "E.W.pre.PCV", "Portugal.pre.PCV",
                                                "Stockholm.pre.PCV", "Sweden.pre.PCV",  "Navajo.post.PCV7"),
                                    labels = c("Alabama pre-PCV","E&W pre-PCV","Portugal pre-PCV",
                                               "Stockholm pre-PCV","Sweden pre-PCV", "Navajo post-PCV7"))

ggplot(adult_overallincidence) + 
  geom_point(aes(x = DS, y = incidence)) + 
  labs(x = "", y = 'Incidence per 100,000ppl') + 
  scale_y_log10() +
  theme_bw() + theme_light() + 
  theme(axis.text.x = element_text(angle = 30, hjust = 1))

# Invasiveness for all DS

locinvchild <- plot.inv.allages.allDS(abs.child) #sero x axis
locinvadult <- plot.inv.allages.allDS(abs.adult)#sero x axis
fig <- ggarrange(locinvchild + ylab('') + ggtitle('Children (< 18 years)'), 
                 locinvadult + ylab('')+ggtitle('Adults (> 18 years)'), 
                 common.legend = TRUE, legend = 'bottom', nrow = 2)
annotate_figure(fig, left = text_grob(expression(Invasiveness~(case~carrier^{-1}~year^{-1})), rot = 90, size = 17)) # pdf 10x16

byDS <- plot.inv.allages.allDS(abs.tot) # pdf dim 5 x 16 for serotype ; pdf dim 5.5 x 16 for DS
bySero <- plot.inv.allages.allDS(abs.tot)
ggarrange(bySero +theme(legend.position = "none"), byDS + theme(legend.position = "none"), 
          nrow=2, labels = c('A', 'B'), common.legend = TRUE, legend = "right") # pdf 8 x 16

# NEW ANALYSIS - NOT USED: log prior vs uniform prior
# abs.child.old <- read.csv('abs.child-new.csv')
# abs.child.old$X <- NULL
# logprivsunif2 <- data.frame(abs.child$Serogroup, 
#                             abs.child$invasiveness, abs.child$invasiveness.low, abs.child$invasiveness.high,
#                             abs.child.old$invasiveness, abs.child.old$invasiveness.low, 
#                             abs.child.old$invasiveness.high)
#                           
# ggplot(logprivsunif2, aes(x = abs.child.invasiveness, y = abs.child.old.invasiveness)) + 
#   geom_point(aes(color = abs.child.Serogroup)) + 
#   geom_errorbar(aes(ymin = abs.child.old.invasiveness.low, ymax = abs.child.old.invasiveness.high)) + 
#   geom_errorbarh(aes(xmin = abs.child.invasiveness.low, xmax = abs.child.invasiveness.high)) + 
#   scale_x_log10() + scale_y_log10() + geom_abline() + 
#   labs(x = 'Log prior', y= 'Uniform prior', colour = 'Sero category')
# 
# 
# logprivsunif2 <- inner_join(abs.child, abs.child.old, by = colnames(abs.child)[1:15])
# ggplot(logprivsunif2) +
#   geom_point(aes(x = invasiveness.x, y = invasiveness.y)) +
#   geom_errorbar(aes(x = invasiveness.x, y = invasiveness.y, ymin = invasiveness.low.y, ymax = invasiveness.high.y)) +
#   geom_errorbarh(aes(x = invasiveness.x, y = invasiveness.y, xmin = invasiveness.low.x, xmax = invasiveness.high.x))+
#   scale_x_log10() + scale_y_log10() +
#   labs(x = 'Uniform prior', y = 'Log prior')


##### Consolidated serotype-specific invasiveness----------------------------------------------------------------------------------------------------------
# Consolidated invasiveness for children in pre-PCV time periods
inv.child.consol <- child.abs.ds %>% #data.frame(map(child.abs.ds, unlist)) %>% 
  filter(!DS %in% c("France.post.PCV7", "Barcelona.post.PCV7", "France.post.PCV13", "Bogota.post.PCV7", "Massachusetts.post.PCV7","Navajo.post.PCV7", 
                    "Netherlands.post.PCV7", "Netherlands.post.PCV10", "Atlanta.post.PCV7")) %>% 
  group_by(Serotype) %>% nest()
nDS <- sapply(inv.child.consol$data, function(x) nrow(x)) # number of datasets for each serotype

excl.sero <- inv.child.consol$Serotype[which(nDS == 1)] # serotypes with only 1 dataset to be excluded for analysis
inv.child.consol <- data.frame(inv.child.consol %>% unnest(cols = c(data)) %>% filter(!Serotype %in% excl.sero))
# eachglobalinvDS <- inv.child.consol %>% select(Serotype, DS)
# write.csv(eachglobalinvDS, 'globalinv_eachDS_children.csv') # Table of which datasets included for each serotype global inv calculation
unique.child.sero <- unique(inv.child.consol$Serotype)
#unique.child.sero <- unique.child.sero[!unique.child.sero %in% c(sero.to.change.low, sero.to.change.high, sero.to.change.higher)]

consol.child2 <- lapply(unique.child.sero, function(x) {calc_invas(inv.child.consol %>% filter(Serotype == x), 
                                                                   invas = seq(0.00000001, 0.001, length.out = 2000))})
consol.child <- do.call("rbind", consol.child2)

# specific serotypes that require manually changing the calc_invas function sequence
sero.to.change.low <- c("14")#, "27")
consol.child.increasedbounds0 <- lapply(sero.to.change.low, function(x) {calc_invas(inv.child.consol %>% filter(Serotype == x), 
                                                                                                                invas = seq(0.00000001, 0.005,length.out = 2000))})
consol.child.increasedbounds0 <- do.call("rbind", consol.child.increasedbounds0)
sero.to.change.high <- c("4", "7F", "18B", "31", "13", "12F")
consol.child.increasedbounds1 <- lapply(sero.to.change.high, function(x) {calc_invas(inv.child.consol %>% filter(Serotype == x), 
                                                                                                                 invas = seq(0.00000001, 0.02, length.out = 2000))})
consol.child.increasedbounds1 <- do.call("rbind", consol.child.increasedbounds1)
sero.to.change.higher <- c("1", "5")
consol.child.increasedbounds2 <- lapply(sero.to.change.higher, function(x) {calc_invas(inv.child.consol %>% filter(Serotype == x), 
                                                                                                                   invas = seq(0.00000001, 0.1, length.out = 2000))})
consol.child.increasedbounds2 <- do.call("rbind", consol.child.increasedbounds2)

consol.child <- consol.child[!(consol.child$Serotype %in% c(sero.to.change.low, sero.to.change.high, sero.to.change.higher)),]
consol.child <- bind_rows(consol.child, consol.child.increasedbounds0, consol.child.increasedbounds1, consol.child.increasedbounds2)

colnames(consol.child) <- c("Serotype", "overall.invasiveness", "overall.invasiveness.low", "overall.invasiveness.high")

consol.child$Serogroup <- NA
consol.child$Serogroup[which(consol.child$Serotype %in% c("4", "6B", "9V", "14", "18C", "19F", "23F"))] <- "VT7"
consol.child$Serogroup[which(consol.child$Serotype %in% c("1", "5", "7F"))] <- "VT10"
consol.child$Serogroup[which(consol.child$Serotype %in% c("3", "6A", "19A"))] <- "VT13"
consol.child$Serogroup[is.na(consol.child$Serogroup)] <- "NVT"

plot.consol.inv(consol.child) # pdf dim 5 x 12
write.csv(consol.child, file = "consol.child-new.csv")

# kendall's tau correlation with colijn OR
or.df<-read.table(file="summarised.invasiveOR_second_revision.tab",header=TRUE)
combined.df<-
  or.df %>%
  dplyr::select(Serotype,InfantOR) %>%
  dplyr::inner_join(consol.child %>% dplyr::select(Serotype,overall.invasiveness) %>% dplyr::distinct())
tau <- cor.test(combined.df$overall.invasiveness, combined.df$InfantOR,  method="kendall") # 0.61, p-value = p-value = 5.982e-08

# Consolidated invasiveness for adults in pre-PCV time periods
inv.adult.consol <- adult.abs.ds %>% #data.frame(map(adult.abs.ds, unlist)) %>% 
  filter(!DS %in% c("Navajo.post.PCV7")) %>% group_by(Serotype) %>% nest()
#eachglobalinvDS_adu <- inv.adult.consol %>% select(Serotype, DS)
#write.csv(eachglobalinvDS_adu, 'globalinv_eachDS_adults.csv')
nDS.adult <- sapply(inv.adult.consol$data, function(x) nrow(x)) # number of datasets for each serotype
excl.sero.adult <- inv.adult.consol$Serotype[which(nDS.adult == 1)] # serotypes with only 1 dataset to be excluded for analysis
inv.adult.consol <- data.frame(inv.adult.consol %>% unnest(cols = c(data)) %>% filter(!Serotype %in% excl.sero.adult))
unique.adult.sero <- unique(inv.adult.consol$Serotype)

consol.adult2 <- lapply(unique.adult.sero, function(x) {calc_invas(inv.adult.consol %>% filter(Serotype == x),
                                                                                               invas = seq(0.00000001, 0.001, length.out = 2000))})
consol.adult <- do.call("rbind", consol.adult2)

increased.bounds.serotypes <- c("1","20", "4", "8", "24F", "17F", "35A")
consol.adult.increasedbounds <- lapply(increased.bounds.serotypes, function(x) {calc_invas(inv.adult.consol %>% filter(Serotype == x),
                                                                                           invas = seq(0.00000001, 0.02, length.out = 2000))})
consol.adult.increasedbounds <- do.call("rbind", consol.adult.increasedbounds)

consol.adult <- consol.adult[!(consol.adult$Serotype %in% increased.bounds.serotypes),]
consol.adult <- bind_rows(consol.adult, consol.adult.increasedbounds)

colnames(consol.adult) <- c("Serotype", "overall.invasiveness", "overall.invasiveness.low", "overall.invasiveness.high")

consol.adult$Serogroup <- NA
consol.adult$Serogroup[which(consol.adult$Serotype %in% c("4", "6B", "9V", "14", "18C", "19F", "23F"))] <- "VT7"
consol.adult$Serogroup[which(consol.adult$Serotype %in% c("1", "5", "7F"))] <- "VT10"
consol.adult$Serogroup[which(consol.adult$Serotype %in% c("3", "6A", "19A"))] <- "VT13"
consol.adult$Serogroup[is.na(consol.adult$Serogroup)] <- "NVT"

plot.consol.inv(consol.adult) # pdf dim 5 x 12
write.csv(consol.adult, file = "consol.adult-new.csv")

# Consolidated invasiveness for both age groups in one graph
consol.child$agegrp <- rep("children", nrow(consol.child))
consol.adult$agegrp <- rep("adults", nrow(consol.adult))
inv.tot.consol <- rbind(consol.child, consol.adult)
plot.consol.inv(inv.tot.consol) # pdf dim 5 x 12.5

# Max and min and ranges
consol.adult$Serotype[which.min(consol.adult$overall.invasiveness)]
min(consol.adult$overall.invasiveness)
consol.adult$Serotype[which.max(consol.adult$overall.invasiveness)]
max(consol.adult$overall.invasiveness)
ggplot(consol.adult) + geom_freqpoly(aes(overall.invasiveness), bins = 150) + coord_cartesian(xlim=c(0,0.0005))

consol.child$Serotype[which.min(consol.child$overall.invasiveness)]
min(consol.child$overall.invasiveness)
consol.child$Serotype[which.max(consol.child$overall.invasiveness)]
max(consol.child$overall.invasiveness)
ggplot(consol.child) + geom_freqpoly(aes(overall.invasiveness), bins = 150) + coord_cartesian(xlim=c(0,0.0005))

consol.child.post$Serotype[which.min(consol.child.post$overall.invasiveness)]
min(consol.child.post$overall.invasiveness)
consol.child.post$Serotype[which.max(consol.child.post$overall.invasiveness)]
max(consol.child.post$overall.invasiveness)
ggplot(consol.child.post) + geom_freqpoly(aes(overall.invasiveness), bins = 150) + coord_cartesian(xlim=c(0,0.0005))


#### Post-vaccination consolidated ----- SUPPLEMENTARY MATERIAL
inv.child.consol.post <- data.frame(map(child.abs.ds, unlist)) %>% 
  filter(DS %in% c("France.post.PCV7", "Barcelona.post.PCV7", "France.post.PCV13", "Bogota.post.PCV7", "Massachusetts.post.PCV7","Navajo.post.PCV7", 
                    "Atlanta.post.PCV7", 'Netherlands.post.PCV7')) %>% group_by(Serotype) %>% nest()
nDS.post <- sapply(inv.child.consol.post$data, function(x) nrow(x)) # number of datasets for each serotype
excl.sero.post <- inv.child.consol.post$Serotype[which(nDS.post == 1)] # serotypes with only 1 dataset to be excluded for analysis
inv.child.consol.post <- data.frame(inv.child.consol.post %>% unnest(cols = c(data)) %>% filter(!Serotype %in% excl.sero.post))
unique.child.sero.post <- unique(inv.child.consol.post$Serotype)

consol.child2.post <- lapply(unique.child.sero.post, function(x) {calc_invas(inv.child.consol.post %>% filter(Serotype == x),
                                                                             invas = seq(0.00000001, 0.001, length.out = 2000))})
consol.child.post <- do.call("rbind", consol.child2.post)

increased.bounds.serotypes.childpost1 <- c("5", "7F", "12F", "25A", "12B")
consol.child.post.increasedbounds1 <- lapply(increased.bounds.serotypes.childpost1, function(x) {calc_invas(inv.child.consol.post %>% filter(Serotype == x),
                                                                                                            invas = seq(0.00000001, 0.02, length.out = 2000))})
consol.child.post.increasedbounds1 <- do.call("rbind", consol.child.post.increasedbounds1)

Sero1.childpost <- calc_invas(inv.child.consol.post %>% filter(Serotype == "1"),
                              invas = seq(0.00000001, 0.04, length.out = 2000))

increased.bounds.serotypes.childpost2 <- c("4", "8", "10B", "33F", "38", "24F", "18C", "18A", "3", "14")
consol.child.post.increasedbounds2 <- lapply(increased.bounds.serotypes.childpost2, function(x) {calc_invas(inv.child.consol.post %>% filter(Serotype == x),
                                                                                                            invas = seq(0.00000001, 0.002, length.out = 2000))})
consol.child.post.increasedbounds2 <- do.call("rbind", consol.child.post.increasedbounds2)

increased.bounds.serotypes.childpost3 <- c("33A", "37")
consol.child.post.increasedbounds3 <- lapply(increased.bounds.serotypes.childpost3, function(x) {calc_invas(inv.child.consol.post %>% filter(Serotype == x),
                                                                                                            invas = seq(0.00000001, 0.1, length.out = 2000))})
consol.child.post.increasedbounds3 <- do.call("rbind", consol.child.post.increasedbounds3)

consol.child.post <- consol.child.post[!(consol.child.post$Serotype %in% c(increased.bounds.serotypes.childpost1, "1", increased.bounds.serotypes.childpost2, increased.bounds.serotypes.childpost3)),]
consol.child.post <- bind_rows(consol.child.post, Sero1.childpost, consol.child.post.increasedbounds1, consol.child.post.increasedbounds2, consol.child.post.increasedbounds3)

colnames(consol.child.post) <- c("Serotype", "overall.invasiveness", "overall.invasiveness.low", "overall.invasiveness.high")#, 'Serogroup')

consol.child.post$Serogroup <- NA
consol.child.post$Serogroup[which(consol.child.post$Serotype %in% c("4", "6B", "9V", "14", "18C", "19F", "23F"))] <- "VT7"
consol.child.post$Serogroup[which(consol.child.post$Serotype %in% c("1", "5", "7F"))] <- "VT10"
consol.child.post$Serogroup[which(consol.child.post$Serotype %in% c("3", "6A", "19A"))] <- "VT13"
consol.child.post$Serogroup[is.na(consol.child.post$Serogroup)] <- "NVT"

plot.consol.inv(consol.child.post) # pdf dim 5 x 12
write.csv(consol.child.post, "consol.child.postPCV.csv")

#### Weinberger's prior vs our likelihood and posterior  ----------------------------------------------------------------------------------------

plot.serolikelihood <- function(df, serotype, agegroup) { # function that plots the likelihood of dataset-specific and Weinberger's posterior 
  
  serotype.df <- df %>% filter(Serotype == serotype)
  dfs <- unique(serotype.df$DS)
  eachdf1 <- lapply(dfs, function(x) get(paste(x, "distrib", agegroup, serotype, sep = ".")))
  eachdf2 <- bind_rows(eachdf1)
  eachdf.new <- list()
  
  for (i in 1:length(eachdf1)) {
    newbinsdev <- seq(0, max(eachdf1[[i]]$bins), length.out = 500)
    inv.cut <- cut(eachdf1[[i]]$bins, newbinsdev)
    dist.sum <- tapply(eachdf1[[i]]$distrib, INDEX = inv.cut, FUN = sum)
    dist.mean <- dist.sum/100000
    dist.mean[is.na(dist.mean)] <- 0
    eachdf.new[[i]] <- data.frame(invas_dev = newbinsdev[-length(newbinsdev)], distrib = dist.mean, DS = dfs[i])
  }
  
  eachdf <- bind_rows(eachdf.new)
  
  sero.invas.dev <- rlnorm(n = 560000, meanlog = inv_priors_df$log.inv.age1[which(inv_priors_df$st == serotype)]-log(1000), 
                           sdlog = inv_priors_df$stddev[which(inv_priors_df$st == serotype)])
  sero.prior <- dlnorm(sero.invas.dev, meanlog = inv_priors_df$log.inv.age1[which(inv_priors_df$st == serotype)]-log(1000), 
                         sdlog = inv_priors_df$stddev[which(inv_priors_df$st == serotype)])
  sero.priordf <- data.frame(invas_dev = sero.invas.dev, distrib = sero.prior*max(eachdf$distrib)/max(sero.prior), DS = "prior")
  eachdf <- rbind(eachdf, sero.priordf)
  eachdf$col <- "Our posterior"
  eachdf$col[which(eachdf$DS == "prior")] <- "Weinberger et al"
  levels(eachdf$DS) <- c(levels(eachdf$DS)[1:(length(levels(eachdf$DS))-1)],"Weinberger et al")
  levels(eachdf$DS) <- gsub("[.]"," ", levels(eachdf$DS))
  levels(eachdf$DS) <- gsub("pre PCV","pre-PCV", levels(eachdf$DS))
  levels(eachdf$DS) <- gsub("post PCV","post-PCV", levels(eachdf$DS))
  levels(eachdf$DS) <- gsub("E W ","E&W ", levels(eachdf$DS))
  
  ggplot(eachdf) + 
    geom_line(aes(x = invas_dev, y = distrib, colour = DS, alpha = factor(col))) + scale_alpha_discrete(range = c(0.3, 1)) +
    scale_y_continuous(labels = scales::scientific) + ggtitle(paste("Serotype", serotype, df$agegrp, "invasiveness", sep = " ")) + 
    labs(x = "Invasiveness", y = "Probability Density") +
    theme_minimal() + theme(legend.position = "none")
}

#
####### children specific invasiveness#

serotypes.with.priors <- unique(child.abs.ds$Serotype)[unique(child.abs.ds$Serotype) %in% inv_priors_df$st] # serotypes that are included in Weinberger's analysis
eachsero <- lapply(serotypes.with.priors, function(x) plot.serolikelihood(df = abs.child, serotype = x, agegroup = "children"))
likelihoodvsprior1 <- grid.arrange(eachsero[[1]], eachsero[[2]], eachsero[[3]], #eachsero[[4]], eachsero[[5]], eachsero[[6]], eachsero[[7]], eachsero[[8]],
                                  #eachsero[[9]], eachsero[[10]], eachsero[[11]], eachsero[[12]], eachsero[[13]], eachsero[[14]], eachsero[[15]], eachsero[[16]],
                                  #
                                  # eachsero[[17]], eachsero[[18]], eachsero[[19]], eachsero[[20]], eachsero[[21]], eachsero[[22]], eachsero[[23]], eachsero[[24]],
                                  # eachsero[[25]], eachsero[[26]], eachsero[[27]], eachsero[[28]], eachsero[[29]], eachsero[[30]], eachsero[[31]], eachsero[[32]],
                                  # 
                                  # eachsero[[33]], eachsero[[34]], eachsero[[35]], eachsero[[36]], eachsero[[37]], eachsero[[38]], eachsero[[39]], eachsero[[40]],
                                  # eachsero[[41]], eachsero[[42]], eachsero[[43]], eachsero[[44]], eachsero[[45]], eachsero[[46]], eachsero[[47]], eachsero[[48]],

                                  # eachsero[[49]], eachsero[[50]], eachsero[[51]], eachsero[[52]], eachsero[[53]], eachsero[[54]], eachsero[[55]], eachsero[[56]],
                                  # eachsero[[57]], eachsero[[58]], eachsero[[59]],
                                  nrow = 1) #pdf dim 8 x 18


## adults specific invasiveness##

serotypes.with.priors.adu <- unique(adult.abs.ds$Serotype)[unique(adult.abs.ds$Serotype) %in% inv_priors_df$st]
eachsero.adu <- lapply(serotypes.with.priors.adu, function(x) plot.serolikelihood(df = abs.adult, serotype = x, agegroup = "adults"))
likelihoodvsprior1.adu <- grid.arrange(eachsero.adu[[1]], eachsero.adu[[2]], eachsero.adu[[3]], eachsero.adu[[4]], eachsero.adu[[5]], eachsero.adu[[6]], 
                                       eachsero.adu[[7]], eachsero.adu[[8]], eachsero.adu[[9]], eachsero.adu[[10]], eachsero.adu[[11]], eachsero.adu[[12]],
                                       eachsero.adu[[13]], eachsero.adu[[14]], eachsero.adu[[15]],
                                       #
                                       # eachsero.adu[[16]], eachsero.adu[[17]], eachsero.adu[[18]], eachsero.adu[[19]], eachsero.adu[[20]], eachsero.adu[[21]],
                                       # eachsero.adu[[22]], eachsero.adu[[23]], eachsero.adu[[24]], eachsero.adu[[25]], eachsero.adu[[26]], eachsero.adu[[27]],
                                       # eachsero.adu[[28]],

                                       # eachsero.adu[[29]], eachsero.adu[[30]], eachsero.adu[[31]], eachsero.adu[[32]], eachsero.adu[[33]], eachsero.adu[[34]],
                                       # eachsero.adu[[35]], eachsero.adu[[36]], eachsero.adu[[37]], eachsero.adu[[38]], eachsero.adu[[39]], eachsero.adu[[40]],
                                       # 
                                       # eachsero.adu[[41]], eachsero.adu[[42]], eachsero.adu[[43]], eachsero.adu[[44]], eachsero.adu[[45]], eachsero.adu[[46]],
                                       # eachsero.adu[[47]], eachsero.adu[[48]], eachsero.adu[[49]],
                                       nrow = 3) #pdf dim 8 x 18


## children consolidated invasiveness ##
consol.child.likelihoodpri <- lapply(unique.child.sero, function(x) {calc_invas(inv.child.consol %>% filter(Serotype == x))}) # plots; uncomment in calc_invas fn

legendgrob <- g_legend(calc_invas(inv.child.consol %>% filter(Serotype == "1"))+ theme(legend.position = "right"))

likelihoodvsprior1.child.consol <- grid.arrange(#consol.child.likelihoodpri[[1]], consol.child.likelihoodpri[[2]], consol.child.likelihoodpri[[3]], 
                                                #consol.child.likelihoodpri[[4]], consol.child.likelihoodpri[[5]], consol.child.likelihoodpri[[6]], 
                                                #consol.child.likelihoodpri[[7]], consol.child.likelihoodpri[[8]], consol.child.likelihoodpri[[9]], 
                                                #consol.child.likelihoodpri[[10]], consol.child.likelihoodpri[[11]], consol.child.likelihoodpri[[12]], 
                                                #consol.child.likelihoodpri[[13]], 
                                                # consol.child.likelihoodpri[[14]], consol.child.likelihoodpri[[15]], consol.child.likelihoodpri[[16]],
                                                # consol.child.likelihoodpri[[17]], consol.child.likelihoodpri[[18]], consol.child.likelihoodpri[[19]],
                                                # consol.child.likelihoodpri[[20]], consol.child.likelihoodpri[[21]], consol.child.likelihoodpri[[22]],
                                                # consol.child.likelihoodpri[[23]], consol.child.likelihoodpri[[24]], consol.child.likelihoodpri[[25]],
                                                # consol.child.likelihoodpri[[26]],
                                                consol.child.likelihoodpri[[27]], consol.child.likelihoodpri[[28]], consol.child.likelihoodpri[[29]],
                                                consol.child.likelihoodpri[[30]], consol.child.likelihoodpri[[31]], consol.child.likelihoodpri[[32]],
                                                consol.child.likelihoodpri[[33]], consol.child.likelihoodpri[[34]], consol.child.likelihoodpri[[35]],
                                                consol.child.likelihoodpri[[36]], consol.child.likelihoodpri[[37]], consol.child.likelihoodpri[[38]],
                                                consol.child.likelihoodpri[[39]],
                                                legendgrob, nrow = 4) #pdf dim 7 x 10


####### Supplementary material for paper
# change calc_invas to return(plot)
sero1consolidatedchildren <- calc_invas(inv.child.consol %>% filter(Serotype == "1"), 
                                        invas = seq(0.00000001, 0.1, length.out = 2000)) + coord_cartesian(xlim = c(0,0.03)) + labs(colour = "")
sero7Fconsolidatedchildren <- calc_invas(inv.child.consol %>% filter(Serotype == "7F"), 
                                         invas = seq(0.00000001, 0.02, length.out = 2000)) + coord_cartesian(xlim = c(0,0.005))
sero12Fconsolidatedchildren <- calc_invas(inv.child.consol %>% filter(Serotype == "12F"), 
                                          invas = seq(0.00000001, 0.02, length.out = 2000)) + coord_cartesian(xlim = c(0,0.004))
#legendgrob <- g_legend(sero1consolidatedchildren + theme(legend.position = "bottom")+ labs(colour = ""))
# consolidated1.7F.12F <- grid.arrange(sero1consolidatedchildren, sero7Fconsolidatedchildren, sero12Fconsolidatedchildren, 
#              bottom = legendgrob, nrow = 1)
consolidated1.7F.12F <- ggarrange(sero1consolidatedchildren, sero7Fconsolidatedchildren, sero12Fconsolidatedchildren, 
                                     common.legend = TRUE, legend = 'bottom', nrow = 1)

specific1.7F.12F <- grid.arrange(eachsero[[1]]+theme_classic()+ theme(legend.position = "right") + 
                                  coord_cartesian(xlim = c(0,0.1))+ guides(alpha = FALSE) + labs(colour = "Dataset") +
                                   ggtitle("Serotype 1"), 
                                 eachsero[[30]]+theme_classic()+ theme(legend.position = "right") + #eachsero[[31]] children
                                   coord_cartesian(xlim = c(0,0.1))+guides(alpha = FALSE) + labs(colour = "Dataset") +
                                   ggtitle("Serotype 7F"), 
                                 eachsero[[38]]+theme_classic() + theme(legend.position = "right") + #eachsero[[39]] children
                                   coord_cartesian(xlim = c(0,0.02))+guides(alpha = FALSE) + labs(colour = "Dataset") +
                                   ggtitle("Serotype 12F"),
                                 nrow = 1) # 6 x 12

supplementarymaterial <- grid.arrange(specific1.7F.12F, consolidated1.7F.12F, nrow = 2) # 8 x 15

## adults consolidated invasiveness ##
consol.adult.likelihoodpri <- lapply(unique.adult.sero, function(x) {calc_invas(inv.adult.consol %>% filter(Serotype == x))}) # plots

legendgrob <- g_legend(calc_invas(inv.adult.consol %>% filter(Serotype == "1"))+ theme(legend.position = "right"))

likelihoodvsprior1.adult.consol <- grid.arrange(consol.adult.likelihoodpri[[1]], consol.adult.likelihoodpri[[2]], consol.adult.likelihoodpri[[3]],
                                                consol.adult.likelihoodpri[[4]], consol.adult.likelihoodpri[[5]], consol.adult.likelihoodpri[[6]],
                                                consol.adult.likelihoodpri[[7]], consol.adult.likelihoodpri[[8]], consol.adult.likelihoodpri[[9]],
                                                consol.adult.likelihoodpri[[10]], consol.adult.likelihoodpri[[11]],
                                                # 
                                                # consol.adult.likelihoodpri[[12]], consol.adult.likelihoodpri[[13]], consol.adult.likelihoodpri[[14]],
                                                # consol.adult.likelihoodpri[[15]], consol.adult.likelihoodpri[[16]], consol.adult.likelihoodpri[[17]],
                                                # consol.adult.likelihoodpri[[18]], consol.adult.likelihoodpri[[19]], consol.adult.likelihoodpri[[20]],
                                                # consol.adult.likelihoodpri[[21]], consol.adult.likelihoodpri[[22]],
                                                # 
                                                # consol.adult.likelihoodpri[[23]], consol.adult.likelihoodpri[[24]], consol.adult.likelihoodpri[[25]],
                                                # consol.adult.likelihoodpri[[26]], consol.adult.likelihoodpri[[27]], consol.adult.likelihoodpri[[28]],
                                                # consol.adult.likelihoodpri[[29]], consol.adult.likelihoodpri[[30]],consol.adult.likelihoodpri[[31]],
                                                # consol.adult.likelihoodpri[[32]], consol.adult.likelihoodpri[[33]],
                                                
                                                legendgrob, nrow = 4) #pdf dim 7 x 10

# Supplementary material adults
sero1consolidatedadults <- calc_invas(inv.adult.consol %>% filter(Serotype == "1"), 
                                      invas = seq(0.00000001, 0.001, length.out = 2000)) + coord_cartesian(xlim = c(0,0.03)) + labs(colour = "")
sero7Fconsolidatedadults <- calc_invas(inv.adult.consol %>% filter(Serotype == "7F"),
                                       invas = seq(0.00000001, 0.001, length.out = 2000)) + coord_cartesian(xlim = c(0,0.005))
sero12Fconsolidatedadults <- calc_invas(inv.adult.consol %>% filter(Serotype == "12F"),
                                        invas = seq(0.00000001, 0.001, length.out = 2000)) + coord_cartesian(xlim = c(0,0.004))
sero3consolidatedadults <- calc_invas(inv.adult.consol %>% filter(Serotype == "3"),
                                      invas = seq(0.00000001, 0.001, length.out = 2000)) + coord_cartesian(xlim = c(0,0.004))

legendgrob <- g_legend(sero1consolidatedadults + theme(legend.position = "right")+ labs(colour = ""))
consolidated1.7F.12F.adu <- grid.arrange(sero1consolidatedadults, sero7Fconsolidatedadults, sero12Fconsolidatedadults, 
                                     right = legendgrob, nrow = 1)

specific1.7F.12F <- grid.arrange(eachsero.adu[[1]]+theme_classic()+ theme(legend.position = "right") + 
                                   coord_cartesian(xlim = c(0,0.5))+ guides(alpha = FALSE) + labs(colour = "Dataset") +
                                   ggtitle("Serotype 1"), 
                                 eachsero.adu[[30]]+theme_classic()+ theme(legend.position = "right") + #eachsero[[31]] children
                                   coord_cartesian(xlim = c(0,0.5))+guides(alpha = FALSE) + labs(colour = "Dataset") +
                                   ggtitle("Serotype 7F"), 
                                 eachsero.adu[[38]]+theme_classic() + theme(legend.position = "right") + #eachsero[[39]] children
                                   coord_cartesian(xlim = c(0,0.2))+guides(alpha = FALSE) + labs(colour = "Dataset") +
                                   ggtitle("Serotype 12F"),
                                 nrow = 1) # 6 x 12

supplementarymaterial <- grid.arrange(specific1.7F.12F, consolidated1.7F.12F, nrow = 2) # 8 x 15

#### Comparison study-specific invasiveness vs global -------------------------------------------------------------------------------------------------------

unique.DS.comparison.child <- unique(inv.child.consol$DS)
prePCV.DS.vs.overall.child <- lapply(unique.DS.comparison.child, function(x) 
  {plot.DS.overall(sero = abs.child %>% filter(DS == toString(x)), consol.child)})


ggarrange(prePCV.DS.vs.overall.child[[1]], prePCV.DS.vs.overall.child[[2]], prePCV.DS.vs.overall.child[[3]],
          prePCV.DS.vs.overall.child[[4]], prePCV.DS.vs.overall.child[[5]], prePCV.DS.vs.overall.child[[6]],
          prePCV.DS.vs.overall.child[[7]], prePCV.DS.vs.overall.child[[8]] + ylab('E&W pre-PCV'),
          prePCV.DS.vs.overall.child[[9]], prePCV.DS.vs.overall.child[[10]], prePCV.DS.vs.overall.child[[11]],
          nrow = 4, ncol = 3, common.legend = TRUE, legend = "bottom") #pdf 9 x 11

# pdf pagess
# mypath.childcomparison.inv <- file.path("Q:","Technical","R","Case-to-carrier","Figures for paper", "Supplementary Fig 4 - Overall vs DS-specific inv", 
#                                         paste("Localvsglobalinv_apr21", ".pdf", sep=""))
# pdf(mypath.childcomparison.inv, width = 7, height = 5)
# print(prePCV.DS.vs.overall.child)
# dev.off()

plot.posteriors <- function(df, serotype, agegroup) { # this fn requires the invasiveness fn to be run and global inv fn
  serotype.df <- df %>% filter(Serotype == serotype)
  dfs <- unique(serotype.df$DS) # unique DS
  eachdf1 <- lapply(dfs, function(x) get(paste(x, "distrib", agegroup, serotype, sep = ".")))
  eachdf2 <- do.call(rbind, eachdf1) #bind_rows(eachdf1)
  eachdf.max <- max(eachdf2$distrib)
  eachdf.norm <- lapply(eachdf1, function(x) x$distrib*eachdf.max/max(x$distrib))
  
  for (i in 1:length(eachdf1)) {
    eachdf1[[i]]$distrib2 <- eachdf.norm[[i]]
    eachdf1[[i]]$DS <- dfs[i]
  }
  
  eachdf <- bind_rows(eachdf2)
  eachdf$DS <- str_replace_all(str_replace_all(eachdf$DS, "[.]", " "), " P", "-P")
  consol.distrib <- get(paste(paste("sero", serotype, sep = ""), "distrib", agegroup, sep = "."))
  consol.distrib$DS <- "Global inv"
  consol.distrib$distrib <- consol.distrib$distrib*eachdf.max/max(consol.distrib$distrib)
  eachdf <- rbind(eachdf, consol.distrib)
  eachdf$col <- "Local"
  eachdf$col[which(eachdf$DS == "Global inv")] <- "Global"
  
  eachdf$DS <- fct_relevel(eachdf$DS, "Global inv", after = Inf)

  ggplot() + 
    geom_line(eachdf, mapping = aes(x = bins, y = distrib, colour = DS, alpha = factor(col))) + scale_alpha_discrete(range = c(1, 0.3)) +
    scale_y_continuous(labels = scales::scientific) + coord_cartesian(xlim = c(0,max(consol.distrib$bins)+0.01)) + 
    ggtitle(paste("Serotype", serotype, sep = " ")) + 
    labs(x = "Invasiveness (cases per carrier yr)", y = "Posterior PDF", alpha = "Model", colour = "Dataset") +
    xlab(expression(Invasiveness~(case~carrier^{-1}~year^{-1}))) +
    guides(colour = guide_legend(order = 1), alpha = guide_legend(order = 2))+
    theme_minimal() 
}

presentation_serotypes <- c('1', '3', '4', '5', '7F', '14')
consolvssingle.plots <- lapply(presentation_serotypes, function(x) {plot.posteriors(df = inv.child.consol, serotype = x, agegroup = "children")}) # unique.child.sero
consol.vs.single <- grid.arrange(consolvssingle.plots[[1]], consolvssingle.plots[[2]], 
  consolvssingle.plots[[3]], 
  consolvssingle.plots[[4]], 
  consolvssingle.plots[[5]],consolvssingle.plots[[6]], 
  #consolvssingle.plots[[7]], consolvssingle.plots[[8]], consolvssingle.plots[[9]], consolvssingle.plots[[10]],
                                 
                                 # consolvssingle.plots[[11]], consolvssingle.plots[[12]], consolvssingle.plots[[13]], consolvssingle.plots[[14]], consolvssingle.plots[[15]],
                                 # consolvssingle.plots[[16]], consolvssingle.plots[[17]], consolvssingle.plots[[18]], consolvssingle.plots[[19]], consolvssingle.plots[[20]],

                                 # consolvssingle.plots[[21]], consolvssingle.plots[[22]], consolvssingle.plots[[23]], consolvssingle.plots[[24]], consolvssingle.plots[[25]],
                                 # consolvssingle.plots[[26]], consolvssingle.plots[[27]], consolvssingle.plots[[28]], consolvssingle.plots[[29]], consolvssingle.plots[[30]],
                                 # 
                                 #consolvssingle.plots[[31]], consolvssingle.plots[[32]], consolvssingle.plots[[33]], consolvssingle.plots[[34]], consolvssingle.plots[[35]],
                                 #consolvssingle.plots[[36]], consolvssingle.plots[[37]], consolvssingle.plots[[38]], consolvssingle.plots[[39]],
                                 
                                 nrow = 2) # pdf dim 5 x 12 #nrow = 4 for all figures
# 3 figures per row = pdf 9 x 15

consolvssingle.plots.adu <- lapply(unique.adult.sero, function(x) {plot.posteriors(df = inv.adult.consol, serotype = x, agegroup = "adults")})

consol.vs.single.adu <- grid.arrange(#consolvssingle.plots.adu[[1]], consolvssingle.plots.adu[[2]], consolvssingle.plots.adu[[3]], consolvssingle.plots.adu[[4]], 
                                     #consolvssingle.plots.adu[[5]], consolvssingle.plots.adu[[6]], consolvssingle.plots.adu[[7]], consolvssingle.plots.adu[[8]],
                                     #consolvssingle.plots.adu[[9]], consolvssingle.plots.adu[[10]], consolvssingle.plots.adu[[11]],
                                 
                                 # consolvssingle.plots.adu[[12]], consolvssingle.plots.adu[[13]], consolvssingle.plots.adu[[14]], consolvssingle.plots.adu[[15]],
                                 # consolvssingle.plots.adu[[16]], consolvssingle.plots.adu[[17]], consolvssingle.plots.adu[[18]], consolvssingle.plots.adu[[19]],
                                 # consolvssingle.plots.adu[[20]], consolvssingle.plots.adu[[21]], consolvssingle.plots.adu[[22]],
                                 
                                 consolvssingle.plots[[23]], consolvssingle.plots[[24]], consolvssingle.plots[[25]], consolvssingle.plots[[26]], consolvssingle.plots[[27]],
                                 consolvssingle.plots[[28]], consolvssingle.plots[[29]], consolvssingle.plots[[30]], consolvssingle.plots[[31]], consolvssingle.plots[[32]],
                                 consolvssingle.plots[[33]],
                                 
                                 nrow = 4) # pdf dim 7 x 11


#### Comparison pre vs post vaccination --------------------------------------------------------------------------------------------------------------------

Atlantaprevspost <- abs.tot[grepl('^Atlanta', abs.tot$DS),]
Atlantaprevspost$DS <- factor(Atlantaprevspost$DS, levels = c("Atlanta.pre.PCV", "Atlanta.post.PCV7"),
                              labels = c("Pre-PCV", "Post-PCV7"))
Atlantaplot <- ggplot(Atlantaprevspost) + 
  geom_point(aes(x = Serotype, y = invasiveness, colour = DS), position = position_dodge(width = 0.5)) + 
  geom_errorbar(aes(x = Serotype, ymin = invasiveness.low, ymax = invasiveness.high, width = 0.01, colour = DS), position = position_dodge(width = 0.5)) + 
  scale_y_continuous(trans = 'log10') + 
 theme_minimal() + theme_bw() + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),#, size = 14),
                                      #axis.text.y = element_text(size = 14),
                                      #legend.text=element_text(size=14),
                                      text = element_text(size = 15)) +   
  scale_color_manual(values = c("Pre-PCV" = '#F8766D', "Post-PCV7" = '#00BFC4')) +
  labs(y = "", title = "Atlanta", colour = "Dataset")

Franceprevspost <- abs.tot[grepl('^France', abs.tot$DS),]
Franceprevspost$DS <- factor(Franceprevspost$DS, levels = c("France.post.PCV7", "France.post.PCV13"),
                             labels = c("Post-PCV7", "Post-PCV13"))
Franceplot <- ggplot(Franceprevspost) + 
  geom_point(aes(x = Serotype, y = invasiveness, colour = DS), position = position_dodge(width = 0.5)) + 
  geom_errorbar(aes(x = Serotype, ymin = invasiveness.low, ymax = invasiveness.high, width = 0.01, colour = DS), position = position_dodge(width = 0.5)) + 
  scale_y_continuous(trans = 'log10') + 
  theme_minimal() + theme_bw() + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),#, size = 14),
                                       #axis.text.y = element_text(size = 14),
                                       #legend.text=element_text(size=14),
                                       text = element_text(size = 15)) +  
  scale_color_manual(values = c("Post-PCV7" = '#00BFC4', "Post-PCV13" = '#7CAE00')) +
  labs(y = "", title = "France", colour = "Dataset")

Bogotaprevspost <- abs.tot[grepl('^Bogota', abs.tot$DS),]
Bogotaprevspost$DS <- factor(Bogotaprevspost$DS, levels = c("Bogota.pre.PCV", "Bogota.post.PCV7"),
                             labels = c("Pre-PCV", "Post-PCV7"))
Bogotaplot <- ggplot(Bogotaprevspost) + 
  geom_point(aes(x = Serotype, y = invasiveness, colour = DS), position = position_dodge(width = 0.5)) + 
  geom_errorbar(aes(x = Serotype, ymin = invasiveness.low, ymax = invasiveness.high, width = 0.01, colour = DS), position = position_dodge(width = 0.5)) + 
  scale_y_continuous(trans = 'log10') + 
   theme_minimal() + theme_bw() + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),#, size = 14),
                                        #axis.text.y = element_text(size = 14),
                                        #legend.text=element_text(size=14),
                                        text = element_text(size = 15)) +  
  scale_color_manual(values = c("Pre-PCV" = '#F8766D', "Post-PCV7" = '#00BFC4')) +
  labs(title = "Bogota", colour = "Dataset", y = "")# +
  #ylab(expression(Invasiveness~(case~carrier^{-1}~year^{-1}))) +

NLprevspost <- abs.tot[grepl('^Netherlands', abs.tot$DS),]
NLprevspost$DS <- factor(NLprevspost$DS, levels = c("Netherlands.pre.PCV", "Netherlands.post.PCV7", "Netherlands.post.PCV10"),
                             labels = c("Pre-PCV", "Post-PCV7", "Post-PCV10"))
NLplot <- ggplot(NLprevspost) + 
  geom_point(aes(x = Serotype, y = invasiveness, colour = DS), position = position_dodge(width = 0.5)) + 
  geom_errorbar(aes(x = Serotype, ymin = invasiveness.low, ymax = invasiveness.high, width = 0.01, colour = DS), position = position_dodge(width = 0.5)) + 
  scale_y_continuous(trans = 'log10') + 
  theme_minimal() + theme_bw() + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),#, size = 14),
                                       #axis.text.y = element_text(size = 14),
                                       #legend.text=element_text(size=14),
                                       text = element_text(size = 15)) +  
  scale_color_manual(values = c("Pre-PCV" = '#F8766D', "Post-PCV7" = '#00BFC4', "Post-PCV10" = '#7CAE00'))+ #'#C77CFF')) +
  labs(title = "Netherlands", colour = "Dataset", y = "")# +
#ylab(expression(Invasiveness~(case~carrier^{-1}~year^{-1}))) +


pdf("Q:/Technical/R/Case-to-carrier/Figures for paper/Fig 5 - Pre vs Post PCV/prevspost_apr21.pdf", width = 12, height = 8)
grid.arrange(Atlantaplot, Bogotaplot, Franceplot, NLplot, nrow = 2,
             left = textGrob(expression(Invasiveness~(case~carrier^{-1}~year^{-1})), rot = 90, vjust = 1,
                             gp = gpar(cex = 1.3))) # pdf dim 7 x 16
dev.off()


#### Comparison children vs adults in each setting ---------------------------------------------------------------------------------------------------------

uniqueDS <- unique(abs.tot$DS)

childrenvsadults_perDS <- function(df) {
  df$agegrp <- factor(df$agegrp, levels = c("children", "adults"), labels = c('Children', 'Adults'))
  title <- gsub('.p', ' p', gsub('.PCV','-PCV',unique(df$DS)))
  ggplot(df) + 
    geom_point(aes(x = Serotype, y = invasiveness, colour = agegrp), position = position_dodge(width = 0.5)) + 
    geom_errorbar(aes(x = Serotype, ymin = invasiveness.low, ymax = invasiveness.high, width = 0.1, colour = agegrp), position = position_dodge(width = 0.5)) +
    labs(title = title, x = "Serotype", y = "Invasiveness (case per carrier yr)", colour = "Age group") + theme_minimal() + theme_bw() + 
    scale_y_continuous(trans = 'log10') + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
                                                text = element_text(size = 17))
}

childrenvsadults_perDS_plots <- lapply(uniqueDS, function(x) childrenvsadults_perDS(abs.tot %>% filter(DS == x)))

mypath.childvsadultsperDS.inv <- file.path("Q:","Technical","R","Case-to-carrier", "childrenvsadultsperDS.pdf", sep="")
pdf(mypath.childvsadultsperDS.inv, width = 15, height = 5)
print(childrenvsadults_perDS_plots)
dev.off()

# childrenadultlegend <- g_legend(childrenvsadults_perDS_plots[[1]] + theme(legend.position = "bottom"))
# grid.arrange(childrenvsadults_perDS_plots[[1]]+theme(legend.position = "none"), childrenvsadults_perDS_plots[[13]]+theme(legend.position = "none"),
#              childrenvsadults_perDS_plots[[12]]+theme(legend.position = "none"), childrenvsadults_perDS_plots[[10]]+theme(legend.position = "none"), 
#              bottom = childrenadultlegend, nrow = 2) # pdf 6 x 20

# includes portugal
legendgr <- get_legend(childrenvsadults_perDS_plots[[1]])
childvsad <- ggarrange(childrenvsadults_perDS_plots[[1]]+theme(legend.position = "none")+rremove('ylab'), 
          childrenvsadults_perDS_plots[[16]]+theme(legend.position = "none")+labs(title = "E&W pre-PCV")+rremove('ylab'),
          childrenvsadults_perDS_plots[[10]]+theme(legend.position = "none")+rremove('ylab'), 
          childrenvsadults_perDS_plots[[15]]+theme(legend.position = "none")+rremove('ylab'), 
          nrow = 2, ncol = 2, common.legend = TRUE, legend = 'bottom')
          #legendgrob = get_legend(childrenvsadults_perDS_plots[[1]]))
annotate_figure(childvsad, 
                left = textGrob(expression(Invasiveness~(case~carrier^{-1}~year^{-1})), rot = 90, 
                                 vjust = 1, hjust = -1, gp = gpar(cex = 1.3))) #pdf 7 x 16

### Serotype profiles: serotypes 4, 14, 1, 5 (suppl figure 5) ----------------------------------------------------------------------------------------------

plot.inv.globallocal <- function(sero) {
  VT7 <- c("4", "14")
  VT10 <- c("1","5")
  sero$DS <- sub(".pre.PCV", ".pre-PCV", sero$DS)
  sero$DS <- sub("\\.", " ",sero$DS)
  ggplot(sero, aes(x = DS, y = invasiveness, group = sero$agegrp)) +
    geom_errorbar(aes(ymin = sero$invasiveness.low, ymax = sero$invasiveness.high, group = sero$agegrp), position = position_dodge(0.75),
                  width = 0.01) + 
    geom_point(aes(color = Serogroup, shape = Invtype), 
               size = 2, position = position_dodge(0.75)) +
    labs(x =  "Dataset", 
         y = "Invasiveness", title = paste("Serotype",sero$Serotype, sep = " ")) +
    scale_color_manual(values= cols) +
    scale_y_continuous(trans = 'log10', labels = function(x) format(x, scientific = TRUE)) +
    theme_bw() + theme_light() +
    theme(axis.text.x = element_text(angle = 30, hjust = 1), legend.position = "none") 
}

abs.tot$DS <- sub("pre.PCV", "pre-PCV", abs.tot$DS)
abs.tot$DS <- sub("post.PCV", "post-PCV", abs.tot$DS)
sero4 <- abs.tot %>% filter(Serotype == '4') %>% filter(agegrp == 'children') %>% filter(vaccinationera == 'pre.PCV') %>% mutate(Invtype = "Local invasiveness")
sero14 <-  abs.tot %>% filter(Serotype == '14') %>% filter(agegrp == 'children') %>% filter(vaccinationera == 'pre.PCV') %>% mutate(Invtype = "Local invasiveness")
sero1 <-  abs.tot %>% filter(Serotype == '1') %>% filter(agegrp == 'children') %>% filter(vaccinationera == 'pre.PCV') %>% mutate(Invtype = "Local invasiveness")
sero5 <-  abs.tot %>% filter(Serotype == '5') %>% filter(agegrp == 'children') %>% filter(vaccinationera == 'pre.PCV') %>% mutate(Invtype = "Local invasiveness")


sero2join <- consol.child %>% filter(Serotype %in% c("4", "14", "1", "5"))
sero2join$X <- NULL
colnames(sero2join) <- c("Serotype", "invasiveness", "invasiveness.low", "invasiveness.high", "Serogroup", "agegrp")
sero2join$vaccinationera <- "pre.PCV"
sero2join$Invtype <- "Global invasiveness"
sero2join$DS <- "Global invasiveness"


sero4 <- full_join(sero4, sero2join %>% filter(Serotype == "4"))
sero14 <- full_join(sero14, sero2join %>% filter(Serotype == "14"))
sero1 <- full_join(sero1, sero2join %>% filter(Serotype == "1"))
sero5 <- full_join(sero5, sero2join %>% filter(Serotype == "5"))

sero4plot <- plot.inv.globallocal(sero4) 
sero14plot <- plot.inv.globallocal(sero14)
sero1plot <- plot.inv.globallocal(sero1)
sero5plot <- plot.inv.globallocal(sero5)

grid.arrange(sero4plot, sero14plot, sero1plot, sero5plot, ncol = 2) #pdf 5x9

# is czech DS lower than other DS?
czswenavj <- abs.tot %>% filter(DS %in% c("Czech.pre.PCV", "Sweden.pre.PCV", 
                                          "Navajo.post.PCV", "Oxford.pre.PCV")) %>% 
  filter(agegrp == "children")

ggplot(czswenavj, aes(x = Serotype, y = invasiveness, group = DS)) + 
  geom_errorbar(aes(ymin = invasiveness.low, ymax = invasiveness.high, group = DS), 
                position = position_dodge(0.75), width = 0.01) + 
  geom_point(aes(colour = Serogroup, shape = DS, group = DS), position = position_dodge(0.75), size = 2) + 
  scale_y_continuous(trans = 'log10') + 
  theme_minimal() + theme(axis.text.x = element_text(angle = 45, hjust = 1)) #pdf 5x15

#### Model Comparison: Bayes Factor  -----------------------------------------------------------------------------------------------------------------------

consol_evidence <- function(df) {
  #set.seed(123)
  n <- nrow(df)
  num_dev <- 5000
  
  invas_dev <- runif(n = num_dev, min = 0.0000001, max = 0.5)
  carr_prev_dev <- sapply(seq(1:n), function(s) {rbeta(n = num_dev, shape1 = df[s,5] + 1, shape2 = df[s,7] - df[s,5] + 1)})
  carr_invas <- sapply(seq(1:n), function(x) invas_dev*carr_prev_dev[,x])
  lambda <- sweep(carr_invas, 2, df[,8]*df[,9], FUN = "*")
  LL <- sapply(seq(1:n), function(s) {dpois(df[s, 6], lambda[,s], log = TRUE)})
  totalLL <- rowSums(LL)
  
  return(exp(totalLL))
}

# serotype 4 consolidated evidence
consol_ev_allDS <- consol_evidence(inv.child.consol %>% filter(Serotype == "4"))
consol_ev <- mean(consol_ev_allDS)
consol_EW_4 <- consol_evidence(inv.child.consol %>% filter(Serotype == "4") %>% filter(DS == "E.W.pre.PCV"))
consol_ev_EW4 <- mean(consol_EW_4)

## function to estimate bayes factor
bayesfxr <- function(serotype, df, agegroup) {
  serotype.df <- df %>% filter(Serotype == serotype)
  dfs <- unique(serotype.df$DS)

  eachdf1 <- lapply(dfs, function(x) get(paste(x, "distrib", agegroup, serotype, sep = ".")))
  # distrib.only <- data.frame(lapply(eachdf1, function(x) x$distrib))
  # common.L <- apply(distrib.only, 1, prod)
  # spec.distrib <- data.frame(bins = eachdf1[[1]]$bins, distrib = common.L)

  if (serotype == "15B/C") {serotype <- "15BC"}
  consol.distrib <- get(paste(paste("sero", serotype, sep = ""), "distrib", agegroup, sep = "."))
  consol.distrib$DS <- "consol"
  colnames(consol.distrib) <- colnames(eachdf1[[1]])
  
  AUC.consol <- AUC(consol.distrib$bins, consol.distrib$distrib, method = "trapezoid")
  #AUC.single.each <- sapply(eachdf1, function(x) AUC(x$bins, x$distrib, method = "trapezoid"))
  #AUC.single.all <- prod(AUC.single.each)
  
  #AUC.consol <- sum(consol.distrib$distrib)/nrow(consol.distrib)
  AUC.single.each <- sapply(eachdf1, function(x) sum(x$distrib)/nrow(x))
  AUC.single.all <- prod(AUC.single.each)
  
  bayesfac <- AUC.consol/AUC.single.all
  #bayesfac <- log(AUC.consol)-log(AUC.single)
  #bayesfac <- exp(bayesfac)
  bayesfac <- log(bayesfac)
  BF <- data.frame(bayesfac = bayesfac, consolevid = log(AUC.consol), specificevid = log(AUC.single.all))
  return(BF)
}

newbayesfun <- function(serotype, df, DS, agegroup) { # function to compare evidence of consol.inv for one DS vs single inv for one DS
  calc_invas(df %>% filter(Serotype == serotype) %>% filter(DS == DS))
  eachdf1 <- get(paste(DS, "distrib", agegroup, serotype, sep = "."))
  consol.distrib <- get(paste(paste("sero", serotype, sep = ""), "distrib", agegroup, sep = "."))
  consol.distrib$DS <- "consol"
  AUC.single <- sum(eachdf1$distrib)/nrow(eachdf1)
  AUC.consol <- sum(consol.distrib$distrib)/nrow(consol.distrib)
  bayesfac <- log(AUC.consol/AUC.single)
  return(bayesfac)
}

bf_levels <- function(val) { # function to add label to BF
  if (val > 100) {
    return("Extreme evidence for global")
  } else if (val < 100 & val > 30) {
    return("Very strong evidence for global")
  } else if (val < 30 & val > 10) {
    return("Strong evidence for global")
  } else if (val < 10 & val > 3) {
    return("Moderate evidence for global")
  } else if (val < 3 & val > 1) {
    return("Anecdotal evidence for global")
  } else if (val == 1) {
    return("No evidence")
  } else if (val < 1 & val > (1/3)) {
    return("Anecdotal evidence for local")
  } else if (val < (1/3) & val > (1/10)) {
    return("Moderate evidence for local")
  } else if (val < (1/10) & val > (1/30)) {
    return("Strong evidence for local")
  } else if (val < (1/30) & val > (1/100)) {
    return("Very strong evidence for local")
  } else {
    return("Extreme evidence for local")
  }
}

### BF for children 
sero15BC.distrib.children <- `sero15B/C.distrib.children`
bayes <- lapply(unique.child.sero, function(x) bayesfxr(serotype = x, df = inv.child.consol, agegroup = "children"))
bayes.df <- bind_rows(bayes)
bayes.df <- cbind(serotype = unique.child.sero, bayes.df)
bayes.df <- cbind(bayes.df, nDS = nDS[nDS > 1])
bayes.df$evidence <- unlist(lapply(bayes.df$bayesfac, FUN = bf_levels))
BFevidence <- bayes.df
# BFevidence <-  read.csv("c2c_evidence_BF_children.csv")
# BFevidence <- arrange(BFevidence, nDS)
# BFevidence <- arrange(BFevidence, Evidence)
# BFevidence$evidence <- factor(BFevidence$evidence, levels = c("Extreme evidence for local",
#                                                                       "Moderate evidence for local",
#                                                                       "Moderate evidence for global", 
#                                                                       "Strong evidence for global", "Very strong evidence for global",
#                                                                       "Extreme evidence for global"))
BFevidence$evidence <- factor(BFevidence$evidence, levels = c("Extreme evidence for local", "Very strong evidence for local",
                                                                      "Strong evidence for local", "Moderate evidence for local",
                                                                      "Anecdotal evidence for local", 
                                                                      "Anecdotal evidence for global" , "Moderate evidence for global", 
                                                                      "Strong evidence for global", "Very strong evidence for global",
                                                                      "Extreme evidence for global"))

colorval <- c("Extreme evidence for local"="#006D2C", 
              "Very strong evidence for local"="#31A354",
              "Strong evidence for local"="#74C476", 
              "Moderate evidence for local"="#BAE4B3",
              "Anecdotal evidence for local"="#EDF8E9", 
              "Anecdotal evidence for global"="#EFF3FF", 
              "Moderate evidence for global"="#BDD7E7", 
              "Strong evidence for global"="#6BAED6", 
              "Very strong evidence for global"="#3182BD",
              "Extreme evidence for global"="#08519C")
bayeschildren <- ggplot(BFevidence,aes(x = reorder(BFevidence$serotype, nDS), y = nDS)) + 
  geom_bar(stat = "identity", aes(fill = evidence)) + 
  #scale_fill_manual(values = c("#006D2C", "#BAE4B3", "#BDD7E7","#6BAED6","#3182BD", "#08519C")) + 
  scale_fill_manual(values = colorval,
                    limits = c("Extreme evidence for local", "Very strong evidence for local",
                               "Strong evidence for local", "Moderate evidence for local",
                               "Anecdotal evidence for local", 
                               "Anecdotal evidence for global" , "Moderate evidence for global", 
                               "Strong evidence for global", "Very strong evidence for global",
                               "Extreme evidence for global")) + 
  theme_minimal() + theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_y_continuous(breaks = function(x) seq(ceiling(x[1]), floor(x[2]), by = 1))+
  labs(x = "Serotype", y = "Number of datasets", fill = "Bayes Factor Label", title = "Evidence: Children") # pdf 4x12

### BF for adults
sero15BC.distrib.adults <- `sero15B/C.distrib.adults`
bayes_adu <- lapply(unique.adult.sero, function(x) bayesfxr(serotype = x, df = inv.adult.consol, agegroup = "adults"))
bayes_adu.df <- bind_rows(bayes_adu)
bayes_adu.df <- cbind(serotype = unique.adult.sero, bayes_adu.df)
bayes_adu.df <- cbind(bayes_adu.df, nDS = nDS.adult[nDS.adult > 1])
bayes_adu.df$evidence <- unlist(lapply(bayes_adu.df$bayesfac, FUN = bf_levels))
BFevidence_adu <- bayes_adu.df
# BFevidence_adu <- arrange(BFevidence_adu, nDS)
# BFevidence_adu <- arrange(BFevidence_adu, Evidence)
BFevidence_adu$evidence <- factor(BFevidence_adu$evidence, levels = c("Extreme evidence for local", "Very strong evidence for local",
                                                                      "Strong evidence for local", "Moderate evidence for local",
                                                                      "Anecdotal evidence for local", 
                                                                      "Anecdotal evidence for global" , "Moderate evidence for global", 
                                                                      "Strong evidence for global", "Very strong evidence for global",
                                                                      "Extreme evidence for global"))

bayesadults <- ggplot(BFevidence_adu,aes(x = reorder(BFevidence_adu$serotype, BFevidence_adu$nDS), y = BFevidence_adu$nDS)) + 
  geom_bar(stat = "identity", aes(fill = BFevidence_adu$evidence)) + 
  scale_fill_manual(values = colorval,
                    limits = c("Extreme evidence for local", "Very strong evidence for local",
                               "Strong evidence for local", "Moderate evidence for local",
                               "Anecdotal evidence for local", 
                               "Anecdotal evidence for global" , "Moderate evidence for global", 
                               "Strong evidence for global", "Very strong evidence for global",
                               "Extreme evidence for global"))+ 
  #c("#006D2C", "#31A354", "#74C476", "#BAE4B3", "#EDF8E9", "#EFF3FF", "#BDD7E7", "#6BAED6", "#3182BD", "#08519C")) + 
  theme_minimal() + theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_y_continuous(breaks = function(x) seq(ceiling(x[1]), floor(x[2]), by = 1))+
  labs(x = "Serotype", y = "Number of datasets", fill = "Bayes Factor Label", title = "Evidence: Adults") # pdf 4 x 12

#grid.arrange(bayeschildren, bayesadults, nrow = 2) # pdf 7 x 12
bf_bothage <- ggarrange(bayeschildren+ylab(''), bayesadults+ylab(''), nrow = 2, legend = "right", common.legend = TRUE) 
annotate_figure(bf_bothage, left = textGrob('Number of datasets', rot = 90))

#### BAYES FACTOR CHECKING CLOSED FORM SOLUTIONS WITH R AND PYTHON SOLUTIONS ----
invpri <- 0.5
numdeviates <- 1000

# serotype 27 local
sero27dummy <-  abs.child %>% filter(Serotype == "27", vaccinationera == "pre.PCV")
nDS27 <- nrow(sero27dummy)
Q27 <- sero27dummy$N*sero27dummy$time.int
prod27 <- prod((sero27dummy$n.swab +1)/(sero27dummy$carriage*Q27))
bayes27 <- log(prod27*(1/invpri^2)) # -15.35229
# serotype 27 consolidated
D.27 <- sum(sero27dummy$disease)
carrdev27 <- lapply(1:nDS27, 
                    function(x) rbeta(n = numdeviates, shape1 = sero27dummy$carriage[x] + 1, shape2 = sero27dummy$n.swab[x] - sero27dummy$carriage[x] + 1))
numerator27a <- lapply(1:nDS27, function(x) (Q27[x]*carrdev27[[x]])^sero27dummy$disease[x])
numerator27b <- numerator27a[[1]]*numerator27a[[2]]+numerator27a[[3]]*numerator27a[[4]]+numerator27a[[5]]*numerator27a[[6]]+numerator27a[[7]]*numerator27a[[8]]+
  numerator27a[[9]]*numerator27a[[10]]

denom27a <- lapply(1:nDS27, function(x) Q27[x]*carrdev27[[x]])
denom27b <- (denom27a[[1]] + denom27a[[2]])^(D.27+1)

complicatedpart27 <- sum(numerator27b/denom27b)/numdeviates

consbayes27 <- (1/invpri)*(factorial(D.27)/prod(factorial(sero27dummy$disease)))*complicatedpart27


# serotype 23B local
sero23Bdummy <-  abs.child %>% filter(Serotype == "23B", vaccinationera == "pre.PCV")
nDS23B <- nrow(sero23Bdummy)
Q23B <- sero23Bdummy$N*sero23Bdummy$time.int
prod23B <- prod((sero23Bdummy$n.swab +1)/(sero23Bdummy$carriage*sero23Bdummy$N*sero23Bdummy$time.int))
bayes23B <- log(prod23B*(1/invpri^2)) # -28.7849
# serotype 23B consolidated
D.23B <- sum(sero23Bdummy$disease)
carrdev23B <- lapply(1:nDS23B, 
                    function(x) rbeta(n = numdeviates, shape1 = sero23Bdummy$carriage[x] + 1, shape2 = sero23Bdummy$n.swab[x] - sero23Bdummy$carriage[x] + 1))
numerator23Ba <- lapply(1:nDS23B, function(x) (Q23B[x]*carrdev23B[[x]])^sero23Bdummy$disease[x])
numerator23Bb <- numerator23Ba[[1]]*numerator23Ba[[2]]

denom23Ba <- lapply(1:nDS23B, function(x) Q23B[x]*carrdev23B[[x]])
denom23Bb <- (denom23Ba[[1]] + denom23Ba[[2]])^(D.23B+1)

complicatedpart23B <- sum(numerator23Bb/denom23Bb)/numdeviates

consbayes23B <- (1/invpri)*(factorial(D.23B)/prod(factorial(sero23Bdummy$disease)))*complicatedpart23B #-20.711

# serotype 6B local
sero6Bdummy <- abs.child %>% filter(Serotype == "6B", vaccinationera == "pre.PCV")
nDS6B <- nrow(sero6Bdummy)
Q6B <- sero6Bdummy$N*sero6Bdummy$time.int
prod6B <- prod((sero6Bdummy$n.swab +1)/(sero6Bdummy$carriage*Q6B))
bayes6B <- log(prod6B*(1/invpri^2)) # -110.7084
# serotype 6B consolidated
D.6B <- sum(sero6Bdummy$disease)
carrdev6B <- lapply(1:nDS6B, 
                    function(x) rbeta(n = numdeviates, shape1 = sero6Bdummy$carriage[x] + 1, shape2 = sero6Bdummy$n.swab[x] - sero6Bdummy$carriage[x] + 1))
numerator6Ba <- lapply(1:nDS6B, function(x) sero6Bdummy$disease[x]*log(carrdev6B[[x]]*Q6B))
numerator6Bb <- numerator6Ba[[1]]+numerator6Ba[[2]]+numerator6Ba[[3]]+numerator6Ba[[4]]+numerator6Ba[[5]]+numerator6Ba[[6]]+numerator6Ba[[7]]+numerator6Ba[[8]]+
  numerator6Ba[[9]]+numerator6Ba[[10]]
denom6Ba <- lapply(1:nDS6B, function(x) Q6B[x]*carrdev6B[[x]])
denom6Bb <- (denom6Ba[[1]] + denom6Ba[[2]]+denom6Ba[[3]] + denom6Ba[[4]]+denom6Ba[[5]] + denom6Ba[[6]]+denom6Ba[[7]] + denom6Ba[[8]]+denom6Ba[[9]] + denom6Ba[[10]])
denom6Bc <- (D.6B+1)*log(denom6Bb)
newfrac6B <- exp(numerator6Bb-denom6Bc)
complicatedpart6B <- sum(newfrac6B)/numdeviates
D6B.fact <- log(factorialZ(D.6B))
proddfact6B <- log(prod(as.bigz(factorial(sero6Bdummy$disease))))
factorialpart6B <- exp(D6B.fact-proddfact6B)
consbayes6B <- log((1/invpri)*factorialpart6B*complicatedpart6B)

# serotype 6A local
sero6Adummy <- abs.child %>% filter(Serotype == "6A", vaccinationera == "pre.PCV")
nDS6A <- nrow(sero6Adummy)
Q6A <- sero6Adummy$N*sero6Adummy$time.int
prod6A <- prod((sero6Adummy$n.swab +1)/(sero6Adummy$carriage*Q6A))
bayes6A <- log(prod6A*(1/invpri^2)) # -107.6759
# serotype 6A consolidated
D.6A <- sum(sero6Adummy$disease)
carrdev6A <- lapply(1:nDS6A, 
                    function(x) rbeta(n = numdeviates, shape1 = sero6Adummy$carriage[x] + 1, shape2 = sero6Adummy$n.swab[x] - sero6Adummy$carriage[x] + 1))
numerator6Aa <- lapply(1:nDS6A, function(x) sero6Adummy$disease[x]*log(carrdev6A[[x]]*Q6A))
numerator6Ab <- numerator6Aa[[1]]+numerator6Aa[[2]]+numerator6Aa[[3]]+numerator6Aa[[4]]+numerator6Aa[[5]]+numerator6Aa[[6]]+numerator6Aa[[7]]+numerator6Aa[[8]]+
  numerator6Aa[[9]]+numerator6Aa[[10]]
denom6Aa <- lapply(1:nDS6A, function(x) Q6A[x]*carrdev6A[[x]])
denom6Ab <- (denom6Aa[[1]] + denom6Aa[[2]]+denom6Aa[[3]] + denom6Aa[[4]]+denom6Aa[[5]] + denom6Aa[[6]]+denom6Aa[[7]] + denom6Aa[[8]]+denom6Aa[[9]] + denom6Aa[[10]])
denom6Ac <- (D.6A+1)*log(denom6Ab)
newfrac6A <- exp(numerator6Ab-denom6Ac)
complicatedpart6A <- sum(newfrac6A)/numdeviates
D6A.fact <- log(factorialZ(D.6A))
proddfact6A <- log(prod(as.bigz(factorial(sero6Adummy$disease))))
factorialpart6A <- exp(D6A.fact-proddfact6A)
consbayes6A <- log((1/invpri)*factorialpart6A*complicatedpart6A) # 119.3598

# serotype 14 local
sero14dummy <- abs.child %>% filter(Serotype == "14", vaccinationera == "pre.PCV")
nDS14 <- nrow(sero14dummy)
Q14 <- sero14dummy$N*sero14dummy$time.int
prod14 <- prod((sero14dummy$n.swab +1)/(sero14dummy$carriage*Q14))
bayes14 <- log(prod14*(1/invpri^2)) # -110.7084
# serotype 14 consolidated
D.14 <- sum(sero14dummy$disease)
carrdev14 <- lapply(1:nDS14, 
                    function(x) rbeta(n = numdeviates, shape1 = sero14dummy$carriage[x] + 1, shape2 = sero14dummy$n.swab[x] - sero14dummy$carriage[x] + 1))
numerator14a <- lapply(1:nDS14, function(x) sero14dummy$disease[x]*log(carrdev14[[x]]*Q14))
numerator14b <- numerator14a[[1]]+numerator14a[[2]]+numerator14a[[3]]+numerator14a[[4]]+numerator14a[[5]]+numerator14a[[6]]+numerator14a[[7]]+numerator14a[[8]]+
  numerator14a[[9]]+numerator14a[[10]]+numerator14a[[11]]
denom14a <- lapply(1:nDS14, function(x) Q14[x]*carrdev14[[x]])
denom14b <- (denom14a[[1]] + denom14a[[2]]+denom14a[[3]] + denom14a[[4]]+denom14a[[5]] + denom14a[[6]]+denom14a[[7]] + denom14a[[8]]+denom14a[[9]] + denom14a[[10]] + 
  denom14a[[11]])
denom14c <- (D.14+1)*log(denom14b)
newfrac14 <- exp(numerator14b-denom14c)
complicatedpart14 <- sum(newfrac14)/numdeviates
D14.fact <- log(factorialZ(D.14))
proddfact14 <- log(prod(as.bigz(factorial(sero14dummy$disease))))
factorialpart14 <- exp(D14.fact-proddfact14)
consbayes14 <- log((1/invpri)*factorialpart14*complicatedpart14)

# serotype 9V local ## CARRIAGE = 0 so closed form solution not possible
sero9Vdummy <- abs.child %>% filter(Serotype == "9V", vaccinationera == "pre.PCV")
nDS9V <- nrow(sero9Vdummy)
Q9V <- sero9Vdummy$N*sero9Vdummy$time.int
prod9V <- prod((sero9Vdummy$n.swab +1)/(sero9Vdummy$carriage*Q9V))
bayes9V <- log(prod9V*(1/invpri^2)) # -110.7084
# serotype 9V consolidated
D.9V <- sum(sero9Vdummy$disease)
carrdev9V <- lapply(1:nDS9V, 
                    function(x) rbeta(n = numdeviates, shape1 = sero9Vdummy$carriage[x] + 1, shape2 = sero9Vdummy$n.swab[x] - sero9Vdummy$carriage[x] + 1))
numerator9Va <- lapply(1:nDS9V, function(x) sero9Vdummy$disease[x]*log(carrdev9V[[x]]*Q9V))
numerator9Vb <- numerator9Va[[1]]+numerator9Va[[2]]+numerator9Va[[3]]+numerator9Va[[4]]+numerator9Va[[5]]+numerator9Va[[6]]+numerator9Va[[7]]+numerator9Va[[8]]+
  numerator9Va[[9]]+numerator9Va[[10]]
denom9Va <- lapply(1:nDS9V, function(x) Q9V[x]*carrdev9V[[x]])
denom9Vb <- (denom9Va[[1]] + denom9Va[[2]]+denom9Va[[3]] + denom9Va[[4]]+denom9Va[[5]] + denom9Va[[6]]+denom9Va[[7]] + denom9Va[[8]]+denom9Va[[9]] + denom9Va[[10]])
denom9Vc <- (D.9V+1)*log(denom9Vb)
newfrac9V <- exp(numerator9Vb-denom9Vc)
complicatedpart9V <- sum(newfrac9V)/numdeviates
D9V.fact <- log(factorialZ(D.9V))
proddfact9V <- log(prod(as.bigz(factorial(sero9Vdummy$disease))))
factorialpart9V <- exp(D9V.fact-proddfact9V)
consbayes9V <- log((1/invpri)*factorialpart9V*complicatedpart9V)


# serotype 22F local
sero22Fdummy <-  abs.child %>% filter(Serotype == "22F", vaccinationera == "pre.PCV")
nDS22F <- nrow(sero22Fdummy)
Q22F <- sero22Fdummy$N*sero22Fdummy$time.int
prod22F <- prod((sero22Fdummy$n.swab +1)/(sero22Fdummy$carriage*sero22Fdummy$N*sero22Fdummy$time.int))
bayes22F <- log(prod22F*(1/invpri^2)) # -18.21582
# serotype 22F consolidated - NOT POSSIBLE BECAUSE CARRIAGE = 0 FOR ONE OF THE DS
D.22F <- sum(sero22Fdummy$disease)
carrdev22F <- lapply(1:nDS22F, 
                     function(x) rbeta(n = numdeviates, shape1 = sero22Fdummy$carriage[x] + 1, shape2 = sero22Fdummy$n.swab[x] - sero22Fdummy$carriage[x] + 1))
numerator22Fa <- lapply(1:nDS22F, function(x) (Q22F[x]*carrdev22F[[x]])^sero22Fdummy$disease[x])
numerator22Fb <- numerator22Fa[[1]]*numerator22Fa[[2]]*numerator22Fa[[3]]

denom22Fa <- lapply(1:nDS22F, function(x) Q22F[x]*carrdev22F[[x]])
denom22Fb <- (denom22Fa[[1]] + denom22Fa[[2]]+denom22Fa[[3]])^(D.22F+1)

complicatedpart22F <- sum(numerator22Fb/denom22Fb)/numdeviates

consbayes22F <- (1/invpri)*(factorial(D.22F)/prod(factorial(sero22Fdummy$disease)))*complicatedpart22F


# #### testing bayes fxr
# testdata <- child.abs.ds %>% filter(Serotype == "6A") %>% filter(DS == "Bogota.pre.PCV")
# test <- invasiveness(testdata)
# BF.df <- bayesfxr.df[(order(bayesfxr.df$invas_dev)),]
# sumlike <- (sum(BF.df$posterior_samp))/nrow(BF.df) # 1.481264e-05
# AUC.test <- AUC(invas_dev, posterior_samp, method = "trapezoid") # 7.180451e-06
# 
# testdata2 <- child.abs.ds %>% filter(Serotype == "6B") %>% filter(DS == "Bogota.pre.PCV")
# test2 <- invasiveness(testdata2)
# BF.df <- bayesfxr.df[(order(bayesfxr.df$invas_dev)),]
# sumlike <- (sum(BF.df$posterior_samp))/nrow(BF.df) # 1.691704e-05
# AUC.test <- AUC(BF.df$invas_dev, BF.df$posterior_samp, method = "trapezoid") # 1.024761e-05
# 
# testdata3 <- child.abs.ds %>% filter(Serotype == "6A") %>% filter(DS == "Atlanta.pre.PCV")
# test3 <- invasiveness(testdata3)
# BF.df <- bayesfxr.df[(order(bayesfxr.df$invas_dev)),]
# sumlike <- (sum(BF.df$posterior_samp))/nrow(BF.df) # 0.0001593778
# AUC.test <- AUC(BF.df$invas_dev, BF.df$posterior_samp, method = "trapezoid") # 6.824021e-05
# 
# testdata4 <- child.abs.ds %>% filter(Serotype == "6A") %>% filter(DS == "Atlanta.pre.PCV")
# test4 <- invasiveness(testdata4)
# BF.df <- bayesfxr.df[(order(bayesfxr.df$invas_dev)),]
# sumlike <- (sum(BF.df$posterior_samp))/nrow(BF.df) # 0.0001593778
# AUC.test <- AUC(BF.df$invas_dev, BF.df$posterior_samp, method = "trapezoid") # 6.824021e-05


#### GPSC data ---------------------------------------------------------------------------------------------------------------------------------------------

GPSC_data <-  read_excel("GPSC_serotype_data.xlsx")
colnames(GPSC_data) <- c("laneID", "Serotype", "GPSC")
GPSC_data$laneID <- NULL
gpsc.unique.sero <- unique(GPSC_data$Serotype)
num.strains <- map(gpsc.unique.sero, function(x) length(unique((GPSC_data %>% filter(Serotype == x))$GPSC)))
consol.strains <- do.call("rbind", num.strains)
consol.gpsc <- data.frame(cbind(gpsc.unique.sero, consol.strains))
colnames(consol.gpsc) <- c("Serotype", "no.strains")
consol.gpsc$no.strains <- as.numeric(as.character(consol.gpsc$no.strains))
consol.gpsc <- consol.gpsc[order(-consol.gpsc$no.strains),]

# child age group
new.child.comb <- dplyr::left_join(consol.gpsc, consol.child) %>% drop_na()
new.child.comb$Serogroup <- factor(new.child.comb$Serogroup, levels = c("VT7", "VT10", "VT13","NVT")) 

ggplot(new.child.comb) + # pdf dim 5 x 10
  geom_smooth(aes(x= no.strains, y = overall.invasiveness, weight = 0.5), method = 'lm', se = F, colour = 'darkgrey', 
              linetype = 'dashed') +
  geom_errorbar(aes(x = no.strains, ymin = overall.invasiveness.low, ymax = overall.invasiveness.high,
                    group = Serogroup, color = Serogroup), 
                width = 0) +
  geom_point(aes(x = no.strains, y = overall.invasiveness, color = Serogroup, group = Serogroup), size = 2) + 
  #geom_text(aes(x = no.strains, y = overall.invasiveness, label = Serotype), 
  #          size = 2, color = "white", fontface = "bold") +
  ggrepel::geom_text_repel(aes(x = no.strains, y = overall.invasiveness, label = Serotype, group = Serogroup, 
                               colour = Serogroup), fontface = "bold",show.legend = FALSE) +
  scale_color_manual(values= cols) + scale_y_continuous(trans = 'log10') +
  labs(x = "No. of strains", y = "Global invasiveness", colour = "Serotype category") + #theme_minimal() +
  theme_classic() + theme_bw() + theme(text = element_text(size = 15)) + 
  stat_cor(aes(x= no.strains, y = overall.invasiveness, label = paste(..rr.label.., ..p.label.., sep = "~`,`~")), label.x = 50)
  #
  #stat_ellipse(aes(x = no.strains, y = overall.invasiveness, color = Serogroup, group = Serogroup, type = "norm"))

# adult age group
new.adult.comb <- dplyr::left_join(consol.gpsc, consol.adult) %>% drop_na()
new.adult.comb$Serogroup <- factor(new.adult.comb$Serogroup, levels = c("VT7", "VT10", "VT13","NVT"))

ggplot(new.adult.comb, aes(x = no.strains, y = overall.invasiveness)) +# pdf dim 5 x 10
  geom_smooth(aes(x= no.strains, y = overall.invasiveness, weight = 0.5), method = 'lm', se = F, colour = 'darkgrey', 
              linetype = 'dashed') +
  geom_errorbar(aes(x = no.strains, ymin = overall.invasiveness.low, ymax = overall.invasiveness.high,
                    group = Serogroup, color = Serogroup), 
                width = 0) +
  geom_point(aes(x = no.strains, y = overall.invasiveness, color = Serogroup, group = Serogroup), size = 2) + 
  #geom_text(aes(x = no.strains, y = overall.invasiveness, label = Serotype), 
  #          size = 2, color = "white", fontface = "bold") +
  ggrepel::geom_text_repel(aes(x = no.strains, y = overall.invasiveness, label = Serotype, group = Serogroup, 
                               colour = Serogroup), fontface = "bold",show.legend = FALSE) +
  scale_color_manual(values= cols) + scale_y_continuous(trans = 'log10') +
  labs(x = "No. of strains", y = "Global invasiveness", colour = "Serotype category") + #theme_minimal() +
  theme_classic() + theme_bw() +theme(text = element_text(size = 15)) +
  stat_cor(aes(x= no.strains, y = overall.invasiveness, label = paste(..rr.label.., ..p.label.., sep = "~`,`~")), label.x = 50)
#stat_ellipse(aes(x = no.strains, y = overall.invasiveness, color = Serogroup, group = Serogroup, type = "norm"))

# CALC GPSC SDI and plot vs invasiveness
# requires import of calc_SimpsDiv from serotype replacement file
GPSC_data <-  read_excel("GPSC_serotype_data.xlsx")
colnames(GPSC_data) <- c("laneID", "Serotype", "GPSC")
GPSC_data$laneID <- NULL
GPSC_div <- map(gpsc.unique.sero, function(x) calc_SimpsDiv((GPSC_data %>% filter(Serotype == x) %>% group_by(GPSC) %>% tally())$n))
GPSC_div_allsero <- do.call("rbind", GPSC_div)
GPSC_div_allsero <- data.frame(cbind(gpsc.unique.sero, GPSC_div_allsero)) %>% drop_na() # NaNs result from serotypes that only had 1 laneID
colnames(GPSC_div_allsero) <- c("Serotype", "SDI", "SDI_low", "SDI_high", "SDI_var")
new_child_GPSC_div <- left_join(GPSC_div_allsero, consol.child) %>% drop_na()
new_child_GPSC_div$Serogroup <- factor(new_child_GPSC_div$Serogroup, levels = c("VT7", "VT10", "VT13","NVT")) 

ggplot(new_child_GPSC_div) + # pdf dim 5 x 10
  geom_point(aes(x = SDI, y = overall.invasiveness, color = Serogroup, group = Serogroup), size = 4) + 
  geom_text(aes(x = SDI, y = overall.invasiveness, label = Serotype), size = 2, color = "white", fontface = "bold") +
  geom_errorbar(aes(x = SDI, ymin = overall.invasiveness.low, ymax = overall.invasiveness.high), width = 0.0001) +
  geom_errorbarh(aes(y = overall.invasiveness, xmin = SDI_low, xmax = SDI_high), height = 0.001)
  scale_color_manual(values= cols) +
  labs(x = "Simpson's Diversity Index", y = "Global Invasiveness")

#### --------------------------NOT USED/COMMENTED OUT-------------------------------------------------------------------------------------------------------

# Consolidated invasiveness for combined age groups in pre-PCV time periods
combined.abs.ds <- rbind(child.abs.ds, adult.abs.ds)
combined.inv.consol <- data.frame(map(combined.abs.ds, unlist)) %>% 
  filter(!DS %in% c("France.post.PCV7", "Barcelona.post.PCV7", "France.post.PCV13", "Bogota.post.PCV7", "Massachusetts.post.PCV7","Navajo.post.PCV7", 
                    "Atlanta.post.PCV7")) %>% group_by(Serotype) %>% nest()
nDS.combo <- sapply(combined.inv.consol$data, function(x) nrow(x)) # number of datasets for each serotype
excl.sero.combo <- combined.inv.consol$Serotype[which(nDS.combo == 1)] # serotypes with only 1 dataset to be excluded for analysis
combined.inv.consol <- data.frame(combined.inv.consol %>% unnest() %>% filter(!Serotype %in% excl.sero.combo))

unique.combined.sero <- unique(combined.inv.consol$Serotype)

consol.combined2 <- lapply(unique.combined.sero, function(x) {calc_invas(combined.inv.consol %>% filter(Serotype == x))})
consol.combined <- do.call("rbind", consol.combined2)
colnames(consol.combined) <- c("Serotype", "overall.invasiveness", "overall.invasiveness.low", "overall.invasiveness.high")

consol.combined$Serogroup <- NA
consol.combined$Serogroup[which(consolidated.combined.inv$Serotype %in% c("4", "6B", "9V", "14", "18C", "19F", "23F"))] <- "VT7"
consol.combined$Serogroup[which(consolidated.combined.inv$Serotype %in% c("1", "5", "7F"))] <- "VT10"
consol.combined$Serogroup[which(consolidated.combined.inv$Serotype %in% c("3", "6A", "19A"))] <- "VT13"
consol.combined$Serogroup[is.na(consolidated.combined.inv$Serogroup)] <- "NVT"

plot.consol.inv(consolidated.combined.inv)

#### ellipse 

#ggplot(abs.child) + 
#  geom_point(aes(x = IPD.inc, y = invasiveness, color = DS, group = DS)) + 
#  stat_ellipse(aes(x = IPD.inc, y = invasiveness, color = DS, group = DS, type = "norm"))

#ggplot(abs.adult) + 
#  geom_point(aes(x = carr.prev, y = IPD.inc, color = DS, group = DS)) + 
#  stat_ellipse(aes(x = carr.prev, y = IPD.inc, color = DS, group = DS, type = "norm"))

plot.carr <- function(sero) { # dataset vs carriage prevalence for *one* serotype & *one* age group
  x.DS <- unlist(sero$DS)
  y.carrprev <- unlist(sero$carr.prev)
  ggplot(sero, aes(x = reorder(x.DS, -y.carrprev), y = y.carrprev)) + 
    geom_point() + 
    geom_errorbar(aes(ymin = unlist(sero$carr.prev.low), ymax = unlist(sero$carr.prev.high)), width = 0.01) +
    labs(x = "Dataset", y = "Carriage Prevalence", title = sero$Serotype[1]) +
    theme_bw() + theme_light()
}

plot.carr.allages <- function(sero) { # dataset vs carriage prevalence for *one* serotype for *both* age groups
  x.DS <- unlist(sero$DS)
  y.carrprev <- unlist(sero$carr.prev)
  ggplot(sero, aes(x = reorder(x.DS, -y.carrprev), y = y.carrprev, color = unlist(sero$agegrp))) + 
    geom_point() + 
    geom_errorbar(aes(ymin = unlist(sero$carr.prev.low), ymax = unlist(sero$carr.prev.high)), width = 0.01) +
    labs(x = "Dataset", y = "Carriage Prevalence", title = sero$Serotype[1], color = "Age Group") +
    theme_bw() + theme_light()
}

plot.dis <- function(sero) { # dataset vs IPD incidence for *one* serotype & *one* age group
  x.DS <- unlist(sero$DS)
  y.IPDinc <- unlist(sero$IPD.inc)
  ggplot(sero, aes(x = reorder(x.DS, -y.IPDinc), y = y.IPDinc)) + 
    geom_point() + 
    geom_errorbar(aes(ymin = unlist(sero$IPD.inc.low), ymax = unlist(sero$IPD.inc.high)), width = 0.01) +
    labs(x = "Dataset", y = "IPD Incidence", title = sero$Serotype[1]) +
    theme_bw() + theme_light()
}

plot.dis.allages <- function(sero) { # dataset vs IPD incidence for *one* serotype for *both* age groups
  x.DS <- unlist(sero$DS)
  y.IPDinc <- unlist(sero$IPD.inc)
  ggplot(sero, aes(x = reorder(x.DS, -y.IPDinc), y = y.IPDinc, color = unlist(sero$agegrp))) + 
    geom_point() + 
    geom_errorbar(aes(ymin = unlist(sero$IPD.inc.low), ymax = unlist(sero$IPD.inc.high)), width = 0.01) +
    labs(x = "Dataset", y = "IPD Incidence", title = sero$Serotype[1], color = "Age Group") +
    theme_bw() + theme_light()
}

plot.invasiveness <- function(sero) { # dataset vs invasiveness for *one* serotype & *one* age group
  x.DS <- unlist(sero$DS)
  y.inv <- unlist(sero$invasiveness)
  ggplot(sero, aes(x = reorder(x.DS, -y.inv), y = y.inv)) + 
  geom_point() + 
  geom_errorbar(aes(ymin = unlist(sero$invasiveness.low), ymax = unlist(sero$invasiveness.high)), width = 0.01) +
  labs(x = "Dataset", y = "Invasiveness", title = sero$Serotype[1]) +
  theme_bw() + theme_light()
}

plot.inv.allages <- function(sero) { # dataset vs invasiveness for *one* serotype for *both* age groups
  x.DS <- unlist(sero$DS)
  y.inv <- unlist(sero$invasiveness)
  sero$agegrp <- factor(sero$agegrp, levels = c("infant","adult"))
  ggplot(sero, aes(x = reorder(x.DS, -y.inv), y = y.inv, color = unlist(sero$agegrp))) + 
    geom_point() + 
    geom_errorbar(aes(ymin = unlist(sero$invasiveness.low), ymax = unlist(sero$invasiveness.high)), width = 0.01) +
    labs(x = "Dataset", y = "Invasiveness", title = sero$Serotype[1], color = "Age Group") +
    theme_bw() + theme_light()
}

plot.inv.allages.DS <- function(sero) { # serotype vs invasiveness for *one* dataset for *both* age groups
  x.Sero <- unlist(sero$Serotype)
  y.inv <- unlist(sero$invasiveness)
  sero$agegrp <- factor(sero$agegrp, levels = c("infant","adult"))
  ggplot(sero, aes(x = reorder(x.Sero, -y.inv), y = y.inv, color = unlist(sero$agegrp))) + 
    geom_point() + 
    geom_errorbar(aes(ymin = unlist(sero$invasiveness.low), ymax = unlist(sero$invasiveness.high)), width = 0.01) +
    labs(x = "Serotype", y = "Invasiveness", title = sero$DS[1], color = "Age Group") +
    theme_bw() + theme_light()
}

# Invasiveness for children only
abs.child <- as.data.frame(t(sapply(1:nrow(child.abs.ds), function(x) invasiveness(child.abs.ds[x,]))))
abs.child.2 <- abs.child %>% dplyr::select(DS, Serotype, invasiveness, invasiveness.low, invasiveness.high)
abs.child.sero <- unlist(unique(abs.child.2$Serotype))
abschildgraphs <- lapply(abs.child.sero, function(x) {plot.invasiveness(abs.child.2 %>% filter(Serotype == x))})
mypath.abschild <- file.path("Q:","Technical","R","Case-to-carrier",paste("Absolute_inv_children", ".pdf", sep=""))
pdf(mypath.abschild, width = 20, height = 5)
print(abschildgraphs)
dev.off()

# Carriage prevalence for children + graphs made into one PDF
abs.child.carr <- abs.child %>% dplyr::select(DS, Serotype, carr.prev, carr.prev.low, carr.prev.high)
abschildgraphs.carr <- lapply(abs.child.sero, function(x) {plot.carr(abs.child.carr %>% filter(Serotype == x))})
mypath.abschild.carr <- file.path("Q:","Technical","R","Case-to-carrier",paste("Absolute_carr_children", ".pdf", sep=""))
pdf(mypath.abschild.carr, width = 20, height = 5)
print(abschildgraphs.carr)
dev.off()

# Disease incidence for children + graphs made into one PDF
abs.child.dis <- abs.child %>% dplyr::select(DS, Serotype, IPD.inc, IPD.inc.low, IPD.inc.high)
abschildgraphs.dis <- lapply(abs.child.sero, function(x) {plot.dis(abs.child.dis %>% filter(Serotype == x))})
mypath.abschild.dis <- file.path("Q:","Technical","R","Case-to-carrier",paste("Absolute_dis_children", ".pdf", sep=""))
pdf(mypath.abschild.dis, width = 20, height = 5)
print(abschildgraphs.dis)
dev.off()

# Invasiveness for adults relative to children's carriage
abs.adult <- as.data.frame(t(sapply(1:nrow(adult.abs.ds), function(x) invasiveness(adult.abs.ds[x,]))))
abs.adult.2 <- abs.adult %>% dplyr::select(DS, Serotype, invasiveness, invasiveness.low, invasiveness.high)
abs.adult.sero <- unlist(unique(abs.adult.2$Serotype))
absadultgraphs <- lapply(abs.adult.sero, function(x) {plot.invasiveness(abs.adult.2 %>% filter(Serotype == x))})
mypath.absadult <- file.path("Q:","Technical","R","Case-to-carrier",paste("Absolute_inv_adults", ".pdf", sep=""))
pdf(mypath.absadult, width = 20, height = 5)
print(absadultgraphs)
dev.off()

# Carriage prevalence for adults + graphs made into one PDF
abs.adult.carr <- abs.adult %>% dplyr::select(DS, Serotype, carr.prev, carr.prev.low, carr.prev.high)
absadultgraphs.carr <- lapply(abs.adult.sero, function(x) {plot.carr(abs.adult.carr %>% filter(Serotype == x))})
mypath.absadult.carr <- file.path("Q:","Technical","R","Case-to-carrier",paste("Absolute_carr_adults", ".pdf", sep=""))
pdf(mypath.absadult.carr, width = 20, height = 5)
print(absadultgraphs.carr)
dev.off()

# Disease incidence for adults + graphs made into one PDF
abs.adult.dis <- abs.adult %>% dplyr::select(DS, Serotype, IPD.inc, IPD.inc.low, IPD.inc.high)
absadultgraphs.dis <- lapply(abs.adult.sero, function(x) {plot.dis(abs.adult.dis %>% filter(Serotype == x))})
mypath.absadult.dis <- file.path("Q:","Technical","R","Case-to-carrier",paste("Absolute_dis_adults", ".pdf", sep=""))
pdf(mypath.absadult.dis, width = 20, height = 5)
print(absadultgraphs.dis)
dev.off()

# Absolute invasiveness of combined children + adults in one graph
abs.child.inv.tot <- abs.child %>% dplyr::select(agegrp, DS, Serotype, invasiveness, invasiveness.low, invasiveness.high, IPD.inc, IPD.inc.low, IPD.inc.high)
abs.adult.inv.tot <- abs.adult %>% dplyr::select(agegrp, DS, Serotype, invasiveness, invasiveness.low, invasiveness.high, IPD.inc, IPD.inc.low, IPD.inc.high)
abs.inv.tot <- bind_rows(abs.child.inv.tot, abs.adult.inv.tot)
abs.sero <- unique(c(abs.child.inv.tot$Serotype, abs.adult.inv.tot$Serotype))
abs.tot.graphs <- lapply(unlist(abs.sero), function(x) {plot.inv.allages(abs.inv.tot %>% filter(Serotype == x))})
mypath.abs.tot.inv <- file.path("Q:","Technical","R","Case-to-carrier",paste("Absolute_inv_tot", ".pdf", sep=""))
pdf(mypath.abs.tot.inv, width = 20, height = 5)
print(abs.tot.graphs)
dev.off()

abs.DS <- unique(c(abs.child.inv.tot$DS, abs.adult.inv.tot$DS))
abs.tot.graphs.DS <- lapply(abs.DS, function(x) {plot.inv.allages.DS(abs.inv.tot %>% filter(DS == x))})
mypath.abs.tot.inv.DS <- file.path("Q:","Technical","R","Case-to-carrier",paste("Absolute_inv_tot_perDS", ".pdf", sep=""))
pdf(mypath.abs.tot.inv.DS, width = 20, height = 5)
print(abs.tot.graphs.DS)
dev.off()

abs.tot.graphs.allDS <- plot.inv.allages.allDS(abs.inv.tot)
mypath.abs.tot.allDS <- file.path("Q:","Technical","R","Case-to-carrier",paste("Absolute_inv_tot_allDS_logscale", ".pdf", sep=""))
pdf(mypath.abs.tot.allDS, width = 25, height = 5)
print(abs.tot.graphs.allDS)
dev.off()

# Absolute IPD incidence of combined children + adults in one graph
abs.child.dis.tot <- abs.child %>% dplyr::select(agegrp, DS, Serotype, IPD.inc, IPD.inc.low, IPD.inc.high)
abs.adult.dis.tot <- abs.adult %>% dplyr::select(agegrp, DS, Serotype, IPD.inc, IPD.inc.low, IPD.inc.high)
abs.dis.tot <- bind_rows(abs.child.dis.tot, abs.adult.dis.tot)
abs.dis.tot.graphs <- lapply(abs.sero, function(x) {plot.dis.allages(abs.dis.tot %>% filter(Serotype == x))})
mypath.abs.tot.dis <- file.path("Q:","Technical","R","Case-to-carrier",paste("Absolute_dis_tot", ".pdf", sep=""))
pdf(mypath.abs.tot.dis, width = 20, height = 5)
print(abs.dis.tot.graphs)
dev.off()

abs.dis.tot.graphs.allDS <- plot.dis.allages.allDS(abs.dis.tot)
mypath.abs.tot.dis.allDS <- file.path("Q:","Technical","R","Case-to-carrier",paste("Absolute_dis_tot_allDS", ".pdf", sep=""))
pdf(mypath.abs.tot.dis.allDS, width = 25, height = 8)
print(abs.dis.tot.graphs.allDS)
dev.off()

something <- child.abs.ds %>% filter(!DS %in% c("France.post.PCV7", "Barcelona.post.PCV7", "France.post.PCV13", "Bogota.post.PCV7", "Massachusetts.post.PCV7",
                                                "Navajo.post.PCV7", "Atlanta.post.PCV7")) %>% mutate(newN = N*time.int) %>% group_by(Serotype, Serogroup) %>% nest()
carriage <- lapply(something$data, function(x) {sum(x$carriage)})
disease <- lapply(something$data, function(x) {sum(x$disease)})
n.swab <- lapply(something$data, function(x) {sum(x$n.swab)})
N <- lapply(something$data, function(x) {sum(x$newN)})
something$carriage <- unlist(carriage)
something$disease <- unlist(disease)
something$n.swab <- unlist(n.swab)
something$N <- unlist(N)
something$time.int <- rep(1, nrow(something))
something$agegrp <- rep("children", nrow(something))
something$DS <- rep(NA, nrow(something))
something$data <- NULL
something.invasiveness <- data.frame(t(sapply(1:nrow(something), function(x) invasiveness(something[x,]))))
ggplot(something.invasiveness, aes(x = reorder(unlist(Serotype), -(unlist(invasiveness))), y = unlist(invasiveness))) + 
  geom_point() + 
  geom_errorbar(aes(ymin = unlist(invasiveness.low), ymax = unlist(invasiveness.high)), width = 0.01) + 
  coord_cartesian(ylim = c(0,0.0035)) +
  labs(x = "Serotype", y = "Invasiveness", title = "Total Consolidated Invasiveness") + theme_bw() + theme_light()

#lapply(something$data, function(x) {mutate(x, new.N=N*time.int)})

# GPSC Data: Number of strains vs Combined Age Groups Invasiveness -- NOT USED / COMMENTED OUT -
#new.comb <- dplyr::left_join(consol.gpsc, consolidated.combined.inv) %>% drop_na()
#new.comb$Serogroup <- factor(new.comb$Serogroup, levels = c("VT7", "VT10", "VT13","NVT")) 

#ggplot(new.comb) + 
#  geom_point(aes(x = no.strains, y = overall.invasiveness*100000, color = Serogroup, group = Serogroup), size = 4) + 
#  geom_text(aes(x = no.strains, y = overall.invasiveness*100000, label = Serotype), size = 2, color = "white", fontface = "bold") +
#  scale_color_manual(values= cols) +
#  labs(x = "Number of strains", y = "Overall Invasiveness")
#stat_ellipse(aes(x = no.strains, y = overall.invasiveness, color = Serogroup, group = Serogroup, type = "norm"))

#### ------------------- Wrong way to calculate invasiveness + overal invasiveness --
#invasiveness_wrong <- function(row_val) { #WRONG INVASIVENESS FUNCTION - NOT USED
  
  set.seed(123)
  
  #row_val <- adult.abs.ds[80,]
  agegroup <- row_val[["agegrp"]]
  DS <- row_val[["DS"]] # dataset/study location
  serogroup <- row_val[["Serogroup"]]
  serotype <- row_val[["Serotype"]]
  n.carriage <- as.numeric(row_val[["carriage"]]) # number of carriers
  n.disease <- as.numeric(row_val[["disease"]]) # number of disease cases
  N <- as.numeric(row_val[["N"]]) # population
  n.swab <- as.numeric(row_val[["n.swab"]]) # number of swabs
  time.int <- as.numeric(row_val[["time.int"]]) # time interval
  
  # carriage prevalence ---
  
  LL_carr <- function(x) { # function computing LL of binomial distribution for carriage prevalence
    prob <- x[1]
    ll_carr <- dbinom(n.carriage, size = n.swab, prob, log = TRUE)
    return(ll_carr)
    }
  
  initial_val_carr <- n.carriage/n.swab
  
  #run MC on LL function for carriage prevalence
  results_carr <- met_gaussian(LL_carr, no_var = 1, ini_value = initial_val_carr, iters = 10000, stepsizes_met = 0.01*initial_val_carr) 
  LL_values_carr <- apply(results_carr,MARGIN=1,FUN=LL_carr) # get a LL value for each iteration
  results_carr_LL <- cbind(results_carr,LL_values_carr)
  
  #pairs(results_carr_LL)
  #hist(results_carr[,1]) # histogram of results
  
  MLE_index_carr <- which.max(LL_values_carr)
  carr_prev <- results_carr[MLE_index_carr,]
  results_carr_LL <- results_carr_LL[order(results_carr_LL[,1]), ]
  carr_prev_low <- results_carr_LL[250]
  carr_prev_high <- results_carr_LL[9750]
  
  # disease incidence ---
  
  set.seed(123)

  LL_dis <- function(x) { # function estimating gamma = mean count rate, with units [cases][person^-1][time^-1] (makes lambda unit = cases)
    gamma <- x[1] # average disease rate
    ll_dis <- dpois(n.disease, gamma*time.int*N, log = TRUE)
    return(ll_dis)
  }
  
  initial_val_dis <- n.disease/(time.int*N)
  
  #run MC on LL function for disease rate
  results_dis <- met_gaussian(LL_dis, no_var = 1, ini_value = initial_val_dis, iters = 10000, stepsizes_met = 0.001*initial_val_dis)
  LL_values_dis <- apply(results_dis, MARGIN=1, FUN=LL_dis)
  results_dis_LL <- cbind(results_dis, LL_values_dis)
  
  #pairs(results_dis_LL)
  #hist(results_dis[,1])
  
  MLE_index_dis <- which.max(LL_values_dis)
  dis_rate <- results_dis[MLE_index_dis,]
  results_dis_LL <- results_dis_LL[order(results_dis_LL[,1]), ]
  dis_rate_low <- results_dis_LL[250]
  dis_rate_high <- results_dis_LL[9750]
  
  # calculate invasiveness ---
  
  inv.df <- data.frame(results_carr, results_dis)
  invas.df <- inv.df %>% transmute(invasiveness = results_dis/results_carr)
  invas.df <- invas.df[order(invas.df[,1]), ]
  invasiveness.low <- invas.df[250]
  invasiveness.high <- invas.df[9750]
  
  invasiveness.mean <- dis_rate/carr_prev
  
  # return function ---

  new_row <- data.frame(row_val, carr.prev = carr_prev, carr.prev.low = carr_prev_low, carr.prev.high = carr_prev_high, 
                        IPD.inc = dis_rate, IPD.inc.low = dis_rate_low, IPD.inc.high = dis_rate_high, 
                        invasiveness = invasiveness.mean, invasiveness.low = invasiveness.low, invasiveness.high = invasiveness.high)
  return(new_row)
#}

#overall_invasiveness <- function(df) { ### NOT USED
  set.seed(123)
  n <- nrow(df)
  LL_inv <- function(x) { # function estimating overall invasiveness for each serotype
    carr_prev <- x[1:n]
    invas <- x[n+1]
    
    if(any(carr_prev < 0.0000001)) return(-Inf)
    if(invas < 0) return(-Inf)

    ll_inv <- sum(sapply(seq(1:n), function(s) {
      dpois(df[s, 6], carr_prev[s]*invas*df[s,9]*df[s,8], log = TRUE) + 
      dbeta(carr_prev[s], shape1 = df[s,5] + 1, shape2 = df[s,7] - df[s,5] + 1, log = TRUE)
    }))
    
    return(ll_inv)
  }
  
  initial_val_carr <- (df$carriage+1)/df$n.swab
  initial_val_invas <- (mean(df$disease)*mean(df$n.swab))/(mean(df$time.int)*mean(df$N)*(mean(df$carriage)+1))
  
  #run MC on LL function for disease rate
  results_invas <- met_gaussian(LL_inv, no_var = n+1, ini_value = c(initial_val_carr, initial_val_invas), 
                                iters = 20000, stepsizes_met = c(0.05*initial_val_carr, 0.05*initial_val_invas))
  LL_values_invas <- apply(results_invas, MARGIN=1, FUN=LL_inv)
  results_invas_LL <- cbind(results_invas, LL_values_invas)
  
  pairs(results_invas_LL)
  #hist(results_invas[,n+1])
  
  MLE_index_invas <- which.max(LL_values_invas)
  overall.invas <- results_invas_LL[MLE_index_invas,n+1]
  results_invas_LL <- results_invas_LL[order(results_invas_LL[,n+1]), ]
  overall.invas.low <- results_invas_LL[250,n+1]
  overall.invas.high <- results_invas_LL[9750,n+1]
  
  inv_row <- data.frame(overall.invasiveness = overall.invas, overall.invasiveness.low = overall.invas.low, overall.invasiveness.high = overall.invas.high)
  
  return(inv_row)#}

  
### Consolidated Invasiveness / WRONG / NOT USED / COMMENTED OUT ---
#inv.child.consol <- data.frame(map(abs.child, unlist)) %>% 
#  filter(!DS %in% c("France.post.PCV7", "Barcelona.post.PCV7", "France.post.PCV13", "Bogota.post.PCV7", "Massachusetts.post.PCV7","Navajo.post.PCV7", 
#                    "Atlanta.post.PCV7")) %>% 
#  mutate(newN = unlist(N)*unlist(time.int)) %>% 
#  mutate(overall.inv.nom.perDS = invasiveness*newN) %>%
#  group_by(Serotype, Serogroup) %>% nest()

#inv.child.consol$overall.inv.denom <- unlist(lapply(inv.child.consol$data, function(x) {sum(x[["newN"]])}))
#inv.child.consol$overall.inv.nom <- unlist(lapply(inv.child.consol$data, function(x) {sum(x[["overall.inv.nom.perDS"]])}))
#inv.child.consol <- inv.child.consol %>% mutate(overall.inv = overall.inv.nom/overall.inv.denom) %>% unnest()

#inv.adult.consol <- data.frame(map(abs.adult, unlist)) %>% 
#  filter(!DS %in% c("France.post.PCV7", "Barcelona.post.PCV7", "France.post.PCV13", "Bogota.post.PCV7", "Massachusetts.post.PCV7","Navajo.post.PCV7", 
#                    "Atlanta.post.PCV7")) %>% 
#  mutate(newN = unlist(N)*unlist(time.int)) %>% 
#  mutate(overall.inv.nom.perDS = invasiveness*newN) %>%
#  mutate(overall.inv.nom.perDS.low = invasiveness.low*newN) %>%
#  mutate(overall.inv.nom.perDS.high = invasiveness.high*newN) %>%
#  group_by(Serotype, Serogroup) %>% nest()

#inv.adult.consol$overall.inv.denom <- unlist(lapply(inv.adult.consol$data, function(x) {sum(x[["newN"]])}))
#inv.adult.consol$overall.inv.nom <- unlist(lapply(inv.adult.consol$data, function(x) {sum(x[["overall.inv.nom.perDS"]])}))
#inv.adult.consol <- inv.adult.consol %>% mutate(overall.inv = overall.inv.nom/overall.inv.denom) %>% unnest()

#plot.consol.inv(inv.adult.consol)
#### Relative Rates---

#### Relative rates data clean up ---

#relative.total.ds <- tot.ds %>% filter(Serotype == "TOTAL")
#relative.total.ds$Serotype <- NULL
#relative.total.ds$Serogroup <- NULL
#colnames(relative.total.ds)[3:4] <- c("Tot.carriage", "Tot.disease")

#rel.ds <- dplyr::left_join(tot.ds, relative.total.ds, by = c("agegrp", "DS"))
#colnames(rel.ds)[7:8] <- c("n.swab", "N")
#rel.ds <- dplyr::left_join(rel.ds, select(population.data, DS, time.int), by = "DS")
#rel.ds$Serotype <- unlist(rel.ds$Serotype)

#child.rel.ds <- rel.ds %>% filter(agegrp == "infant") %>% filter(!is.na(time.int)) %>% filter(Serotype != "TOTAL") # children's datasets
#adult.rel.ds <- rel.ds %>% filter(agegrp == "adult") %>% filter(!is.na(time.int)) %>% filter(Serotype != "TOTAL") # adults' dataset

# Consolidate all serotype data to get one invasiveness measure for each serotype
#something <- child.rel.ds %>% filter(!DS %in% c("France.post.PCV7", "Barcelona.post.PCV7", "France.post.PCV13", "Bogota.post.PCV7", "Massachusetts.post.PCV7",
#                                             "Navajo.post.PCV7", "Atlanta.post.PCV7")) %>% group_by(Serotype, Serogroup) %>% nest()
#carriage <- lapply(something$data, function(x) {sum(x$carriage)})
#disease <- lapply(something$data, function(x) {sum(x$disease)})
#n.swab <- lapply(something$data, function(x) {sum(x$n.swab)})
#N <- lapply(something$data, function(x) {sum(x$N)})
#something$carriage <- unlist(carriage)
#something$disease <- unlist(disease)
#something$n.swab <- unlist(n.swab)
#something$N <- unlist(N)
#something$time.int <- rep(1, nrow(something))
#something$agegrp <- rep("children", nrow(something))
#something$DS <- rep(NA, nrow(something))
#something$data <- NULL
#something.invasiveness <- data.frame(t(sapply(1:nrow(something), function(x) invasiveness(something[x,]))))
#ggplot(something.invasiveness, aes(x = reorder(unlist(Serotype), -(unlist(invasiveness))), y = unlist(invasiveness))) + 
#  geom_point() + 
#  geom_errorbar(aes(ymin = unlist(invasiveness.low), ymax = unlist(invasiveness.high)), width = 0.01) + 
#  labs(x = "Serotype", y = "Invasiveness", title = "Total Consolidated Invasiveness") + theme_bw() + theme_light()

# Invasiveness for children + graphs made into one PDF
#rel.child <- data.frame(t(sapply(1:nrow(child.rel.ds), function(x) invasiveness(child.rel.ds[x,]))))
#rel.child.2 <- rel.child %>% dplyr::select(DS, Serotype, invasiveness, invasiveness.low, invasiveness.high)
#rel.child.sero <- unlist(unique(rel.child.2$Serotype))
#relchildgraphs <- lapply(rel.child.sero, function(x) {plot.invasiveness(rel.child.2 %>% filter(Serotype == x))})
#mypath.relchild <- file.path("Q:","Technical","R","Case-to-carrier",paste("Relative_inv_children", ".pdf", sep=""))
#pdf(mypath.relchild, width = 20, height = 5)
#print(relchildgraphs)
#dev.off()

# Carriage prevalence for children + graphs made into one PDF
#rel.child.carr <- rel.child %>% dplyr::select(DS, Serotype, carr.prev, carr.prev.low, carr.prev.high)
#relchildgraphs.carr <- lapply(rel.child.sero, function(x) {plot.carr(rel.child.carr %>% filter(Serotype == x))})
#mypath.relchild.carr <- file.path("Q:","Technical","R","Case-to-carrier",paste("Relative_carr_children", ".pdf", sep=""))
#pdf(mypath.relchild.carr, width = 20, height = 5)
#print(relchildgraphs.carr)
#dev.off()

# Disease incidence for children + graphs made into one PDF
#rel.child.dis <- rel.child %>% dplyr::select(DS, Serotype, IPD.inc, IPD.inc.low, IPD.inc.high)
#relchildgraphs.dis <- lapply(rel.child.sero, function(x) {plot.dis(rel.child.dis %>% filter(Serotype == x))})
#mypath.relchild.dis <- file.path("Q:","Technical","R","Case-to-carrier",paste("Relative_dis_children", ".pdf", sep=""))
#pdf(mypath.relchild.dis, width = 20, height = 5)
#print(relchildgraphs.dis)
#dev.off()

# Invasiveness for adults + graphs made into one PDF
#rel.adult <- as.data.frame(t(sapply(1:nrow(adult.rel.ds), function(x) invasiveness(adult.rel.ds[x,]))))
#rel.adult.2 <- rel.adult %>% dplyr::select(DS, Serotype, invasiveness, invasiveness.low, invasiveness.high)
#rel.adult.sero <- unlist(unique(rel.adult.2$Serotype))
#reladultgraphs <- lapply(rel.adult.sero, function(x) {plot.invasiveness(rel.adult.2 %>% filter(Serotype == x))})
#mypath.reladult <- file.path("Q:","Technical","R","Case-to-carrier",paste("Relative_inv_adults", ".pdf", sep=""))
#pdf(mypath.reladult, width = 20, height = 5)
#print(reladultgraphs)
#dev.off()

# Carriage prevalence for adults + graphs made into one PDF
#rel.adult.carr <- rel.adult %>% dplyr::select(DS, Serotype, carr.prev, carr.prev.low, carr.prev.high)
#reladultgraphs.carr <- lapply(rel.adult.sero, function(x) {plot.carr(rel.adult.carr %>% filter(Serotype == x))})
#mypath.reladult.carr <- file.path("Q:","Technical","R","Case-to-carrier",paste("Relative_carr_adults", ".pdf", sep=""))
#pdf(mypath.reladult.carr, width = 20, height = 5)
#print(reladultgraphs.carr)
#dev.off()

# Disease incidence for adults + graphs made into one PDF
#rel.adult.dis <- rel.adult %>% dplyr::select(DS, Serotype, IPD.inc, IPD.inc.low, IPD.inc.high)
#reladultgraphs.dis <- lapply(rel.adult.sero, function(x) {plot.dis(rel.adult.dis %>% filter(Serotype == x))})
#mypath.reladult.dis <- file.path("Q:","Technical","R","Case-to-carrier",paste("Relative_dis_adults", ".pdf", sep=""))
#pdf(mypath.reladult.dis, width = 20, height = 5)
#print(reladultgraphs.dis)
#dev.off()

# Relative invasiveness of combined children + adults in one graph
#rel.child.inv.tot <- rel.child %>% dplyr::select(agegrp, DS, Serotype, invasiveness, invasiveness.low, invasiveness.high, IPD.inc, IPD.inc.low, IPD.inc.high)
#rel.adult.inv.tot <- rel.adult %>% dplyr::select(agegrp, DS, Serotype, invasiveness, invasiveness.low, invasiveness.high, IPD.inc, IPD.inc.low, IPD.inc.high)
#rel.inv.tot <- bind_rows(rel.child.inv.tot, rel.adult.inv.tot)
#rel.sero <- unique(c(rel.child.sero, rel.adult.sero))
#rel.tot.graphs <- lapply(rel.sero, function(x) {plot.inv.allages(rel.inv.tot %>% filter(Serotype == x))})
#mypath.rel.tot.inv <- file.path("Q:","Technical","R","Case-to-carrier",paste("Relative_inv_tot_JointFit", ".pdf", sep=""))
#pdf(mypath.rel.tot.inv, width = 20, height = 5)
#print(rel.tot.graphs)
#dev.off()

#rel.DS <- unique(c(unlist(unique(rel.adult$DS)), unlist(unique(rel.child$DS))))
#plot.inv.allages.allDS(rel.inv.tot)

# Relative IPD incidence of combined children + adults in one graph
#rel.child.dis.tot <- rel.child %>% dplyr::select(agegrp, DS, Serotype, IPD.inc, IPD.inc.low, IPD.inc.high)
#rel.adult.dis.tot <- rel.adult %>% dplyr::select(agegrp, DS, Serotype, IPD.inc, IPD.inc.low, IPD.inc.high)
#rel.dis.tot <- bind_rows(rel.child.dis.tot, rel.adult.dis.tot)
#rel.dis.tot.graphs <- lapply(rel.sero, function(x) {plot.dis.allages(rel.dis.tot %>% filter(Serotype == x))})
#mypath.rel.tot.dis <- file.path("Q:","Technical","R","Case-to-carrier",paste("Relative_dis_tot", ".pdf", sep=""))
#pdf(mypath.rel.tot.dis, width = 20, height = 5)
#print(rel.dis.tot.graphs)
#dev.off()

# Relative = invasiveness
#plot.inv.allages.allDS.rel <- function(sero) {
#  x.Sero <- unlist(sero$Serotype)
#  y.inv <- unlist(sero$invasiveness)
#  IPD.inc <- unlist(sero$IPD.inc)
#  sero$agegrp <- factor(sero$agegrp, levels = c("infant","adult"))
#  ggplot(sero, aes(x = reorder(x.Sero, -unlist(IPD.inc)), y = y.inv, shape = unlist(sero$agegrp), color = unlist(sero$DS))) + 
#    geom_point() + 
#    geom_errorbar(aes(ymin = unlist(sero$invasiveness.low), ymax = unlist(sero$invasiveness.high)), width = 0.01) +
#    labs(x = "Serotype", y = "Invasiveness", title = "Serotype-specific invasiveness", color = "Dataset", shape = "Age Group") +
#    scale_color_manual(values=c("#4e0e46","#5c5200","#3a0a5e","#cc283c", "#956100","#980060","#cd1971","#cb0f55","#6159cf","#aaa90a","#00954a","#00c79f","#2ce38d",
#                                "#86e8ab","#91da5b","#32b2ff","#ff9872","#ff8089","#ffc771",
#                                "#6c7fd7", "#573789","#466e20", "#ce9c2d", "#43c29e", "#36dee6")) +
    #scale_y_continuous(trans = 'log10', labels = scales::comma) +
#    theme_bw() + theme_light()
#}



#### Supplementary/unused code in data cleanup function that sits right before return function ---


  #### Plot case vs carrier raw numbers of each serotype (each point represents a dataset) [COMMENTED OUT] ---
  #case2carrier.plot <- function(index) {
  #  data <- nst.initial.data.2$data[index]
  #  c2c <- ggplot(as.data.frame(data), aes(x = carriage, y = disease, color = DS)) +
  #    geom_point() + theme_bw() +
  #    labs(x = "No. of carriage isolates", y = "No. of disease isolates", title = nst.initial.data.2$Serotype[index])
  #}
  
  #g <- lapply(seq(nst.initial.data.2$data), function(x) {case2carrier.plot(x)}) COMMENTED OUT
  
  #### make pdf files for each serogroup [COMMENTED OUT] ---
  
  #VT7_file <- file.path("Q:","Technical","R", "Case-to-carrier",paste(data.name, "VT7",".pdf", sep=""))
  #pdf(file = VT7_file, width = 15, height = 8)
  #multiplot(plotlist = g[which(nst.initial.data.2$Serogroup == "VT7")], cols = 2) #plot all VT7 in one graph
  #dev.off()
  
  #VT10_file <- file.path("Q:","Technical","R", "Case-to-carrier",paste(data.name, "VT10",".pdf", sep=""))
  #pdf(file = VT10_file, width = 15, height = 8)
  #multiplot(plotlist = g[which(nst.initial.data.2$Serogroup == "VT10")], cols = 2) #plot all VT10 in one graph
  #dev.off()
  
  #VT13_file <- file.path("Q:","Technical","R", "Case-to-carrier",paste(data.name, "VT13",".pdf", sep=""))
  #pdf(file = VT13_file, width = 15, height = 8)
  #multiplot(plotlist = g[which(nst.initial.data.2$Serogroup == "VT13")], cols = 2) #plot all VT13 in one graph
  #dev.off()
  
  #NVT_file <- file.path("Q:","Technical","R", "Case-to-carrier",paste(data.name, "NVT",".pdf", sep=""))
  #pdf(file = NVT_file, width = 10, height = 8)
  #print(g[which(nst.initial.data.2$Serogroup == "NVT")])
  #dev.off()
  
  #nst.initial.data.ds <- nst.initial.data.2 %>% unnest() %>% group_by(DS) %>% nest() # nest dataframe according to each dataset
  
  #### add a column with crude case:carrier ratio [COMMENTED OUT] ---
  
  #c2c.ratio <- function(data) {
  #  data <- data %>% dplyr::filter(carriage != 0 | disease != 0) %>%  # eliminate serotypes with no cases in disease or carriage
  #    dplyr::mutate(c2c.ratio = disease/carriage) # add a case to carrier ratio
  #}
  
  #nst.initial.data.ds$data <- lapply(nst.initial.data.ds$data, function(x) clean.data(x))
  
  #cols_serogroup <- c("VT7" = "#70AD47", "VT10" = "#4472C4", "VT13" = "#ED7D31", "NVT" = "black")
  
  #### plot each serotype's case:carrier ratio (one plot per dataset) [COMMENTED OUT] ---
  
  #ds.specific.plot <- function(index) {
  #  data <- nst.initial.data.ds$data[index]
  #  ds.spec <- ggplot(as.data.frame(data), aes(x = reorder(Serotype, -c2c.ratio), y = c2c.ratio, color = Serogroup)) +
  #    geom_point() + theme_bw() +
  #    labs(x = "Serotype", y = "Case-to-carrier ratio", title = nst.initial.data.ds$DS[index]) + 
  #    geom_hline(yintercept = 1, linetype = "dashed") + # horizontal line to indicate where cases > carriers (y > 1) vs cases < carriers (y < 1)
  #    scale_color_manual(breaks = c("VT7", "VT10", "VT13", "NVT"), name = "Serogroup", values = cols_serogroup)
  #}
  
  #d <- lapply(seq(nst.initial.data.ds$data), function(x) {ds.specific.plot(x)})
  
  #### make pdfs for plots [COMMENTED OUT]  ---
  
  #for (n in 1:length(nst.initial.data.ds$DS)){
  #  filename <- file.path("Q:","Technical","R", "Case-to-carrier",paste(data.name, nst.initial.data.ds$DS[n], ".pdf", sep=""))
  #  pdf(filename, width = 15, height = 7)
  #  print(d[n])
  #  dev.off()
  #}



#### Heatmaps and matrices - NOT USED/COMMENTED OUT ---


#clean.invasiveness.data <- function(dataset) { # function to make a matrix for the heatmap COMMENTED OUT
#  ds <- dataset %>% dplyr::select(DS, Serotype, invasiveness, invasiveness.low, invasiveness.high)
#  ds.spread <- ds %>% tidyr::spread(Serotype, invasiveness)
#  row.names(ds.spread) <- ds.spread$DS
#  ds.spread <- ds.spread[,-1]
#  ds.spread[ds.spread == "NULL"] <- NA
#  ds.spread[ds.spread == "Inf"] <- max(unlist(ds.spread[ds.spread != Inf | ds.spread != NA])) + 1
#  ds.matrix <- as.matrix(ds.spread)
#  ds.matrix <- apply(ds.matrix, 1, FUN = as.numeric)
#  rownames(ds.matrix) <- colnames(ds.spread)
#  return(ds.matrix)
#}

# Heatmap for children
#rel.child.matrix <- clean.invasiveness.data(rel.child)
#heatmap(rel.child.matrix, Rowv = NA, Colv = NA, scale = "row", col = colorRampPalette(brewer.pal(8, "Blues"))(25))

# Heatmap for adults
#rel.adult.matrix <- clean.invasiveness.data(rel.adult)
#heatmap(rel.adult.matrix, Rowv = NA, Colv = NA, scale = "row", col = colorRampPalette(brewer.pal(8, "Blues"))(25))

# Heatmap for children only
#abs.child.matrix <- clean.invasiveness.data(abs.child)
#heatmap(abs.child.matrix, Rowv = NA, Colv = NA, scale = "row", col = colorRampPalette(brewer.pal(8, "Reds"))(25))

# Heatmap for adult disease vs children carriage
#abs.adult.matrix <- clean.invasiveness.data(abs.adult)
#heatmap(abs.adult.matrix, Rowv = NA, Colv = NA, scale = "row", col = colorRampPalette(brewer.pal(8, "Blues"))(25))