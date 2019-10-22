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

infant <- read.csv("Q:/Technical/R/Case-to-carrier/pneumo_invasiveness_inf.csv")
adult <- read.csv("Q:/Technical/R/Case-to-carrier/pneumo_invasiveness_adu.csv")

inv_priors_df <- read.csv("Q:/Technical/R/Case-to-carrier/Weinberger_invpriors.csv")
inv_priors_df <- inv_priors_df %>% mutate(stddev = sqrt(1/log.inv.prec.age1))

#### Rename column names of each dataset -------------------------------------------------------------------------------------------------------------------
colnames(infant)[2:24] <- paste("carriage_",colnames(infant)[2:24], sep = "")
colnames(infant)[26:48] <- paste("disease_",colnames(infant)[26:48], sep = "")
colnames(infant)[26:48] <- substr(colnames(infant)[26:48],1, nchar(colnames(infant)[26:48])-2)

colnames(adult)[2:8] <- paste("carriage_",colnames(adult)[2:8], sep = "")
colnames(adult)[10:16] <- paste("disease_",colnames(adult)[10:16], sep = "")
colnames(adult)[10:16] <- substr(colnames(adult)[10:16],1, nchar(colnames(adult)[10:16])-2)

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

#### Make population dataframe with study time interval, n.swab, population numbers ------------------------------------------------------------------------
# sources of population numbers available upon request #

DS <- c("Alabama.pre.PCV", "Alaska.pre.PCV", "Barcelona.post.PCV7", "Bogota.post.PCV7", "Bogota.pre.PCV", "Caracas.pre.PCV", "Czech.pre.PCV","France.post.PCV7",
        "Goroka.pre.PCV", "Kenya.pre.PCV", "Navajo.post.PCV7", "Oxford.pre.PCV", "Sweden.pre.PCV", "E.W.pre.PCV", "France.post.PCV13", "Massachusetts.post.PCV7", 
        "Atlanta.post.PCV7", "Atlanta.pre.PCV", "Finland.pre.PCV", "Greece.pre.PCV", "Ontario.pre.PCV", "Iceland.pre.PCV", "Morocco.pre.PCV", "Portugal.pre.PCV",
        "Stockholm.pre.PCV")
time.int <- c(3, 4, 4, 1, 4, 1, 9, 6,
              6, 4/2, 2, 3, 8, 1, 6, 3,
              1, 1, 5, 5, 1, 9, 1, 2,
              1)
n.swab <- c(827, NA, 209, 246, 197, 1004, 425, 1212, 
            2844, NA, 6541, 501, 550, 3752, 1212, 2969,
            451, 231, 3290, 464, 1139, NA, 200, 769,
            611)
N.children <- c(19316, NA, 228000, 357200, 357200, 146125, 784863, 1592988,
                3957, NA, 30558, 3778636, 376838, 3091000, 1597519, 820000,
                298831, 204680, 318083, NA, 580507, NA, 212566, NA,
                NA)
N.adults <- c(232373, NA, NA, NA, NA, NA, NA, NA, 
              NA, NA, 108789, NA, 7066020, 48702414, NA, NA,
              NA, NA, NA, NA, NA, NA, NA, 10356117,
              2004152)

population.data <- data.frame(DS, n.swab, N.children, N.adults, time.int)

#### Absolute rates data clean up --------------------------------------------------------------------------------------------------------------------------

## Data cleaning for children carriage vs children disease

child_only <- population.data %>% filter(!is.na(N.children)) %>% dplyr::select(-4) # keep only popln data for which children have data on n.swabs and N.children
colnames(child_only)[3] <- "N" # change column name of imported DF so it's easier to append population values onto main dataframe

keepers <- paste(c(as.character(child_only$DS)), collapse = '|') # reformat dataset names that we are keeping 

child_df <- tot.ds %>% group_by(agegrp, DS) %>% nest() %>%
  filter(!grepl("adult", agegrp)) %>% # remove adult age group
  filter(grepl(keepers, DS)) %>% # keep only the keepers
  unnest()

child.abs.ds <- dplyr::left_join(child_df, child_only, by = "DS") %>% filter(Serotype != "TOTAL") # append dataframes together + remove total, working dataframe

## Data cleaning for children carriage vs adults disease

adult_only <- population.data %>% filter(!is.na(N.adults)) %>% dplyr::select(-3)
colnames(adult_only)[3] <- "N"

keepers.adult <- paste(c(as.character(adult_only$DS)), collapse = '|')

adult_df <- tot.ds %>% group_by(agegrp, DS) %>% nest() %>%
  filter(!grepl("children", agegrp)) %>%
  filter(grepl(keepers.adult, DS)) %>% 
  unnest

adult.abs.ds <- dplyr::left_join(adult_df, adult_only, by = "DS") %>% filter(Serotype != "TOTAL") # append dataframes together + remove total, working dataframe

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
  
  invas_dev <- runif(n = num.dev, min = 0.0000001, max = 0.5)
  
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
  
  #### bins & credible intervals ----
  bin.df <- data.frame(bins = bins[-length(bins)], distrib = dist.mean, DS = rep(DS, nrow(dist.mean)))
  invas_cp <- cumsum(bin.df$distrib/sum(bin.df$distrib))
  invas.low <- bin.df$bins[which(invas_cp >= 0.025)[1]] #sum(bin.df$bins[which(invas_cp >= 0.025)[1]],bin.df$bins[which(invas_cp >= 0.025)[2]])/2
  invas.high <- sum(bin.df$bins[which(invas_cp >= 0.975)[1]],bin.df$bins[which(invas_cp >= 0.975)[2]])/2
  
  # to plot vs weinberger distributions; and for model comparison
  bin.df.new <- data.frame(bins = invas_dev, distrib = posterior_samp) # new bin.df
  bin.df.new <- bin.df.new[order(bin.df.new$bins),]
  df.name <- paste(DS, "distrib", agegroup, serotype, sep = ".")
  assign(df.name, bin.df.new, envir = .GlobalEnv)
   
  
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

calc_invas <- function(df) { # function that returns a single invasiveness with credible intervals for a serotype
                             # given all the datasets that contain it
  invas <- 
        #seq(0.00000001, 0.001, # children pre and postPCV; adults
        #seq(0.00000001, 0.002, # children prePCV sero 14, 27;
                                # children postPCV sero 4, 8, 10B, 33F, 38, 24F, 18C, 18A, 3, 14
        #seq(0.00000001, 0.02, #0.01, # children prePCV sero 4, 7F, 18B, 13, 31, 12F; 
                               # adults sero 1, 4, 8, 20, 24F, 17F, 35A; 
                               # children postPCV sero 5, 7F, 12F, 25A, 12B
        seq(0.00000001, 0.1, #children prePCV sero 1, 5
        #seq(0.00000001, 0.04, # children postPCV sero 1
        length.out = 2000)
  likelihood <- sapply(invas, function(s) overall_inv_2(df, s))
  #plot(invas, LL, main = df$Serotype[1])
  
  Serotype <- as.character(unique(df$Serotype))
  
  # Weinberger prior
  #if (Serotype %in% inv_priors_df$st) {
  #  prior_invas <- dlnorm(invas, meanlog = inv_priors_df$log.inv.age1[which(inv_priors_df$st == Serotype)]-log(1000), 
  #                               sdlog = inv_priors_df$stddev[which(inv_priors_df$st == Serotype)])
  #} else {
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
  #likelihooddf <- data.frame(invas = invas, likelihood = likelihood, distribution = "Our results")
  #weinbergerdf <- data.frame(invas = invas, likelihood = prior_invas*max(likelihood)/max(prior_invas), distribution = "Weinberger et al")
  #consol.df <- bind_rows(likelihooddf, weinbergerdf)
  #ggplot(consol.df) + geom_line(aes(x = invas, y = likelihood, colour = distribution)) +
  #  coord_cartesian(xlim = c(0,0.005)) + ggtitle(Serotype) + theme_classic() + labs(x = "Invasiveness", y = "Probability Density") + 
  #  scale_y_continuous(labels = scales::scientific) + theme(legend.position = "none")
  
  
  bin.df <- data.frame(bins = invas, distrib = likelihood)
  df.name <- paste(paste("sero", df$Serotype[1], sep = ""), "distrib", df$agegrp[1], sep = ".")
  assign(df.name, bin.df, envir = .GlobalEnv)
  
  return(inv_row)
}

#### -------------------------------------------------------------------------------------------------------------------------------------------------------
#### CODE FOR PLOTS   --------------------------------------------------------------------------------------------------------------------------------------
#### -------------------------------------------------------------------------------------------------------------------------------------------------------
cols <- c("VT7" = "#70AD47", "VT10"= "#4472C4", "VT13" = "#ED7D31", "NVT" = "black")

# comment out either [geom_boxplot()] OR [geom_point() + geom_errorbar()] depending on what you need to see
plot.inv.allages.allDS <- function(sero) {
  VT7 <- c("4","6B","9V", "14", "18C", "19F", "23F")
  VT10 <- c("1","5","7F")
  VT13 <- c("3", "6A", "19A")
  NVT <- as.character(unique(sero$Serotype[!sero$Serotype %in% c(VT7, VT10, VT13)]))
  
  x.sero <- sero$DS
            #sero$Serotype
  y.inv <- sero$invasiveness
  sero$Serogroup <- factor(sero$Serogroup, levels = c("VT7", "VT10", "VT13", "NVT"))
  sero$agegrp <- factor(sero$agegrp, levels = c("children","adults"))
  # order x-axis DS by vaccination periods
  # comment out if x-axis is serotypes
  x.sero <- factor(x.sero, levels = c("Alabama.pre.PCV", "Atlanta.pre.PCV", "Bogota.pre.PCV","Caracas.pre.PCV","Czech.pre.PCV","E.W.pre.PCV","Finland.pre.PCV",
                                       "Goroka.pre.PCV","Morocco.pre.PCV","Ontario.pre.PCV","Oxford.pre.PCV","Portugal.pre.PCV","Stockholm.pre.PCV","Sweden.pre.PCV",
                                       "Atlanta.post.PCV7", "Barcelona.post.PCV7", "Bogota.post.PCV7","France.post.PCV7","Massachusetts.post.PCV7","Navajo.post.PCV7",
                                       "France.post.PCV13"))
  
  # order serotypes in x-axis from VT7, VT10, VT13 to NVT 
  # comment out if x-axis is datasets (DS)
  #x.sero <- factor(x.sero, levels = c(VT7, VT10, VT13, NVT))
  
  ggplot(sero, aes(x = x.sero, y = y.inv, group = sero$agegrp)) +
   geom_errorbar(aes(ymin = sero$invasiveness.low, ymax = sero$invasiveness.high, group = sero$agegrp), position = position_dodge(0.75),
                  width = 0.01) + 
    geom_point(aes(color = sero$Serogroup, shape = sero$agegrp, group = sero$agegrp), 
               size = 2, position = position_dodge(0.75)) +
    labs(x =  "Dataset", #"Serotype", #
         y = "Invasiveness", title = "Dataset-Specific Invasiveness", #"Serotype-specific invasiveness",#
         color = "Vaccine Category", 
         shape = "Age Group") +
    guides(color = guide_legend(order = 1), shape = guide_legend(order = 2)) + 
    scale_color_manual(values= cols) +
    scale_y_continuous(trans = 'log10') +
    theme_bw() + theme_light() +
    theme(axis.text.x = element_text(angle = 30, hjust = 1)) # DS angle = 30; serotype angle = 45
}

plot.carr.allages.allDS <- function(sero) { # plots carriage prevalence
  x.sero <- #sero$DS
    sero$Serotype
  y.carrprev <- sero$carr.prev
  sero$Serogroup <- factor(sero$Serogroup, levels = c("VT7", "VT10", "VT13", "NVT"))
  sero$vaccinationera <- factor(sero$vaccinationera, levels = c("pre.PCV", "post.PCV7", "post.PCV13"))

  ggplot(sero, aes(x = reorder(x.sero, -y.carrprev), y = y.carrprev)) + 
    geom_errorbar(aes(ymin = sero$carr.prev.low, ymax = sero$carr.prev.high), width = 0.01) +
    geom_point(aes(color = sero$Serogroup, shape = sero$vaccinationera)) + 
    labs(x = "Serotype", y = "Carriage Prevalence", title = "Carriage Prevalence", color = "Vaccine Category",
         shape = "Vaccination Era") +
    scale_color_manual(values=c("#70AD47", "#4472C4", "#ED7D31", "black")) + # PCV7 green, PCV10 blue, PCV13 orange, NVT black
    #scale_y_continuous(trans = 'log10', labels = scales::comma) +
    theme_bw() + theme_light() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) 
}

plot.dis.allages.allDS <- function(sero) { # plots all datasets' IPD incidence for all serotypes in a graph
  x.sero <- sero$DS # change to sero$DS for dataset plot
    #sero$Serotype
  
  x.Nt <- sero$N*sero$time.int
  
  # comment out if x-axis is Serotypes
  x.sero <- factor(x.sero, levels = c("Alabama.pre.PCV", "Atlanta.pre.PCV", "Bogota.pre.PCV","Caracas.pre.PCV","Czech.pre.PCV","E.W.pre.PCV","Finland.pre.PCV", 
                                        "Goroka.pre.PCV","Morocco.pre.PCV","Ontario.pre.PCV","Oxford.pre.PCV","Portugal.pre.PCV","Stockholm.pre.PCV","Sweden.pre.PCV",  
                                        "Atlanta.post.PCV7", "Barcelona.post.PCV7", "Bogota.post.PCV7","France.post.PCV7","Massachusetts.post.PCV7","Navajo.post.PCV7",
                                        "France.post.PCV13"))
  
  y.IPDinc <- (sero$lambda/x.Nt)*100000
  sero$agegrp <- factor(sero$agegrp, levels = c("children","adults"))
  sero$Serogroup <- factor(sero$Serogroup, levels = c("VT7", "VT10", "VT13", "NVT"))
  ggplot(sero, aes(x = x.sero, y = y.IPDinc)) + 
    geom_errorbar(aes(ymin = (sero$lambda.low/x.Nt)*100000, ymax = (sero$lambda.high/x.Nt)*100000,
                      group = sero$agegrp), width = 0.01, position = position_dodge(0.75)) +
    geom_point(aes(color = sero$Serogroup, shape = sero$agegrp, group = sero$agegrp), #size = 5, 
               position = position_dodge(0.75)) +
    labs(x = "Dataset", y = "IPD Incidence per 100,000 ppl", title = "IPD Incidence", color = "Vaccine Category", shape = "Age Group") +
    scale_color_manual(values= cols) + 
    scale_y_continuous(trans = 'log10', labels = scales::comma) +
    theme_bw() + theme_light() +
    theme(axis.text.x = element_text(angle = 30, hjust = 1))
}

plot.consol.inv <- function(sero) { # plots the overall/consolidated invasiveness of all serotypes in one graph
  x.sero <- sero$Serotype
  y.inv <- sero$overall.invasiveness
  sero$agegrp <- factor(sero$agegrp, levels = c("children","adults"))
  sero$Serogroup <- factor(sero$Serogroup, levels = c("VT7", "VT10", "VT13","NVT"))
  
  ggplot(sero, aes(x = reorder(x.sero, -y.inv), y = y.inv)) +
    geom_errorbar(aes(ymin = sero$overall.invasiveness.low, ymax = sero$overall.invasiveness.high, group = sero$agegrp), 
                  width = 0.01, position = position_dodge(0.75)) +
    geom_point(aes(color = sero$Serogroup, shape = sero$agegrp, group = sero$agegrp), 
               position = position_dodge(0.75), size = 2) +
    labs(x = "Serotype", y = "Invasiveness", title = "Overall Invasiveness", color = "Vaccine Category", shape = "Age Group") +
    scale_color_manual(values= cols) +
    scale_y_continuous(trans = 'log10') +
    theme_bw() + theme_light()
}

plot.DS.overall <- function(sero, consolidated.inv) { # plots one DS invasiveness with overall/consolidated invasiveness
  
  #### scatterplot w overall invasiveness vs DS-specific invasiveness
  new.sero <- sero %>% select(DS, Serogroup, Serotype, invasiveness, invasiveness.low, invasiveness.high)
  new.sero <- data.frame(lapply(new.sero, FUN = unlist))
  
  consolidated.inv.dat <- consolidated.inv %>% filter(Serotype %in% new.sero$Serotype)
  consolidated.inv.dat$DS <- rep("overall", nrow(consolidated.inv.dat))
  consolidated.inv.dat$agegrp <- NULL
  consolidated.inv.dat <- data.frame(lapply(consolidated.inv.dat, FUN = unlist))
  colnames(consolidated.inv.dat) <- c("Serotype", "invasiveness", "invasiveness.low", "invasiveness.high", "Serogroup", "DS")
  
  new.df <- dplyr::left_join(new.sero, consolidated.inv.dat, by = c("Serotype", "Serogroup"))
  new.df$Serogroup <- factor(new.df$Serogroup, levels = c("VT7", "VT10", "VT13", "NVT"))
  
  ggplot(data = new.df, aes(x = invasiveness.x, y = invasiveness.y, color = Serogroup)) + 
    geom_point() + # (size = 5.5) + if geom_text included 
    geom_abline(slope = 1, intercept = 0) +
    geom_errorbar(aes(ymin = invasiveness.low.y, ymax = invasiveness.high.y), width = 0.0001) +
    geom_errorbarh(aes(xmin = invasiveness.low.x, xmax = invasiveness.high.x), height = 0.001) +
    scale_y_continuous(trans = 'log10') +
    scale_x_continuous(trans = 'log10') +
    scale_color_manual(values= cols) +
    #geom_text(aes(label = Serotype), size = 3, color = "white", fontface = "bold") +
    labs(x = "Overall invasiveness", y = paste(sero$DS, "invasiveness", sep = " "), title = paste("Overall vs", sero$DS, "invasiveness", sep = " ")) +
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

#### -------------------------------------------------------------------------------------------------------------------------------------------------------
#### -------------------------------------------------------------------------------------------------------------------------------------------------------
#### Invasiveness Rates ------------------------------------------------------------------------------------------------------------------------------------
#### -------------------------------------------------------------------------------------------------------------------------------------------------------
#### -------------------------------------------------------------------------------------------------------------------------------------------------------
##### Sero and DS specific invasiveness ---------
# Invasiveness for children only
abs.child <- as.data.frame(t(sapply(1:nrow(child.abs.ds), function(x) invasiveness(child.abs.ds[x,]))))
abs.child <- data.frame(map(abs.child, unlist))

write.csv(abs.child, "abs.child-new.csv")

# Invasiveness for adults relative to children's carriage
abs.adult <- as.data.frame(t(sapply(1:nrow(adult.abs.ds), function(x) invasiveness(adult.abs.ds[x,]))))
abs.adult <- data.frame(map(abs.adult, unlist))

write.csv(abs.adult, "abs.adult-new.csv")

# Invasiveness for children and adults in same DF
abs.tot <- bind_rows(abs.child, abs.adult)
abs.tot$agegrp <- as.factor(abs.tot$agegrp)
abs.tot$DS <- as.factor(abs.tot$DS)
abs.tot$Serotype <- as.factor(abs.tot$Serotype)
abs.tot$vaccinationera <- as.factor(sub("^.*\\.p","p", abs.tot$DS))
abs.child <- abs.tot %>% filter(agegrp == "children")
abs.adult <- abs.tot %>% filter(agegrp == "adults")

# Absolute carriage prevalence of children for all DS in one graph
plot.carr.allages.allDS(abs.child)

# Absolute disease incidence for all DS 
plot.dis.allages.allDS(abs.tot)

# Absolute invasiveness for all DS

plot.inv.allages.allDS(abs.child)
plot.inv.allages.allDS(abs.adult)
plot.inv.allages.allDS(abs.tot) # pdf dim 4 x 16 for serotype ; pdf dim 4 x 20 for DS

##### Consolidated serotype-specific invasiveness ----
# Consolidated invasiveness for children in pre-PCV time periods
inv.child.consol <- data.frame(map(child.abs.ds, unlist)) %>% 
  filter(!DS %in% c("France.post.PCV7", "Barcelona.post.PCV7", "France.post.PCV13", "Bogota.post.PCV7", "Massachusetts.post.PCV7","Navajo.post.PCV7", 
                                                                             "Atlanta.post.PCV7")) %>% group_by(Serotype) %>% nest()
nDS <- sapply(inv.child.consol$data, function(x) nrow(x)) # number of datasets for each serotype
excl.sero <- inv.child.consol$Serotype[which(nDS == 1)] # serotypes with only 1 dataset to be excluded for analysis
inv.child.consol <- data.frame(inv.child.consol %>% unnest() %>% filter(!Serotype %in% excl.sero))
unique.child.sero <- unique(inv.child.consol$Serotype)
#unique.child.sero <- unique.child.sero[!unique.child.sero %in% c(sero.to.change.low, sero.to.change.high, sero.to.change.higher)]

consol.child2 <- lapply(unique.child.sero, function(x) {calc_invas(inv.child.consol %>% filter(Serotype == x))})
consol.child <- do.call("rbind", consol.child2)

sero.to.change.low <- c("14", "27")
consol.child.increasedbounds0 <- lapply(sero.to.change.low, function(x) {calc_invas(inv.child.consol %>% filter(Serotype == x))})
consol.child.increasedbounds0 <- do.call("rbind", consol.child.increasedbounds0)
sero.to.change.high <- c("4", "7F", "18B", "31", "13", "12F")
consol.child.increasedbounds1 <- lapply(sero.to.change.high, function(x) {calc_invas(inv.child.consol %>% filter(Serotype == x))})
consol.child.increasedbounds1 <- do.call("rbind", consol.child.increasedbounds1)
sero.to.change.higher <- c("1", "5")
consol.child.increasedbounds2 <- lapply(sero.to.change.higher, function(x) {calc_invas(inv.child.consol %>% filter(Serotype == x))})
consol.child.increasedbounds2 <- do.call("rbind", consol.child.increasedbounds2)

consol.child <- consol.child[!(consol.child$Serotype %in% c(sero.to.change.low, sero.to.change.high, sero.to.change.higher)),]
consol.child <- bind_rows(consol.child, consol.child.increasedbounds0, consol.child.increasedbounds1, consol.child.increasedbounds2)

colnames(consol.child) <- c("Serotype", "overall.invasiveness", "overall.invasiveness.low", "overall.invasiveness.high")

consol.child$Serogroup <- NA
consol.child$Serogroup[which(consol.child$Serotype %in% c("4", "6B", "9V", "14", "18C", "19F", "23F"))] <- "VT7"
consol.child$Serogroup[which(consol.child$Serotype %in% c("1", "5", "7F"))] <- "VT10"
consol.child$Serogroup[which(consol.child$Serotype %in% c("3", "6A", "19A"))] <- "VT13"
consol.child$Serogroup[is.na(consol.child$Serogroup)] <- "NVT"

plot.consol.inv(consol.child) # pdf dim 5 x 15
write.csv(consol.child, file = "consol.child-new.csv")

# Consolidated invasiveness for adults in pre-PCV time periods
inv.adult.consol <- data.frame(map(adult.abs.ds, unlist)) %>% 
  filter(!DS %in% c("Navajo.post.PCV7")) %>% group_by(Serotype) %>% nest()
nDS.adult <- sapply(inv.adult.consol$data, function(x) nrow(x)) # number of datasets for each serotype
excl.sero.adult <- inv.adult.consol$Serotype[which(nDS.adult == 1)] # serotypes with only 1 dataset to be excluded for analysis
inv.adult.consol <- data.frame(inv.adult.consol %>% unnest() %>% filter(!Serotype %in% excl.sero.adult))
unique.adult.sero <- unique(inv.adult.consol$Serotype)

consol.adult2 <- lapply(unique.adult.sero, function(x) {calc_invas(inv.adult.consol %>% filter(Serotype == x))})
consol.adult <- do.call("rbind", consol.adult2)

increased.bounds.serotypes <- c("1","20", "4", "8", "24F", "17F", "35A")
consol.adult.increasedbounds <- lapply(increased.bounds.serotypes, function(x) {calc_invas(inv.adult.consol %>% filter(Serotype == x))})
consol.adult.increasedbounds <- do.call("rbind", consol.adult.increasedbounds)

consol.adult <- consol.adult[!(consol.adult$Serotype %in% increased.bounds.serotypes),]
consol.adult <- bind_rows(consol.adult, consol.adult.increasedbounds)

colnames(consol.adult) <- c("Serotype", "overall.invasiveness", "overall.invasiveness.low", "overall.invasiveness.high")

consol.adult$Serogroup <- NA
consol.adult$Serogroup[which(consol.adult$Serotype %in% c("4", "6B", "9V", "14", "18C", "19F", "23F"))] <- "VT7"
consol.adult$Serogroup[which(consol.adult$Serotype %in% c("1", "5", "7F"))] <- "VT10"
consol.adult$Serogroup[which(consol.adult$Serotype %in% c("3", "6A", "19A"))] <- "VT13"
consol.adult$Serogroup[is.na(consol.adult$Serogroup)] <- "NVT"

plot.consol.inv(consol.adult) # pdf dim 5 x 15
write.csv(consol.adult, file = "consol.adult-new.csv")

# Consolidated invasiveness for both age groups in one graph
consol.child$agegrp <- rep("children", nrow(consol.child))
consol.adult$agegrp <- rep("adults", nrow(consol.adult))
inv.tot.consol <- rbind(consol.child, consol.adult)
plot.consol.inv(inv.tot.consol) # pdf dim 5 x 15

#### Post-vaccination consolidated ----- SUPPLEMENTARY MATERIAL
inv.child.consol.post <- data.frame(map(child.abs.ds, unlist)) %>% 
  filter(DS %in% c("France.post.PCV7", "Barcelona.post.PCV7", "France.post.PCV13", "Bogota.post.PCV7", "Massachusetts.post.PCV7","Navajo.post.PCV7", 
                    "Atlanta.post.PCV7")) %>% group_by(Serotype) %>% nest()
nDS.post <- sapply(inv.child.consol.post$data, function(x) nrow(x)) # number of datasets for each serotype
excl.sero.post <- inv.child.consol.post$Serotype[which(nDS.post == 1)] # serotypes with only 1 dataset to be excluded for analysis
inv.child.consol.post <- data.frame(inv.child.consol.post %>% unnest() %>% filter(!Serotype %in% excl.sero.post))
unique.child.sero.post <- unique(inv.child.consol.post$Serotype)

consol.child2.post <- lapply(unique.child.sero.post, function(x) {calc_invas(inv.child.consol.post %>% filter(Serotype == x))})
consol.child.post <- do.call("rbind", consol.child2.post)

increased.bounds.serotypes.childpost1 <- c("5", "7F", "12F", "25A", "12B")
consol.child.post.increasedbounds1 <- lapply(increased.bounds.serotypes.childpost1, function(x) {calc_invas(inv.child.consol.post %>% filter(Serotype == x))})
consol.child.post.increasedbounds1 <- do.call("rbind", consol.child.post.increasedbounds1)

Sero1.childpost <- calc_invas(inv.child.consol.post %>% filter(Serotype == "1"))

increased.bounds.serotypes.childpost2 <- c("4", "8", "10B", "33F", "38", "24F", "18C", "18A", "3", "14")
consol.child.post.increasedbounds2 <- lapply(increased.bounds.serotypes.childpost2, function(x) {calc_invas(inv.child.consol.post %>% filter(Serotype == x))})
consol.child.post.increasedbounds2 <- do.call("rbind", consol.child.post.increasedbounds2)

consol.child.post <- consol.child.post[!(consol.child.post$Serotype %in% c(increased.bounds.serotypes.childpost1, "1", increased.bounds.serotypes.childpost2)),]
consol.child.post <- bind_rows(consol.child.post, Sero1.childpost, consol.child.post.increasedbounds1, consol.child.post.increasedbounds2)

colnames(consol.child.post) <- c("Serotype", "overall.invasiveness", "overall.invasiveness.low", "overall.invasiveness.high")

consol.child.post$Serogroup <- NA
consol.child.post$Serogroup[which(consol.child.post$Serotype %in% c("4", "6B", "9V", "14", "18C", "19F", "23F"))] <- "VT7"
consol.child.post$Serogroup[which(consol.child.post$Serotype %in% c("1", "5", "7F"))] <- "VT10"
consol.child.post$Serogroup[which(consol.child.post$Serotype %in% c("3", "6A", "19A"))] <- "VT13"
consol.child.post$Serogroup[is.na(consol.child.post$Serogroup)] <- "NVT"

plot.consol.inv(consol.child.post) # pdf dim 5 x 15
write.csv(consol.child.post, "consol.child.postPCV-priors.csv")

#### Debugging: Weinberger's prior vs likelihood and posterior  --------------------------------------------------------------------------------------------

plot.serolikelihood <- function(df, serotype, agegroup) {
  serotype.df <- df %>% filter(Serotype == serotype)
  dfs <- unique(serotype.df$DS)
  eachdf1 <- lapply(dfs, function(x) get(paste(x, "distrib", agegroup, serotype, sep = ".")))
  eachdf2 <- bind_rows(eachdf1)
  eachdf.max <- max(eachdf2$distrib)
  eachdf.norm <- lapply(eachdf1, function(x) x$distrib*eachdf.max/max(x$distrib))
  
  for (i in 1:length(eachdf1)) {
    eachdf1[[i]]$distrib <- eachdf.norm[[i]]
  }
  
  eachdf <- bind_rows(eachdf2)
  
  sero.invas.dev <- rlnorm(n = 560000, meanlog = inv_priors_df$log.inv.age1[which(inv_priors_df$st == serotype)]-log(1000), 
                           sdlog = inv_priors_df$stddev[which(inv_priors_df$st == serotype)])
  sero.prior <- dlnorm(sero.invas.dev, meanlog = inv_priors_df$log.inv.age1[which(inv_priors_df$st == serotype)]-log(1000), 
                         sdlog = inv_priors_df$stddev[which(inv_priors_df$st == serotype)])
  sero.priordf <- data.frame(bins = sero.invas.dev, distrib = sero.prior*max(eachdf$distrib)/max(sero.prior), DS = "prior")
  eachdf <- rbind(eachdf, sero.priordf)
  eachdf$col <- "Our posterior"
  eachdf$col[which(eachdf$DS == "prior")] <- "Weinberger et al"
  
  ggplot() + 
    geom_line(eachdf, mapping = aes(x = bins, y = distrib, colour = DS, alpha = factor(col))) + scale_alpha_discrete(range = c(0.3, 1)) +
    scale_y_continuous(labels = scales::scientific) + coord_cartesian(xlim = c(0,0.01)) + ggtitle(paste("Serotype", serotype, df$agegrp, "invasiveness", sep = " ")) + 
    theme_minimal() + theme(legend.position = "none")
}

## children single invasiveness##

serotypes.with.priors <- unique(child.abs.ds$Serotype)[unique(child.abs.ds$Serotype) %in% inv_priors_df$st]
eachsero <- lapply(serotypes.with.priors, function(x) plot.serolikelihood(df = abs.child, serotype = x, agegroup = "children"))
likelihoodvsprior1 <- grid.arrange(eachsero[[1]], eachsero[[2]], eachsero[[3]], eachsero[[4]], eachsero[[5]], eachsero[[6]], eachsero[[7]], eachsero[[8]], 
                                  eachsero[[9]], eachsero[[10]], eachsero[[11]], eachsero[[12]], eachsero[[13]], eachsero[[14]], eachsero[[15]], eachsero[[16]],
                                  #
                                  # eachsero[[17]], eachsero[[18]], eachsero[[19]], eachsero[[20]], eachsero[[21]], eachsero[[22]], eachsero[[23]], eachsero[[24]],
                                  # eachsero[[25]], eachsero[[26]], eachsero[[27]], eachsero[[28]], eachsero[[29]], eachsero[[30]], eachsero[[31]], eachsero[[32]],
                                  # 
                                  # eachsero[[33]], eachsero[[34]], eachsero[[35]], eachsero[[36]], eachsero[[37]], eachsero[[38]], eachsero[[39]], eachsero[[40]],
                                  # eachsero[[41]], eachsero[[42]], eachsero[[43]], eachsero[[44]], eachsero[[45]], eachsero[[46]], eachsero[[47]], eachsero[[48]],

                                  # eachsero[[49]], eachsero[[50]], eachsero[[51]], eachsero[[52]], eachsero[[53]], eachsero[[54]], eachsero[[55]], eachsero[[56]],
                                  # eachsero[[57]], eachsero[[58]], eachsero[[59]],
                                  nrow = 4) #pdf dim 8 x 18

## adults single invasiveness##

serotypes.with.priors.adu <- unique(adult.abs.ds$Serotype)[unique(adult.abs.ds$Serotype) %in% inv_priors_df$st]
eachsero.adu <- lapply(serotypes.with.priors.adu, function(x) plot.serolikelihood(df = abs.adult, serotype = x, agegroup = "adults"))
likelihoodvsprior1.adu <- grid.arrange(#eachsero.adu[[1]], eachsero.adu[[2]], eachsero.adu[[3]], eachsero.adu[[4]], eachsero.adu[[5]], eachsero.adu[[6]], 
                                       # eachsero.adu[[7]], eachsero.adu[[8]], eachsero.adu[[9]], eachsero.adu[[10]], eachsero.adu[[11]], eachsero.adu[[12]], 
                                       # eachsero.adu[[13]], eachsero.adu[[14]], eachsero.adu[[15]], 
                                       #
                                       # eachsero.adu[[16]], eachsero.adu[[17]], eachsero.adu[[18]], eachsero.adu[[19]], eachsero.adu[[20]], eachsero.adu[[21]],
                                       # eachsero.adu[[22]], eachsero.adu[[23]], eachsero.adu[[24]], eachsero.adu[[25]], eachsero.adu[[26]], eachsero.adu[[27]],
                                       # eachsero.adu[[28]],

                                       # eachsero.adu[[29]], eachsero.adu[[30]], eachsero.adu[[31]], eachsero.adu[[32]], eachsero.adu[[33]], eachsero.adu[[34]],
                                       # eachsero.adu[[35]], eachsero.adu[[36]], eachsero.adu[[37]], eachsero.adu[[38]], eachsero.adu[[39]], eachsero.adu[[40]],
                                       # 
                                       eachsero.adu[[41]], eachsero.adu[[42]], eachsero.adu[[43]], eachsero.adu[[44]], eachsero.adu[[45]], eachsero.adu[[46]],
                                       eachsero.adu[[47]], eachsero.adu[[48]], eachsero.adu[[49]],
                                       nrow = 4) #pdf dim 8 x 18


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

#### Comparison each location invasiveness w overall -------------------------------------------------------------------------------------------------------

unique.DS.comparison.child <- unique(inv.child.consol$DS)
prePCV.DS.vs.overall.child <- lapply(unique.DS.comparison.child, function(x) 
  {plot.DS.overall(sero = abs.child %>% filter(DS == x), consol.child)})

mypath.childcomparison.inv <- file.path("Q:","Technical","R","Case-to-carrier","Figures for paper", paste("comp inv overall", ".pdf", sep=""))
pdf(mypath.childcomparison.inv, width = 7, height = 5)
print(prePCV.DS.vs.overall.child)
dev.off()

plot.posteriors <- function(df, serotype, agegroup) {
  serotype.df <- df %>% filter(Serotype == serotype)
  dfs <- unique(serotype.df$DS)
  eachdf1 <- lapply(dfs, function(x) get(paste(x, "distrib", agegroup, serotype, sep = ".")))
  eachdf2 <- do.call(rbind, eachdf1) #bind_rows(eachdf1)
  eachdf.max <- max(eachdf2$distrib)
  eachdf.norm <- lapply(eachdf1, function(x) x$distrib*eachdf.max/max(x$distrib))
  
  for (i in 1:length(eachdf1)) {
    eachdf1[[i]]$distrib <- eachdf.norm[[i]]
  }
  
  eachdf <- bind_rows(eachdf2)
  consol.distrib <- get(paste(paste("sero", serotype, sep = ""), "distrib", agegroup, sep = "."))
  consol.distrib$DS <- "consol"
  consol.distrib$distrib <- consol.distrib$distrib*eachdf.max/max(consol.distrib$distrib)
  eachdf <- rbind(eachdf, consol.distrib)
  eachdf$col <- "single"
  eachdf$col[which(eachdf$DS == "consol")] <- "consol"

  ggplot() + 
    geom_line(eachdf, mapping = aes(x = bins, y = distrib, colour = DS, alpha = factor(col))) + scale_alpha_discrete(range = c(1, 0.3)) +
    scale_y_continuous(labels = scales::scientific) + coord_cartesian(xlim = c(0,max(consol.distrib$bins)+0.01)) + 
    ggtitle(paste("Serotype", serotype, df$agegrp, "invasiveness", sep = " ")) + 
    theme_minimal() + theme(legend.position = "none")
}

consolvssingle.plots <- lapply(unique.child.sero, function(x) {plot.posteriors(df = inv.child.consol, serotype = x, agegroup = "children")})
consol.vs.single <- grid.arrange(#consolvssingle.plots[[1]], consolvssingle.plots[[2]], consolvssingle.plots[[3]], consolvssingle.plots[[4]], consolvssingle.plots[[5]],
                                 #consolvssingle.plots[[6]], consolvssingle.plots[[7]], consolvssingle.plots[[8]], consolvssingle.plots[[9]], consolvssingle.plots[[10]],
                                 
                                 # consolvssingle.plots[[11]], consolvssingle.plots[[12]], consolvssingle.plots[[13]], consolvssingle.plots[[14]], consolvssingle.plots[[15]],
                                 # consolvssingle.plots[[16]], consolvssingle.plots[[17]], consolvssingle.plots[[18]], consolvssingle.plots[[19]], consolvssingle.plots[[20]],

                                 # consolvssingle.plots[[21]], consolvssingle.plots[[22]], consolvssingle.plots[[23]], consolvssingle.plots[[24]], consolvssingle.plots[[25]],
                                 # consolvssingle.plots[[26]], consolvssingle.plots[[27]], consolvssingle.plots[[28]], consolvssingle.plots[[29]], consolvssingle.plots[[30]],
                                 # 
                                 consolvssingle.plots[[31]], consolvssingle.plots[[32]], consolvssingle.plots[[33]], consolvssingle.plots[[34]], consolvssingle.plots[[35]],
                                 consolvssingle.plots[[36]], consolvssingle.plots[[37]], consolvssingle.plots[[38]], consolvssingle.plots[[39]],
                                 
                                 nrow = 4) # pdf dim 7 x 11

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

consol_ev_allDS <- consol_evidence(inv.child.consol %>% filter(Serotype == "4"))
consol_ev <- mean(consol_ev_allDS)

consol_EW_4 <- consol_evidence(inv.child.consol %>% filter(Serotype == "4") %>% filter(DS == "E.W.pre.PCV"))
consol_ev_EW4 <- mean(consol_EW_4)

calc_invas <- function(df) { # function that returns a single invasiveness with credible intervals for a serotype
  # given all the datasets that contain it
  invas <- 
    #seq(0.00000001, 0.001, # children pre and postPCV; adults
    #seq(0.00000001, 0.002, # children prePCV sero 14, 27;
    # children postPCV sero 4, 8, 10B, 33F, 38, 24F, 18C, 18A, 3, 14
    seq(0.00000001, 0.02, #0.01, # children prePCV sero 4, 7F, 18B, 13, 31, 12F; 
        # adults sero 1, 4, 8, 20, 24F, 17F, 35A; 
        # children postPCV sero 5, 7F, 12F, 25A, 12B
        #seq(0.00000001, 0.1, #children prePCV sero 1, 5
        #seq(0.00000001, 0.04, # children postPCV sero 1
        length.out = 2000)
  likelihood <- sapply(invas, function(s) overall_inv_2(df, s))
  #plot(invas, LL, main = df$Serotype[1])
  
  Serotype <- as.character(unique(df$Serotype))
  
  # Weinberger prior
  #if (Serotype %in% inv_priors_df$st) {
  #  prior_invas <- dlnorm(invas, meanlog = inv_priors_df$log.inv.age1[which(inv_priors_df$st == Serotype)]-log(1000), 
  #                               sdlog = inv_priors_df$stddev[which(inv_priors_df$st == Serotype)])
  #} else {
  prior_invas <- dunif(invas, min = 0.00001, max = 0.5)
}

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
  
  # AUC.consol <- AUC(consol.distrib$bins, consol.distrib$distrib, method = "trapezoid")
  # AUC.single.each <- sapply(eachdf1, function(x) AUC(x$bins, x$distrib, method = "trapezoid"))
  # AUC.single.all <- prod(AUC.single.each)
  
  AUC.consol <- sum(consol.distrib$distrib)/nrow(consol.distrib)
  AUC.single.each <- sapply(eachdf1, function(x) sum(x$distrib)/nrow(x))
  AUC.single.all <- prod(AUC.single.each)
  
  bayesfac <- AUC.consol/AUC.single.all
  #bayesfac <- log(AUC.consol)-log(AUC.single)
  #bayesfac <- exp(bayesfac)
  bayesfac <- log(bayesfac)
  return(bayesfac)
}

newbayesfun <- function(serotype, df, DS, agegroup) { # function to compare evidence of consol.inv for one DS vs single inv for one DS
  calc_invas(df %>% filter(Serotype == serotype) %>% filter(DS == DS))
  eachdf1 <- get(paste(DS, "distrib", agegroup, serotype, sep = "."))
  consol.distrib <- get(paste(paste("sero", serotype, sep = ""), "distrib", agegroup, sep = "."))
  consol.distrib$DS <- "consol"
  AUC.single <- sum(eachdf1$distrib)/nrow(eachdf1)
  AUC.consol <- sum(consol.distrib$distrib)/nrow(consol.distrib)
  bayesfac <- AUC.consol/AUC.single
  return(bayesfac)
}

sero15BC.distrib.children <- `sero15B/C.distrib.children`
bayes <- lapply(unique.child.sero, function(x) bayesfxr(serotype = x, df = inv.child.consol, agegroup = "children"))
bayes.df <- data.frame(serotype = unique.child.sero, bayes.fxr = unlist(bayes), agegroup = "children", nDS = nDS[nDS > 1])
seroDS.df <- data.frame(serotype = unique(inv.child.consol$Serotype), nDS = nDS[nDS > 1])
seroDS.df <- seroDS.df[order(seroDS.df$nDS),]

sero15BC.distrib.adults <- `sero15B/C.distrib.adults`
bayes.adults <- lapply(unique.adult.sero, function(x) bayesfxr(serotype = x, df = inv.adult.consol, agegroup = "adults")) 
bayes.adults.df <- data.frame(serotype = unique.adult.sero, bayes.fxr = unlist(bayes.adults), agegroup = "adults", nDS = nDS.adult[nDS.adult > 1])
seroDS.adults.df <- data.frame(serotype = unique(inv.adult.consol$Serotype), nDS = nDS.adult[nDS.adult > 1])
seroDS.adults.df <- seroDS.adults.df[order(seroDS.adults.df$nDS),]

bayesfxr.all <- rbind(melt(bayesfxr.df), melt(bayesfxr.df.adults))
bayesfxr.all <- bayesfxr.all[!is.na(bayesfxr.all$value),]
colnames(bayesfxr.all) <- c("agegroup", "serotype", "bayes.fxr")
write.csv(bayesfxr.all, "bayesfactor.csv")

# bayes.conditions <- function(bf){
# if (bf == 1) {"No evidence"
# } else if (bf > 1 && bf < 3) {"Consolidated: anecdotal evidence"
# } else if (bf > 3 && bf < 10) {"Consolidated: moderate evidence"
# } else if (bf > 10 && bf < 30) {"Consolidated: strong evidence"
# } else if (bf > 30 && bf < 100) {"Consolidated: very strong evidence"
# } else if (bf > 100) {"Consolidated: moderate evidence"
# } else if (bf > (1/3) && bf < 1) {"Single: anecdotal evidence"
# } else if (bf > (1/10) && bf < (1/3)) {"Single: moderate evidence"
# } else if (bf > (1/30) && bf < (1/10)) {"Single: strong evidence"
# } else if (bf > (1/100) && bf < (1/30)) {"Single: very strong evidence"
# } else if (bf < (1/100)) {"Single: extreme evidence"
# }
#   return(label)
# }
# 
# bayesfxr.all$label <- sapply(bayesfxr.all$bayes.fxr, function(x) bayes.conditions(x))

#### testing bayes fxr
testdata <- child.abs.ds %>% filter(Serotype == "6A") %>% filter(DS == "Bogota.pre.PCV")
test <- invasiveness(testdata)
BF.df <- bayesfxr.df[(order(bayesfxr.df$invas_dev)),]
sumlike <- (sum(BF.df$posterior_samp))/nrow(BF.df) # 1.481264e-05
AUC.test <- AUC(invas_dev, posterior_samp, method = "trapezoid") # 7.180451e-06

testdata2 <- child.abs.ds %>% filter(Serotype == "6B") %>% filter(DS == "Bogota.pre.PCV")
test2 <- invasiveness(testdata2)
BF.df <- bayesfxr.df[(order(bayesfxr.df$invas_dev)),]
sumlike <- (sum(BF.df$posterior_samp))/nrow(BF.df) # 1.691704e-05
AUC.test <- AUC(BF.df$invas_dev, BF.df$posterior_samp, method = "trapezoid") # 1.024761e-05

testdata3 <- child.abs.ds %>% filter(Serotype == "6A") %>% filter(DS == "Atlanta.pre.PCV")
test3 <- invasiveness(testdata3)
BF.df <- bayesfxr.df[(order(bayesfxr.df$invas_dev)),]
sumlike <- (sum(BF.df$posterior_samp))/nrow(BF.df) # 0.0001593778
AUC.test <- AUC(BF.df$invas_dev, BF.df$posterior_samp, method = "trapezoid") # 6.824021e-05

testdata4 <- child.abs.ds %>% filter(Serotype == "6A") %>% filter(DS == "Atlanta.pre.PCV")
test4 <- invasiveness(testdata4)
BF.df <- bayesfxr.df[(order(bayesfxr.df$invas_dev)),]
sumlike <- (sum(BF.df$posterior_samp))/nrow(BF.df) # 0.0001593778
AUC.test <- AUC(BF.df$invas_dev, BF.df$posterior_samp, method = "trapezoid") # 6.824021e-05

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
  geom_point(aes(x = no.strains, y = overall.invasiveness, color = Serogroup, group = Serogroup), size = 4) + 
  geom_text(aes(x = no.strains, y = overall.invasiveness, label = Serotype), size = 2, color = "white", fontface = "bold") +
  scale_color_manual(values= cols) +
  labs(x = "Number of strains", y = "Overall Invasiveness")
  #stat_ellipse(aes(x = no.strains, y = overall.invasiveness, color = Serogroup, group = Serogroup, type = "norm"))

# adult age group
new.adult.comb <- dplyr::left_join(consol.gpsc, consol.adult) %>% drop_na()
new.adult.comb$Serogroup <- factor(new.adult.comb$Serogroup, levels = c("VT7", "VT10", "VT13","NVT")) 

ggplot(new.adult.comb) +# pdf dim 5 x 10
  geom_point(aes(x = no.strains, y = overall.invasiveness, color = Serogroup, group = Serogroup), size = 4) + 
  geom_text(aes(x = no.strains, y = overall.invasiveness, label = Serotype), size = 2.2, color = "white", fontface = "bold") +
  scale_color_manual(values= cols) +
  labs(x = "Number of strains", y = "Overall Invasiveness")
#stat_ellipse(aes(x = no.strains, y = overall.invasiveness, color = Serogroup, group = Serogroup, type = "norm"))

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