library(ggplot2)
library(reshape2)
library(plyr)
library(car)
library(afex)
source('./helper_functions.R')

#A series of functions to run and plot linear mixed models
#The general structure is a script (with prefix print_) that melts the data into long dataformat for "curOutcomes" runs (using the get_ helper scripts) and plots the statistical models. 

get_aov_car_1w_results <- function(df, withinVar, covs = '', idvar = 'subNum', dv = 'value'){
  #linear mixed model with 1 within subject factor
  #df = dataframe
  #covs = list of the column names of the covariates to include in model 
  #withinVar = the variable that contains repeates (i.e. time, occasion, session etc.)
  
  #besure to set sum contrasts before running
  #linear model generation
  
  mylmer <- aov_car(eval(parse(text=paste0(dv, ' ~ ', ifelse(all(covs== ''), '', paste0(covs, sep=' + ', collapse='')), withinVar, ' + Error(', idvar,'/', withinVar, ')', collapse=''))), df)
  model_text <- paste0(withinVar, ' p= ', round(mylmer$anova_table[withinVar, 'Pr(>F)'], 3))
  data.frame(model_text = model_text)
  
}


get_ezAnova_1w_results <- function(df, withinVar, covs = '', idvar = 'subNum', dv = 'value'){
  #linear mixed model with 1 within subject factor
  #df = dataframe
  #covs = list of the column names of the covariates to include in model 
  #withinVar = the variable that contains repeates (i.e. time, occasion, session etc.)
  
  #besure to set sum contrasts before running
  #linear model generation
  mylmer <- eval(parse(text=paste0('ezANOVA(data=df, dv=', dv, ', wid=', idvar, ', within=', withinVar, ', type=3)')))
  model_text <- paste0(withinVar, ' p= ', round(mylmer$ANOVA['p'], 3))
  data.frame(model_text = model_text)
  
}


get_lmer_1w_results <- function(df, withinVar, covs = '', idvar = 'subNum'){
  #linear mixed model with 1 within subject factor
  #df = dataframe
  #covs = list of the column names of the covariates to include in model 
  #withinVar = the variable that contains repeates (i.e. time, occasion, session etc.)
  
  #besure to set sum contrasts before running
  #linear model generation
  mylmer <- mixed(eval(parse(text=paste0('value ~ ', ifelse(all(covs== ''), '', paste0(covs, sep=' + ', collapse='')), withinVar, ' + (1|', idvar, ')', collapse=''))), df)
  model_text <- paste0(withinVar, ' p= ', round(mylmer$anova_table[withinVar, 'Pr(>F)'], 3))
  data.frame(model_text = model_text)
  
}

get_lmer_1w1b_results <- function(df, withinVar, betweenVar, covs = '', idvar = 'subNum', y ='value'){
  #linear mixed model with 1 within subject factor
  #df = dataframe
  #withinVar = the variable that contains repeates (i.e. time, occasion, session etc.)
  #betweenVar = the variable that contains the between subject factor
  #covs = list of the column names of the covariates to include in model 
  
  #besure to set sum contrasts before running
  #linear model generation
  #Of the form y ~ withinVar * betweenVar + covs + (1|idvar)
  mylmer <- mixed(eval(parse(text=paste0(y, ' ~ ', ifelse(all(covs== ''), '', paste0(covs, sep=' + ', collapse='')), withinVar, '*', betweenVar, ' + (1|', idvar, ')', collapse=''))), df)
  
  model_text <- paste0('MEw= ', round(mylmer$anova_table[withinVar, 'Pr(>F)'], 3), ' MEb= ', round(mylmer$anova_table[betweenVar, 'Pr(>F)'], 3), ' Int= ', round(mylmer$anova_table[paste0(withinVar, ':', betweenVar), 'Pr(>F)'], 3))
  data.frame(model_text = model_text)
}


print_multiPlot_1w_oneway_lmer <- function(withinVar, curOutputs, data, ncol=3, scales="free_y", covs='', outlierThresh = 10, plotTitle = '', log_values =F, idvar = 'subNum'){
  set_sum_contrasts() #need this in order to properly use anova function 
  
  if(!idvar %in% colnames(data)) stop(paste0('idvar not found in dataframe'))
  
  #Collapse into long data frame to be able to run everything in parallel
  if(covs == ""){
    longDF <- melt(data[, c(idvar, withinVar, curOutputs)], id.vars = c(idvar, withinVar))
  }else{
    longDF <- melt(data[, c(idvar, withinVar, curOutputs, covs)], id.vars = c(idvar, withinVar, covs) )
  }
  
  #Determine outlier 
  longDF <- ddply(longDF, c(withinVar, 'variable'), here(mutate), lowSD_thresh = low_outlier_thresh(value, thresh  = outlierThresh), highSD_thresh = high_outlier_thresh(value, thresh = outlierThresh), outlier = ifelse(value < lowSD_thresh | value > highSD_thresh, T, F))
  
  #Make plotting variables
  summaryStats <- ddply(longDF[!longDF$outlier & !is.na(longDF$outlier), ], c('variable', withinVar), summarise, 
                        Mean = mean(value), 
                        lowSEM = low_sem(value), 
                        highSEM = high_sem(value))
  
  #get lmer method stats stats
  lmer.results <-  ddply(longDF[!longDF$outlier & !is.na(longDF$outlier), ], .(variable),  function(x) get_lmer_1w_results(x, withinVar, covs = covs, idvar = idvar)) 
  
  summaryStats <- merge.data.frame(summaryStats, lmer.results, by = 'variable', all.x = T)
  
  #graph line plot
  limits <- aes( ymax = lowSEM, ymin=highSEM)
  ggplot(summaryStats, aes_string(x=withinVar, y='Mean', group = 'variable')) + 
    geom_point(position = position_dodge(.1)) + 
    geom_line(position = position_dodge(.1)) + 
    geom_errorbar( position=position_dodge(.9), limits, width=0.25, color='black') + 
    facet_wrap(~variable , scales = scales, ncol=ncol) + 
    
    geom_text(aes(label=model_text) , data=summaryStats[summaryStats[, withinVar] ==levels(summaryStats[, withinVar])[1], ], x=-Inf, y=-Inf, hjust=-.02, vjust=-.3, size=4) +
    theme_bw() + 
    theme(legend.position="none") + 
    ggtitle(plotTitle)
  
}

print_multiPlot_1w1b_twoway_lmer <- function(withinVar, betweenVar, curOutputs, data, ncol=3, scales="free_y", covs='', outlierThresh = 10, plotTitle = '', log_values =F, idvar = 'subNum'){
  set_sum_contrasts() #need this in order to properly use anova function 
  
  if(!idvar %in% colnames(data)) stop(paste0('idvar not found in dataframe'))
  
  #Collapse into long data frame to be able to run everything in parallel
  if(covs == ""){
    longDF <- melt(data[, c(idvar, withinVar, betweenVar, curOutputs)], id.vars = c(idvar, withinVar, betweenVar))
  }else{
    longDF <- melt(data[, c(idvar, withinVar, betweenVar, curOutputs, covs)], id.vars = c(idvar, withinVar, betweenVar, covs) )
  }
  
  #Determine outlier 
  longDF <- ddply(longDF, c(withinVar, betweenVar, 'variable'), here(mutate), 
                  lowSD_thresh = low_outlier_thresh(value, thresh =  outlierThresh), 
                  highSD_thresh = high_outlier_thresh(value, thresh = outlierThresh), 
                  outlier = ifelse(value < lowSD_thresh | value > highSD_thresh, T, F))
  
  #Make plotting variables
  summaryStats <- ddply(longDF[!longDF$outlier & !is.na(longDF$outlier), ], c('variable', withinVar, betweenVar), summarise, 
                        Mean = mean(value), 
                        lowSEM = low_sem(value), 
                        highSEM = high_sem(value))
  
  #get lmer method stats stats
  lmer.results <-  ddply(longDF[!longDF$outlier & !is.na(longDF$outlier), ], .(variable),  function(x) get_lmer_1w1b_results(x, withinVar, betweenVar, covs = covs, idvar = idvar))
  
  summaryStats <- merge.data.frame(summaryStats, lmer.results, by = 'variable', all.x = T)
  
  limits <- aes( ymax = lowSEM, ymin=highSEM)
  
  #graph line plot
  ggplot(summaryStats, aes_string(x=withinVar, y='Mean', group = betweenVar, color = betweenVar)) + 
    geom_point(position = position_dodge(.1)) + 
    geom_line(position = position_dodge(.1)) + 
    geom_errorbar( position=position_dodge(.1), limits, width=0.25) + 
    facet_wrap(~variable , scales = scales, ncol=ncol) + 
    geom_text(aes(label=model_text) , data=summaryStats[summaryStats[, withinVar] ==levels(summaryStats[, withinVar])[1], ], x=-Inf, y=-Inf, hjust=-.02, vjust=-.3, size=4) +
    scale_color_manual(values=c('darkred', 'darkblue', 'gold')) +
    theme_bw() + 
    #theme(legend.position="none") + 
    ggtitle(plotTitle)
  
}
