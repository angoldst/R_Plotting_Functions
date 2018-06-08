

oneway_collapse <- function(splitCol, meanCol, data){
  #helper function for generating the barplots 
  value <- as.data.frame(as.list(aggregate(data[, meanCol] ~ data[, splitCol], FUN = function(x) c(mean(x), lowsd(x), highsd(x)), simplify=T)))
  
  colnames(value) <- c('Splitby', 'Mean',  'LowSD', 'HighSD')
  value[, 'Splitby'] <- as.factor(value[, 'Splitby'])
  return(value)
}

print_multiPlot_1w_onewayanova_outlier <- function(withinVar, curOutputs, data, ncol=3, scales="free_y", covs='', thresh = 3, plotTitle = '', log_values =F){
  options(contrasts=c('contr.sum','contr.poly'))
  newDF <- data.frame(output = character(0),   withinVar = numeric(0), Mean = numeric(0), 'LowSD' = numeric(0), 'HighSD' = numeric(0), 'curText' = character(0), 'colCol' = character(0))
  
  for (i in seq_along(curOutputs)){
    var <- curOutputs[i]
    #curCol = getColor(var, curOutputs, colors[1:length(curOutputs)])
    
    lowSD_thresh <- mean(data[, var], na.rm=T) - thresh*sd(data[, var], na.rm=T)
    highSD_thresh <- mean(data[, var], na.rm=T) + thresh*sd(data[, var], na.rm=T)
    
    tempDF <- cbind(output = var,  oneway_collapse(withinVar, var, data[data[, var]<highSD_thresh & data[, var]>lowSD_thresh, ]))
    curAnova <- Anova(lmer(eval(parse(text=paste0(var, ' ~ ',  withinVar, ' + (1|subNum) ' , covs))), data[data[, var]<highSD_thresh & data[, var]>lowSD_thresh, ]), type=3, test.statistic = 'F')
    
    if(log_values == T){
      curAnova <- Anova(lmer(eval(parse(text=paste0('log(', var, ') ~ ',  withinVar, ' + (1|subNum) ' , covs))), data[data[, var]<highSD_thresh & data[, var]>lowSD_thresh, ]), type=3, test.statistic = 'F')
      
    }else{
      curAnova <- Anova(lmer(eval(parse(text=paste0(var, ' ~ ',  withinVar, ' + (1|subNum) ' , covs))), data[data[, var]<highSD_thresh & data[, var]>lowSD_thresh, ]), type=3, test.statistic = 'F')
      
    }
    sig <-  ifelse(curAnova[withinVar, 'Pr(>F)']<=0.05, ifelse(curAnova[withinVar, 'Pr(>F)']<=0.01, ifelse(curAnova[withinVar, 'Pr(>F)']<=0.001, '***', '**'), '*'), '')
    tempDF$curText <- paste0('p=', round(curAnova[withinVar, 'Pr(>F)'], 3), sig)
    newDF <- rbind(newDF, tempDF)
  }  
  
  #newDF <- newDF[order(newDF$colCol), ] #Only use if you're doing paired plots e.g. left and right regions
  longDF <- melt(data, id.vars = colnames(data)[!(colnames(data) %in% curOutputs)], curOutputs)
  newDF <- merge.data.frame(newDF, longDF, by.x=c('output', 'Splitby'), by.y=c('variable', withinVar))
  
  colors <- colorRampPalette(c('red', 'yellow'))(length(curOutputs))
  colorList <- ddply(newDF, .(output), summarize, color = getColor(output, curOutputs = curOutputs, colors = colors))
  limits <- aes( ymax = LowSD, ymin=HighSD)
  ggplot(newDF, aes_string(x='Splitby', y='Mean', fill='output', group='Splitby')) + 
    geom_bar(position="dodge", stat="identity", color='black', aes_string(alpha='Splitby')) + ylab('') + 
    geom_errorbar( position=position_dodge(.9), limits, width=0.25, color='black') + 
    geom_point(position=position_jitterdodge(),  aes_string(y='value', fill='output'),shape=21, color = 'black') +
    #geom_line(aes_string(y='value',)) + 
    facet_wrap(~output , scales = scales, ncol=ncol) + 
    
    geom_text(aes(label=curText) , data=newDF[newDF[, 'Splitby'] ==levels(newDF[, 'Splitby'])[1] & newDF[, 'subNum'] == newDF$subNum[1], ], x=-Inf, y=-Inf, hjust=-.02, vjust=-.3, size=4) +
    scale_fill_manual(values=colorList$color, guide="none") +
    scale_alpha_discrete(range = c(1, .3, .6)) + 
    theme_bw() + 
    theme(legend.position="none") + 
    ggtitle(plotTitle)
  
}

getColor <- function(x, curOutputs,  colors ){
  colors <- colors[1:length(curOutputs)]
  if(any(grepl(paste0('^', x, '$'), curOutputs))){
    color <- colors[grepl(paste0('^', x, '$'), curOutputs)]
  }else{
    
    color <- 'white'}
  color
}