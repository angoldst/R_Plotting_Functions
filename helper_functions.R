#Identify outlier by a threshold of standard deviations
low_outlier_thresh <- function(x, thresh = 3){mean(x, na.rm=T) - thresh * sd(x, na.rm=T)}
high_outlier_thresh <- function(x, thresh = 3){mean(x, na.rm=T) + thresh * sd(x, na.rm=T)}

#calculate standard error of the mean and +/- mean for plotting
sem <- function(x){sqrt(var(x,na.rm=TRUE)/length(na.omit(x)))}
low_sem <- function(x){return(mean(x, na.rm=T)-sem(x))}
high_sem <- function(x){return(mean(x, na.rm=T)+sem(x))}
