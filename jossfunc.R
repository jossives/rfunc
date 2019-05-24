## Summarizes data.
## Gives count, mean, standard deviation, standard error of the mean, and confidence interval (default 95%).
##   data: a data frame.
##   measurevar: the name of a column that contains the variable to be summariezed
##   groupvars: a vector containing names of columns that contain grouping variables
##   na.rm: a boolean that indicates whether to ignore NA's
##   conf.interval: the percent range of the confidence interval (default is 95%)
summarySE <- function(data=NULL, measurevar, groupvars=NULL, na.rm=FALSE,
                      conf.interval=.95, .drop=TRUE) {
  require(plyr)
  
  # New version of length which can handle NA's: if na.rm==T, don't count them
  length2 <- function (x, na.rm=FALSE) {
    if (na.rm) sum(!is.na(x))
    else       length(x)
  }
  
  # This does the summary. For each group's data frame, return a vector with
  # N, mean, and sd
  datac <- ddply(data, groupvars, .drop=.drop,
                 .fun = function(xx, col) {
                   c(N    = length2(xx[[col]], na.rm=na.rm),
                     mean = mean   (xx[[col]], na.rm=na.rm),
                     median = median   (xx[[col]], na.rm=na.rm),
                     sd   = sd     (xx[[col]], na.rm=na.rm)
                   )
                 },
                 measurevar
  )


  
  datac$se <- datac$sd / sqrt(datac$N)  # Calculate standard error of the mean
  
  # Confidence interval multiplier for standard error
  # Calculate t-statistic for confidence interval: 
  # e.g., if conf.interval is .95, use .975 (above/below), and use df=N-1
  ciMult <- qt(conf.interval/2 + .5, datac$N-1)
  datac$ci <- datac$se * ciMult
  
  # Calculate binomial error, for use if data is binomial or binary
  datac$binomial.error <- sqrt(datac$mean * (1-datac$mean) / datac$N)
  
  # Rename the "mean" column    
  datac <- rename(datac, c("mean" = measurevar))
  
  return(datac)
}

## fisher.table
fisher.table <- function (cell.1.1, cell.1.2, cell.2.1, cell.2.2){
  ctable = table(c("control","treatment"),c(0,1))
  ctable[1,1] = cell.1.1
  ctable[1,2] = cell.1.2
  ctable[2,1] = cell.2.1
  ctable[2,2] = cell.2.2
  ctable          
  test.result <- fisher.test(ctable)
  return(test.result)
}

## fisher.2vector
# 2 vectors of 0s and 1s are passed
# * The first is control (top row)
# * The second is treatement (bottom row)
fisher.2vector <- function (top.row,bottom.row){
  ctable = table(c("control","treatment"),c(0,1))
  ctable[1,1] = length(which(top.row == 0)) 
  ctable[1,2] = length(which(top.row == 1)) 
  ctable[2,1] = length(which(bottom.row == 0))
  ctable[2,2] = length(which(bottom.row == 1)) 
  ctable          
  test.result <- fisher.test(ctable)
  return(test.result)
}

## report.2vector
# do a fisher exact test and report values in some sort of useful way
# At this point it is specific to this dataframe
report.2vector <- function (local.data,QNUM1,QNUM2){
    ftest.1.pre <- fisher.2vector(
        local.data$QCORR[local.data$PREPOST == "PRE" & local.data$QNUM==QNUM1 & local.data$TREATMENT == 0],
        local.data$QCORR[local.data$PREPOST == "PRE" & local.data$QNUM==QNUM1 & local.data$TREATMENT == 1]
    )

    ftest.1.post <- fisher.2vector(
        local.data$QCORR[local.data$PREPOST == "POST" & local.data$QNUM==QNUM1 & local.data$TREATMENT == 0],
        local.data$QCORR[local.data$PREPOST == "POST" & local.data$QNUM==QNUM1 & local.data$TREATMENT == 1]
    )

    ftest.2.pre <- fisher.2vector(
        local.data$QCORR[local.data$PREPOST == "PRE" & local.data$QNUM==QNUM2 & local.data$TREATMENT == 0],
        local.data$QCORR[local.data$PREPOST == "PRE" & local.data$QNUM==QNUM2 & local.data$TREATMENT == 1]
    )

    ftest.2.post <- fisher.2vector(
        local.data$QCORR[local.data$PREPOST == "POST" & local.data$QNUM==QNUM2 & local.data$TREATMENT == 0],
        local.data$QCORR[local.data$PREPOST == "POST" & local.data$QNUM==QNUM2 & local.data$TREATMENT == 1]
    )

    cat(ftest.1.pre$p.value, ftest.1.pre[[3]][[1]],"\n")
    cat(ftest.1.post$p.value, ftest.1.post[[3]][[1]],"\n")
    cat(ftest.2.pre$p.value, ftest.2.pre[[3]][[1]],"\n")
    cat(ftest.2.post$p.value, ftest.2.post[[3]][[1]],"\n")

}

## compact.in
compact.in <- function (local.data,QNUM1,QNUM2){

    mydata.pair <- subset(local.data, QNUM == QNUM1 | QNUM == QNUM2)

    qpair <- summarySE(mydata.pair, measurevar="QCORR",
                              groupvars=c("QNUM","TREATMENT","POSTFLAG","PREPOST"))
                              
    options(repr.plot.width=5, repr.plot.height=3)
    
    limits <- aes(ymax = QCORR + binomial.error, ymin = QCORR - binomial.error)

    ggplot(qpair, aes(x=PREPOST, y=QCORR, group=TREATMENT, colour=TREATMENT, shape=TREATMENT) ) + 
        geom_line() + 
        geom_errorbar(limits, width=0.25) +
        facet_wrap(~QNUM)  
   
}


## plot.crtsplit.byQ
plot.crtsplit.byQ <- function (local.data,QNUM1,QNUM2,local.PREPOST){

    mydata.pair <- subset(local.data, QNUM == QNUM1 | QNUM == QNUM2)
    mydata.crt <- subset(mydata.pair, CRTSPLIT == 0 | CRTSPLIT == 1)
    mydata.crtsplit <- subset(mydata.crt, PREPOST = local.PREPOST)

    qpair <- summarySE(mydata.crtsplit, measurevar="QCORR",
                              groupvars=c("QNUM","TREATMENT","CRTSPLIT"))
                              
    options(repr.plot.width=5, repr.plot.height=3)
    
    limits <- aes(ymax = QCORR + binomial.error, ymin = QCORR - binomial.error)

    ggplot(qpair, aes(x=CRTSPLIT, y=QCORR, group=TREATMENT, colour=TREATMENT, shape=TREATMENT) ) + 
        geom_line() + 
        geom_errorbar(limits, width=0.25) +
        facet_wrap(~QNUM)  
   
}

## plot.crtsplit.byQ
plot.crtsplit.byQ.alt <- function (local.data,QNUM1,QNUM2){

    mydata.pair <- subset(local.data, QNUM == QNUM1 | QNUM == QNUM2)
    mydata.crtsplit <- subset(mydata.pair, CRTSPLIT == 0 | CRTSPLIT == 1)

    qpair <- summarySE(mydata.crtsplit, measurevar="QCORR",
                              groupvars=c("QNUM","TREATMENT","CRTSPLIT","PREPOST"))
                              
    options(repr.plot.width=5, repr.plot.height=3)
    
    limits <- aes(ymax = QCORR + binomial.error, ymin = QCORR - binomial.error)

    options(repr.plot.width=8, repr.plot.height=2.5)

    ggplot(qpair, aes(x=CRTSPLIT, y=QCORR, group=TREATMENT, colour=TREATMENT, shape=TREATMENT) ) + 
        geom_line() + 
        geom_errorbar(limits, width=0.25) +
        facet_wrap(~ PREPOST * QNUM, ncol = 4 ,labeller = label_wrap_gen(multi_line=FALSE))  
        # cyl~am+vs, labeller = label_wrap_gen(multi_line=FALSE)
   
}

## 
cohens.d.from.odds.simple <- function (odds){
# d = log_odds * root(3) / pi
    return(log(odds)*sqrt(3)/pi)   

}

