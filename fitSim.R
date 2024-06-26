# Fit all sims in a given folder
fitSim <- function( estFolder="mod1-uncor",
                    simFolder="sim_N1",
                    nSim=100,
                    nParallelCores=8 )
{

  suppressWarnings(dir.create(paste("simOutputs",simFolder,estFolder,sep="/")))

  mods <- 1:nSim

  cl <- makeCluster(nParallelCores)

  parSapply( cl=cl, X=mods, FUN=fitIndividualSim, simFolder=simFolder, estFolder=estFolder )

  stopCluster(cl)

#  processSimEst( ctrl=ctrl,
#                 simFolder=simFolder,
#                 estFolder=estFolder )

}

fitIndividualSim <- function( i, simFolder, estFolder )
{
  source("initRR.R")
  ctlFile <- paste(estFolder,"/estControlFile.txt",sep="")
  load(paste("simOutputs/",simFolder,"/sim",i,".Rdata",sep=""))
  rpt <- fitRR( ctlFile=ctlFile, folder=".", simData=sim$simData, saveRun=FALSE )
  save(rpt,file=paste("simOutputs/",simFolder,"/",estFolder,"/fit",i,".Rdata",sep=""))
}

processSimEst <- function( ctrl=NULL,
                           simFolder="sim_N1",
                           estFolder="mod1-uncor",
                           nSim=100 )
{
  if( is.null(ctrl) )
  {
    controlTable  <- .readParFile( "estControlFile.txt" )
    ctrl <- .createList( controlTable )
  }

  converge <- numeric(nSim)

  # Count number of sims in folder
  simDir <- paste('simOutputs/',simFolder,sep="")
  estDir <- paste('simOutputs/',simFolder,"/",estFolder,sep="")
  simFiles <- list.files(simDir, pattern='sim')
  fitFiles <- list.files(estDir, pattern='fit')

  if(length(fitFiles) != length(simFiles))
    cat('WARNING: # sims =',length(simFiles),
                 '| # fits =',length(fitFiles),'\n')

  # get fit nums
  fitNums <- grep(fitFiles,pattern='fit')
  for( i in fitNums )
  {
    if( is.na(i) ) browser()
    load(paste(simDir,"/sim",i,".Rdata",sep="")) # obj name = 'sim'
    load(paste(estDir,"/fit",i,".Rdata",sep="")) # obj name = 'rpt'
    cvg <- substr(rpt$opt$message,1,4)
    if(length(cvg)>0)
      converge[i] <- cvg %in% c("sing","rela")
    #converge[i] <- !rpt$opt$convergence 
    #In TMB, 0 = convergence, 1=non-convergence

  }

  # Select only reps that converged for summary stats
  fitNums <- fitNums[converge==1]

  # Create objects to store sub-stock information for operating model (i.e. 'sim') and estimates (i.e. 'fit') for sub-stock run sizes
  stocks <- ctrl$stocks
  yrs <- ctrl$initYear:ctrl$lastYear
  nReps <- length(fitNums) # i index in array
  nStocks <- length(stocks)
  nYrs    <- length(yrs)

  om    <-list()
  om$mu_ist    <- array(NA, dim=c(nReps,nStocks,nYrs))
  om$N_ist     <- array(NA, dim=c(nReps,nStocks,nYrs))
  om$filename  <- vector(length=nReps)

  est   <- list()
  est$mu_ist   <- array(NA, dim=c(nReps,nStocks,nYrs))
  est$N_ist    <- array(NA, dim=c(nReps,nStocks,nYrs))
  est$filename <- vector(length=nReps)
  est$resids_its <- array(NA, dim=c(nReps,nYrs,nStocks))
  est$relErr_its <- array(NA, dim=c(nReps,nYrs,nStocks))
  est$muRelErr_its <- array(NA, dim=c(nReps,nYrs,nStocks))

  stats <- list()
  stats$RMSE_is   <- array(NA, dim=c(nReps,nStocks))
  stats$MRE_is    <- array(NA, dim=c(nReps,nStocks))
  stats$muMRE_is  <- array(NA, dim=c(nReps,nStocks))
  stats$CV_is     <- array(NA, dim=c(nReps,nStocks))
  stats$stocks    <- stocks

  # Now pull out stats for runs that converged
  for( i in 1:length(fitNums) )
  {
   if( is.na(i) ) browser() 
    # populate objects from simulation
    load(paste(simDir,"/sim",fitNums[i],".Rdata",sep="")) # obj name = 'sim'
    
    omMu_st <- sim$om$mu_st
    omN_st <- apply(sim$om$N_sdt, c(1,3), sum, na.rm=T)
    om$mu_ist[i,,] <- omMu_st
    om$N_ist[i,,] <- omN_st
    om$filename <- paste("sim",fitNums[i],".Rdata",sep="")

    # populate objects from simulation
    load(paste(estDir,"/fit",fitNums[i],".Rdata",sep="")) # obj name = 'rpt'

    estMu_st <- rpt$mu_st
    estN_st <- apply(rpt$N_dst, c(2,3), sum, na.rm=T)
    est$mu_ist[i,,] <- estMu_st
    est$N_ist[i,,] <- estN_st
    est$filename <- paste("fit",fitNums[i],".Rdata",sep="")

    # Calculate residuals for each stock
    for (s in 1:nStocks)
    {
      muResids <- omMu_st[s,]-estMu_st[s,]
      muRelErr <- muResids/omMu_st[s,]

      resids <- omN_st[s,]-estN_st[s,]
      relErr <- resids/omN_st[s,]

      est$resids_its[i,,s] <- resids
      est$relErr_its[i,,s] <- relErr
      est$muRelErr_its[i,,s] <- muRelErr
      
      # Root-mean squared error (i.e. accuracy)
      MSE <- mean(resids^2)
      stats$RMSE_is[i,s] <- sqrt(MSE)

      # Median relative error (i.e. bias)
      stats$MRE_is[i,s] <- median(relErr)
      stats$muMRE_is[i,s] <- median(muRelErr)

      # CV of residuals (i.e. precision)
      se <- sd(resids)/sqrt(length(resids))
      stats$CV_is[i,s] <- sd(resids)
      #stats$CV_is[i,s] <- se/mean(abs(resids))

    }  

  }  

  perf <- list(reps  = fitNums,
               om    = om,
               est   = est,
               stats = stats)

  filename <- file.path(estDir,'perf.Rdata')
  save(perf, file=filename)
  
  # Generate Plots
  #plotResids(perf, simDir=simDir)


}

plotResids <- function(obj, simDir)
{

  # Calculate means & quantiles across stocks
  stocks <- obj$stats$stocks
  RMSE_is <- obj$stats$RMSE_is
  meanRMSE_s <- apply(RMSE_is,2,mean)
  upRMSE_s  <- apply(RMSE_is,2,quantile,0.90)
  lowRMSE_s  <- apply(RMSE_is,2,quantile,0.1)

  CV_is <- obj$stats$CV_is
  meanCV_s <- apply(CV_is,2,mean)
  upCV_s  <- apply(CV_is,2,quantile,0.90)
  lowCV_s  <- apply(CV_is,2,quantile,0.1)

  MRE_is <- obj$stats$MRE_is
  meanMRE_s <- apply(MRE_is,2,mean)
  upMRE_s  <- apply(MRE_is,2,quantile,0.90)
  lowMRE_s  <- apply(MRE_is,2,quantile,0.1)

  clrs <- brewer.pal(length(stocks),'Dark2')

  # PLOTS

  filename <- file.path(simDir,'perfStats.pdf')
  pdf(file=filename)

  # 3 Panel Plot
  par(mfrow=c(3,1), mgp=c(2.5,0.6,0),
      mar=c(2,4,0.2,0), tck=-0.03)
  
  # Plot 1 - RMSE
  YLIM <- c(min(lowRMSE_s), max(upRMSE_s))
  plot(x=1:length(stocks), y=meanRMSE_s,las=1,
         type='n',xlab='', ylab='RMSE', ylim=YLIM)
  segments(x0=1:length(stocks), lwd=1.5,
           y0=lowRMSE_s, y1=upRMSE_s,
           col=clrs)
  points(x=1:length(stocks), y=meanRMSE_s,
           pch=1, col=clrs, cex=1, lwd=1.5)

  # Plot 2 - CV
  YLIM <- c(min(lowCV_s), max(upCV_s))
  plot(x=1:length(stocks), y=meanCV_s,las=1,
         type='n',xlab='', ylab='CV', ylim=YLIM)
  segments(x0=1:length(stocks), lwd=1.5,
           y0=lowCV_s, y1=upCV_s,
           col=clrs)
  points(x=1:length(stocks), y=meanCV_s,
           pch=1, col=clrs, cex=1, lwd=1.5)

  # Plot 3 - MRE
  YLIM <- c(min(lowMRE_s), max(upMRE_s))
  plot(x=1:length(stocks), y=meanMRE_s,las=1,
         type='n',xlab='', ylab='MRE', ylim=YLIM)
  abline(h=0,lty=3)
  segments(x0=1:length(stocks), lwd=1.5,
           y0=lowMRE_s, y1=upMRE_s,
           col=clrs)
  points(x=1:length(stocks), y=meanMRE_s,
           pch=1, col=clrs, cex=1, lwd=1.5)

  dev.off()

}  



