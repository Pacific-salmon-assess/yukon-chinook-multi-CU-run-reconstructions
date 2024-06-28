# Fit run reconstruction
fitRR <- function( ctlFile="estControlFile.txt",
                   arrivSD=NULL,
                   folder="fittedMod",
                   simData=NULL,
                   saveRun=TRUE,
                   runRW=1 )
{
  # Initialize
  source("initRR.R")

  # Read in control file
  controlTable  <- .readParFile( ctlFile )
  # Create control list
  ctrl <- .createList( controlTable )
  init  <- ctrl$inits

  # Load abundance indices and stock composition data
  load("data/chinookYkData.Rdata")

  # Dimensions
  years  <- ctrl$initYear:ctrl$lastYear
  days   <- chinookYkData$days
  gears  <- chinookYkData$gears
  stocks <- ctrl$stocks
  nT     <- length(years)
  nD     <- length(days)
  nG     <- length(gears)
  nS     <- length(stocks)

  # Create directory for saving plots and report
  suppressWarnings(dir.create(folder))

  # DATA -------------------------------------------------------------------- #

  # Observed numbers by stock, day, year and gear
  n0_sdtg <- chinookYkData$n_pdtg[ , ,as.character(years), ]
  n_sdtg <- n0_sdtg
  n_sdtg[ , , ,1] <- n0_sdtg[ , , ,2]
  n_sdtg[ , , ,2] <- n0_sdtg[ , , ,1]
  # Abundance indices by day, year and gear
  E0_dtg  <- chinookYkData$E_dtg[ ,as.character(years), ]
  E_dtg <- E0_dtg
  E_dtg[ , ,1] <- E0_dtg[ , ,2]
  E_dtg[ , ,2] <- E0_dtg[ , ,1]
  # Mark-recapture run size indices for fish wheel years only (<2005)
  borderDat <- read.csv("data/border-passage.csv") %>%
               filter( year %in% years )
  I_t <- borderDat$mean
  mrCV_t <- borderDat$cv
  # Don't fit to mark-recapture index in sonar years
  I_t[years>=2005] <- NA
  mrCV_t[years>=2005] <- NA
  

  # Set NAs before first obs and after last obs to 0
  for( t in 1:nT )
  {
    for( g in 1:nG )
    {
      E <- E_dtg[ ,t,g]
      if( sum(!is.na(E)) > 0 )
      {
        E_dtg[is.na(E),t,g] <- 0
      }
    }
  }

  # Create TMB data object
  data <- list( n_sdtg    = n_sdtg,
                E_dtg     = E_dtg,
                I_t       = I_t,
                day_d     = days,
                mrCV_t    = mrCV_t,
                runRW     = runRW,    # Run size random walk switch (0=off, 1=on)
                runSD     = 1 )       # Run size random walk std dev

  # Simulated data
  if( !is.null(simData) )
    data <- simData

  # PARAMETERS -------------------------------------------------------------- #

  # Initial conditions
  runSize_st <- matrix( data=init$runSize_s, nrow=nS, ncol=nT )
  initMu_s   <- init$arrivMu_s
  sigma_s    <- init$arrivSD_s
  arrivErr_st <- matrix( data=0, nrow=nS, ncol=nT-1 )
  qE_sg      <- matrix( data=1, nrow=nS, ncol=nG )
  qE_sg[ ,2] <- 0.05
  qI_s       <- rep(1,nS)
  errSD_s    <- init$errSD_s
  obsErrSD_g <- rep(0.1,nG)
  cor_ss     <- matrix( data=0, nrow=nS, ncol=nS )
  diag(cor_ss) <- 1 # Correlation matrix must have 1 on diagonal

  lnDisp_tg  <- matrix( data=log(1e-4), nrow=nT, ncol=nG )
  #lnDisp_tg[4,2]  <- log(1)
  #lnDisp_tg[22,2] <- log(0.5)
  #lnDisp_tg[23,2] <- log(0.5)
  lnDisp_tg[ ,2] <- log(0.02)
  lnDisp_tg[17,2] <- log(2)
  #lnDisp_tg[5,2] <- log(1e-5)
  #lnDisp_tg[21,2] <- log(1e-5)
  #lnDisp_tg[22,2] <- log(0.5)
  #lnDisp_tg[23,2] <- log(0.5)

  # Create TMB parameter object
  pars <- list( lnRunSize_st = log(runSize_st),
                lnArrivMu_s  = log(initMu_s),
                #lnArrivSD    = log(mean(sigma_s)),
                lnArrivSD_s  = log(sigma_s),
                arrivErr_st  = arrivErr_st,
                #lnErrSD      = log(mean(errSD_s)),
                lnErrSD_s    = log(errSD_s),
                logitCor_ss  = logit(cor_ss,lb=-1,ub=1),
                lnqE_sg      = log(qE_sg),
                lnqI_s       = log(qI_s),
                lnDisp_tg    = lnDisp_tg )

  if(!is.null(arrivSD))
    pars$lnArrivSD_s <- rep(log(arrivSD),nS)

  # MAP --------------------------------------------------------------------- #

  qEmap_sg <- matrix( data=1:(nS*nG), nrow=nS, ncol=nG )
  qEmap_sg[ ,1] <- NA
  qEmap_sg[ ,2] <- ctrl$map$qFishWheel_s

  corMap_ss <- matrix( data=1:(nS*nS), nrow=nS, ncol=nS )
  # Always fix diagonal
  diag(corMap_ss) <- NA
  # Optional - define map for lower triangle
  # Let's fix the lower triangle (except diag) at a single value
  for( s in 2:nS )
  {
    if( ctrl$map$corType=="uncor" )
      corMap_ss[s,1:(s-1)] <- NA
    else if( ctrl$map$corType=="single" )
      corMap_ss[s,1:(s-1)] <- nS^2+1
  }
  # Mirror upper & lower triangle
  corMap_ss <- mirrorMatrix(corMap_ss)

  dispMap <- NA*lnDisp_tg

  mapArrivErr_st <- matrix( data=1:(nS*(nT-1)), nrow=nS, ncol=nT-1 )
  mapArrivErr_st[ ,which(apply(n_sdtg,3,sum,na.rm=TRUE)==0)-1] <- NA

  map <- list( lnArrivMu_s = as.factor(ctrl$map$arrivMu_s),
               lnArrivSD_s = as.factor(ctrl$map$arrivSD_s),
               arrivErr_st = as.factor(mapArrivErr_st),
               lnErrSD_s   = as.factor(ctrl$map$errSD_s),
               logitCor_ss = as.factor(corMap_ss),
               lnqE_sg     = as.factor(qEmap_sg),
               lnqI_s      = as.factor(NA*ctrl$map$qI_s),
               lnDisp_tg   = as.factor(dispMap) )

  # BUILD AND OPTIMIZE OBJECTIVE FUNCTION ----------------------------------- #

  # Load the DLL
  if("yukonChinookRunRecon" %in% names(getLoadedDLLs()))
    dyn.unload(dynlib("yukonChinookRunRecon"))          # unlink the C++ code if already linked
  # Compile the model
  compile("yukonChinookRunRecon.cpp", flags = "-g")
  # Load the DLL
  dyn.load(dynlib("yukonChinookRunRecon"))          # Dynamically link the C++ code

  # Build objective function
  obj <- MakeADFun( data       = data,
                    parameters = pars,
                    map        = map,
                    DLL        = "yukonChinookRunRecon",
                    random     = NULL )

  # Set bounds
  low <- obj$par*0-Inf
  upp <- obj$par*0+Inf

  # Optimization controlsarrivSD
  optCtrl <- list(  eval.max = ctrl$maxFunEval, 
                    iter.max = ctrl$maxIterations )

  # Optimize
  #sink("sink.txt")
  opt <- try( nlminb( start     = obj$par,
                      objective = obj$fn,
                      gradient  = obj$gr,
                      lower     = low,
                      upper     = upp,
                      control   = optCtrl ) )
  
  #sink()

  rptFE <- obj$report()

  if( ctrl$randEffects )
  {
    # Build objective function
    obj <- MakeADFun( data       = data,
                      parameters = rptFE[names(pars)],
                      map        = map,
                      DLL        = "yukonChinookRunRecon",
                      random     = "arrivErr_st" )

    # Optimize
    sink("sink.txt")
    opt <- try( nlminb( start     = obj$par,
                        objective = obj$fn,
                        gradient  = obj$gr,
                        lower     = low,
                        upper     = upp,
                        control   = optCtrl ) )
    sink()
  }

  #hes <- round(solve(obj$he(opt$par)),10)
  #samps <- mvtnorm::rmvnorm(7, opt$par, hes)

  sdrep <- NULL
  errGrad <- NULL

  if( mode(opt)=="character" )
  {
    rpt <- obj$report()
    rpt$opt$convergence <- 1
  }
  else
  {
    # Retrieve optimized parameters and gradients
    par <- data.frame( par  = names(opt$par),
                       val  = opt$par,
                       grad = as.numeric(obj$gr()) )
  
    # Calculate standard errors via delta method
    sink("sink.txt")
    sdobj <- sdreport( obj )
    sdrpt <- summary( sdobj )
    sink()
    if( mode(sdrpt)!="character" )
    {
      colnames(sdrpt) <- c("val","se")
    
      sdrpt <- as.data.frame(sdrpt) %>%
             mutate( par = rownames(sdrpt),
                     lCI = val - qnorm(.95)*se,
                     uCI = val + qnorm(.95)*se ) %>%
             dplyr::select( par, val, se, lCI, uCI )
    }


    # Build report object
    rpt <- obj$report()
    rpt$opt    <- opt
    rpt$years  <- years
    rpt$gears  <- gears
    rpt$stocks <- stocks
    rpt$par    <- par
    rpt$sdrpt  <- sdrpt
    #rpt$errGrad_st <- errGrad_st

    nObs <- sum(!is.na(data$n_sdtg[1, , , ])) +
            sum(!is.na(data$E_dtg)) + 
            length(data$I_t)
    nPar <- length(opt$par)
    nll  <- opt$objective
    aic  <- 2*nPar + 2*nll + 2*nPar*(nPar+1)/(nObs-nPar-1)
    rpt$aic <- aic


#  library(tmbstan)
#  options(mc.cores = 3)
#  mcinit <- list()
#  for( i in 1:3 )
#    mcinit[[i]] <- rnorm(n=length(obj$par),mean=obj$par,sd=1e-5)
#  
#  fit <- tmbstan( obj = obj,
#                  chains = 3,
#                  iter = 1e3,
#                  init = mcinit )
#  rpt$fit <- fit
  

    if( saveRun )
    {
      plotAll(rpt=rpt,folder=folder)
      save( rpt, file=paste(folder,"/rpt.Rdata",sep="") )
      system( paste("cp ",ctlFile," ",folder,"/estControlFile.txt",sep="") )
    }

  }

  rpt

}