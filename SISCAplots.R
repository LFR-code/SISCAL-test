# Plotting functions for SISCAL TMB model


# plotPriorPosteriorSteep()
# Prior and posterior densities of stock-recruit steepness
plotPriorPosteriorSteep <- function( repList = reports )
{
  if(is.null(repList$posts))
  {
    message("\n\nError: Posteriors missing from report object!\n\n")
    return()
  }

  steepPriorPars <- repList$pars$rSteepBetaPrior

  x <- seq(from = 0.001, to = 1, length.out = 100 )
  priorDens <- dbeta(x, shape1 = steepPriorPars[1], shape2 = steepPriorPars[2] )

  steepPost <- repList$posts$h_ip
  steepPostHist <- hist(  steepPost, 
                          breaks = "FD", 
                          plot = FALSE)

  plot( x = c(0.5,1), y = c(0,max(priorDens,steepPostHist$density)), 
        type = "n", las = 1, xaxs = "i",
        xlab = "Steepness", ylab = "Density")
    hist(steepPost, add = TRUE, probability = TRUE,
          col = "grey60")    
    lines( x = x, y = priorDens, col = "black", lwd = 3, lty = 2)
} # END plotPriorPosteriorSteep()

# plotPriorPosteriorM()
# Prior and posterior densities for sex-specific natural
# mortality in the first year and averaged over time
plotPriorPosteriorM <- function( repList = reports )
{
  if(is.null(repList$posts))
  {
    message("\n\nError: Posteriors missing from report object!\n\n")
    return()
  }

  MpriorPars <- repList$ctlList$hypo$initMprior

  x <- seq( from = log(0.135), 
            to = log(0.155), length.out = 100 )
  logdens <- dnorm( x, mean = log(MpriorPars[1]), sd = MpriorPars[2] )

  postM_ix  <- repList$posts$M_iaxpt[,1,,1,1]
  postM0_ix <- repList$posts$M0_ixp[,,1]

  x <- exp(x)

  # Histogram init M
  histM_male    <- hist(postM_ix[,1], plot = FALSE)
  histM_female  <- hist(postM_ix[,2], plot = FALSE)
  

  # Histogram M0
  histM0_male    <- hist(postM0_ix[,1], plot = FALSE)
  histM0_female  <- hist(postM0_ix[,2], plot = FALSE)  


  par(mfrow = c(2,1), oma = c(3,4,1,1), mar = c(2,1,2,1) )
  plot( x = range(x), y = c(0, max(logdens,histM_male$density,histM_female$density)),
        las = 1, type = "n", xaxs = "i",
        xlab = "Natural mortality (/yr) in 1970",
        ylab = "Density" )
    mtext( side = 3, text = "Initial M", font = 2)

    hist( postM_ix[,1], add = TRUE, freq = FALSE,
          probability = TRUE, col = "steelblue")
    hist( postM_ix[,2], add = TRUE, freq = FALSE,
          probability = TRUE, col = "salmon")
    
    lines(x = x, y = logdens, lwd = 3, col = "black", lty = 2)

  plot( x = range(x), y = c(0, max(logdens,histM0_male$density,histM0_female$density)),
        las = 1, type = "n", xaxs = "i",
        xlab = "Time-averaged natural mortality (/yr) ",
        ylab = "Density" )
    mtext( side = 3, text = "Time-averaged M0", font = 2)
    hist( postM0_ix[,1], add = TRUE, freq = FALSE,
          probability = TRUE, col = "steelblue")
    hist( postM0_ix[,2], add = TRUE, freq = FALSE,
          probability = TRUE, col = "salmon")
    
    lines(x = x, y = logdens, lwd = 3, col = "black", lty = 2)

    mtext( side = 1, outer = TRUE, text = "Natural Mortality (/yr)", line = 1)
    mtext( side = 2, outer = TRUE, text = "Density", line = 3)

} # END plotPriorPosteriorM()


# Plot of MCMC diagnostics (hist of Rhat, ESS)
plotMCdiagnostics <- function( fitObj = fit )
{
  stanfit <- fitObj$stanfit

  sumStanfit <- summary(stanfit)$summary
  Rhat <- sumStanfit[,"Rhat"]
  nEff <- sumStanfit[,"n_eff"]
  cv   <- sumStanfit[,"sd"]/abs(sumStanfit[,"mean"])

  par(mfrow = c(2,2), mar = c(2,2,2,2), oma = c(3,3,1,1) )

  hist( Rhat, xlab = "Rhat", las = 1, main = "" )
  box()
  abline(v = 1.02, lty= 2, col = "red" )
  mtext(side = 1, text = "Rhat", line = 2, font = 2)

  hist( nEff, xlab = "Effective Sample Size", las = 1, main = "")
  box()
  abline( v = 400, lty = 2, col = "red" )
  mtext(side = 1, text = "Effective Sample Size", line = 2, font = 2)

  hist( cv, xlab = "Coefficient of variation",
        las = 1, main = "", breaks = "FD", xlim = c(0,20) )
  mtext(side = 1, text = "CV(theta)", line = 2, font = 2)
  box()

  mtext( side = 3, outer = TRUE, text = "Posterior Chain Diagnostics",
          font = 2, line = -1)

}

# Compare ponded fish in a fleet between two model fits.
comparePondedFish <- function(  fit1 = "fit_parBatWCVI_predBatch4",
                                fit2 = "fit_parBatWCVI_predBatch5",
                                groupFolder = "",
                                fleetIdx = 12,
                                pIdx = 1 )
{
  # First load the first fit
  .loadFit(fit1, groupFolder = groupFolder )
  rpt1 <- reports

  .loadFit(fit2, groupFolder = groupFolder )
  rpt2 <- reports

  reports <- NULL
  gc()

  # Now pull ponded fish
  pondC1_t <- rpt1$repOpt$pondC_pgt[pIdx,fleetIdx,]
  pondC2_t <- rpt2$repOpt$pondC_pgt[pIdx,fleetIdx,]

}

# plotLegalB()
# Time series of legal biomass (above size limit)
plotTotLegalB <- function( repList = reports )
{
  repObj  <- repList$repOpt

  nT <- repObj$nT
  nP <- repObj$nP

  fYear <- repList$fYear

  isPosts <- !is.null(repList$posts)

  SB_pt     <- repObj$SB_pt
  legB_pt   <- repObj$legB_pt  
  totB_pt   <- repObj$B_pt
  yrSizeLim <- repList$data$fYearSizeLim_g[1]

  pRel_axg <- repObj$pRel_axg
  N_axpt   <- repObj$N_axpt
  wt_ax    <- repObj$wt_ax

  if(isPosts)
  {
    # Pull posterior legal, spawning and total B
    posts <- repList$posts

    # SB
    SB_qpt <- apply(X = posts$SB_ipt, MARGIN = 2:3,
                    FUN = quantile, probs = c(0.025,0.5,0.975) )
    SB_pt  <- apply(X = posts$SB_ipt, MARGIN = 2:3,
                    FUN = mean )


    # SB
    legB_qpt <- apply(X = posts$legB_ipt, MARGIN = 2:3,
                    FUN = quantile, probs = c(0.025,0.5,0.975) )
    legB_pt  <- apply(X = posts$legB_ipt, MARGIN = 2:3,
                    FUN = mean )

    totB_pt   <- apply( X = posts$B_ipt, MARGIN = 2:3,
                        FUN = mean )

  }

  


  chkLegB_pt <- legB_pt
  chkTotB_pt <- legB_pt

  for( p in 1:nP )
  {
    for(t in 1:nT)
    {
      if(t > yrSizeLim)
        chkLegB_pt[p,t] <- sum( N_axpt[,,p,t] * (1 - pRel_axg[,,1]) * wt_ax )
      else chkLegB_pt[p,t] <- sum( N_axpt[,,p,t] * wt_ax )
      
      chkTotB_pt[p,t] <- sum( N_axpt[,,p,t] *  wt_ax )
    }
  }

  # browser()

  yrs <- seq(from = fYear, by = 1, length.out = nT)

  par(  mfrow = c(nP,1),
        mar = c(.1,1,.1,1),
        oma = c(3,3,1,1 ) )

  for(p in 1:nP)
  {
    plot(x = range(yrs), y = c(0,max(totB_pt[p,],legB_pt[p,],na.rm = T)),
          axes = FALSE, type = "n" )
      mfg <- par("mfg")
      if(mfg[1] == mfg[3])
        axis(side = 1)
      if(mfg[2] == 1 )
        axis(side = 2, las = 1 )

      grid()
      box()

      if(isPosts)
      {
        spawnPolyCol <- scales::alpha("red", .5)
        legalPolyCol <- scales::alpha("steelblue", .5)
        polygon(  x = c(yrs,rev(yrs)), 
                  y = c(SB_qpt[1,p,1:nT],SB_qpt[3,p,nT:1]),
                  col = spawnPolyCol, border = NA )

        polygon(  x = c(yrs,rev(yrs)), 
                  y = c(legB_qpt[1,p,1:nT],legB_qpt[3,p,nT:1]),
                  col = legalPolyCol, border = NA )

      }

      lines(x = yrs, y = totB_pt[p,1:nT], lwd = 2, col = "black" )
      lines(x = yrs, y = legB_pt[p,1:nT], lwd = 2, col = "steelblue" )
      lines(x = yrs, y = SB_pt[p,1:nT], lwd = 2, col = "red" )
      abline( v = yrs[yrSizeLim+1], lty = 2 )

      if(p == 1)
        legend( x = "topleft",
                legend = c("Total","Legal","Spawning"),
                bty = "n",
                lty = 1, lwd = 2, col = c("black","steelblue","red") )
  }

  mtext(side = 1, outer = TRUE, text = "Year", line = 2)
  mtext(side = 2, outer = TRUE, text = "Biomass (kt)", line = 2)


} # END plotLegalB()

# plotTotMortAtAge()
# Total continuous mortality at age
plotTotMortAtAge <- function( repList = reports,
                              nCols = 2,
                              pIdx = 1,
                              cBrewPal="Paired")
{
  repObj  <- repList$repOpt
  ctsFs   <- approxCtsF( repList = repList )

  F_apgt    <- ctsFs$F_apgt
  nA        <- repObj$nA
  nT        <- repObj$nT
  M_apt     <- repObj$M_apt
  Z_apt     <- repObj$Z_apt
  appZ_apt  <- ctsFs$appZ_apt

  nRows   <- ceiling(nA/nCols)

  fYear   <- repList$fYear
  lYear   <- repList$lYear
  yrs     <- fYear:lYear

  fishingGears <- which(repList$data$fleetType_g > 0)
  gearNames <- names(fishingGears)
  nFishingGears <- length(fishingGears)
  gearCols <- brewer.pal(nFishingGears,cBrewPal)

  if( nFishingGears>12)
    gearCols <- c(gearCols,'grey')


  maxMort <- max(Z_apt[,pIdx,],F_apgt[,pIdx,,], M_apt[,pIdx,])

  par(  mfcol = c(nRows,nCols), 
        oma = c(3,3,3,1),
        mar = c(.1,.5,.1,.5))

    for( a in 1:nA )
    {
      baseMort <- rep(0,nT)
      plot(x = range(yrs), y = c(0,1.05*maxMort),
            type = "n", axes = FALSE )
      mfg <- par("mfg")
      text(x = fYear + 2, y = 0.95*maxMort,
            labels = paste0("Age ",a), font = 2, cex = 1 )
      
      if(mfg[1] == mfg[3])
        axis(side = 1)
      if( mfg[2] == 1)
        axis( side = 2, las = 1 )
      grid()
      box()

      polygon(  x = c(yrs,rev(yrs)),
                y = c(baseMort,M_apt[a,pIdx,nT:1]),
                col = "grey" )
      
      baseMort <- M_apt[a,pIdx,1:nT]
      for( g in 1:nFishingGears )
      {
        gIdx <- fishingGears[g]

        polygon(  x = c(yrs,rev(yrs)),
                  y = c(baseMort,baseMort[nT:1] + F_apgt[a,pIdx,gIdx,nT:1]),
                  col = gearCols[g] )
        baseMort <- baseMort + F_apgt[a,pIdx,gIdx,1:nT]
      }

      # lines(x = yrs, y=  appZ_apt[a,pIdx,1:nT], lwd = 3 )
      
    }

    mtext( side = 2, outer = TRUE, text = "Mortality (/yr)", line = 1.5 )

    par(fig=c(0, 1, 0, 1), oma=c(0, 0, 0, 0), mar=c(0, 0, 0, 0), new=TRUE)
    plot( x = 0, y = 0, type = "n", axes = FALSE, bty = "n")
    legend( x = "top", bty = "n",# xpd = TRUE,
            horiz = TRUE, 
            legend = c("Total","Residual",gearNames),
            col = c("black","grey",gearCols),
            pch = 15, pt.cex = 1.5 )


} # END plotTotMortAtAge

# plotNatMortAtAge()
# Total continuous natural mortality at age (i.e., basal mortality + predation mortality)
plotNatMortAtAge <- function( repList = reports,
                              nCols = 2,
                              pIdx = 1,
                              cBrewPal='Paired')
{
  repObj  <- repList$repOpt
  ctsFs   <- approxCtsF( repList = repList )

  F_apgt    <- ctsFs$F_apgt
  nA        <- repObj$nA
  nT        <- repObj$nT
  M_apt     <- repObj$M_apt
  Z_apt     <- repObj$Z_apt
  appZ_apt  <- ctsFs$appZ_apt

  nRows   <- ceiling(nA/nCols)

  fYear   <- repList$fYear
  lYear   <- repList$lYear
  yrs     <- fYear:lYear

  fishingGears <- repList$data$fleetType_g 
  gearNames <- names(fishingGears)
  predG     <- which(gearNames %in% c("hakeLt50", "hakeGt50", "HS", "SSL", "HB","hbW","hbS"))
  # predG     <- which(gearNames %in% c("hakeLt50", "hakeGt50", "HS", "SSL", "HB","hbW"))


  nPredG <- length(predG)
  gearCols <- brewer.pal(nPredG,cBrewPal)

  if(nPredG>12)
    gearCols <- c(gearCols,'grey')


  predM_apgt <- F_apgt[,,predG,,drop=FALSE]

  F_apt <- apply(predM_apgt, FUN=sum, MAR=c(1,2,4), drop=FALSE)

  aggM_apt <- F_apt[,pIdx,] + M_apt[,pIdx,1:nT]
  maxMort <- max(aggM_apt)



  par(  mfcol = c(nRows,nCols), 
        oma = c(3,3,3,1),
        mar = c(.1,.5,.1,.5))

    for( a in 1:nA )
    {
      baseMort <- rep(0,nT)
      plot(x = range(yrs), y = c(0,1.05*maxMort),
            type = "n", axes = FALSE )
      mfg <- par("mfg")
      text(x = fYear + 2, y = 0.95*maxMort,
            labels = paste0("Age ",a), font = 2, cex = 1 )
      
      if(mfg[1] == mfg[3])
        axis(side = 1)
      if( mfg[2] == 1)
        axis( side = 2, las = 1 )
      grid()
      box()

      polygon(  x = c(yrs,rev(yrs)),
                y = c(baseMort,M_apt[a,pIdx,nT:1]),
                col = "grey" )
      
      baseMort <- M_apt[a,pIdx,1:nT]

      if(nPredG>0)
      {
        for( g in 1:nPredG )
        {
          gIdx <- predG[g]

          polygon(  x = c(yrs,rev(yrs)),
                    y = c(baseMort,baseMort[nT:1] + F_apgt[a,pIdx,gIdx,nT:1]),
                    col = gearCols[g] )
          baseMort <- baseMort + F_apgt[a,pIdx,gIdx,1:nT]

        }      
      }  

      # lines(x = yrs, y=  appZ_apt[a,pIdx,1:nT], lwd = 3 )
      
    }

    mtext( side = 2, outer = TRUE, text = "Natural Mortality (/yr)", line = 1.5 )

    par(fig=c(0, 1, 0, 1), oma=c(0, 0, 0, 0), mar=c(0, 0, 0, 0), new=TRUE)
    plot( x = 0, y = 0, type = "n", axes = FALSE, bty = "n")
    if(nPredG>0)
      legend( x = "top", bty = "n",# xpd = TRUE,
            horiz = TRUE, 
            legend = c("Basal M",gearNames[predG]),
            col = c("grey50",gearCols),
            pch = 15, pt.cex = 1.5 )
    else
      legend( x = "top", bty = "n",# xpd = TRUE,
            horiz = TRUE, 
            legend = c("Basal M"),
            col = c("grey50"),
            pch = 15, pt.cex = 1.5 )


} # END plotTotMortAtAge


# compareAggPosts()
# Reads in two fits and compares their posterior
# time series estimates at the aggregate level
compareAggPosts <- function(  fit1 = "fit_primaryMS",
                              fit2 = "fit_primaryAgg",
                              posts1 = NULL,
                              posts2 = NULL,
                              groupFolder = "",
                              baseDir = NULL )
{
  if(is.null(baseDir))
    baseDir = "Outputs/fits"
  # Load each fit
  if(is.null(posts1))
  {
    .loadFit(fit1, groupFolder= groupFolder, baseDir = baseDir,
              quiet = TRUE)
    if(reports$repOpt$nP > 1)
      posts1 <- makeAggPosts(reports)
    else posts1 <- reports$posts
  }

  if(is.null(posts2))
  {
    .loadFit(fit2, groupFolder= groupFolder, baseDir = baseDir,
              quiet = TRUE)
    if(reports$repOpt$nP > 1)
      posts2 <- makeAggPosts(reports)
    else posts2 <- reports$posts  
  }

  nDraws1 <- dim(posts1$B0_ip)[1]
  nDraws2 <- dim(posts2$B0_ip)[1]

  for( i in 1:max(nDraws1,nDraws2))
  {
    if( i <= nDraws1 )
      posts1$SB_ipt[i,,] <- posts1$SB_ipt[i,,]/posts1$B0_ip[i,]

    if( i <= nDraws2 )
      posts2$SB_ipt[i,,] <- posts2$SB_ipt[i,,]/posts2$B0_ip[i,]
  }


  nT <- dim(posts1$SB_ipt)[3]-1
  fYear <- 1951
  yrs <- seq(from = fYear, by = 1, length.out = nT )



  par( mfrow = c(3,1), oma = c(4,5,1,1), mar = c(.1,.1,.1,.1) )
  # Plot biomass
  maxB <- max(posts1$SB_ipt, posts2$SB_ipt, na.rm = T)
  B01   <- mean(posts1$B0_ip, na.rm = T)
  B02   <- mean(posts2$B0_ip, na.rm = T)
  polyCols <- scales::alpha(c("grey30","salmon"),alpha = .5)
  plot( x = range(yrs), y = c(0,maxB),
        type = "n", axes = FALSE)
    grid()
    axis(side = 2, las = 1 )
    mtext( side = 2, text = expression( paste("Depletion ", B[t]/B[0],"")), line = 3.5 )
    box()
    B1_qt <- apply(X = posts1$SB_ipt[,,1:nT], FUN = quantile, MARGIN = 2,
                    probs = c(0.025, 0.5, 0.975), na.rm = T )
    
    muB1_t <- apply(X = posts1$SB_ipt[,,1:nT], FUN = mean, MARGIN = 2 )
    B2_qt <- apply(X = posts2$SB_ipt[,,1:nT], FUN = quantile, MARGIN = 2,
                    probs = c(0.025, 0.5, 0.975), na.rm = T )
    muB2_t <- apply(X = posts2$SB_ipt[,,1:nT], FUN = mean, MARGIN = 2 )

    polygon(x = c(yrs, rev(yrs)), y = c(B1_qt[1,1:nT],rev(B1_qt[3,1:nT])),
            col = polyCols[1], border = NA )
    lines( x = yrs, y = muB1_t, lwd = 3, lty = 1 )

    polygon(x = c(yrs, rev(yrs)), y = c(B2_qt[1,1:nT],rev(B2_qt[3,1:nT])),
            col = polyCols[2], border = NA )
    lines( x = yrs, y = muB2_t, lwd = 3, lty = 2, col = "red" )
    abline(h = 1, lty = 2, lwd = 0.8 )
    legend( x = "topright",
            legend = c("SISCAH","SCAH"),
            col = c("black","red"),
            pch = 22,
            pt.lwd = 0, cex = 1.5,
            pt.cex = 3, bty = "n",
            pt.bg = polyCols,
            lwd = 3, lty = c(1,2))

  # Now recruitment
  # Do in jittered dot/line plots
  maxR  <- max(posts1$R_ipt, posts2$R_ipt, na.rm = T )
  R01   <- mean(posts1$R0_ip)
  R02   <- mean(posts2$R0_ip)
  logmaxR <- log(maxR, base = 10)
  yAxisLabs <- c(1,10,100,1000,10000)
  yAxisLocs <- log(yAxisLabs,base = 10)
  
  plot( x = range(yrs), y = c(1,logmaxR),
         axes = FALSE, type = "n")
    axis( side = 2, at = yAxisLocs, labels = yAxisLabs, las = 1 )
    grid()
    box()
    mtext( side = 2, text = "Recruitment (1e6)", line = 3.5 )
    R1_ipt <- posts1$R_ipt
    R2_ipt <- posts2$R_ipt

    R1_qt <- apply(X = log(R1_ipt[,,1:nT],base = 10), FUN = quantile, MARGIN = 2,
                    probs = c(0.025, 0.5, 0.975), na.rm = T )
    muR1_t <- apply(X = log(R1_ipt[,,1:nT],base = 10), FUN = mean, MARGIN = 2 )
    R2_qt <- apply(X = log(R2_ipt[,,1:nT],base = 10), FUN = quantile, MARGIN = 2,
                    probs = c(0.025, 0.5, 0.975), na.rm = T )
    muR2_t <- apply(X = log(R2_ipt[,,1:nT],base = 10), FUN = mean, MARGIN = 2 )

    abline(h = log(R01,base = 10), lty = 2, col = "black", lwd = 0.8)
    abline(h = log(R02,base = 10), lty = 2, col = "red", lwd = 0.8)

    segments( x0 = yrs - .3, y0 = R1_qt[1,], y1 = R1_qt[3,],
              col = "grey30", lwd = 1.5)
    points( x = yrs-.3, y = muR1_t, lwd = 3, lty = 1, 
            pch = 16, cex = 1  )

    segments( x0 = yrs + .3, y0 = R2_qt[1,], y1 = R2_qt[3,],
              col = "salmon", lty = 2, lwd = 1.5 )
    points( x = yrs + .3, y = muR2_t, lwd = 3, lty = 1, col = "red",
            pch = 16, cex = 1  )   

    legend( x = "bottomright",
            legend = c("SISCAH","SCAH"),
            col = c("black","red"),
            pch = 16,
            cex = 1.5,
            pt.cex = 1.75, bty = "n",
            lwd = 1.5, lty = c(1,2)) 


  # Now mortality
  M1_it <- posts1$M_iapt[,2,1,]
  M2_it <- posts2$M_iapt[,2,1,]
  maxM <- max(M1_it,M2_it, na.rm = T)
  plot( x = range(yrs), y= c(0,maxM), axes = FALSE, type = "n" )
    axis( side = 2, las = 1)
    axis( side = 1 )
    grid()
    box()
    mtext( side = 2, text = "Natural Mortality (/yr)", line = 3.5)
    M1_qt <- apply(X = M1_it[,1:nT], FUN = quantile, MARGIN = 2,
                    probs = c(0.025, 0.5, 0.975), na.rm = T )
    
    muM1_t <- apply(X = M1_it[,1:nT], FUN = mean, MARGIN = 2 )
    M2_qt <- apply(X = M2_it[,1:nT], FUN = quantile, MARGIN = 2,
                    probs = c(0.025, 0.5, 0.975), na.rm = T )
    muM2_t <- apply(X = M2_it[,1:nT], FUN = mean, MARGIN = 2 )

    polygon(x = c(yrs, rev(yrs)), y = c(M1_qt[1,],rev(M1_qt[3,])),
            col = scales::alpha("grey30",.5), border = NA )
    lines( x = yrs, y = muM1_t, lwd = 3, lty = 1 )

    polygon(x = c(yrs, rev(yrs)), y = c(M2_qt[1,],rev(M2_qt[3,])),
            col = scales::alpha("salmon",.5), border = NA )
    lines( x = yrs, y = muM2_t, lwd = 3, lty = 2, col = "red" )

    legend( x = "bottomright",
            legend = c("SISCAH","SCAH"),
            col = c("black","red"),
            pch = 22,
            pt.lwd = 0, cex = 1.5,
            pt.cex = 3, bty = "n",
            pt.bg = polyCols,
            lwd = 3, lty = c(1,2))


  mtext( side = 1, text = "Year", outer = TRUE, line = 2.5)
} # END compareAggPosts

# plotDDmort()
# Plots density dependent mortality relationship
# from SISCA.
plotDDmort <- function( repObj )
{
  # Pull mortality and biomass
  M_apt <- repObj$M_apt
  B_pt  <- repObj$B_pt

  juveMage <- repObj$juveMage

  # eqbm values
  totB0_p <- repObj$totB0_p
  M0_p    <- repObj$M0_p

  # parameters
  M_p     <- repObj$M_p
  m1_p    <- repObj$m1_p
  nP      <- repObj$nP

  D_pt <- B_pt
  for( p in 1:nP )
    D_pt[p,] <- B_pt[p,]/totB0_p[p]

  Dseq <- seq(from = 0, to = max(D_pt), length.out = 1000 )

  M_pb <- array(0, dim = c(nP,length(Dseq)) )
  for( p in 1:nP )
    M_pb[p,] <- M_p[p] + exp(-m1_p[p] * Dseq )

  par(mfrow = c(nP,1), mar = c(.1,1,.1,1), oma = c(3,4,1,1) )
  for( p in 1:nP )
  {
    plot( x = c(0,max(D_pt)), y = c(0,max(M_pb[p,])), type = "n",
          axes = FALSE, xaxs = "i", yaxs = "i" )
      mfg <- par("mfg")
      if(mfg[1] == mfg[3])
        axis(side = 1)
      axis( side = 2, las = 1)
      grid()
      box()
      lines( x = Dseq, y = M_pb[p,], lwd = 2, col = "salmon" )
      points( x = D_pt[p,], y = M_apt[juveMage+1,p,], pch = 16 )
      segments(x0 = 1, y0 = 0, y1 = M0_p[p], lty = 2 )
      segments(x0 = 0, x1 = 1, y0 = M0_p[p], lty = 2 )
      
      # rmtext( txt = )
  }
  mtext( side = 1, text = "Age 2+ biomass depletion", line = 2, outer = TRUE)
  mtext( side = 2, text = "Age 2+ mortality (/yr)", line = 2, outer = TRUE)

} # END plotDDmort()

# compareAggMLEs()
# Reads in two fits and compares their posterior
# time series estimates at the aggregate level
compareAggMLEs <- function( fit1 = "fit_condMS3_SOG",
                            fit2 = "fit_SOG_RWM",
                            modelNames = NULL,
                            groupFolder = "",
                            baseDir = NULL )
{
  if(is.null(modelNames))
    modelNames <- c(fit1,fit2)

  if(is.null(baseDir))
    baseDir = "Outputs/fits"
  # Load each fit
  .loadFit(fit1, groupFolder= groupFolder, baseDir = baseDir,
              quiet = TRUE)
  fit1 <- reports

  .loadFit(fit2, groupFolder= groupFolder, baseDir = baseDir,
              quiet = TRUE)
  fit2 <- reports

  

    
  nT <- dim(fit1$repOpt$SB_pt)[2]-1
  fYear <- 1951
  yrs <- seq(from = fYear, by = 1, length.out = nT )

  # Aggregate Bt

  SB1_t <- apply(X = fit1$repOpt$SB_pt, FUN = sum, MARGIN = 2)
  SB2_t <- apply(X = fit2$repOpt$SB_pt, FUN = sum, MARGIN = 2)


  par( mfrow = c(3,1), oma = c(4,5,1,1), mar = c(.1,.1,.1,.1) )
  # Plot biomass
  maxB <- max(SB1_t, SB2_t, na.rm = T)
  B01   <- sum(fit1$repOpt$B0_p)
  B02   <- sum(fit2$repOpt$B0_p)
  # polyCols <- scales::alpha(c("grey30","salmon"),alpha = .5)
  plot( x = range(yrs), y = c(0,maxB),
        type = "n", axes = FALSE)
    grid()
    axis(side = 2, las = 1 )
    mtext( side = 2, text = "Spawning Biomass (kt)", line = 3 )
    box()
    lines( x = yrs, y = SB1_t[1:nT], lwd = 3, lty = 1, col = "red" )
    abline( h = B01, col = "red", lwd = 2, lty = 2)
    lines( x = yrs, y = SB2_t[1:nT], lwd = 3, lty = 1, col = "blue" )
    abline( h = B02, col = "blue", lwd = 2, lty = 2)

  
    legend( x = "topright", bty = "n",
            legend = c(modelNames,"Unfished biomass"),
            col = c("red","blue","black"),
            lwd = c(3,3,2), lty = c(1,1,2))

  # Now recruitment
  # Do in jittered dot/line plots
  R1_t  <- apply(X = fit1$repOpt$R_pt, FUN = sum, MARGIN = 2)
  R2_t  <- apply(X = fit2$repOpt$R_pt, FUN = sum, MARGIN = 2)
  maxR  <- max(R1_t, R2_t, na.rm = T )
  R01   <- sum(fit1$repOpt$R0_p)
  R02   <- sum(fit2$repOpt$R0_p)
  
  plot( x = range(yrs), y = c(1,maxR),
         axes = FALSE, type = "n")
    axis( side = 2, las = 1 )
    grid()
    box()
    mtext( side = 2, text = "Recruitment (1e6)", line = 3 )
    
    abline(h = R01, lty = 2, col = "red", lwd = 2)
    abline(h = R02, lty = 2, col = "blue", lwd = 2)

    
    lines( x = yrs, y = R1_t[1:nT], lwd = 3, lty = 1, 
            col = "red" )

    lines( x = yrs, y = R2_t[1:nT], lwd = 3, lty = 1, col = "blue"  )   

    legend( x = "topright", bty = "n",
            legend = c(modelNames,"Unfished Recruitment"),
            col = c("red","blue","black"),
            lwd = c(3,3,2), lty = c(1,1,2))

    # legend( x = "bottomright",
    #         legend = c("SISCAH","SCAH"),
    #         col = c("black","red"),
    #         pch = 16,
    #         cex = 1.5,
    #         pt.cex = 1.75, bty = "n",
    #         lwd = 1.5, lty = c(1,2)) 


  # Now mortality
  M1_t  <- fit1$repOpt$M_apt[2,1,1:nT]
  M2_t  <- fit2$repOpt$M_apt[2,1,1:nT]
  maxM <- max(M1_t,M2_t, na.rm = T)
  plot( x = range(yrs), y= c(0,maxM), axes = FALSE, type = "n" )
    axis( side = 2, las = 1)
    axis( side = 1 )
    grid()
    box()
    mtext( side = 2, text = "Natural Mortality (/yr)", line = 3)
    lines( x = yrs, y = M1_t, lwd = 3, lty = 1, col = "red" )
    lines( x = yrs, y = M2_t, lwd = 3, lty = 1, col = "blue" )

    # legend( x = "bottomright",
    #         legend = c("No Predators","Predators"),
    #         col = c("red","blue"),
    #         pch = 22,
    #         lwd = 3 )


  mtext( side = 1, text = "Year", outer = TRUE, line = 2.5)
} # END compareAggPosts

# plotProcObsErrors()
plotProcObsErrors <- function(repList)
{
  repObj    <- repList$repOpt
  posts     <- repList$posts
  nT        <- repObj$nT
  nP        <- repObj$nP
  fYear     <- repList$fYear
  yrs       <- seq(from = fYear, by = 1, length.out = nT)

  # Pull first/last rec dev years
  tFirstRecDev_p  <- repList$data$firstRecDev_p
  tLastRecDev_p   <- repList$data$lastRecDev_p

  # errors
  omegaM_pt <- apply(X = posts$omegaM_ipt, FUN = mean, MARGIN = c(2,3))
  SRdevs_pt <- apply(X = posts$SRdevs_ipt, FUN = mean, MARGIN = c(2,3))
  zComb_pt  <- apply(X = posts$zComb_ipt, FUN = mean, MARGIN = c(2,3))

  # Remove zeroes
  omegaM_pt[omegaM_pt == 0] <- NA
  SRdevs_pt[SRdevs_pt == 0] <- NA

  # Clear rec devs before/after first and last rec dev
  for( p in 1:nP )
  {
    SRdevs_pt[p,1:tFirstRecDev_p[p]] <- NA
    SRdevs_pt[p,tLastRecDev_p[p]:nT] <- NA
  }
  # zComb_pt[zComb_pt == 0] <- NA

  # Get standard errors
  tauComb_pt  <- repObj$tauComb_pt
  sigmaR      <- repObj$sigmaR
  sigmaM      <- repObj$sigmaM

  stockCols   <- c("darkgreen","salmon","steelblue")
  stockPch    <- 15 + 1:nP
  xJitter     <- c(-.25,0,+.25)
  stockLty    <- c(1,2,4)


  par( mfrow = c(2,1), mar = c(.5,.5,.5,.5), 
        oma = c(5,5,2,2) )
  # Mortality
  plot( x = range(yrs), y = c(-3,3),
        type = "n", axes = FALSE  )
    grid()
    axis(side = 2, las = 1)
    box()
    mtext( side = 2, text = expression(omega[M]), line = 2 )
    abline( h = 0, lty = 2, lwd = 1 )
    for( p in 1:nP )
    {
      lines( x = xJitter[p] + yrs[2:nT], y = omegaM_pt[p,1:(nT-1)],
              col = stockCols[p], lty = stockLty[p], lwd = 2 )
      # points( x = xJitter[p] + yrs[2:nT], y = omegaM_pt[p,1:(nT-1)],
      #         col = stockCols[p], pch = stockPch[p], cex = 2 )
    }
    legend( x = "topright", bty = "n",
            legend = c("C/S","JP/S","Lou"),
            lty = stockLty,
            lwd = 2, cex = 1.5,
            col = stockCols )

  # Recruitment
  plot( x = range(yrs), y = c(-3,3),
        type = "n", axes = FALSE  )
    grid()
    axis(side = 2, las = 1)
    axis( side = 1 )
    box()
    mtext( side = 2, text = expression(omega[R]), line = 2 )
    abline( h = 0, lty = 2, lwd = 1 )
    for( p in 1:nP )
    {
      lines( x = xJitter[p] + yrs, y = SRdevs_pt[p,1:nT],
              col = stockCols[p], lty = stockLty[p], lwd = 2 )
      # points( x = xJitter[p] + yrs, y = SRdevs_pt[p,1:nT],
      #         col = stockCols[p], pch = stockPch[p], cex = 2 )
    }

  # # Obs errors
  # plot( x = range(yrs), y = c(-3,3),
  #       type = "n", axes = FALSE  )
  #   grid()
  #   axis(side = 2, las = 1)
  
  #   box()
  #   mtext( side = 2, text = expression(z[p][t]), line = 2 )
  #   abline( h = 0, lty = 2, lwd = 1 )
  #   for( p in 1:nP )
  #   {
  #     lines( x = xJitter[p] + yrs, y = zComb_pt[p,],
  #             col = stockCols[p], lty = stockLty[p], lwd = 1.5 )
  #     points( x = xJitter[p] + yrs, y = zComb_pt[p,],
  #             col = stockCols[p], pch = stockPch[p], cex = 2 )
  #   }

  mtext( side = 1, text = "Year", line = 3.5 )


}


# rmtext()
# Refactored procedure to plot right hand inner
# margin mtext with the bottom towards the middle
# of the plot
rmtext <- function( line = 1, 
                    txt = "Sample", 
                    font = 1,
                    cex = 1,
                    outer = FALSE,
                    yadj = .5)
{
  corners <- par("usr") #Gets the four corners of plot area (x1, x2, y1, y2)
  if( outer )
    par(xpd = NA) #Draw outside the figure region
  if( !outer )
    par(xpd = TRUE)

  xRange <- corners[2] - corners[1]


  text( x = corners[2]+line*xRange, 
        y = yadj * sum(corners[3:4]), 
        labels = txt, srt = 270,
        font = font, cex = cex )
  par(xpd = FALSE)
} # END rmtext()

# Plot multi-panel time series plots
plotMulti <- function(  ts = c("Nt","Bt","Ft"),
                        repList = reports,
                        labcex = 0.8, heading = "dimLabs" )
{

  nP <- repList$repOpt$nP
  par(mfcol = c(length(ts),nP), oma = c(3,4,3,5), mar = c(0.2,3,0.2,2),
        cex.lab = labcex )

  argList <- list(  repList = repList,  
                    noPar = TRUE, 
                    labcex = labcex )

  for( pIdx in 1:nP)  
  {
    argList$pIdx = pIdx
    for( tsName in ts )
    {
      if( tsName %in% c("SBt","SBtIdx") )
        plotArgs <- c(argList, list(heading = heading) )
      else plotArgs <- argList
      plotCall <- paste("plot",tsName,sep = "")
      do.call(  what = eval(plotCall), 
                args = plotArgs, 
                envir = .GlobalEnv )
    }
  }
}



# Plot biomass
plotSBt <- function(  repList = reports,  
                      noPar = FALSE, 
                      pIdx = 1,
                      labcex = .8,
                      heading = NULL )
{
  report    <- repList$repOpt
  initYear  <- repList$fYear


  # Pull stuff from report
  SBt       <- report$SB_pt[pIdx,]
  Cgt       <- report$totC_pgt[pIdx,,]

  B0        <- signif(report$B0_p[pIdx],3)
  M0_x      <- signif(report$M0_xp[,pIdx],3)

  plotCI <- FALSE

  if(!is.null(repList$posts))
  {
    plotCI <- TRUE
    tInitModel_p <- repList$repOpt$tInitModel_p
    SB_it <- repList$posts$SB_ipt[,pIdx,]

    SB_qt <- apply(X = SB_it, FUN = quantile,
                    MARGIN = c(2), probs = c(0.025, 0.5, 0.975),
                    na.rm = TRUE )

    SBt <- apply(X = SB_it, FUN = mean, MARGIN = c(2))

    B0 <- signif(mean(repList$posts$B0_ip[,pIdx]),3)
    M  <- signif(mean(repList$posts$M_iapt[,1,pIdx,tInitModel_p[pIdx]+ 1]),3)

  } else {
    SB_qt <- 0
  }

  SBt[SBt == 0] <- NA

  # Sum catch
  Ct        <- apply(X = Cgt, FUN = sum, MARGIN = 2)

  stockLabs <- dimnames(report$SB_pt)[[1]][pIdx]

  if( !is.null(heading))
  {
    if( heading == "dimLabs")
      stockLabel <- stockLabs
    else
      stockLabel <- try(heading[stockLabs])

    if( class(stockLabel) == "try-error" ) 
      stockLabel <- stockLabs
  }

  
  # Get number of time steps
  nT    <- report$nT
  years <- seq(from = initYear, length = nT+1, by = 1)
  vertLines <- seq(from = 1950, to = max(years), by = 10)

  depRange  <- range(SBt/B0, na.rm = T)
  maxDep    <- ceiling(depRange[2])

  depAxis   <- round(seq(0, maxDep, length.out = 5),3)

  # Set up plotting area
  if(!noPar)
    par(mfrow = c(1,1), mar = c(.5,3,.5,3), oma = c(3,1,3,1) )

  # Now plot recruitments
  plot( x = range(years), y = c(0,max(SBt,SB_qt, na.rm =T) ),
        type = "n", xlab = "", ylab = "",
        las = 1, axes = FALSE )
    mfg <- par("mfg")
    # Add SSB axis
    axis(side = 2, las =1 )
    # Add depletion axis
    axis( side = 4, las = 1, labels = depAxis, at = B0 * depAxis)

    if( mfg[1] == mfg[3] )
      axis(side = 1 )

    box()
    grid()
    rect( xleft = years[1:nT] - .3, xright = years[1:nT] + .3,
          ybottom = 0, ytop = Ct, border = NA, col = "grey60" )
    # Plot Biomass
    if( plotCI )
    {
      polyCol <- scales::alpha("red", .5)
      polygon(  x = c(years[1:nT],rev(years[1:nT])),
                y = c(SB_qt[1,1:nT], rev(SB_qt[3,1:nT])),
                col = polyCol, border = NA )
    }
    lines( x = years[1:nT], y = SBt[1:nT], lwd = 2, col = "red" )
    # abline( v = vertLines, lwd = .8, lty = 2, col = "grey80")
    abline( h = B0, lty = 3, lwd = .8, col = "red" )
    points(x = years[(nT+1)], y = SBt[(nT+1)], col = "black", bg = "red", pch = 21)
    
    if( mfg[2] == 1 )
      mtext(side = 2, text = "Spawning \n Biomass (kt)", line = 3.5, cex = labcex)
    if( mfg[2] == mfg[4] )
      mtext(side = 4, text = "Depletion", line = 3, cex = labcex)
    if( mfg[1] == 1 & !is.null(heading) )
      mtext(side = 3, text = stockLabel, line = 1, font = 2, cex = labcex )

    panLab( x = 0.7, y = 0.8, txt = paste( "B0 = ",  B0, sep = "") )
    panLab( x = 0.7, y = 0.7, txt = paste( "Mm = ",  M_x[1], sep = "") )
    panLab( x = 0.7, y = 0.6, txt = paste( "Mf = ",  M_x[2], sep = "") )
} # END plotSBt()

# Plot biomass
plotSBtIdx <- function( repList = reports,  
                        noPar = FALSE, 
                        pIdx = 1,
                        labcex = 1,
                        heading = NULL,
                        plotCt = TRUE )
{
  initYear <- repList$fYear
  report   <- repList$repOpt

  # Pull stuff from report
  SBt       <- report$SB_pt[pIdx,]
  Cgt       <- report$totC_pgt[pIdx,,]
  Igt       <- report$I_pgt[pIdx,,]
  qg        <- report$qhat_pg[pIdx,]
  rI_gt     <- repList$data$rI_pgt[pIdx,,]
  vulnBgt   <- report$vulnB_pgt[pIdx,,]
  vulnNgt   <- report$vulnN_pgt[pIdx,,]

  combIdx <- FALSE

  if( any(repList$data$whichCombIdx_g > 0) )
  {
    combIdx <- TRUE
    # Replace q values with qComb_g
    qg <- report$qComb_pg[pIdx,]
  }
  
  vulnBgt[ vulnBgt == 0 ] <- NA
  vulnNgt[ vulnNgt == 0 ] <- NA

  # Get survey type and whether
  # idx is relative or absolute
  surveyType  <- report$survType_g
  indexType   <- report$indexType_g
  calcIndex   <- report$calcIndex_g
  surveyGears <- which(calcIndex == 1)

  # Get number of time steps/gear types
  nT    <- report$nT
  nG    <- report$nG

  B0        <- signif(report$B0_p[pIdx],3)
  M0_x      <- signif(report$M0_xp[,pIdx],3)

  plotCI <- FALSE

  if(!is.null(repList$posts))
  {
    plotCI <- TRUE
    tInitModel_p <- repList$repOpt$tInitModel_p
    SB_it <- repList$posts$SB_ipt[,pIdx,]

    SB_qt <- apply(X = SB_it, FUN = quantile,
                    MARGIN = c(2), probs = c(0.025, 0.5, 0.975),
                    na.rm = TRUE )

    SBt <- apply(X = SB_it, FUN = mean, MARGIN = c(2))

    B0 <- signif(mean(repList$posts$B0_ip[,pIdx]),3)
    M0 <- signif(mean(repList$posts$M0_ip[,pIdx]),3)

    q_ig  <- repList$posts$q_ipg[,pIdx,]
    qg   <- apply( X = q_ig, FUN = mean, MARGIN = 2)

  } else {
    SB_qt <- array(0, dim = c(3,nT+1))
  }



  SBt[ SBt == 0 ] <- NA
  SB_qt[SB_qt == 0] <- NA

  # Sum catch
  Ct        <- apply(X = Cgt, FUN = sum, MARGIN = 2)

  stockLabs <- dimnames(report$SB_pt)[[1]][pIdx]

  if( !is.null(heading))
  {
    if( heading == "dimLabs")
      stockLabel <- stockLabs
    else
      stockLabel <- try(heading[stockLabs])

    if( class(stockLabel) == "try-error" ) 
      stockLabel <- stockLabs
  }


  Igt[Igt<0] <- NA

  scaledIndices <- Igt

  for( g in surveyGears )
  {
    posIdx <- which(Igt[g,1:nT] > 0)

    if( surveyType[g] == 1 )
      vulnTS <- rI_gt[g,1:nT] * vulnNgt[g,1:nT]
    
    if( surveyType[g] == 0 )
      vulnTS <- rI_gt[g,1:nT] * vulnBgt[g,1:nT]

    if( surveyType[g] == 2 )
      vulnTS <- rI_gt[g,1:nT] * SBt[1:nT]
    
    # Compute scaled indices
    scaledIndices[g,posIdx] <- Igt[g,posIdx] / qg[g]  * SBt[posIdx] / vulnTS[posIdx] 

  }

  tInitModel <- report$tInitModel_p[pIdx] + 1

  scaledCombIdx <- repList$data$combI_pt[pIdx,]
  scaledCombIdx[scaledCombIdx < 0 ] <- NA

  if( combIdx )
  {
    
    probPosCombIdx <- report$probPosCombIdx_pt[pIdx,]
    qComb_t <- report$qComb_pt[pIdx,]
    scaledCombIdx <- scaledCombIdx / qComb_t / probPosCombIdx
  }


  # Create x axis label vector and vertical lines
  # for easy multipanel plotting
  years <- seq(from = initYear, length = nT+1, by = 1)
  vertLines <- seq(from = 1950, to = max(years), by = 10)

  cols <- brewer.pal(nG,"Paired")

  depRange  <- range(SBt/B0, na.rm = T)
  maxDep    <- ceiling(depRange[2])

  depAxis   <- round(seq(0,maxDep, length.out = 5),2)

  # Set up plotting area
  if(!noPar)
    par(mfrow = c(1,1), mar = c(.5,3,.5,3), oma = c(3,1,3,1) )

  yMax <- max(SBt,scaledIndices, na.rm =T)

  # Now plot recruitments
  plot( x = range(years), y = c(0,yMax ),
        type = "n", xlab = "", ylab = "",
        las = 1, axes = FALSE )
    mfg <- par("mfg")
    # Add SSB axis
    axis(side = 2, las =1 )
    # Add depletion axis
    axis( side = 4, las = 1, labels = depAxis, at = B0 * depAxis)

    if( mfg[1] == mfg[3] )
      axis(side = 1 )

    box()
    grid()
    # Plot data
    if( combIdx )
    {
      posIdx  <- which(scaledCombIdx > 0)
      zeroIdx <- which(scaledCombIdx == 0)
      zeroIdx <- zeroIdx[zeroIdx >= tInitModel ]
      points( x = years[posIdx], y = scaledCombIdx[posIdx],
              col = "grey50", pch = 16)

      # points( x = years[zeroIdx], y = scaledCombIdx[zeroIdx],
      #         col = "grey40", pch = 1)
      axis( side = 3, at = years[zeroIdx], labels = FALSE,
            col.ticks = "grey40", tck = .02, lwd.ticks = 3 )
    }

    # Plot Biomass
    if( plotCI )
    {
      polyCol <- scales::alpha("red", .5)
      polygon(  x = c(years[1:nT],rev(years[1:nT])),
                y = c(SB_qt[1,1:nT], rev(SB_qt[3,1:nT])),
                col = polyCol, border = NA )
    }
    abline( h = B0, lty = 3, lwd = .8, col = "red" )
    lines( x = years[1:nT], y = SBt[1:nT], lwd = 2, col = "red" )
    if(plotCt)
      rect( xleft = years[1:nT] - .3, xright = years[1:nT] + .3,
            ybottom = 0, ytop = Ct, border = NA, col = "grey60" )
    # points(x = years[(nT+1)], y = SBt[(nT+1)], col = "black", bg = "red", pch = 21)
    
    if( mfg[2] == 1 )
      mtext(side = 2, text = "Spawning\nBiomass (kt)", line = 3.5, cex = labcex)
    if( mfg[2] == mfg[4] )
      rmtext(txt = "Depletion", line = .1, cex = 1.5, outer = TRUE)
    if( mfg[1] == 1 & !is.null(heading) )
      mtext(side = 3, text = stockLabel, line = 1, font = 2, cex = labcex )

    for( g in surveyGears )
    {
      posIdx  <- which(scaledIndices[g,] > 0)
      zeroIdx <- which(scaledIndices[g,] == 0)
      points( x = years[posIdx], y = scaledIndices[g,posIdx],
            col = alpha(cols[g],.5), pch = 16 )
      points( x = years[zeroIdx], y = scaledIndices[g,zeroIdx],
            col = alpha(cols[g],.5), pch = 1 )
    }

    
    text( x = years[nT - 10], y = c(0.9,.8,.7,.6,.5) * yMax,
          labels = c( paste( "B0 = ",  B0, sep = ""),
                      paste( "M0 = ", paste(M0_x,collapse = ","), sep = ""),
                      paste( "qRV = ",  signif(qg[5],2), sep = ""),
                      paste( "qHSf = ",  signif(qg[6],2), sep = ""),
                      paste( "qHSr = ",  signif(qg[7],2), sep = "")),

          cex = .9 )
} # END plotSBtIdx



# Plot fishing mortality
plotNt <- function( report = repFE, 
                    initYear = fYear, 
                    noPar = FALSE, pIdx = 1 )
{
  # Pull stuff from report
  Nat       <- report$N_atp[,,pIdx]
  Nt        <- colSums(Nat, na.rm = T)
  
  # Get number of time steps
  nT    <- report$nT
  years <- seq(from = initYear, length = nT+1, by = 1)
  vertLines <- seq(from = 1950, to = max(years), by = 10)

  # Set up plotting area
  if(!noPar)
    par(mfrow = c(1,1), mar = c(.5,3,.5,3), oma = c(3,1,3,1) )

  # Now plot recruitments
  plot( x = range(years), y = c(0,max(Nt, na.rm =T) ),
        type = "n", xlab = "", ylab = "",
        las = 1, axes = FALSE )
    mfg <- par("mfg")
    axis(side = 2, las =1 )
    if( mfg[1] == mfg[3] )
      axis(side = 1 )
    box()
    # Plot recruitment
    abline( v = vertLines, lwd = .8, lty = 2, col = "grey80")
    lines( x = years, y = Nt, lwd = 2, col = "grey40" )
    if(mfg[2] == 1)
      mtext(side = 2, text = "Numbers (1e6)", line = 3)
}


# Plot fishing mortality
plotFtg <- function(  repList = reports,  
                      noPar = FALSE, 
                      pIdx = 1,
                      labcex = 1,
                      cBrewPal='Paired' )
{
  report    <- repList$repOpt
  initYear  <- repList$fYear

  # Pull stuff from report
  U_gt     <- report$U_pgt[pIdx,,]
  gearIDs <- dimnames(report$F_pgt)[[2]]
  nG      <- report$nG
  nT      <- report$nT
  lUmsy_p <- report$refPts$FmsyRefPts$lUmsy_p

  legB_t  <- report$legB_pt[pIdx,]

  # Pull catch
  C_gt   <- report$C_pgt[pIdx,,]
  C_t    <- apply( X = C_gt, FUN = sum, MARGIN = 2)

  if(!is.null(repList$posts))
  {
    tInitModel_p <- repList$repOpt$tInitModel_p
    U_igt <- repList$posts$U_ipgt[,pIdx,,]
    U_gt  <- apply( X = U_igt, FUN = mean,
                    MARGIN = c(2,3), na.rm = T)
    lUmsy_p <- apply(X = repList$posts$lUmsy_ip, FUN = mean, MARGIN = 2)

    legB_t <- apply(  X = repList$posts$legB_ipt[,pIdx,],
                      FUN = mean, MARGIN = 2 )
  } 

  U_t    <- C_t/legB_t
  commGears <- which(repList$ctlList$data$fleetType[gearIDs] != 0 )

  minTime_g <- rep(0,nG)
  maxTime_g <- rep(nT,nG)
  for( g in commGears )
  {
    minTime_g[g] <- min(which(U_gt[g,] > 0),na.rm = T)
    maxTime_g[g] <- max(which(U_gt[g,] > 0),na.rm = T)
  }
  minTime_g[!is.finite(minTime_g)] <- 1
  maxTime_g[!is.finite(maxTime_g)] <- nT


  cols    <- brewer.pal( length(commGears), cBrewPal )

  years <- seq(from = initYear, length.out = nT+1, by = 1)
  vertLines <- seq(from = 1950, to = max(years), by = 10)

  # Set up plotting area
  if(!noPar)
    par(mfrow = c(1,1), mar = c(.5,3,.5,3), oma = c(3,1,3,1) )

  # Now plot F series
  plot( x = range(years), y = c(0, min(1,2*max(U_gt,U_t,na.rm = T)) ),
        type = "n", xlab = "", ylab = "",
        las = 1, axes = FALSE, yaxs = "i" )
    mfg <- par("mfg")
    axis(side = 2, las =1 )
    if( mfg[1] == mfg[3] )
      axis(side = 1 )
    box()
    grid()

    # Plot harvest rates by gear
    for(gIdx in 1:length(commGears))
    {
      g <- commGears[gIdx]
      gYrs <-  max(minTime_g[g]-1,1):min(maxTime_g[g]+1,nT)
      lines( x = years[gYrs], y = U_gt[g,gYrs], lwd = 2, col = cols[gIdx] )
    }
    # Now add total legal HR
    lines(  x = years[1:nT],y = U_t, col = "grey30",
            lwd = 3, lty = 2 )

    abline( h = lUmsy_p[pIdx],
            lty = 3, lwd = 2, col = "red")

    
    if(mfg[2] == 1)
    {
      legend( x = "topright", col = cols,
              legend = gearIDs[commGears], lwd = 2, bty = "n")
      mtext(side = 2, text = "Legal\nHarvest Rate", line = 3.5, cex = labcex )
    }

}


# Plot natural mortality
plotMt <- function( repList = reports, 
                    noPar = FALSE, 
                    pIdx = 1,
                    labcex = .8)
{
  report    <- repList$repOpt 
  initYear  <- repList$fYear

  nX <- report$nX
  nA <- report$nA

  # Pull stuff from report
  M_xt      <- report$M_axpt[nA,,pIdx,]
  Mbar      <- report$M
  M_x       <- report$M_xp[pIdx]
  M0_x      <- report$M0_xp[,pIdx]
  nT        <- report$nT
  Mjuve     <- report$Mjuve_p[pIdx]
  juveMage  <- report$juveMage

  # Pull info on predation mortality
  Fgt     <- report$F_pgt[pIdx,,,drop = TRUE]
  gearIDs <- dimnames(report$F_pgt)[[2]]
  nG      <- report$nG
  nT      <- report$nT

  predG   <- which(gearIDs %in% c("hakeLt50", "hakeGt50", "HS", "SSL", "HB"))

  M_xt[M_xt == 0] <- NA

  plotCI <- FALSE

  
  if(!is.null(repList$posts))
  {
    plotCI <- TRUE
    tInitModel_p  <- repList$repOpt$tInitModel_p
    M_ixt         <- repList$posts$M_iaxpt[,nA,,pIdx,]
    
    Mjuve <- NA
    if(juveMage > 0)
    {
      Mjuve_i       <- repList$posts$M_iapt[,juveMage,pIdx,tInitModel_p[pIdx]+1]
      Mjuve         <- mean(Mjuve_i)
    }

    M_qxt <- apply(  X = M_ixt, FUN = quantile,
                      MARGIN = c(2:3), probs = c(0.025, 0.5, 0.975),
                      na.rm = TRUE )

    M_xt  <- apply(X = M_ixt, FUN = mean, MARGIN = c(2:3))

    M0_x <- apply(X = repList$posts$M0_ixp[,,pIdx], FUN = mean, MARGIN = 2)

  } else {
    M_qxt <- NA
  }

  # if(addPredM & length(predG)>0 )
  # {

  #   predM   <- apply(Fgt[predG,], FUN=sum, MAR=c(2))
  #   basalM  <- Mat[juveMage+1,1,1:nT]
  #   Mt    <- predM + basalM

  #   YLIM <- c(0,max(Mt, na.rm =T) )      
  # }  else
  # YLIM <- c(0,max(Mat,M_qapt, na.rm =T) )
 

  years     <- seq(from = initYear, length.out = nT+1, by = 1)
  vertLines <- seq(from = 1950, to = max(years), by = 10)



  # Set up plotting area
  if(!noPar)
    par(mfrow = c(1,1), mar = c(.5,3,.5,3), oma = c(3,1,3,1) )

  # Now plot F series
  plot( x = range(years), y = c(0,max(M_xt,M_qxt, na.rm =T) ),
        type = "n", xlab = "", ylab = "", las = 1, axes = FALSE )
    mfg <- par("mfg")
    axis(side = 2, las =1 )
    if( mfg[1] == mfg[3] )
      axis(side = 1 )
    box()
    grid()
    # Plot credibility intervals
    polyCol <- scales::alpha("salmon", .5)
    if( plotCI )
    {
      polygon(  x = c(years[1:nT],rev(years[1:nT])),
                y = c(M_qxt[1,nX,1:nT], rev(M_qxt[3,nX,1:nT])),
                col = polyCol, border = NA )
    }
    
    # abline( h = M, lty = 2, lwd = 2, col = "salmon")
    lines(  x = years[1:nT], y = M_xt[nX,1:nT], col = "salmon",
              lwd = 3  )
    abline( h = M0_x[nX], lty = 2, lwd = 3, col = "grey50")


    if(mfg[2] == 1)
    {
      mtext(side = 2, text = "Natural\nMortality (/yr)", line = 3.5, cex = labcex)
      legend( x = "topleft", 
              lwd = c(2,2,2,NA),
              lty = c(1,2,2,NA),
              pch = c(22,NA,NA),
              pt.bg = c(polyCol,NA,NA), 
              pt.lwd = c(0,NA,NA),
              pt.cex = c(1.5,NA,NA),
              cex = .8,
              col = c("salmon","grey50"),
              legend = c("Mt","M0"), bty = "n")
    }
}

# plotSlg()
# Selectivity-at-length for each gear
plotSlg <- function(repList = reports )
{
  report <- repList$repOpt
  # pull selectivity from report
  Saxg       <- report$sel_axg
  Saxpgt     <- report$sel_axpgt
  
  # Pull length-at-age from reports
  len_ax     <- report$len_ax

  gearLabs <- dimnames(Saxg)[[3]]
  stockLabels <- dimnames(Saxpgt)[[3]]

  # Update gear cols so they match the commGears
  commGears <- 1:length(gearLabs)

  # model dimensions
  nA      <- report$nA
  nG      <- report$nG
  nT      <- report$nT
  nP      <- report$nP
  nX      <- report$nX


  gearCols <- RColorBrewer::brewer.pal(length(commGears), "Paired")

  maxL <- max(len_ax)

  par(mfrow = c(nG,nP), mar = c(.1,1,.1,1), oma = c(3,3,1,1))
  for(g in commGears)
  {
    plot(x = c(0,maxL), y = c(0,1), type = "n", axes = FALSE )
      mfg <- par("mfg")
      if(mfg[1] == mfg[3])
        axis(side = 1 )
      if(mfg[2] == 1)
        axis(side = 2, las = 1)
      box()
      lines(  x = len_ax[,nX], y = Saxg[,nX,g],
              col = gearCols[g], lwd = 3 )
      if(mfg[2] == mfg[4])
        rmtext( txt = gearLabs[g], outer = TRUE, 
                line = 2 )

  }
  mtext(side = 2, outer = TRUE, text = "Selectivity-at-length", line = 2 )
  mtext(side = 1, outer = TRUE, text = "Length (cm)", line = 2 )

 
} # END plot Slg


# Selectivity at age for each gear
plotSagComp <- function(  repList = list(repFE),
                          gIdx  = 1:3,
                          noPar = TRUE,
                          legend = FALSE,
                          gearLabels = NULL )
{
  SagList <- vector(mode = "list", length = length(repList) )

  nReps <- length(repList)

  # model dimensions
  nA      <- repList[[1]]$repOpt$nA
  nG      <- min(repList$repOpt$nG,length(gIdx))

  legNames <- character(length = 0)

  for( rIdx in 1:nReps )
  {
    SagList[[rIdx]] <- repList[[rIdx]]$repOpt$sel_ag[,gIdx, drop = FALSE]
    legNames <- c(legNames,repList[[rIdx]]$ctlList$ctrl$dataScenarioName)
  }

  gearLabs <- dimnames(SagList[[1]])[[2]]

  if(is.null(gearLabels))
    names(gearLabs) <- gearLabs
  else names(gearLabs) <- gearLabels

  cols <- alpha(brewer.pal( n = nReps, "Set1" ), .8)


  # Make plotting window
  if(!noPar)
    par(  mfrow = c(nG,1), 
          oma = c(3,4,.1,.1), 
          mar = c(0,1,0,1) )

  # loop over gears
  for( g in 1:nG )
  {
    plot( x = c(1,nA), y = c(0,1), type = "n",
          axes = FALSE, xlab = "", ylab = "" )
      mfg <- par("mfg")
      if(mfg[1] == mfg[3])
        axis(side = 1)
      if( mfg[2] == 1)
        axis( side = 2, las =1 )
      box()
      for( rIdx in 1:nReps)
      {
        lines( x = 1:nA, y = SagList[[rIdx]][,g], col = cols[rIdx],
                lwd = 3 )
        points( x = 1:nA, y = SagList[[rIdx]][,g], col = cols[rIdx],
                pch = 16 )
      }
      if( mfg[2] == mfg[4])
        mtext( side = 4, text = names(gearLabs)[g], line = 2 )
  }
  if( legend )
    legend( x = "bottomright", bty = "n",
            legend = legNames, col = cols, lwd = 2 )

  if( !noPar )
  {
    mtext( side = 1, text = "Age", outer = T, line = 1 )
    mtext( side = 2, text = "Selectivity", outer = T, line = 2 )
  }

}


# plotC_pgt()
# Breaks catch out into stacked bars
# by gear type. Shows SOK catch as total ponded
# fish, with a red border for dead ponded
# fish.
plotC_pgt <- function(  repList = reports,
                        yrRange = 1951:2019,
                        cBrewPal='Paired' )
{
  repObj    <- repList$repOpt 
  fYear     <- repList$fYear
  # Pull stuff from report
  C_pgt     <- repObj$totC_pgt
  gearLabs  <- dimnames(repObj$F_pgt)[[2]]
  nG        <- repObj$nG
  nT        <- repObj$nT
  nP        <- repObj$nP

  stockLabs <- dimnames(repObj$F_pgt)[[1]]

  gearCols  <- brewer.pal( nG, 'Paired')

  if(nG>12)
    gearCols <- c(gearCols,'grey')

  fleetType   <- repObj$fleetType_g
  postPondM_g <- repObj$postPondM_g

  gears <- which(fleetType > 0)
  nComm <- length(gears)

  years <- seq(from = fYear, length = nT, by = 1)
  vertLines <- seq(from = 1950, to = max(years), by = 10)

  plotdx <- which(years %in% yrRange)

  Csum_pt <- apply( X = C_pgt, FUN = sum, MARGIN = c(1,3), na.rm = T )

  # Now loop over p
  par( mfrow = c(nP,1), mar = c(0.5,.5,0,0), oma = c(3,3,3,3))

  for( p in 1:nP )
  {
    plot( x = range(years[plotdx]), y = c(0,max(Csum_pt[p,plotdx])),
          type = "n", axes = FALSE )
      mfg <- par("mfg")
      if( mfg[1] == mfg[3] )
        axis( side = 1 )
      axis( side = 2, las = 1)
      box()
      grid()
      mtext( side = 4, text = stockLabs[p], font = 2 )

      for( g in gears )
      {
        if( g == 1 )
          baseC <- rep(0,nT)
        if( g > 1 )
          baseC <- apply( X = C_pgt[p,1:(g-1),,drop = FALSE], FUN = sum,
                          MARGIN = 3 )

        # if SOK, add shading
        if( fleetType[g] == 2 )
        {
          catVec <- (1 - exp(-postPondM_g[g]))*C_pgt[p,g,]
        } else catVec <- C_pgt[p,g,]

        rect( xleft = years - .3, xright = years + .3,
              ybottom = baseC, ytop = baseC + catVec,
              col = gearCols[g], border = NA )

        
        

        if( p == 1 )
          legend( x = "topleft", bty = "n",
                  pt.bg = c(gearCols[gears]),
                  pch = 22,
                  col = c(rep(NA,nComm)),
                  pt.lwd = c(rep(0,nComm)),
                  legend = c(gearLabs[gears]))

      }
  }

  mtext( side = 1, outer = TRUE, text = "Year", line = 2)
  mtext( side = 2, outer = TRUE, text = "Catch (kt)", line = 2 )

} # END plotC_pgt()

# Selectivity at age for each gear
plotSag <- function(  repList = reports, fleets = NULL )
{
  report <- repList$repOpt
  # pull selectivity from report
  Saxg       <- report$sel_axg
  Saxpgt     <- report$sel_axpgt

  gearLabs <- dimnames(Saxg)[[3]]

  stockLabels <- dimnames(Saxpgt)[[3]]

  # Update gear cols so they match the commGears
  commGears <- 1:length(gearLabs)

  # model dimensions
  nA      <- report$nA
  nG      <- report$nG
  nT      <- report$nT
  nP      <- report$nP
  nX      <- report$nX

  if(is.null(fleets))
    fleets <- commGears

  gearCols <- RColorBrewer::brewer.pal(length(commGears), "Paired")

  # Make plotting window
  par(  mfcol = c(length(fleets),nX), 
        oma = c(3,4,.1,.1), 
        mar = c(0.1,2,0.1,2) )
  # loop over gears
  for( x in 1:nX )
    for( gIdx in 1:length(fleets) )
    {
      g <- fleets[gIdx]
      plot( x = c(1,nA), y = c(0,1), type = "n",
            axes = FALSE, xlab = "", ylab = "" )
        mfg <- par("mfg")
        if(mfg[1] == mfg[3])
          axis(side = 1)
        axis( side = 2, las =1 )
        box()

        
        panLab( x = 0.1, y = 0.9, txt = gearLabs[g], font = 2 )
        
        for( p in 1:nP )
        {
          for( t in 1:nT)
            lines( x = 1:nA, y = Saxpgt[,x,p,g,t],
                    col = "grey80", lwd = .8 )

          lines(  x = 1:nA, y = Saxpgt[,x,p,g,1], 
                  col = "grey30", lwd = 1,
                  lty = p + 1 )
        }
        
        lines( x = 1:nA, y = Saxg[,x,g], col = gearCols[gIdx],
                lwd = 3 )

        if( g == 1 )
        {
          legend( x = "bottomright", bty = "n",
                  legend = stockLabels, col = "black", lwd = 2,
                  lty = 1:nP + 1 )
        }

    }
  
  mtext( side = 1, text = "Age", outer = T, line = 2 )
  mtext( side = 2, text = "Selectivity-at-age", outer = T, line = 2 )

} # END plotSag()

# plotIdxFits()
plotIdxFits <- function(  repList = reports,
                          labcex = .8 )
{
  report    <- repList$repOpt 
  initYear  <- repList$fYear

  # Pull vulnerable biomass and indices
  vulnBpgt      <- report$vulnB_pgt
  vulnNpgt      <- report$vulnN_pgt
  SBpt          <- report$SB_pt
  Ipgt          <- report$I_pgt
  qpg           <- report$qhat_pg
  taupg         <- report$tauObs_pg

  # Mixed indices
  mIgt          <- report$mI_gt
  qg            <- report$qhat_g
  taug          <- report$tauObs_g


  # Get survey type and whether
  # idx is relative or absolute
  surveyType  <- report$survType_g
  indexType   <- report$indexType_g
  calcIndex   <- report$calcIndex_g

  # Model dimensions
  nP          <- report$nP
  nG          <- report$nG
  nT          <- report$nT
  surveyGears <- which(calcIndex == 1)

  gearLabels  <- dimnames(vulnBpgt)[[2]]
  stockLabels <- dimnames(vulnBpgt)[[1]]

  # replace missing indices with NAs
  Ipgt[ Ipgt < 0 ] <- NA
  mIgt[ mIgt < 0 ] <- NA

  # Can we combine Ipgt and mItg into one, with 
  # an extra slice for the mixed stock if nP > 1,
  # or just overwrite the single stock gear
  # entries if nP == 1

  if( nP > 1 )
  {
    # Combine the indices
    combinedI_pgt <- array( NA, dim = c(nP+1,nG,nT))
    combinedI_pgt[1:nP,,] <- Ipgt
    combinedI_pgt[nP+1,,] <- mIgt

    Ipgt <- combinedI_pgt

    # now combine the qs
    comb_qpg <- array(NA, dim = c(nP+1,nG))
    comb_qpg[1:nP,] <- qpg
    comb_qpg[nP+1,] <- qg

    qpg <- comb_qpg

    # And the taus
    comb_taupg <- array(NA, dim = c(nP+1,nG))
    comb_taupg[1:nP,] <- taupg
    comb_taupg[nP+1,] <- taug

    taupg <- comb_taupg

    stockLabels <- c( stockLabels, "Mixed stock" )
  }

  if(nP == 1)
  {
    # Loop over gear types, merge the two together
    for( g in surveyGears )
    {
      if( all( is.na(Ipgt[,g,]) ) & any(!is.na(mIgt[g,])) )
      {
        Ipgt[,g,] <- mIgt[g,]
        qpg[,g]   <- qg[g]
        taupg[,g] <- taug[g]
      }
    }    

  }

  cols <- brewer.pal(nG,"Paired")

  if( nP > 1 )
    nSeries <- nP + length(surveyGears) - 1
  else nSeries <- 1

  if(length(surveyGears) > 1)
    par(  mfcol = c(length(surveyGears)+1,nSeries), 
          oma = c(3,4,1,1),
          mar = c(0,1.5,0,1.5) )
  if( length(surveyGears) == 1 )
    par(  mfrow = c(nSeries,1), 
          oma = c(3,4,1,1),
          mar = c(0,1.5,0,1.5) )

  years <- seq(initYear, by = 1, length = nT + 1 )
  vertLines <- seq(from = 1950, to = max(years), by = 10)

  for( sIdx in 1:nSeries)
  {
    if( sIdx == nP + 1)
      pIdx <- 1:nP
    else pIdx <- sIdx
    for( g in surveyGears )
    {
      if( surveyType[g] == 1 )
      {
        vulnTS <- vulnNpgt[pIdx,g,,drop = FALSE]
        vulnTS <- apply( X = vulnTS, FUN = sum, MARGIN = 3)
        yLabel  <- "Vulnerable Numbers (1e6)"
      }
      if( surveyType[g] == 0 )
      {
        vulnTS <- vulnBpgt[pIdx,g,1:nT,drop = FALSE]
        vulnTS <- apply( X = vulnTS, FUN = sum, MARGIN = 3)
        yLabel  <- "Vulnerable Biomass (kt)"
      }
      if( surveyType[g] == 2 )
      {
        vulnTS <- SBpt[pIdx,1:nT,drop = FALSE]
        vulnTS <- apply( X = vulnTS, FUN = sum, MARGIN = 2)
        yLabel  <- "Spawning Biomass (kt)"
      }



      yRange <- range(vulnTS,Ipgt[sIdx,g,]/qpg[sIdx,g], na.rm = T )

      plot( x = range(years), y = yRange,
            xlab = "", ylab = "", type = "n",
            axes = F )
        mfg <- par("mfg")
        if(mfg[1] == mfg[3])
          axis(side = 1)
        if(mfg[1] == 1)
          mtext( side = 3, text = stockLabels[sIdx], font = 2,
                  line = 1, cex = labcex )
        if(mfg[2] == 1)
        {
          axis(side = 2, las = 1)
          mtext(  side = 2, text = yLabel, outer = F, 
                  line = 3, cex = labcex)
        }
        grid()
        box()
        if( mfg[2] == mfg[4] )
          mtext( side = 4, text = gearLabels[g], font = 2, cex = 1,
                  line = 1 )
        # Only plot lines/points if there is data
        if( all(is.na(Ipgt[sIdx,g,]) ) )
          next
        points( x = years[1:nT], y = Ipgt[sIdx,g,]/qpg[sIdx,g], 
                col = cols[g], pch = 16 )
        lines( x = years[1:nT], y = vulnTS, lwd = 2, col = "grey40" )
        panLab( x = 0.8, y = 0.8, 
                txt = paste( "q = ", signif(qpg[sIdx,g],2), sep = "") )
        panLab( x = 0.8, y = 0.7, 
                txt = paste( "tau = ", signif(taupg[sIdx,g],2), sep = "") )
        panLab( x = 0.1, y = 0.9, txt = stockLabels[sIdx] )
        
    }

    # Now plot residuals
    plotItResids(noPar = TRUE,pIdx = pIdx, repList = repList ) 
  }


  mtext(  side = 1, text = "Year", outer = T, line = 2 )
}

# Vulnerable biomass
plotMixedIdxFits <- function( repList = reports )
{
  report    <- repList$repOpt 
  initYear  <- repList$fYear

  # Pull vulnerable biomass and indices
  vulnBtgp      <- report$vulnB_tgp
  vulnNtgp      <- report$vulnN_tgp
  SBpt          <- report$SB_pt
  Itgp          <- report$I_tgp
  qpg           <- report$qhat_gp
  taugp         <- report$tauObs_gp

  # Mixed indices
  mItg          <- report$mI_tg
  qg            <- report$qhat_g
  taug          <- report$tauObs_g

  # Get survey type and whether
  # idx is relative or absolute
  surveyType  <- report$survType_g
  indexType   <- report$indexType_g
  calcIndex   <- report$calcIndex_g

  # Model dimensions
  nP          <- report$nP
  nG          <- report$nG
  nT          <- report$nT
  surveyGears <- which(calcIndex == 1)

  gearLabels  <- dimnames(vulnBtgp)[[2]]
  stockLabels <- dimnames(vulnBtgp)[[3]]

  # replace missing indices with NAs
  mItg[ mItg < 0 ] <- NA
  cols <- brewer.pal(nG,"Paired")

  if(length(surveyGears) > 1)
    par(  mfcol = c(length(surveyGears),1), 
          oma = c(3,4,1,1),
          mar = c(0,1.5,0,1.5) )
  if( length(surveyGears) == 1 )
    par(  mfrow = c(1,1), 
          oma = c(3,4,1,1),
          mar = c(0,1.5,0,1.5) )

  years <- seq(initYear, by = 1, length = nT + 1 )
  vertLines <- seq(from = 1950, to = max(years), by = 10)

  for( g in surveyGears )
  {
    if( surveyType[g] == 1 )
    {
      vulnTS  <- vulnNtgp[,g,1:nP,drop = FALSE]
      vulnTS  <- apply( X = vulnTS, FUN = sum, MARGIN = 1 )
      yLabel  <- "Vulnerable Numbers (1e6)"
    }
    if( surveyType[g] == 0 )
    {
      vulnTS  <- vulnBtgp[,g,1:nP,drop = FALSE]
      vulnTS  <- apply( X = vulnTS, FUN = sum, MARGIN = 1 )
      yLabel  <- "Vulnerable Biomass (kt)"
    }
    if( surveyType[g] == 2 )
    {
      vulnTS  <- SBpt[1:nP,1:nT,drop = FALSE]
      vulnTS  <- apply( X = vulnTS, FUN = sum, MARGIN = 2 )
      yLabel  <- "Spawning Biomass (kt)"
    }

    yRange <- range(vulnTS,mItg[,g]/qg[g], na.rm = T )

    plot( x = range(years), y = yRange,
          xlab = "", ylab = "", type = "n",
          axes = F )
      mfg <- par("mfg")
      if(mfg[1] == mfg[3])
        axis(side = 1)
      if(mfg[1] == 1)
        mtext( side = 3, text = "Mixed", font = 2, cex = 1,
                line = 1 )
      axis(side = 2, las = 1)
      box()
      if( mfg[2] == mfg[4] )
        mtext( side = 4, text = gearLabels[g], font = 2, cex = 1,
                line = 1 )
      # Only plot lines/points if there is data
      if( all(is.na(mItg[,g]) ) )
        next
      points( x = years[1:nT], y = mItg[,g]/qg[g], 
              col = cols[g], pch = 16 )
      lines( x = years[1:nT], y = vulnTS, lwd = 2, col = "grey40" )
      panLab( x = 0.8, y = 0.8, 
              txt = paste( "q = ", signif(qg[g],2), sep = "") )
      panLab( x = 0.8, y = 0.7, 
              txt = paste( "tau = ", signif(taug[g],2), sep = "") )
      panLab( x = 0.1, y = 0.9, txt = "Mixed stock" )
      
  }
  mtext(  side = 1, text = "Year", outer = T, line = 2 )
  mtext(  side = 2, text = yLabel, outer = T, 
          line = 2)
}

# Plot stock-recruitment curve
plotSR <- function( repList = reports )
{
  report    <- repList$repOpt 
  initYear  <- repList$fYear

  # Pull stuff from report
  Rpt       <- report$R_pt
  Bpt       <- report$SB_pt
  omegaRtp  <- report$omegaR_tp
  sigmaR    <- report$sigmaR
  R0_p      <- report$R0_p
  B0_p      <- report$B0_p
  reca_p    <- report$reca_p
  recb_p    <- report$recb_p
  h_p       <- round(report$rSteepness_p,2)
  h         <- round(report$rSteepness,2)
  phi_p     <- report$phi_p
  
  # Get model dimensions
  nT      <- report$nT
  nP      <- report$nP

  # Get max depletion level
  D_pt    <- Bpt
  for( p in 1:nP )
    D_pt[p,] <- D_pt[p,] / B0_p[p]

  # stock names
  stockNames <- dimnames(Rpt)[[1]]

  cols    <- brewer.pal(n = max(nP,3), "Set1") 

  SB <- seq(0,1.2*max(Bpt,B0_p),length = 1000 )
  # Recruitment
  R_p <-  matrix(0, nrow = nP, ncol = length(SB) )
  for( p in 1:nP) 
    R_p[p,] <- reca_p[p] * SB / (1 + recb_p[p]*SB)
  # 20% level
  B20_p <- 0.2*B0_p
  R20_p <- reca_p * B20_p / (1 + recb_p*B20_p)

  par( mfrow = c(nP,1), mar = c(1,1,1.5,1), oma = c(3,3,1,1))

  for( p in 1:nP )
  {
    plot( x = c(0,max(Bpt[p,],B0_p[p])), y = c(0,1.2 * max(R_p[p,],Rpt[p,],na.rm = T ) ),
          type = "n", las = 1, xlab = "",
          ylab = "")
      mfg <- par("mfg")
      lines( x = SB, y = R_p[p,], lwd = 3, col = cols[p] )
      points( x = Bpt[p,1:nT], y = Rpt[p,2:(nT+1)], pch = 16, col = cols[p] )
      mtext( side = 3, text = stockNames[p], font = 2)
      # Plot B0,R0
      segments( x0 = B0_p[p], y0 = 0, y1 = R0_p[p],
                lty = 2, lwd = 2, col = cols[p] )
      segments( x0 = 0, x1 = B0_p[p], y0 = R0_p[p], y1 = R0_p[p],
                lty = 2, lwd = 2, col = cols[p] )
      # Plot B20,R20
      segments( x0 = B20_p[p], x1 = B20_p[p], y0 = 0, y1 = R20_p[p],
                lty = 2, lwd = 2, col = cols[p] )
      segments( x0 = 0, x1 = B20_p[p], y0 = R20_p[p], y1 = R20_p[p],
                lty = 2, lwd = 2, col = cols[p] )
      # Label with steepness
      text( x = 1.1, y = R0_p[p]*1.3, labels = paste("h = ", h_p[p], sep = "") )
  }
  


  mtext( side = 1, text = "Spawning Stock Biomass", line = 1.5, outer = TRUE )
  mtext( side = 2, text = "Age-1 Recruits (1e6)", line = 2, outer =  TRUE )

}

# recruitments
plotRt <- function( repList = reports,
                    noPar = FALSE, 
                    pIdx = 1, 
                    labcex = .8,
                    plotLog = FALSE ,
                    heading = "dimLabs" )
{
  report    <- repList$repOpt 
  initYear  <- repList$fYear

  # Pull stuff from report
  Rt      <- report$R_pt[pIdx,]
  omegaRt <- report$omegaR_pt[pIdx,]
  sigmaR  <- report$sigmaR
  R0      <- report$R0_p[pIdx]

  ageData_axgt <- repList$data$A_axpgt[,,pIdx,,]
  ageData_axgt[ageData_axgt<0] <- NA
  ageData_agt <- apply(X = ageData_axgt, FUN = sum, MARGIN = c(1,3,4),na.rm = T)

  colRec_t <- rep("grey40", length(Rt))

  checkPos <- function( x )
  {
    if(any(x > 0))
      return(1)
    else return(0)
  }

  posAgeIdx_t <- apply( X = ageData_agt, FUN = checkPos, MARGIN = 3 )

  colRec_t[posAgeIdx_t == 0] <- "white"
  


  plotCI <- FALSE

  if(!is.null(repList$posts))
  {
    plotCI <- TRUE
    tInitModel_p <- repList$repOpt$tInitModel_p
    R_it <- repList$posts$R_ipt[,pIdx,]

    R_qt <- apply(X = R_it, FUN = quantile,
                    MARGIN = c(2), probs = c(0.025, 0.5, 0.975),
                    na.rm = TRUE )

    Rt <- apply(X = R_it, FUN = mean, MARGIN = c(2))

    R0 <- round(mean(repList$posts$R0_ip[,pIdx]),2)

  } else {
    R_qt <- 0
  }

  R0lab <- round(R0,2)

  Rt[Rt == 0] <- NA
  R_qt[R_qt == 0] <- NA

  if(plotLog)
  {
    R0 <- log(R0, base = 10)
    R_qt <- log(R_qt, base = 10)
    Rt <- log(Rt, base = 10)

  }

  # stock labels
  stockLabs <- dimnames(report$R_pt)[[1]]


  if( !is.null(heading))
  {
    if( heading == "dimLabs")
      stockLabel <- stockLabs
    else
      stockLabel <- try(heading[stockLabs])

    if( class(stockLabel) == "try-error" ) 
      stockLabel <- stockLabs
  }

  # Get number of time steps
  nT    <- report$nT
  years <- seq(from = initYear, length = nT+1, by = 1)

  vertLines <- seq(from = 1950, to = max(years), by = 10)

  # Set up plotting area
  if(!noPar)
    par(mfrow = c(1,1), mar = c(.5,3,.5,3), oma = c(3,1,3,1) )

  # Now plot recruitments
  plot( x = range(years), y = c(0,max(Rt,R_qt, na.rm =T) ),
        type = "n", xlab = "", ylab = "",
        las = 1, axes = FALSE )
    mfg <- par("mfg")
    if(plotLog)
    {
      maxY <- ceiling(max(Rt,R_qt, na.rm =T))
      yTicks <- 0:maxY
      yLabs <- 10^yTicks
    }
    if( mfg[1] == mfg[3])
      axis(side = 1)
    if(!plotLog)
      axis(side = 2, las = 1  )
    if(plotLog)
      axis(side = 2, las = 1, at = yTicks, labels = yLabs )
    box()
    grid()
    # Plot recruitment
    # lines( x = years[1:nT], y = Rt[1:nT], lwd = 2, col = "grey40" )
    abline( h = R0, lty = 2, lwd = 2)
    if( plotCI )
      segments( x0 = years[1:nT],
                y0 = R_qt[1,1:nT],
                y1 = R_qt[3,1:nT],
                lwd = 1, col = "grey40" )
    points( x = years[1:nT], y = Rt[1:nT], lwd = 2, bg = "grey40",
            pch = 21, col = "grey40" )

    axis( side = 3, at = years[posAgeIdx_t == 0], labels = FALSE,
          col.ticks = "grey40", tck = .02, lwd.ticks = 3 )
    
    panLab( x = 0.8, y = 0.1, txt = paste("R0 = ", R0lab, sep = "") )
    panLab( x = 0.8, y = 0.2, txt = "|  Miss Ages")
    if( mfg[1] == mfg[3])
      mtext( side = 1, text = "Year", line = 2)
    if( mfg[2] == 1 )
      mtext(side = 2, text = "Recruits\n(millions)", line = 3.5, cex = labcex)
    if( mfg[1] == 1 & !is.null(heading) )
      mtext(side = 3, text = stockLabel[pIdx], line = 1, font = 2, cex = labcex )

}

# recruitments, adjusted for brood year
plotRtResids <- function( repList = reports,
                          noPar = FALSE, 
                          pIdx = 1, labcex = .8 )
{
  report    <- repList$repOpt 
  initYear  <- repList$fYear

  # Pull stuff from report
  Rt      <- report$R_pt[pIdx,]
  omegaRt <- report$SRdevs_pt[pIdx,]
  sigmaR  <- report$sigmaR
  R0      <- report$R0_p[pIdx]

  omegaRt <- omegaRt / sigmaR

  # Get number of time steps
  nT    <- report$nT
  years <- seq(from = initYear, length = nT+1, by = 1)

  # Adjust by brood year
  vertLines <- seq(from = 1950, to = max(years), by = 10)

  # Set up plotting area
  if(!noPar)
    par(mfrow = c(1,1), mar = c(.5,3,.5,3), oma = c(3,1,3,1) )

  # Now plot recruitment resids
  plot( x = range(years), y = range(omegaRt,omegaRt-sigmaR^2/2, na.rm =T),
        type = "n", xlab = "", ylab = "", axes = FALSE,
        las = 1 )
    # Add axis label if on the left margin
    mfg <- par("mfg")
    if( mfg[2] == 1 )
      mtext(side = 2, text = "Recruitment std.\nlog-residuals", line = 3, cex = labcex)  
    if( mfg[1] == mfg[3] )
    {
      axis( side = 1 )
      mtext( side = 1, text = "Year", line = 2 )
    }
    # Add y axis
    axis( side = 2, las = 1 )
    box()
    # Plot recruitment - update to include SEs later
    abline( h = 0, lty = 2, lwd = 1)
    grid()
    abline( h = mean(omegaRt,na.rm =T), lty = 3, lwd = 2, col = "red")
    points( x = years[1:(nT)], y = omegaRt[1:nT], col = "grey40", pch = 16 )
    panLab( x = 0.8, y = 0.8, txt = paste("sigmaR = ", round(sigmaR,2), sep = "" ) )

    
    panLegend(  x = 0.1, y = 0.95, 
                bty = "n",
                legTxt = c("Mean Resid"),
                lty = c(3),
                lwd = c(2),
                col = c("red") )
}

# # Plot age fits
# plotAgeFitYrs <- function(  report = repFE,
#                             initYear = fYear,
#                             gearLabels = gearLabs,
#                             pIdx = 1 )
# {
#   # Pull predicted and observed ages
#   predAge <- report$predPA_apgt[,,,pIdx]
#   obsAge  <- report$A_apgt[,,,pIdx]

#   # Pull model dims
#   nG      <- report$nG
#   nT      <- report$nT
#   nA      <- report$nA

  

#   # Make years vector
#   years   <- seq(initYear, length = nT+1, by = 1)

#   # Now, we want to loop over gear types now, 
#   # and record the gIdxes for which there are
#   # observations
#   ageGears <- c()
#   gearTimes <- vector(mode = "list", length = nG)
#   for( gIdx in 1:nG )
#   {
#     if( any(obsAge[1,,gIdx] >= 0) )
#     {
#       ageGears <- c(ageGears,gIdx)
#       gearTimes[[gIdx]] <- which(obsAge[1,,gIdx] >= 0)
#     }
#   }

#   # Make colours vector
#   cols    <- brewer.pal( n = length(ageGears), "Paired" )

#   # ok, ageGears are the ones we want to plot,
#   # and gearTimes is the list of time indices
#   for( aIdx in 1:length(ageGears) )
#   {
#     gIdx <- ageGears[aIdx]
#     dev.new()
#     times <- gearTimes[[gIdx]]
#     # Count the number of age observations
#     # there are, and make the plotting window
#     nObs <- length(times)
#     nCols <- round(sqrt(nObs))
#     nRows <- ceiling(nObs/nCols)

#     par(  mfcol = c(nRows,nCols), 
#           mar = c(1,1,1,1),
#           oma = c(3,3,3,3) )

#     ageObs <- obsAge[,,gIdx]
#     # ageObs <- sweep( x = ageObs, FUN = "/", MARGIN = c(2), STATS = "sum")

#     for( tIdx in times )
#     { 
#       # get age obs and preds
#       ageObs.t  <- ageObs[,tIdx]/sum(ageObs[,tIdx])
#       agePred <- predAge[,tIdx,gIdx]

#       plot( x = c(1,nA), y = c(0,max(ageObs.t,agePred,na.rm = T) ),
#             xlab = "", ylab = "", type = "n", las = 1 )
#         rect( xleft = 1:nA - .3, xright = 1:nA + .3,
#               ybottom = 0, ytop = ageObs.t,
#               col = "grey40", border = NA )
#         lines(  x = 1:nA, y = agePred, lwd = 1,
#                 col = cols[aIdx] )
#         points(  x = 1:nA, y = agePred,
#                 col = cols[aIdx], pch = 21 )
#         panLab( x=.5, y = .95, txt = years[tIdx] )

#     }
#     mtext( side = 1, outer = T, text = "Age", line = 2 )
#     mtext( side = 2, outer = T, text = "Proportion", line = 2 )
#   }

# }

# Plot age fits averaged over time and stock
plotAgeFitAggAvg <- function( repList = reports )
{
  report    <- repList$repOpt 
  initYear  <- repList$fYear

  # Pull predicted and observed ages
  predAge <- report$predPA_apgt
  obsAge  <- report$A_apgt

  # Pull mixed ages, add these to the area-specific
  # ones
  predAgeMixed <- report$predMixedPA_agt
  obsAgeMixed  <- report$mA_agt

  predAgeMixed[predAgeMixed < 0] <- NA
  obsAgeMixed[obsAgeMixed < 0] <- NA

  minAge_g <- repList$data$minAge_g

  gearLabs  <- dimnames(obsAge)[[3]]

  # replace missing entries with NAs
  predAge[predAge < 0] <- NA
  obsAge[obsAge < 0] <- NA

  # Get time steps
  nT <- dim(predAge)[4]
  nG <- dim(predAge)[3]
  nP <- dim(predAge)[2]

  # Multiply the yearly predicted proportions by the
  # yearly sample sizes
  for( g in 1:nG )
    for( p in 1:nP )
      for( t in 1:nT )
      {
        thisSamp <- sum( obsAge[,p,g,t], na.rm = T )
        if( thisSamp < 0 )
          thisSamp <- 0

        predAge[,p,g,t] <- thisSamp * predAge[,p,g,t]
      }

  
  
  # Average over time
  predAge <- apply( X = predAge, FUN = sum, MARGIN = c(1,3), na.rm = T )
  obsAge <- apply( X = obsAge, FUN = sum, MARGIN = c(1,3), na.rm = T )    

  obsAgeMixed <- apply(X = obsAgeMixed, FUN = sum, MARGIN = c(1,2), na.rm = T )
  predAgeMixed <- apply(X = predAgeMixed, FUN = sum, MARGIN = c(1,2), na.rm = T )

  predAge <- predAge + predAgeMixed
  obsAge <- obsAge + obsAgeMixed


  
  # Pull model dims
  nG      <- report$nG
  nA      <- report$nA
  nT      <- report$nT
  nP      <- report$nP
  minPA   <- report$minPropAge


  # Now, we want to loop over gear types now, 
  # and record the gIdxes for which there are
  # observations
  ageGears <- c()
  gearTimes <- vector(mode = "list", length = nG)
  for( gIdx in 1:nG )
    if( sum(obsAge[,gIdx]) > 0  )
      ageGears <- c(ageGears,gIdx)

  # Make colours vector
  cols    <- brewer.pal( n = nG, "Paired" )


  par(  mfcol = c(length(ageGears),1), 
        mar = c(1,2,1,2),
        oma = c(3,3,3,3) )

  # ok, ageGears are the ones we want to plot,
  for( aIdx in 1:length(ageGears) )
  { 
    gIdx <- ageGears[aIdx]
    # get age obs and preds
    ageObs  <- obsAge[,gIdx]
    nSamp   <- sum(ageObs)
    ageObs  <- ageObs / sum(ageObs)
    agePred <- predAge[,gIdx]
    agePred <- agePred / sum( agePred )

    plot( x = c(1,nA), y = c(0,1.1*max(ageObs,agePred,na.rm = T) ),
          xlab = "", ylab = "", type = "n", las = 1 )
      legend( x = "topright",
              bty = "n",
              legend = paste("N = ", nSamp, sep = "" ) )
      rect( xleft = minAge_g[gIdx]:nA - .3, xright = minAge_g[gIdx]:nA + .3,
            ybottom = 0, ytop = ageObs[minAge_g[gIdx]:nA],
            col = "grey40", border = NA )
      lines(  x = minAge_g[gIdx]:nA, y = agePred[minAge_g[gIdx]:nA], lwd = 3,
              col = cols[gIdx] )
      points(  x = minAge_g[gIdx]:nA, y = agePred[minAge_g[gIdx]:nA],
              col = cols[gIdx], pch = 16, cex = 1.5 )
      abline( h = minPA, lty = 2, lwd = .8 )
      
      
      # panLab( x=.5, y = .95, txt = gearLabels[gIdx] )
      mfg <- par("mfg")
      if(mfg[1] == 1 )
        mtext( side = 3, text = "Aggregate SAR", line = 1, font = 2)

      if( mfg[2] == mfg[4] )
        rmtext( txt = gearLabs[gIdx], line = 0.01, outer = TRUE, cex = 1.5,
                font = 2)
      
  }
  
  mtext( side = 1, outer = T, text = "Age", line = 2 )
  mtext( side = 2, outer = T, text = "Proportion", line = 2 )

}

# Plot age fits averaged over time
plotAgeFitAvg <- function(  repList = reports )
{
  report    <- repList$repOpt 
  initYear  <- repList$fYear
  gearLabels<- repList$gearLabs

  # Pull predicted and observed ages
  predAge <- report$predPA_apgt
  obsAge  <- report$A_apgt


  minAge_g <- repList$data$minAge_g

  gearLabs  <- dimnames(obsAge)[[3]]
  stockLabs <- dimnames(obsAge)[[2]]


  # replace missing entries with NAs
  predAge[predAge < 0] <- NA
  obsAge[obsAge < 0] <- NA


  # Get time steps
  nT <- dim(predAge)[4]
  nG <- dim(predAge)[3]
  nP <- dim(predAge)[2]

  # Multiply the yearly predicted proportions by the
  # yearly sample sizes
  for( g in 1:nG )
    for( p in 1:nP )
      for( t in 1:nT )
      {
        thisSamp <- sum( obsAge[,p,g,t], na.rm = T )
        if( thisSamp < 0 )
          thisSamp <- 0

        predAge[,p,g,t] <- thisSamp * predAge[,p,g,t]
      }


  
  # Average over time
  predAge <- apply( X = predAge, FUN = sum, MARGIN = c(1,2,3), na.rm = T )
  obsAge <- apply( X = obsAge, FUN = sum, MARGIN = c(1,2,3), na.rm = T )    

  # Pull model dims
  nG      <- report$nG
  nA      <- report$nA
  nT      <- report$nT
  nP      <- report$nP
  minPA   <- report$minPropAge

  # Make years vector
  years   <- seq(initYear, length = nT+1, by = 1)

  # Now, we want to loop over gear types now, 
  # and record the gIdxes for which there are
  # observations
  ageGears <- c()
  gearTimes <- vector(mode = "list", length = nG)
  for( gIdx in 1:nG )
    if( sum(obsAge[,,gIdx]) > 0  )
      ageGears <- c(ageGears,gIdx)

  # Make colours vector
  cols    <- brewer.pal( n = nG, "Paired" )


  par(  mfcol = c(length(ageGears),nP), 
        mar = c(1,2,1,2),
        oma = c(3,3,3,3) )

  # ok, ageGears are the ones we want to plot,
  for( pIdx in 1:nP)
    for( aIdx in 1:length(ageGears) )
    { 
      gIdx <- ageGears[aIdx]
      # get age obs and preds
      ageObs  <- obsAge[,pIdx,gIdx]
      nSamp   <- sum(ageObs)
      ageObs  <- ageObs / sum(ageObs)
      agePred <- predAge[,pIdx,gIdx]
      agePred <- agePred / sum( agePred )

      plot( x = c(1,nA), y = c(0,1.1*max(ageObs,agePred,na.rm = T) ),
            xlab = "", ylab = "", type = "n", las = 1 )
        legend( x = "topright",
                bty = "n",
                legend = paste("N = ", nSamp, sep = "" ) )
        rect( xleft = minAge_g[gIdx]:nA - .3, xright = minAge_g[gIdx]:nA + .3,
              ybottom = 0, ytop = ageObs[minAge_g[gIdx]:nA],
              col = "grey40", border = NA )
        lines(  x = minAge_g[gIdx]:nA, y = agePred[minAge_g[gIdx]:nA], lwd = 3,
                col = cols[gIdx] )
        points(  x = minAge_g[gIdx]:nA, y = agePred[minAge_g[gIdx]:nA],
                col = cols[gIdx], pch = 16, cex = 1.5 )
        abline( h = minPA, lty = 2, lwd = .8 )
        
        
        # panLab( x=.5, y = .95, txt = gearLabels[gIdx] )
        mfg <- par("mfg")
        if(mfg[1] == 1 )
          mtext( side = 3, text = stockLabs[pIdx], line = 1, font = 2)

        if( mfg[2] == mfg[4] )
        {
          corners <- par("usr") #Gets the four corners of plot area (x1, x2, y1, y2)
          par(xpd = TRUE) #Draw outside plot area
          text( x = corners[2]+0.5, 
                y = mean(corners[3:4]), 
                labels = gearLabs[gIdx], srt = 270,
                font = 2, cex = 1.5 )
          par(xpd = FALSE)
        }
    }
  
  mtext( side = 1, outer = T, text = "Age", line = 2 )
  mtext( side = 2, outer = T, text = "Proportion", line = 2 )

}

# Plot length fits averaged over time
plotLengthFitAvg <- function(  repList = reports )
{
  report    <- repList$repOpt 
  initYear  <- repList$fYear
  gearLabels<- repList$gearLabs

  # Pull predicted and observed ages
  pLen_lxpgt <- report$predPL_lxpgt
  oLen_lxpgt <- report$L_lxpgt

  minLenBin_g <- repList$data$minLenBin_g
  maxLenBin_g <- repList$data$maxLenBin_g

  gearLabs  <- dimnames(oLen_lxpgt)[[4]]
  stockLabs <- dimnames(oLen_lxpgt)[[3]]
  lenBins_l <- dimnames(oLen_lxpgt)[[1]]
  sexLabs_x <- dimnames(oLen_lxpgt)[[2]]

  # replace missing entries with NAs
  pLen_lxpgt[pLen_lxpgt < 0] <- NA
  oLen_lxpgt[oLen_lxpgt < 0] <- NA


  # Get time steps
  nT    <- dim(pLen_lxpgt)[5]
  nG    <- dim(pLen_lxpgt)[4]
  nP    <- dim(pLen_lxpgt)[3]
  nLenX <- dim(pLen_lxpgt)[2]
  nL    <- dim(pLen_lxpgt)[1]

  # Multiply the yearly predicted proportions by the
  # yearly sample sizes
  for( x in 1:nLenX)
    for( g in 1:nG )
      for( p in 1:nP )
        for( t in 1:nT )
        {
          thisSamp <- sum( oLen_lxpgt[,x,p,g,t], na.rm = T )
          if( thisSamp < 0 )
            thisSamp <- 0

          pLen_lxpgt[,x,p,g,t] <- thisSamp * pLen_lxpgt[,x,p,g,t]
        }


  
  # Average over time
  pLen_lxpg <- apply( X = pLen_lxpgt, FUN = sum, MARGIN = c(1:4), na.rm = T )
  oLen_lxpg <- apply( X = oLen_lxpgt, FUN = sum, MARGIN = c(1:4), na.rm = T )    

  # Pull model dims
  nG      <- report$nG
  nX      <- report$nX
  nT      <- report$nT
  nP      <- report$nP
  minPL   <- report$minPropLen

  # Make years vector
  years   <- seq(initYear, length = nT+1, by = 1)

  # Now, we want to loop over gear types now, 
  # and record the gIdxes for which there are
  # observations
  lenGears <- c()
  gearTimes <- vector(mode = "list", length = nG)
  for( gIdx in 1:nG )
    if( sum(oLen_lxpg[,,,gIdx]) > 0  )
      lenGears <- c(lenGears,gIdx)

  # Make colours vector
  cols    <- brewer.pal( n = nG, "Paired" )

  par(  mfcol = c(length(lenGears),nP * nLenX), 
        mar = c(.1,2,.1,2),
        oma = c(3,3,3,3) )

  xAxisTicks <- seq(from = 1, to = nL, by = 10)
  xAxisLabs  <- lenBins_l[xAxisTicks]

  # ok, ageGears are the ones we want to plot,
  for( pIdx in 1:nP)
    for( x in 1:nLenX)
    for( lIdx in 1:length(lenGears) )
      { 
        gIdx <- lenGears[lIdx]
        # get age obs and preds
        oLen_l <- oLen_lxpg[,x,pIdx,gIdx]
        nSamp  <- sum(oLen_l)
        pLen_l <- pLen_lxpg[,x,pIdx,gIdx]
        oLen_l <- oLen_l/nSamp
        pLen_l <- pLen_l/sum(pLen_l)
        
        plot( x = c(1,nL), y = c(0,1.1*max(oLen_l,pLen_l,na.rm = T) ),
            xlab = "", ylab = "", type = "n", las = 1,
            axes = FALSE )
          mfg <- par("mfg")
          axis( side = 2, las = 1)
          if(mfg[1] == mfg[3])
            axis(side = 1, at = xAxisTicks, labels = xAxisLabs )
          box()
          legend( x = "topright",
                  bty = "n",
                  legend = paste("N = ", round(nSamp), sep = "" ) )
            lenBinIdx <- minLenBin_g[gIdx]:maxLenBin_g[gIdx]

          rect( xleft = lenBinIdx - .3, xright = lenBinIdx + .3,
                ybottom = 0, ytop = oLen_l[lenBinIdx],
                col = "steelblue", border = NA )  
          lines(  x = lenBinIdx, y = pLen_l[lenBinIdx], lwd = 3, col = cols[gIdx] )
          points( x = lenBinIdx, y = pLen_l[lenBinIdx],
                  col = cols[gIdx], pch = 16, cex = 1.5 )
          abline( h = minPL, lty = 2, lwd = .8 )

          
          if(mfg[1] == 1 )
            mtext( side = 3, text = paste0(stockLabs[pIdx],":",sexLabs_x[x]), line = 1, font = 2)


          if( mfg[2] == mfg[4] )
          {
            corners <- par("usr") #Gets the four corners of plot area (x1, x2, y1, y2)
            par(xpd = TRUE) #Draw outside plot area
            text( x = corners[2]+4, 
                  y = mean(corners[3:4]), 
                  labels = gearLabs[gIdx], srt = 270,
                  font = 2, cex = 1 )
            par(xpd = FALSE)
          }

        }
  
  mtext( side = 1, outer = T, text = "Length", line = 2 )
  mtext( side = 2, outer = T, text = "Proportion", line = 2 )

}

# Plot age fits averaged over time
plotMixedAgeFitAvg <- function( repList = reports )
{
  report    <- repList$repOpt 
  initYear  <- repList$fYear
  gearLabels<- repList$gearLabs

  # Pull predicted and observed ages
  predAge <- report$predMixedPA_agt
  obsAge  <- report$mA_agt

  gearLabs  <- dimnames(obsAge)[[2]]

  # replace missing entries with NAs
  predAge[predAge < 0] <- NA
  obsAge[obsAge < 0] <- NA
  
  # Average over time
  predAge <- apply( X = predAge, FUN = mean, MARGIN = c(1,2), na.rm = T )
  obsAge <- apply( X = obsAge, FUN = mean, MARGIN = c(1,2), na.rm = T )    

  # Pull model dims
  nG      <- report$nG
  nA      <- report$nA
  nT      <- report$nT
  minPA   <- report$minPropAge

  # Make years vector
  years   <- seq(initYear, length = nT+1, by = 1)

  # Now, we want to loop over gear types now, 
  # and record the gIdxes for which there are
  # observations
  ageGears <- c()
  gearTimes <- vector(mode = "list", length = nG)
  for( gIdx in 1:nG )
    if( any(!is.na(obsAge[1,gIdx]) ) )
      ageGears <- c(ageGears,gIdx)

  # Make colours vector
  cols    <- brewer.pal( n = nG, "Paired" )


  par(  mfcol = c(length(ageGears),1), 
        mar = c(1,1,1,1),
        oma = c(3,3,3,3) )

  # ok, ageGears are the ones we want to plot,
  for( aIdx in 1:length(ageGears) )
  { 
    gIdx <- ageGears[aIdx]
    # get age obs and preds
    ageObs  <- obsAge[,gIdx]
    ageObs  <- ageObs / sum(ageObs)
    agePred <- predAge[,gIdx]

    plot( x = c(1,nA), y = c(0,max(ageObs,agePred,na.rm = T) ),
          xlab = "", ylab = "", type = "n", las = 1 )
      rect( xleft = 1:nA - .3, xright = 1:nA + .3,
            ybottom = 0, ytop = ageObs,
            col = "grey40", border = NA )
      lines(  x = 1:nA, y = agePred, lwd = 3,
              col = cols[aIdx] )
      points(  x = 1:nA, y = agePred,
              col = cols[aIdx], pch = 16, cex = 1.5 )
      abline( h = minPA, lty = 2, lwd = .8 )
      panLab( x=.5, y = .95, txt = gearLabels[gIdx] )
      mfg <- par("mfg")
      if(mfg[1] == 1)
        mtext( side = 3, text = "Mixed stock", line = 1, font = 2)
  }
  
  mtext( side = 1, outer = T, text = "Age", line = 2 )
  mtext( side = 2, outer = T, text = "Proportion", line = 2 )

}



# Plot Biomass, catch and indices on the same
# axes
plotItResids <- function( repList =reports,
                          noPar = FALSE,
                          pIdx = 1, labcex = .8 )
{
  report    <- repList$repOpt 
  initYear  <- repList$fYear
  gearLabels<- repList$gearLabs
  # Pull stuff from report
  Bt     <- report$SB_pt[pIdx,]
  Igt    <- report$I_pgt[pIdx,,]
  qg     <- report$qhat_pg[pIdx,]
  rIgt   <- repList$data$rI_pgt[pIdx,,]

  # Need to pull vuln biomass/numbers
  vulnBgt <- report$vulnB_pgt[pIdx,,]
  vulnNgt <- report$vulnN_pgt[pIdx,,]

  # Pull probability of positive indices
  # probPosIdx_gt <- report$probPosIdx_pgt[pIdx,,]

  # Get survey type and whether
  # idx is relative or absolute
  surveyType  <- report$survType_g      # abs or relative?
  indexType   <- report$indexType_g     # biomass or numbers
  calcIndex   <- report$calcIndex_g     # calc index at all?

  # Model dimensions
  nG          <- report$nG
  nT          <- report$nT
  surveyGears <- which(calcIndex == 1)

  # Get survey names
  surveyNames <- dimnames(vulnBgt)[[1]][surveyGears]
  stockLabel  <- dimnames(report$I_pgt)[[1]]

  if( stockLabel == "Aggregate" )
  {
    # We need to pull the aggregated stock
    # indices in at this point
    mIgt    <- report$mI_gt
    mqg     <- report$qhat_g

    mIgt[mIgt < 0] <- 0

    for( g in surveyGears )
    {
      if( all(Igt[g,] < 0) & any(mIgt[g,] > 0))
      {
        Igt[g,] <- mIgt[g,]
        qg[g]   <- mqg[g]
      }
    }
  }

  Igt[Igt<=0] <- NA

  if(!is.null(repList$posts))
  {
    tInitModel_p <- repList$repOpt$tInitModel_p
    SB_it <- repList$posts$SB_ipt[,pIdx,]

    Bt <- apply(X = SB_it, FUN = mean, MARGIN = c(2))

    q_ig  <- repList$posts$q_ipg[,pIdx,]
    qg   <- apply( X = q_ig, FUN = mean, MARGIN = 2)

    # probPosIdx_igt <- repList$posts$probPosIdx_ipgt[,pIdx,,]
    # probPosIdx_gt  <- apply( X = probPosIdx_igt, FUN = mean, MARGIN = c(2,3))

  } 


  # Scale indices by q
  IgtScaled <- Igt
  for( g in 1:nG)
  {
    IgtScaled[g,] <- Igt[g,] / qg[g] 
  }

  IgtScaled[IgtScaled < 0] <- NA

  cols <- brewer.pal(nG, name = "Paired")

  # Get number of time steps
  years <- seq(from = initYear, length = nT, by = 1)
  vertLines <- seq(from = 1950, to = max(years), by = 10)

  tauObs_g <- report$tauObs_pg[pIdx,]

  resids <- IgtScaled

  # Calc residuals
  for( g in 1:nG )
  {
    posIdx <- which(IgtScaled[g,] > 0)

    # Check index type, calcs resids
    if( surveyType[g] == 0 )
      resids[g,posIdx] <- log(IgtScaled[g,posIdx] / vulnBgt[g,posIdx] ) / tauObs_g[g]
    if( surveyType[g] == 1 )
      resids[g,posIdx] <- log(IgtScaled[g,posIdx] / vulnNgt[g,posIdx] ) / tauObs_g[g]
    if( surveyType[g] == 2 )
      resids[g,posIdx] <- log(IgtScaled[g,posIdx] / Bt[posIdx] ) / tauObs_g[g]

  }   

  combResids <- rep(NA,nT)
  # if( any(repList$data$whichCombIdx_g == 1) )
  #   combResids <- -1 * report$zComb_pt[pIdx,]


  # Now set up plotting window
  if( !noPar )
    par( mfrow = c(1, 1), mar = c(2,3,1,3), oma = c(3,3,1,1) )

  plot( x = range(years), y = range(resids,combResids,na.rm=T), xlab = "",
          ylab = "", type = "n", las = 1 )
    mfg <- par("mfg")
    # show 0 line
    abline( h = 0, lty = 3, lwd = .8 )
    abline( v = vertLines, lwd = .8, lty = 2, col = "grey80")
    # Now show the resids
    for(g in surveyGears)
    {
      residTimes <- which(!is.na(resids[g,]))
      if(length(residTimes) > 0)
      {
        points( x = years[residTimes], y = resids[g,residTimes], pch = 16,
                col = cols[g] )
        # Fit a regression
        dat <- data.frame(x = years[residTimes], y = resids[g,residTimes])
        residModel <- lm( y~x, data = dat )
        dat$pred <- predict.lm(object = residModel, newdata = dat)

        pVal <- round(summary(residModel)$coefficients[2,4],2)

        lines( x = dat$x, y = dat$pred,
                col = cols[g], lwd = 2 )
        text( x = dat$x[4], y = min(dat$y) + .3, label = paste("p = ", pVal, sep = ""),
              col = cols[g], font = 2 )

        meanResid <- mean(resids[g,residTimes])
        segments(  x0 = years[residTimes[1]],
                  x1 = years[residTimes[length(residTimes)]], 
                  y0 = meanResid, col = cols[g], lty = 2, lwd = 2 ) 
      }
    }
    combResids[combResids == 0] <- NA
    combResidTimes <- which(!is.na(combResids))

    if(length(combResidTimes) > 0)
    {
      points( x = years[combResidTimes], y = combResids[combResidTimes], pch = 16,
                col = "grey40" )
        # Fit a regression
        dat <- data.frame(x = years[combResidTimes], y = combResids[combResidTimes])
        residModel <- lm( y~x, data = dat )
        dat$pred <- predict.lm(object = residModel, newdata = dat)

        pVal <- round(summary(residModel)$coefficients[2,4],2)

        lines( x = dat$x, y = dat$pred,
                col = "grey40", lwd = 2 )
        text( x = dat$x[7], y = min(dat$y + .3), label = paste("p = ", pVal, sep = ""),
              col = "grey40", font = 2 )

        meanResid <- mean(combResids[combResidTimes])
        segments(  x0 = years[combResidTimes[1]],
                  x1 = years[combResidTimes[length(combResidTimes)]], 
                  y0 = meanResid, col = cols[g], lty = 2, lwd = 2 )

        surveyNames <- c(surveyNames,"Spawn Index")
    }
    if(mfg[2] == 1 )
      legend( x = "topleft",
              legend = surveyNames,
              pch = 16,
              col = c(cols[which(calcIndex == 1)],"grey40") )
    if( mfg[2] == 1)
      mtext(  side = 2, outer = F,
              text = paste("Std. log-residuals", sep = ""),
              line = 3, cex = labcex )
}

# Plot standardised residuals for each gear
plotCompResids <- function( repList = reports,
                            pIdx = 1 )
{
  repObj    <- repList$repOpt 
  initYear  <- repList$fYear
  gearLabels<- repList$gearLabs

  # Pull age data and predictions
  nL            <- repObj$nL
  resids_lxgt   <- repObj$lenResids_lxpgt[,,pIdx,,]
  yLab          <- "Length (cm)"
  tau2Len_xg    <- repObj$tau2Len_xpg[,pIdx,]
  lenBins_l     <- repList$data$lenBinMids_l

  yAxisTicks    <- seq(from = lenBins_l[1], to = lenBins_l[nL],
                        by = 10 )

  # Take sqrt for stdizing reisds
  tauLen_xg      <- sqrt(tau2Len_xg)

  
  # Pull gear names
  gearLabs  <- dimnames( repObj$predPL_lxpgt )[[4]]
  stockName <- dimnames( repObj$predPL_lxpgt )[[3]][pIdx]
  nG        <- dim(resids_lxgt)[3]
  nT        <- dim(resids_lxgt)[4]
  nLenX     <- dim(resids_lxgt)[2]

  yrs <- seq( initYear, by = 1, length.out = nT )

  lenGears  <- c()

  for( g in 1:nG )
    if( any( resids_lxgt[,,g,] != 0 ) )
      lenGears <- c(lenGears,g)

  # Count age gears
  nLenGears <- length(lenGears)

  resids_lxgt <- resids_lxgt/max(resids_lxgt,na.rm = T)
  

  par(  mfcol = c(nLenGears,nLenX), 
        mar = c(1,.5,1,.5),
        oma = c(3,3,3,2) )

  # resids_lgt <- resids_lgt / max(abs(range(resids_lgt)))

  for( x in 1:nLenX)
    for( g in lenGears )
    {
      # Plot window
      plot( x = range(yrs), y = c(0,nL+1),
            axes = FALSE, type = "n", 
            xlab = "", ylab = "" )
        grid()

        mfg <- par("mfg")
        if( mfg[2] == 1 )
          axis( side = 2, las = 1, at = yAxisTicks, labels = lenBins_l[yAxisTicks] )
        if( mfg[1] == mfg[3] )
          axis( side = 1 )
        
        
        box()


        # Gear label
        if(mfg[2] == mfg[4])
          rmtext( txt = gearLabs[g],
                  outer = TRUE, 
                  font = 2, line = 1 )

        # Loop over years, plot circles
        for( t in 1:nT )
        {
          if( any(!is.finite(resids_lxgt[,x,g,t])) )
            next

          if( all( resids_lxgt[,x,g,t] == 0 )  )
            next

          res_l <- resids_lxgt[,x,g,t]/tauLen_xg[x,g]

          cols  <- rep( "black", nL )
          cols[res_l < 0] <- "red"
          cols[res_l == 0 ] <- NA

          radians <- seq(0,2*pi, length.out = 100)
          for( l in 1:nL)
            if( res_l[l] != 0 )
              lines(  x = yrs[t] + abs(res_l[l]) * cos(radians),
                      y = l + abs(res_l[l]) * sin(radians),
                      col = cols[l], lwd = 1 )


        }
        legend( x = "topright", bty = "n",
                legend = paste("tau = ", round(tauLen_xg[x,g],2), sep = "") )


      
    }
  mtext( side = 2, text = "Length (cm)", outer = T, line = 2)
  mtext( side = 1, text = "Year", outer = T, line = 2)

} #END plotCompResids()

# Plot comp fits
plotCompFitYrs <- function( repList = reports,
                            comps = "len",
                            save = FALSE,
                            savePath = "plotFitYrs",
                            pIdx = 1,
                            xIdx = 1,
                            tc = FALSE )
{
  repObj    <- repList$repOpt 
  initYear  <- repList$fYear
  gearLabels<- repList$gearLabs
  sexLabs   <- dimnames(repList$data$L_lxpgt)[[2]]

  # Pull predicted and observed ages
  if( comps == "age" )
  {
    max       <- repObj$nA
    pred_xgt  <- repObj$predPA_apgt[1:max,pIdx,,]
    obs_xgt   <- repObj$A_apgt[1:max,pIdx,,]  
    xLab      <- "Age"
    yLab      <- "Age"

    xAxisTicks <- seq(from = 1, to = max, by = 5)
    xAxisLabs  <- xAxisTicks

    minX      <- repList$data$minAge

    # browser()

    if( tc )
    {
      pred_xgt  <- repObj$tcPred_apgt[1:max,pIdx,,]
      obs_xgt   <- repObj$tcComps_apgt[1:max,pIdx,,]
    }
  }

  # Pull predicted and observed ages
  if( comps == "len" )
  {
    max       <- repObj$nL
    pred_xgt  <- repObj$predPL_lxpgt[1:max,xIdx,pIdx,,]
    obs_xgt   <- repObj$L_lxpgt[1:max,xIdx,pIdx,,]  
    xLab      <- "Length (cm)"
    yLab      <- "Length"

    lenBins_l <- repList$data$lenBinMids_l

    xAxisTicks <- seq(from = 1, to = max, by = 5)
    xAxisLabs  <- lenBins_l[xAxisTicks]

    minX      <- 1

    # browser()
  }
  
  dimNames    <- dimnames(pred_xgt)
  compNames   <- dimNames[[1]]
  yearNames   <- dimNames[[3]]
  gearNames   <- dimNames[[2]]
  stockName   <- dimnames(repObj$L_lxpgt)[[3]][pIdx]

  # Pull model dims
  nG      <- repObj$nG
  nT      <- repObj$nT  
  minPA   <- repObj$minPropAge

  # Make colours vector
  cols    <- brewer.pal( n = nG, "Paired" )

  # Make years vector
  years   <- seq(initYear, length = nT+1, by = 1)

  # Now, we want to loop over gear types now, 
  # and record the gIdxes for which there are
  # observations
  obsGears <- c()
  gearTimes <- vector(mode = "list", length = nG)
  for( gIdx in 1:nG )
  {
    if( any(obs_xgt[1,gIdx,] >= 0) )
    {
      obsGears <- c(obsGears,gIdx)
      gearTimes[[gIdx]] <- which(obs_xgt[1,gIdx,] >= 0)
    }
  }

  # ok, obsGears are the ones we want to plot,
  # and gearTimes is the list of time indices
  for( gIdx in obsGears )
  {
    if(save)
    {
      gearPath <- paste(savePath,gearNames[gIdx],".png",sep = "")
      png(  gearPath, 
            width = 11, height = 8.5,
            units = "in", res = 300 )
    }

    times <- gearTimes[[gIdx]]
    # Count the number of age observations
    # there are, and make the plotting window
    nObs <- length(times)
    nCols <- round(sqrt(nObs))
    nRows <- ceiling(nObs/nCols)

    par(  mfcol = c(nRows,nCols), 
          mar = c(0,0.5,0,0.5),
          oma = c(3,3,3,3) )

    yMax <- max(pred_xgt, na.rm = T ) 

    for( tIdx in times )
    { 
      Nobs <- sum(obs_xgt[,gIdx,tIdx])
      # get age obs and preds
      compObsProp_x  <- obs_xgt[,gIdx,tIdx]/sum(obs_xgt[,gIdx,tIdx])
      compPred_x     <- pred_xgt[,gIdx,tIdx]



      plot( x = c(1,max), y = c(0,yMax ),
            xlab = "", ylab = "", type = "n", las = 1,
            axes = FALSE )
        box()
        grid()
        mfg <- par("mfg")
        if( mfg[1] == mfg[3] )
          axis( side = 1, at = xAxisTicks, labels = xAxisLabs )

        if( mfg[2] == 1 )
          axis( side = 2, las = 1 )

        rect( xleft = minX:max - .3, xright = minX:max + .3,
              ybottom = 0, ytop = compObsProp_x[minX:max],
              col = "grey40", border = NA )
        lines(  x = minX:max, y = compPred_x[minX:max], lwd = 1,
                col = cols[gIdx] )
        points(  x = minX:max, y = compPred_x[minX:max],
                col = cols[gIdx], pch = 16 )
        abline( h = minPA, lty = 2, lwd = .8 )
        panLab( x=.5, y = .95, txt = years[tIdx], font = 2 )
        legend( x = "topright",
                legend = paste("N = ", round(Nobs), sep = "" ),
                bty = "n" )

    }
    mtext(  side = 3, 
            outer = T, 
            text = paste0(stockName, ":", gearNames[gIdx],":",sexLabs[xIdx]), 
            line = 2, font = 2)
    mtext(  side = 1, outer = T, text = xLab, line = 2 )
    mtext(  side = 2, outer = T, text = paste("Proportion-at-",yLab,sep=""), 
            line = 2 )
    if(save)
      dev.off()
  }

} # END plotCompFitYrs()

# plotProbRelease()
# Probability-at-age of releasing fish
plotProbRelease <- function( repList = reports )
{
  repObj <- repList$repOpt

  nA <- repObj$nA
  nP <- repObj$nP
  nX <- repObj$nX
  nG <- repObj$nG

  gearLabs <- repList$gearLabs

  pRel_axg <- repObj$pRel_axg

  

  plot( x = c(1,nA), y = c(0,1), 
        xlab = "Age", ylab = "Probability of release",
        las = 1, type = "n" )
    for( x in 1:nX )
        lines(x = 1:nA, y = pRel_axg[,x,1],
              lty = x, lwd = 3, col = "grey35" )

    legend( x = "topright",
            col = "grey35",
            bty = "n",
            legend = c("Males","Females"),
            lwd = 3, lty = 1:2 )
} # END plotProbRelease

# plotRelCatch()
# Total released biomass
plotRelCatch <- function( repList = reports )
{
  repObj <- repList$repOpt

  nA <- repObj$nA
  nP <- repObj$nP
  nX <- repObj$nX
  nG <- repObj$nG
  nT <- repObj$nT

  gearLabs <- repList$gearLabs

  relC_pgt <- repObj$relC_pgt

  fYear <- repList$fYear
  yrs <- seq(from = fYear, by = 1, length.out = nT)
  gearCols <- RColorBrewer::brewer.pal(nG,"Paired")
  commGears <- which(repObj$fleetType_g > 0)

  nComm <- length(commGears)

  par(mfrow = c(nP,1),
      mar = c(.1,1,.1,1),
      oma = c(3,3,1,1) )
  for( p in 1:nP )
  {
    plot( x = range(yrs), y = range(relC_pgt[p,,]),
          type = "n", axes = FALSE )
      mfg <- par("mfg")
      if(mfg[1] == mfg[3] )
        axis( side = 1 )
      if(mfg[2] == 1 )
        axis( side = 2, las = 1 )
      box()
      for(g in commGears)  
      {
        if( all(relC_pgt[p,g,] == 0))
          next
        lines( x = yrs, y = relC_pgt[p,g,], lwd = 2, col = gearCols[g] )
      }

    if( p == 1 )
      legend( x = "topleft", bty = "n",
              col = c(gearCols[commGears]),
              lwd = 2,
              legend = c(gearLabs[commGears]))
  }

  mtext(side = 1, text = "Year", outer = TRUE, line = 2 )
  mtext(side = 2, text = "Released catch (kt)", outer = TRUE, line = 2 )

} # END plotRelCatch


# plotProbLenAge()
# Probabilities (roughly) for length-at-age distributions
# for each sex
plotProbLenAge <- function( repList = reports )
{
  repObj <- repList$repOpt

  nA <- repObj$nA
  nP <- repObj$nP
  nX <- repObj$nX
  nG <- repObj$nG
  nT <- repObj$nT
  
  # Prob of length-at-age
  probLenAge_alx  <- repObj$probLenAge_alx
  len_ax          <- repObj$len_ax
  sizeLim_g       <- repList$data$sizeLim_g

  lenBins_l <- repList$data$lenBinMids_l
  maxL      <- max(lenBins_l)

  sexCols <- c("steelblue","salmon")

  par(mfrow = c(nX,1), mar = c(.1,1,.1,1),
      oma = c(3,3,1,1) )
  for( x in 1:nX )
  {
    plot( x = c(0,maxL), y = c(0,max(probLenAge_alx)),
          type = "n", las = 1, axes = FALSE)
      mfg <- par("mfg")
      if(mfg[1] == mfg[3] )
        axis(side = 1 )
      axis(side = 2, las = 1 )
      box()
      abline(v = sizeLim_g, lty = 2, col = "grey40", lwd =1)
      for( a in 1:nA )
      {
         lines(x = lenBins_l, y = probLenAge_alx[a,,x],
               lwd = 1, col = sexCols[x] )

         text( x = len_ax[a,x], y = 1,
               labels = paste0(a), cex = .5 )

      }
      
  }
  mtext(side = 1, text = "Length (cm)", outer = TRUE, line= 2)
  mtext(side = 2, text = "Probability", outer = TRUE, line= 2)


 
} # END plotPobLenAge()

plotCompLikelihoodError <- function(  repList = reports)
{
  repObj    <- repList$repOpt 
  initYear  <- repList$fYear
  gearLabels<- repList$gearLabs
  # Pull tail compressed compositions and
  # predictions
  tcComps_apgt  <- repObj$tcComps_apgt
  tcPred_apgt   <- repObj$tcPred_apgt 

  nP <- repObj$nP
  nG <- repObj$nG
  nT <- repObj$nT

  stockLabs <- dimnames(tcComps_apgt)[[2]]
  gearLabs <- dimnames(tcComps_apgt)[[3]]

  # Pull annual weights and geometric means
  gmObs_pgt   <- repObj$gmObs_pgt
  gmPred_pgt  <- repObj$gmPred_pgt
  ageWt_pgt   <- repObj$ageWt_pgt

  # Pull correlation matrices
  Corr_gaa    <- repObj$Corr_gaa

  # Pull internally calculated component
  intrnlAgeLikCpt_pg <- repObj$intrnlAgeLikCpt_pg

  # Pull etaSumSq
  modelEtaSumSq_pg  <- repObj$etaSumSq_pg

  # Pull model tau2Age
  tau2Age_pg        <- repObj$tau2Age_pg

  # And full likelihood function value
  ageObsNLL_pg <- repObj$ageObsNLL_pg

  # Calculate mean sample size
  meanA_apgt <- repObj$tcComps_apgt
  meanA_apgt[meanA_apgt < 0] <- NA
  meanA_pgt <- apply( X = meanA_apgt, FUN = sum, MARGIN = c(2,3,4), na.rm = T)
  meanA_pgt[meanA_pgt == 0] <- NA
  meanA_pg <- apply( X = meanA_pgt, FUN = mean, MARGIN = c(1,2), na.rm = T )
  meanA_pg[is.nan(meanA_pg)] <- 0

  # Now, we need to loop over stock, gear and time,
  # plot the errors in the annual weights, 
  # and annotate with the difference between
  # the etaSumSq values and the internally
  # calculated likelihood component.
  wtErr_pgt     <- array( NA, dim = c(nP, nG, nT) )  # On the log scale
  nBins_pgt     <- array( 0,  dim = c(nP, nG, nT) )  # number of bins
  gmObsErr_pgt  <- array( NA, dim = c(nP, nG, nT) )  # error in geometric mean calc for obs
  gmPredErr_pgt <- array( NA, dim = c(nP, nG, nT) )  # error in geometric mean calc for preds
  etaSumSq_pgt  <- array( 0,  dim = c(nP, nG, nT) )  # sum of squared resids
  intNLLCpt_pgt <- array( 0,  dim = c(nP, nG, nT) )  # Internal likelihood component

  # Differences in likelihood values
  diffIntCpt_pg   <- array( 0, dim = c(nP, nG) )
  diffEtaSumSq_pg <- array( 0, dim = c(nP, nG) )
  diffTotLike_pg  <- array( 0, dim = c(nP, nG) )
  tau_pg          <- array( 0, dim = c(nP, nG) )
  difftau_pg      <- array( 0, dim = c(nP, nG) )


  for( g in 1:nG )
    for( p in 1:nP )
    {
      if( any( tcComps_apgt[,p,g,] > 0 ) )
      {
        for( t in 1:nT )
        {
          if( all(tcComps_apgt[,p,g,t] < 0 ) )
            next

          thisSampSize <- sum( tcComps_apgt[,p,g,t] )
          thisWt <- sqrt(meanA_pg[p,g] / thisSampSize)

          wtErr_pgt[p,g,t] <- log(thisWt) - log(ageWt_pgt[p,g,t])
          # Count the number of bin
          posBins <- which(tcComps_apgt[,p,g,t] > 0)
          nBins   <- length(posBins)
          nBins_pgt[p,g,t] <- nBins - 1

          # Calculate the covariance matrix
          Kmat <- matrix( 0, nrow = nBins - 1, ncol = nBins )
          Fmat <- matrix( 0, nrow = nBins - 1, ncol = nBins )
          Hmat <- matrix( 1, nrow = nBins - 1, ncol = nBins - 1)
          Kmat[1:(nBins - 1),1:(nBins-1)] <- diag(1, nBins - 1 )
          Fmat[1:(nBins - 1),1:(nBins-1)] <- diag(1, nBins - 1 )
          Kmat[,nBins] <- -1
          Fmat[,nBins] <- -1
          diag(Hmat) <- 2

          Corr    <- Corr_gaa[g,posBins,posBins]
          Vcheck  <- Kmat %*% Corr %*% t(Kmat)
          Vinv    <- ginv(Vcheck)
          logdetV <- log(det(Vcheck))
          Hinv    <- ginv(Hmat)

          Gamma   <- t(Fmat) %*% Hinv %*% Vcheck %*% Hmat %*% Fmat


          # Get obs and pred vecs without zeroes
          obsVec  <- tcComps_apgt[posBins,p,g,t]
          obsVec  <- obsVec / sum(obsVec)
          predVec <- tcPred_apgt[posBins,p,g,t]

          gmObsErr_pgt[p,g,t]   <-  (prod(obsVec)^(1/nBins))  - gmObs_pgt[p,g,t]
          gmPredErr_pgt[p,g,t]  <-  (prod(predVec)^(1/nBins)) - gmPred_pgt[p,g,t]

          obsY  <- obsVec[-nBins] / obsVec[nBins]
          predY <- predVec[-nBins] / predVec[nBins]

          # Calculate SSR
          etaSumSq_pgt[p,g,t] <- (t( log(obsY) - log(predY) ) %*% Vinv %*% ( log(obsY) - log(predY) )) / thisWt^2

          # Calculate internal component 
          intNLLCpt_pgt[p,g,t] <- 0.5 * logdetV + (nBins - 1) * log(thisWt);
        } 

        totBins <- sum(nBins_pgt[p,g,])

        # Compute internal likelihood component
        internalLike        <- sum( intNLLCpt_pgt[p,g,] )
        diffIntCpt_pg[p,g]  <- intrnlAgeLikCpt_pg[p,g] - internalLike
        
        # Compute sigma
        tau_pg[p,g]         <- sqrt(sum( etaSumSq_pgt[p,g,] ) / totBins )
        difftau_pg[p,g]     <- tau_pg[p,g] - sqrt(tau2Age_pg[p,g])

        # total Likelihood
        totLike             <- log(tau_pg[p,g]) * totBins + 
                                internalLike + 
                                0.5 * sum(etaSumSq_pgt[p,g,]) / (tau_pg[p,g]^2)

        diffEtaSumSq_pg[p,g] <- sum( etaSumSq_pgt[p,g,] ) - modelEtaSumSq_pg[p,g]
        diffTotLike_pg[p,g]  <- totLike - ageObsNLL_pg[p,g]

      }

    }
  
  years <- fYear:(fYear + nT - 1)

  errorCols <- RColorBrewer::brewer.pal(3, "Set2" )

  par(  mfrow = c(3,nP),
        mar = c(1,1,1,1),
        oma = c(3,4,2,3) )

  for( g in 1:3 )
    for( p in 1:nP )
    {
      plot( x = range(years), 
            y = range(wtErr_pgt,gmObsErr_pgt,gmPredErr_pgt,na.rm = T),
            type = "n", xlab = "", ylab = "",
            axes = FALSE )
        mfg <- par("mfg")
        if( mfg[1] == mfg[3] )
          axis(side = 1 )
        if( mfg[2] == 1 )
          axis( side = 2, las = 1 )
        if( mfg[1] == 1 )
          mtext( side = 3, text = stockLabs[p] )

        if( mfg[2] == mfg[4] )
          mtext( side = 4, text = gearLabs[g] )
        box()
        # Plot points for each error seriess
        points( x = years, y = wtErr_pgt[p,g,],
                pch = 1, col = errorCols[1] )
        points( x = years, y = gmObsErr_pgt[p,g,],
                pch = 1, col = errorCols[2] )
        points( x = years, y = gmPredErr_pgt[p,g,],
                pch = 1, col = errorCols[3] )

        text( x = years[nT - 10], y = c(.5,.6,.7),
              labels = c( paste("diffNLL = ", round(diffTotLike_pg[p,g],2)),
                          paste("diffSSR = ", round(diffEtaSumSq_pg[p,g],2)),
                          paste("diffInt = ", round(diffIntCpt_pg[p,g],2)) ) )
    }

  mtext( side = 1, outer = TRUE, text = "Year" )
  mtext( side = 2, outer = TRUE, text = "Error" )

}

# Plot comp fits
plotMixedCompFitYrs <- function(  repList = reports,
                                  comps = "age",
                                  save = FALSE,
                                  fitLines = TRUE,
                                  savePath = "plotFitYrs" )
{
  repObj    <- repList$repOpt 
  initYear  <- repList$fYear
  gearLabels<- repList$gearLabs

  # Pull predicted and observed ages
  if( comps == "age" )
  {
    max       <- repObj$nA
    pred_xgt  <- repObj$predMixedPA_agt[1:max,,, drop = FALSE]
    obs_xgt   <- repObj$mA_agt[1:max,,, drop = FALSE]  
    xLab      <- "Age"
  }
  

  dimNames  <- dimnames(obs_xgt)
  compNames <- dimNames[[1]]
  yearNames <- dimNames[[3]]
  gearNames <- dimNames[[2]]

  # Pull model dims
  nG      <- repObj$nG
  nT      <- repObj$nT  
  minPA   <- repObj$minPropAge

  # Make colours vector
  cols    <- brewer.pal( n = nG, "Paired" )

  # Make years vector
  years   <- seq(initYear, length = nT+1, by = 1)

  # Now, we want to loop over gear types now, 
  # and record the gIdxes for which there are
  # observations
  obsGears <- c()
  gearTimes <- vector(mode = "list", length = nG)
  for( gIdx in 1:nG )
  {
    if( any(obs_xgt[1,gIdx,] >= 0) )
    {
      obsGears <- c(obsGears,gIdx)
      gearTimes[[gIdx]] <- which(obs_xgt[1,gIdx,] >= 0)
    }
  }

  # ok, obsGears are the ones we want to plot,
  # and gearTimes is the list of time indices
  for( gIdx in obsGears )
  {
    if(!save)
      dev.new()

    if(save)
    {
      gearPath <- paste(savePath,gearNames[gIdx],".png",sep = "")
      png(  gearPath, 
            width = 11, height = 8.5,
            units = "in", res = 300 )
    }

    times <- gearTimes[[gIdx]]
    # Count the number of age observations
    # there are, and make the plotting window
    nObs <- length(times)
    nCols <- round(sqrt(nObs))
    nRows <- ceiling(nObs/nCols)

    par(  mfcol = c(nRows,nCols), 
          mar = c(1,2,1,2),
          oma = c(3,3,3,3) )

    for( tIdx in times )
    { 
      Nobs <- sum(obs_xgt[,gIdx,tIdx])
      # get age obs and preds
      compObsProp_x  <- obs_xgt[,gIdx,tIdx]/sum(obs_xgt[,gIdx,tIdx])
      compPred_x     <- pred_xgt[,gIdx,tIdx]



      plot( x = c(1,max), y = c(0,max(compObsProp_x,compPred_x,na.rm = T) ),
            xlab = "", ylab = "", type = "n", las = 1 )
        rect( xleft = 1:max - .3, xright = 1:max + .3,
              ybottom = 0, ytop = compObsProp_x,
              col = "grey40", border = NA )
        if( fitLines)
        {
          lines(  x = 1:max, y = compPred_x, lwd = 1,
                  col = cols[gIdx] )
          points(  x = 1:max, y = compPred_x,
                  col = cols[gIdx], pch = 16 )
        }
        
        abline( h = minPA, lty = 2, lwd = .8 )
        panLab( x=.5, y = .95, txt = years[tIdx] )
        panLab( x=.1, y = .95, txt = paste( "N = ", Nobs, sep = "" ) )

    }
    mtext( side = 3, outer = T, text = gearNames[gIdx], line = 2, font = 2)
    mtext(  side = 1, outer = T, text = xLab, line = 2 )
    mtext(  side = 2, outer = T, text = paste("Proportion-at-",xLab,sep=""), 
            line = 2 )
    if(save)
      dev.off()
  }

} # END plotCompFitYrs()


# plotKobePhase()
# Phase plot of U/Umsy vs B/Bmsy (Kobe plot)
plotKobePhase <- function(  repList = reports,
                            pIdx = 1 )
{
  repObj <- repList$repOpt
  # Pull reference points object
  if(is.null(repObj$refPts))
    repObj <- calcRefPts(repObj)

  refPoints <- repObj$refPts
  refCurves <- refPoints$refCurves

  stockNames <- dimnames(repObj$SB_pt)[[1]]

  nT <- repObj$nT

  # pull equilibria
  F           <- refCurves$F
  Yeq_f       <- refCurves$lYeq_pf[pIdx,]
  Ueq_f       <- refCurves$lUeq_pf[pIdx,]
  Beq_f       <- refCurves$SBeq_pf[pIdx,]
  legBeq_f    <- refCurves$lBeq_pf[pIdx,]
  Fmsy        <- refPoints$FmsyRefPts$Fmsy_p[pIdx]
  Umsy        <- refPoints$FmsyRefPts$lUmsy_p[pIdx]
  YeqFmsy     <- refPoints$FmsyRefPts$lYeqFmsy_p[pIdx]
  BeqFmsy     <- refPoints$FmsyRefPts$SBeqFmsy_p[pIdx]
  lBeqFmsy    <- refPoints$FmsyRefPts$lBeqFmsy_p[pIdx]
  B0          <- repObj$B0_p[pIdx]

  # Pull time series
  SB_t  <- repObj$SB_pt[pIdx,]
  U_t   <- repObj$U_pt[pIdx,]

  relB_t    <- SB_t / BeqFmsy
  relU_t    <- U_t / Umsy
  
  relBeq_f  <- Beq_f/BeqFmsy
  relUeq_f  <- Ueq_f/Umsy

  maxFidx <- max(which(Yeq_f > 0))
  maxF <- F[maxFidx]
  maxU <- Ueq_f[maxFidx]


  plot(x = c(0,3), y = c(0,3), type = "n",
        las = 1, ylab = "U / Umsy", xlab = "B / Bmsy" )
    grid()
    # Make the eqbm line
    lines(y = relUeq_f[1:maxFidx], x = relBeq_f[1:maxFidx], 
          col = "grey80", lwd = 2 )
    # Show quadrant borders
    abline( v = 1, lwd = 1, lty = 2)
    abline( h = 1, lwd = 1, lty = 2)

    arrows( y0 = relU_t[1:(nT-1)], y1 = relU_t[2:nT],
            x0 = relB_t[1:(nT-1)], x1 = relB_t[2:nT],
            length = .1, angle = 20,
            col = "grey30" )

    text(y = relU_t[1]+.1,x = relB_t[1]-.1, labels = "1970")
    text(y = relU_t[nT]+.1,x = relB_t[nT]+.1, labels = "2020")


}


# RefPtsU
# Plot equilibrium yield and spawning biomass curves as a function
# of legal harvest rate. 
plotRefPtsU <- function(  repList = reports,
                          pIdx = 1 )
{
  repObj <- repList$repOpt
  # Pull reference points object
  if(is.null(repObj$refPts))
    repObj <- calcRefPts(repObj)

  refPoints <- repObj$refPts
  refCurves <- refPoints$refCurves

  stockNames <- dimnames(repObj$SB_pt)[[1]]

  isPosts <- !is.null(repList$posts)

  # pull equilibrium yields
  F           <- refCurves$F
  Yeq_f       <- refCurves$lYeq_pf[pIdx,]
  Ueq_f       <- refCurves$lUeq_pf[pIdx,]
  Beq_f       <- refCurves$SBeq_pf[pIdx,]
  lBeq_f      <- refCurves$lBeq_pf[pIdx,]
  Fmsy        <- refPoints$FmsyRefPts$Fmsy_p[pIdx]
  Umsy        <- refPoints$FmsyRefPts$lUmsy_p[pIdx]
  YeqFmsy     <- refPoints$FmsyRefPts$lYeqFmsy_p[pIdx]
  BeqFmsy     <- refPoints$FmsyRefPts$SBeqFmsy_p[pIdx]
  lBeqFmsy    <- refPoints$FmsyRefPts$lBeqFmsy_p[pIdx]
  B0          <- repObj$B0_p[pIdx]

  # Set dummies in case posts don't exist
  YeqFmsy_q <- NA


  if(isPosts)
  {
    B0        <- mean(repList$posts$B0_ip[,pIdx])
    B0_q      <- quantile(  repList$posts$B0_ip[,pIdx], 
                            probs = c(0.025,0.5,0.975) )

    # Yield
    lYeq_if <- repList$posts$lYeq_ipf[,pIdx,]
    lYeq_if[lYeq_if < 0] <- 0
    Yeq_qf    <- apply( X = lYeq_if, FUN = quantile,
                        MARGIN = 2, probs = c(0.025,0.5,0.975), na.rm = T )

    Yeq_f     <- apply( X = repList$posts$lYeq_ipf[,pIdx,], FUN = mean,
                        MARGIN = 2 )

    # Just one mean harvest rate to avoid plotting confusion
    lUeq_if <- repList$posts$lUeq_ipf[,pIdx,]
    lUeq_if[lYeq_if < 0] <- NA
    Ueq_f     <- apply( X = lUeq_if, FUN = mean,
                        MARGIN = 2,na.rm = T )

    # Biomass
    Beq_qf    <- apply( X = repList$posts$SBeq_ipf[,pIdx,], FUN = quantile,
                        MARGIN = 2, probs = c(0.025,0.5,0.975) )
    Beq_f     <- apply( X = repList$posts$SBeq_ipf[,pIdx,], FUN = mean,
                        MARGIN = 2 )

    # Legal Biomass
    lBeq_if   <- repList$posts$lBeq_ipf[,pIdx,]
    lBeq_if[lBeq_if < 0] <- 0
    lBeq_qf   <- apply( X = lBeq_if, FUN = quantile,
                        MARGIN = 2, probs = c(0.025,0.5,0.975) )
    lBeq_f    <- apply( X = repList$posts$lBeq_ipf[,pIdx,], FUN = mean,
                        MARGIN = 2 )

    # Umsy
    Umsy_q    <- quantile(  repList$posts$lUmsy_ip[,pIdx],
                            probs = c(0.025,0.5,0.975) )
    Umsy      <- mean( repList$posts$lUmsy_ip[,pIdx])

    # MSY
    YeqFmsy_q <- quantile(  repList$posts$lMSY_ip[,pIdx],
                            probs = c(0.025,0.5,0.975) )
    YeqFmsy   <- mean( repList$posts$lMSY_ip[,pIdx] )


    # Bmsy
    BeqFmsy_q <- quantile(  repList$posts$SBmsy_ip[,pIdx],
                            probs = c(0.025,0.5,0.975) )
    BeqFmsy   <- mean( repList$posts$SBmsy_ip[,pIdx])

    # Legal Bmsy
    lBeqFmsy_q <- quantile(  repList$posts$lBmsy_ip[,pIdx],
                            probs = c(0.025,0.5,0.975) )
    lBeqFmsy   <- mean( repList$posts$lBmsy_ip[,pIdx])

  }

  maxFidx <- max(which(Yeq_f > 0))
  maxF <- F[maxFidx]
  maxU <- Ueq_f[maxFidx]



  
  par(mfrow = c(2,1), mar = c(.1,.1,.1,.1), oma = c(4,4,1,1))

  plot( x = c(0,maxU), y = c(0,max(YeqFmsy,YeqFmsy_q,na.rm = T)), type = "n", 
        xlab = "", ylab = "",
        axes = FALSE )
    axis(side = 2, las = 1)
    grid()
    box()
    mtext(side = 2, text = "Legal yield (kt)", line = 3 )
    if(isPosts)
    {

      polygon(  x = c(Ueq_f[1:maxFidx],rev(Ueq_f[1:maxFidx])),
                y = c(Yeq_qf[1,1:maxFidx],Yeq_qf[3,maxFidx:1]),
                border = NA, col = "grey60" )

      abline( v = Umsy_q[c(1,3)], lty = 2, lwd = 1, col = "grey75" ) 
    }
    lines(x = Ueq_f[1:maxFidx], y = Yeq_f[1:maxFidx],
            col = "grey30", lwd = 2 )
    points( x = Umsy, y = YeqFmsy, pch = 16, col = "grey30",
            cex = 2 )
    # if( isPosts )
    #   segments( x0 = Umsy_q[1], x1 = Umsy_q[3], lty = 1, lwd = 1,
    #             y0 = YeqFmsy )

    legend( x = "topright", bty = "n",
            legend = c( paste0("legUmsy = ", round(Umsy,2)),
                        paste0("legMSY = ",  round(YeqFmsy,2)) ))


  plot( x = c(0,maxU), y = c(0,max(B0,lBeq_f,na.rm = T)), type = "n", 
        xlab = "Legal Harvest Rate", ylab = "Spawning Biomass (kt)",
        axes = FALSE  )
    axis(side = 1)
    axis(side = 2, las = 1)
    grid()
    box()
    mtext(side = 1, text = "Legal harvest rate", line = 2)
    mtext(side = 2, text = "Biomass (kt)", line = 3)
    if(isPosts)
    {
      legalPolyCol <- scales::alpha("steelblue", alpha = 0.5)

      spawnPolyCol <- scales::alpha("red", alpha = 0.5)

      # Plot legal biomass poly
      polygon(  x = c(Ueq_f[1:maxFidx],Ueq_f[maxFidx:1]),
                y = c(lBeq_qf[1,1:maxFidx],lBeq_qf[3,maxFidx:1]),
                border = NA, col = legalPolyCol )
      # Ploy spawning biomass polygon
      polygon(  x = c(Ueq_f[1:maxFidx],Ueq_f[maxFidx:1]),
                y = c(Beq_qf[1,1:maxFidx],Beq_qf[3,maxFidx:1]),
                border = NA, col = spawnPolyCol )      

      abline( v = Umsy_q[c(1,3)], lty = 2, lwd = 1, col = "grey75" ) 
    }
    lines(x = Ueq_f[1:maxFidx], y = Beq_f[1:maxFidx],
            col = "red", lwd = 2 )
    lines(x = Ueq_f[1:maxFidx], y = lBeq_f[1:maxFidx],
            col = "steelblue", lwd = 2 )
    points( x = Umsy, y = BeqFmsy, pch = 16, col = "red",
            cex = 2 )
    points( x = Umsy, y = lBeqFmsy, pch = 16, col = "steelblue",
            cex = 2 )
    legend( x = "topright", bty = "n",
            legend = c( "Spawning biomass",
                        paste0("SBmsy = ", round(BeqFmsy,2)),
                        "Legal biomass",
                        paste0("legBmsy = ",  round(lBeqFmsy,2)) ),
            lwd = c(2,NA,2,NA), 
            pch = c(NA,16,NA,16),
            col = c("red","red","steelblue","steelblue") )


} # END RefPtsU()

# plotAvgAgeFitComparison()
# Comparison of age fits, averaged over time, 
# for ISCAM and SISCA.
plotAvgAgeFitComparison <- function(  repObj      = reports$repOpt,
                                      iscamRep    = ISCAMrep,
                                      fYearSISCA  = reports$fYear,
                                      lYearSISCA  = reports$lYear )
{
  # Pull predicted and observed ages
  predAge_apgt <- repObj$predPA_apgt
  obsAge_apgt  <- repObj$A_apgt

  gearLabs  <- dimnames(obsAge_apgt)[[3]]

  # replace missing entries with NAs
  predAge_apgt[predAge_apgt < 0] <- NA
  obsAge_apgt[obsAge_apgt < 0] <- NA

  fYearISCAM <- iscamRep$syr
  lYearISCAM <- iscamRep$nyr


  # Get time steps
  nT <- dim(predAge_apgt)[4]
  nG <- dim(predAge_apgt)[3]
  nP <- dim(predAge_apgt)[2]


  # Multiply the yearly predicted proportions by the
  # yearly sample sizes
  for( g in 1:nG )
    for( p in 1:nP )
      for( t in 1:nT )
      {
        thisSamp <- sum( obsAge_apgt[,p,g,t], na.rm = T )
        if( thisSamp < 0 )
          thisSamp <- 0

        predAge_apgt[,p,g,t] <- thisSamp * predAge_apgt[,p,g,t]
      }

  
  
  # Average over time
  predAge_ag <- apply( X = predAge_apgt, FUN = sum, MARGIN = c(1,3), na.rm = T )
  obsAge_ag <- apply( X = obsAge_apgt, FUN = sum, MARGIN = c(1,3), na.rm = T )    

  
  # Pull model dims
  nG      <- repObj$nG
  nA      <- repObj$nA
  nT      <- repObj$nT
  nP      <- repObj$nP
  minPA   <- repObj$minPropAge

  # Now pull ISCAM ages
  iscamObsAge   <- iscamRep$d3_A
  iscamPredAge  <- iscamRep$A_hat

  # Gotta rearrange
  iscamObsAge_agt   <- array( NA, dim = c(nA,nG,nT) )
  iscamPredAge_agt  <- array( NA, dim = c(nA,nG,nT) )

  for( k in 1:nrow(iscamObsAge) )
  {
    iscamYr <- iscamObsAge[k,1]

    if( iscamYr < fYearSISCA | iscamYr > lYearSISCA )
      next

    gearIdx <- iscamObsAge[k,2]
    yrIdx   <- iscamYr - fYearSISCA + 1

    ageObs <- iscamObsAge[k,6:14]
    sampSize <- sum(ageObs)

    iscamObsAge_agt[2:10,gearIdx,yrIdx] <- ageObs
    iscamPredAge_agt[2:10,gearIdx,yrIdx] <- sampSize * iscamPredAge[k,]

  }

  iscamObsAge_ag  <- apply( X = iscamObsAge_agt, FUN = sum, MARGIN = c(1,2), na.rm = T)
  iscamPredAge_ag <- apply( X = iscamPredAge_agt, FUN = sum, MARGIN = c(1,2), na.rm = T)

  # Now, we want to loop over gear types now, 
  # and record the gIdxes for which there are
  # observations
  ageGears <- c()
  gearTimes <- vector(mode = "list", length = nG)
  for( gIdx in 1:nG )
    if( sum(obsAge_ag[,gIdx]) > 0  )
      ageGears <- c(ageGears,gIdx)

  # Make colours vector
  cols    <- RColorBrewer::brewer.pal( n = nG, "Dark2" )


  par(  mfrow = c(length(ageGears),2), 
        mar = c(1,2,1,2),
        oma = c(3,3,3,3) )

  # ok, ageGears are the ones we want to plot,
  for( aIdx in 1:length(ageGears) )
  { 
    gIdx <- ageGears[aIdx]
    # get age obs and preds for SISCA
    ageObs  <- obsAge_ag[,gIdx]
    nSamp   <- sum(ageObs)
    ageObs  <- ageObs / sum(ageObs)
    agePred <- predAge_ag[,gIdx]
    agePred <- agePred / sum( agePred )

    plot( x = c(1,nA), y = c(0,1.1*max(ageObs,agePred,na.rm = T) ),
          xlab = "", ylab = "", type = "n", las = 1 )
      legend( x = "topright",
              bty = "n",
              legend = paste("N = ", nSamp, sep = "" ) )
      rect( xleft = 3:nA - .3, xright = 3:nA + .3,
            ybottom = 0, ytop = ageObs[3:nA],
            col = "grey40", border = NA )
      lines(  x = 3:nA, y = agePred[3:nA], lwd = 3,
              col = cols[gIdx] )
      points(  x = 3:nA, y = agePred[3:nA],
              col = cols[gIdx], pch = 16, cex = 1.5 )
      abline( h = minPA, lty = 2, lwd = .8 )
      
      
      # panLab( x=.5, y = .95, txt = gearLabels[gIdx] )
      mfg <- par("mfg")
      if(mfg[1] == 1 )
        mtext( side = 3, text = "SISCAH", line = 1, font = 2)

    # get age obs and preds for ISCAM
    ageObs  <- iscamObsAge_ag[,gIdx]
    nSamp   <- sum(ageObs)
    ageObs  <- ageObs / sum(ageObs)
    agePred <- iscamPredAge_ag[,gIdx]
    agePred <- agePred / sum( agePred )

    plot( x = c(1,nA), y = c(0,1.1*max(ageObs,agePred,na.rm = T) ),
          xlab = "", ylab = "", type = "n", las = 1 )
      legend( x = "topright",
              bty = "n",
              legend = paste("N = ", nSamp, sep = "" ) )
      rect( xleft = 2:nA - .3, xright = 2:nA + .3,
            ybottom = 0, ytop = ageObs[2:nA],
            col = "grey40", border = NA )
      lines(  x = 2:nA, y = agePred[2:nA], lwd = 3,
              col = cols[gIdx] )
      points(  x = 2:nA, y = agePred[2:nA],
              col = cols[gIdx], pch = 16, cex = 1.5 )
      abline( h = 0.02, lty = 2, lwd = .8 )
      
      
      # panLab( x=.5, y = .95, txt = gearLabels[gIdx] )
      mfg <- par("mfg")
      if(mfg[1] == 1 )
        mtext( side = 3, text = "ISCAM", line = 1, font = 2)


      if( mfg[2] == mfg[4] )
      {
        corners <- par("usr") #Gets the four corners of plot area (x1, x2, y1, y2)
        par(xpd = TRUE) #Draw outside plot area
        text( x = corners[2]+0.5, 
              y = mean(corners[3:4]), 
              labels = gearLabs[gIdx], srt = 270,
              font = 2, cex = 1.5 )
        par(xpd = FALSE)
      }
  }
  
  mtext( side = 1, outer = T, text = "Age", line = 2 )
  mtext( side = 2, outer = T, text = "Proportion", line = 2 )

} # END plotAvgAgeFitComparison()


# plotRetroSBt()
# Takes in a report list, and 
# plots the spawning biomass retrospective
# plots
plotRetroSBt <- function( reportList = repList,
                          dep = FALSE,
                          relative = FALSE,
                          refIdx = 1,
                          legendText = NULL  )
{

  nFits     <- length(reportList)
  nTs       <- numeric( length = nFits )
  fYears    <- numeric( length = nFits )
  lYears    <- numeric( length = nFits )


  # Get number of stocks
  nP <- reportList[[1]]$repOpt$nP

  # Calculate Mohn's rho
  # Number of peels is nFits - 1
  Bpeel_pf  <- array( NA, dim = c(nP,nFits) )
  Bref_pf   <- array( NA, dim = c(nP,nFits) )

  stockNames <- reportList[[1]]$stock

  # Get nTs and model start/end years
  for( k in 1:nFits )
  {
    fYears[k] <-  reportList[[k]]$fYear
    lYears[k] <-  reportList[[k]]$lYear
  }
  nTs    <- lYears - fYears + 1
  maxT <- max(nTs)

  minfYear <- min(fYears)
  maxlYear <- max(lYears)

  if( all(fYears == minfYear ))
    retroType <- "end"
  if( all(lYears == maxlYear ))
    retroType <- "start"

  # Make an array for holding biomass
  bio_kpt <- array(NA, dim = c(nFits,nP,maxT) )

  for( k in 1:nFits )
  {
    thisfYear <- fYears[k]
    thislYear <- lYears[k]

    tdx <- thisfYear:thislYear - minfYear + 1

    bio_kpt[k,,tdx] <- reportList[[k]]$repOpt$SB_pt[,1:nTs[k]]

    # if(retroType == "end")
    # {
    #   Bpeel_pf[,k] <- reportList[[k]]$repOpt$SB_pt[,nTs[k]]
    #   Bref_pf[,k]  <- reportList[[refIdx]]$repOpt$SB_pt[,nTs[k]]
    # }

    # if(retroType == "start")
    # {
    #   Bpeel_pf[,k] <- reportList[[k]]$repOpt$SB_pt[,1]
    #   Bref_pf[,k]  <- reportList[[refIdx]]$repOpt$SB_pt[,tdx[1]]
    # }

    # Make depletion if asked for
    if( dep )
    {
      for( p in 1:nP)
        bio_kpt[k,p,tdx] <- bio_kpt[k,p,tdx] / reportList[[k]]$repOpt$B0_p[p]
    }
  }

  # mohnBerr_pf <- (Bpeel_pf - Bref_pf) / Bref_pf
  # mohnBerr_pf <- mohnBerr_pf[,-refIdx]

  # rho_p       <- round(apply( X = mohnBerr_pf, FUN = mean, MARGIN = 1, na.rm = T),3)

  # Replace zeroes with NAs
  bio_kpt[bio_kpt == 0 ] <- NA

  # Pull reference out
  refB_pt <- array(NA, dim = c(nP,maxT))
  refB_pt[1:nP,] <- bio_kpt[refIdx,1:nP,]

  if( relative )
  {
    for( k in 1:nFits )
      bio_kpt[k,,] <- bio_kpt[k,,] / refB_pt
  
    refB_pt[1:nP,] <- refB_pt[1:nP,] / refB_pt[1:nP,]
  }


  # And remove from biomasss
  bio_kpt <- bio_kpt[-refIdx,,,drop = FALSE]

  # Pull out colours

  caseCols <- scales::viridis_pal(alpha = 1, 
                                  begin = 0,
                                  end = 1,
                                  direction = -1)(nFits-1)

  par( mfrow = c( nP,1 ), mar = c(.1,1,.1,2),
        oma = c(3,4,2,2) )

  for( p in 1:nP )
  {
    plot( x = c(minfYear,maxlYear),
          y = c(0,max(1,bio_kpt[,p,],refB_pt[p,],na.rm = T)),
          type = "n", axes = FALSE )
      mfg <- par("mfg")
      if( mfg[1] == mfg[3] )
        axis( side = 1 )

      if( mfg[2] == mfg[4] )
      {
        rmtext( txt = stockNames[p], outer = TRUE,
                font = 2, line = .02, cex = 1.5)
      }

      axis( side = 2, las = 1 )
      box()
      grid()


      # rhoLab <- bquote(rho == .(rho_p[p]) )

      # text( x = 1960, y = 1 * max(bio_kpt[,p,],na.rm = T),
      #       labels = rhoLab )

      # Now do the reference model
      lines(  x = minfYear:maxlYear, y = refB_pt[p,],
              col = "black", lwd = 3 )

      # And retro lines
      for( k in 1:(nFits-1) )
      {
        lines( x = minfYear:maxlYear, y = bio_kpt[k,p,],
                col = caseCols[k], lwd = 2 )
      }

      
  }

  if(!is.null(legendText))
  {
    # legOrder <- order(legendText)
    # legCols <- character(length = nFits)
    # legCols[refIdx] <- "black"
    # legCols[-refIdx] <- caseCols
    legOrder <- order(legendText,decreasing =T)

    legend( x = "topleft", bty = "n",
            legend = legendText[legOrder],
            lty = 1, lwd = 2,
            col = c('black',caseCols[legOrder[-1]-1]))
  }

  yLab <- "Spawning Biomass (kt)"
  if( dep )
    yLab <- expression( paste("Spawning Biomass Depletion (", SB[t]/SB[0], ")", sep = "") )

  mtext( outer = TRUE, line = 2, "Year", side = 1 )
  mtext( outer = TRUE, line = 2, yLab,   side = 2 )


} # END plotRetroSBt()

# plotBatchMt()
# Takes in a report list, and 
# plots the spawning biomass retrospective
# plots
plotBatchMt <- function(  reportList = repList,
                          refIdx = 1,
                          relative = FALSE,
                          legendText = NULL,
                          aIdx = 2  )
{

  nFits     <- length(reportList)
  nTs       <- numeric( length = nFits )
  fYears    <- numeric( length = nFits )
  lYears    <- numeric( length = nFits )


  # Get number of stocks
  nP <- reportList[[1]]$repOpt$nP

  # Calculate Mohn's rho
  # Number of peels is nFits - 1
  Mpeel_pf  <- array( NA, dim = c(nP,nFits) )
  Mref_pf   <- array( NA, dim = c(nP,nFits) )

  stockNames <- reportList[[1]]$stock

  # Get nTs and model start/end years
  for( k in 1:nFits )
  {
    fYears[k] <-  reportList[[k]]$fYear
    lYears[k] <-  reportList[[k]]$lYear
  }
  nTs    <- lYears - fYears + 1
  maxT <- max(nTs)

  minfYear <- min(fYears)
  maxlYear <- max(lYears)

  if( all(fYears == minfYear ))
    retroType <- "end"
  if( all(lYears == maxlYear ))
    retroType <- "start"

  # Make an array for holding biomass
  M_kpt <- array(NA, dim = c(nFits,nP,maxT) )

  for( k in 1:nFits )
  {
    thisfYear <- fYears[k]
    thislYear <- lYears[k]

    tdx <- thisfYear:thislYear - minfYear + 1

    M_kpt[k,,tdx] <- reportList[[k]]$repOpt$M_apt[aIdx,,1:nTs[k]]

   
  }

  # Replace zeroes with NAs
  M_kpt[M_kpt == 0 ] <- NA

  # Pull reference out
  refM_pt <- array(NA, dim = c(nP,maxT))
  refM_pt[1:nP,] <- M_kpt[refIdx,1:nP,]

  if( relative )
  {
    for( k in 1:nFits )
      M_kpt[k,,] <- M_kpt[k,,] / refM_pt
  
    refM_pt[1:nP,] <- refM_pt[1:nP,] / refM_pt[1:nP,]
  }


  # And remove from biomasss
  M_kpt <- M_kpt[-refIdx,,,drop = FALSE]

  # Pull out colours

  caseCols <- scales::viridis_pal(alpha = 1, 
                                  begin = 0,
                                  end = 1,
                                  direction = -1)(nFits-1)

  par( mfrow = c( nP,1 ), mar = c(.1,1,.1,2),
        oma = c(3,4,2,2) )

  for( p in 1:nP )
  {
    plot( x = c(minfYear,maxlYear),
          y = c(0,max(1,M_kpt[,p,],refM_pt[p,],na.rm = T)),
          type = "n", axes = FALSE )
      mfg <- par("mfg")
      if( mfg[1] == mfg[3] )
        axis( side = 1 )

      if( mfg[2] == mfg[4] )
      {
        rmtext( txt = stockNames[p], outer = TRUE,
                font = 2, line = .02, cex = 1.5)
      }

      axis( side = 2, las = 1 )
      box()
      grid()

      # Now do the reference model
      lines(  x = minfYear:maxlYear, y = refM_pt[p,],
              col = "black", lwd = 3 )

      # And retro lines
      for( k in 1:(nFits-1) )
      {
        lines( x = minfYear:maxlYear, y = M_kpt[k,p,],
                col = caseCols[k], lwd = 2 )
      }

      
  }

  if(!is.null(legendText))
  {
    # legOrder <- order(legendText)
    # legCols <- character(length = nFits)
    # legCols[refIdx] <- "black"
    # legCols[-refIdx] <- caseCols

    idx       <- 1:nFits
    nonRefIdx <- idx[!idx%in%refIdx]
    legOrder <- c(refIdx,nonRefIdx)

    legend( x = "bottomleft", bty = "n",
            legend = legendText[legOrder],
            lty = 1, lwd = 2,
            col = c('black',caseCols))

  }

  yLab <- "Natural Mortaliyt (/yr)"
  if( relative )
    yLab <- "Natural mortality relative to reference case (unitless)"

  mtext( outer = TRUE, line = 2, "Year", side = 1 )
  mtext( outer = TRUE, line = 2, yLab,   side = 2 )


} # END plotBatchMt()


# plotRetroRt()
# Takes in a report list, and 
# plots the spawning biomass retrospective
# plots
plotRetroRt <- function(  reportList = repList,
                          dep = FALSE,
                          relative = FALSE,
                          refIdx = 1,
                          legendText = NULL,
                          cBrewPal='Dark2' )
{

  nFits     <- length(reportList)
  nTs       <- numeric( length = nFits )
  fYears    <- numeric( length = nFits )
  lYears    <- numeric( length = nFits )

  # Get number of stocks
  nP <- reportList[[1]]$repOpt$nP

  stockNames <- reportList[[1]]$stock


  # Get nTs and model start/end years
  for( k in 1:nFits )
  {
    fYears[k] <-  reportList[[k]]$fYear
    lYears[k] <-  reportList[[k]]$lYear
  }
  nTs    <- lYears - fYears + 1
  maxT <- max(nTs)

  minfYear <- min(fYears)
  maxlYear <- max(lYears)

  # Make an array for holding biomass
  R_kpt <- array(NA, dim = c(nFits,nP,maxT) )

  for( k in 1:nFits )
  {
    thisfYear <- fYears[k]
    thislYear <- lYears[k]

    tdx <- thisfYear:thislYear - minfYear + 1


    R_kpt[k,,tdx] <- reportList[[k]]$repOpt$R_pt[,1:nTs[k]]

    # Make depletion if asked for
    if( dep )
    {
      for( p in 1:nP)
        R_kpt[k,p,tdx] <- R_kpt[k,p,tdx] / reportList[[k]]$repOpt$R0_p[p]
    }    
  }

  R_kpt[R_kpt == 0 ] <- NA

  refR_pt <- array(NA, dim = c(nP,maxT))
  refR_pt[1:nP,] <- R_kpt[refIdx,1:nP,]

  if( relative )
  {
    

    for( k in 1:nFits )
      R_kpt[k,,] <- R_kpt[k,,] / refR_pt

  }
  refR_kpt  <- R_kpt[refIdx,,,drop = FALSE]
  R_kpt     <- R_kpt[-refIdx,,,drop = FALSE]

  # Pull out colours
  # caseCols <- RColorBrewer::brewer.pal(nFits,cBrewPal)

  caseCols <- scales::viridis_pal(alpha = 1, 
                                  begin = 0,
                                  end = 1,
                                  direction = -1)(nFits-1)


  par( mfrow = c( nP,1 ), mar = c(.1,1,.1,2),
        oma = c(3,4,2,2) )

  for( p in 1:nP )
  {
    plot( x = c(minfYear,maxlYear),
          y = c(0,max(1,R_kpt[,p,],refR_pt[p,],na.rm = TRUE)),
          type = "n", axes = FALSE )
      mfg <- par("mfg")
      if( mfg[1] == mfg[3] )
        axis( side = 1 )
      if( mfg[2] == mfg[4] )
      {
        corners <- par("usr") #Gets the four corners of plot area (x1, x2, y1, y2)
        par(xpd = TRUE) #Draw outside plot area
        text( x = corners[2]+1, 
              y = mean(corners[3:4]), 
              labels = stockNames[p], srt = 270,
              font = 2, cex = 1.5 )
        par(xpd = FALSE)
      }

      axis( side = 2, las = 1 )
      box()
      grid()


      # Now add lines for other models
      for( k in 1:(nFits-1) )   
          lines( x = minfYear:maxlYear, y = R_kpt[k,p,],
                col = caseCols[k], lwd = 1.5 )

      # Plot reference case first
      refColr <- adjustcolor('black',alpha.f=0.7)
      lines(  x = minfYear:maxlYear, y = refR_kpt[1,p,], 
              lwd = 3, col=refColr)


  }

  if(!is.null(legendText))
  {
    legOrder <- order(legendText,decreasing =T)

    legend( x = "topleft", bty = "n",
            legend = legendText[legOrder],
            lty = 1, lwd = 2,
            col = c('black',caseCols[legOrder[-1]-1]))
  }

  yLab <- "Recruitment (1e6)"
  if( dep )
    yLab <- expression( paste("Relative Recruitment  (", R[t]/R[0], ")", sep = "") )
  if( relative )
    yLab <- expression(paste( "Recruitment relative to reference model, (", R[t]/Rref[t] ,")" ))

  mtext( outer = TRUE, line = 2, "Year", side = 1 )
  mtext( outer = TRUE, line = 2, yLab,   side = 2 )


} # END plotRetroRt

# plotACF()
# Plots autocorrelation estimates for
# an array of time series.
plotACF <- function(arr_pt, max.lag = 10, maxy = .5)
{
  nP <- dim(arr_pt)[1]

  # Now plot indicators of +ve or 0 idx
  stockCols <- scales::viridis_pal(option = "C",alpha = 1,
                                    begin = 0, end = .5)(3)

  calcACF <- function(idx, arr_pt)
  {
    xVec <- arr_pt[idx,]

    ACF <- acf(xVec, na.action = na.pass, plot = FALSE)
    ACF
  }

  acfList <- lapply( X = 1:nP, FUN = calcACF, arr_pt = arr_pt)
  xJitter <- seq( from = -.15, to = .15, length.out = nP)

  plot( x = c(0, max.lag), y = c(-maxy, maxy),
        xlab = "", ylab = "", axes = FALSE, type = "n" )
  mfg <- par("mfg")
  if( mfg[1] == mfg[3])
    axis( side = 1, at = 0:max.lag)

  if( mfg[2] == 1)
    axis( side = 2, las = 1 )

  if( mfg[1] == 1 )
    legend( x = "topright",
            legend = c("C/S","JP/S","Lou"),
            lty = 1, lwd = 2,
            col = stockCols, bty = "n")

  box()
  grid()
  
  abline(h = 0, lty = 2, col = "grey30")

  for( p in 1:nP )
  {
    segments( x0 = acfList[[p]]$lag + xJitter[p],
              y0 = 0, y1 = acfList[[p]]$acf,
              col = stockCols[p], lwd = 2 )
  }
  # abline(v = .5 + c(0:10), lty = 1, lwd = 1)
} # END plotACF()

plotDeltaLNFits <- function( reports = reports )
{
  # Get reports and fYear
  repObj <- reports$repOpt
  fYear  <- reports$fYear
  datObj <- reports$data

  # Get dimension
  nP <- repObj$nP
  nG <- repObj$nG
  nT <- repObj$nT

  stock <- reports$stock



  # Get indicator vector for delta model
  deltaIdx_pg <- datObj$deltaIdx_pg
  deltaIdx_g  <- apply( X = deltaIdx_pg, FUN = sum, MARGIN = 2)
  deltaIndVar <- datObj$deltaIndVar

  # Pull biomass, probability of positive indices,
  B0_p              <- repObj$B0_p
  SB_pt             <- repObj$SB_pt
  probPosIdx_pgt    <- repObj$probPosIdx_pgt
  I_pgt             <- repObj$I_pgt

  probPosCombIdx_pt <- repObj$probPosCombIdx_pt


  maxBt <- max(SB_pt[,38:69], na.rm = T)
  BtSeq <- seq(from = 0, to = maxBt, length.out = 100)



  # Compute all probability lines for the posterior
  indI_pgt <- I_pgt
  indI_pgt[I_pgt > 0] <- 1
  indI_pgt[I_pgt < 0] <- NA


  if( any(datObj$combI_pt >= 0) )
    indI_pgt[indI_pgt == 0] <- NA

  indCombI_pt <- datObj$combI_pt
  indCombI_pt[indCombI_pt > 0] <- 1
  indCombI_pt[indCombI_pt < 0] <- NA

  # get model parameters too
  meanProbPosIdx_pg <- repObj$meanProbPosIdx_pg
  SDProbPosIdx_pg   <- repObj$SDProbPosIdx_pg

  if(!is.null(reports$posts))
  {
    probPosIdx_ipgt <- reports$posts$probPosIdx_ipgt
    nPosts <- dim(probPosIdx_ipgt)[1]
    probPosIdx_pgt  <- apply( X = probPosIdx_ipgt, FUN = mean, MARGIN = c(2,3,4))

    meanProbPosIdx_ipg  <- reports$posts$meanProbPosIdx_ipg
    SDProbPosIdx_ipg    <- reports$posts$SDProbPosIdx_ipg

    

    meanProbPosIdx_pg <- apply( X = reports$posts$meanProbPosIdx_ipg, FUN = mean, MARGIN = c(2,3))
    SDProbPosIdx_pg   <- apply( X = reports$posts$SDProbPosIdx_ipg  , FUN = mean, MARGIN = c(2,3))

    modelProbPos_igb <- array( 0, dim = c(nPosts,nG,100) )
    for(i in 1:nPosts )
    {
      for( g in 1:nG )
        modelProbPos_igb[i,g,] <- 1 / (1 + exp( - (meanProbPosIdx_ipg[i,1,g] + SDProbPosIdx_ipg[i,1,g] * BtSeq)))
    }

    modelProbPos_qgb <- apply( X = modelProbPos_igb, FUN = quantile, na.rm = T,
                                MARGIN = c(2,3), probs = c(0.025, 0.5, 0.975) )

    meanModelProb_gb <- apply( X = modelProbPos_igb, FUN = mean,
                                MARGIN = c(2,3), na.rm = T )
  } 

  probPosIdx_pgt[I_pgt < 0] <- NA




  par(  mar = c(.1,.1,.1,.1),
        oma = c(3,5.5,1,3) )


  fleetCols <- RColorBrewer::brewer.pal( nG, "Dark2" )
  fleetPch  <- 1:nG

  # Now plot indicators of +ve or 0 idx
  stockCols <- scales::viridis_pal(option = "C",alpha = 1,
                                    begin = 0, end = .5)(nP)

  plot( x = c(0,30),y = c(0,1),
        type = "n", axes = FALSE, xaxs = "i" ) 
    mfg <- par("mfg")
    if(mfg[1] == mfg[3])
    {
      axis( side = 1 )
      mtext(  side = 1, text = "Biomass",
              line = 2)
    }

    axis( side = 2, las = 1)

    box()
    grid()
    deltaFleets <- which(deltaIdx_g > 0)
    nDeltaFleets <- length(deltaFleets)

    # Loop over gears, 
    for( g in deltaFleets)
    {
      # plot model of detection prob
      if(is.null(reports$posts))
      {
        indVar <- BtSeq
        if(deltaIndVar == 2)
          indVar <- BtSeq / B0_p[p]
        modelProb <- 1 / (1 + exp( - (meanProbPosIdx_pg[1,g] + SDProbPosIdx_pg[1,g] * indVar)))
      } else {
        modelProb <- meanModelProb_gb[g,]
        polygon(  x = c(BtSeq,rev(BtSeq)),
                  y = c( modelProbPos_qgb[1,g,], rev(modelProbPos_qgb[3,g,]) ),
                  col = "grey70",
                  border = NA )
      }

      lines(  x = BtSeq, y = modelProb,
              col = fleetCols[g], lwd = 2 )

      for( p in 1:nP )
      {
        # Get ticks if using split indices
        if( any(!is.na(indI_pgt)))
        {
          zeroTicks <- SB_pt[p,which(indI_pgt[p,g,] == 0)] 
          oneTicks <- SB_pt[p,which(indI_pgt[p,g,] == 1)] 
        }

        if( any(!is.na(indCombI_pt)))
        {
          zeroTicks <- SB_pt[p,which(indCombI_pt[p,] == 0)] 
          oneTicks <- SB_pt[p,which(indCombI_pt[p,] == 1)]  
        }

        axis( side = 1, at = zeroTicks, labels = FALSE,
              tck = .02, col.ticks = stockCols[p])
        axis( side = 3, at = oneTicks, labels = FALSE,
              tck = .02, col.ticks = stockCols[p])
       
      }

      
    }

      
    legend( x = "bottomright",bty = "n",
            legend = c(stock,"Dive"),
            col = c(stockCols,fleetCols[deltaFleets]),
            pch = c(NA,NA,NA,22),
            pt.bg = "grey70",
            pt.lwd = 0,
            lty = 1,
            lwd = 1,
            pt.cex = 1.5 )

  mtext( side = 2, text = "Probability of\ndetecting spawn",
         outer = TRUE, line = 3, font = 2 )

}

plotIratio <- function( repList = reports,
                        gIdx = c(4,5) )
{
  datObj <- repList$data

  I_pgt     <- datObj$I_pgt
  rI_pgt    <- datObj$rI_pgt
  combI_pt  <- datObj$combI_pt

  stockLabs <- dimnames(I_pgt)[[1]]
  gearLabs  <- dimnames(I_pgt)[[2]]
  nT        <- dim(I_pgt)[3]
  nP        <-  dim(I_pgt)[1]

  yrs <- 1951:2019

  combI_pt[combI_pt <= 0] <- NA

  fleetCols <- brewer.pal(n = 6, "Dark2")

  par( mfrow = c(nP,1), oma = c(4,4,2,3), mar = c(.1,1,.1,1) )

  for( p in 1:nP)
  {
    maxIdx <- max(combI_pt[p,],na.rm = T)

    ratTicks <- seq(from = 0, to = maxIdx, length.out = 5)
    ratLabs  <- (0:4)/4
    plotRat <- rI_pgt[p,4,] * maxIdx
    plot( x = range(yrs), y = c(0,maxIdx),
          axes = FALSE, type = "n" )
      mfg <- par("mfg")
      if( mfg[1] == mfg[3] )
        axis( side = 1 )
      axis( side = 2, las = 1 )
      box()
      grid()
      # Surface survey
      rect( xleft = yrs - .3,
            xright = yrs + .3,
            ybottom = 0,
            ytop = combI_pt[p,] * rI_pgt[p,4,],
            col = fleetCols[4], border = NA )
      # Dive survey
      rect( xleft = yrs - .3,
            xright = yrs + .3,
            ybottom = combI_pt[p,] * rI_pgt[p,4,],
            ytop = combI_pt[p,],
            col = fleetCols[5], border = NA )
      if( p == 1 )
        legend( x = "topright", bty = "n",
                legend = c("Surface","Dive"),
                col = fleetCols[4:5],
                pch = 22,
                pt.bg = fleetCols[4:5],
                cex = 2)

    rmtext( txt = stockLabs[p], font = 2, line = 0.02, cex = 1.5, outer = TRUE )
  }
  mtext( side = 1, outer = TRUE, text =  "Year", line = 2)
  mtext( side = 2, outer = TRUE, text =  "Spawn Index (kt)", line = 2)
}

# plotJuveMLikeProfile()
plotJuveMLikeProfile <- function( groupFolder = "juveM_likeProf",
                                  prefix = "juveM" )
{ 
  # Load batch reports
  batchReports <- .loadBatch( groupFolder = groupFolder,
                            baseDir = "Outputs/fits",
                            prefix = prefix )
  # Load the base OM
  .loadFit("fit_MSbase")

  baseJuveM <- reports$repOpt$M_apt[1,2,1]


  baseModelLike <-  sum(reports$repOpt$obsIdxNLL_pg) +
                    sum(reports$repOpt$obsIdxDeltaNLL_pg) +
                    sum(reports$repOpt$obsCombIdxNLL_p) +
                    sum(reports$repOpt$obsCombIdxDeltaNLL_p) +
                    sum(reports$repOpt$ageObsNLL_pg) +
                    sum(reports$repOpt$obsMixedIdxNLL_g) +
                    sum(reports$repOpt$ageObsNLL_g)



  reports$repOpt$objFun
  baseB0    <- sum(reports$repOpt$B0_p)

  nModels <- length(batchReports)
  # Now get x value and likelihood function values
  juveM_m     <- numeric(length = nModels)
  modelLike_m <- numeric(length = nModels)
  B0_m        <- numeric(length = nModels)

  for( k in 1:nModels)
  {
    rep_k <- batchReports[[k]]$repOpt
    juveM_m[k]      <-  rep_k$M_apt[1,2,1]
    modelLike_m[k]  <-  sum(rep_k$obsIdxNLL_pg) +
                        sum(rep_k$obsIdxDeltaNLL_pg) +
                        sum(rep_k$obsCombIdxNLL_p) +
                        sum(rep_k$obsCombIdxDeltaNLL_p) +
                        sum(rep_k$ageObsNLL_pg) +
                        sum(rep_k$obsMixedIdxNLL_g) +
                        sum(rep_k$ageObsNLL_g)
    B0_m[k]         <-  sum(rep_k$B0_p)
  }

  modelLike_m <- modelLike_m[order(juveM_m)]
  B0_m <- B0_m[order(juveM_m)]
  juveM_m <- juveM_m[order(juveM_m)]


  minLikeIdx <- which.min(modelLike_m)

  likeDiff <- modelLike_m - baseModelLike

  # Fit a spline to likeDiff
  likeDiffSpline <- splinefun(  x = juveM_m, y = modelLike_m - baseModelLike )

  bounds <- c(juveM_m[minLikeIdx],2)
  zeroLikeDiff  <- uniroot(f = likeDiffSpline, interval = bounds)$root
  minLikeM      <- uniroot( f = likeDiffSpline, interval = c(0.6,2), deriv = 1)$root
  
  # Now make a spline for plotting
  likeSpline <- splinefun(  x = juveM_m, y = modelLike_m )
  plotM <- seq(from = 0.6, to = 2, by = 0.01)
  plotLike <- likeSpline(plotM)

  par(mfrow = c(2,1), mar = c(.1,.1,.1,.1), oma = c(3,5,2,2))
  plot( x = c(0.6,2),
        y = range(modelLike_m),
        las = 1, type = "n", xlab = "Age-1 M (/yr)",
        ylab = "Model minimum negative log-posterior",
        axes = FALSE)
    lines( x = plotM, y = plotLike, lty = 1, lwd = 3)
    points( x = baseJuveM, y  = baseModelLike, pch = 24, bg = "green", col = "black", cex = 2)
    points( x = minLikeM, y = likeSpline(minLikeM), 
            pch = 21, bg = "red", cex = 2 )
    abline( h = baseModelLike, lty = 2, lwd = 2, col = "steelblue")
    axis( side = 2, las = 1)
    box()
    points( x = zeroLikeDiff, y = likeSpline(zeroLikeDiff), 
            pch = 21, bg = "steelblue", lwd = 1, cex = 2 )
    mtext( side = 2, text = "Model minimum negative log-likelihood", line = 3.5)

    legend( x = "topleft", bty = "n",
            legend = c( "Base OM, Age-1 M = 0.72",
                        paste("Min NLL, Age-1M =",round(minLikeM,2)),
                        paste("Age-1 M =", round(zeroLikeDiff,2)),
                        "Base OM NLL"),
            pch = c(24,21,21, NA),
            lwd = c(NA,NA,NA,2),
            lty = c(NA,NA,NA,2),
            col = c("black","black","black","steelblue"),
            pt.bg = c("green","red","steelblue",NA),
            pt.lwd = 1,
            cex = 1,
            pt.cex = 2 )

  B0spline <- splinefun(  x = juveM_m, y = B0_m )

  plot( x = c(0.6,2),
        y = range(B0_m),
        las = 1, type = "n", xlab = "Age-1 M (/yr)",
        ylab = "Aggregate unfished biomass (kt)", axes = FALSE)
    axis( side = 1 )
    axis( side = 2, las = 1)
    box()
    lines( x = juveM_m, y = B0_m, lty = 1, lwd = 3)
    points( x = baseJuveM, y  = baseB0, pch = 24, bg = "green", cex = 2)
    points( x = minLikeM, y = B0spline(minLikeM), pch = 21, bg = "red", cex = 2 )
    points( x = zeroLikeDiff, y = B0spline(zeroLikeDiff), pch = 21, 
            bg = "steelblue", lwd = 1, cex = 2 )
    mtext( side = 2, text = "Aggregate unfished biomass (kt)", line = 3.5)
    mtext( side = 1, text = "Age-1 Natural Mortality (/yr)", line = 2)



} # END plotJuveMLikeProfile

# plotSelComparison()
# Comparison of SISCAN and previous model selectivity
# estimates
plotSelComparison <- function( repObj = reports$repOpt,
                                prevRep = SCALrep  )
{
  nP <- repObj$nP
  nT <- repObj$nT
  nG <- repObj$nG 
  nA <- repObj$nA 
  nX <- repObj$nX 

  prevG <- 6

  sel_axpgt <- repObj$sel_axpgt[,,,,1:nT,drop = FALSE]

  prevSel_axpgt <- array(NA, dim = dim(sel_axpgt))

  gearCols <- RColorBrewer::brewer.pal(nG,"Paired")

  for(g in 1:nG )
  {
    if( g > prevG )
      next
    for( x in 1:nX )
    {
      if(x == 1)
        selLab <- paste0("sel_gta_m",g)
      if(x == 2)
        selLab <- paste0("sel_gta_f",g)


      if(!is.null(prevRep[[selLab]]))
      {
        prevT <- dim(prevRep[[selLab]])[1]
        if(prevT > nT)
          prevT <- nT
        prevSel_axpgt[,x,1,g,1:prevT] <- t(prevRep[[selLab]][1:prevT,])
      }
    }
  }

  par(mfcol = c(nG,nX), mar = c(.1,1,.1,1),
      oma = c(4,4,2,2) )
    for( x in 1:nX )
      for( g in 1:nG )
      {
        plot(x = c(1,nA), y = c(0,1), type = "n",
              axes = FALSE )
          mfg <- par("mfg")
          if( mfg[1] == mfg[3])
            axis(side = 1)
          if( mfg[2] == 1)
            axis(side = 2, las = 1)
          grid()
          box()

          for( t in 1:nT )
          {
            lines(x = 1:nA, y = sel_axpgt[,x,1,g,t],
                    lwd = 0.8, col = "grey80" )
            lines(x = 1:nA, y = prevSel_axpgt[,x,1,g,t],
                    lwd = 0.8, col = "black", lty = 2 )
          }

          lines(  x = 1:nA, y = sel_axpgt[,x,1,g,1],
                  col = gearCols[g], lwd = 3 )
          lines(  x = 1:nA, y = prevSel_axpgt[,x,1,g,1],
                  col = "black", lty = 2, lwd = 3 )
      }
}

# PlotAggBtRtComparison()
# Comparison of SISCA and ISCAM Bt and age-2 Rt at the major stock
# level.
plotBtRtComparison <- function( repObj      = reports$repOpt,
                                datObj      = reports$data,
                                prevRep     = SCALrep,
                                fYearSISCA  = reports$fYear,
                                lYearSISCA  = reports$lYear )
{

  nP <- repObj$nP
  nT <- repObj$nT
  nG <- repObj$nG

  if( any( datObj$combI_pt > 0 ) )
  {
    combIdx <- TRUE
    qComb_pt <- repObj$qComb_pt
  } else combIdx <- FALSE


  # fYear <- repList$fYear
  yrs <- seq(from = fYearSISCA, by = 1, length.out = nT)

  # Get SBt, sum
  SB_pt <- repObj$SB_pt[1:nP,,drop = FALSE]
  SB_t  <- apply( X = SB_pt, FUN = sum, MARGIN = 2, na.rm = T)

  prevB0   <- prevRep$SSB0/1e3
  siscaB0  <- sum(repObj$B0_p)

  prevR0   <- prevRep$R0/1e6
  siscaR0   <- sum(repObj$R0_p)

  # Get age-2 numbers at age
  R_pt        <- repObj$R_pt

  R_pt[R_pt == 0] <- NA
  R_t   <- apply(X = R_pt, FUN = sum, MARGIN = 2, na.rm = T)


  # Get iscam SBt and age-2 recruitment
  prevSB_t  <- prevRep$SSBt/1e3
  Nt1_m     <- prevRep$Nta_m[,1]
  Nt1_f     <- prevRep$Nta_m[,2]
  prevR_t   <- Nt1_m/1e6 + Nt1_f/1e6    

  nTprev <- min(length(prevSB_t),length(prevR_t))


  # Now plot
  par(  mfrow = c(2,1), mar = c(0.1,1,.1,1),
        oma = c(3,4,1,1) )
  # First, SB
  plot( x = c(fYearSISCA,lYearSISCA),
        y = c(0, max(SB_t, prevSB_t, prevB0, siscaB0)), type = "n",
        axes = FALSE, xlab = "", ylab = "" )
    axis( side = 2, las = 1)
    box()
    grid()
    mtext( side = 2, text=  "Spawning Biomass (kt)", line = 3)

    # Plot SISCAH biomass
    lines(  x = yrs, y = SB_t[1:nT], col = "steelblue", lwd = 3)
  
    lines( x = yrs[1:nTprev], y = prevSB_t[1:nTprev],
            lwd = 3, col = "red" )

    legText <- c( "SCAL Bt",
                  "SCAL B0",
                  "SISCAH Bt",
                  "SISCAH B0" )

    legend( x = "topright", bty = "n",
            legend = legText,
            col = c("red","red","steelblue","steelblue" ),
            pt.bg = c(NA,NA,NA,NA,NA,NA,"grey40"),
            lwd = c(3,3,3,3,NA,NA,NA),
            pch = c(NA,NA,NA,NA,4,5,21),
            lty = c(1,2,1,3,NA,NA,NA),
            pt.cex = 1.5)

    abline( h = c(prevB0,siscaB0),
          lty = c(2,3), lwd = 2, col = c("red","steelblue") )

  plot( x = c(fYearSISCA,lYearSISCA),
        y = c(0, max(R_t, prevR_t)), type = "n",
        axes = FALSE, xlab = "", ylab = "" )
    axis( side = 2, las = 1)
    axis( side = 1 )
    mtext( side = 2, text=  "Age-1 Recruitment (1e6)", line = 3)
    mtext( side = 1, text=  "Year", line = 2)
    box()
    grid()
    legend( x = "topright", bty = "n",
            legend = c( "SISCAL Rt", "SISCAL R0", 
                        "SCAL Rt", "SCAL R0"),
            col = c("steelblue","steelblue","red","red"),
            lwd = c(3,3,3,3),
            lty = c(1,3,1,2))
    lines( x = yrs, y = R_t[1:nT],
            lwd = 3, col = "steelblue" )
    lines( x = yrs[1:nTprev], y = prevR_t[1:nTprev],
            lwd = 3, col = "red" )
    abline( h = c(prevR0,siscaR0),
            lty = c(2,3), lwd = 2, col = c("red","steelblue") )

} # END plotAggBtRtComparison()

# plotPerCapSP
# Per Capita surplus production plotted against
# biomass depletion
plotPerCapSP <- function( repList = reports,
                          minAge = NULL,
                          maxY = 2,
                          minYear = 1951  )
{
  repObj <- repList$repOpt

  # Pull biomass and catch
  B_apt   <- repObj$B_apt
  C_pgt   <- repObj$C_pgt
  nA      <- repObj$nA
  nP      <- repObj$nP
  nT      <- repObj$nT
  initT_p <- repObj$tInitModel_p
  B0_p    <- repObj$B0_p
  R0_p    <- repObj$R0_p

  fYear <- repList$fYear
  tInit <- minYear - fYear + 1


  if( !is.null(minAge))
  {
    B0minAge_p <- B0_p

    for( p in 1:nP )
    {
      eqbmN_a <- R0_p[p] * repObj$initSurv_ap[,p]
      B0minAge_p[p] <- sum(eqbmN_a[minAge:nA] * repObj$meanWt_ap[minAge:nA,p])
    }

    B0_p <- B0minAge_p
  }

  # Get biomass and catch
  B_pt    <- repObj$SB_pt[,1:nT]
  C_pt    <- apply( X = C_pgt, FUN = sum, MARGIN = c(1,3) )

  if( !is.null(minAge))
    B_pt    <- apply( X = B_apt[minAge:nA,,1:nT], FUN = sum, MARGIN = c(2,3) )
  
  # Calculate surplus production
  P_pt     <- array(NA, dim = c(nP,nT))

  for( p in 1:nP)
    for( t in (initT_p[p]+1):(nT-1) )
      P_pt[p,t] <- B_pt[p,t+1] - B_pt[p,t] + C_pt[p,t]



  # Now do the aggregate
  # B_t <- apply( X = B_apt[minAge:nA,,1:nT], FUN = sum, MARGIN = c(3))
  B_t <- apply( X = B_pt[,1:nT], FUN = sum, MARGIN = c(2))
  C_t <- apply( X = C_pgt, FUN = sum, MARGIN = c(3))

  P_t   <- array(NA, dim = c(nT))
  P_t[1:(nT-1)] <- B_t[2:nT] - B_t[1:(nT - 1)] + C_t[1:(nT-1)]



  ptCols <- scales::viridis_pal(option = "C",alpha = 1)(nT - tInit+1)

  names(ptCols) <- minYear:2019

  stockNames <- repList$stock

  pcP_pt <- P_pt/B_pt
  pcP_t  <- P_t/B_t

  depB_pt <- B_pt
  for( p in 1:nP )
    depB_pt[p,] <- B_pt[p,]/B0_p[p]
  depB_t  <- B_t/sum(B0_p)

  if( is.null(maxY) )
    maxY <- max(pcP_pt,pcP_t,na.rm =T)

  minY <- min(pcP_pt,pcP_t,na.rm =T)

  par(  mfrow = c(2,2),
        mar = c(1,.1,1,.1),
        oma = c(4,4,1,1) )
  for( p in 1:nP )
  {
    plot( x = c(0,max(B_pt[p,tInit:nT],na.rm = T)),
          y = c(minY,maxY),
          axes = FALSE, type = "n", las = 1 )
      mfg <- par("mfg")
      # if( mfg[1] == mfg[3])
        axis(side = 1)
      if( mfg[2] == 1)
        axis(side = 2, las = 1)
      box()
      abline(h = 0, lty = 2 )
      
      arrows( x0 = B_pt[p,tInit:(nT-1)],
              x1 = B_pt[p,(tInit+1):(nT)],
              y0 = P_pt[p,tInit:(nT-1)]/B_pt[p,tInit:(nT-1)],
              y1 = P_pt[p,(tInit+1):(nT)]/B_pt[p,(tInit+1):(nT)],
              col = "grey75", lwd = .8, length = .1 )
      points( x = B_pt[p,tInit:nT],
              y = P_pt[p,tInit:nT]/B_pt[p,tInit:nT],
              bg = ptCols, pch = 21 )
      abline( v = .3*B0_p[p], lty = 3, lwd = .8)
      abline( v = B0_p[p], lty = 3, lwd = .8, col = "red")

      legend( x = "topleft",
              legend = stockNames[p], bty = "n")

  }

  plot( x = c(0,max(B_t[tInit:nT],na.rm = T)),
          y = c(minY,maxY),
          axes = FALSE, type = "n", las = 1 )
    mfg <- par("mfg")
    # if( mfg[1] == mfg[3])
      axis(side = 1)
    if( mfg[2] == 1)
      axis(side = 2, las = 1)
    box()
    abline(h = 0, lty = 2 )
    arrows( x0 = B_t[tInit:(nT-1)],
            x1 = B_t[(tInit+1):(nT)],
            y0 = P_t[tInit:(nT-1)]/B_t[tInit:(nT-1)],
            y1 = P_t[(tInit+1):(nT)]/B_t[(tInit+1):(nT)],
            col = "grey75", lwd = .8, length = .1 )
    points( x = B_t[tInit:nT],
            y = P_t[tInit:nT]/B_t[tInit:nT],
            bg = ptCols, pch = 21 )

    abline( v = .3 * sum(B0_p), lty = 3, lwd = .8)
    abline( v = sum(B0_p), lty = 3, lwd = .8, col = "red")

    legend( x = "topright",
            bty = "n",
            legend = c( minYear,
                        "2019"),
            pch = 21, pt.bg = ptCols[c(1,nT-tInit+1)] )
    legend( x = "topleft",
            legend = "Aggregate",
            bty = "n")

    xLab <- "Spawning Biomass (kt)"

    if( !is.null(minAge))
      xLab <- paste( minAge,"+ Biomass (kt)",sep = "")

    mtext(side = 1, outer = TRUE, text = xLab,
          line = 2.5)
    mtext(side = 2, outer = TRUE, text = "Per Capita Surplus Production",
          line = 2.5)
} # END plotPerCapSP()

# plotRecPerSpawner()
# Per capita recruits per spawner, as a function
# of biomass
plotRecPerSpawner <- function(  repList = reports,
                                maxY = NULL  )
{
  repObj <- repList$repOpt

  # Pull biomass and catch
  SB_pt   <- repObj$SB_pt
  R_pt    <- repObj$R_pt
  C_pgt   <- repObj$C_pgt
  nA      <- repObj$nA
  nP      <- repObj$nP
  nT      <- repObj$nT
  initT_p <- repObj$tInitModel_p
  B0_p    <- repObj$B0_p

  B_pt    <- repObj$SB_pt[,1:nT]

  # Calculate surplus production
  RPS_pt     <- array(NA, dim = c(nP,nT))

  for( p in 1:nP)
    for( t in (initT_p[p]+1):(nT-1) )
      RPS_pt[p,t] <- R_pt[p,t+1]/B_pt[p,t]



  # Now do the aggregate
  # B_t <- apply( X = B_apt[minAge:nA,,1:nT], FUN = sum, MARGIN = c(3))
  B_t <- apply( X = B_pt[,1:nT], FUN = sum, MARGIN = c(2))
  R_t <- apply( X = R_pt, FUN = sum, MARGIN = c(2))

  RPS_t   <- array(NA, dim = c(nT))
  RPS_t[1:(nT-1)] <- R_t[2:nT]/ B_t[1:(nT - 1)]

  DecadeCols <- scales::viridis_pal(option = "C",alpha = 1)(7)

  ptCols <- c(  rep(DecadeCols[1],10),
                rep(DecadeCols[2],10),
                rep(DecadeCols[3],10),
                rep(DecadeCols[4],10),
                rep(DecadeCols[5],10),
                rep(DecadeCols[6],10),
                rep(DecadeCols[7],9) )

  names(ptCols) <- 1951:2019

  stockNames <- repList$stock

  if( is.null(maxY) )
    maxY <- max(RPS_pt,RPS_t,na.rm =T)

  minY <- 0

  par(  mfrow = c(2,2),
        mar = c(1,1.5,1,1),
        oma = c(4,3,1,1) )
  for( p in 1:nP )
  {
    plot( x = c(0,max(B_pt[p,],na.rm = T)),
          y = c(0,max(RPS_pt[p,],na.rm = T)),
          axes = FALSE, type = "n", las = 1 )
      mfg <- par("mfg")
      # if( mfg[1] == mfg[3])
        axis(side = 1)
      # if( mfg[2] == 1)
        axis(side = 2, las = 1)
      box()
      points( x = B_pt[p,],
              y = RPS_pt[p,],
              bg = ptCols, pch = 21 )
      abline( v = .3*B0_p[p], lty = 3, lwd = .8)
      abline( v = B0_p[p], lty = 3, lwd = .8, col = "red")

      legend( x = "topleft",
              legend = stockNames[p], bty = "n")

  }

  plot( x = c(0,max(B_t,na.rm = T)),
          y = c(0,max(RPS_t,na.rm = T)),
          axes = FALSE, type = "n", las = 1 )
    mfg <- par("mfg")
    # if( mfg[1] == mfg[3])
      axis(side = 1)
    # if( mfg[2] == 1)
      axis(side = 2, las = 1)
    box()
    points( x = B_t,
            y = RPS_t,
            bg = ptCols, pch = 21 )

    abline( v = .3 * sum(B0_p), lty = 3, lwd = .8)
    abline( v = sum(B0_p), lty = 3, lwd = .8, col = "red")

    legend( x = "topright",
            bty = "n",
            legend = c( "1951 - 1960",
                        "1961 - 1970",
                        "1971 - 1980",
                        "1981 - 1990",
                        "1991 - 2000",
                        "2001 - 2010",
                        "2011 - 2019"),
            pch = 21, pt.bg = DecadeCols )
    legend( x = "topleft",
            legend = "Aggregate",
            bty = "n")

    mtext(side = 1, outer = TRUE, text = "Spawning Biomass (kt)",
          line = 1.5)
    mtext(side = 2, outer = TRUE, text = "Recruits per Spawning Stock Biomass (1e6/kt)",
          line = 1.5)
} # END plotRecPerSpawner()

# plotMvsB()
# Plots M as a function of beginning year biomass in each 
# year
plotMvsB <- function( repList=reports,
                      minAge = 3 )
{
  repObj <- repList$repOpt

  # Pull biomass and catch
  B_apt   <- repObj$B_apt
  M_pt    <- repObj$M_apt[minAge,,]
  nA      <- repObj$nA
  nP      <- repObj$nP
  nT      <- repObj$nT
  initT_p <- repObj$tInitModel_p
  B0_p    <- repObj$B0_p
  R0_p    <- repObj$R0_p

  stockNames <- repList$stock

  B_pt    <- apply( X = B_apt[minAge:nA,,], FUN = sum, MARGIN = c(2,3) )

  B_t     <- apply( X = B_pt, FUN = sum, MARGIN = 2 )

  B0minAge_p <- B0_p

  for( p in 1:nP )
  {
    eqbmN_a <- R0_p[p] * repObj$initSurv_ap[,p]
    B0minAge_p[p] <- sum(eqbmN_a[minAge:nA] * repObj$meanWt_ap[minAge:nA,p])
  }

  B0_p <- B0minAge_p

  M_t     <- numeric( length = nT )
  for( t in 1:nT )
    M_t[t] <- sum(M_pt[,t] * B0_p)/sum(B0_p)

  depB_pt <- B_pt
  for( p in 1:nP )
    depB_pt[p,] <- depB_pt[p,] / B0_p[p]
  
  depB_t  <- B_t / sum(B0_p)

  depB_pt[depB_pt == 0] <- NA
  M_pt[M_pt == 0] <- NA

  stockCols <- brewer.pal(n = 4, "Set2")

  plot( x = c(0,max(depB_pt,depB_t,na.rm = T)),
        y = c(0,max(M_pt,M_t,na.rm = T)),
        axes = FALSE, type = "n", las = 1,
        xlab = "", ylab = "" )
    mfg <- par("mfg")
    axis(side = 1)
    axis(side = 2, las = 1)
    box()
    for(p in 1:nP )
    {
      points( x = depB_pt[p,],
              y = M_pt[p,],
              bg = alpha(stockCols[p],alpha = .4), pch = 21, col = NA )
      lines(loess.smooth( x = depB_pt[p,],
                          y = M_pt[p,], span = 1),
            col = stockCols[p], lwd = 2)
    }
    points( x = depB_t[1:nT],
              y = M_t[1:nT],
              bg = alpha(stockCols[4],alpha = .4), pch = 21, col = NA )
    lines(  loess.smooth( x = depB_t[1:nT],
                          y = M_t[1:nT], span = 1),
            col = stockCols[4], lwd = 2)


    abline( v = .3, lty = 3, lwd = .8)
    abline( v = 1, lty = 3, lwd = .8, col = "red")


    legend( x = "topright",
            bty = "n",
            legend = c(stockNames,"Aggregate"),
            pch = 21, pt.bg = stockCols, pt.lwd = 0,
            lty = 1, col = stockCols, lwd = 2 )


    mtext(side = 1, outer = FALSE, text = "3+ Biomass Depletion",
          line = 2)
    mtext(side = 2, outer = FALSE, text = "Natural Mortality (yr)",
          line = 2.5)

}# END plotMvsB()

# plotPropFem()
# Visualises the observed and expected proportion
# females for each gear type, year, and length bin
plotPropFem <- function(  repList = reports,
                          pIdx = 1, gIdx = 6 )
{
  datObj <- repList$data
  repObj <- repList$repOpt

  # Pull data and predictions
  pF_lt     <- datObj$pF_lpgt[,pIdx,gIdx,]
  predPF_lt <- repObj$predPF_lpgt[,pIdx,gIdx,]

  fYear <- repList$fYear


  isPos <- function(x)
  {
    ind <- any(x > 0)
    ind
  }

  yrDat   <- which(apply(X = pF_lt, FUN = isPos, MARGIN = 2) )
  nYrDat  <- length(yrDat) 
  

  # Get first and last length bins
  firstLenBin_lxt  <- repObj$firstLenBin_lxpgt[,,pIdx,gIdx,]
  lastLenBin_lxt   <- repObj$lastLenBin_lxpgt[,,pIdx,gIdx,]

  nL <- repObj$nL
  nT <- repObj$nT
  nP <- repObj$nP

  yrs <- seq(from = fYear, by = 1, length.out = nT)

  nCols <- ceiling(sqrt(nYrDat)) + 1
  nRows <- floor(sqrt(nYrDat))

  lenBinMids_l <- datObj$lenBinMids_l

  par(  mfcol = c(nRows,nCols), 
        mar = c(.1,.1,.1,.1),
        oma = c(3,3,1,3) )

  for( t in 1:nYrDat )
  {
    yIdx <- yrDat[t]
    yr <- yrs[yIdx]
    plot( x = c(0,max(lenBinMids_l)),
          y = c(0,1),
          type = "n", axes = FALSE )
      mfg <- par("mfg")
      if(mfg[1] == mfg[3])
        axis(side =1)
      if(mfg[2] == 1)
        axis(side =2, las = 1 )

      if(mfg[2] == mfg[4] | nYrDat - t < nRows )
        axis( side = 4, las = 1 )
      grid()
      box()

      rect( xleft   = lenBinMids_l - .3, 
            xright  = lenBinMids_l + .3, 
            ybottom = rep(0,nL),
            ytop    = predPF_lt[,yIdx], col = "grey75", border = NA )

      
      datIdx <- which(pF_lt[,yIdx] > 0)
      rect( xleft   = lenBinMids_l[datIdx] - .3, 
            xright  = lenBinMids_l[datIdx] + .3, 
            ybottom = rep(0,length(datIdx)),
            ytop    = predPF_lt[datIdx,yIdx], col = NA, border = "black", lwd=.8 )

      points( x = lenBinMids_l[datIdx], y = pF_lt[datIdx,yIdx],
              col = "darkgreen", pch = 1 )      

      legend(x = "topleft", legend = as.character(yr), bty = "n" )
  }

}



# plotStockDataSummary()
# A visual summary of the data
# available for a given stock
plotStockDataSummary <- function( dataList,
                                  fYear = 1951 )
{
  # Pull data
  I_pgt     <- dataList$I_pgt
  A_apgt    <- dataList$A_apgt
  C_pgt     <- dataList$C_pgt
  combI_pt  <- dataList$combI_pt
  # Mixed catch
  mC_gt  <- dataList$mC_gt

  # First, get dimension
  nT <- dim(I_pgt)[3]
  nP <- dim(I_pgt)[1]
  yrs <- seq( from = fYear, by = 1, length.out = nT+4)

  # Replace -1s with NA
  A_apgt[A_apgt < 0] <- NA
  I_pgt[I_pgt < 0] <- 0
  combI_pt[combI_pt < 0] <- NA

  gearCols <- RColorBrewer::brewer.pal(6, "Dark2")

  A_pgt <- apply( X = A_apgt, FUN = sum, MARGIN = c(2,3,4),
                  na.rm = TRUE )
  A_gt <- apply( X = A_pgt, FUN = sum, MARGIN = c(2,3),
                  na.rm = TRUE )
  sampSize_pg <- apply( X = A_pgt, FUN = sum, MARGIN = c(1,2) )
  sampSize_g  <- apply( X = A_gt, FUN = sum, MARGIN = c(1) )
  
  pchA_pgt <- A_pgt
  pchA_pgt[A_pgt < 200 ] <- 1
  pchA_pgt[A_pgt >= 200 ] <- 16
  pchA_pgt[A_pgt == 0 ] <- NA

  pchA_gt <- A_gt
  pchA_gt[A_gt < 200 ] <- 1
  pchA_gt[A_gt >= 200 ] <- 16
  pchA_gt[A_gt == 0 ] <- NA

  pchI_pgt <- I_pgt
  pchI_pgt[ I_pgt > 0 ] <- 16
  pchI_pgt[ I_pgt == 0 ] <- NA

  pchC_pgt <- C_pgt
  pchC_pgt[ C_pgt > 0 ] <- 16
  pchC_pgt[ C_pgt == 0 ] <- NA

  pchCombI_pt <- combI_pt
  pchCombI_pt[ combI_pt > 0] <- 16
  pchCombI_pt[ combI_pt == 0] <- 1


  I_gt <- apply( X = I_pgt, FUN = sum, MARGIN = c(2,3))
  pchI_gt <- I_gt
  pchI_gt[ I_gt > 0 ] <- 16
  pchI_gt[ I_gt == 0 ] <- NA  

  C_gt <- apply( X = C_pgt, FUN = sum, MARGIN = c(2,3))
  C_gt <- C_gt + mC_gt
  pchC_gt <- C_gt
  pchC_gt[ C_gt > 0 ] <- 16
  pchC_gt[ C_gt == 0 ] <- NA 

  combI_t <- apply( X = combI_pt, FUN = sum, MARGIN = c(2), na.rm = T)
  pchCombI_t <- combI_t
  pchCombI_t[combI_t > 0] <- 16


  ageGears <- 1:3
  catGears <- c(1:3,6)
  idxGears <- 4:5

  # Make a plotting area

  stockLab <- c("C/S","JP/S","Lou")

  yLabs <- c("Ages", "Catch", "Spawn Idx" )

  gearLabs <- c("Reduction","SeineRoe","Gillnet","Surface","Dive","SOK")

  if( any(!is.na(combI_pt)))
    gearLabs <- c("Reduction","SeineRoe","Gillnet","Spawn Idx","SOK")

  par(  mfrow = c(nP+1,1), mar = c(0.5,1,.5,1),
        oma = c(7,6,3,4), xpd = FALSE )
  # Now plot aggregate/mixed data
  plot( x = range(yrs), y = c(0.5,3.5),
          type = "n", axes = FALSE )
    axis( side = 2, las = 1, at = 1:3,
          labels = yLabs)
    mtext( side = 3, text = "SISCA Data Summary", font = 2)

    box()
    grid()  
    abline(h = c(1:2) + .5, lwd = 2)
    rmtext( txt = "Aggregate", cex = 1.5,
            line = .05, outer = TRUE, font = 2 )
    # Plot age comps
    points( x = yrs[1:nT], y = rep(1,nT)-.2,
            pch = pchA_gt[1,],
            col = gearCols[1] )
    points( x = yrs[1:nT], y = rep(1,nT),
            pch = pchA_gt[2,],
            col = gearCols[2] )
    points( x = yrs[1:nT], y = rep(1,nT)+.2,
            pch = pchA_gt[3,],
            col = gearCols[3] )
    # Plot Catch
    points( x = yrs[1:nT], y = rep(2,nT)-.3,
            pch = pchC_gt[1,],
            col = gearCols[1] )
    points( x = yrs[1:nT], y = rep(2,nT)-.1,
            pch = pchC_gt[2,],
            col = gearCols[2] )
    points( x = yrs[1:nT], y = rep(2,nT)+.1,
            pch = pchC_gt[3,],
            col = gearCols[3] )
    points( x = yrs[1:nT], y = rep(2,nT)+.3,
            pch = pchC_gt[6,],
            col = gearCols[6] )


    # Plot Indices
    points( x = yrs[1:nT], y = rep(3,nT)-.2,
            pch = pchI_gt[4,],
            col = gearCols[4] )
    points( x = yrs[1:nT], y = rep(3,nT)+.2,
            pch = pchI_gt[5,],
            col = gearCols[5] )
    points( x = yrs[1:nT], y = rep(3,nT),
            pch = pchCombI_t,
            col = "grey40" )



    text( x = yrs[nT]+4, y = 1 + c(-.2,0,.2), 
          labels = sampSize_g[1:3], cex = .8 )
  for( p in 1:nP )
  {
    plot( x = range(yrs), y = c(0.5,3.5),
          type = "n", axes = FALSE )

      axis( side = 2, las = 1, at = 1:3,
            labels = yLabs)
      box()
      grid()
      rmtext( txt = stockLab[p], cex = 1.5,
              font = 2, line = .05, outer = TRUE )

      # Plot age comps
      points( x = yrs[1:nT], y = rep(1,nT)-.2,
              pch = pchA_pgt[p,1,],
              col = gearCols[1] )
      points( x = yrs[1:nT], y = rep(1,nT),
              pch = pchA_pgt[p,2,],
              col = gearCols[2] )
      points( x = yrs[1:nT], y = rep(1,nT) + .2,
              pch = pchA_pgt[p,3,],
              col = gearCols[3] )
      # Plot catch
      points( x = yrs[1:nT], y = rep(2,nT)-.3,
              pch = pchC_pgt[p,1,],
              col = gearCols[1] )
      points( x = yrs[1:nT], y = rep(2,nT)-.1,
              pch = pchC_pgt[p,2,],
              col = gearCols[2] )
      points( x = yrs[1:nT], y = rep(2,nT)+.1,
              pch = pchC_pgt[p,3,],
              col = gearCols[3] )
      points( x = yrs[1:nT], y = rep(2,nT)+.3,
              pch = pchC_pgt[p,6,],
              col = gearCols[6] )

      # Plot indices
      points( x = yrs[1:nT], y = rep(3,nT)-.2,
              pch = pchI_pgt[p,4,],
              col = gearCols[4] )
      points( x = yrs[1:nT], y = rep(3,nT)+.2,
              pch = pchI_pgt[p,5,],
              col = gearCols[5] )
      points( x = yrs[1:nT], y = rep(3,nT),
              pch = pchCombI_pt[p,],
              col = "grey40" )
      abline(h = c(1:2)+.5, lwd = 2)

  
      text( x = yrs[nT]+4, y = 1 + c(-.2,0,.2), 
            labels = sampSize_pg[p,1:3], cex = .8 )

  }
    axis( side = 1)
    mtext( side = 1, text = "Year", font = 2, outer = TRUE,
            line = 2)

    if( any(!is.na(combI_pt)))
      gearCols <- c(gearCols[1:3],"grey40",gearCols[6])

    par(xpd=TRUE, oma = c(0,0,0,0))
    # plot( x = c(0,1), y = c(0,1), add = TRUE)
    legend( x = "bottom", bty = "n",
            # inset = c(0,-1), xpd = TRUE,
            legend = gearLabs, 
            pch = 16, col = gearCols,
            horiz = TRUE, cex = 1.5, pt.cex = 2)
}


plotSimilarityMetrics <- function(  sensTable, 
                                    winX = c(.3,.6),
                                    winY = c(0,.32) )
{
  AREB0range      <-  c(0,max(sensTable$AREB0))
  AREstatusrange  <-  c(0,max(sensTable$AREstatus))
  MAREBtrange     <-  c(0,max(sensTable$MAREBt))

  nSens <- length(unique(sensTable$Sensitivity))

  cols <- scales::viridis_pal(option = "C",alpha = .8)(nSens)


  sensNames <-  sensTable %>%
                group_by(Sensitivity) %>%
                summarise( nModels = n() ) %>%
                ungroup() %>%
                mutate( sensCol = cols )

  sensTable <-  sensTable %>%
                left_join( sensNames, by = "Sensitivity" ) %>%
                mutate( cex = .8 + 4*( AREB0 - min(AREB0) )/( max(AREB0) - min(AREB0)) ,
                        pch = ifelse(pdHess, 21, 22) )

  simMetrics <- as.matrix(sensTable[,c("AREstatus","AREB0","MAREBt")])
  simMetricCorrelation <- cor(simMetrics)


  par( mfrow = c(2,1), mar = c(1,2,1,2), oma =c(3,3,1,1) )
  
  plot( x = MAREBtrange, y = AREstatusrange,
        type = "n", las = 1 )
    
    points( x = sensTable$MAREBt,
            y = sensTable$AREstatus,
            cex = sensTable$cex,
            bg = sensTable$sensCol,
            pch = sensTable$pch,
            col = "black" )
    points( x = baseCase$MAREBt, y = baseCase$AREstatus,
            col = "black", pch = 2, cex = 1.5)
    # Make a window for the zoomed in plot
    segments( x0 = winX,
              y0 = rep(winY[1],2), y1 = rep(winY[2],2),
              col = "black", lwd = 1)

    segments( x0 = rep(winX[1],2), x1 = rep(winX[2],2),
              y0 = winY,
              col = "black", lwd = 1)
                
    plot( x = winX, y = winY,
          type = "n", las = 1 )
      # mtext( side = 2, text = "ARE(B2019/B0)")
      points( x = sensTable$MAREBt,
              y = sensTable$AREstatus,
              cex = sensTable$cex,
              bg = sensTable$sensCol,
              pch = sensTable$pch,
              col = "black" )
      points( x = baseCase$MAREBt, y = baseCase$AREstatus,
              col = "black", pch = 2, cex = 1.5)
      legend( x = "topright", bty = "n",
            legend = c(sensNames$Sensitivity[-nSens],"Base Model"),
            pt.bg = c(sensNames$sensCol[-nSens],NA),
            pch = c(rep(21,nSens-1),2),
            col = c(rep("black",nSens-1),"black"),
            pt.cex = 1.5 )
  
  mtext( side = 2, text = expression(rho(B[2019]/B[0])), outer = TRUE,
          line = 1.5) 
  mtext( side = 1, text = expression(rho(B[t])),
          outer = TRUE, line = 2 )
}

plotSimMetricCorrelation <- function(  sensTable )
{
  require(corrplot)
  simMetrics <- as.matrix(sensTable[,c("AREstatus","AREB0","MAREBt")])
  colnames(simMetrics) <- c("rho(Status)","rho(B0)","rho(Bt)")
  # colnames(simMetrics) <- c(  expression( rho(B[2019]/B[0])) ,
  #                             expression( rho(B[0]) ),
  #                             expression( rho(B[t]) ) )
  simMetricCorrelation <- cor(simMetrics)

  Dmetrics <- sensTable %>%
              dplyr::select(  Equal,
                              Status = StatusDom,
                              B0 = B0dom,
                              Bt = BtDom,
                              ageSE = meanAgeSE,
                              idxSE = meanIdxSE ) %>%
              as.matrix()


  DmetricCorrelation <- cor(Dmetrics)

  par( mfrow = c(2,1), mar = c(0,0,0,0), oma = c(0.,0.,0,0) )

  corrplot.mixed( simMetricCorrelation, lower.col = "black", number.cex = 1 )
  corrplot.mixed( DmetricCorrelation, lower.col = "black", number.cex = 1 )
}

# plotSpawnIdxSP()
# Plot of spawn indices and surplus production,
# reproducing Kronlund et al Figure 22 - 29
plotSpawnIdxSP <- function( repList = reports,
                            stock = "Agg",
                            lowSSBquant = .3 )
{
  # First, make SP tables
  SPtables <- makeSPtables(repList, save = FALSE)

  # Pull table for stock
  spTable <- SPtables[[stock]]

  repObj <- repList$repOpt

  # Pull biomass and catch
  nA      <- repObj$nA
  nP      <- repObj$nP
  nT      <- repObj$nT

  fYear   <- repList$fYear
  yrs     <- spTable$Year

  SSBquant  <- quantile(spTable$SSB, probs = lowSSBquant, na.rm = T)
  spTable <- spTable %>% mutate( ptCol = ifelse( SSB <= SSBquant, "grey40","white") )


  par(mfrow = c(3,1), oma = c(3,4,2,2), mar = c(.1,.1,.1,.1))

  # Plot spawn indices
  plot( x = range(yrs), y = range(spTable$spawnIdx, na.rm = T),
        type = "n", axes = FALSE )
    axis( side = 2, las = 1 )
    grid()
    box()
    mtext( side = 2, line = 2.5, text = "Spawn Index (kt)" )
    lines( x = yrs, y = spTable$spawnIdx, lwd = 2 )
    points( x = yrs, y = spTable$spawnIdx, pch = 21,
            bg = spTable$ptCol, cex = 1.5 )
    lines( x = yrs[-(1:2)], y = pracma::movavg(x = spTable$spawnIdx, n = 3)[-(1:2)],
            col = "steelblue", lwd = 3 )
    abline(v = 1987.5, lty = 2, lwd = 1)
    legend( x = "topright",
            legend = c( paste("SSB < ", 100*lowSSBquant, "th percentile", sep = ""),
                        "3 Year Moving Avg"),
            pch = c(21,NA),
            col = c("black","steelblue"),
            pt.bg = c("grey40",NA),
            lty = c(NA,1),
            lwd = c(NA,3) )

  plot( x = range(yrs), y = range(spTable$SurplusProd, na.rm = T),
        type = "n", axes = FALSE )
    axis( side = 2, las = 1 )
    grid()
    box()
    mtext( side = 2, line = 2.5, text = "Surplus Production (kt)" )
    abline( h = 0, lty = 2, lwd = .8)
    lines( x = yrs, y = spTable$SurplusProd, lwd = 2 )
    points( x = yrs, y = spTable$SurplusProd, pch = 21,
            bg = spTable$ptCol, cex = 1.5 )
    lines( x = yrs[-(1:2)], y = pracma::movavg(x = spTable$SurplusProd, n = 3)[-(1:2)],
            col = "steelblue", lwd = 3 )
    abline(v = 1987.5, lty = 2, lwd = 1)

  plot( x = range(yrs), y = range(spTable$perCapSurplusProd, na.rm = T),
        type = "n", axes = FALSE )
    axis( side = 2, las = 1 )
    axis( side = 1)
    grid()
    box()
    abline( h = 0, lty = 2, lwd = .8)
    mtext( side = 2, line = 2.5, text = "Surplus Production Rate" )
    lines( x = yrs, y = spTable$perCapSurplusProd, lwd = 2 )
    points( x = yrs, y = spTable$perCapSurplusProd, pch = 21,
            bg = spTable$ptCol, cex = 1.5 )
    lines(  x = yrs[-(1:2)], 
            y = pracma::movavg(x = spTable$perCapSurplusProd, n = 3)[-(1:2)],
            col = "steelblue", lwd = 3 )
    abline(v = 1987.5, lty = 2, lwd = 1)


} # END plotSpawnIdxSP

# plotSPphasePlots()
# Plots SP phase plots similar to Kronlund et al
plotSPphasePlots <- function( repList = reports,
                              stock = "Agg",
                              depLines = c(.1,.25,.3,.5,.6,1.0),
                              depCols  = c("red","grey30") )
{
  # First, make SP tables
  SPtables <- makeSPtables(repList, save = FALSE)

  # Pull table for stock
  spTable <- SPtables[[stock]]

  repObj <- repList$repOpt

  # Pull biomass and catch
  nA      <- repObj$nA
  nP      <- repObj$nP
  nT      <- repObj$nT

  if( stock == "Agg" )
    B0 <- sum(repObj$B0_p)
  else
    B0 <- repObj$B0_p[which(names(SPtables == stock))]

  fYear   <- repList$fYear
  yrs     <- spTable$Year

  surfYrs <- 1951:1987
  diveYrs <- 1988:2019

  nSurf <- length(surfYrs)
  nDive <- length(diveYrs)

  surfTable <- spTable %>% filter(Year %in% surfYrs)
  diveTable <- spTable %>% filter(Year %in% diveYrs)

  depLinePalette  <- colorRampPalette(depCols)
  depLineCols     <- depLinePalette(length(depLines)) 

  ptCols <- scales::viridis_pal(option = "C",alpha = 1)(nT)
  names(ptCols) <- fYear:2019

  # 2x2 plot, 
  # left column is 1951-1987, right is 1988+
  # top row is abs SP, bottom is per cap
  par(  mfrow = c(2,2),
        mar = c(.1,1.5,.1,1.5),
        oma = c(4,2,1,1) )
  
  plot( x = c(0,max(surfTable$SSB)),
        y = range(surfTable$SurplusProd),
        axes = FALSE, type = "n", las = 1 )
    axis(side = 2, las = 1)
    box()
    abline(h = 0, lty = 2 )
    abline(v = depLines * B0, col = depLineCols, lty = 2 )
    mtext( side = 2, text = "Surplus Production (kt)", line = 2)
    # Put in arrows first

    arrows( x0 = surfTable$SSB[1:(nSurf-1)],
            x1 = surfTable$SSB[2:(nSurf)],
            y0 = surfTable$SurplusProd[1:(nSurf-1)],
            y1 = surfTable$SurplusProd[2:(nSurf)],
            col = "grey75", lwd = .8, length = .1 )
    points( x = surfTable$SSB,
            y = surfTable$SurplusProd,
            bg = ptCols[1:nSurf], pch = 21 )
    text( x = surfTable$SSB[c(1,nSurf)] ,
          y = surfTable$SurplusProd[c(1,nSurf)]+ c(-5,5),
          labels = surfTable$Year[c(1,nSurf)] )

    legend( x = "topright",
            bty = "n",
            legend = c( "1951",
                        "1987"),
            pch = 21, pt.bg = ptCols[c(1,nSurf)] )
    mtext( side = 3, text="1951-1987", font = 2)

  plot( x = c(0,max(diveTable$SSB)),
        y = range(diveTable$SurplusProd, na.rm = T),
        axes = FALSE, type = "n", las = 1 )
    axis(side = 2, las = 1)
    box()
    abline(h = 0, lty = 2 )
    abline(v = depLines * B0, col = depLineCols, lty = 2 )
    # Put in arrows first

    arrows( x0 = diveTable$SSB[1:(nDive-1)],
            x1 = diveTable$SSB[2:(nDive)],
            y0 = diveTable$SurplusProd[1:(nDive-1)],
            y1 = diveTable$SurplusProd[2:(nDive)],
            col = "grey75", lwd = .8, length = .1 )
    points( x = diveTable$SSB,
            y = diveTable$SurplusProd,
            bg = ptCols[(nSurf+1):nT], pch = 21 )
    text( x = diveTable$SSB[c(1,nDive-1)] + c(-1,0),
          y = diveTable$SurplusProd[c(1,nDive-1)] + c(0,0),
          labels = diveTable$Year[c(1,nDive-1)] )

    legend( x = "topright",
            bty = "n",
            legend = c( "1988",
                        "2018"),
            pch = 21, pt.bg = ptCols[c(nSurf+1,nT-1)] )
    mtext( side = 3, text="1988-2019", font = 2)


  plot( x = c(0,max(surfTable$SSB)),
        y = range(surfTable$perCapSurplusProd, na.rm = T),
        axes = FALSE, type = "n", las = 1 )
    axis(side = 1)
    axis(side = 2, las = 1)
    box()
    abline(h = 0, lty = 2 )
    abline(v = depLines * B0, col = depLineCols, lty = 2 )
    mtext( side = 2, text = "Surplus Production Rate", line = 2)
    # Put in arrows first

    arrows( x0 = surfTable$SSB[1:(nSurf-1)],
            x1 = surfTable$SSB[2:(nSurf)],
            y0 = surfTable$perCapSurplusProd[1:(nSurf-1)],
            y1 = surfTable$perCapSurplusProd[2:(nSurf)],
            col = "grey75", lwd = .8, length = .1 )
    points( x = surfTable$SSB,
            y = surfTable$perCapSurplusProd,
            bg = ptCols[1:nSurf], pch = 21 )
    text( x = surfTable$SSB[c(1,nSurf)],
          y = surfTable$perCapSurplusProd[c(1,nSurf)] + c(.1,.1),
          labels = surfTable$Year[c(1,nSurf)] )

    legend( x = "topright",
            bty = "n",
            legend = paste(depLines,"B0",sep = ""),
            lty = 2, col = depLineCols )

  plot( x = c(0,max(diveTable$SSB)),
        y = range(diveTable$perCapSurplusProd, na.rm = T),
        axes = FALSE, type = "n", las = 1 )
    axis(side = 1)
    axis(side = 2, las = 1)
    box()
    abline(h = 0, lty = 2 )
    abline(v = depLines * B0, col = depLineCols, lty = 2 )
    # Put in arrows first

    arrows( x0 = diveTable$SSB[1:(nDive-1)],
            x1 = diveTable$SSB[2:(nDive)],
            y0 = diveTable$perCapSurplusProd[1:(nDive-1)],
            y1 = diveTable$perCapSurplusProd[2:(nDive)],
            col = "grey75", lwd = .8, length = .1 )
    points( x = diveTable$SSB,
            y = diveTable$perCapSurplusProd,
            bg = ptCols[(nSurf+1):nT], pch = 21 )
    text( x = diveTable$SSB[c(1,nDive-1)],
          y = diveTable$perCapSurplusProd[c(1,nDive-1)] + .1,
          labels = diveTable$Year[c(1,nDive-1)] )



    mtext(side = 1, outer = TRUE, text = "Spawning Stock Biomass (kt)",
          line = 2.5)

}
