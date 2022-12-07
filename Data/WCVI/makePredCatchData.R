# Read length comps
hakeList <- readRDS("wcviHakeConsumption.rds")
load("MMpredMod.RData")

# We need to assign fleet numbers
# fleetNums <- c( hakeLt50 = 7, hakeGt50 = 8, HS = 9, SSL = 10, HB = 11, GW = 12 )
fleetNums <- c( hakeLt50 = 7, hakeGt50 = 8, HS = 9, SSL = 10, hbW = 11, hbS=12, GW = 13 )
nPreds <- length(fleetNums)
catchType <- c(rep(1,nPreds-1),2)

yrs <- 1951:2019
nT <- length(yrs)


C_tp <- array(0, dim = c(nT,nPreds))

gwCat <- 1.5/1e3

predCatchTable <- matrix( NA, nrow = nT * nPreds, ncol = 6 )
colnames(predCatchTable) <- c("Year","Gear","Area","Type","Value","Stock")
predCatchTable <- as.data.frame(predCatchTable)
predCatchTable$Stock <- "Agg"



for( k in 1:nPreds )
{
  pred <- names(fleetNums)[k]

  if( pred == "hakeLt50" )
  {
    C_t <- hakeList$C_lt[1,]
    hakeYrs <- names(C_t)
    addC_t <- rep(0,length = nT)
    names(addC_t) <- yrs
    addC_t[hakeYrs] <- C_t
    addC_t[!yrs %in% as.numeric(hakeYrs)] <- C_t[1]
    C_t <- addC_t
  }

  if( pred == "hakeGt50" )
  {
    C_t <- hakeList$C_lt[2,]
    hakeYrs <- names(C_t)
    addC_t <- rep(0,length = nT)
    names(addC_t) <- yrs
    addC_t[hakeYrs] <- C_t
    addC_t[!yrs %in% as.numeric(hakeYrs)] <- C_t[1]
    C_t <- addC_t
  }

  if(pred %in% c("hbS","hbW","SSL","HS"))
  {  
    if( pred == "hbS" )
      mmPredIdx <- which(predMod$ctrl$predNames_p=='humpbackS')

    if( pred == "hbW" )
      mmPredIdx <- which(predMod$ctrl$predNames_p=='humpbackW')

    if( pred == "SSL" )
      mmPredIdx <- which(predMod$ctrl$predNames_p=='stellarSeaLion')

    if( pred == "HS" )
      mmPredIdx <- which(predMod$ctrl$predNames_p=='harbourSeal')
  
  C_t <- predMod$C_tp[2:70,mmPredIdx]/1e3
  }

  if( pred == "GW" )
    C_t <- rep(gwCat,nT)

  rowIdx <- (k-1) * nT + 1:nT

  predCatchTable[rowIdx,"Year"] <- yrs
  predCatchTable[rowIdx,"Gear"] <- fleetNums[k]
  predCatchTable[rowIdx,"Area"] <- 1
  predCatchTable[rowIdx,"Type"] <- catchType[k]
  predCatchTable[rowIdx,"Value"] <- round(C_t,3)

  C_tp[,k] <- round(C_t,3)

}

write.csv(predCatchTable, file = "catchDataMixed.csv")

predC_t <- apply(X = C_tp, FUN = sum, MARGIN = 1)

# Now plot stacked bar plot of catch
predCols <- RColorBrewer::brewer.pal(nPreds,"Dark2")
plot(x = range(yrs), y = c(0,max(predC_t)), las = 1,
      type = "n", xlab = "Year", ylab = "Predator Consumption (kt)" )
lines(x = yrs, y = predC_t, col = "black", lwd = 3)
for( p in 1:nPreds )
{
  lines( x = yrs, y = C_tp[,p], col = predCols[p], lwd = 2)
  # rect( xleft= yrs - .3, xright = yrs + .3, 
  #       ybottom = baseC_t, ytop = baseC_t + C_tp[,p],
  #       col = predCols[p], border = NA )
  
}
legend( x = "topleft", bty = "n",
        legend = c("Total",names(fleetNums)),
        col = c("black",predCols),
        lwd = c(3,rep(2,nPreds)) )
