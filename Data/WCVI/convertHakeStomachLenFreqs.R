# Read length comps
hakeDietLenComps <- read.csv("hakeDietLenComps.csv")

# Get inverse age length key
invAgeLenKey_al <- read.csv("./wcviBio/HerringInvAgeLenKey.csv", row.names = 1)
invAgeLenKey_al <- as.matrix(invAgeLenKey_al[,-ncol(invAgeLenKey_al)])


hakeDietAgeComps <- array(0, dim = c(nrow(hakeDietLenComps),11) )
colnames(hakeDietAgeComps) <- c("Year",paste0("a",1:10))
hakeDietAgeComps <- as.data.frame(hakeDietAgeComps)
hakeDietAgeComps$Year <- hakeDietLenComps$Year

lenCompsMat <- as.matrix(hakeDietLenComps[,-1])

hakeDietAgeComps[,2:11] <- t(invAgeLenKey_al %*% t(lenCompsMat))
for( rIdx in 1:5 )
  hakeDietAgeComps[rIdx,2:11] <- hakeDietAgeComps[rIdx,2:11] / sum(hakeDietAgeComps[rIdx,2:11])


write.csv(hakeDietAgeComps, file = "hakeDietAgeComps.csv")