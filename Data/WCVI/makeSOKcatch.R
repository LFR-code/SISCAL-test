# Read length comps
library(tidyverse)

area23 <- read.csv("SOKWCVI_Area23.csv") %>%
          mutate( AreaName = "Area23",
                  Area = 3 )

area24 <- read.csv("SOKWCVI_Area24.csv") %>%
          mutate( AreaName = "Area24",
                  Area = 2 )
area25 <- read.csv("SOKWCVI_Area25.csv") %>%
          mutate( AreaName = "Area25",
                  Area = 1 )


allAreas <- area25 %>%
            rbind( area24 ) %>%
            rbind( area23 ) %>%
            mutate( gearID = 6,
                    Harvest = Harvest/1e6,
                    Type = 2 )



aggSOB <- allAreas %>%
          group_by( Year ) %>%
          summarise(  Value = sum(Harvest),
                      Gear = mean(gearID),
                      Area = min(Area),
                      Type = mean(Type) ) %>%
          mutate(Stock = "Agg" )

catchData <- read.csv("catchData.csv", header = TRUE)

catchData <- rbind(catchData,aggSOB)

write.csv(catchData, file = "catchDataWSOK.csv")