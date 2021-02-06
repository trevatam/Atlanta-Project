########################################################################
###
### ATLANTA PROJECT - POISSON MODELS & CAR MODELS
### LAST UPDATED FEB 5, 2021
### RUN SUCCESSFULLY IN R VERSION 4.0.0
###
########################################################################

rm(list = ls())
setwd("C:/Users/treva/Documents/Penn/Project_MA paper/MA_Publish/")

library(data.table)
library(dplyr) #data manipulation
library(rgdal) # needed for readOGR, and working with projections
library(ggplot2)
library(classInt)
library(car)
library(spdep)
library(RColorBrewer)
library(StatDA)
library(rgdal)
library(maptools)
library(MASS)
library(lmtest)
library(sp)
library(spatstat)
library(nortest)
library(sphet)
library(spatialreg)
library(stargazer)
library(CARBayes)
library(xtable)

########################################################################
###
### Load Data: Atlanta MSA Shapefile and project
###
########################################################################
set.seed(2020)

#load Georgia Tracts
GAtracts <- readOGR(dsn=getwd(), layer="GA_tracts")
#GAtracts <- st_read("GA_tracts/GA_tracts.shp")
#Counties in Atlanta MSA (except Morgan County - added after 2010)
CountyFIPS <- c('013', '015', '035', '045', '057', '063', '067', '077', '085', '089', '097', '113', '117', '121', '135', '143', '149', '151', '159', '171', '199', '217', '223', '227', '231', '247', '255', '297')
GAtracts <- st_as_sf(GAtracts)
ATLtracts <- GAtracts %>%
  filter(COUNTYFP10 %in% CountyFIPS)
ATLtracts <- as_Spatial(ATLtracts)

class(ATLtracts) # it should be a Spatial Polygons Data Frame

# check the projection information (will be the .prj shapefile)
is.projected(ATLtracts)
# project the data
ATLtracts <- spTransform(ATLtracts, CRS("+init=EPSG:26967"))
# the projections use the ESPG codes. They are listed on http://spatialreference.org/ref/
proj4string(ATLtracts)
prj <- proj4string(ATLtracts) # save the projection file


##### load data files from "Code_data" ## ORIGINAL DATA

data.long <- read.csv("data.long.csv")
data.long <- data.long[,-1]
data.long$year1 <- data.long$year
data.long$year <- ifelse(data.long$year1 == 2010, 1, 0) # only data from 2000 and 2010 because of lag
data.wide <- read.csv("data.wide.csv") # 936 tracts
data.wide <- data.wide[,-1]

#### make integer variables ####

data.long$WH <- round(data.long$WH)
data.long$BL <- round(data.long$BL)
data.long$ASN <- round(data.long$ASN)
data.long$LTN <- round(data.long$LTN)
data.long$POP <- round(data.long$POP)

data.long$lag.WH <- round(data.long$lag.WH)
data.long$lag.BL <- round(data.long$lag.BL)
data.long$lag.ASN <- round(data.long$lag.ASN)
data.long$lag.LTN <- round(data.long$lag.LTN)
data.long$lag.POP <- round(data.long$lag.POP)

##### merge data

wide.sp <- sp::merge(ATLtracts, data.wide, by.x="GEOID10", by.y="GEOID", all.x=FALSE)
writeOGR(wide.sp, dsn = "C:/Users/treva/Documents/Penn/Project_MA paper/MA_Publish/Final_Shpfiles", layer = "fulldata", driver = "ESRI Shapefile" )

long.sp <- sp::merge(ATLtracts, data.long, by.x="GEOID10", by.y="GEOID", all.x=FALSE, duplicateGeoms = TRUE)

#### data for GIS ######

dat1990 <- data.long %>%
  filter(year1 == 2000) %>%
  dplyr::select(GEOID, GEOQNAME, TRACTCE,
                lag.WH, lag.BL, lag.ASN, lag.LTN, 
                lag.WHPCT, lag.BLPCT, lag.ASNPCT, lag.LTNPCT,
                lag.HHINC, lag.B20YR, lag.B20YRPCT) %>%
  rename(WH = lag.WH, BL = lag.BL, ASN = lag.ASN, LTN = lag.LTN, 
         WH.PCT = lag.WHPCT, BL.PCT = lag.BLPCT, ASN.PCT = lag.ASNPCT, LTN.PCT = lag.LTNPCT,
         newHHINC = lag.HHINC, B20YR = lag.B20YR, B20YR.PCT = lag.B20YRPCT)
dat2000 <- data.long %>%
  filter(year1 == 2000) %>%
  dplyr::select(GEOID, GEOQNAME, TRACTCE,
                WH, BL, ASN, LTN,
                WH.PCT, BL.PCT, ASN.PCT, LTN.PCT,
                newHHINC, B20YR, B20YR.PCT)
dat2010 <- data.long %>%
  filter(year1 == 2010) %>%
  dplyr::select(GEOID, GEOQNAME, TRACTCE,
                WH, BL, ASN, LTN,
                WH.PCT, BL.PCT, ASN.PCT, LTN.PCT,
                newHHINC, B20YR, B20YR.PCT)

dat1990$year <- 1990
dat2000$year <- 2000
dat2010$year <- 2010

dat.long <- bind_rows(dat1990, dat2000, dat2010)
dat.long <- as.data.table(dat.long)

cols <- colnames(dat.long)
remove <- c("GEOID", "GEOQNAME", "TRACTCE","year")
cols_r <- cols[! cols %in% remove]
dat.wide <- dcast(dat.long, GEOID+GEOQNAME+TRACTCE~year, value.var=c(cols_r))
dat.wide$GEOID <- as.factor(dat.wide$GEOID)

dat.sp <- sp::merge(ATLtracts, dat.wide, by.x="GEOID10", by.y="GEOID", all.x=FALSE)
library(raster)
#shapefile(dat.sp, filename = "Final_Shpfiles/Maps/datafull.shp")

####### create segregation index for MSA #####

library(seg)
seg1990 <- dat.sp[c("WH_1990", "BL_1990", "ASN_1990", "LTN_1990")]
seg2000 <- dat.sp[c("WH_2000", "BL_2000", "ASN_2000", "LTN_2000")]
seg2010 <- dat.sp[c("WH_2010", "BL_1990", "ASN_2010", "LTN_2010")]

ATLseg <- function(x){ #x has to be some kind of dataset with spatial
  seg.dat.wb <- x[c(1,2)]
  seg.dat.wa <- x[c(1,3)]
  seg.dat.wl <- x[c(1,4)]
  seg.dat.ba <- x[c(2,3)]
  seg.dat.bl <- x[c(2,4)]
  seg.dat.la <- x[c(3,4)]
  out.wb <- spseg(seg.dat.wb, method = "all")
  out.wa <- spseg(seg.dat.wa, method = "all")
  out.wl <- spseg(seg.dat.wl, method = "all")
  out.ba <- spseg(seg.dat.ba, method = "all")
  out.bl <- spseg(seg.dat.bl, method = "all")
  out.la <- spseg(seg.dat.la, method = "all")
  
  list("wH-BL dissim" = out.wb@d,
       "WH-ASN dissim" = out.wa@d,
       "WH-LTN dissim" = out.wl@d,
       "BL-ASN dissim" = out.ba@d,
       "BL-LTN dissim" = out.bl@d,
       "ASN-LTN dissim" = out.la@d)
}

a <- unlist(ATLseg(seg1990))
b <- unlist(ATLseg(seg2000))
c <- unlist(ATLseg(seg2010))

x <- cbind(a,b,c)
colnames(x) <- c("1990", "2000", "2010")
rownames(x) <- c("White-Black", "White-Asian", "White-Latinx", "Black-Asian", "Black-Latinx", "Asian-Latinx")
print(xtable(x, 
             digits = 2,
             caption = "Segregation Indices in Atlanta, 1990-2010"),
      include.rownames = T,
      include.colnames = T,
      booktabs = T,
      caption.placement = "top",
      auto=T)

########################################################################
###
### Initial Data Set up, weights setting, and EDA
###
########################################################################

### separate year data sets ###

data1 <- data.long %>%
  filter(year1 ==2000)
data1 <- sp::merge(ATLtracts, data1, by.x="GEOID10", by.y="GEOID", all.x=FALSE)

data2 <- data.long %>%
  filter(year1 ==2010)
data2 <- sp::merge(ATLtracts, data2, by.x="GEOID10", by.y="GEOID", all.x=FALSE)

######### CREATE A WEIGHTS MATRIX for Global Moran's I:#########

nb <- poly2nb(wide.sp)
## Queen's Continguity
nb.Q1 <- nb2listw(nb, zero.policy = F) # first order, queens contiguity
W.qmat <- nb2mat(nb, style = "B") # weights matrix

########## Look at correlation #####################

library(GGally)
model1 <- glm(WH ~ lag.BLPCT + lag.ASNPCT + lag.LTNPCT + lag.MHSPCT + lag.HHINC + lag.OWNPCT + lag.MULTIPCT + lag.VACANTPCT + B20YR.PCT + lag.MEDVAL, family="poisson", data = data1)
ggpairs(model1)

model2 <- glm(ASN ~ offset(log(lag.ASN+1)) + lag.BLPCT + lag.ASNPCT + lag.LTNPCT + lag.MHSPCT + lag.HHINC + lag.OWNPCT + lag.MULTIPCT + lag.VACANTPCT + B20YR.PCT+ lag.MEDVAL, family="poisson", data = data2)
ggpairs(model2)

#################################################################
###
### POISSON STOCK MODELS
###
#################################################################

############################
# white population stock reg
############################

# Global Moran's I
moran.test(wide.sp$WH_2000, nb.Q1, zero.policy=T) # 0.62
moran.test(wide.sp$WH_2010, nb.Q1, zero.policy=T) # 0.64

stockmod.w1 <- glm(WH ~ year + lag.BLPCT + lag.ASNPCT + lag.LTNPCT + lag.MHSPCT + lag.HHINC, data=data.long, family = poisson)

stockmod.w2 <- glm(WH ~ year + lag.BLPCT + lag.ASNPCT + lag.LTNPCT + lag.MHSPCT + lag.HHINC + lag.MULTIPCT + lag.VACANTPCT + B20YR.PCT, data=data.long, family = poisson)

stockmod.w3  <- glm(WH ~ year*lag.BLPCT + lag.ASNPCT*year + lag.LTNPCT*year + lag.MHSPCT*year + lag.HHINC*year + lag.MULTIPCT*year + lag.VACANTPCT*year + B20YR.PCT*year, data=data.long, family = poisson)

summary(stockmod.w3)
exp(stockmod.w3[[1]])

#split up into 2 years to get Moran's I

stockmod.w3.t1 <- glm(WH ~ lag.BLPCT + lag.ASNPCT + lag.LTNPCT + lag.MHSPCT + lag.HHINC + lag.MULTIPCT + lag.VACANTPCT + B20YR.PCT, data=data1, family = poisson)

stockmod.w3.t2 <- glm(WH ~ lag.BLPCT + lag.ASNPCT + lag.LTNPCT + lag.MHSPCT + lag.HHINC + lag.MULTIPCT + lag.VACANTPCT + B20YR.PCT, data=data2, family = poisson)

# Best Model Fit (AIC)
stockmod.w1[[11]] # 973,783
stockmod.w2[[11]] # 911,705
stockmod.w3[[11]] # 879,878.2

#split up by year
stockmod.w3.t1[[11]] # 421,486.7
stockmod.w3.t2[[11]] # 451,039.8

### Moran's I on residuals per each year
# 2000
yhat1 <- stockmod.w3.t1$fitted.values
y1 <- data1$WH
plot(yhat1, y1)
lines(lowess(y1 ~ yhat1), lwd=2, col=4)
abline(lm(y1 ~ yhat1), lwd=2, col=2)
# 2010
yhat2 <- stockmod.w3.t2$fitted.values
y2 <- data2$WH
plot(yhat2, y2)
lines(lowess(y2 ~ yhat2), lwd=2, col=4)
abline(lm(y2 ~ yhat2), lwd=2, col=2)

# Run the Moran test on the residuals:
moran.test(stockmod.w3.t1$residuals, nb.Q1, alternative="two.sided") # 0.368
moran.test(stockmod.w3.t2$residuals, nb.Q1, alternative="two.sided") # 0.369

########################
# Black population stock
########################

#Global Moran's I
moran.test(wide.sp$BL_2000, nb.Q1, zero.policy=T) # 0.77
moran.test(wide.sp$BL_2010, nb.Q1, zero.policy=T) # 0.63

stockmod.b1 <- glm(BL ~ year + lag.BLPCT + lag.ASNPCT + lag.LTNPCT + lag.MHSPCT + lag.HHINC, data=data.long, family = poisson)

stockmod.b2 <- glm(BL ~ year + lag.BLPCT + lag.ASNPCT + lag.LTNPCT + lag.MHSPCT + lag.HHINC + lag.MULTIPCT + lag.VACANTPCT + B20YR.PCT, data=data.long, family = poisson)

stockmod.b3 <- glm(BL ~ year*lag.BLPCT + lag.ASNPCT*year + lag.LTNPCT*year + lag.MHSPCT*year + lag.HHINC*year + lag.MULTIPCT*year + lag.VACANTPCT*year + B20YR.PCT*year + lag.MEDVAL*year, data=data.long, family = poisson)

#split up into 2 years to get Moran's I

stockmod.b3.t1 <- glm(BL ~ lag.BLPCT + lag.ASNPCT + lag.LTNPCT + lag.MHSPCT + lag.HHINC + lag.MULTIPCT + lag.VACANTPCT + B20YR.PCT, data=data1, family = poisson)

stockmod.b3.t2 <- glm(BL ~ lag.BLPCT + lag.ASNPCT + lag.LTNPCT + lag.MHSPCT + lag.HHINC + lag.MULTIPCT + lag.VACANTPCT + B20YR.PCT, data=data2, family = poisson)

# Best Model Fit (AIC)
reg1.b[[11]] # 1584409
reg2.b[[11]] # 1373926
reg3.b[[11]] # 1222036

#split up by year
reg3.b.1[[11]] # 525331.9
reg3.b.2[[11]] # 696704.3

### Moran's I on residuals per each year
# 2000
yhat1 <- stockmod.b3.t1$fitted.values
y1 <- data1$BL
plot(yhat1, y1)
lines(lowess(y1 ~ yhat1), lwd=2, col=4)
abline(lm(y1 ~ yhat1), lwd=2, col=2)
# 2010
yhat2 <- stockmod.b3.t2$fitted.values
y2 <- data2$BL
plot(yhat2, y2)
lines(lowess(y2 ~ yhat2), lwd=2, col=4)
abline(lm(y2 ~ yhat2), lwd=2, col=2)

# Run the Moran test on the residuals:
moran.test(stockmod.b3.t1$residuals, nb.Q1, alternative="two.sided") # 0.58
moran.test(stockmod.b3.t2$residuals, nb.Q1, alternative="two.sided") # 0.47

########################
# Asian population stock
########################

# Global Moran's I
moran.test(wide.sp$ASN_2000, nb.Q1, zero.policy=T) # 0.67
moran.test(wide.sp$ASN_2010, nb.Q1, zero.policy=T) # 0.56

stockmod.a1 <- glm(ASN ~ year + lag.BLPCT + lag.ASNPCT + lag.LTNPCT + lag.MHSPCT + lag.HHINC, data=data.long, family = poisson)

stockmod.a2 <- glm(ASN ~ year + lag.BLPCT + lag.ASNPCT + lag.LTNPCT + lag.MHSPCT + lag.HHINC + lag.MULTIPCT + lag.VACANTPCT + B20YR.PCT, data=data.long, family = poisson)

stockmod.a3 <- glm(ASN ~ year*lag.BLPCT + lag.ASNPCT*year + lag.LTNPCT*year + lag.MHSPCT*year + lag.HHINC*year + lag.MULTIPCT*year + lag.VACANTPCT*year + B20YR.PCT*year + lag.MEDVAL*year, data=data.long, family = poisson)

#split up into 2 years to get Moran's I

stockmod.a3.t1<- glm(ASN ~ lag.BLPCT + lag.ASNPCT + lag.LTNPCT + lag.MHSPCT + lag.HHINC + lag.MULTIPCT + lag.VACANTPCT + B20YR.PCT, data=data1, family = poisson)

stockmod.a3.t2 <- glm(ASN ~ lag.BLPCT + lag.ASNPCT + lag.LTNPCT + lag.MHSPCT + lag.HHINC + lag.MULTIPCT + lag.VACANTPCT + B20YR.PCT, data=data2, family = poisson)

# Best Model Fit (AIC)
reg1.a[[11]] # 253,019.9
reg2.a[[11]] # 229,248.3
reg3.a[[11]] # 217,109.1

#split up by year
reg3.a.1[[11]] # 83,590.56
reg3.a.2[[11]] # 133,518.5

### Moran's I on residuals per each year
# 2000
yhat1 <- stockmod.a3.t1$fitted.values
y1 <- data1$ASN
plot(yhat1, y1)
lines(lowess(y1 ~ yhat1), lwd=2, col=4)
abline(lm(y1 ~ yhat1), lwd=2, col=2)
# 2010
yhat2 <- stockmod.a3.t2$fitted.values
y2 <- data2$ASN
plot(yhat2, y2)
lines(lowess(y2 ~ yhat2), lwd=2, col=4)
abline(lm(y2 ~ yhat2), lwd=2, col=2)

# Run the Moran test on the residuals:
moran.test(stockmod.a3.t1$residuals, nb.Q1, alternative="two.sided") # 0.127
moran.test(stockmod.a3.t2$residuals, nb.Q1, alternative="two.sided") # 0.331

########################
# Latinx population stock
########################

# Global Moran's I
moran.test(wide.sp$LTN_2000, nb.Q1, zero.policy=T) # 0.57
moran.test(wide.sp$LTN_2010, nb.Q1, zero.policy=T) # 0.60

stockmod.l1 <- glm(LTN ~ year + lag.BLPCT + lag.ASNPCT + lag.LTNPCT + lag.MHSPCT + lag.HHINC, data=data.long, family = poisson)

stockmod.l2 <- glm(LTN ~ year + lag.BLPCT + lag.ASNPCT + lag.LTNPCT + lag.MHSPCT + lag.HHINC + lag.MULTIPCT + lag.VACANTPCT + B20YR.PCT, data=data.long, family = poisson)

stockmod.l3 <- glm(LTN ~ year*lag.BLPCT + lag.ASNPCT*year + lag.LTNPCT*year + lag.MHSPCT*year + lag.HHINC*year + lag.MULTIPCT*year + lag.VACANTPCT*year + B20YR.PCT*year + lag.MEDVAL*year, data=data.long, family = poisson)

summary(reg3.l)
exp(reg3.l[[1]])

#split up into 2 years to get Moran's I

stockmod.l3.t1 <- glm(LTN ~ lag.BLPCT + lag.ASNPCT + lag.LTNPCT + lag.MHSPCT + lag.HHINC + lag.MULTIPCT + lag.VACANTPCT + B20YR.PCT, data=data1, family = poisson)

stockmod.l3.t2 <- glm(LTN ~ lag.BLPCT + lag.ASNPCT + lag.LTNPCT + lag.MHSPCT + lag.HHINC + lag.MULTIPCT + lag.VACANTPCT + B20YR.PCT, data=data2, family = poisson)

# Best Model Fit (AIC)
reg1.l[[11]] # 581,678.3
reg2.l[[11]] # 519,416.4
reg3.l[[11]] # 452,323.7

#split up by year
reg3.l.1[[11]] # 179,833.1
reg3.l.2[[11]] # 272,490.5

### Moran's I on residuals per each year
# 2000
yhat1 <- stockmod.l3.t1$fitted.values
y1 <- data1$LTN
plot(yhat1, y1)
lines(lowess(y1 ~ yhat1), lwd=2, col=4)
abline(lm(y1 ~ yhat1), lwd=2, col=2)
# 2010
yhat2 <- stockmod.l3.t2$fitted.values
y2 <- data2$LTN
plot(yhat2, y2)
lines(lowess(y2 ~ yhat2), lwd=2, col=4)
abline(lm(y2 ~ yhat2), lwd=2, col=2)

# Run the Moran test on the residuals:
moran.test(stockmod.l3.t1$residuals, nb.Q1, alternative="two.sided") # 0.382
moran.test(stockmod.l3.t2$residuals, nb.Q1, alternative="two.sided") # 0.389

########################
###
### FLOW MODELS
###
########################


########################
# white population flow
########################

# Global Moran's I
moran.test(wide.sp$WH_2000, nb.Q1, zero.policy=T) # 0.62
moran.test(wide.sp$WH_2010, nb.Q1, zero.policy=T) # 0.64

flowmod.w1 <- glm(WH ~ year + lag.BLPCT + lag.ASNPCT + lag.LTNPCT + lag.MHSPCT + lag.HHINC, offset = log(lag.WH+1), data=data.long, family = poisson)

flowmod.w2 <- glm(WH ~ year + lag.BLPCT + lag.ASNPCT + lag.LTNPCT + lag.MHSPCT + lag.HHINC + lag.MULTIPCT + lag.VACANTPCT + B20YR.PCT, offset = log(lag.WH+1), data=data.long, family = poisson)

flowmod.w3 <- glm(WH ~ year*lag.BLPCT + lag.ASNPCT*year + lag.LTNPCT*year + lag.MHSPCT*year + lag.HHINC*year + lag.MULTIPCT*year + lag.VACANTPCT*year + B20YR.PCT*year + lag.MEDVAL*year, offset = log(lag.WH+1), data=data.long, family = poisson)

summary(flowreg3.w)
exp(flowreg3.w[[1]])

#split up into 2 years to get Moran's I

flowmod.w3.t1 <- glm(WH ~ lag.BLPCT + lag.ASNPCT + lag.LTNPCT + lag.MHSPCT + lag.HHINC + lag.MULTIPCT + lag.VACANTPCT + B20YR.PCT, offset = log(lag.WH+1), data=data1, family = poisson)

flowmod.w3.t2 <- glm(WH ~ lag.BLPCT + lag.ASNPCT + lag.LTNPCT + lag.MHSPCT + lag.HHINC + lag.MULTIPCT + lag.VACANTPCT + B20YR.PCT, data=data2, offset = log(lag.WH+1), family = poisson)

# Best Model Fit (AIC)
flowreg1.w[[11]] # 747237.5
flowreg2.w[[11]] # 593401.5
flowreg3.w[[11]] # 568253.2

#split up by year
flowreg3.w.1[[11]] # 334408.4
flowreg3.w.2[[11]] # 233848.2

### Moran's I on residuals per each year
# 2000
yhat1 <- flowmod.w3.t1$fitted.values
y1 <- data1$WH
plot(yhat1, y1)
lines(lowess(y1 ~ yhat1), lwd=2, col=4)
abline(lm(y1 ~ yhat1), lwd=2, col=2)
# 2010
yhat2 <- flowmod.w3.t2$fitted.values
y2 <- data2$WH
plot(yhat2, y2)
lines(lowess(y2 ~ yhat2), lwd=2, col=4)
abline(lm(y2 ~ yhat2), lwd=2, col=2)

# Run the Moran test on the residuals:
moran.test(flowmod.w3.t1$residuals, nb.Q1, alternative="two.sided") # 0.388
moran.test(flowmod.w3.t2$residuals, nb.Q1, alternative="two.sided") # 0.248


########################
# Black population flow
########################

#Global Moran's I
moran.test(wide.sp$BL_2000, nb.Q1, zero.policy=T) # 0.77
moran.test(wide.sp$BL_2010, nb.Q1, zero.policy=T) # 0.63

flowmod.b1 <- glm(BL ~ year + lag.BLPCT + lag.ASNPCT + lag.LTNPCT + lag.MHSPCT + lag.HHINC, offset = log(lag.BL+1), data=data.long, family = poisson)

flowmod.b2 <- glm(BL ~ year + lag.BLPCT + lag.ASNPCT + lag.LTNPCT + lag.MHSPCT + lag.HHINC + lag.MULTIPCT + lag.VACANTPCT + B20YR.PCT, offset = log(lag.BL+1), data=data.long, family = poisson)

flowmod.b3 <- glm(BL ~ year*lag.BLPCT + lag.ASNPCT*year + lag.LTNPCT*year + lag.MHSPCT*year + lag.HHINC*year + lag.MULTIPCT*year + lag.VACANTPCT*year + B20YR.PCT*year + lag.MEDVAL*year, offset = log(lag.BL+1), data=data.long, family = poisson)

summary(flowreg3.b)
exp(flowreg3.b[[1]])

#split up into 2 years to get Moran's I

flowmod.b3.t1 <- glm(BL ~ lag.BLPCT + lag.ASNPCT + lag.LTNPCT + lag.MHSPCT + lag.HHINC + lag.MULTIPCT + lag.VACANTPCT + B20YR.PCT, offset = log(lag.BL+1), data=data1, family = poisson)

flowmod.b3.t2 <- glm(BL ~ lag.BLPCT + lag.ASNPCT + lag.LTNPCT + lag.MHSPCT + lag.HHINC + lag.MULTIPCT + lag.VACANTPCT + B20YR.PCT, data=data2, offset = log(lag.BL+1), family = poisson)

# Best Model Fit (AIC)
flowreg1.b[[11]] # 788397
flowreg2.b[[11]] # 610637.6
flowreg3.b[[11]] # 575442.7

#split up by year
flowreg3.b.1[[11]] # 285565
flowreg3.b.2[[11]] # 289954.5

### Moran's I on residuals per each year
# 2000
yhat1 <- flowmod.b3.t1$fitted.values
y1 <- data1$BL
plot(yhat1, y1)
lines(lowess(y1 ~ yhat1), lwd=2, col=4)
abline(lm(y1 ~ yhat1), lwd=2, col=2)
# 2010
yhat2 <- flowmod.b3.t2$fitted.values
y2 <- data2$BL
plot(yhat2, y2)
lines(lowess(y2 ~ yhat2), lwd=2, col=4)
abline(lm(y2 ~ yhat2), lwd=2, col=2)

# Run the Moran test on the residuals:
moran.test(flowmod.b3.t1$residuals, nb.Q1, alternative="two.sided") # 0.267
moran.test(flowmod.b3.t2$residuals, nb.Q1, alternative="two.sided") # 0.334

########################
# Asian population flow
########################

# Global Moran's I
moran.test(wide.sp$ASN_2000, nb.Q1, zero.policy=T) # 0.67
moran.test(wide.sp$ASN_2010, nb.Q1, zero.policy=T) # 0.56

flowmod.a1 <- glm(ASN ~ year + lag.BLPCT + lag.ASNPCT + lag.LTNPCT + lag.MHSPCT + lag.HHINC, offset = log(lag.ASN+1), data=data.long, family = poisson)

flowmod.a2 <- glm(ASN ~ year + lag.BLPCT + lag.ASNPCT + lag.LTNPCT + lag.MHSPCT + lag.HHINC + lag.MULTIPCT + lag.VACANTPCT + B20YR.PCT, offset = log(lag.ASN+1), data=data.long, family = poisson)

flowmod.a3 <- glm(ASN ~ year*lag.BLPCT + lag.ASNPCT*year + lag.LTNPCT*year + lag.MHSPCT*year + lag.HHINC*year + lag.MULTIPCT*year + lag.VACANTPCT*year + B20YR.PCT*year + lag.MEDVAL*year, offset = log(lag.ASN+1), data=data.long, family = poisson)

summary(flowreg3.a)
exp(flowreg3.a[[1]])

#split up into 2 years to get Moran's I

flowmod.a3.t1 <- glm(ASN ~ lag.BLPCT + lag.ASNPCT + lag.LTNPCT + lag.MHSPCT + lag.HHINC + lag.MULTIPCT + lag.VACANTPCT + B20YR.PCT, offset = log(lag.ASN+1), data=data1, family = poisson)

flowmod.a3.t2 <- glm(ASN ~ lag.BLPCT + lag.ASNPCT + lag.LTNPCT + lag.MHSPCT + lag.HHINC + lag.MULTIPCT + lag.VACANTPCT + B20YR.PCT, data=data2, offset = log(lag.ASN+1), family = poisson)

# Best Model Fit (AIC)
flowreg1.a[[11]] # 190531.9
flowreg2.a[[11]] # 159126.3
flowreg3.a[[11]] # 156956.7

#split up by year
flowreg3.a.1[[11]] # 64401.39
flowreg3.a.2[[11]] # 92603.85

### Moran's I on residuals per each year
# 2000
yhat1 <- flowmod.a3.t1$fitted.values
y1 <- data1$ASN
plot(yhat1, y1)
lines(lowess(y1 ~ yhat1), lwd=2, col=4)
abline(lm(y1 ~ yhat1), lwd=2, col=2)
# 2010
yhat2 <- flowmod.a3.t2$fitted.values
y2 <- data2$ASN
plot(yhat2, y2)
lines(lowess(y2 ~ yhat2), lwd=2, col=4)
abline(lm(y2 ~ yhat2), lwd=2, col=2)

# Run the Moran test on the residuals:
moran.test(flowmod.a3.t1$residuals, nb.Q1, alternative="two.sided") # 0.061
moran.test(flowmod.a3.t2$residuals, nb.Q1, alternative="two.sided") # 0.222

#########################
# Latinx population flow
#########################

# Global Moran's I
moran.test(wide.sp$LTN_2000, nb.Q1, zero.policy=T) # 0.57
moran.test(wide.sp$LTN_2010, nb.Q1, zero.policy=T) # 0.60

flowmod.l1 <- glm(LTN ~ year + lag.BLPCT + lag.ASNPCT + lag.LTNPCT + lag.MHSPCT + lag.HHINC, offset = log(lag.LTN+1), data=data.long, family = poisson)

flowmod.l2 <- glm(LTN ~ year + lag.BLPCT + lag.ASNPCT + lag.LTNPCT + lag.MHSPCT + lag.HHINC + lag.MULTIPCT + lag.VACANTPCT + B20YR.PCT, offset = log(lag.LTN+1), data=data.long, family = poisson)

flowmod.l3 <- glm(LTN ~ year*lag.BLPCT + lag.ASNPCT*year + lag.LTNPCT*year + lag.MHSPCT*year + lag.HHINC*year + lag.MULTIPCT*year + lag.VACANTPCT*year + B20YR.PCT*year + lag.MEDVAL*year, offset = log(lag.LTN+1), data=data.long, family = poisson)

#split up into 2 years to get Moran's I

flowmod.l3.t1 <- glm(LTN ~ lag.BLPCT + lag.ASNPCT + lag.LTNPCT + lag.MHSPCT + lag.HHINC + lag.MULTIPCT + lag.VACANTPCT + B20YR.PCT, offset = log(lag.LTN+1), data=data1, family = poisson)

flowmod.l3.t2 <- glm(LTN ~ lag.BLPCT + lag.ASNPCT + lag.LTNPCT + lag.MHSPCT + lag.HHINC + lag.MULTIPCT + lag.VACANTPCT + B20YR.PCT, data=data2, offset = log(lag.LTN+1), family = poisson)

# Best Model Fit (AIC)
flowreg1.l[[11]] # 271537.1
flowreg2.l[[11]] # 255221.2
flowreg3.l[[11]] # 229778.7

#split up by year
flowmod.l3.t1[[11]] # 116856.5
flowmod.l3.t2[[11]] # 112863.2

### Moran's I on residuals per each year
# 2000
yhat1 <- flowmod.l3.t1$fitted.values
y1 <- data1$LTN
plot(yhat1, y1)
lines(lowess(y1 ~ yhat1), lwd=2, col=4)
abline(lm(y1 ~ yhat1), lwd=2, col=2)
# 2010
yhat2 <- flowmod.l3.t2$fitted.values
y2 <- data2$LTN
plot(yhat2, y2)
lines(lowess(y2 ~ yhat2), lwd=2, col=4)
abline(lm(y2 ~ yhat2), lwd=2, col=2)

# Run the Moran test on the residuals:
moran.test(flowmod.l3.t1$residuals, nb.Q1, alternative="two.sided") # 0.174
moran.test(flowmod.l3.t2$residuals, nb.Q1, alternative="two.sided") # 0.174

####################################################

# Fixing weights matrix, matching names

####################################################

#data1$wstock1.resid <- residuals(reg1.w)
#residuals.temp <- data1
#W.nb <- poly2nb(residuals.temp, row.names = residuals.temp@data$GEOID10) 
#W <- nb2mat(W.nb, style = "B")
#W.list <- nb2listw(W.nb, style = "B")
#lookup <- data.frame(GEOID10 = residuals.temp@data$GEOID10,
#                     spatialorder = 1:nrow(residuals.temp@data))
#data.temp <- merge(x=data1, y=lookup, by="GEOID10" )
#dat.ordered <- arrange(data.temp@data, spatialorder)

#spatmod1.w.s <- S.CARbym(formula = WH ~ lag.MHSPCT, family="poisson", data = dat.ordered, W = W, burnin=200000, n.sample=2200000, thin = 1000)
#print(spatmod1.w.s) #median and interquartile range

##############################################################
###
### Run CAR models
###
##############################################################

####### Formulas #######

## got rid of MedVal

# stock
stockform.w <- paste("WH ~ lag.BLPCT + lag.ASNPCT + lag.LTNPCT + lag.MHSPCT + lag.HHINC + lag.MULTIPCT + lag.VACANTPCT + B20YR.PCT")
stockform.b <- paste("BL ~ lag.BLPCT + lag.ASNPCT + lag.LTNPCT + lag.MHSPCT + lag.HHINC + lag.MULTIPCT + lag.VACANTPCT + B20YR.PCT")
stockform.a <- paste("ASN ~ lag.BLPCT + lag.ASNPCT + lag.LTNPCT + lag.MHSPCT + lag.HHINC + lag.MULTIPCT + lag.VACANTPCT + B20YR.PCT")
stockform.l <- paste("LTN ~ lag.BLPCT + lag.ASNPCT + lag.LTNPCT + lag.MHSPCT + lag.HHINC + lag.MULTIPCT + lag.VACANTPCT + B20YR.PCT")

# flow
flowform.w <- paste("WH ~ offset(log(lag.WH+1)) + lag.BLPCT + lag.ASNPCT + lag.LTNPCT + lag.MHSPCT + lag.HHINC + lag.MULTIPCT + lag.VACANTPCT + B20YR.PCT + lag.MEDVAL")
flowform.b <- paste("BL ~ offset(log(lag.BL+1)) + lag.BLPCT + lag.ASNPCT + lag.LTNPCT + lag.MHSPCT + lag.HHINC + lag.MULTIPCT + lag.VACANTPCT + B20YR.PCT + lag.MEDVAL")
flowform.a <- paste("ASN ~ offset(log(lag.ASN+1)) + lag.BLPCT + lag.ASNPCT + lag.LTNPCT + lag.MHSPCT + lag.HHINC + lag.MULTIPCT + lag.VACANTPCT + B20YR.PCT + lag.MEDVAL")
flowform.l <- paste("LTN ~ offset(log(lag.LTN+1)) + lag.BLPCT + lag.ASNPCT + lag.LTNPCT + lag.MHSPCT + lag.HHINC + lag.MULTIPCT + lag.VACANTPCT + B20YR.PCT + lag.MEDVAL")


####### BYM MODELS #######

######################
### stock white models
######################

set.seed(838)
CARstock.w1 <- S.CARbym(WH ~ lag.BLPCT + lag.ASNPCT + lag.LTNPCT + lag.MHSPCT + lag.HHINC + lag.MULTIPCT + lag.VACANTPCT + B20YR.PCT, family="poisson", data = data1, W = W.qmat, burnin = 400000, n.sample = 4400000, thin = 2000)
print(CARstock.w1) #median and interquartile range
plot(CARstock.w1$samples$sigma2) #iterate more if it doesn't flatten out

resid.CARstock.w1 <- residuals(CARstock.w1, type="response")
moran.mc(x=resid.CARstock.w1, nb.Q1, nsim=999) # Moran's I: -0.119

set.seed(838)
CARstock.w2 <- S.CARbym(WH ~ lag.BLPCT + lag.ASNPCT + lag.LTNPCT + lag.MHSPCT + lag.HHINC + lag.MULTIPCT + lag.VACANTPCT + B20YR.PCT, family="poisson", data = data2, W = W.qmat, burnin = 400000, n.sample = 4400000, thin = 2000)
print(CARstock.w2) #median and interquartile range
plot(CARstock.w2$samples$sigma2) #iterate more if it doesn't flatten out

resid.CARstock.w2 <- residuals(CARstock.w2, type="response")
moran.mc(x=resid.CARstock.w2, nb.Q1, nsim=999) # Moran's I: -0.100

######################
### stock Black models
######################

set.seed(838)
CARstock.b1 <- S.CARbym(BL ~ lag.BLPCT + lag.ASNPCT + lag.LTNPCT + lag.MHSPCT + lag.HHINC + lag.MULTIPCT + lag.VACANTPCT + B20YR.PCT, family="poisson", data = data1, W = W.qmat, burnin = 600000, n.sample = 6600000, thin = 3000)
print(CARstock.b1) #median and interquartile range
plot(CARstock.b1$samples$sigma2) #iterate more if it doesn't flatten out

resid.CARstock.b1 <- residuals(CARstock.b1, type="response")
moran.mc(x=resid.CARstock.b1, nb.Q1, nsim=999) # Moran's I: -0.119

set.seed(838)
CARstock.b2 <- S.CARbym(BL ~ lag.BLPCT + lag.ASNPCT + lag.LTNPCT + lag.MHSPCT + lag.HHINC + lag.MULTIPCT + lag.VACANTPCT + B20YR.PCT, family="poisson", data = data2, W = W.qmat, burnin = 800000, n.sample = 8800000, thin = 4000)
print(CARstock.b2) #median and interquartile range
plot(CARstock.b2$samples$sigma2) #iterate more if it doesn't flatten out

resid.CARstock.b2 <- residuals(CARstock.b2, type="response")
moran.mc(x=resid.CARstock.b2, nb.Q1, nsim=999) # Moran's I: -0.100

######################
### stock Asian models
######################

set.seed(838)
CARstock.a1 <- S.CARbym(ASN ~ lag.BLPCT + lag.ASNPCT + lag.LTNPCT + lag.MHSPCT + lag.HHINC + lag.MULTIPCT + lag.VACANTPCT + B20YR.PCT, family="poisson", data = data1, W = W.qmat, burnin = 400000, n.sample = 4400000, thin = 2000)
print(CARstock.a1) #median and interquartile range
plot(CARstock.a1$samples$sigma2) #iterate more if it doesn't flatten out

resid.CARstock.a1 <- residuals(CARstock.a1, type="response")
moran.mc(x=resid.CARstock.a1, nb.Q1, nsim=999) # Moran's I: -0.119

set.seed(838)
CARstock.a2 <- S.CARbym(ASN ~ lag.BLPCT + lag.ASNPCT + lag.LTNPCT + lag.MHSPCT + lag.HHINC + lag.MULTIPCT + lag.VACANTPCT + B20YR.PCT, family="poisson", data = data2, W = W.qmat, burnin = 400000, n.sample = 4400000, thin = 2000)
print(CARstock.a2) #median and interquartile range
plot(CARstock.a2$samples$sigma2) #iterate more if it doesn't flatten out

resid.CARstock.a2 <- residuals(CARstock.a2, type="response")
moran.mc(x=resid.CARstock.a2, nb.Q1, nsim=999) # Moran's I: -0.100

######################
### stock Latinx models
######################

set.seed(838)
CARstock.l1 <- S.CARbym(LTN ~ lag.BLPCT + lag.ASNPCT + lag.LTNPCT + lag.MHSPCT + lag.HHINC + lag.MULTIPCT + lag.VACANTPCT + B20YR.PCT, family="poisson", data = data1, W = W.qmat, burnin = 400000, n.sample = 4400000, thin = 2000)
print(CARstock.l1) #median and interquartile range
plot(CARstock.l1$samples$sigma2) #iterate more if it doesn't flatten out

resid.CARstock.l1 <- residuals(CARstock.l1, type="response")
moran.mc(x=resid.CARstock.l1, nb.Q1, nsim=999) # Moran's I: -0.119

set.seed(838)
CARstock.l2 <- S.CARbym(LTN ~ lag.BLPCT + lag.ASNPCT + lag.LTNPCT + lag.MHSPCT + lag.HHINC + lag.MULTIPCT + lag.VACANTPCT + B20YR.PCT, family="poisson", data = data2, W = W.qmat, burnin = 400000, n.sample = 4400000, thin = 2000)
print(CARstock.l2) #median and interquartile range
plot(CARstock.l2$samples$sigma2) #iterate more if it doesn't flatten out

resid.CARstock.l2 <- residuals(CARstock.l2, type="response")
moran.mc(x=resid.CARstock.l2, nb.Q1, nsim=999) # Moran's I: -0.100

######################
### flow white models
######################

set.seed(838)
CARflow.w1 <- S.CARbym(WH ~ offset(log(lag.WH+1)) + lag.BLPCT + lag.ASNPCT + lag.LTNPCT + lag.MHSPCT + lag.HHINC + lag.MULTIPCT + lag.VACANTPCT + B20YR.PCT, family="poisson", data = data1, W = W.qmat, burnin = 400000, n.sample = 4400000, thin = 2000)
print(CARflow.w1) #median and interquartile range
plot(CARflow.w1$samples$sigma2) #iterate more if it doesn't flatten out

resid.CARflow.w1 <- residuals(CARflow.w1, type="response")
moran.mc(x=resid.CARflow.w1, nb.Q1, nsim=999) # Moran's I: -0.119

set.seed(838)
CARflow.w2 <- S.CARbym(WH ~ offset(log(lag.WH+1)) + lag.BLPCT + lag.ASNPCT + lag.LTNPCT + lag.MHSPCT + lag.HHINC + lag.MULTIPCT + lag.VACANTPCT + B20YR.PCT, family="poisson", data = data2, W = W.qmat, burnin = 400000, n.sample = 4400000, thin = 2000)
print(CARflow.w2) #median and interquartile range
plot(CARflow.w2$samples$sigma2) #iterate more if it doesn't flatten out

resid.CARflow.w2 <- residuals(CARflow.w2, type="response")
moran.mc(x=resid.CARflow.w2, nb.Q1, nsim=999) # Moran's I: -0.100

######################
### flow Black models
######################

set.seed(838)
CARflow.b1 <- S.CARbym(BL ~ offset(log(lag.BL+1)) + lag.BLPCT + lag.ASNPCT + lag.LTNPCT + lag.MHSPCT + lag.HHINC + lag.MULTIPCT + lag.VACANTPCT + B20YR.PCT, family="poisson", data = data1, W = W.qmat, burnin = 600000, n.sample = 6600000, thin = 3000)
print(CARflow.b1) #median and interquartile range
plot(CARflow.b1$samples$sigma2) #iterate more if it doesn't flatten out

resid.CARflow.b1 <- residuals(CARflow.b1, type="response")
moran.mc(x=resid.CARflow.b1, nb.Q1, nsim=999) # Moran's I: -0.119

set.seed(838)
CARflow.b2 <- S.CARbym(BL ~ offset(log(lag.BL+1)) + lag.BLPCT + lag.ASNPCT + lag.LTNPCT + lag.MHSPCT + lag.HHINC + lag.MULTIPCT + lag.VACANTPCT + B20YR.PCT, family="poisson", data = data2, W = W.qmat, burnin = 600000, n.sample = 6600000, thin = 3000)
print(CARflow.b2) #median and interquartile range
plot(CARflow.b2$samples$sigma2) #iterate more if it doesn't flatten out

resid.CARflow.b2 <- residuals(CARflow.b2, type="response")
moran.mc(x=resid.CARflow.b2, nb.Q1, nsim=999) # Moran's I: -0.100

######################
### flow Asian models
######################

set.seed(838)
CARflow.a1 <- S.CARbym(ASN ~ offset(log(lag.ASN+1)) + lag.BLPCT + lag.ASNPCT + lag.LTNPCT + lag.MHSPCT + lag.HHINC + lag.MULTIPCT + lag.VACANTPCT + B20YR.PCT, family="poisson", data = data1, W = W.qmat, burnin = 400000, n.sample = 4400000, thin = 2000)
print(CARflow.a1) #median and interquartile range
plot(CARflow.a1$samples$sigma2) #iterate more if it doesn't flatten out

resid.CARflow.a1 <- residuals(CARflow.a1, type="response")
moran.mc(x=resid.CARflow.a1, nb.Q1, nsim=999) # Moran's I: -0.119

set.seed(838)
CARflow.a2 <- S.CARbym(ASN ~ offset(log(lag.ASN+1)) + lag.BLPCT + lag.ASNPCT + lag.LTNPCT + lag.MHSPCT + lag.HHINC + lag.MULTIPCT + lag.VACANTPCT + B20YR.PCT, family="poisson", data = data2, W = W.qmat, burnin = 400000, n.sample = 4400000, thin = 2000)
print(CARflow.a2) #median and interquartile range
plot(CARflow.a2$samples$sigma2) #iterate more if it doesn't flatten out

resid.CARflow.a2 <- residuals(CARflow.a2, type="response")
moran.mc(x=resid.CARflow.a2, nb.Q1, nsim=999) # Moran's I: -0.100

######################
### flow Latinx models
######################

set.seed(838)
CARflow.l1 <- S.CARbym(LTN ~ offset(log(lag.LTN+1)) + lag.BLPCT + lag.ASNPCT + lag.LTNPCT + lag.MHSPCT + lag.HHINC + lag.MULTIPCT + lag.VACANTPCT + B20YR.PCT, family="poisson", data = data1, W = W.qmat, burnin = 400000, n.sample = 4400000, thin = 2000)
print(CARflow.l1) #median and interquartile range
plot(CARflow.l1$samples$sigma2) #iterate more if it doesn't flatten out

resid.CARflow.l1 <- residuals(CARflow.l1, type="response")
moran.mc(x=resid.CARflow.l1, nb.Q1, nsim=999) # Moran's I: -0.119

set.seed(838)
CARflow.l2 <- S.CARbym(LTN ~ offset(log(lag.LTN+1)) + lag.BLPCT + lag.ASNPCT + lag.LTNPCT + lag.MHSPCT + lag.HHINC + lag.MULTIPCT + lag.VACANTPCT + B20YR.PCT, family="poisson", data = data2, W = W.qmat, burnin = 400000, n.sample = 4400000, thin = 2000)
print(CARflow.l2) #median and interquartile range
plot(CARflow.l2$samples$sigma2) #iterate more if it doesn't flatten out

resid.CARflow.l2 <- residuals(CARflow.l2, type="response")
moran.mc(x=resid.CARflow.l2, nb.Q1, nsim=999) # Moran's I: -0.100

###########################################
###
### TABLES
###
###########################################

# White Stock
stargazer(stockmod.w1, stockmod.w2, stockmod.w3,
          title = "Poisson Model: white stock", 
          align = T, 
          dep.var.labels="Absolute number of white residents", 
          covariate.labels=c("year", "Percent Black", "Percent Asian", "Percent Latinx", "Percent More than HS edu", "Median HouseHold Income", "Percent Multi-family Housing", "Vacancy Rate", "Percent New Construction", "Median Housing Value", "year x Percent Black", "year x Percent Asian", "year x Percent Latinx", "year x Percent More than HS edu", "year x Median HouseHold Income", "year x Percent Multi-family Housing", "year x Vacancy Rate", "year x Percent New Construction", "year x Median Housing Value"),
          model.numbers = F,
          column.labels = c("model 1", "model 2", "model 3"),
          no.space = T, 
          omit.stat=c("ser","f"))

# Black Stock
stargazer(stockmod.b1, stockmod.b2, stockmod.b3,
          title = "Poisson Model: Black stock", 
          align = T, 
          dep.var.labels="Absolute number of Black residents", 
          covariate.labels=c("year", "Percent Black", "Percent Asian", "Percent Latinx", "Percent More than HS edu", "Median HouseHold Income", "Percent Multi-family Housing", "Vacancy Rate", "Percent New Construction", "Median Housing Value", "year x Percent Black", "year x Percent Asian", "year x Percent Latinx", "year x Percent More than HS edu", "year x Median HouseHold Income", "year x Percent Multi-family Housing", "year x Vacancy Rate", "year x Percent New Construction", "year x Median Housing Value"),
          model.numbers = F,
          column.labels = c("model 1", "model 2", "model 3"),
          no.space = T, 
          omit.stat=c("ser","f"))

# Asian Stock
stargazer(stockmod.a1, stockmod.a2, stockmod.a3,
          title = "Poisson Model: Asian stock", 
          align = T, 
          dep.var.labels="Absolute number of Asian residents", 
          covariate.labels=c("year", "Percent Black", "Percent Asian", "Percent Latinx", "Percent More than HS edu", "Median HouseHold Income",  "Percent Multi-family Housing", "Vacancy Rate", "Percent New Construction", "Median Housing Value", "year x Percent Black", "year x Percent Asian", "year x Percent Latinx", "year x Percent More than HS edu", "year x Median HouseHold Income", "year x Percent Multi-family Housing", "year x Vacancy Rate", "year x Percent New Construction", "year x Median Housing Value"),
          model.numbers = F,
          column.labels = c("model 1", "model 2", "model 3"),
          no.space = T, 
          omit.stat=c("ser","f"))

# Latinx Stock
stargazer(stockmod.l1, stockmod.l2, stockmod.l3,
          title = "Poisson Model: Latinx stock", 
          align = T, 
          dep.var.labels="Absolute number of Latinx residents", 
          covariate.labels=c("year", "Percent Black", "Percent Asian", "Percent Latinx", "Percent More than HS edu", "Median HouseHold Income",  "Percent Multi-family Housing", "Vacancy Rate", "Percent New Construction", "Median Housing Value", "year x Percent Black", "year x Percent Asian", "year x Percent Latinx", "year x Percent More than HS edu", "year x Median HouseHold Income", "year x Percent Multi-family Housing", "year x Vacancy Rate", "year x Percent New Construction", "year x Median Housing Value"),
          model.numbers = F,
          column.labels = c("model 1", "model 2", "model 3"),
          no.space = T, 
          omit.stat=c("ser","f"))

# White Flow
stargazer(flowmod.w1, flowmod.w2, flowmod.w3,
          title = "Poisson Model: white flow", 
          align = T, 
          dep.var.labels="White Growth", 
          covariate.labels=c("year", "Percent Black", "Percent Asian", "Percent Latinx", "Percent More than HS edu", "Median HouseHold Income", "Ownership Rate", "Percent Multi-family Housing", "Vacancy Rate", "Percent New Construction", "Median Housing Value", "year x Percent Black", "year x Percent Asian", "year x Percent Latinx", "year x Percent More than HS edu", "year x Median HouseHold Income", "year x Ownership Rate", "year x Percent Multi-family Housing", "year x Vacancy Rate", "year x Percent New Construction", "year x Median Housing Value"),
          model.numbers = F,
          column.labels = c("model 1", "model 2", "model 3"),
          no.space = T, 
          omit.stat=c("ser","f"))

# Black Flow
stargazer(flowmod.b1, flowmod.b2, flowmod.b3,
          title = "Poisson Model: Black flow", 
          align = T, 
          dep.var.labels="Black Growth", 
          covariate.labels=c("year", "Percent Black", "Percent Asian", "Percent Latinx", "Percent More than HS edu", "Median HouseHold Income", "Ownership Rate", "Percent Multi-family Housing", "Vacancy Rate", "Percent New Construction", "Median Housing Value", "year x Percent Black", "year x Percent Asian", "year x Percent Latinx", "year x Percent More than HS edu", "year x Median HouseHold Income", "year x Ownership Rate", "year x Percent Multi-family Housing", "year x Vacancy Rate", "year x Percent New Construction", "year x Median Housing Value"),
          model.numbers = F,
          column.labels = c("model 1", "model 2", "model 3"),
          no.space = T, 
          omit.stat=c("ser","f"))

# Asian Flow
stargazer(flowmod.a1, flowmod.a2, flowmod.a3,
          title = "Poisson Model: Asian flow", 
          align = T, 
          dep.var.labels="Asian Growth", 
          covariate.labels=c("year", "Percent Black", "Percent Asian", "Percent Latinx", "Percent More than HS edu", "Median HouseHold Income", "Ownership Rate", "Percent Multi-family Housing", "Vacancy Rate", "Percent New Construction", "Median Housing Value", "year x Percent Black", "year x Percent Asian", "year x Percent Latinx", "year x Percent More than HS edu", "year x Median HouseHold Income", "year x Ownership Rate", "year x Percent Multi-family Housing", "year x Vacancy Rate", "year x Percent New Construction", "year x Median Housing Value"),
          model.numbers = F,
          column.labels = c("model 1", "model 2", "model 3"),
          no.space = T, 
          omit.stat=c("ser","f"))

# Latinx Flow
stargazer(flowmod.l1, flowmod.l2, flowmod.l3,
          title = "Poisson Model: Latinx flow", 
          align = T, 
          dep.var.labels="Latinx Growth", 
          covariate.labels=c("year", "Percent Black", "Percent Asian", "Percent Latinx", "Percent More than HS edu", "Median HouseHold Income", "Ownership Rate", "Percent Multi-family Housing", "Vacancy Rate", "Percent New Construction", "Median Housing Value", "year x Percent Black", "year x Percent Asian", "year x Percent Latinx", "year x Percent More than HS edu", "year x Median HouseHold Income", "year x Ownership Rate", "year x Percent Multi-family Housing", "year x Vacancy Rate", "year x Percent New Construction", "year x Median Housing Value"),
          model.numbers = F,
          column.labels = c("model 1", "model 2", "model 3"),
          no.space = T, 
          omit.stat=c("ser","f"))

##### Moran's I Table #####

# Global Moran's I
moran.test(wide.sp$WH_2000, nb.Q1, zero.policy=T) # 0.6235654988 
moran.test(wide.sp$WH_2010, nb.Q1, zero.policy=T) # 0.6417130624
moran.test(wide.sp$BL_2000, nb.Q1, zero.policy=T) # 0.7671965904 
moran.test(wide.sp$BL_2010, nb.Q1, zero.policy=T) # 0.6290144370 
moran.test(wide.sp$ASN_2000, nb.Q1, zero.policy=T) # 0.6703690583
moran.test(wide.sp$ASN_2010, nb.Q1, zero.policy=T) # 0.5571060951 
moran.test(wide.sp$LTN_2000, nb.Q1, zero.policy=T) # 0.5651845957 
moran.test(wide.sp$LTN_2010, nb.Q1, zero.policy=T) # 0.5988483937
globalM <- c("0.62^{***}", "0.64^{***}", "0.77^{***}", "0.63^{***}", "0.67^{***}", "0.56^{***}", "0.57^{***}", "0.60^{***}")

# Moran's Test on Stock Model 2 Residuals
moran.test(stockmod.w3.t1$residuals, nb.Q1, alternative="two.sided") # 0.3680828576 -s #multi: 0.3799546862
moran.test(stockmod.w3.t2$residuals, nb.Q1, alternative="two.sided") # 0.3689001253 -s 
moran.test(stockmod.b3.t1$residuals, nb.Q1, alternative="two.sided") # 0.3694403302 -s
moran.test(stockmod.b3.t2$residuals, nb.Q1, alternative="two.sided") # 0.06665351 -s
moran.test(stockmod.a3.t1$residuals, nb.Q1, alternative="two.sided") # 0.1214328013 -s
moran.test(stockmod.a3.t2$residuals, nb.Q1, alternative="two.sided") # 0.3333329077 -s
moran.test(stockmod.l3.t1$residuals, nb.Q1, alternative="two.sided") # 0.2853849105 -s
moran.test(stockmod.l3.t2$residuals, nb.Q1, alternative="two.sided") # 0.3586348686 -s
stockmod.M <- c("0.37^{***}", "0.37^{***}", "0.37^{***}", "0.07^{***}", "0.12^{***}", "0.33^{***}", "0.29^{***}", "0.36^{***}")

# Moran's Test on Flow Model 2 Residuals
moran.test(flowmod.w3.t1$residuals, nb.Q1, alternative="two.sided") # 0.3890259815 -s 
moran.test(flowmod.w3.t2$residuals, nb.Q1, alternative="two.sided") # 0.2387115489 -s
moran.test(flowmod.b3.t1$residuals, nb.Q1, alternative="two.sided") # 0.2566704498 -s 
moran.test(flowmod.b3.t2$residuals, nb.Q1, alternative="two.sided") # 0.3275759082 -s
moran.test(flowmod.a3.t1$residuals, nb.Q1, alternative="two.sided") # 0.05554205 -s
moran.test(flowmod.a3.t2$residuals, nb.Q1, alternative="two.sided") # 0.2266638311 -s
moran.test(flowmod.l3.t1$residuals, nb.Q1, alternative="two.sided") # 0.1631298307 -s 
moran.test(flowmod.l3.t2$residuals, nb.Q1, alternative="two.sided") # 0.1630738069 -s
flowmod.M <- c("0.39^{***}", "0.24^{***}", "0.26^{***}", "0.33^{***}", "0.05^{***}", "0.23^{***}", "0.16^{***}", "0.16^{***}")

# Moran's Test on Stock Car Model Residuals
moran.mc(x=resid.CARstock.w1, nb.Q1, nsim=999) # -0.12366 - ns
moran.mc(x=resid.CARstock.w2, nb.Q1, nsim=999) # -0.10908 - ns
moran.mc(x=resid.CARstock.b1, nb.Q1, nsim=999) # -0.12137 - ns
moran.mc(x=resid.CARstock.b2, nb.Q1, nsim=999) # -0.13544 - ns
moran.mc(x=resid.CARstock.a1, nb.Q1, nsim=999) # -0.13535 - ns
moran.mc(x=resid.CARstock.a2, nb.Q1, nsim=999) # -0.133 - ns
moran.mc(x=resid.CARstock.l1, nb.Q1, nsim=999) # -0.13743 - ns
moran.mc(x=resid.CARstock.l2, nb.Q1, nsim=999) # -0.12266 - ns
stockcar.M <- c("-0.12", "-0.11", "-0.12", "-0.14", "-0.14", "-0.13", "-0.14", "-0.12")

# Moran's Test on Flow Car Model Residuals
moran.mc(x=resid.CARflow.w1, nb.Q1, nsim=999) # -0.10811 - ns
moran.mc(x=resid.CARflow.w2, nb.Q1, nsim=999) # -0.14624 - ns
moran.mc(x=resid.CARflow.b1, nb.Q1, nsim=999) # -0.11227 - ns
moran.mc(x=resid.CARflow.b2, nb.Q1, nsim=999) # -0.099583 - ns
moran.mc(x=resid.CARflow.a1, nb.Q1, nsim=999) # -0.079507 - ns
moran.mc(x=resid.CARflow.a2, nb.Q1, nsim=999) # -0.092725 - ns
moran.mc(x=resid.CARflow.l1, nb.Q1, nsim=999) # -0.11903 - ns
moran.mc(x=resid.CARflow.l2, nb.Q1, nsim=999) # -0.06961 - ns
flowcar.M <- c("-0.11", "-0.15", "-0.11", "-0.10", "-0.08", "-0.09", "-0.12", "-0.07")

year <- c(2000,2010,2000,2010,2000,2010,2000,2010)
group <- c("White", "White", "Black", "Black", "Asian", "Asian", "Latinx", "Latinx")

moran.tab <- cbind(group, year, globalM, stockmod.M, stockcar.M, flowmod.M,  flowcar.M)
moran.tab <- as.data.table(moran.tab)

moran.tab1 <- moran.tab %>% filter(year == 2000)
moran.tab1 <- moran.tab1[,-2]
colnames(moran.tab1) <- c(" ", "Global Moran", "Poisson Stock", "CAR Stock", "Poisson Flow", "CAR Flow")
moran.tab2 <- moran.tab %>% filter(year == 2010)
moran.tab2 <- moran.tab2[,-2]
colnames(moran.tab2) <- colnames(moran.tab1)

# 2000 Moran table
print(xtable(moran.tab1, caption= "2000 Moran's I Comparison: Global, Poisson Models, CAR models"),
      include.rownames = F,
      booktabs = T,
      caption.placement = "top",
      auto=T,
      sanitize.text.function = function(moran.tab1){moran.tab1})

# 2010 Moran table
print(xtable(moran.tab2, caption= "2010 Moran's I Comparison: Global, Poisson Models, CAR models"),
      include.rownames = F,
      booktabs = T,
      caption.placement = "top",
      auto=T,
      sanitize.text.function = function(moran.tab2){moran.tab2})

# full Moran years
moran.df <- dcast(moran.tab, group~year, value.var=c("globalM", "stockmod.M", "stockcar.M", "flowmod.M", "flowcar.M"))
moran.df <- moran.df[c(4,2,1,3),]
colnames(moran.df) <- c(" ", "Global Moran", "Global Moran 2010", "Poisson Stock: Moran", "Poisson Stock: Moran 2010", "CAR Stock: Moran", "CAR Stock: Moran 2010", "Poisson Flow: Moran", "Poisson Flow: Moran 2010", "CAR Flow: Moran", "CAR Flow: Moran 2010")

print(xtable(moran.df, caption= "Moran's I Comparison: Global, Poisson Models, CAR models"),
      include.rownames = F,
      booktabs = T,
      caption.placement = "top",
      auto=T,
      sanitize.text.function = function(moran.df){moran.df})

######## CAR Model Table ############

library(numform)

# CAR STOCK SEGMENTS
CAR.w.s1 <- as.data.frame(CARstock.w1[[1]])
CAR.w.s1 <- CAR.w.s1[,c(1,2,3)]
CAR.w.s1$CI <- paste("[", f_num(CAR.w.s1$`2.5%`, 3), ", ", f_num(CAR.w.s1$`97.5%`,3), "]", sep="")
tabCARstock.w1 <- CAR.w.s1[,c(1,4)]

CAR.w.s2 <- as.data.frame(CARstock.w2[[1]])
CAR.w.s2 <- CAR.w.s2[,c(1,2,3)]
CAR.w.s2$CI <- paste("[", f_num(CAR.w.s2$`2.5%`,3), ", ", f_num(CAR.w.s2$`97.5%`,3), "]", sep="")
tabCARstock.w2 <- CAR.w.s2[,c(1,4)]

CAR.b.s1 <- as.data.frame(CARstock.b1[[1]])
CAR.b.s1 <- CAR.b.s1[,c(1,2,3)]
CAR.b.s1$CI <- paste("[", f_num(CAR.b.s1$`2.5%`,3), ", ", f_num(CAR.b.s1$`97.5%`,3) ,"]", sep="")
tabCARstock.b1 <- CAR.b.s1[,c(1,4)]

CAR.b.s2 <- as.data.frame(CARstock.b2[[1]])
CAR.b.s2 <- CAR.b.s2[,c(1,2,3)]
CAR.b.s2$CI <- paste("[", f_num(CAR.b.s2$`2.5%`,3), ", ", f_num(CAR.b.s2$`97.5%`,3), "]", sep="")
tabCARstock.b2 <- CAR.b.s2[,c(1,4)]

CAR.a.s1 <- as.data.frame(CARstock.a1[[1]])
CAR.a.s1 <- CAR.a.s1[,c(1,2,3)]
CAR.a.s1$CI <- paste("[", f_num(CAR.a.s1$`2.5%`,3), ", ", f_num(CAR.a.s1$`97.5%`,3), "]", sep="")
tabCARstock.a1 <- CAR.a.s1[,c(1,4)]

CAR.a.s2 <- as.data.frame(CARstock.a2[[1]])
CAR.a.s2 <- CAR.a.s2[,c(1,2,3)]
CAR.a.s2$CI <- paste("[", f_num(CAR.a.s2$`2.5%`,3), ", ", f_num(CAR.a.s2$`97.5%`,3), "]", sep="")
tabCARstock.a2 <- CAR.a.s2[,c(1,4)]

CAR.l.s1 <- as.data.frame(CARstock.l1[[1]])
CAR.l.s1 <- CAR.l.s1[,c(1,2,3)]
CAR.l.s1$CI <- paste("[", f_num(CAR.l.s1$`2.5%`,3), ", ", f_num(CAR.l.s1$`97.5%`,3), "]", sep="")
tabCARstock.l1 <- CAR.l.s1[,c(1,4)]

CAR.l.s2 <- as.data.frame(CARstock.l2[[1]])
CAR.l.s2 <- CAR.l.s2[,c(1,2,3)]
CAR.l.s2$CI <- paste("[", f_num(CAR.l.s2$`2.5%`,3), ", ", f_num(CAR.l.s2$`97.5%`,3), "]", sep="")
tabCARstock.l2 <- CAR.l.s2[,c(1,4)]

# CAR FLOW SEGMENTS
CAR.w.f1 <- as.data.frame(CARflow.w1[[1]])
CAR.w.f1 <- CAR.w.f1[,c(1,2,3)]
CAR.w.f1$CI <- paste("[", f_num(CAR.w.f1$`2.5%`,3), ", ", f_num(CAR.w.f1$`97.5%`,3), "]", sep="")
tabCARflow.w1 <- CAR.w.f1[,c(1,4)]

CAR.w.f2 <- as.data.frame(CARflow.w2[[1]])
CAR.w.f2 <- CAR.w.f2[,c(1,2,3)]
CAR.w.f2$CI <- paste("[", f_num(CAR.w.f2$`2.5%`,3), ", ", f_num(CAR.w.f2$`97.5%`,3), "]", sep="")
tabCARflow.w2 <- CAR.w.f2[,c(1,4)]

CAR.b.f1 <- as.data.frame(CARflow.b1[[1]])
CAR.b.f1 <- CAR.b.f1[,c(1,2,3)]
CAR.b.f1$CI <- paste("[", f_num(CAR.b.f1$`2.5%`,3), ", ", f_num(CAR.b.f1$`97.5%`,3), "]", sep="")
tabCARflow.b1 <- CAR.b.f1[,c(1,4)]

CAR.b.f2 <- as.data.frame(CARflow.b2[[1]])
CAR.b.f2 <- CAR.b.f2[,c(1,2,3)]
CAR.b.f2$CI <- paste("[", f_num(CAR.b.f2$`2.5%`,3), ", ", f_num(CAR.b.f2$`97.5%`,3), "]", sep="")
tabCARflow.b2 <- CAR.b.f2[,c(1,4)]

CAR.a.f1 <- as.data.frame(CARflow.a1[[1]])
CAR.a.f1 <- CAR.a.f1[,c(1,2,3)]
CAR.a.f1$CI <- paste("[", f_num(CAR.a.f1$`2.5%`,3), ", ", f_num(CAR.a.f1$`97.5%`,3), "]", sep="")
tabCARflow.a1 <- CAR.a.f1[,c(1,4)]

CAR.a.f2 <- as.data.frame(CARflow.a2[[1]])
CAR.a.f2 <- CAR.a.f2[,c(1,2,3)]
CAR.a.f2$CI <- paste("[", f_num(CAR.a.f2$`2.5%`,3), ", ", f_num(CAR.a.f2$`97.5%`,3), "]", sep="")
tabCARflow.a2 <- CAR.a.f2[,c(1,4)]

CAR.l.f1 <- as.data.frame(CARflow.l1[[1]])
CAR.l.f1 <- CAR.l.f1[,c(1,2,3)]
CAR.l.f1$CI <- paste("[", f_num(CAR.l.f1$`2.5%`,3), ", ", f_num(CAR.l.f1$`97.5%`,3), "]", sep="")
tabCARflow.l1 <- CAR.l.f1[,c(1,4)]

CAR.l.f2 <- as.data.frame(CARflow.l2[[1]])
CAR.l.f2 <- CAR.l.f2[,c(1,2,3)]
CAR.l.f2$CI <- paste("[", f_num(CAR.l.f2$`2.5%`,3), ", ", f_num(CAR.l.f2$`97.5%`,3), "]", sep="")
tabCARflow.l2 <- CAR.l.f2[,c(1,4)]

# CAR TABLES
CAR.mat.w <- cbind(tabCARstock.w1, tabCARstock.w2, tabCARflow.w1, tabCARflow.w2) # 2000, no group names
CAR.mat.b <- cbind(tabCARstock.b1, tabCARstock.b2, tabCARflow.b1, tabCARflow.b2)
CAR.mat.a <- cbind(tabCARstock.a1, tabCARstock.a2, tabCARflow.a1, tabCARflow.a2)
CAR.mat.l <- cbind(tabCARstock.l1, tabCARstock.l2, tabCARflow.l1, tabCARflow.l2)

rownames(CAR.mat.w) <- c("intercept", "% Black", "% Asian", "% Latinx", "% More than HS", "Med. HHINC", "Multi-Family Rate", "Vacancy Rate", "% New Construction","Tau", "Sigma")
rownames(CAR.mat.b) <- c("intercept", "% Black", "% Asian", "% Latinx", "% More than HS", "Med. HHINC", "Multi-Family Rate", "Vacancy Rate", "% New Construction","Tau", "Sigma")
rownames(CAR.mat.a) <- c("intercept", "% Black", "% Asian", "% Latinx", "% More than HS", "Med. HHINC", "Multi-Family Rate", "Vacancy Rate", "% New Construction","Tau", "Sigma")
rownames(CAR.mat.l) <- c("intercept", "% Black", "% Asian", "% Latinx", "% More than HS", "Med. HHINC", "Multi-Family Rate", "Vacancy Rate", "% New Construction","Tau", "Sigma")

addtorow <- list()
addtorow$pos <- list(0,0,0)
addtorow$command <- c("& \\multicolumn{2}{c}{Stock Model: 2000} & \\multicolumn{2}{c}{Stock Model: 2010} & \\multicolumn{2}{c}{Flow Model: 2000} & \\multicolumn{2}{c}{Flow Model: 2010} \\\\\n", " \\cline{2-3} \\cline{4-5} \\cline{6-7} \\cline{8-9} \\\\\n", "& \\multicolumn{1}{c}{Coeff.} & \\multicolumn{1}{c}{95\\% Int} & \\multicolumn{1}{c}{Coeff.} & \\multicolumn{1}{c}{95\\% Int} & \\multicolumn{1}{c}{Coeff.} & \\multicolumn{1}{c}{95\\% Int} & \\multicolumn{1}{c}{Coeff.} & \\multicolumn{1}{c}{95\\% Int} \\\\\n ")

print(xtable(CAR.mat.w, 
             digits = 3,
             caption = "CAR Models for White population in 2000 and 2010"),
      add.to.row = addtorow,
      include.rownames = T,
      include.colnames = F,
      booktabs = T,
      caption.placement = "top",
      auto=T)

print(xtable(CAR.mat.b, 
             digits = 3,
             caption = "CAR Models for Black population in 2000 and 2010"),
      add.to.row = addtorow,
      include.rownames = T,
      include.colnames = F,
      booktabs = T,
      caption.placement = "top",
      auto=T)

print(xtable(CAR.mat.a, 
             digits = 3,
             caption = "CAR Models for Asian population in 2000 and 2010"),
      add.to.row = addtorow,
      include.rownames = T,
      include.colnames = F,
      booktabs = T,
      caption.placement = "top",
      auto=T)

print(xtable(CAR.mat.l, 
             digits = 3,
             caption = "CAR Models for Latinx population in 2000 and 2010"),
      add.to.row = addtorow,
      include.rownames = T,
      include.colnames = F,
      booktabs = T,
      caption.placement = "top",
      auto=T)

