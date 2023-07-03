#----------------------------------------------------------------------------------
#script name:  Distance_distribution_dynamics.R
#written by :  Original script by Barbara (2015) and Ivar Olaf Peereboom
#              restructured, revised and further developed by A.Stavang (2019).
#              Developed for urban hydrology by Thomas Skaugen 2022
#
# Model input:    1) DEM: Digital elevation model of catchment
#                 2) catchment boundary: .shp-file
#                 3) River network shapefile
#                 4) Landuse shapefie. buildings -impermeable areas
#                 5) Landuse shapefile roads  -impermeable areas
#
# Model output:   1) Hypsographic curve
#                 2) Distance distribution statistics for Permeable areas (Soils) 
#                 3) Distance distribution statistics for Impermeable areas (Rooftops and roads)
#                 4) Distance distribution statistics for river network. 
#                 5) Parameterfile for DDDUrban
# 
#----------------------------------------------------------------------------------
rm(list=ls()) # removes all from global environment

#Packages
library(raster)
library(rgdal)
library(rgeos)
library(igraph)
library(maptools)
library (DMMF) 

# Input files
station <- "123.38"
res <- 500 # CA = 500 m2

# UTM 33 coordinates for outlet of catchment 
xoutlet <- 271457# coordinates of the outlet
youtlet <- 7038169# coordinates of the outlet

# lat long coordinates for ca centroid of catchment (for Energy Balance calculations) 
lat <-63.39615
long <- 10.42417

DEM <- raster("\\\\nve.no\\fil\\h\\HV\\Bruker\\aktv\\DDDUrban\\Risvollan\\Risvollan_burned.tif")  # 1 m. UTM zone N33. Import DEM (digital elevation map) file. 
catchment <- readOGR(dsn = path.expand("\\\\nve.no\\fil\\h\\HV\\Bruker\\aktv\\DDDUrban\\Risvollan\\Grunnlagsdata\\ths\\Risvollan_basin500.shp"))
landuse1 <-  readOGR(dsn = path.expand("\\\\nve.no\\fil\\h\\HV\\Bruker\\aktv\\DDDUrban\\Risvollan\\Grunnlagsdata\\ths\\Risvollan_Bygning.shp")) #hustak
landuse2 <-  readOGR(dsn = path.expand("\\\\nve.no\\fil\\h\\HV\\Bruker\\aktv\\DDDUrban\\Risvollan\\Grunnlagsdata\\ths\\Risvollan_AR5.shp")) #veier
river.shapefile <- readOGR(dsn = path.expand("\\\\nve.no\\fil\\h\\HV\\Bruker\\aktv\\DDDUrban\\Risvollan\\Grunnlagsdata\\ths\\Risvollan_stream500.shp")) #

plotfile <- paste("\\\\nve.no\\fil\\h\\HB\\HB-modellering\\DDDtestbenk\\DDDUrban\\DDplots\\",station,"_",res,".png",sep="")
#---------------------------------------------------------------------------------#
#  (A)                            HYPSOGRAPHIC CURVE                              #
#---------------------------------------------------------------------------------#

# 1. Adjust coordinate systems. The two files need to be on the same coordinate system. Use DEM's crs as standard. 
#    To check: use crs- and/or extent-command on both files. 

# 2. Masking the DEM with the catchment boundary. 
#    Crop a raster (DEM) to vector extent of the catchment using the mask-function. 
#    Catchment boundary contains all Norwegian HBV-felt. 

catchment_boundary <-  catchment 
DEM_masked <- mask(DEM, catchment_boundary) #Digital elevation map for "Felt". 

# 3. Calculate hypsographic curve, mean elevation and area (m2)

#   a) calculate values from Hypsographic curve
#   a00  	Hypsografic curve  minimum value
#   a01		Hypsografic curve  elev first areal 10% 
#   a02		Hypsografic curve  elev second areal 10%
#   a03		Hypsografic curve  elev 3 areal 10%
#   a04		Hypsografic curve  elev 4 areal 10%
#   a05		Hypsografic curve  elev 5 areal 10%
#   a06		Hypsografic curve  elev 6 areal 10%
#   a07		Hypsografic curve  elev 7 areal 10%
#   a08		Hypsografic curve  elev 8 areal 10%
#   a09	  Hypsografic curve  elev 9 areal 10%
#   a10		Hypsografic curve  maximum value

a00  <- round(minValue(DEM_masked))
a10  <- round(maxValue(DEM_masked))
hypso <- quantile(DEM_masked,probs = c(0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1),names=F)
a01<- round(hypso[2]);
a02<- round(hypso[3]);
a03<- round(hypso[4]);
a04<- round(hypso[5]);
a05<- round(hypso[6]);
a06<- round(hypso[7]);
a07<- round(hypso[8]);
a08<- round(hypso[9]);
a09<- round(hypso[10]);

hfelt <- round(cellStats(DEM_masked,mean))

#   c) area ## Area in m2
options(scipen=999) # disable scientific notation
cells.data  <- length(which(is.na(values(DEM_masked))== F)) # cells with DEM data
area <- cells.data * xres(DEM_masked) * yres(DEM_masked)
print(paste("area is[m2]:",area ))
#---------------------------------------------------------------------------------#
#  (B)               DISTANCE DISTRIBUTION AREA TYPE                              #
#---------------------------------------------------------------------------------#
#     Need to calculate the distance distribution from each catchment characteristic (CC) to point in river. 

# Make catchment raster
#   Set ut a raster template. The pixel size, extent of raster layer will follow this template.

r <- raster(extent(DEM_masked), ncol=ncol(DEM_masked), nrow=nrow(DEM_masked)) #DEM_masked has 10x10 grid #original fra Anne
extent_felt <-rasterize(catchment_boundary, r)#           #rastersize/grid the catchment with cell size = r. 


# Make raster of rivernet
#   Make sure files are same coordinate system
river.shapefile <- spTransform(river.shapefile,crs(DEM)) 
grd_rivernet_felt <- mask(extent_felt, river.shapefile)

grd_rivernet_felt <- reclassify(grd_rivernet_felt, matrix(c(1,2,1), byrow=TRUE)) #Makes all values = river = 1 # de to første i c(1,2,1) skal være intervallet verdiene er i 

# Calculate the euclidian distance from river raster to points in the catchment. 
grd_distance <- distance(grd_rivernet_felt, method = "euclidian") #Calculates the distance from river network to points in catchment

# Mask the DD with the shape of catchment
grd_distance_masked <- mask(grd_distance, catchment_boundary) #This is the DD for point in catchment to river.

# Make the Distance rasters for area type
# Make sure the crs is the same. 
landuse1_read <- spTransform(landuse1, crs(DEM))        # Bygninger sjekk attributt tabell i QGIS for å se hva som finnes...
hustak <- landuse1_read[landuse1_read$omradeid=="5001",]#
landuse2_read <- spTransform(landuse2, crs(DEM))        ## 12 er samferdsel, check the "attributt" table in QGIS
veg <-landuse2_read[landuse2_read$arealtype=="12",]


Soilgrd_distance1 <- mask(grd_distance_masked, hustak)
tullhustak <-length(which(is.na(values(Soilgrd_distance1))== F)) # number of cells with rooftop data

Soilgrd_distance2<- mask(grd_distance_masked, veg)
tullveg <-length(which(is.na(values(Soilgrd_distance2))== F)) # number of cells with roads data

IPgrd_distance <- merge(Soilgrd_distance1,Soilgrd_distance2)
plot(IPgrd_distance)

husveg <- which(values(IPgrd_distance)>0.0) #position for roads and roofs
Pgrd_distance <-grd_distance_masked
Pgrd_distance[husveg] <- NA
plot(Pgrd_distance)


# Parameters to paramfile
# 1. Areal fraction (of total area)
# 2. Distance distribution: max distance, mean distance
# 3. zero distance areal fraction.

#--- 1.AREA FRACTION
# Permeable areas
cells_data_P  <- length(which(values(Pgrd_distance) >0.0)) # count of cells = SOILS dvs  Permeable fraction
areaP <- cells_data_P * xres(Soilgrd_distance1) * yres(Soilgrd_distance1) # area of SOILS
cells_data_IP <- length(husveg)
areaIP <- cells_data_IP * xres(Soilgrd_distance1) * yres(Soilgrd_distance1) # area of SOILS
Pfrac <- areaP/area # Fraction of soils of total area m2
IPfrac <- 1-Pfrac  #areaIP/area # Fraction of soils of total area m2. Needs to add up to 1 anyway
print(Pfrac)
print(IPfrac)


#--- 2. MAX, MIN and MEAN 
# Permeable area
Pmax <-0.0
Pmin <-0.0
Pmid <-0.0
print(Pmax <- Pgrd_distance@data@max) #saves the maximum distance
print(Pmin <- Pgrd_distance@data@min) #saves the minimum distance
print(Pmid <- cellStats(Pgrd_distance, stat = 'mean', na.rm = TRUE))

# Impermeable area
IPmax <-0.0
IPmin <-0.0
IPmid <-0.0
print(IPmax <- IPgrd_distance@data@max) #saves the maximum distance
print(IPmin <- IPgrd_distance@data@min) #saves the minimum distance
print(IPmid <- cellStats(IPgrd_distance, stat = 'mean', na.rm = TRUE))

freqDD <-freq(Pgrd_distance, useNA ="no")
#must make akkumulated area
akkarea <- vector("numeric", length(freqDD[,2]))
frac <- vector("numeric", length(freqDD[,2]))
for(i in 1 :length(freqDD[,2]))
{
  akkarea[i] <- sum(freqDD[,2][1:i])/sum(freqDD[,2])
  frac[i] <- (1-exp(-((freqDD[,1][i])/(Pmid))))
}
windows(18,16)
par(mfrow=c(1,1))

plot(freqDD[,1],akkarea,type = "l", xlab="Distance [m]", ylab="Fraction", main= paste("Distance distribution ",station,"_CA",res," m2",sep=""))
lines(freqDD[,1],frac, col="red")
text (Pmid,0.1,paste("Pmean =",round(Pmid,2)))

print(P_SD <- cellStats(Pgrd_distance, stat = 'sd', na.rm = TRUE)) # In arcGIS = 243 vs 239 = R
print(P_SUM <- cellStats(Pgrd_distance, stat = 'sum', na.rm =TRUE)) #In arcGIS = 15910387 vs 14406273 = R


#---- 3. ZERO AREAL FRACTION
#Permeable
zP <- length(which(Pgrd_distance@data@values < 1)) #Each cell is 1x1m. 
area_zP <- zP * xres(Soilgrd_distance1) * yres(Soilgrd_distance1)
print(Pz <- area_zP/areaP)

#Impermeable
zIP <- length(which(IPgrd_distance@data@values < 1)) #Each cell is 1x1m. 
area_zIP <- zIP * xres(Soilgrd_distance1) * yres(Soilgrd_distance1)
print(IPz <- area_zIP/areaIP)

#-----4. DISTANCE DISTRIBUTION RIVER NETWORK 
#---------------------------------------------------------------------------------#
#  (C)           DISTANCE DISTRIBUTION RIVER NETWORK                              #
#---------------------------------------------------------------------------------#
# Cumulative distribution functions of distances between the outlet and points in the river network
outlet.catchment <- c(xoutlet,youtlet) # coordinates of outlet from maps UTM33. Risvollan

outlet <- extract(grd_rivernet_felt,matrix(outlet.catchment,nrow=1),method='simple',cellnumbers=TRUE) # this gives the cell number, by default only provides raster cell value
outlet.coor <- xyFromCell(grd_rivernet_felt,outlet[1]) # Coordinates of the outletpoint. 
grd_rivernet_felt[outlet[1]] <- 2 

river.cells <- Which(grd_rivernet_felt==1,cells = TRUE)   
print(paste("Lengde elv=", length(river.cells), " x 1 [m] if cells eq 1 X 1 m"))

river.cells.coor<- xyFromCell(grd_rivernet_felt,river.cells)

#calculate distance from river points to outlet
river.distance.outlet     <- matrix(NA,nrow=length(river.cells),ncol=4)
river.distance.outlet[,1] <- river.cells 
river.distance.outlet[,2:3] <- river.cells.coor
for(i in 1:length(river.cells)){
  pst  <- rbind(river.cells.coor[i,],outlet.coor)
  pd   <- pointDistance(pst,lonlat=FALSE)
  river.distance.outlet[i,4] <- unique(pd[pd>0])
}

#Parameters to Paramfile.
maxFL <-0.0
midFl <-0.0
stdFL <-0.0

maxFL <- max(river.distance.outlet[,4])
midFl <- mean(river.distance.outlet[,4])
stdFL <- sd(river.distance.outlet[,4])

#---------------------------------------------------------------------------------#
#              CREATE PARAMETER FILE                                              #
#---------------------------------------------------------------------------------#

param.file <- paste("\\\\nve.no\\fil\\h\\HB\\HB-modellering\\DDDtestbenk\\DDD_Urban\\Parameters\\123.38\\Par_123.38_10min.csv",sep= "") #original param file, brukes som Mal
nyparam.file <- paste("\\\\nve.no\\fil\\h\\HB\\HB-modellering\\DDDtestbenk\\DDD_Urban\\Parameters\\",station,"\\Par_",station,"_QGIS_",res,".csv",sep="") 

prm <- read.csv(param.file,sep=";", header=FALSE)#original parameters
print(hfelt	<- ((a01+a02+a03+a04+a05+a06+a07+a08+a09+a10)/10))

prmny <- prm
prmny$V2[2] <- a00
prmny$V2[3] <- a01
prmny$V2[4] <- a02
prmny$V2[5] <- a03
prmny$V2[6] <- a04
prmny$V2[7] <- a05
prmny$V2[8] <- a06
prmny$V2[9] <- a07
prmny$V2[10] <- a08
prmny$V2[11] <- a09
prmny$V2[12] <- a10
#Parameters 13-30 as template file, needs perhaps editing!!!
prmny$V2[31] <- area
#Parameters 32-35 as template file, needs perhaps editing!!!
prmny$V2[36] <- Pfrac
prmny$V2[37] <- IPfrac
#Parameter 38 as template file, needs perhaps editing!!!
prmny$V2[39] <- Pmax  
prmny$V2[40] <- IPmax
#Parameter 41 as template file, needs perhaps editing!!!
prmny$V2[42] <- Pmid
prmny$V2[43] <- IPmid
#Parameter 44 as template file, needs perhaps editing!!!
prmny$V2[45] <- Pz
prmny$V2[46] <- IPz
#Parameter 47 as template file, needs perhaps editing!!!
prmny$V2[48] <- midFl # mean river length  
prmny$V2[49] <- stdFL # std river length
prmny$V2[50] <- maxFL # max river length
#Parameters 51-56 as template file, needs perhaps editing!!!

write.table(prmny, file=nyparam.file, row.names = FALSE,col.names = FALSE,sep=";")

#savePlot(filename = plotfile, type = c("png"), device = dev.cur(), restoreConsole = TRUE)



