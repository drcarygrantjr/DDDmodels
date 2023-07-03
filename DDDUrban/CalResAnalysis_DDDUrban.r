#Scriptet skal sejekke området for GshInt og GscInt samt cerlerities for Eline sin kalibrering 7.8.2019

rm(list=ls())    # removes old variables

pathinn <-"\\\\nve.no\\fil\\h\\HB\\HB-modellering\\DDDtestbenk\\DDD_Urban\\Dottyplots\\"
path <-"\\\\nve.no\\fil\\h\\HB\\HB-modellering\\DDDtestbenk\\DDD_Urban\\"
st_file <- paste(path,"inndata\\St_list.txt",sep="")

st_list <- scan(st_file, character(0))
catchment <- st_list[6] # 1 er Ila, 2 er GK 3 er Mini GK, 4 er RIS, 5 er Vestli, 6 er Sandsli
TR <- "10min"  #tidsoppløsning
DRN <-1
plotfile <- paste("Dottyplots_Urban_",TR,"_",catchment,"_SL4.tif",sep="")

r2fil <- paste(path,"utdata\\r2fil_",catchment,"_SL4.csv",sep="")
#Stasjon	pro	pkorr	skorr	TX	TS	CX	Cglac	cea	rv	gtcel	hmid	a0	d	MAD	area		bogfrac	Gshape	Gscale	GshInt	GscInt	midDL	midFl	glacfrac

lowlim <- 0.0
TR <-"10min"

calres <- read.csv(r2fil,header = FALSE,sep= ";") #reads from metric file
#     u,        pro,         TX,      PCritFlux,    OFP,        GshInt,   GscInt,
#     u,        pro,         TX,      PCritFlux,    OFP,        GshInt,   GscInt, Persons 
navn <- c("NSE","KGE","BIAS","u","pro", "TX","PCritFlux","OFP","GShInt","GscInt")
names(calres)[1:10] <-navn[1:10]

#tull <- which(calres[,4] ==0) # da er verdien NA
#if(length(tull)>0) calres[tull,2:6] <-0.0  

parampath <-"\\\\nve.no\\fil\\h\\HB\\HB-modellering\\DDDtestbenk\\DDD_Urban\\Parameters\\"
param.file <-paste(parampath,"Best_par_",catchment,"_10min_SL4.csv",sep="")# I første omgang har vi bare ptq filer for 24 timer med elevation    #REDIGER !!
par <- read.csv(param.file,header = FALSE,sep= ";")

lengde <- length(calres$pro)
fart <- vector("numeric", lengde)
for (i in 1:lengde) fart[i] <- round(calres$GShInt[i]*calres$GscInt[i]*par$V2[42]/par$V2[28] ,6)

windows(16,16)
par(mfrow=c(4,3))

bestv <- max(calres$KGE)
loc1 <- which(calres$KGE==bestv)#finner plasseringene til høyeste KGE 
maxNSE <- max(calres$NSE[loc1]) #highest NSE amongst the highest KGE
loc2 <- which(calres$NSE[loc1]==maxNSE) # finds the location of the highest NSE amongst the highest KGE

best <- loc1[loc2] # Beste kombinasjona av beste KGE og den beste NSE blant beste KGE

plot(calres$u,calres$KGE,ylim=c(lowlim,1.1), ylab ="KGE",xlab="u", pch=1, main = catchment)
points(calres$u[best],calres$KGE[best], pch=19,cex=1.5, col="red")
text(calres$u[best],0.99, round(calres$u[best],3),cex=1.2)

plot(calres$pro,calres$KGE,ylim=c(lowlim,1.1), ylab ="KGE",xlab="pro", pch=1, main = paste("max KGE=",round(max(calres$KGE),3)))
points(calres$pro[best],calres$KGE[best], pch=19, cex=1.5,col="red")
text(calres$pro[best],0.99, round(calres$pro[best],3),cex=1.2)

plot(calres$OFP,calres$KGE,ylim=c(lowlim,1.1), ylab ="KGE",xlab="OFP", pch=1)
points(calres$OFP[best],calres$KGE[best], pch=19, cex=1.5,col="red")
text(calres$OFP[best],0.99, round(calres$OFP[best],3),cex=1.2)

plot(calres$GShInt,calres$KGE,ylim=c(lowlim,1.1), ylab ="KGE",xlab="Gshape", pch=1)
points(calres$GShInt[best],calres$KGE[best], pch=19,cex=1.5, col="red")
text(calres$GShInt[best],0.99, round(calres$GShInt[best],3),cex=1.2)

plot(calres$GscInt,calres$KGE,ylim=c(lowlim,1.1),xlim=c(0.0,0.004), ylab ="KGE",xlab="Gscale", pch=1)
points(calres$GscInt[best],calres$KGE[best], pch=19,cex=1.5, col="red")
text(calres$GscInt[best],0.99, round(calres$GscInt[best],3),cex=1.2)

plot(calres$PCritFlux,calres$KGE,ylim=c(lowlim,1.1), ylab ="KGE",xlab="PCritFlux", pch=1)
points(calres$PCritFlux[best],calres$KGE[best], pch=19, cex=1.5,col="red")
text(calres$PCritFlux[best],0.99, round(calres$PCritFlux[best],4),cex=1.2)

plot(calres$TX,calres$KGE,ylim=c(lowlim,1.1), ylab ="KGE",xlab="TX", pch=1)
points(calres$TX[best],calres$KGE[best], pch=19, cex=1.5,col="red")
text(calres$TX[best],0.99, round(calres$TX[best],3),cex=1.2)

plot(fart,calres$KGE,ylim=c(lowlim,1.1),xlim=c(0.0,0.0005), ylab ="KGE",xlab="Velocity", pch=1)
points(fart[best],calres$KGE[best], pch=19, cex=1.5,col="red")
text(fart[best]*0.5,0.99, round(fart[best],6),cex=1.2)

#plot(calres$Persons,calres$KGE,ylim=c(lowlim,1.1), ylab ="KGE",xlab="persons", pch=1)
#points(calres$Persons[best],calres$KGE[best], pch=19, cex=1.5,col="red")
#text(calres$Persons[best]*0.5,0.99, round(calres$Persons[best],6),cex=1.2)

print(bestv)
print(maxNSE)
print(calres$BIAS[best])
print(fart[best])

savePlot(filename = paste(pathinn,plotfile, sep=""), type = c("tif"), device = dev.cur(), restoreConsole = TRUE)
