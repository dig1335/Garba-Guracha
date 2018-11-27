setwd("E:/R (GG)")

library(rioja)
library(grid)
library(gridExtra)
library(dplyr)
library(vegan)
library(analogue)
library(ggpalaeo)

### MAKE DICTIONARIES FOR HABITAT TYPES TO MAKE % TYPE PLOTS LATER ###

dictionary<-read.table("dictionary.csv", header=TRUE, sep=",", stringsAsFactors = FALSE)


### conc. data seperatley to have conc data where there is no other data ###
# to show the very low conc. 

concdepth <- diatom1$"depth"
concage <- diatom1$"AgeCal"
conc <- diatom1$"valves/g dry sediment"
#reduce to make conc. more readale on axis ticks
conc <- conc/1000
conc[is.na(conc <- conc)] <- 0
conc <- data.frame(conc)

## remove empty count rows by using "total counted" 0 values
## change to NA to be able to filter by this column thereby keeping
## depth for each row that is kept for plotting
## Only for conc. data do we need ALL rows
diatom[diatom == 0] <- NA

completeFun <- function(diatom, desiredCols) {
  completeVec <- complete.cases(diatom[, desiredCols])
  return(diatom[completeVec, ])
}

diatom <- completeFun(diatom, "Total counted")
#######


# all data

diatom <- read.csv("Garba Guracha core count.csv", header=TRUE, sep=",", check.names=FALSE)

# bind age chronology to diatom dataset and filter only ages 
# needed for the depths you have

core.age <- read.csv("gg_age.csv", header=TRUE, check.names = F)
diatom <- read.csv("Garba Guracha core count.csv", header=TRUE, sep=",", check.names=F)
age.diatom <- core.age %>% filter(depth %in% diatom$depth)
age.diatom <- age.diatom [, "median"]
age.cal.diatom<-setNames(data.frame(age.diatom), c("AgeCal"))
diatom <-cbind(age.cal.diatom, diatom)

depth <- diatom$"depth"
age <- diatom$"AgeCal"
conc <- diatom$"valves/g dry sediment"/1000
conc[is.na(conc <- conc)] <- 0
conc <- data.frame(conc)

## change any NA to 0 to be processed properly
diatom[is.na(diatom <- diatom)] <- 0

# extract depths for plotting later and remove unnec. columns
diatom <- diatom[-1:-4]

#prcurve
di.pcca <- prcurve(diatom, method = "ca", trace = T, vary = T, penalty = 1.4, plotit = T)

#extract scores to plot
prc <- scores(di.pcca)

#normalised 0-1 and then by variance explained as in Bennion
prcnorm <- (prc-min(prc))/(max(prc)-min(prc))
prcnorm <- varExpl(di.pcca)*prcnorm

# Remove diatom taxa where total abundanace is less than 5%
## convert to % ##
diatom <- diatom/rowSums(diatom)*100
mxc <- apply(diatom, 2, max)
diatom <- diatom[, mxc>5]

## bind prc scores to the diatom dataset
diatom <- cbind(prc, diatom)

## remove PRC to plot seperatley 
prcsep <- diatom[1]
diatom <- diatom[-1]

#add a blank column for lithology (prcsep[2]) and zones (prcsep[3])
conc["depth"] <- depth
prcsep["Lith"] <- 0
prcsep["Zones"] <- 0
conc <- data.frame(conc)

# use for a future dendrogram from constrained cluster analysis
diss <- dist(sqrt(diatom/100)^2)
clust <- chclust(diss, method="coniss")

# broken stick model suggests significant zones
bstick(clust)

# plot limits for axes always a certain distance from eachother based on first plot
# reduces need to change them manually in the code 
# just change relative proportions here...
# Total must </= 1 to be plotted properly

Depthaxis <- 0.15
Firstplotend <- Depthaxis + 0.05
Ageaxis <- Depthaxis - 0.058
Agesaxisend <- Depthaxis - 0.02
Secondplotend <- Firstplotend + 0.55
Thirdplotend <- Secondplotend + 0.07
Fourthplotend <- Thirdplotend + 0.07
Fifthplotend <- (0.999 - Fourthplotend) + Fourthplotend

# remember to export from here if you dont use pdf() function!
windows(width=14, height=10)

blank <- strat.plot(prcsep[2], 
                    yvar = depth, y.rev=TRUE,
                    xLeft=Depthaxis, xRight=Firstplotend, yBottom = 0.075, yTop = 0.8, x.pc.lab=T, 
                    y.tks=c(0, 200, 400, 600, 800, 1000, 1200, 1400, 1600),
                    col.line = "NA", col.poly.line = "NA", col.bar = "NA", ylabel="Depth (cm)",
                    y.axis=T, srt.xlabel=0, cex.xlabel=1, xSpace=0.01)

x <- strat.plot(diatom[1:23], yvar=depth, y.axis=F,
                y.rev=TRUE, xLeft=Firstplotend, xRight = Secondplotend, yBottom = 0.075, yTop = 0.8,
                plot.line=FALSE, plot.poly=FALSE, plot.bar=TRUE,
                col.bar="blue", lwd.bar=2, scale.percent=TRUE, 
                cex.xlabel=0.9, srt.xlabel=45, add=T,
                xSpace=0.01, x.pc.lab=TRUE, x.pc.inc=20, x.pc.omit0=TRUE,
                las=2)

prc <- strat.plot(prcsep[1], 
           yvar = depth, y.rev=TRUE, x.pc.lab=TRUE,x.pc.inc=100, scale.percent=T,
           xLeft = Secondplotend, xRight=Thirdplotend, yBottom = 0.075, yTop = 0.8, x.pc.omit0=F,
           y.axis=FALSE, add=TRUE, srt.xlabel=0, cex.xlabel=1)

concplot <- strat.plot(conc[1], 
           yvar = depth, y.rev=TRUE, x.pc.lab=TRUE, x.pc.inc=100, scale.percent=T,
           xLeft = Thirdplotend, xRight=Fourthplotend, yBottom = 0.075, yTop = 0.8, x.pc.omit0=F,
           y.axis=FALSE, add=TRUE, srt.xlabel=0, cex.xlabel=1, xSpace=0.01, min.width=400,
           plot.line=T, plot.poly=FALSE, plot.bar=FALSE, plot.symb=TRUE, symb.pch=16, symb.cex=0.6, 
           col.line="Black")

zones <- strat.plot(prcsep[3], 
           yvar = depth, y.rev=TRUE, 
           xLeft = Fourthplotend, xRight=Fifthplotend, yBottom = 0.075, yTop = 0.8, x.pc.lab=T,
           col.line = "NA", col.poly.line = "Black", col.bar = "NA",
           y.axis=FALSE, add=TRUE, srt.xlabel=0, cex.xlabel=1, xSpace=0.01, 
           clust=clust, clust.width=0.05)

addClustZone(x, clust, 4, col="red", lty = 2)
addClustZone(prc, clust, 4, col="red", lty = 2)
addClustZone(concplot, clust, 4, col="red", lty = 2)
addClustZone(zones, clust, 4, col="red", lty = 2)

# from ggpalaeo

secondary.scale(yvar = depth, xLeft=Ageaxis, cex.ylabel2 = 1, yvar2 = age, n = 75, y.rev = TRUE, ylabel2 = "")

dev.off()

fix(secondary.scale)