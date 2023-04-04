#Data and code for the manuscript "Prebiotic chiral transfer from self-aminoacylating ribozymes may favor either handedness" (2023) by Kenchel et al.

#Program to takes data from variable "chiralityImportData" and fits it to kinetic model
#Additional code to generate all Figures for the manuscript is included
#Source data is stored in the variable "chiralityImportData"

#The ie of all ribozymes will be discoverable as the variable "OutputResults"
#Some graphical elements bug out depending on the size of the computer screen -- reason unknown.


############ Begin program ##########

allDataForFit <- chiralityImportData[c("Sequence.name","molarBFOconcentration","F","F_L","F_D")]

#Vector of sequence names
seqNames <- c("S-1A.1-a", "S-1A.1-n", "S-1B.1-a", "S-1B.2-a", "S-1B.3-a", "S-1C.1-a", "S-2.1-a", "S-2.1-t", "S-2.2-a", "S-3.1-a", "P4-1", "P4-2", "P4-3", "P2-4", "P1-5") 

Lcolor <- "red"
Dcolor <- "blue"
Fcolor <- "purple"

nBootReps <- 1000 #number of bootstrap replicates
BFOconc <- (1:1100)/1000000
confInterval <- 0.95

halfConfInt <- confInterval/2
confUpper <- ceiling((0.5 + halfConfInt)*nBootReps)
confLower <- ceiling((0.5 - halfConfInt)*nBootReps)
confMedian <- round(0.5*nBootReps)

#Set parameter guesses. These values seem to work
Aguess <- 0.5
kguess <- 50
kLguess <- 60
kDguess <- 10
t <- 1200
alpha <- 0.038

#Turn the data table into a list, with each sequence's data set a separate list entry
nSeqs <- length(seqNames)

dataList <- list()
for(i in 1:nSeqs){
  rowNums <- which(allDataForFit[,"Sequence.name"] == seqNames[i])
  dataList[[seqNames[i]]] <- allDataForFit[rowNums,]
}

#Just trying something
finalLowerPercent <- numeric(nSeqs)

#Set graphing parameters
par(mfrow=c(5,3), mar=c(2,2,1,1))

#Create vectors to store parameter fits
ks <- array(dim = c(nSeqs, nBootReps))
As <- array(dim = c(nSeqs, nBootReps))
kLs <- array(dim = c(nSeqs, nBootReps))
kDs <- array(dim = c(nSeqs, nBootReps))
kDkLs <- array(dim = c(nSeqs, nBootReps))
kfks <- array(dim = c(nSeqs, nBootReps))
ksSorted <- array(dim = c(nSeqs, nBootReps))
AsSorted <- array(dim = c(nSeqs, nBootReps))
kLsSorted <- array(dim = c(nSeqs, nBootReps))
kDsSorted <- array(dim = c(nSeqs, nBootReps))
kDkLSorted <- array(dim = c(nSeqs, nBootReps))
kfksSorted <- array(dim = c(nSeqs, nBootReps))
kDkLMedian <- numeric(nSeqs)
kDkLUpper <- numeric(nSeqs)
kDkLLower <- numeric(nSeqs)
kfksMedian <- numeric(nSeqs)
kfksUpper <- numeric(nSeqs)
kfksLower <- numeric(nSeqs)

#Loop through list of data tables and fit to model
for(i in c(1:9, 11:15, 10)){
  print(i)
  ktemp <- numeric(nBootReps)
  Atemp <- numeric(nBootReps)
  kLtemp <- numeric(nBootReps)
  kDtemp <- numeric(nBootReps)
  kDkLtemp <- numeric(nBootReps)
  kfkstemp <- numeric(nBootReps)
  
  allFCurves <- array(dim = c(nBootReps, length(BFOconc)))
  allLCurves <- array(dim = c(nBootReps, length(BFOconc)))
  allDCurves <- array(dim = c(nBootReps, length(BFOconc)))

  #Take a dataset out of the list
  dataToFit <- dataList[[i]]
  concs <- unique(dataToFit[,"molarBFOconcentration"])
  nConcs <- length(concs)
  nReps <- numeric(nConcs)
  for(j in 1:nConcs){
    nReps[j] <- sum(concs[j] == as.numeric(dataToFit[,"molarBFOconcentration"]))
  }

  for(j in 1:nBootReps){
    
    bootData <- NULL
    for(k in 1:nConcs){
      concRows <- which(as.numeric(dataToFit[,"molarBFOconcentration"]) == concs[k])
      if(length(concRows) > 1){
        bootData <- rbind(bootData, dataToFit[sample(concRows, size = nReps[k], replace = TRUE),])
      } else{
        bootData <- rbind(bootData, dataToFit[concRows,])
      }
    }
    
    #First fit to original kinetic model to get k and A
    kAfit <- nls(
      F ~ A*(1 - exp(-1*k*molarBFOconcentration*t*alpha)),
      data = bootData, 
      start = list(A = Aguess, k = kguess),
      lower = c(0,0),
      upper = c(1,1000),
      trace = FALSE,
      algorithm = "port"
    )
    
    #Extract k and A from fit object
    A <- coef(kAfit)[1]
    k <- coef(kAfit)[2]
    
    #Fit to separate models for L and D using values obtained for k and A
    kLfit <- nls(
      F_L ~ (A*((kL)/(2*k)))*(1 - exp(-1*(2*k)*molarBFOconcentration*0.5*t*alpha)), 
      data = bootData, 
      start = list(kL = k),
      trace = FALSE,
      lower = 0,
      upper = 1000,
      algorithm = "port"
    )
    
    kDfit <- nls(
      F_D ~ (A*((kD)/(2*k)))*(1 - exp(-1*(2*k)*molarBFOconcentration*0.5*t*alpha)), 
      data = bootData, 
      start = list(kD = k),
      trace = FALSE,
      lower = 0,
      upper = 1000,
      algorithm = "port"
    )
    
    kL <- coef(kLfit)
    kD <- coef(kDfit)
    
    ktemp[j] <- k
    Atemp[j] <- A
    kLtemp[j] <- kL
    kDtemp[j] <- kD
    kDkLtemp[j] <- kD/kL
    kfkstemp[j] <- max(c(kL, kD))/min(c(kL, kD))
    
    allFCurves[j,] <- A*(1 - exp(-1*k*BFOconc*t*alpha))
    allLCurves[j,] <- A*(kL/(kL + kD))*(1 - exp(-1*(kL + kD)*0.5*BFOconc*t*alpha))
    allDCurves[j,] <- A*(kD/(kL + kD))*(1 - exp(-1*(kL + kD)*0.5*BFOconc*t*alpha))
    
  }
  
  #store fits
  ks[i,] <- ktemp
  As[i,] <- Atemp
  kLs[i,] <- kLtemp
  kDs[i,] <- kDtemp
  kDkLs[i,] <- kDkLtemp
  kfks[i,] <- kfkstemp
  
  kstemp <- sort(ktemp)
  Astemp <- sort(Atemp)
  kLstemp <- sort(kLtemp)
  kDstemp <- sort(kDtemp)
  kDkLtemp <- sort(kDkLtemp)
  kfksstemp <- sort(kfkstemp)
  
  ksSorted[i,] <- kstemp
  AsSorted[i,] <- Astemp
  kLsSorted[i,] <- kLstemp
  kDsSorted[i,] <- kDstemp
  kDkLSorted[i,] <- kDkLtemp
  kfksSorted[i,] <- kfksstemp
  
  kDkLMedian[i] <- median(kDkLtemp)
  kDkLUpper[i] <- kDkLtemp[confUpper]
  kDkLLower[i] <- kDkLtemp[confLower]
  
  kfksMedian[i] <- median(kfksstemp)
  kfksUpper[i] <- kfksstemp[confUpper]
  kfksLower[i] <- kfksstemp[confLower]
  
  kMedian <- median(kstemp)
  kUpper <- kstemp[confUpper]
  kLower <- kstemp[confLower]
  AMedian <- median(Astemp)
  AUpper <- Astemp[confUpper]
  ALower <- Astemp[confLower]
  kLMedian <- median(kLstemp)
  kLUpper <- kLstemp[confUpper]
  kLLower <- kLstemp[confLower]
  kDMedian <- median(kDstemp)
  kDUpper <- kDstemp[confUpper]
  kDLower <- kDstemp[confLower]
  
  for(j in 1:length(BFOconc)){
    allFCurves[,j] <- sort(allFCurves[,j])
    allLCurves[,j] <- sort(allLCurves[,j])
    allDCurves[,j] <- sort(allDCurves[,j])
  }

    
    #Generate curves based on fit parameters
    fitCurveF <- allFCurves[ceiling(nBootReps/2),]
    fitCurveFUpper <- allFCurves[confUpper,]
    fitCurveFLower <- allFCurves[confLower,]
    fitCurveL <- allLCurves[ceiling(nBootReps/2),]
    fitCurveLUpper <- allLCurves[confUpper,]
    fitCurveLLower <- allLCurves[confLower,]
    fitCurveD <- allDCurves[ceiling(nBootReps/2),]
    fitCurveDUpper <- allDCurves[confUpper,]
    fitCurveDLower <- allDCurves[confLower,]
    
    finalLowerPercent[i] <- min(fitCurveL[length(fitCurveL)],fitCurveD[length(fitCurveD)])/fitCurveF[length(fitCurveF)]
    
    #jitter clustered points
    jitteredXcoords <- dataToFit[,"molarBFOconcentration"]*(10^6)
    spread <- runif(length(jitteredXcoords), -20, 20)
    jitteredXcoords <- jitteredXcoords + spread
    
    #Plot data and curves
    if(i == 14){
    plot(NULL, xlim = c(0, 1000), ylim = c(0,1), cex.main = 1.3,
         xlab = "", ylab = "", main = seqNames[i],
         yaxp = c(0, 1, 4), xaxp = c(0, 1000, 4), cex.axis = 1.3)
    } else if(sum(i == c(1,4,7)) > 0){
      plot(NULL, xlim = c(0, 1000), ylim = c(0,1),
           xlab = "", ylab = "", main = seqNames[i], cex.main = 1.3,
           yaxp = c(0, 1, 4), xaxp = c(0, 1000, 4), xaxt = 'n', cex.axis = 1.3)
      Axis(side=1, labels=FALSE, at = c(0, 250, 500, 750, 1000))
    } else if(sum(i == c(2,3,5,6,8,9)) > 0){
      #par(mar=c(1,1,1,1))
      plot(NULL, xlim = c(0, 1000), ylim = c(0,1), cex.main = 1.3,
           xlab = "", ylab = "", main = seqNames[i],
           yaxp = c(0, 1, 4), xaxp = c(0, 1000, 4),
           xaxt = 'n', yaxt = 'n', cex.axis = 1.3)
      Axis(side=1, labels=FALSE, at = c(0, 250, 500, 750, 1000))
      Axis(side=2, labels=FALSE, at = c(0, 0.25, 0.5, 0.75, 1))
    } else if(i == 11){
      plot(NULL, xlim = c(0, 1000), ylim = c(0,0.2), cex.main = 1.3,
           xlab = "", ylab = "", main = seqNames[i],
           yaxp = c(0, 0.2, 4), xaxp = c(0, 1000, 4),
           xaxt = 'n', cex.axis = 1.3)
      Axis(side=1, labels=FALSE, at = c(0, 250, 500, 750, 1000))
    } else if(sum(i == c(12, 13)) > 0){
      plot(NULL, xlim = c(0, 1000), ylim = c(0,0.2), cex.main = 1.3,
           xlab = "", ylab = "", main = seqNames[i],
           yaxp = c(0, 0.2, 4), xaxp = c(0, 1000, 4),
           xaxt = 'n', yaxt = 'n', cex.axis = 1.3)
      Axis(side=1, labels=FALSE, at = c(0, 250, 500, 750, 1000))
      Axis(side=2, labels=FALSE, at = c(0, 0.05, 0.1, 0.15, 0.2))
    } else if(sum(i == c(10, 15)) > 0){
      #par(mar=c(2,0,1,1))
      plot(NULL, xlim = c(0, 1000), ylim = c(0,1), cex.main = 1.3,
           xlab = "", ylab = "", main = seqNames[i],
           yaxp = c(0, 1, 4), xaxp = c(0, 1000, 4),
           yaxt = 'n', cex.axis = 1.3)
      Axis(side=2, labels=FALSE, at = c(0, 0.25, 0.5, 0.75, 1))
    }
    polygon(c(BFOconc*1000000, rev(BFOconc)*1000000), c(fitCurveFLower, rev(fitCurveFUpper)),
            col = adjustcolor("purple", alpha.f = 0.5), border = NA)
    polygon(c(BFOconc*1000000, rev(BFOconc)*1000000), c(fitCurveLLower, rev(fitCurveLUpper)),
            col = adjustcolor("red", alpha.f = 0.5), border = NA)
    polygon(c(BFOconc*1000000, rev(BFOconc)*1000000), c(fitCurveDLower, rev(fitCurveDUpper)),
            col = adjustcolor("blue", alpha.f = 0.5), border = NA)
    lines(x = BFOconc*1000000, y = fitCurveF, col = Fcolor)
    lines(x = BFOconc*1000000, y = fitCurveL, col = Lcolor)
    lines(x = BFOconc*1000000, y = fitCurveD, col = Dcolor)
    points(x = jitteredXcoords, y = dataToFit[,"F"], bg = Fcolor, pch = 21)
    points(x = jitteredXcoords, y = dataToFit[,"F_L"], bg = Lcolor, pch = 21)
    points(x = jitteredXcoords, y = dataToFit[,"F_D"], bg = Dcolor, pch = 21)
    if(i == 1){
      legend("topleft", legend = c("L-Phe", "D-Phe", "Both"), bg = "transparent", bty = 'n',
             pch = 21, pt.bg = c(Lcolor,Dcolor,Fcolor), cex = 1.3)
    }

}

#Reset graphing parameters
par(mfrow=c(1,1), mar=c(5.1,4.1,4.1,2.1))

#not currently used
kMedians <- ksSorted[,confMedian]
kUppers <- ksSorted[,confUpper]
kLowers <- ksSorted[,confLower]
AMedians <- AsSorted[,confMedian]
AUppers <- AsSorted[,confUpper]
ALowers <- AsSorted[,confLower]

#make data frame for Outputresults table
OutputResults <- data.frame(Ribozyme = seqNames, i_e = log2(kDkLMedian))
print(OutputResults)


########## Horizontal dotplot of ie estimates ##########

motifs = c(1, 1, 1, 1, 1, 1, 2, 2, 2, 3, 2, 2, 2, 2, 2)
dotplotData <- data.frame(ribozyme=seqNames,ie=log2(kDkLMedian),
                          ieUpper=log2(kDkLUpper),ieLower=log2(kDkLLower),motif=motifs)
dotplotData <- rbind(dotplotData, dotplotData[10,])
dotplotData <- dotplotData[c(1:9, 11:16),]

dotplotColorPalette = colorRampPalette(c("red","white","blue"), space = ("rgb"))(400)
dotplotColors <- dotplotColorPalette[ceiling((dotplotData$ie + 4)*400/8)]

errorBarY <- c(1, 4:11, 14:19)

dotchart(x=rev(dotplotData$ie), labels = rev(dotplotData$ribozyme),
         pch = 21, bg = rev(dotplotColors), cex = 2.5,
         groups = rev(dotplotData$motif))
invisible(sapply(1:15, function(i) {
  segments(rev(dotplotData$ieLower)[i], errorBarY[i], rev(dotplotData$ieUpper)[i], errorBarY[i],
           lwd = 2)
}))


########## Histogram of ie ##########

par(mar=c(6.1,6.1,2.1,2.1))
kDkLs <- kDs[,round(ncol(kDs)/2)]/kLs[,round(ncol(kLs)/2)]
enantioselectivityIndex <- log2(kDkLs)
hist(enantioselectivityIndex,
     xlab = "",
     ylab = "Frequency",
     main = "",
     breaks = 8,
     col = colorRampPalette(c("red","white","blue"), space = ("rgb"))(8),
     cex.axis = 3,
     cex.lab = 3)
title(xlab = expression(paste("Enantioselectivity index (", log2(frac("k"["D"], "k"["L"])), ")")),
      line = 5, cex.lab = 3)
par(mar=c(5.1,4.1,4.1,2.1))



########### Make plot of genetic distance vs. enantiomeric preference ###########

par(mar=c(6.1,6.1,2.1,2.1))
nSeqs <- 15
nBootReps <- 1000
confInterval <- 0.95
ieDiffs <- matrix(ncol = nSeqs, nrow = nSeqs)
ieRatios <- matrix(ncol = nSeqs, nrow = nSeqs)
kDkLRatios <- matrix(ncol = nSeqs, nrow = nSeqs)
ieDiff <- numeric(sum(1:(nSeqs-1)))
colorCode <- character(sum(1:(nSeqs-1)))
editDist <- numeric(sum(1:(nSeqs-1)))

motif1indices <- c(1,2,3,4,5,6)
motif2indices <- c(7,8,9,11,12,13,14,15)
motif3indices <- c(10)
motifIndices <- list(motif1indices,motif2indices,motif3indices)
motifIdentities <- numeric(2)

for(i in 1:nSeqs){
  for(j in 1:nSeqs){
    ieDiffs[i, j] <- abs(enantioselectivityIndex[i] - enantioselectivityIndex[j])
    ieRatios[i, j] <- enantioselectivityIndex[i]/enantioselectivityIndex[j]
    kDkLRatios[i, j] <- kDkLs[i]/kDkLs[j]
  }
}

#make vectors of ie difference and edit distnance for all unique pairs
count <- 0
for(i in 1:(nSeqs-1)){
  for(j in (i+1):nSeqs){
    for(k in 1:3){
      if(sum(i==motifIndices[[k]])>0){
        motifIdentities[1] <- k
      }
      if(sum(j==motifIndices[[k]])>0){
        motifIdentities[2] <- k
      }
    }
    
    count <- count + 1
    ieDiff[count] <- abs(enantioselectivityIndex[i] - enantioselectivityIndex[j])
    editDist[count] <- distances[i,j]
    
    #Color-code points according to motifs
    if(sum(motifIdentities == 1) == 2){
      colorCode[count] <- "blue"
    } else if(sum(motifIdentities == 2) == 2){
      colorCode[count] <- "red"
    } else if(sum(motifIdentities == 3) == 2){
      colorCode[count] <- "white"
    } else if(sum(motifIdentities == 1)==1 & sum(motifIdentities == 2)==1){
      colorCode[count] <- "purple"
    } else if(sum(motifIdentities == 1)==1 & sum(motifIdentities == 3)==1){
      colorCode[count] <- "green"
    } else if(sum(motifIdentities == 2)==1 & sum(motifIdentities == 3)==1){
      colorCode[count] <- "orange"
    }
  }
}
meanieDiffs <- numeric(nSeqs + 1)
sdieDiffs <- numeric(nSeqs + 1)
meanieRatios <- numeric(nSeqs + 1)
for(i in 0:nSeqs){
  meanieDiffs[i + 1] <- mean(ieDiffs[which(distances == i)])
  sdieDiffs[i + 1] <- sd(ieDiffs[which(distances == i)])
  meanieRatios[i + 1] <- mean(ieRatios[which(distances == i)])
}

ievsDistance <- data.frame(editDistance = editDist, 
                           ieDifference = ieDiff, 
                           color = colorCode)

#figure out which edit distances have 3 or more points
errorBarDists <- numeric()
for(i in 1:nSeqs){
  nPoints <- sum(editDist == i)
  if(nPoints >= 3){
    errorBarDists <- c(errorBarDists, i)
  }
}
nErrorBars <- length(errorBarDists)

#calculate SEM for those edit distances
errorBarSEM <- numeric(nErrorBars)
for(i in 1:nErrorBars){
  errorData <- ieDiff[which(editDist == errorBarDists[i])]
  errorBarSEM[i] <- sd(errorData)/sqrt(length(errorData))
}

#determine points to jitter
pointsToJitter <- logical(length(editDist))
for(i in 1:length(ieDiff)){
  editDistSubst <- editDist
  editDistSubst[i] <- 0
  possibleClusters <- which(editDistSubst == editDist[i])
  if(sum(abs(ieDiff[i] - ieDiff[possibleClusters]) < 0.2) > 0){
    pointsToJitter[i] <- TRUE
  }
}
editDistJittered <- editDist*(!pointsToJitter) + (editDist + rnorm(length(editDist), 0, 0.1))*pointsToJitter
plot(x = editDistJittered,
     y = ievsDistance[, "ieDifference"],
     xlim = c(0,15),
     xlab = "Edit distance (nt)",
     ylab = expression(paste("Difference in enantioselectivity (", Delta, "i"["e"], ")")),
     xaxp = c(0, 15, 3),
     cex.axis = 2,
     cex.lab = 2,
     cex = 2,
     pch = 21,
     bg = as.character(ievsDistance[, "color"])
)
axis(side = 1, at = c(0:15), labels = FALSE)
par(mar=c(5.1,4.1,4.1,2.1))




########## L-RNA vs D-RNA horizontal barplot ##########

par(mar=c(5.1,8.1,2.1,2.1))
D1ApercenteeD <- c(0.7500, 0.6604, 0.7957)
L1ApercenteeD <- c(-0.8279, -0.8701, -0.8239)
Lmean <- mean(L1ApercenteeD)
Dmean <- mean(D1ApercenteeD)
Lerrbar <- 2*sd(L1ApercenteeD)
Derrbar <- 2*sd(D1ApercenteeD)
Llower <- Lmean - Lerrbar
Lupper <- Lmean + Lerrbar
Dlower <- Dmean - Derrbar
Dupper <- Dmean + Derrbar
meanPercenteeD <- c(Lmean, Dmean)
lowerErrbar <- c(Llower, Dlower)
upperErrbar <- c(Lupper, Dupper)
barplotColors <- c(adjustcolor("red", alpha.f = abs(Lmean)),
                   adjustcolor("blue", alpha.f = abs(Dmean)))
savedBarplot <- barplot(meanPercenteeD, xlim = c(-1, 1), horiz = TRUE,
                        col = barplotColors,
                        ylab = "",
                        xlab = expression(paste("(F"["D"], " - ", "F"["L"], ")", "/", "(F"["D"], " + ", "F"["L"], ")")),
                        names.arg = c("L", "D"),
                        cex.axis = 2,
                        cex.lab = 2,
                        cex.names = 2)
arrows(lowerErrbar, savedBarplot, upperErrbar, savedBarplot,
       angle = 90, code = 3, length = 0.05, lwd = 2)
points(c(D1ApercenteeD,L1ApercenteeD), y=c(1.8,1.8,1.8,0.6,0.6,0.6), pch = 16, cex=1, lwd=1.3)
title(ylab = "S-1A.1-a\nRNA handedness", cex.lab = 2, line = 4)
par(mar=c(5.1,4.1,4.1,2.1))
