## FDA-Final-Project: Data on primary biliary cirrhosis
## UC3M, January 2019
## Authors: Carreras García, Danae and  
##          Hernández Roig, Harold A. and 
##          Orozco Gámez, Víctor 


## --- Load & Transform Data ---
library('survival')

# all.data <- pbcseq
# # Some basics on the dataset & what we are using:
# (***) id:	case number
# (***) status:	status at endpoint, 0/1/2 for censored, transplant, dead
# (***) day:	number of days between enrollment and this visit date
# (***) bili:	serum bilirunbin (mg/dl)
# (***) futime:	number of days between registration and the earlier of death,
#             transplantion, or study analysis in July, 1986


d.bili <- data.frame(id = pbcseq$id, status = pbcseq$status,
                    day = pbcseq$day, l.bili = log(pbcseq$bili), futime = pbcseq$futime)
d.bili <- d.bili[!(d.bili$futime <= 910 & d.bili$status == 2),] # retain those who survived first 910 days
d.bili <- d.bili[d.bili$day <= 910,]   # covariates only from first 910 days
summary(d.bili$day)

# Reformulate "status": long-lived = 1 (in) & short-lived = 0
d.bili$status[d.bili$status == 0] = 1 # censored and transplant
d.bili$status[d.bili$status == 2] = 0 # dead labeled as 0

# d.bili <- d.bili[d.bili$status > 0,]
# d.bili$status[d.bili$status == 2] = 0

d.bili$status = as.factor(d.bili$status)
summary(d.bili)

#Plot individual growth curves
library(ggplot2)
p <- ggplot(d.bili, aes(x=day, y=l.bili, colour=status, group=id))+ 
  geom_line(alpha=.5) + geom_point() +
  labs(x = "Days") + labs(y = "log( serum bilirunbin )")
p

# Saving plot:
ggsave("curvas.png", plot = last_plot(), device = NULL, path = NULL,
       scale = 1, width = 12, height = 7, units = c("in"),
       dpi = 300, limitsize = TRUE)

# Convert long-formatted data into wide
w.bili <- reshape(d.bili, v.names="l.bili", timevar="day",
               idvar="id", direction="wide")
summary(w.bili$status)

## --- PACE ---

library(fdapace)
library(plot3D)
in.bili <- MakeFPCAInputs(IDs = d.bili$id, 
                           tVec=d.bili$day, 
                           yVec=d.bili$l.bili)
in.bili0 <- MakeFPCAInputs(IDs = d.bili$id[d.bili$status == 0], 
                           tVec=d.bili$day[d.bili$status == 0], 
                           yVec=d.bili$l.bili[d.bili$status == 0])
in.bili1 <- MakeFPCAInputs(IDs = d.bili$id[d.bili$status == 1], 
                           tVec=d.bili$day[d.bili$status == 1], 
                           yVec=d.bili$l.bili[d.bili$status == 1])

# PCA analysis via PACE for different subsets of status:

# all subjects
fpcaObjbili <- FPCA(in.bili$Ly, in.bili$Lt, 
                    optns = list(dataType='Sparse', 
                                 kernel = 'epan' , userBwMu = 370, userBwCov = 370, methodBwCov = 'GCV' ))
plot(fpcaObjbili)
CreateCovPlot(fpcaObjbili, covPlotType = "Smoothed", isInteractive = TRUE,
              colSpectrum = c('blue', 'red'))
# Check cum. Functional Variance Explained:
fpcaObjbili[["cumFVE"]] # 2 explains more than 97.9% and 3 explains 99.2% !!!

# just short-lived = 0
fpcaObjbili0 <- FPCA(in.bili0$Ly, in.bili0$Lt,
                     optns = list(dataType='Sparse', kernel = 'epan' , userBwMu = 370, userBwCov = 370 ))
plot(fpcaObjbili0)

fpcaObjbili1 <- FPCA(in.bili1$Ly, in.bili1$Lt,
                     optns = list(dataType='Sparse', kernel = 'epan' , userBwMu = 370, userBwCov = 370 ))
plot(fpcaObjbili1)

# Plotting all smoothed-means
muss <- data.frame(mu.a = fpcaObjbili$mu, mu0 = fpcaObjbili0$mu, 
                   mu1 = fpcaObjbili1$mu, days = fpcaObjbili$workGrid)
p <- ggplot() +
  geom_line(data=muss, aes(x = days, y = mu.a, color = "all patients", linetype = "all patients")) +
  geom_line(data=muss, aes(x = days, y = mu1, color = "long-lived",linetype = "long-lived")) +
  geom_line(data=muss, aes(x = days, y = mu0, color = "short-lived", linetype = "short-lived"))
p +  labs(color = "Eigenfunctions", linetype = "Means for...", x = "Days" , y = "Smoothed Means") + guides(color = FALSE)

# Saving plot:
ggsave("medias.png", plot = last_plot(), device = NULL, path = NULL,
       scale = 1, width = 12, height = 7, units = c("in"),
       dpi = 300, limitsize = TRUE)

# Plotting first K = 3 smoothed eigenfunctions
phi3 <- data.frame(phi = fpcaObjbili$phi[,1:3], days = fpcaObjbili$workGrid)
p <- ggplot() +
  geom_line(data=phi3, aes(x = days, y = phi.1, color = "phi_1", linetype = "phi_1")) +
  geom_line(data=phi3, aes(x = days, y = phi.2, color = "phi_2", linetype = "phi_2")) +
  geom_line(data=phi3, aes(x = days, y = phi.3, color = "phi_3", linetype = "phi_3"))
p +  labs(color = "Eigenfunctions", linetype = "Eigenfunctions", x = "Days", y = " ") + guides(color = FALSE)

# Saving plot
ggsave("phis3.png", plot = last_plot(), device = NULL, path = NULL,
       scale = 1, width = 12, height = 7, units = c("in"),
       dpi = 300, limitsize = TRUE)

# library(ks)
# CreateOutliersPlot(fpcaObjbili, optns = list(K = 2, variant = 'KDE'))

# Just checking Fitted-Curves:
par(mfrow=c(3,2),oma = c(0, 0, 2, 0))
CreatePathPlot(fpcaObjbili, subset = c(5), K = 2, showObs = TRUE) ; grid()
CreatePathPlot(fpcaObjbili, subset = c(29), K = 2, showObs = TRUE) ; grid()
CreatePathPlot(fpcaObjbili, subset = c(150), K = 2, showObs = TRUE) ; grid()
CreatePathPlot(fpcaObjbili, subset = c(71), K = 2, showObs = TRUE) ; grid()
CreatePathPlot(fpcaObjbili, subset = c(111), K = 2, showObs = TRUE) ; grid()
CreatePathPlot(fpcaObjbili, subset = c(70), K = 2, showObs = TRUE) ; grid()
mtext("Predicted Log(bilirubin) Trajectories (K = 2)", outer = TRUE, cex = 1)

# FPCS 1 vs. FPCS2: there's overlapping!

xi.12 <- data.frame(fpcs = fpcaObjbili$xiEst[,1:2], status = w.bili$status)

p <- ggplot(data=xi.12, aes(x = fpcs.1, y = fpcs.2, group=status)) +
  geom_point(aes(shape=status, color=status))
p + xlab("1st FPCS")  + ylab("2nd FPCS")

# Saving plot
ggsave("fpcs1vs2.png", plot = last_plot(), device = NULL, path = NULL,
       scale = 1, width = 12, height = 7, units = c("in"),
       dpi = 300, limitsize = TRUE)

# Clustering the curves
library(EMCluster)
A <- FClust(in.bili$Ly, in.bili$Lt, optnsFPCA = list(methodMuCovEst = 'smooth', userBwCov = 370, FVEthreshold = 0.95), k = 2)
# The Neg-Entropy Criterion can be found as: A$clusterObj@bestResult@criterionValue 
CreatePathPlot( fpcaObjbili, K=2, showObs=FALSE, lty=1, col= A$cluster, xlab = 'Days', ylab = 'log( Serum Bilirubin )')
grid()

## ---- Classifying/Discriminating  ----

# Split data into Train/Test 70/30 %
set.seed(123) # for reproducibility
train_ind<- sample(1:dim(xi.12)[1] , 0.7*nrow(xi.12))
train <- xi.12[train_ind,]
test <- xi.12[-train_ind,]

# Fitting the binomial model:
model <- glm(status~.,family=binomial(link='logit'),data=train)
summary(model) # just to check the estimated parameters

anova(model, test = "Chisq") 
# estiamte for 2nd comp. parameter is no signif. (but we keep it)
# How to interpret it:
# The difference between the null deviance and the residual deviance shows how our 
# model is doing against the null model (a model with only the intercept). The wider 
# this gap, the better.

# Accuracy 
library(ROCR)
p <- predict(model, newdata=test, type="response")
pr <- prediction(p, test$status)
prf <- performance(pr, measure = "tpr", x.measure = "fpr")
plot(prf)

auc <- performance(pr, measure = "auc")
auc <- auc@y.values[[1]]
auc # ~ 0.7548077 is not that good! Although it depends on the data splitting...
