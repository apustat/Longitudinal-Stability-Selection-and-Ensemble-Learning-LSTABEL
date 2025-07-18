library(dplyr) #for data processing
library(glmmLasso) #for lasso
library(PGEE) #for penalized gee
library(glmtoolbox) #for gee
library(made4)
library(foreach) #for parallel computing
library(parallel) #for parallel computing
library(doParallel) #for parallel computing
library(Matrix)
library(lme4)
library(CVXR) #for optimization
library(REEMtree) #for random effect model
library(nlme) #for linear mixed model
library(lmmen)
library(GMMBoost) #for boosting
library(dplyr)

df=read.csv("C:/Users/....../final.data.adas13.csv")

continuous_vars <- df[,-c(1,2)]  # Select all numeric (continuous) covariates
#binary_var <- dat1[, (p+1)]  # Select the binary variable
standardized_continuous_vars <- scale(continuous_vars)

df1=cbind(standardized_continuous_vars, df[,c(1,2)])

N=length(unique(df1$RID))
p=dim(df1)[2]-3  #number of covariates
B=50
threshold=0.70
set.seed(6578)
train=sample(unique(df1$RID), N*0.6, replace=FALSE)
df.train=subset(df1,RID %in% train) #training data
df.test=subset(df1,!(RID %in% train)) #testing data
lam <- seq(0,500, length=30)
BIC_vec <- rep(NA, length(lam))
AIC_vec <- rep(NA, length(lam))

preds=colnames(df1)[1:p] #excluding RID, VISCODE, ADAS13
fm <- as.formula(paste("ADAS13~", paste(preds, collapse = "+")))

set.seed(100)
for (j in 1:length(lam)){
  print(paste("Iteration ", j, sep=""))
  glm1 <- glmmLasso(fm, rnd = list(RID=~1),
                    family = gaussian(link = "identity"), data = df.train, lambda=lam[j],
                    switch.NR=TRUE, final.re = FALSE, control=list(print.iter=TRUE))
  BIC_vec[j]<-glm1$bic
  AIC_vec[j]<-glm1$aic
  
}

best.lam=lam[which.min(BIC_vec)]
               
############################################################################
mat.coef.lasso1=matrix(NA, p, B); mat.coef.lasso2=matrix(NA, p, B)
for (j in 1:B){
  set.seed(1000+j)
  lasso.subsam=sample(c(train), length(train)*0.5, replace=FALSE)
  df.lasso.sub.sam1=subset(df.train, RID %in% lasso.subsam)
  df.lasso.sub.sam2=subset(df.train,!(RID %in% lasso.subsam))
  mod.lasso1 <- glmmLasso(fm, rnd = list(RID=~1),
                          family = gaussian(link = "identity"), data = df.lasso.sub.sam1, lambda=best.lam,
                          switch.NR=TRUE, final.re = TRUE, control=list(print.iter=TRUE))
  coef.lasso1=mod.lasso1$coefficients[-1]
  mat.coef.lasso1[,j]=coef.lasso1
  
  mod.lasso2 <- glmmLasso(fm, rnd = list(RID=~1),
                          family = gaussian(link = "identity"), data = df.lasso.sub.sam2, lambda=best.lam,
                          switch.NR=TRUE, final.re = TRUE, control=list(print.iter=TRUE))
  coef.lasso2=mod.lasso2$coefficients[-1]
  mat.coef.lasso2[,j]=coef.lasso2
} 

coef.lasso.cont=cbind(mat.coef.lasso1,mat.coef.lasso2)
coef.lasso.dis=ifelse(coef.lasso.cont !=0 ,1,0)
freq.lasso=apply(coef.lasso.dis, 1, sum)
final.set=which(freq.lasso>(2*B*threshold))

########################superlearner cross-validation#############################
set.seed(6578)
nfold=5

# Create a vector of unique subject IDs
subject_ids <- unique(df.train$RID)

# Randomly shuffle the subject IDs and divide into folds
folds <- cut(seq(1, length(subject_ids)), breaks = nfold, labels = FALSE)
rand_folds=sample(folds)

# Create a dataframe with subject IDs and their assigned fold
subject_folds <- data.frame(RID = subject_ids, Fold = rand_folds)

# Merge the fold information back into the original dataset
df.train1 <- df.train %>%
  left_join(subject_folds, by = "RID")

#Create subsets based on the new subset variable
subset1 <- subset(df.train1, Fold == 1)[, -(p+4)]
subset2 <- subset(df.train1, Fold == 2)[, -(p+4)]
subset3 <- subset(df.train1, Fold == 3)[, -(p+4)]
subset4 <- subset(df.train1, Fold == 4)[, -(p+4)]
subset5 <- subset(df.train1, Fold == 5)[, -(p+4)]

fold1=unique(subset1$RID)
fold2=unique(subset2$RID)
fold3=unique(subset3$RID)
fold4=unique(subset4$RID)
fold5=unique(subset5$RID)

all.fold <- list(fold1,fold2,fold3,fold4,fold5)
combinations <- combn(all.fold, (nfold-1), simplify = FALSE) # Generate all combinations of 4 vectors
combined.index <- lapply(combinations, function(comb) do.call(c, comb)) # Perform the combination operation (e.g., concatenate)

###############################################################################
###############linear mixed model##############################

###picking the last visit of outcome ADAS13 only
true.y.test <- df.test %>%
  group_by(RID) %>%
  arrange(VISCODE, .by_group = TRUE) %>%  # Ensure data is ordered by VISCODE within each group
  slice_tail(n = 1) %>%                   # Select the last row for each group
  ungroup() %>%
  pull(ADAS13)
#####################################################
##################GLMM###########

fm.glm <- as.formula(paste("ADAS13  ~", paste(colnames(df.train)[final.set], collapse = " + "), "+ (1 | RID)"))
pred.glm.vec=vector("numeric")
for (l in 1:nfold){
  mod.glm <- lmer(fm.glm, data = subset(df.train, RID %in% c(combined.index[[l]])))
  ndf.all.rm=subset(df.train, !(RID %in% combined.index[[l]])) #testing sets with all repeated measures
  df_last_measure <- ndf.all.rm %>%
    group_by(RID) %>% 
    slice_max(VISCODE) %>%  # Pick the row with the maximum VISCODE (i.e., last visit)
    ungroup()
  ndf=df_last_measure[, c(final.set, p+2)]
  pred.glm=predict(mod.glm, newdata=ndf, re.form=~(1|RID), allow.new.levels = TRUE) #for test set
  pred.glm.vec=c(pred.glm.vec, pred.glm)
}

######################################################################
###################generalized estimating equations -gee##############
fm.gee <- as.formula(paste("ADAS13 ~", paste(colnames(df.train)[final.set], collapse = " + ")))
pred.gee.vec=vector("numeric")
for (l in 1:nfold){
  mod.gee=glmgee(fm.gee, data= subset(df.train, RID %in% c(combined.index[[l]])), id=RID, family=gaussian(), #doesn't have method=ML or REML options as it's not likelihood based approach
                 corstr="AR-M-dependent(1)",tol = 10^-4, start= rep(0, (length(final.set)+1)))
  ndf.all.rm=subset(df.train, !(RID %in% combined.index[[l]])) #testing sets with all repeated measures
  df_last_measure <- ndf.all.rm %>%
    group_by(RID) %>% 
    slice_max(VISCODE) %>%  # Pick the row with the maximum VISCODE (i.e., last visit)
    ungroup()
  ndf=df_last_measure[, c(final.set, p+2)]
  pred.gee=predict(mod.gee, newdata=ndf, se.fit = FALSE, varest="robust")
  pred.gee.vec=c(pred.gee.vec, pred.gee)
} 

#######################Randomforest##############################
fm.rf <- as.formula(paste("ADAS13 ~", paste(colnames(df.train)[final.set], collapse = " + ")))
pred.rf.vec=vector("numeric")
true.y.vec=c()  #keeping record of true y for cross-validation only
for (l in 1:nfold){
  mod.rf=REEMtree(fm.rf, data=subset(df.train, RID %in% c(combined.index[[l]])), random=~1|RID, correlation=corAR1(), method="ML")
  ndf.all.rm=subset(df.train, !(RID %in% combined.index[[l]])) #testing sets with all repeated measures
  df_last_measure <- ndf.all.rm %>%
    group_by(RID) %>% 
    slice_max(VISCODE) %>%  # Pick the row with the maximum VISCODE (i.e., last visit)
    ungroup()
  ndf.rf=data.frame(df_last_measure[, c("ADAS13",paste(colnames(df.train)[final.set]),"RID")])
  true.y=ndf.rf$ADAS13
  true.y.vec=c(true.y.vec, true.y)
  pred.rf=predict(mod.rf,newdata=ndf.rf, id=ndf.rf$RID, EstimateRandomEffects=TRUE)
  pred.rf.vec=c(pred.rf.vec, pred.rf)
}



####################################################################
####################Superlearner####################################
pred.mat=data.frame(pred.glm.vec, pred.gee.vec, pred.rf.vec)
names(pred.mat)=c("GLM", "GEE", "RF")

w <- Variable(dim(pred.mat)[2])
objective <- Minimize(sum((true.y.vec - pred.mat %*% w)^2))

C<-matrix(rep(1,dim(pred.mat)[2]), nrow=1)
constraint1 <- w[1]>=0
constraint2 <- w[2]>=0
constraint3 <- w[3]>=0
constraint4 <- C%*%w == 1
problem <- Problem(objective, constraints = list(constraint1, constraint2, constraint3,  constraint4))
result <- solve(problem)
wght = round(result$getValue(w),3)


######################Prediction on testing set########################
mod.glm.test <- lmer(fm.glm, data = df.train)
df_last_measure_test <- df.test %>%
  group_by(RID) %>% 
  slice_max(VISCODE) %>%  # Pick the row with the maximum VISCODE (i.e., last visit)
  ungroup()
ndf.test=df_last_measure_test[, c(final.set, p+2)]
pred.glm.test=predict(mod.glm.test, newdata=ndf.test, re.form=~(1|RID), allow.new.levels = TRUE) #for test set
pmse.glm.test=mean((true.y.test-pred.glm.test)^2) #PMSE from GLMM
mape.glm.test=mean(abs(true.y.test-pred.glm.test)) #MAPE from GLMM
rpmse.glm.test=sqrt(pmse.glm.test) #RPMSE from GLMM

mod.gee.test=glmgee(fm.gee, data=df.train, id=RID, family=gaussian(), #doesn't have method=ML or REML options as it's not likelihood based approach
                    corstr="AR-M-dependent(1)",tol = 10^-4, start= rep(0, (length(final.set)+1)))
pred.gee.test=predict(mod.gee.test, newdata=ndf.test, se.fit = FALSE, varest="robust")
pmse.gee.test=mean((true.y.test-pred.gee.test)^2) #PMSE from GEE
mape.gee.test=mean(abs(true.y.test-pred.gee.test)) #MAPE from GEE
rpmse.gee.test=sqrt(pmse.gee.test) #RPMSE from GEE

mod.rf.test<-REEMtree(fm.rf, data=df.train, random=~1|RID, correlation=corAR1(), method="ML")
newdata.test=data.frame(df_last_measure_test[, c("ADAS13",paste(colnames(df.train)[final.set]),"RID")])
pred.rf.test=predict(mod.rf.test,newdata=newdata.test, id=newdata.test$RID, EstimateRandomEffects=TRUE)
pmse.rf.test=mean((true.y.test-pred.rf.test)^2) #PMSE from RF
mape.rf.test=mean(abs(true.y.test-pred.rf.test)) #MAPE from RF
rpmse.rf.test=sqrt(pmse.rf.test) #RPMSE from RF

pred.mat.test=data.frame(pred.glm.test, pred.gee.test, pred.rf.test)
pred.el.test=as.matrix(pred.mat.test)%*%wght
pmse.el.test=mean((true.y.test-pred.el.test)^2) #PMSE from ensemble learning
mape.el.test=mean(abs(true.y.test-pred.el.test)) #MAPE from ensemble learning
rpmse.el.test=sqrt(pmse.el.test) #RPMSE from ensemble learning


###########################################################################################
#####################GLMMLASSO######################################
###########################################################################################
lasso.mod <- glmmLasso(fm, rnd = list(RID=~1),
                       family = gaussian(link="identity"), data = df.train, lambda=best.lam,
                       switch.NR=TRUE, final.re = TRUE, control=list(print.iter=TRUE))
final.set.lasso=which(lasso.mod$coefficients[-1] !=0) #selected variables by GLMMLASSO

pred.test.lasso=predict(lasso.mod, data.frame(df_last_measure_test))
pmse.test.lasso=mean((true.y.test-pred.test.lasso)^2) #PMSE from GLMMLASSO
mape.test.lasso=mean(abs(true.y.test-pred.test.lasso)) # MAPE from GLMMLASSO
rpmse.test.lasso=sqrt(pmse.test.lasso) #RPMSE from GLMMLASSO

###########################################################################################
#####################GLMMLASSO with ensemble learning######################################
###########################################################################################

##################GLMM###########
fm.glm.lasso <- as.formula(paste("ADAS13  ~", paste(colnames(df.train)[final.set.lasso], collapse = " + "), "+ (1 | RID)"))
pred.glm.vec.lasso=vector("numeric")
for (l in 1:nfold){
  mod.glm.lasso <- lmer(fm.glm.lasso, data = subset(df.train, RID %in% c(combined.index[[l]])))
  ndf.all.rm.lasso=subset(df.train, !(RID %in% combined.index[[l]])) #testing sets with all repeated measures
  df_last_measure_lasso <- ndf.all.rm.lasso %>%
    group_by(RID) %>% 
    slice_max(VISCODE) %>%  # Pick the row with the maximum VISCODE (i.e., last visit)
    ungroup()
  ndf.lasso=df_last_measure_lasso[, c(final.set.lasso, p+2)]
  pred.glm.lasso=predict(mod.glm.lasso, newdata=ndf.lasso, re.form=~(1|RID), allow.new.levels = TRUE) #for test set
  pred.glm.vec.lasso=c(pred.glm.vec.lasso, pred.glm.lasso)
}

###################generalized estimating equations -gee##############
fm.gee.lasso <- as.formula(paste("ADAS13 ~", paste(colnames(df.train)[final.set.lasso], collapse = " + ")))
pred.gee.vec.lasso=vector("numeric")
for (l in 1:nfold){
  mod.gee.lasso=glmgee(fm.gee.lasso, data= subset(df.train, RID %in% c(combined.index[[l]])), id=RID, family=gaussian(), #doesn't have method=ML or REML options as it's not likelihood based approach
                 corstr="AR-M-dependent(1)",tol = 10^-4, start= rep(0, (length(final.set.lasso)+1)))
  ndf.all.rm.lasso=subset(df.train, !(RID %in% combined.index[[l]])) #testing sets with all repeated measures
  df_last_measure_lasso <- ndf.all.rm.lasso %>%
    group_by(RID) %>% 
    slice_max(VISCODE) %>%  # Pick the row with the maximum VISCODE (i.e., last visit)
    ungroup()
  ndf.lasso=df_last_measure_lasso[, c(final.set.lasso, p+2)]
  pred.gee.lasso=predict(mod.gee.lasso, newdata=ndf.lasso, se.fit = FALSE, varest="robust")
  pred.gee.vec.lasso=c(pred.gee.vec.lasso, pred.gee.lasso)
} 

#######################Randomforest##############################
fm.rf.lasso <- as.formula(paste("ADAS13 ~", paste(colnames(df.train)[final.set.lasso], collapse = " + ")))
pred.rf.vec.lasso=vector("numeric")
for (l in 1:nfold){
  mod.rf.lasso=REEMtree(fm.rf.lasso, data=subset(df.train, RID %in% c(combined.index[[l]])), random=~1|RID, correlation=corAR1(), method="ML")
  ndf.all.rm.lasso=subset(df.train, !(RID %in% combined.index[[l]])) #testing sets with all repeated measures
  df_last_measure_lasso <- ndf.all.rm.lasso %>%
    group_by(RID) %>% 
    slice_max(VISCODE) %>%  # Pick the row with the maximum VISCODE (i.e., last visit)
    ungroup()
  ndf.rf.lasso=data.frame(df_last_measure_lasso[, c("ADAS13",paste(colnames(df.train)[final.set.lasso]),"RID")])
  pred.rf.lasso=predict(mod.rf.lasso,newdata=ndf.rf.lasso, id=ndf.rf.lasso$RID, EstimateRandomEffects=TRUE)
  pred.rf.vec.lasso=c(pred.rf.vec.lasso, pred.rf.lasso)
}

###########################################################################
#####################Ensemble learning####################################
pred.mat.lasso=data.frame(pred.glm.vec.lasso, pred.gee.vec.lasso, pred.rf.vec.lasso)
names(pred.mat.lasso)=c("GLMM", "GEE", "RF")

w.lasso <- Variable(dim(pred.mat.lasso)[2])
objective.lasso <- Minimize(sum((true.y.vec - pred.mat.lasso %*% w.lasso)^2))

C.lasso<-matrix(rep(1,dim(pred.mat.lasso)[2]), nrow=1)
constraint1.lasso <- w.lasso[1]>=0
constraint2.lasso <- w.lasso[2]>=0
constraint3.lasso <- w.lasso[3]>=0
constraint4.lasso <- C.lasso%*%w.lasso == 1
problem.lasso <- Problem(objective.lasso, constraints = list(constraint1.lasso,  constraint2.lasso,  constraint3.lasso, constraint4.lasso))
result.lasso <- solve(problem.lasso)
wght.lasso = round(result.lasso$getValue(w.lasso),3)

mod.glm.test.lasso = lmer(fm.glm.lasso, data = df.train)
ndf.test.lasso=df_last_measure_test[, c(final.set.lasso, p+2)] #p+2 is for RID
pred.glm.test.lasso=predict(mod.glm.test.lasso, newdata=ndf.test.lasso, re.form=~(1|RID), allow.new.levels = TRUE) #for test set
pmse.glm.test.lasso=mean((true.y.test-pred.glm.test.lasso)^2) #PMSE from GLMML
mape.glm.test.lasso=mean(abs(true.y.test-pred.glm.test.lasso)) #MAPE from GLMML
rpmse.glm.test.lasso=sqrt(pmse.glm.test.lasso) #RPMSE from GLMM

mod.gee.test.lasso=glmgee(fm.gee.lasso, data=df.train, id=RID, family=gaussian(), #doesn't have method=ML or REML options as it's not likelihood based approach
                    corstr="AR-M-dependent(1)",tol = 10^-4, start= rep(0, (length(final.set.lasso)+1)))
pred.gee.test.lasso=predict(mod.gee.test.lasso, newdata=ndf.test.lasso, se.fit = FALSE, varest="robust")
pmse.gee.test.lasso=mean((true.y.test-pred.gee.test.lasso)^2) #PMSE from GEE
mape.gee.test.lasso=mean(abs(true.y.test-pred.gee.test.lasso)) #MAPE from GEE
rpmse.gee.test.lasso=sqrt(pmse.gee.test.lasso) #RPMSE from GEE

mod.rf.test.lasso = REEMtree(fm.rf.lasso, data=df.train, random=~1|RID, correlation=corAR1(), method="ML")
newdata.test.lasso=data.frame(df_last_measure_test[, c("ADAS13",paste(colnames(df.train)[final.set.lasso]),"RID")])
pred.rf.test.lasso=predict(mod.rf.test.lasso,newdata=newdata.test.lasso, id=newdata.test.lasso$RID, EstimateRandomEffects=TRUE)
pmse.rf.test.lasso=mean((true.y.test-pred.rf.test.lasso)^2) #PMSE from RF
mape.rf.test.lasso=mean(abs(true.y.test-pred.rf.test.lasso)) #MAPE from RF
rpmse.rf.test.lasso=sqrt(pmse.rf.test.lasso) #RPMSE from RF

pred.mat.test.lasso=data.frame(pred.glm.test.lasso, pred.gee.test.lasso, pred.rf.test.lasso)
pred.el.test.lasso=as.matrix(pred.mat.test.lasso)%*%wght.lasso
pmse.el.test.lasso=mean((true.y.test-pred.el.test.lasso)^2) #PMSE from ensemble learning
mape.el.test.lasso=mean(abs(true.y.test-pred.el.test.lasso)) #MAPE from ensemble learning
rpmse.el.test.lasso=sqrt(pmse.el.test.lasso) #RPMSE from ensemble learning

#########################################################################
#####################Ensemble learning with all biomarkers################
#########################################################################

##################GLMM###########
all.fm.glm <- as.formula(paste("ADAS13  ~", paste(colnames(df.train)[1:p], collapse = " + "), "+ (1 | RID)"))
all.pred.glm.vec=vector("numeric")
for (l in 1:nfold){
  all.mod.glm <- lmer(all.fm.glm, data = subset(df.train, RID %in% c(combined.index[[l]])))
  all.ndf.all.rm=subset(df.train, !(RID %in% combined.index[[l]])) #testing sets with all repeated measures
  df_last_measure <- all.ndf.all.rm %>%
    group_by(RID) %>% 
    slice_max(VISCODE) %>%  # Pick the row with the maximum VISCODE (i.e., last visit)
    ungroup()
  all.ndf=df_last_measure[, c(1:p, p+2)]
  all.pred.glm=predict(all.mod.glm, newdata=all.ndf, re.form=~(1|RID), allow.new.levels = TRUE) #for test set
  all.pred.glm.vec=c(all.pred.glm.vec, all.pred.glm)
}

###################generalized estimating equations -gee##############
all.fm.gee <- as.formula(paste("ADAS13 ~", paste(colnames(df.train)[1:p], collapse = " + ")))
all.pred.gee.vec=vector("numeric")
for (l in 1:nfold){
  all.mod.gee=glmgee(all.fm.gee, data= subset(df.train, RID %in% c(combined.index[[l]])), id=RID, family=gaussian(), #doesn't have method=ML or REML options as it's not likelihood based approach
                       corstr="AR-M-dependent(1)",tol = 10^-4, start= rep(0, p+1))
  all.ndf.all.rm=subset(df.train, !(RID %in% combined.index[[l]])) #testing sets with all repeated measures
  df_last_measure <- all.ndf.all.rm %>%
    group_by(RID) %>% 
    slice_max(VISCODE) %>%  # Pick the row with the maximum VISCODE (i.e., last visit)
    ungroup()
  all.ndf=df_last_measure
  all.pred.gee=predict(all.mod.gee, newdata=all.ndf, se.fit = FALSE, varest="robust")
  all.pred.gee.vec=c(all.pred.gee.vec, all.pred.gee)
} 

#######################Randomforest##############################
all.fm.rf <- as.formula(paste("ADAS13 ~", paste(colnames(df.train)[1:p], collapse = " + ")))
all.pred.rf.vec=vector("numeric")
for (l in 1:nfold){
  all.mod.rf=REEMtree(all.fm.rf, data=subset(df.train, RID %in% c(combined.index[[l]])), random=~1|RID, correlation=corAR1(), method="ML")
  all.ndf.all.rm=subset(df.train, !(RID %in% combined.index[[l]])) #testing sets with all repeated measures
  df_last_measure <- all.ndf.all.rm %>%
    group_by(RID) %>% 
    slice_max(VISCODE) %>%  # Pick the row with the maximum VISCODE (i.e., last visit)
    ungroup()
  ndf.rf=data.frame(df_last_measure[, c("ADAS13",paste(colnames(df.train)[1:p]),"RID")])
  all.pred.rf=predict(all.mod.rf,newdata=ndf.rf, id=ndf.rf$RID, EstimateRandomEffects=TRUE)
  all.pred.rf.vec=c(all.pred.rf.vec, all.pred.rf)
}

###########################################################################
#####################Ensemble learning####################################
all.pred.mat=data.frame(all.pred.glm.vec, all.pred.gee.vec, all.pred.rf.vec)
names(all.pred.mat)=c("GLMM", "GEE", "RF")
all.w <- Variable(dim(all.pred.mat)[2])
all.objective <- Minimize(sum((true.y.vec - all.pred.mat %*% all.w)^2))

all.C<-matrix(rep(1,dim(all.pred.mat)[2]), nrow=1)
all.constraint1 <- all.w[1]>=0

all.constraint2 <- all.w[2]>=0
all.constraint3 <- all.w[3]>=0
all.constraint4 <- all.C%*%all.w == 1
all.problem <- Problem(all.objective, constraints = list(all.constraint1,  all.constraint2,  all.constraint3, all.constraint4))
all.result <- solve(all.problem)
all.wght = round(all.result$getValue(all.w),3)

######################Prediction on testing set########################
all.mod.glm.test <- lmer(all.fm.glm, data = df.train)
df_last_measure_test <- df.test %>%
  group_by(RID) %>% 
  slice_max(VISCODE) %>%  # Pick the row with the maximum VISCODE (i.e., last visit)
  ungroup()
all.ndf.test=df_last_measure_test[, c(1:p, p+2)]
all.pred.glm.test=predict(all.mod.glm.test, newdata=all.ndf.test, re.form=~(1|RID), allow.new.levels = TRUE) #for test set
all.pmse.glm.test=mean((true.y.test-all.pred.glm.test)^2)
all.mape.glm.test=mean(abs(true.y.test-all.pred.glm.test))
all.rpmse.glm.test=sqrt(all.pmse.glm.test)

all.mod.gee.test=glmgee(all.fm.gee, data=df.train, id=RID, family=gaussian(), #doesn't have method=ML or REML options as it's not likelihood based approach
                    corstr="AR-M-dependent(1)",tol = 10^-4, start= rep(0, p+1))
all.pred.gee.test=predict(all.mod.gee.test, newdata=all.ndf.test, se.fit = FALSE, varest="robust")
all.pmse.gee.test=mean((true.y.test-all.pred.gee.test)^2)
all.mape.gee.test=mean(abs(true.y.test-all.pred.gee.test))
all.rpmse.gee.test=sqrt(all.pmse.gee.test)

all.mod.rf.test<-REEMtree(all.fm.rf, data=df.train, random=~1|RID, correlation=corAR1(), method="ML")
all.newdata.test=data.frame(df_last_measure_test[, c("ADAS13",paste(colnames(df.train)[1:p]),"RID")])
all.pred.rf.test=predict(all.mod.rf.test,newdata=all.newdata.test, id=all.newdata.test$RID, EstimateRandomEffects=TRUE)
all.pmse.rf.test=mean((true.y.test-all.pred.rf.test)^2)
all.mape.rf.test=mean(abs(true.y.test-all.pred.rf.test))
all.rpmse.rf.test=sqrt(all.pmse.rf.test)

all.pred.mat.test=data.frame(all.pred.glm.test, all.pred.gee.test, all.pred.rf.test)
all.pred.el.test=as.matrix(all.pred.mat.test)%*%all.wght
all.pmse.el.test=mean((true.y.test-all.pred.el.test)^2)
all.mape.el.test=mean(abs(true.y.test-all.pred.el.test))
all.rpmse.el.test=sqrt(all.pmse.el.test)


#######################################################################
##############Adding synthetic biomarker################################
########################################################################
library(dplyr) #for data processing
library(glmmLasso) #for lasso
library(PGEE) #for penalized gee
library(glmtoolbox) #for gee
library(made4)
library(foreach) #for parallel computing
library(parallel) #for parallel computing
library(doParallel) #for parallel computing
library(Matrix)
library(lme4)
library(CVXR) #for optimization
library(REEMtree) #for random effect model
library(nlme) #for linear mixed model
library(lmmen)
library(GMMBoost) #for boosting
library(dplyr)

df=read.csv("C:/Users/.../final.data.adas13.csv")

continuous_vars <- df[,-c(1,2)]  # Select all numeric (continuous) covariates
standardized_continuous_vars <- scale(continuous_vars)

df1=cbind(standardized_continuous_vars, df[,c(1,2)])

p1=20  #20 additional noise variables
Sigma_cov <- 0.6^abs(outer(1:p1, 1:p1, "-"))

n=dim(df1)[1] # total number of observations
X <- mvrnorm(n, mu = rep(0, p1), Sigma = Sigma_cov)
colnames(X)= paste0("X", 1:p1)

df2=cbind(X, df1)

N=length(unique(df2$RID))
p=dim(df2)[2]-3  #number of covariates
B=50
threshold=0.70
set.seed(6578)
train=sample(unique(df2$RID), N*0.6, replace=FALSE)
df.train=subset(df2,RID %in% train) #training data
df.test=subset(df2,!(RID %in% train)) #testing data

##########################################################################
lam <- seq(0,40, length=20)
BIC_vec <- rep(NA, length(lam))
AIC_vec <- rep(NA, length(lam))

preds=colnames(df2)[1:p] #excluding RID, VISCODE, ADAS13
fm <- as.formula(paste("ADAS13~", paste(preds, collapse = "+")))

set.seed(200)
for (j in 1:length(lam)){
  print(paste("Iteration ", j, sep=""))
  glm1 <- glmmLasso(fm, rnd = list(RID=~1),
                    family = gaussian(link = "identity"), data = df.train, lambda=lam[j],
                    switch.NR=TRUE, final.re = FALSE, control=list(print.iter=TRUE))
  BIC_vec[j]<-glm1$bic
  AIC_vec[j]<-glm1$aic
  
}

best.lam=lam[which.min(BIC_vec)]


mat.coef.lasso1=matrix(NA, p ,B); mat.coef.lasso2=matrix(NA, p ,B)
for (j in 1:B){
  set.seed(1000+j)
  lasso.subsam=sample(c(train), length(train)*0.5, replace=FALSE)
  df.lasso.sub.sam1=subset(df.train, RID %in% lasso.subsam)
  df.lasso.sub.sam2=subset(df.train,!(RID %in% lasso.subsam))
  mod.lasso1 <- glmmLasso(fm, rnd = list(RID=~1),
                          family = gaussian(link = "identity"), data = df.lasso.sub.sam1, lambda=best.lam,
                          switch.NR=TRUE, final.re = TRUE, control=list(print.iter=TRUE))
  coef.lasso1=mod.lasso1$coefficients[-1]
  mat.coef.lasso1[,j]=coef.lasso1
  
  mod.lasso2 <- glmmLasso(fm, rnd = list(RID=~1),
                          family = gaussian(link = "identity"), data = df.lasso.sub.sam2, lambda=best.lam,
                          switch.NR=TRUE, final.re = TRUE, control=list(print.iter=TRUE))
  coef.lasso2=mod.lasso2$coefficients[-1]
  mat.coef.lasso2[,j]=coef.lasso2
} 

coef.lasso.cont=cbind(mat.coef.lasso1,mat.coef.lasso2)
coef.lasso.dis=ifelse(coef.lasso.cont !=0 ,1,0)
freq.lasso=apply(coef.lasso.dis, 1, sum)
final.set=which(freq.lasso>(2*B*threshold)) #variables selected bu LSTABEL


lasso.mod <- glmmLasso(fm, rnd = list(RID=~1),
                       family = gaussian(link="identity"), data = df.train, lambda=best.lam,
                       switch.NR=TRUE, final.re = TRUE, control=list(print.iter=TRUE))
final.set.lasso=which(lasso.mod$coefficients[-1] !=0) #variables selected bu GLMMLASSO

