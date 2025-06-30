library(mvtnorm)
library(glmmLasso)
library(PGEE)
library(lme4)
library(pROC)
library(glmtoolbox) #for gee
library(made4)
library(foreach)
library(parallel)
library(doParallel)
library(Matrix)
library(nlme) #for linear mixed model
library(REEMtree) #for random forest
library(GMMBoost) #for boosting
library(CVXR) #for optimization

#N, sample size
#p, number of covariates
#rho=0 if X is independent, rho=0.5 if X is correlated
#threshold=0.70 (fixed), threshold probability for stability selection
#iter=100 (fixed), number of iteration
#B=50 (fixed), 50*2, total 100 subsamples
#train.pcntg=0.6 (fixed), percentage of data in the training set
#best.lam, tuning parameter for traditional LASSO

#sigma_u=1.2 and 0.8
#Based on the AIC/BIC values following values are obtained for lambda which is denoted by best.lam
#N=600, p=200, rho_cv=0, rho_rm=0.7, sigma2_u=1.2, best.lam=130
#N=600, p=800, rho_cv=0, rho_rm=0.7, sigma2_u=1.2, best.lam=190
#N=1000, p=500, rho_cv=0, rho_rm=0.7, sigma2_u=1.2, best.lam=210
#N=1000, p=1500, rho_cv=0, rho_rm=0.7, sigma2_u=1.2, best.lam=250
#N=600, p=200, rho_cv=0, rho_rm=0.7, sigma2_u=0.8, best.lam=130
#N=600, p=800, rho_cv=0, rho_rm=0.7, sigma2_u=0.8, best.lam=190
#N=1000, p=500, rho_cv=0, rho_rm=0.7, sigma2_u=0.8, best.lam=210
#N=1000, p=1500, rho_cv=0, rho_rm=0.7, sigma2_u=0.8, best.lam=250

#N=600, p=200, rho_cv=0.6, rho_rm=0.7, sigma2_u=1.2, best.lam=120
#N=600, p=800, rho_cv=0.6, rho_rm=0.7, sigma2_u=1.2, best.lam=120
#N=1000, p=500, rho_cv=0.6, rho_rm=0.7, sigma2_u=1.2, best.lam=200
#N=1000, p=1500, rho_cv=0.6, rho_rm=0.7, sigma2_u=1.2, best.lam=220
#N=600, p=200, rho_cv=0.6, rho_rm=0.7, sigma2_u=0.8, best.lam=120
#N=600, p=800, rho_cv=0.6, rho_rm=0.7, sigma2_u=0.8, best.lam=120
#N=1000, p=500, rho_cv=0.6, rho_rm=0.7, sigma2_u=0.8, best.lam=200
#N=1000, p=1500, rho_cv=0.6, rho_rm=0.7, sigma2_u=0.8, best.lam=220


t1=Sys.time()
glmm_wrapper=function(N, p, m, rho_cov, rho_rm, sigma2_u, threshold, iter, B, train.pcntg, best.lam){
  ########################################Cross-validation#################################
  id.vect <- rep(1:N, each = m)
  
  cl <- parallel::makeCluster(detectCores()-1, setup_strategy = "sequential")
  # Activate cluster for foreach library
  registerDoParallel(cl)
  
  unregister_dopar <- function() {
    env <- foreach:::.foreachGlobals
    rm(list=ls(name=env), pos=env)
  }
  
  stabel_iter=foreach(i = 1:iter, .packages = c("glmmLasso", "PGEE", "mvtnorm", "lme4", "glmtoolbox", "made4", "nlme", "REEMtree", "CVXR")) %dopar%{
    set.seed(12345*i)
    
    Sigma_cov <- rho_cov^abs(outer(1:p, 1:p, "-")) # AR(1) parameter for covariates
    
    # Generate a matrix for all subjects at once, with each row corresponding to a covariate measurement
    X <- mvrnorm(N * m, mu = rep(0, p), Sigma = Sigma_cov)
    
    b = matrix(0,p,1)
    b[1:8] = c(-1,1,-1,1,-1,1,-1,1)#fixed effect, p*1
    Z=matrix(0, nrow=N*m, ncol=N) #n*Nq=1800*600, only random intercept is considered
    for (col_index in seq_len(N)) {
      # Calculate the row range for each column based on m
      start_row <- (col_index - 1) * m + 1
      end_row <- col_index * m
      row_range <- seq(start_row, end_row)
      
      Z[row_range, col_index] <- 1
    }
    u=rnorm(N, 0, sqrt(sigma2_u)) # random effect, which was denoted by b in the paper
    
    
    # covariance matrix of error
    R=rho_rm^abs(outer(1:m, 1:m, "-")) 
    e = rmvnorm(N, mean = rep(0,m), R)
    
    y.vect=X%*%b+Z%*%u+as.vector(t(e))
    
    df = data.frame(X, id.vect,y.vect)
    colnames(df)=c(paste("X",1:length(b),sep=""), "id","y")
    
    train=sample(1:N, N*0.6, replace=FALSE) #60% subjects in the training data
    
    df.train=subset(df,id %in% train) #trainning data frame
    df.test=subset(df,!(id %in% train)) #testing data frame
    
    fm <- as.formula(paste("y ~", paste0("X", 1:p, collapse = "+"))) #formula using all the variables
    mat.coef.lasso1=matrix(NA,p,B); mat.coef.lasso2=matrix(NA,p,B)
    for (j in 1:B){
      set.seed(2000+j)
      lasso.subsam=sample(c(train), length(train)*0.5, replace=FALSE) #randomly taking 50% subjects from the training set
      df.lasso.sub.sam1=subset(df.train,id %in% lasso.subsam)
      df.lasso.sub.sam2=subset(df.train,!(id %in% lasso.subsam))
      mod.lasso1 <- glmmLasso(fm, rnd = list(id=~1),
                              family = gaussian(link="identity"), data = df.lasso.sub.sam1, lambda=best.lam,
                              switch.NR=TRUE, final.re = TRUE, control=list(print.iter=TRUE))
      coef.lasso1=mod.lasso1$coefficients[-1]
      mat.coef.lasso1[,j]=coef.lasso1
      
      mod.lasso2 <- glmmLasso(fm, rnd = list(id=~1),
                              family = gaussian(link="identity"), data = df.lasso.sub.sam2, lambda=best.lam,
                              switch.NR=TRUE, final.re = TRUE, control=list(print.iter=TRUE))
      coef.lasso2=mod.lasso2$coefficients[-1]
      mat.coef.lasso2[,j]=coef.lasso2
    } 
    
    coef.lasso.cont=cbind(mat.coef.lasso1,mat.coef.lasso2)
    coef.lasso.dis=ifelse(coef.lasso.cont !=0 ,1,0)
    freq.lasso=apply(coef.lasso.dis, 1, sum)
    final.set= which(freq.lasso > threshold*2*B) #final set of selected variables by STABEL
    
    ###False positive and false negative (by LSTABEL)
    true.b=which(b != 0)
    tp.stabel = length(intersect(final.set, true.b)) #True Positives
    fp.stabel = length(setdiff(final.set, true.b)) #False Positives
    fn.stabel = length(setdiff(true.b, final.set)) #False Negatives
    tn.stabel = length(setdiff(1:p, union(true.b, final.set))) #True negatives
    
    TPR.stabel= tp.stabel/(tp.stabel+fn.stabel) #proportion of true covariates selected in the final model.
    FDR.stabel= fp.stabel/(fp.stabel+tp.stabel) #the number of noise variables selected as a proportion of the total number of variables selected in the final model.
    FDER.stabel=fp.stabel/(fp.stabel+tn.stabel) #the number of noise variables selected in the final model as a proportion of the total noise variables available for selection
    
    ########################building ensemble learning through cross-validation using the selected variables (final.set)#############################
    set.seed(1234)
    nfold=5 #5 fold cross varlidation for longitudinal super learner
    subject_ids <- unique(df.train$id) # Create a subject-to-subset mapping
    n_subjects <- length(subject_ids)
    subset_assignment <- rep(1:nfold, length.out = n_subjects)
    
    subset_assignment <- sample(subset_assignment) # Shuffle the assignment to randomize
    subject_to_subset <- data.frame(id = subject_ids, subset = subset_assignment) # Create a data frame to map subject_id to subset
    df.train.cv <- merge(df.train, subject_to_subset, by = "id") # Merge this mapping with the training data
    covariate_cols <- grep("^X", names(df.train.cv), value = TRUE)
    new_order <- c(covariate_cols, "id", "y", "subset") # Define the new column order
    df.train_reorganized <- df.train.cv[, new_order]
    
    #Create subsets based on the new subset variable
    subset1 <- subset(df.train_reorganized, subset == 1)[, -(p+3)]
    subset2 <- subset(df.train_reorganized, subset == 2)[, -(p+3)]
    subset3 <- subset(df.train_reorganized, subset == 3)[, -(p+3)]
    subset4 <- subset(df.train_reorganized, subset == 4)[, -(p+3)]
    subset5 <- subset(df.train_reorganized, subset == 5)[, -(p+3)]
    
    fold1=unique(subset1$id) #id belongs to fold 1
    fold2=unique(subset2$id)
    fold3=unique(subset3$id)
    fold4=unique(subset4$id)
    fold5=unique(subset5$id)
    
    all.fold <- list(fold1,fold2,fold3,fold4,fold5)
    combinations <- combn(all.fold, (nfold-1), simplify = FALSE) # Generate all combinations of 4 vectors
    combined.index <- lapply(combinations, function(comb) do.call(c, comb)) # Perform the combination operation (e.g., concatenate)
    
    
    #####################################################
    ##################GLMM for ensemble learning###########
    
    fm.glm <- as.formula(paste("y ~", paste0("X", final.set, collapse = " + "), "+ (1 | id)"))
    pred.glm.vec=vector("numeric")
    for (l in 1:nfold){
      mod.glm <- lmer(fm.glm, data = subset(df.train, id %in% c(combined.index[[l]])))
      ndf.all.rm=subset(df.train, !(id %in% combined.index[[l]])) #testing sets with all repeated measures
      ndf=ndf.all.rm[seq(m, dim(ndf.all.rm)[1], by=m), c(final.set, p+1)]
      pred.glm=predict(mod.glm, newdata=ndf, re.form=~(1|id), allow.new.levels = TRUE) #for test set
      pred.glm.vec=c(pred.glm.vec, pred.glm)
    }
    
    
    #################################################################
    ##generalized estimating equations for ensebmle learning#########
    
    fm.gee <- as.formula(paste("y ~", paste(colnames(df.train)[final.set], collapse = " + ")))
    pred.gee.vec=vector("numeric")
    for (l in 1:nfold){
      mod.gee=glmgee(fm.gee, data= subset(df.train, id %in% c(combined.index[[l]])), id=id, family=gaussian(), #doesn't have method=ML or REML options as it's not likelihood based approach
                     corstr="AR-M-dependent(1)",tol = 10^-4, start= rep(0, (length(final.set)+1)))
      ndf.all.rm=subset(df.train, !(id %in% combined.index[[l]])) #testing sets with all repeated measures
      ndf=ndf.all.rm[seq(m, dim(ndf.all.rm)[1], by=m), c(paste(colnames(df.train)[final.set]),"id")]
      pred.gee=predict(mod.gee, newdata=ndf, se.fit = FALSE, varest="robust")
      pred.gee.vec=c(pred.gee.vec, pred.gee)
    }  
    
    # #######################Randomforest for ensemble learning##############################
    fm.rf <- as.formula(paste("y ~", paste(colnames(df.train)[final.set], collapse = " + ")))
    pred.rf.vec=vector("numeric")
    true.y.vec=c()  #keeping record of true y for cross-validation only
    for (l in 1:nfold){
      mod.rf=REEMtree(fm.rf, data=subset(df.train, id %in% c(combined.index[[l]])), random=~1|id, correlation=corAR1(), method="ML")
      ndf.all.rm=subset(df.train, !(id %in% combined.index[[l]])) #testing sets with all repeated measures
      ndf.rf=ndf.all.rm[seq(m, dim(ndf.all.rm)[1], by=m), c("y", paste(colnames(df.train)[final.set]),"id")] #the true y (for superlearner) is recorded from this new dataframe(cross-validation)
      true.y=ndf.rf$y
      true.y.vec=c(true.y.vec, true.y)
      pred.rf=predict(mod.rf,newdata=ndf.rf, id=ndf.rf$id, EstimateRandomEffects=TRUE)
      pred.rf.vec=c(pred.rf.vec, pred.rf)
    }
    
    ###########################################################################
    #####################Ensemble learning####################################
    pred.mat=data.frame(pred.glm.vec, pred.gee.vec, pred.rf.vec)
    names(pred.mat)=c("GLMM", "GEE", "RF")
    w <- Variable(dim(pred.mat)[2])
    objective <- Minimize(sum((true.y.vec - pred.mat %*% w)^2))
    
    C<-matrix(rep(1,dim(pred.mat)[2]), nrow=1)
    constraint1 <- w[1]>=0
    
    constraint2 <- w[2]>=0
    constraint3 <- w[3]>=0
    constraint4 <- C%*%w == 1
    problem <- Problem(objective, constraints = list(constraint1,  constraint2,  constraint3, constraint4))
    result <- solve(problem)
    wght = round(result$getValue(w),3) #estimated weight for ensemble learning
    
    ######################Prediction on testing set########################
    fm.glm <- as.formula(paste("y ~", paste0("X", final.set, collapse = " + "), "+ (1 | id)"))
    mod.glm.test <- lmer(fm.glm, data = df.train, REML=FALSE)
    ndf.test=df.test[seq(m, dim(df.test)[1], by=m),c(final.set, p+1)] #p+1 is for id
    pred.glm.test=predict(mod.glm.test, newdata=ndf.test, re.form=~(1|id), allow.new.levels = TRUE) #for test set
    pmse.glm.test=mean((df.test$y[seq(m, dim(df.test)[1], by=m)]-pred.glm.test)^2) #PMSE from GLMM
    mape.glm.test=mean(abs(df.test$y[seq(m, dim(df.test)[1], by=m)]-pred.glm.test)) # MAPE from GLMM
    rpmse.glm.test=sqrt(pmse.glm.test) #RPMSE from GLMM
    
    fm.gee <- as.formula(paste("y ~", paste(colnames(df.train)[final.set], collapse = " + ")))
    mod.gee.test=glmgee(fm.gee, data=df.train, id=id, family=gaussian(), #doesn't have method=ML or REML options as it's not likelihood based approach
                        corstr="AR-M-dependent(1)",tol = 10^-4, start= rep(0, (length(final.set)+1)))
    pred.gee.test=predict(mod.gee.test, newdata=df.test[seq(m, dim(df.test)[1], by=m), c(paste(colnames(df.train)[final.set]),"id")], se.fit = FALSE, varest="robust")
    pmse.gee.test=mean((df.test$y[seq(m, dim(df.test)[1], by=m)]-pred.gee.test)^2) #PMSE from GEE
    mape.gee.test=mean(abs(df.test$y[seq(m, dim(df.test)[1], by=m)]-pred.gee.test)) # MAPE from GEE
    rpmse.gee.test=sqrt(pmse.gee.test) #RPMSE from GEE
    
    fm.rf<- as.formula(paste("y ~", paste(colnames(df.train)[final.set], collapse = " + ")))
    mod.rf.test<-REEMtree(fm.rf, data=df.train, random=~1|id, correlation=corAR1(), method="ML", ErrorTolerance=0.001, tree.control=rpart.control(cp=0.001))
    newdata.test=df.test[seq(m, dim(df.test)[1], by=m), c("y",paste(colnames(df.train)[final.set]),"id")]
    pred.rf.test=predict(mod.rf.test,newdata=newdata.test, id=newdata.test$id, EstimateRandomEffects=TRUE)
    pmse.rf.test=mean((df.test$y[seq(m, dim(df.test)[1], by=m)]-pred.rf.test)^2) #PMSE from RF
    mape.rf.test=mean(abs(df.test$y[seq(m, dim(df.test)[1], by=m)]-pred.rf.test)) # MAPE from RF
    rpmse.rf.test=sqrt(pmse.rf.test) #RPMSE from RF
    
    pred.mat.test=data.frame(pred.glm.test, pred.gee.test, pred.rf.test)
    pred.el.test=as.matrix(pred.mat.test)%*%wght
    pmse.el.test=mean((df.test$y[seq(m, dim(df.test)[1], by=m)]-pred.el.test)^2) #PMSE from ensemble learning
    mape.el.test=mean(abs(df.test$y[seq(m, dim(df.test)[1], by=m)]-pred.el.test)) # MAPE from ensemble learning
    rpmse.el.test=sqrt(pmse.el.test) #RPMSE from ensemble learning
    
    #####################################################################################
    #####################GLMMLASSO#######################################################
    #####################################################################################
    lasso.mod <- glmmLasso(fm, rnd = list(id=~1),
                           family = gaussian(link="identity"), data = df.train, lambda=best.lam,
                           switch.NR=TRUE, final.re = TRUE, control=list(print.iter=TRUE))
    final.set.lasso=which(lasso.mod$coefficients[-1] !=0) #selected variables by GLMMLASSO
    
    pred.test.lasso=predict(lasso.mod, df.test[seq(m, dim(df.test)[1], by=m),]) 
    pmse.test.lasso=mean((df.test$y[seq(m, dim(df.test)[1], by=m)]-pred.test.lasso)^2) #PMSE from GLMMLASSO
    mape.test.lasso=mean(abs(df.test$y[seq(m, dim(df.test)[1], by=m)]-pred.test.lasso))# MAPE from GLMMLASSO
    rpmse.test.lasso=sqrt(pmse.test.lasso) #RPMSE from GLMMLASSO
    
    ###########################################################################################
    #####################GLMMLASSO with ensemble learning######################################
    ###########################################################################################
    
    ##################GLMM###########
    fm.glm.lasso <- as.formula(paste("y ~", paste0("X", final.set.lasso, collapse = " + "), "+ (1 | id)"))
    pred.glm.vec.lasso=vector("numeric")
    for (l in 1:nfold){
      mod.glm.lasso <- lmer(fm.glm.lasso, data = subset(df.train, id %in% c(combined.index[[l]])))
      ndf.all.rm.lasso=subset(df.train, !(id %in% combined.index[[l]])) #testing sets with all repeated measures
      ndf.lasso=ndf.all.rm.lasso[seq(m, dim(ndf.all.rm.lasso)[1], by=m), c(final.set.lasso, p+1)]
      pred.glm.lasso=predict(mod.glm.lasso, newdata=ndf.lasso, re.form=~(1|id), allow.new.levels = TRUE) #for test set
      pred.glm.vec.lasso=c(pred.glm.vec.lasso, pred.glm.lasso)
    }
    
    
    ##generalized estimating equations
    fm.gee.lasso <- as.formula(paste("y ~", paste(colnames(df.train)[final.set.lasso], collapse = " + ")))
    pred.gee.vec.lasso=vector("numeric")
    for (l in 1:nfold){
      mod.gee.lasso=glmgee(fm.gee.lasso, data= subset(df.train, id %in% c(combined.index[[l]])), id=id, family=gaussian(), #doesn't have method=ML or REML options as it's not likelihood based approach
                           corstr="AR-M-dependent(1)",tol = 10^-4, start= rep(0, (length(final.set.lasso)+1)))
      ndf.all.rm.lasso=subset(df.train, !(id %in% combined.index[[l]])) #testing sets with all repeated measures
      ndf.lasso=ndf.all.rm.lasso[seq(m, dim(ndf.all.rm.lasso)[1], by=m), c(paste(colnames(df.train)[final.set.lasso]),"id")]
      pred.gee.lasso=predict(mod.gee.lasso, newdata=ndf.lasso, se.fit = FALSE, varest="robust")
      pred.gee.vec.lasso=c(pred.gee.vec.lasso, pred.gee.lasso)
    }  
    
    ########################Randomforest##############################
    fm.rf.lasso <- as.formula(paste("y ~", paste(colnames(df.train)[final.set.lasso], collapse = " + ")))
    pred.rf.vec.lasso=vector("numeric")
    true.y.vec=c()  #keeping record of true y for cross-validation only
    for (l in 1:nfold){
      mod.rf.lasso=REEMtree(fm.rf.lasso, data=subset(df.train, id %in% c(combined.index[[l]])), random=~1|id, correlation=corAR1(), method="ML")
      ndf.all.rm.lasso=subset(df.train, !(id %in% combined.index[[l]])) #testing sets with all repeated measures
      ndf.rf.lasso=ndf.all.rm.lasso[seq(m, dim(ndf.all.rm.lasso)[1], by=m), c("y", paste(colnames(df.train)[final.set.lasso]),"id")] #the true y (for superlearner) is recorded from this new dataframe(cross-validation)
      true.y=ndf.rf.lasso$y
      true.y.vec=c(true.y.vec, true.y)
      pred.rf.lasso=predict(mod.rf.lasso,newdata=ndf.rf.lasso, id=ndf.rf.lasso$id, EstimateRandomEffects=TRUE)
      pred.rf.vec.lasso=c(pred.rf.vec.lasso, pred.rf.lasso)
    }
    
    ###########################################################################
    ##########Ensemble learning-variables are selected by GLMMLASSO############
    
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
    wght.lasso = round(result.lasso$getValue(w.lasso),3) #estimated weight for ensemble learning
    
    ######################Prediction on testing set########################
    fm.glm <- as.formula(paste("y ~", paste0("X", final.set.lasso, collapse = " + "), "+ (1 | id)"))
    mod.glm.test.lasso <- lmer(fm.glm.lasso, data = df.train)
    ndf.test.lasso=df.test[seq(m, dim(df.test)[1], by=m),c(final.set.lasso, p+1)] #p+1 is for id
    pred.glm.test.lasso=predict(mod.glm.test.lasso, newdata=ndf.test.lasso, re.form=~(1|id), allow.new.levels = TRUE) #for test set
    pmse.glm.test.lasso=mean((df.test$y[seq(m, dim(df.test)[1], by=m)]-pred.glm.test.lasso)^2) #PMSE from GLMM
    mape.glm.test.lasso=mean(abs(df.test$y[seq(m, dim(df.test)[1], by=m)]-pred.glm.test.lasso)) #MAPE from GLMM
    rpmse.glm.test.lasso=sqrt(pmse.glm.test.lasso) #RPMSE from GLMM
    
    fm.gee <- as.formula(paste("y ~", paste(colnames(df.train)[final.set.lasso], collapse = " + ")))
    mod.gee.test.lasso=glmgee(fm.gee.lasso, data=df.train, id=id, family=gaussian(), #doesn't have method=ML or REML options as it's not likelihood based approach
                              corstr="AR-M-dependent(1)",tol = 10^-4, start= rep(0, (length(final.set.lasso)+1)))
    pred.gee.test.lasso=predict(mod.gee.test.lasso, newdata=df.test[seq(m, dim(df.test)[1], by=m), c(paste(colnames(df.train)[final.set.lasso]),"id")], se.fit = FALSE, varest="robust")
    pmse.gee.test.lasso=mean((df.test$y[seq(m, dim(df.test)[1], by=m)]-pred.gee.test.lasso)^2) #PMSE from GEE
    mape.gee.test.lasso=mean(abs(df.test$y[seq(m, dim(df.test)[1], by=m)]-pred.gee.test.lasso)) #MAPE from GEE
    rpmse.gee.test.lasso=sqrt(pmse.gee.test.lasso) #RPMSE from GEE
    
    fm.rf<- as.formula(paste("y ~", paste(colnames(df.train)[final.set.lasso], collapse = " + ")))
    mod.rf.test.lasso<-REEMtree(fm.rf.lasso, data=df.train, random=~1|id, correlation=corAR1(), method="ML")
    newdata.test.lasso=df.test[seq(m, dim(df.test)[1], by=m), c("y",paste(colnames(df.train)[final.set.lasso]),"id")]
    pred.rf.test.lasso=predict(mod.rf.test.lasso, newdata=newdata.test.lasso, id=newdata.test.lasso$id, EstimateRandomEffects=TRUE)
    pmse.rf.test.lasso=mean((df.test$y[seq(m, dim(df.test)[1], by=m)]-pred.rf.test.lasso)^2) #PMSE from RF
    mape.rf.test.lasso=mean(abs(df.test$y[seq(m, dim(df.test)[1], by=m)]-pred.rf.test.lasso)) #MAPE from RF
    rpmse.rf.test.lasso=sqrt(pmse.rf.test.lasso) #RPMSE from RF
    
    pred.mat.test.lasso=data.frame(pred.glm.test.lasso, pred.gee.test.lasso, pred.rf.test.lasso)
    pred.el.test.lasso=as.matrix(pred.mat.test.lasso)%*%wght.lasso
    pmse.el.test.lasso=mean((df.test$y[seq(m, dim(df.test)[1], by=m)]-pred.el.test.lasso)^2) #PMSE from ensemble learning
    mape.el.test.lasso=mean(abs(df.test$y[seq(m, dim(df.test)[1], by=m)]-pred.el.test.lasso)) #MAPE from ensemble learning
    rpmse.el.test.lasso=sqrt(pmse.el.test.lasso) #RPMSE from ensemble learning
    
    ###False positive and false negative (by GLMMLASSO)
    tp.lasso = length(intersect(final.set.lasso, true.b)) #True Positives
    fp.lasso = length(setdiff(final.set.lasso, true.b)) #False Positives
    fn.lasso = length(setdiff(true.b, final.set.lasso)) #False Negatives
    tn.lasso = length(setdiff(1:p, union(true.b, final.set.lasso))) #True negatives
    
    TPR.lasso = tp.lasso /(tp.lasso +fn.lasso ) #proportion of true covariates selected in the final model.
    FDR.lasso = fp.lasso /(fp.lasso +tp.lasso ) #the number of noise variables selected as a proportion of the total number of variables selected in the final model.
    FDER.lasso =fp.lasso /(fp.lasso +tn.lasso ) #the number of noise variables selected in the final model as a proportion of the total noise variables available for selection
    
    #####combining all results
    length.stabel=length(final.set)
    length.lasso=length(final.set.lasso)
    
    c(length.stabel, TPR.stabel, FDR.stabel, FDER.stabel, 
      pmse.glm.test, mape.glm.test, rpmse.glm.test, 
      pmse.gee.test, mape.gee.test, rpmse.gee.test, 
      pmse.rf.test, mape.rf.test, rpmse.rf.test, 
      pmse.el.test, mape.el.test, rpmse.el.test,
      length.lasso,TPR.lasso, FDR.lasso, FDER.lasso, 
      pmse.test.lasso, mape.test.lasso, rpmse.test.lasso,
      pmse.el.test.lasso, mape.el.test.lasso, rpmse.el.test.lasso)
    
  }
  
  pred.sum=as.data.frame(do.call(rbind, stabel_iter)) #making a dataframe using all results
  
}

result=glmm_wrapper(N=600, p=200, m=3, rho_cov=0, rho_rm=0.7, sigma2_u=1.2, threshold=0.70, iter=10, B=50, train.pcntg=0.6, best.lam=130)
Sys.time()-t1


colnames(result) <- c("length.stabel","TPR.stabel", "FDR.stabel", "FDER.stabel", "pmse.glm", "mape.glm", "rpmse.glm", "pmse.gee", "mape.gee", "rpmse.gee", "pmse.rf", "mape.rf", "rpmse.rf", "pmse.el", "mape.el", "rpmse.el",
                      "length.lasso","TPR.lasso", "FDR.lasso", "FDER.lasso", "pmse.lasso", "mape.lasso", "rpmse.lasso", "pmse.lasso.el", "mape.lasso.el", "rpmse.lasso.el")
result1=as.data.frame(lapply(result, as.numeric))
apply(result1, 2, mean)

result2=as.data.frame(lapply(result1[,c(-1,-17)], as.numeric))
