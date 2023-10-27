#Set repository to install packages

#output file
OUT_file <- "/OUTPUT/FILE/DIR/"
# upload the data
variab_info<-"df_lmer_variability_4m_wetness.csv"
table_variab<-read.table(variab_info,header=TRUE,sep=",",fill=TRUE,na.strings=c(""," ","NA"))

#upload calib parameters
calib_info<-"/Users/ebruni/Desktop/DOTTORATO/SCRIPT_MODELLI/PROVE_CALIB_GLOBAL/calib_params_mm_4m.csv"
table_calib<-read.table(calib_info,header=TRUE,sep=",",fill=TRUE,na.strings=c(""," ","NA"))

#Remove Century
table_calib<- subset(table_calib, select = -c(AB_MSfrac,BE_MSfrac,Q10,T_ref))

sort.table_calib <- with(table_calib,  table_calib[order(X) , ])
all_tog<-cbind(table_variab,sort.table_calib)

#download necessary library to do lme
library(nlme)
library(MASS)

jpeg(filename=paste0(OUT_file,"calib_parameters.jpg"),
    width=1000, height=2000,res=300)

#Get param names
X_list <-colnames(table_calib)
X_list <-X_list[X_list !="X"]
X_list_expr<-c(expression("(a)  " ~ T[param] ~ " ("~degree~C~")"),
               expression("(b)  " ~ k[1] ~ " ("~{yr}^{-1}~ ")"),
               expression("(c)  " ~ k[2] ~ " ("~{yr}^{-1}~ ")"),
               "(d)   r",
               expression("(e)  " ~ k[0]~ " ("~{yr}^{-1}~ ")"))

sum_rsq<-0
sum_adj_rsq<-0
#If plotting all parameters
n_rows = dim(table_calib)[2]-1
#If plotting only 3 params
#n_rows = 3
par(mfcol = c(n_rows, 1))
par(mar=c(4.5, 5, 2, 3) + 0.1,cex.main=1.5,cex.lab=1.5,family = "sans")
it<-1
it_plot<-1
for(i in X_list){
  print("running")
  lm_outlier<-lm(get(i) ~Initial_SOC+Temp+Prec+PET+Litter_in+Clay+Carbonate+pH+Soil.C.N+Wetness_class2, data=all_tog, na.action=na.omit)
  #print(anova(lm_outlier))
  
  remove_outliers<-FALSE
  if(remove_outliers==TRUE){
    out_text<-"Removed outliers"
    #Outliers test
    cooksd <- cooks.distance(lm_outlier)
    all_tog1<-all_tog
    #Eliminate outlier
    all_tog1<-all_tog1[c(cooksd<4*mean(cooksd, na.rm=T)),]
    #Redo lm without outlier
    lm<-lm(get(i) ~Initial_SOC+Temp+Prec+PET+Litter_in+Clay+Carbonate+pH+Soil.C.N+Wetness_class2, data=all_tog1, na.action=na.omit)
    
  }else{
    out_text<-"With outliers"
    lm<-lm_outlier
    remove(lm_outlier)
    all_tog1<-all_tog
  }

  #Do the stepwise selection
  #BIC (Bayesian)
  #model_tot_step <-stepAIC(lm, direction=c("both"),trace = 0, k=log(nrow(all_tog)))
  #AIC (Akaike)->penalizes less the more complex models
  model_tot_step <-stepAIC(lm, direction=c("both"),trace = 0, k=2) 
  summary_aic<-summary(model_tot_step)

  #Estimated coefficients
  list_coeff <- summary_aic$coefficients[,"Estimate"]
  assign(paste0("name",it), data.frame(list_coeff))
  
  print("####################")
  print("Parameter predicted:")
  print(i)
  print("####################")
  print(summary_aic)
  r_sq<-summary_aic$r.squared
  adj_r_sq<-summary_aic$adj.r.squared
  sum_rsq<-sum_rsq+r_sq
  sum_adj_rsq<-sum_adj_rsq+adj_r_sq
  print("Multiple R-squared")
  print(r_sq)
  print("Adjusted R-squared")
  print(adj_r_sq)
  
  #To see how standard deviation of coefficients are calculated
  #https://stat.ethz.ch/pipermail/r-help/2006-September/113115.html
  
  table_predict<-cbind(predict(model_tot_step),all_tog1[[i]])
  #table_predict<-cbind(predict(lm),all_tog1[[i]])
  #plot_pred<-TRUE
  #if(plot_pred==TRUE){
  #all params
  if(i=="T_param" | i=="k1" | i=='k2' | i=='r' | i=='k0'){
    print("is plotting")
  #3params only
  #if(i=="T_param" | i=="r" | i=='k0'){  
    #Plot (save inches: (2;9.2) if 5 params)
    #Plot (save inches: (3;7.2) if 3 params)
    plot(table_predict[,2],table_predict[,1],xlim=c(min(table_predict),max(table_predict)),
         ylim=c(min(table_predict),max(table_predict)),main=X_list_expr[it_plot],adj=0,xlab=" ",ylab=" ",pch=19)
    abline(a=0,b=1)
    
    text_pos<-mean(c(min(table_predict),max(table_predict)))
    #text(text_pos-0.2*(max(table_predict)-min(table_predict)),text_pos+0.1*(max(table_predict)-min(table_predict)),paste(expression(R^2 ~ = ~),round(r_sq,digits=2))))
    text(text_pos-0.35*(max(table_predict)-min(table_predict)),text_pos+0.21*(max(table_predict)-min(table_predict)),expression(R^2),cex=1.2)
    text(text_pos-0.20*(max(table_predict)-min(table_predict)),text_pos+0.2*(max(table_predict)-min(table_predict)),paste0(" = ",round(r_sq,digits=2)),cex=1.2)

    if(it_plot==5){
      title(xlab="Parameters calibrated on the LTEs")
    }
    #if(it_plot==3){
    #  title(ylab="Predicted parameter")
    #}
    it_plot=it_plot+1
  }
  it<-it+1
}
#print(table_predict)

mtext("Statistically calibrated parameters",side=2,line=-2,outer=TRUE,cex=1.)

dev.off()

mean_rsq<-sum_rsq/length(X_list)
mean_adj_rsq<-sum_adj_rsq/length(X_list)

print(mean_rsq)
print(mean_adj_rsq)

df_list <-list("T_param"=name1,"k1"=name2,"k2"=name3,"r"=name4,
                            "k0"=name5) 
#Save functions
loc_outputfiles<-'/Users/ebruni/Desktop/DOTTORATO/SCRIPT_MODELLI/PROVE_CALIB_GLOBAL/CALIB_PARAM/V8_wetness2/'
write.table( data.frame(name1), paste0(loc_outputfiles,'T_param_calib_func.csv'),sep=',')
write.table( data.frame(name2), paste0(loc_outputfiles,'k1_calib_func.csv'),sep=',')
write.table( data.frame(name3), paste0(loc_outputfiles,'k2_calib_func.csv'),sep=',')
write.table( data.frame(name4), paste0(loc_outputfiles,'r_calib_func.csv'),sep=',')
write.table( data.frame(name5), paste0(loc_outputfiles,'k0_calib_func.csv'),sep=',')
#write.table( data.frame(name6), paste0(loc_outputfiles,'AB_MSfrac_calib_func.csv'),sep=',')
#write.table( data.frame(name7), paste0(loc_outputfiles,'BE_MSfrac_calib_func.csv'),sep=',')
#write.table( data.frame(name8), paste0(loc_outputfiles,'Q10_calib_func.csv'),sep=',')
#write.table( data.frame(name9), paste0(loc_outputfiles,'T_ref_calib_func.csv'),sep=',')
#write.table( data.frame(name10), paste0(loc_outputfiles,'av_calib_func.csv'),sep=',')
#write.table( data.frame(name11), paste0(loc_outputfiles,'ak_calib_func.csv'),sep=',')
#write.table( data.frame(name12), paste0(loc_outputfiles,'fMET_calib_func.csv'),sep=',')
#write.table( data.frame(name13), paste0(loc_outputfiles,'eact_lb_calib_func.csv'),sep=',')
#write.table( data.frame(name14), paste0(loc_outputfiles,'eact_pl_calib_func.csv'),sep=',')
#write.table( data.frame(name15), paste0(loc_outputfiles,'kaff_pl_calib_func.csv'),sep=',')


library(caret)
library(dplyr)

###############################
#Cross validation with split sampling
################################

MPESS_df <- data.frame(matrix(ncol = 2, nrow = length(X_list)))
colnames(MPESS_df) <- c('Param', 'MPESS')
j<-1
for(param_i in X_list){
  set.seed(123) 
  print("####################")
  print(param_i)
  print("####################")
  
  #splits data into p/1-p, X times
  split_index <- createDataPartition(all_tog[,param_i], p = 0.7, list = FALSE,times=100) 
  #print(split_index)
  #Initialize output dataframe
  error_df <- data.frame(matrix(ncol = 2, nrow = ncol(split_index)))
  colnames(error_df) <- c(paste0('test_error',as.character(param_i)), 'fold')
  
  for(i in 1:nrow(error_df)){
    #i=1
    # use ith column of split_index to create feature and target training/test sets
    features_train <- all_tog[ split_index[,i], !(names(all_tog) %in% c('X'))] 
    features_test  <- all_tog[-split_index[,i], !(names(all_tog) %in% c('X'))]
    target_train <- all_tog[ split_index[,i], as.character(param_i)]
    target_test <- all_tog[-split_index[,i], as.character(param_i)]
    
    # Fit the model and predict
    lm_fit <-lm(get(param_i) ~Initial_SOC+Temp+Prec+PET+Litter_in+Clay+Carbonate+pH+Soil.C.N+Wetness_class2, data=features_train, na.action=na.omit)
    model_fit_step <-stepAIC(lm_fit, direction=c("both"),trace = 0, k=2) 
    summary_aic_fit<-summary(model_fit_step)
    lm_pred <- predict(model_fit_step, features_test, type = 'response' )
    
    # Calculate error and store it
    #error <- mean(ifelse(target_test != lm_pred, 1, 0))
    error <- mean(abs((lm_pred-target_test)/target_test))
    #print("Predicted")
    #print(lm_pred)
    #print("Real")
    #print(target_test)
    error_df[i,paste0('test_error',as.character(param_i))] <- error
    error_df[i, 'fold'] <- i
    
  }
  print(error_df)
  mean_error_pred = mean(error_df$test_error)
  print(paste("Mean error over 100 tests for param",param_i))
  print(mean_error_pred)

  MPESS_df[j,"Param"]<-param_i
  MPESS_df[j,"MPESS"]<-mean_error_pred
  j<-j+1
}


rownames(MPESS_df)<-X_list
print(MPESS_df)

###############################
#Leave one out cross validation
################################
#RothC
ctrl <- trainControl(method = "LOOCV")
#fit a regression model and use LOOCV to evaluate performance
model_rothc <- train(T_param~Initial_SOC+Temp+Prec+PET+Litter_in+Clay+Carbonate+pH+Soil.C.N+Wetness_class2, data = all_tog, method = "lmStepAIC", trControl = ctrl)
#summary(model_rothc)
print(model_rothc)

#AMG
#fit a regression model and use LOOCV to evaluate performance
model_amg <- train(k0~Initial_SOC+Temp+Prec+PET+Litter_in+Clay+Carbonate+pH+Soil.C.N+Wetness_class2, data = all_tog, method = "lmStepAIC", trControl = ctrl)
print(model_amg)

#ICBM
#fit a regression model and use LOOCV to evaluate performance
model_ICBM_k1 <- train(k1~Initial_SOC+Temp+Prec+PET+Litter_in+Clay+Carbonate+pH+Soil.C.N+Wetness_class2, data = all_tog, method = "lmStepAIC", trControl = ctrl)
print(model_ICBM_k1)
model_ICBM_k2 <- train(k2~Initial_SOC+Temp+Prec+PET+Litter_in+Clay+Carbonate+pH+Soil.C.N+Wetness_class2, data = all_tog, method = "lmStepAIC", trControl = ctrl)
print(model_ICBM_k2)
model_ICBM_r <- train(r~Initial_SOC+Temp+Prec+PET+Litter_in+Clay+Carbonate+pH+Soil.C.N+Wetness_class2, data = all_tog, method = "lmStepAIC", trControl = ctrl)
print(model_ICBM_r)

err_res<-c("RMSE","Rsquared","MAE")
df_LOOCV <- rbind(model_rothc$results[err_res],model_amg$results[err_res],
                  model_ICBM_k1$results[err_res],model_ICBM_k2$results[err_res],model_ICBM_r$results[err_res])
row_ord <-c("T_param","k0","k1","k2","r")
row.names(df_LOOCV) <- row_ord

means_param <- c(mean(all_tog[[row_ord[1]]]),mean(all_tog[[row_ord[2]]]),mean(all_tog[[row_ord[3]]]),mean(all_tog[[row_ord[4]]]),mean(all_tog[[row_ord[5]]]))

#Relative root mean squared error
df_LOOCV$rRMSE<-df_LOOCV$RMSE/means_param

#If you want to add the cross val with split sampling
table_fin <- merge(df_LOOCV,MPESS_df,by='row.names', all = TRUE)
#Save
write.table(table_fin,file = paste0(OUT_file,"wetness2_table_LOOCV.csv"),sep=",", row.names = F)




            
