###Note this is mainly a play program, used for testing and various views and sortings

library(pROC)
library(ROCR)
##pROC
data(aSAH)
myroc <- roc(aSAH$outcome, aSAH$s100b)
auc(myroc)
roc(outcome ~ s100b, aSAH, smooth=TRUE) 
plot(myroc, auc.polygon=TRUE)
plot(roc(outcome ~ s100b, aSAH, smooth=TRUE), auc.polygon=TRUE)

plot(myroc)

##ROCR
labels <- ifelse(aSAH$outcome == "Poor", 1, 0)
pred <- prediction(aSAH$s100b, labels)
perf <- performance(pred, "tpr", "fpr")
plot(perf)
abline(a=0, b= 1)

#They both have diagonal lines

####
data(ROCR.simple)

pred <- prediction(ROCR.simple$predictions,ROCR.simple$labels)
sn = slotNames(pred)
sapply(sn, function(x) length(slot(pred, x)))

sapply(sn, function(x) class(slot(pred, x)))
roc.perf = performance(pred, measure = "tpr", x.measure = "fpr")
par(pty='s') #makes the plot square
plot(roc.perf,main="ROC curve (ROCR)")
abline(a=0, b= 1)

roc.proc <- roc(ROCR.simple$labels, ROCR.simple$predictions)
plot(roc.proc, main="ROC curve (pROC)",legacy.axes=TRUE) #legacy axes makes it 1-spec instead of spec

#They're both step functions
N<-101
Np<-51
Nn<-50
K<-5
Nk<-20
Nnk<-10
Npk<-10
TPk<-10 
TNk<-10
SpecialFoldTP<-9
SpecialFoldTN<-9
TP<-(K-1)*TPk+SpecialFoldTP
TN<-(K-1)*TNk+SpecialFoldTN

MortenP<-(1-1/Np)*(1-1/(Nn))
MortenA<-(1/K)*((K-1)+(1-(1/(Np/K)))*(1-(1/(Nn/K))))
MelanieP<-(1/2)*((TN/Nn)+(TP/Np))
MelanieA<- (K-1)*(Nk/N)*(1/2)*(TNk/Nnk+TPk/Npk)+(Nk/N)*(1/2)*(SpecialFoldTP/Npk+SpecialFoldTN/Nnk)

(1/2)*((TN/Nn)+(49/51))
(K-1)*(Nk/N)*(1/2)*(TNk/Nnk+TPk/Npk)+(Nk/N)*(1/2)*(SpecialFoldTP/11+SpecialFoldTN/Nnk)

MortenP
MortenA
MelanieP
MelanieA


## Play data set
dl<-copy(kloop_BIGBIG) #[ K == 10 ] & Stat=="Averaged" Metric == "Precision" & type == "stratified" &
dl<-dl[!startsWith(dataset, "clean")]
dl<-dl[Metric != "AUC ROCR"]
dl[, prev :=
     fcase(
       startsWith(dataset, "sim_low"), "low",
       startsWith(dataset, "sim_mid"), "mid",
       startsWith(dataset, "sim_bal"), "bal",
       default = NA_character_)]
dl[, prev_num := factor(
  fcase(
    prev == "low",  0.01,
    prev == "mid",  0.10,
    prev == "bal",  0.50,
    default = NA_real_
  ),
  levels = c(0.01, 0.10, 0.50),
  labels = c("0.01", "0.1", "0.5")
)]
dl[,abs_diff := abs(diff_est)]
dl[, size_label := {
  # sort unique N within this prev level
  Ns <- sort(unique(N))
  # assign labels in order
  factor(N,
         levels = Ns,
         labels = size_labels[seq_along(Ns)])
}, by = prev]
dl[mc_dt, on = .(Metric, prev, size_label), MC := i.MC]
dl[, est.err := Value - MC]
#dl[,absest.err := abs(Value - MC)]

w <- setdiff(names(dl), c("Stat","Value","Var","Std", "est.err"))
dl <- dcast(
  dl,
  formula = as.formula(paste(paste(w, collapse = " + "), "~ Stat")),
  value.var = c("Value","Var","Std", "est.err") #,"absest.err"
)

w <- c("Value_Averaged", "Var_Averaged", "Std_Averaged", "Value_Pooled", "Var_Pooled", "Std_Pooled", "diff_est", "diff_var", "abs_diff", "est.err_Averaged","est.err_Pooled") #,"absest.err.1_Averaged","absest.err.1_Pooled"
dl <- dl[, {
  out <- list(n_seeds = .N)
  for (col in w) {
    x <- get(col)
    out[[paste0(col, "_mean")]]   <- mean(x, na.rm = TRUE)
    #    out[[paste0(col, "_sd")]]     <- sd(x, na.rm = TRUE)
    out[[paste0(col, "_se")]]     <- sd(x, na.rm = TRUE) / sqrt(sum(!is.na(x)))
    #    out[[paste0(col, "_median")]] <- median(x, na.rm = TRUE)
    #    out[[paste0(col, "_q025")]]   <- as.numeric(quantile(x, 0.025, na.rm = TRUE, names = FALSE))
    #    out[[paste0(col, "_q975")]]   <- as.numeric(quantile(x, 0.975, na.rm = TRUE, names = FALSE))
  }
  out
}, by =  c("Metric", "type", "N", "prev", "prev_num", "size_label", "dataset", "MC")]
dl[,MSE_average:= Var_Averaged_mean+est.err_Averaged_mean^2]
dl[,MSE_pool := Var_Pooled_mean+est.err_Pooled_mean^2]
#Summary tables
w<-dl[,c("Metric", "type", "N", "prev", "prev_num", "size_label","Var_Averaged_mean", "Var_Pooled_mean", "Value_Averaged_se", "Value_Pooled_se","est.err_Averaged_mean", "est.err_Pooled_mean","MSE_average","MSE_pool","absest.err.1_Averaged_mean","absest.err.1_Pooled_mean")]
w<-w[, ':='
     (
       Var_Averaged_mean=round(Var_Averaged_mean,6), 
       Var_Pooled_mean=round(Var_Pooled_mean,6) ,
       est.err_Averaged_mean = round(est.err_Averaged_mean,6),
       est.err_Pooled_mean = round(est.err_Pooled_mean,6),
       MSE_average = round(MSE_average,6),
       MSE_pool = round(MSE_pool,6)
     )]
w<-w[,diff_bias := est.err_Averaged_mean - est.err_Pooled_mean]
w<-w[,absdiff_bias := absest.err.1._Averaged_mean - absest.err.1._Pooled_mean]
w<-w[,diff_mse := MSE_average - MSE_pool]
w<-w[, ':='
     (
       Value_Averaged_mean=round(Value_Averaged_mean,6), 
       Value_Pooled_mean=round(Value_Pooled_mean,6) ,
       diff_est_mean=round(diff_est_mean,6),
       diff_est_q025=round(diff_est_q025,6),
       diff_est_q975=round(diff_est_q975,6),
       est.err_Averaged_mean = round(est.err_Averaged_mean,6),
       est.err_Pooled_mean = round(est.err_Pooled_mean,6)
     )]
w<-dl[,c("Metric", "type", "N", "prev", "prev_num", "size_label","Value_Averaged_mean","Value_Pooled_mean","diff_est_mean","abs_diff_mean", "Value_Averaged_se", "Value_Pooled_se")]
w<-w[,':='( 
  Value_Averaged_mean = round(Value_Averaged_mean,6),
  Value_Pooled_mean = round(Value_Pooled_mean,6),
  diff_est_mean = round(diff_est_mean,6),
  abs_diff_mean = round(abs_diff_mean,6),
  Value_Averaged_se=round(Value_Averaged_se,6), 
  Value_Pooled_se=round(Value_Pooled_se,6)
)]

dl <- melt(
  dl,
  measure.vars  = c("Value_Averaged_se", "Value_Pooled_se"),
  variable.name = "se_kind",
  value.name    = "se"
)

dl[, se_kind := fifelse(se_kind == "Value_Averaged_se", "Average", "Pool")]
dl<-dl[,c("Metric", "type", "N", "prev", "prev_num", "size_label", "se_kind", "se")]
#MSE
dl <- melt(
  dl,
  measure.vars  = c("MSE_average", "MSE_pool"),
  variable.name = "MSE_kind",
  value.name    = "MSE"
)

dl[, MSE_kind := fifelse(MSE_kind == "MSE_average", "Average", "Pool")]
dl[, size_label := factor(size_label, levels = c("Smallest", "Small"  ,  "Medium"  , "Large"  ,  "Largest"))]
dl[, col_id := interaction(MSE_kind, prev_num, drop = TRUE)]
w<- c("Average.0.01"=paste0("Average, ", "\u03C0","=0.01"),"Pool.0.01"=paste0("Pool, ", "\u03C0","=0.01"), "Average.0.1"=paste0("Average, ", "\u03C0","=0.1"), "Pool.0.1"=paste0("Pool, ", "\u03C0","=0.1") ,"Average.0.5"=paste0("Average, ", "\u03C0","=0.5") ,"Pool.0.5"=paste0("Pool, ", "\u03C0","=0.5"))
