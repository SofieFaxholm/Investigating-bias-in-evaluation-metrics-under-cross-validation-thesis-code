library(data.table)
library(ggplot2)
library(splitTools)
library(pROC)
library(ROCR)
######################################################
######################################################
######################################################
##Data 
######################################################
######################################################
######################################################
wdbc <- read.csv("~/Documents/Matematik/Speciale/SpecialeR/breast+cancer+wisconsin+diagnostic/wdbc.data", header=FALSE)
setDT(wdbc)
wdbc<-wdbc[,!c("V1")] 
wdbc$V2<-as.factor(wdbc$V2)
levels(wdbc$V2)[1]
labs <- levels(wdbc$V2)
wdbcmodel <- glm(V2~. , data=wdbc,family = binomial(link = "logit"))
summary(wdbcmodel)

ggplot(wdbc, aes_string(x = wdbc$V2)) +
  geom_bar(fill = c("B" = "lightgreen", "M" = "lightskyblue"), color = "black") +
  scale_x_discrete(labels = c("Negative", "Positive")) +
  labs(x = "Diagnosis", y = "Count") +
  theme_bw()


library(ggcorrplot)
corr <- cor(wdbc, use = "pairwise.complete.obs")
ggcorrplot(
  corr,
  type = "lower",       
  lab = FALSE,             
  outline.col = "white",   
  show.diag = TRUE,
  colors = c("darksalmon", "white", "navyblue")
) +
  ggplot2::theme_bw() +
  ggplot2::theme(
    axis.text.x  = ggplot2::element_blank(),  
    axis.ticks.x = ggplot2::element_blank(),
    axis.title.x = ggplot2::element_blank(),
    axis.title.y = ggplot2::element_blank() 
  )


######################################################
######################################################
######################################################
### Simulated data sets
######################################################
######################################################
######################################################
sim_fun <- function(n, target_prev) {
  dt_sim_reduced <- test
  W_cov          <- setdiff(names(dt_sim_reduced), "V2")
  if (is.null(target_prev)) {

    lev <- levels(test$V2)
    target_prev <- mean(test$V2 == lev[2])
  }
  Ms <- round(n * target_prev) # otherwise the "exact" prevalence is not guaranteed
  Bs <- n - Ms
  total  <- Bs + Ms
  n_rows <- nrow(dt_sim_reduced)
  lev    <- levels(dt_sim_reduced$V2)   # e.g. c("B","M")
  
  out_list <- vector("list", total)
  y_vec    <- character(total)
  
  i       <- 1L
  count_B <- 0L
  count_M <- 0L
  
  while (count_B < Bs || count_M < Ms) {
    idx   <- sample.int(n_rows, 1L)
    W_row <- dt_sim_reduced[idx, ..W_cov]
    W_row[, "(Intercept)" := 1]
    
    p <- predict(testmodel, newdata = W_row, type = "response")
    y <- rbinom(1L, 1L, p)
    
    if (y == 1L && count_M < Ms) {
      cls     <- lev[2]   # e.g. "M"
      count_M <- count_M + 1L
    } else if (y == 0L && count_B < Bs) {
      cls     <- lev[1]   # e.g. "B"
      count_B <- count_B + 1L
    } else {
      next
    }
    
    W_row[, "(Intercept)" := NULL]
    out_list[[i]] <- W_row
    y_vec[i]      <- cls
    i             <- i + 1L
  }
  
  sim_data <- rbindlist(out_list)
  sim_data[, V2 := factor(y_vec, levels = lev)]
  
  sim_data[]
}

size_labels <- c("Smallest", "Small", "Medium", "Large", "Largest") 
lowpi_Ns <- c(10007,19997,40013, 79997, 160033)
midpi_Ns <- c(997,2003,4003,7993, 16007)
balpi_Ns <- c(197,397,797,1607, 3203)

#simulates the different data sets (individually)
for (i in seq_along(lowpi_Ns)){
  N<-lowpi_Ns[i]
  assign(
    paste0("sim_low", i),
    sim_fun(n = N, target_prev = 0.01)
  )
}
for (i in seq_along(midpi_Ns)){
  N<-midpi_Ns[i]
  assign(
    paste0("sim_mid", i),
    sim_fun(n = N, target_prev = 0.1)
  )
}
for (i in seq_along(balpi_Ns)){
  N<-balpi_Ns[i]
  assign(
    paste0("sim_bal", i),
    sim_fun(n = N, target_prev = 0.5)
  )
}

cleanlowpiN <-c(10000,20000,40000,80000, 160000)
cleanmidpiN <-c(1000,2000,4000,8000,16000)
cleanbalpiN <-c(200,400,800,1600,3200)
for (i in seq_along(cleanlowpiN)){
  N<-cleanlowpiN[i]
  assign(
    paste0("clean_sim_low", i),
    sim_fun(n = N, target_prev = 0.01)
  )
}
for (i in seq_along(cleanmidpiN)){
  N<-cleanmidpiN[i]
  assign(
    paste0("clean_sim_mid", i),
    sim_fun(n = N, target_prev = 0.1)
  )
}
for (i in seq_along(cleanbalpiN)){
  N<-cleanbalpiN[i]
  assign(
    paste0("clean_sim_bal", i),
    sim_fun(n = N, target_prev = 0.5)
  )
}
######################################################
######################################################
######################################################
## Monte Carlo data sets and estimates 
######################################################
######################################################
######################################################

#MCsim_low <- sim_fun(n = 1e6, target_prev = 0.01)
#MCsim_mid <- sim_fun(n = 1e6, target_prev = 0.1)
#MCsim_bal <- sim_fun(n = 1e6, target_prev = 0.5)

# This way avoids having three 1e6 sized data sets (didn't have enough memory)
mc <- function(fit, n, target_prev, test, testmodel, threshold, labs, chunk=50000) {
  dt <- as.data.table(test)
  W  <- setdiff(names(dt), "V2")
  lev <- levels(dt$V2)
  
  Ms <- round(n * target_prev); Bs <- n - Ms
  pred <- numeric(n); y1 <- logical(n)
  m <- 0L; cM <- 0L; cB <- 0L; N0 <- nrow(dt)
  
  while (m < n) {
    idx <- sample.int(N0, chunk, replace=TRUE)
    X <- dt[idx, ..W]
    X2 <- copy(X); X2[, "(Intercept)" := 1]
    yprop <- rbinom(chunk, 1L, as.numeric(predict(testmodel, X2, type="response")))
    
    keep <- logical(chunk)
    for (j in 1:chunk) {
      if (yprop[j]==1L) { if (cM < Ms) { keep[j]=TRUE; cM=cM+1L } }
      else             { if (cB < Bs) { keep[j]=TRUE; cB=cB+1L } }
      if (cM>=Ms && cB>=Bs) break
    }
    if (!any(keep)) next
    
    Xk <- X[keep]
    yk <- yprop[keep]
    pk <- as.numeric(predict(fit, Xk, type="response"))
    
    take <- min(length(pk), n-m)
    pred[(m+1):(m+take)] <- pk[1:take]
    y1[(m+1):(m+take)]   <- (yk[1:take]==1L)
    m <- m + take
  }
  
  cls <- factor(ifelse(pred > threshold, labs[2], labs[1]), levels=labs)
  tru <- factor(ifelse(y1, labs[2], labs[1]), levels=labs)
  
  TP <- sum(cls==labs[2] & tru==labs[2])
  FP <- sum(cls==labs[2] & tru==labs[1])
  FN <- sum(cls==labs[1] & tru==labs[2])
  TN <- sum(cls==labs[1] & tru==labs[1])
  
  data.table(
    Metric = c("err","spec","prec","auc"),
    MC = c(
      (FP+FN)/(TP+TN+FP+FN),
      TN/(TN+FP),
      TP/(TP+FP),
      as.numeric(pROC::auc(tru, pred, levels=labs))
    )
  )
}
#this estimates the MC estimates from the 15 data sets 
sims_dtLIST <- list(
  low = list(sim_low1, sim_low2, sim_low3, sim_low4, sim_low5),
  mid = list(sim_mid1, sim_mid2, sim_mid3, sim_mid4, sim_mid5),
  bal = list(sim_bal1, sim_bal2, sim_bal3, sim_bal4, sim_bal5)
)
mc_15 <- function(sims_dtLIST, test, testmodel, threshold, labs, n_mc=1e6, chunk=50000) {
  out <- list(); k <- 1L
  
  for (prev in names(sims_dtLIST)) {
    target_prev <- if (prev=="low") 0.01 else if (prev=="mid") 0.10 else 0.50
    dts <- sims_dtLIST[[prev]]  # list of 5 datasets
    
    for (i in 1:5) {
      dt  <- as.data.table(dts[[i]])
      fit <- glm(V2 ~ ., data = dt, family = binomial())
      
      tmp <- mc(fit, n_mc, target_prev, test, testmodel, threshold, labs, chunk)
      tmp[, `:=`(prev = prev, size_label = size_labels[i])]
      
      out[[k]] <- tmp
      k <- k + 1L
    }
  }
  
  rbindlist(out)
}


mc_dt <- mc_15(sims_dtLIST, test, testmodel, threshold, labs, n_mc=1e6, chunk=50000)

w<- c(
  err  = "Misclassification error",
  spec = "Specificity",
  prec = "Precision",
  auc  = "AUC pROC"
)
mc_dt<-mc_dt[, Metric := w[Metric]]
