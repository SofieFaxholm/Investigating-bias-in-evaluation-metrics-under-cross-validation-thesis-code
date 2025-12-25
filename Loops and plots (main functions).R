
seeds <- sample.int(.Machine$integer.max, 200)
K <- 10
Ks <- c(5,10)
threshold <- 0.5

res <- c(
  paste0("clean_sim_bal", 1:5), 
  paste0("clean_sim_mid", 1:5),
  paste0("clean_sim_low", 1:5),
  paste0("sim_bal", 1:5),
  paste0("sim_mid", 1:5),
  paste0("sim_low", 1:5)
)

## MAIN FUNCTION - does cv, fits logistic, evaluates etc
cv_types_metric_comp <- function(data, K, seed, type = c("stratified","basic"), threshold) {
  type <- match.arg(type)
  
  labs <- levels(data$V2)            
  n    <- nrow(data)
  
  set.seed(seed)
  rFolds <- create_folds(data$V2, k = K, type = type, invert = TRUE, seed = seed)
  
  averageMetric <- data.table(
    a.err = rep(NA_real_, K),
    a.spec = NA_real_,
    a.prec = NA_real_,
    a.auc_pROC  = NA_real_,
    a.auc_ROCR = NA_real_,
    weight = NA_real_,
    varerr = NA_real_,
    varspec =NA_real_,
    varprec = NA_real_,
    varauc_pROC = NA_real_
  ) 
  
  
  poolMetric <- data.table(
    prob = rep(NA_real_, n),
    pred = factor(rep(NA_character_, n), levels = labs),
    fold = rep(NA_integer_, n)
  )
  
  for (i in seq_along(rFolds)) {
    te <- rFolds[[i]]
    tr <- setdiff(seq_len(n), te)
    
    testset  <- data[te, ]
    trainset <- data[tr, ]
    fold <- i
    
    foldfit     <- glm(V2 ~ ., data = trainset, family = binomial(link = "logit"))
    predictions <- predict(foldfit, newdata = testset, type = "response")
    
    classifier <- factor(ifelse(predictions > threshold, labs[2], labs[1]), levels = labs)
    
    true <- testset$V2
    TP <- sum(classifier == labs[2] & true == labs[2])
    FP <- sum(classifier == labs[2] & true == labs[1])
    FN <- sum(classifier == labs[1] & true == labs[2])
    TN <- sum(classifier == labs[1] & true == labs[1])
    
    err  <- (FP + FN) / (TP + TN + FP + FN)
    spec <- TN / (TN + FP)
    prec <- TP / (TP + FP)
    auc_pROC <- as.numeric(pROC::auc(response = true, predictor = predictions, levels = labs))
    pred_rocr <- ROCR::prediction(predictions, true, label.ordering = labs)
    #   auc_ROCR  <- as.numeric(ROCR::performance(pred_rocr, "auc")@y.values[[1]]) 
    auc_ROCR <- NA_real_ #I removed all ROCRs for faster running times 
    
    var.err <- (err*(1-err)*length(te))/(length(te)-1)
    var.spec <- spec * (1 - spec) * (TN + FP)/ (TN + FP -1)
    var.prec <- prec * (1 - prec)*(TP + FP) / (TP + FP-1)
    var.auc_pROC <-pROC::var(roc(response = true, predictor = predictions, levels = labs))
    #    var.auc_ROCR <- var(ROCR::performance(pred_rocr, "auc")@y.values[[1]]) #ROCR package doesn't have a variance function
    
    averageMetric[i, `:=`(
      a.err  = err,
      a.spec = spec,
      a.prec = prec,
      a.auc_pROC = auc_pROC,
      a.auc_ROCR = auc_ROCR,
      weight = length(te),
      varerr = var.err,
      varspec = var.spec,
      varprec = var.prec,
      varauc_pROC = var.auc_pROC
      #      varauc_ROCR = var.auc_ROCR
    )]
    
    poolMetric[te, `:=`(prob = predictions, pred = classifier, fold = fold)]
  }
  
  averageERR  <- weighted.mean(averageMetric$a.err,  averageMetric$weight)
  averageSPEC <- weighted.mean(averageMetric$a.spec, averageMetric$weight)
  averagePREC <- weighted.mean(averageMetric$a.prec, averageMetric$weight)
  averageAUC_pROC  <- weighted.mean(averageMetric$a.auc_pROC,  averageMetric$weight)
  #  averageAUC_ROCR  <- weighted.mean(averageMetric$a.auc_ROCR,  averageMetric$weight)
  averageAUC_ROCR <- NA_real_ 
  
  poolMetric[, true := data$V2]
  pTP <- sum(poolMetric$pred == labs[2] & poolMetric$true == labs[2])
  pFP <- sum(poolMetric$pred == labs[2] & poolMetric$true == labs[1])
  pFN <- sum(poolMetric$pred == labs[1] & poolMetric$true == labs[2])
  pTN <- sum(poolMetric$pred == labs[1] & poolMetric$true == labs[1])
  
  poolERR  <- (pFP + pFN) / (pTP + pTN + pFP + pFN)
  poolSPEC <- pTN / (pTN + pFP)
  poolPREC <- pTP / (pTP + pFP)
  poolAUC_pROC  <- as.numeric(pROC::auc(response = poolMetric$true, predictor = poolMetric$prob, levels = labs))
  pred_rocr_pool <- ROCR::prediction(poolMetric$prob, poolMetric$true, label.ordering = labs)
  #  poolAUC_ROCR  <- as.numeric(ROCR::performance(pred_rocr_pool, "auc")@y.values[[1]])
  poolAUC_ROCR <- NA_real_
  
  w <- averageMetric$weight/sum(averageMetric$weight)
  var.a.err <- round(sum(averageMetric$varerr*(w^2)),6)
  var.a.spec <- round(sum(averageMetric$varspec*(w^2)),6)
  var.a.prec <- round(sum(averageMetric$varprec*(w^2)),6)
  var.a.auc_pROC <- sum(averageMetric$varauc_pROC*(w^2))
  #  var.a.auc_ROCR <- sum(averageMetric$varauc_ROCR*(w^2))
  
  var.p.err <- round(poolERR * (1 - poolERR) *(pTP + pTN + pFP + pFN)/ (pTP + pTN + pFP + pFN - 1),6)
  var.p.spec <- round(poolSPEC * (1 - poolSPEC)*(pTN + pFP) / (pTN + pFP - 1),6)
  var.p.prec <- round(poolPREC * (1 - poolPREC)*(pTP + pFP) / (pTP + pFP - 1),6)
  var.p.auc_pROC <-pROC::var(roc(response = poolMetric$true,
                                 predictor = poolMetric$prob,
                                 levels   = labs))
  #  var.p.auc_ROCR <- var(ROCR::performance(pred_rocr_pool, "auc")@y.values[[1]]) 
  
  std.a.err <- round(sqrt(var.a.err),5)
  std.a.spec <- round(sqrt(var.a.spec),5)
  std.a.prec <- round(sqrt(var.a.prec),5)
  std.a.auc_pROC <- round(sqrt(var.a.auc_pROC),5)
  #  std.a.auc_ROCR <- round(sqrt(var.a.auc_ROCR),5)
  
  std.p.err <- round(sqrt(var.p.err),5)
  std.p.spec <- round(sqrt(var.p.spec),5)
  std.p.prec <- round(sqrt(var.p.prec),5)
  std.p.auc_pROC <-round(sqrt(var.p.auc_pROC),5)
  #  std.p.auc_ROCR <- round(sqrt(var.p.auc_ROCR),5)
  
  compare <- data.table(
    Metric   = factor(c("Misclassification error","Specificity","Precision","AUC pROC", "AUC ROCR"),
                      levels = c("Misclassification error","Specificity","Precision","AUC pROC", "AUC ROCR" )),
    Averaged = c(averageERR, averageSPEC, averagePREC, averageAUC_pROC, averageAUC_ROCR),
    Pooled   = c(poolERR,    poolSPEC,    poolPREC,    poolAUC_pROC, poolAUC_ROCR),
    Averaged_var = c(var.a.err, var.a.spec, var.a.prec, var.a.auc_pROC, NA_real_),
    Pooled_var = c(var.p.err, var.p.spec, var.p.prec, var.p.auc_pROC, NA_real_),
    Averaged_std = c(std.a.err, std.a.spec, std.a.prec, std.a.auc_pROC, NA_real_),
    Pooled_std = c(std.p.err, std.p.spec, std.p.prec, std.p.auc_pROC, NA_real_) 
  )
  
  list(
    by_fold = averageMetric[],
    pooled  = poolMetric[, .(true, prob, pred)],
    summary = compare[]
  )
}

######################################################
######################################################
######################################################
## Loops over seeds
loop_seeds <- function(data, seeds, K, type = c("stratified","basic"), threshold) {
  type <- match.arg(type)
  rbindlist(lapply(seeds, function(s) {
    res <- cv_types_metric_comp(data, K = K, seed = s, type = type, threshold = threshold)
    dt  <- as.data.table(res$summary)                        
    val_long <- melt(dt, id.vars = "Metric", measure.vars = c("Averaged", "Pooled"), variable.name = "Stat", value.name = "Value")
    var_long <- melt(dt, id.vars = "Metric", measure.vars = c("Averaged_var", "Pooled_var"), variable.name = "Stat_var", value.name = "Var")
    std_long <- melt(dt, id.vars = "Metric", measure.vars = c("Averaged_std", "Pooled_std"), variable.name = "Stat_std", value.name = "Std")
    out <- cbind(val_long, Var = var_long$Var, Std = std_long$Std)
    out[ , `:=`(seed = s, K = K, type = type)]
  }), use.names = TRUE)
}

seedloop_basic <- loop_seeds(test, seeds, K,"basic", threshold)
seedloop_strat <- loop_seeds(test, seeds, K,"stratified", threshold)
seedloop_combined <- rbindlist(list(seedloop_basic, seedloop_strat), use.names = TRUE)

seedloop_combined[, Value := round(Value, 6)]
seedloop_combined[, type := factor(type, levels = c("basic","stratified"))]
seedloop_combined[, seed := factor(seed)]
seedloop_combined[, K := factor(K)]
seedloop_combined[, diff := {
  a <- Value[Stat == "Averaged"][1]
  p <- Value[Stat == "Pooled"][1]
  ifelse(Stat == "Averaged", a - p, p - a)
}, by = .(Metric, type, seed, K)]

seedloop_vardt <- dcast(
  seedloop_combined,
  Metric + type + K + seed ~ Stat,
  value.var = c("Value", "Var")  # creates Value_Averaged, Value_Pooled, Var_Averaged, Var_Pooled
)
seedloop_vardt[, `:=`(
  diff_mean     = Value_Averaged - Value_Pooled, 
  diff_var_est  = Var_Averaged   - Var_Pooled   
)]
seedloop_vardt <- seedloop_vardt[Metric != "AUC ROCR"]
vals_long <- seedloop_vardt[, {
  m  <- mean(diff_var_est)
  s  <- sd(diff_var_est)
  n  <- .N
  .(
    mean_diff_var   = m,
    median_diff_var = median(diff_var_est),
    sd_diff_var     = s,
    se_diff_var     = s / sqrt(n),
    # relative difference: mean diff relative to average of the two methods' mean variances
    rel_mean_diff_var = m / mean(c(mean(Var_Averaged), mean(Var_Pooled))),
    # Wilcoxon signed-rank test for median difference = 0
    p_signed_rank = if (all(diff_var_est == 0)) NA_real_ else
      wilcox.test(diff_var_est, mu = 0, alternative = "two.sided", exact = FALSE)$p.value
  )
}, by = .(Metric, type, K)]


### Loops over seeds, just like above, but also spits out Scompare
# Note name is to limit number of things in environment
loop_seeds_wCompare <- function(data, seeds, K, threshold) {
  seedloop_basic <- loop_seeds(data, seeds, K, "basic",      threshold)
  seedloop_strat <- loop_seeds(data, seeds, K, "stratified", threshold)
  seedloop_combined <- rbindlist(list(seedloop_basic, seedloop_strat), use.names = TRUE)
  
  seedloop_combined[, Value := round(Value, 6)]
  seedloop_combined[, type := factor(type, levels = c("basic","stratified"))]
  seedloop_combined[, seed := factor(seed)]
  seedloop_combined[, K    := factor(K)]
  seedloop_combined[, diff := {
    a <- Value[Stat == "Averaged"][1]
    p <- Value[Stat == "Pooled"][1]
    ifelse(Stat == "Averaged", a - p, p - a)
  }, by = .(Metric, type, seed, K)]
  
  seedloop_vardt <- data.table::dcast(
    seedloop_combined,
    Metric + type + K + seed ~ Stat,
    value.var = c("Value", "Var")  
  )
  
  seedloop_vardt[, `:=`(
    diff_mean    = Value_Averaged - Value_Pooled,
    diff_var_est = Var_Averaged   - Var_Pooled
  )]
  
  seedloop_summary <- seedloop_vardt[
    , .(
      Value_Averaged = mean(Value_Averaged),
      Value_Pooled   = mean(Value_Pooled),
      Var_Averaged   = mean(Var_Averaged),
      Var_Pooled     = mean(Var_Pooled),
      diff_mean      = mean(diff_mean),
      diff_var_est   = mean(diff_var_est)
    ),
    by = .(Metric, type, K)
  ]
  
  seedloop_summary[, N := nrow(data)]
  
  list(
    seedloop_combined = seedloop_combined,
    Scompare    = seedloop_summary
  )
}


## Histogram of raw estimates
plot_hist_type_metric <- function(dt, metric_name,
                                  cv_type = c("basic","stratified"),
                                  bins = NULL, binwidth = NULL,
                                  use_density = FALSE, metric_v, prev_v, sizelabel) {
  cv_type <- match.arg(cv_type)
  
  # filter rows for this metric & type; keep both Averaged and Pooled
  d <- dt[Metric == metric_name & type == cv_type & Stat %in% c("Averaged","Pooled")]
  d[, Stat := factor(Stat, levels = c("Averaged","Pooled"))]
  
  vline_x <- mc_dt[Metric == metric_v & prev == prev_v & size_label == sizelabel, MC][1]
  
  # default binwidth if neither bins nor binwidth supplied
  if (is.null(bins) && is.null(binwidth)) {
    bw <- 2 * IQR(d$Value) / (length(d$Value)^(1/3))
    if (!is.finite(bw) || bw <= 0) {
      rng <- diff(range(d$Value))
      bw <- ifelse(is.finite(rng) && rng > 0, rng/30, 0.1)
    }
    binwidth <- bw
  }
  
  title_type <- if (cv_type == "basic") "Random" else "Stratified"
  title_metric_name <- if(metric_name == "AUC pROC") "AUC" else metric_name
  
  p <- ggplot(d, aes(x = Value, fill = Stat))
  
  if (use_density) {
    p <- p + {
      if (is.null(bins))
        geom_histogram(aes(y = after_stat(density)),
                       position = "identity", alpha = 0.45,
                       binwidth = binwidth, color = NA)
      else
        geom_histogram(aes(y = after_stat(density)),
                       position = "identity", alpha = 0.45,
                       bins = bins, color = NA)
    } +
      labs(y = "Density")
  } else {
    p <- p + {
      if (is.null(bins))
        geom_histogram(position = "identity", alpha = 0.45,
                       binwidth = binwidth, color = NA)
      else
        geom_histogram(position = "identity", alpha = 0.45,
                       bins = bins, color = NA)
    } +
      labs(y = "Count")
  }
  p +
    geom_vline(xintercept = vline_x, linetype = "dashed", color = "skyblue4") +
    scale_fill_manual(values = c(Averaged = "seagreen", Pooled = "orange"), name = NULL) +
    labs(title = paste0(title_metric_name, " — ", title_type), x = "Estimate") +
    theme_bw(base_size = 13) +
    theme(plot.title = element_text(hjust = 0.5))
}
plot_hist_type_metric(seedloop_combined, "Misclassification error", cv_type = "basic", bins = 50, use_density = TRUE)

## Histogram of differences
# Note in data set "wdbc" and "sim" column is "diff" in simulated sets (clean..,sim_low etc) it's "diff_est"
plot_hist_type_metric_diff <- function(dt, metric_name,
                                       cv_type = c("basic","stratified"),
                                       bins = NULL, binwidth = NULL,
                                       use_density = FALSE) {
  cv_type <- match.arg(cv_type)
  
  # keep only Averaged rows so diff = Averaged - Pooled
  d <- dt[Metric == metric_name & type == cv_type & Stat == "Averaged"]
  
  # default binwidth if neither bins nor binwidth supplied
  xvar <- d$diff
  if (is.null(bins) && is.null(binwidth)) {
    bw <- 2 * IQR(xvar, na.rm = TRUE) / (length(xvar)^(1/3))
    if (!is.finite(bw) || bw <= 0) {
      rng <- diff(range(xvar, na.rm = TRUE))
      bw <- ifelse(is.finite(rng) && rng > 0, rng/30, 0.1)
    }
    binwidth <- bw
  }
  
  title_type <- if (cv_type == "basic") "Random" else "Stratified"
  title_metric_name <- if (metric_name == "AUC pROC") "AUC" else metric_name
  
  p <- ggplot(d, aes(x = diff))
  
  if (use_density) {
    if (is.null(bins)) {
      p <- p + geom_histogram(
        aes(y = after_stat(density)),
        position = "identity",
        alpha = 0.7,
        binwidth = binwidth,
        color = NA,
        fill = "#4e79a7"
      )
    } else {
      p <- p + geom_histogram(
        aes(y = after_stat(density)),
        position = "identity",
        alpha = 0.7,
        bins = bins,
        color = NA,
        fill = "#4e79a7"
      )
    }
    p <- p + labs(y = "Density")
  } else {
    if (is.null(bins)) {
      p <- p + geom_histogram(
        position = "identity",
        alpha = 0.7,
        binwidth = binwidth,
        color = NA,
        fill = "#4e79a7"
      )
    } else {
      p <- p + geom_histogram(
        position = "identity",
        alpha = 0.7,
        bins = bins,
        color = NA,
        fill = "#4e79a7"
      )
    }
    p <- p + labs(y = "Count")
  }
  
  p +
    geom_vline(xintercept = 0, linetype = "dashed") +
    labs(
      title = paste0(title_metric_name, " — ", title_type),
      x = "Averaged - Pooled"
    ) +
    theme_bw(base_size = 13) +
    theme(plot.title = element_text(hjust = 0.5))
}
plot_hist_type_metric_diff(seedloop_combined, "Misclassification error", cv_type = "basic", use_density = TRUE)

######################################################
######################################################
######################################################
#### Loop over K and seeds
loop_Ks <- function(data, Ks, seeds, type = c("stratified","basic"), threshold) {
  type <- match.arg(type)
  rbindlist(lapply(Ks, function(K) {
    loop_seeds(data = data, seeds = seeds, K = K, type = type, threshold = threshold)
  }), use.names = TRUE)
}

kloop_basic <- loop_Ks(sim, Ks, seeds,"basic", threshold)
kloop_strat <- loop_Ks(sim, Ks, seeds, "stratified", threshold)
kloop_combined <- rbindlist(list(kloop_basic, kloop_strat), use.names = TRUE)

kloop_combined[, Value := round(Value, 6)]
kloop_combined[, type := factor(type, levels = c("basic","stratified"))]
kloop_combined[, seed := factor(seed)]
kloop_combined[, K := factor(K)]
kloop_combined[, diff := {
  a <- Value[Stat == "Averaged"][1]
  p <- Value[Stat == "Pooled"][1]
  ifelse(Stat == "Averaged", a - p, p - a)
}, by = .(Metric, type, seed, K)]
#std errors (recall that errors is the difference between seeds)
kloop_combinedSUM <- kloop_combined[
  , .(
    n             = sum(Stat == "Averaged"),
    mean          = round(mean(diff[Stat == "Averaged"]), 5),
    se            = round(sd(diff[Stat == "Averaged"]) / sqrt(sum(Stat == "Averaged")), 5),
    mean_average  = round(mean(Value[Stat == "Averaged"]), 5),
    mean_pooled   = round(mean(Value[Stat == "Pooled"]), 5)
  ), by = .(Metric, K, type)]

View(kloop_combinedSUM[K==10])

######################################################
######################################################
######################################################
#### Loop all 30 data sets over K and seeds
process_K_dataset <- function(data, Ks, seeds, threshold, data_name = NULL) {
  if (is.null(data_name)) data_name <- deparse(substitute(data))
  N_val <- nrow(data)
  
  kloop_basic <- loop_Ks(data, Ks, seeds, "basic",      threshold)
  kloop_strat <- loop_Ks(data, Ks, seeds, "stratified", threshold)
  
  kloop_combined <- rbindlist(list(kloop_basic, kloop_strat), use.names = TRUE)

  kloop_combined[, Value := round(Value, 6)]
  kloop_combined[, type := factor(type, levels = c("basic","stratified"))]
  kloop_combined[, seed := factor(seed)]
  kloop_combined[, K    := factor(K)]
  kloop_combined[, diff_est := {
    a <- Value[Stat == "Averaged"][1]
    p <- Value[Stat == "Pooled"][1]
    a - p
  }, by = .(Metric, type, seed, K)]
  
  kloop_combined[, diff_var := {
    va <- Var[Stat == "Averaged"][1]
    vp <- Var[Stat == "Pooled"][1]
    va - vp
  }, by = .(Metric, type, seed, K)]
  

  kloop_combined[, `:=`(N = N_val, dataset = data_name)]
  kloop_combinedSUM <- kloop_combined[
    , .(
      N   = N_val,
      est = mean(Value),
      var = mean(Var),
      std = mean(Std)
    ),
    by = .(Metric, K, type, Stat)
  ]
  kloop_combinedSUM[
    , `:=`(
      mean_diff_est = est[Stat == "Averaged"] - est[Stat == "Pooled"],
      mean_diff_var = var[Stat == "Averaged"] - var[Stat == "Pooled"]
    ),
    by = .(Metric, K, type)
  ]
  
  kloop_combinedSUM[, dataset := data_name]
  
  list(
    combined     = kloop_combined,
    combined_sum = kloop_combinedSUM
  )
}
#runs all 30 sim sets (carefull before running +10h running time)
kloop_results <- setNames(
  lapply(res, function(nm) {
    process_K_dataset(
      data      = get(nm),      
      Ks        = K,
      seeds     = seeds,
      threshold = threshold,
      data_name = nm            
    )
  }),
  res
)

kloop_combLIST    <- lapply(kloop_results, `[[`, "combined")
kloop_combSumLIST <- lapply(kloop_results, `[[`, "combined_sum")
names(kloop_combLIST)    <- res
names(kloop_combSumLIST) <- res

#one big summary table across all 30
kloop_sumBIG <- rbindlist(kloop_combSumLIST, use.names = TRUE)
options(scipen = 999)
kloop_sumBIG[, mean_diff_est := round(mean_diff_est, 7)]
kloop_sumBIG[, mean_diff_var := round(mean_diff_var, 7)]
kloop_sumBIG <- kloop_sumBIG[     #this one for the summaries
  Metric == "Specificity" &
    K == 10 &
    type == "stratified" &
    Stat == "Averaged",
  .(dataset, est, var, mean_diff_est, mean_diff_var)
]

kloop_BIGBIG <- rbindlist(kloop_combLIST, use.names = TRUE)

#for each prevalence
dt_low <- rbindlist(list(kloop_combLIST$sim_low1,kloop_combLIST$sim_low2,kloop_combLIST$sim_low3,kloop_combLIST$sim_low4,kloop_combLIST$sim_low5), use.names = TRUE, fill = TRUE)
dt_mid <- rbindlist(list(kloop_combLIST$sim_mid1,kloop_combLIST$sim_mid2,kloop_combLIST$sim_mid3,kloop_combLIST$sim_mid4,kloop_combLIST$sim_mid5), use.names = TRUE, fill = TRUE)
dt_bal <- rbindlist(list(kloop_combLIST$sim_bal1,kloop_combLIST$sim_bal2,kloop_combLIST$sim_bal3,kloop_combLIST$sim_bal4,kloop_combLIST$sim_bal5), use.names = TRUE, fill = TRUE)

dt_cleanlow <- rbindlist(list(kloop_combLIST$clean_sim_low1,kloop_combLIST$clean_sim_low2,kloop_combLIST$clean_sim_low3,kloop_combLIST$clean_sim_low4,kloop_combLIST$clean_sim_low5), use.names = TRUE, fill = TRUE)
dt_cleanmid <- rbindlist(list(kloop_combLIST$clean_sim_mid1,kloop_combLIST$clean_sim_mid2,kloop_combLIST$clean_sim_mid3,kloop_combLIST$clean_sim_mid4,kloop_combLIST$clean_sim_mid5), use.names = TRUE, fill = TRUE)
dt_cleanbal <- rbindlist(list(kloop_combLIST$clean_sim_bal1,kloop_combLIST$clean_sim_bal2,kloop_combLIST$clean_sim_bal3,kloop_combLIST$clean_sim_bal4,kloop_combLIST$clean_sim_bal5), use.names = TRUE, fill = TRUE)


#### Plots

##Boxplot multiple Ks one data set
plot_box <- function(dt, metric_name) {
  d <- dt[Metric == metric_name & Stat == "Averaged" & !is.na(diff) & K != 20 & K != 3 & K!=5]
  d[, type_label := factor(ifelse(type == "basic", "Random", "Stratified"),
                           levels = c("Random","Stratified"))]
  d[, K := factor(K)]  
  
  Ks <- levels(d$K)
  x_levels <- as.vector(rbind(paste0(Ks, "\nRandom"),
                              paste0(Ks, "\nStratified")))
  d[, xlab := factor(paste0(K, "\n", type_label), levels = x_levels)]
  
  ggplot(d, aes(x = xlab, y = diff, fill = type_label, color = type_label)) +
    geom_hline(yintercept = 0, linetype = "dashed") +
    geom_boxplot(outlier.shape = NA, alpha = 0.25, width = 0.75, color = NA) +
    geom_jitter(width = 0, height = 0, size = 2, alpha = 0.7) +
    stat_summary(aes(color = type_label), fun = mean, geom = "crossbar",
                 position = position_dodge(width = 0.6)) + 
    # show only the K values under each pair
    scale_x_discrete(labels = rep(Ks, each = 2)) +
    scale_fill_manual(values = c(Random = "#d1495b", Stratified = "#2e86ab"), guide = "none") +
    scale_color_manual(values = c(Random = "#d1495b", Stratified = "#2e86ab"), name = "CV type") +
    labs(title = metric_name, x = NULL, y = "Difference") +
    theme_bw(base_size = 13) +
    theme(axis.text.x = element_text(lineheight = 0.9),
          panel.grid.major.x = element_blank(),
          plot.title = element_text(hjust = 0.5))
}

metrics <- unique(kloop_combined$Metric)
plots_box <- lapply(metrics, function(m) plot_box(seedloop_combined, m))
names(plots_box) <- metrics
print(plots_box[["Misclassification error"]])
print(plots_box[["Specificity"]])
print(plots_box[["Precision"]])
print(plots_box[["AUC pROC"]])

##Boxplots over N instead of K, for a single splitting type
plot_box_by_N <- function(dt, metric_name, K_val, cv_type = c("basic", "stratified")) {
  cv_type <- match.arg(cv_type)
  
  type_label <- if (cv_type == "basic") "Random" else "Stratified"
  type_color <- if (cv_type == "basic") "#d1495b" else "#2e86ab"
  
  setDT(dt)
  
  d <- dt[
    Metric == metric_name &
      Stat == "Averaged" &
      !is.na(diff_est) &
      K == K_val &
      type == cv_type
  ]
  
  d[, N := factor(N, levels = sort(unique(N)))]
  
  ggplot(d, aes(x = N, y = diff_est)) +
    geom_hline(yintercept = 0, linetype = "dashed") +
    geom_boxplot(alpha = 0.25,
                 width = 0.6,
                 fill = type_color,
                 colour = NA) +
    geom_jitter(width = 0, height = 0, size = 2, alpha = 0.7) +
    stat_summary(
      fun = mean,
      geom = "crossbar",
      width = 0.4,
      colour = type_color
    ) +
    #    coord_cartesian(ylim = c(-0.01, 0.01)) +   #set specific y-axis when they're annoying
    labs(
      title = paste0(metric_name, " - Balanced" ),
      x = "N",
      y = "Difference"
    ) +
    theme_bw(base_size = 13) +
    theme(
      panel.grid.major.x = element_blank(),
      plot.title = element_text(hjust = 0.5)
    )
}
#Random  metric_name,
plots_byN <- setNames(
  lapply(metrics, function(m)
    plot_box_by_N(dt_bal,
                  metric_name = m,
                  K_val = 10,
                  cv_type = "basic")),
  metrics
)

#print(plots_byN[["Misclassification error"]])
print(plots_byN[["Specificity"]])
print(plots_byN[["Precision"]])
print(plots_byN[["AUC pROC"]])
#Stratified
plots_byN <- setNames(
  lapply(metrics, function(m)
    plot_box_by_N(dt_bal,
                  metric_name = m,
                  K_val = 10,
                  cv_type = "stratified")),
  metrics
)
#print(plots_byN[["Misclassification error"]])
print(plots_byN[["Specificity"]])
print(plots_byN[["Precision"]])
print(plots_byN[["AUC pROC"]])


#### Wilcoxon
### Wilcoxon
wilc<-kloop_combined[Stat == "Averaged" & K==10]
wilcoxon <- wilc[
  , {
    w <- wilcox.test(diff, mu = 0, alternative = "two.sided", exact = FALSE)
    data.table(n = .N, median_diff = median(diff), p = round(w$p.value,5))
  },
  by = .(Metric, type, K)
][, p_adj := p.adjust(p, method = "BH")]  # multiple-testing correction (optional)
wilcoxon[]

# just one test
#Basic
wilcA<-data[Metric=="Precision" & type=="basic" & K==10 & Stat=="Averaged"]
wilcP<-data[Metric=="Precision" & type=="basic" & K==10 & Stat=="Pooled"]
wilcox.test(wilcA$Value,wilcP$Value,alternative = "two.sided",paired=TRUE)

#Stratified
wilcA<-data[Metric=="Precision" & type=="stratified" & K==10 & Stat=="Averaged"]
wilcP<-data[Metric=="Precision" & type=="stratified" & K==10 & Stat=="Pooled"]
wilcox.test(wilcA$Value,wilcP$Value,alternative = "two.sided",paired=TRUE)

# Does all the tests
wilcox_pvals <- function(data,
                         value_col = c("Value", "Var"),
                         metrics   = c("Misclassification error",
                                       "Specificity",
                                       "Precision",
                                       "AUC pROC"),
                         types     = c("basic", "stratified"),
                         K_value   = 10) {
  dt <- as.data.table(data)
  value_col <- match.arg(value_col)
  K_char <- as.character(K_value)
  
  sub <- dt[
    Metric %in% metrics &
      type   %in% types &
      K      == K_char &
      Stat   %in% c("Averaged", "Pooled")
  ]
  
  res <- sub[, {
    vals <- get(value_col)
    wilcA <- vals[Stat == "Averaged"]
    wilcP <- vals[Stat == "Pooled"]
    
    if (length(wilcA) == 0L || length(wilcP) == 0L) {
      .(column    = value_col,
        n_pairs   = NA_integer_,
        statistic = NA_real_,
        p_value   = NA_real_)
    } else {
      test <- wilcox.test(wilcA, wilcP,
                          alternative = "two.sided",
                          paired      = TRUE,
                          exact       = FALSE)
      .(column    = value_col,
        n_pairs   = length(wilcA),
        p_value   = round(unname(test$p.value), 7))
    }
  }, by = .(dataset, N, Metric, type, K)]
  res[]
}

wilcox <- wilcox_pvals(kloop_BIGBIG, value_col = "Value")

## Cleans data a bit, so we can plot p-values againt N
wilcox <- wilcox[!startsWith(dataset, "clean")]
wilcox[, prev :=
         fcase(
           startsWith(dataset, "sim_low"), "low",
           startsWith(dataset, "sim_mid"), "mid",
           startsWith(dataset, "sim_bal"), "bal",
           default = NA_character_)]
wilcox<-wilcox[Metric == "Precision" & type == "stratified" ] #& prev=="bal"
#wilcox<-wilcox[,p_value:=round(p_value, 9)]
## Plot p_values over N (need to subset the wilcox data set first to specify a metric, type and prevalence)
ggplot(wilcox, aes(x = N, y = p_value)) +
  geom_line(color="steelblue4") +
  geom_point(color="steelblue2") +
  geom_hline(yintercept = 0.05, linetype = "dashed") +
  labs(x = "Sample size",
       y = "p-value",
       title = "Precision" ) +
  theme_bw(base_size = 13) +
  theme(plot.title = element_text(hjust = 0.5))



wilcox[, prev_num := factor(
  fcase(
    prev == "low",  0.01,
    prev == "mid",  0.10,
    prev == "bal",  0.50,
    default = NA_real_
  ),
  levels = c(0.01, 0.10, 0.50),
  labels = c("0.01", "0.1", "0.5")
)]

size_labels <- c("Smallest", "Small", "Medium", "Large", "Largest")
wilcox[, size_label := {
  Ns <- sort(unique(N))
  factor(N,
         levels = Ns,
         labels = size_labels[seq_along(Ns)])
}, by = prev]

#plot imbalance against p-value (one line for N)
ggplot(wilcox, aes(x = prev_num, y = p_value,
                   color = size_label, group = size_label)) +
  geom_point(size = 3) +
  geom_line()+
  labs(y = "p-value",
       color = "N scenario") +
  theme_bw(base_size = 13)

#plot imbalance against effect size (one line for N)
ggplot(dl, aes(x = prev_num, y = mean_diff_est,
               color = size_label, group = size_label)) +
  geom_point(size = 3) +
  geom_line()+
  labs(y = "effect size",
       color = "N scenario") +
  theme_bw(base_size = 13)



######################################################
### Variance plots
dl <- dcast(
  dt_mid[Stat %in% c("Pooled","Averaged")],
  Metric + type + N + K + seed + dataset ~ Stat,
  value.var = "Var"
)
dl <- dl[!is.na(Pooled) & !is.na(Averaged) & Pooled > 0 & Averaged > 0]
dl[, `:=`(
  ratio_var = Pooled / Averaged
)]
dl[, relative_diff := (Pooled-Averaged)/Averaged]
dl[,diff := Averaged-Pooled]

##Boxplots over N instead of K, for a single splitting type - Variance
plot_box_by_N <- function(dt, metric_name, K_val, cv_type = c("basic", "stratified")) {
  cv_type <- match.arg(cv_type)

  type_label <- if (cv_type == "basic") "Random" else "Stratified"
  type_color <- if (cv_type == "basic") "#d1495b" else "#2e86ab"
  
  setDT(dt)
  d <- dt[
    Metric == metric_name &
      !is.na(relative_diff) &
      K == K_val &
      type == cv_type 
  ]

  d[, N := factor(N, levels = sort(unique(N)))]
  
  ggplot(d, aes(x = N, y = diff)) +
    #    geom_hline(yintercept = 0, linetype = "dashed") +
    geom_boxplot(alpha = 0.25,
                 width = 0.6,
                 fill = type_color,
                 colour = NA) +
    geom_jitter(width = 0, height = 0, size = 2, alpha = 0.7) +
    geom_hline(yintercept = 0, linetype = "dashed") +
    stat_summary(
      fun = mean,
      geom = "crossbar",
      width = 0.4,
      colour = type_color
    ) +
    #    coord_cartesian(ylim = c(-0.01, 0.01)) +   #sets specific y-axis
    labs(
      title = paste0("Severely imbalanced"),
      x = "N",
      y = "Difference"
    ) +
    theme_bw(base_size = 13) +
    theme(
      panel.grid.major.x = element_blank(),
      plot.title = element_text(hjust = 0.5)
    )
}
#Random
plots_byN <- setNames(
  lapply(metrics, function(m)
    plot_box_by_N(dl,
                  metric_name = m,
                  K_val = 10,
                  cv_type = "basic")),
  metrics
)

print(plots_byN[["Misclassification error"]])
print(plots_byN[["Specificity"]])
print(plots_byN[["Precision"]])
print(plots_byN[["AUC pROC"]])
#Stratified
plots_byN <- setNames(
  lapply(metrics, function(m)
    plot_box_by_N(dl,
                  metric_name = m,
                  K_val = 10,
                  cv_type = "stratified")),
  metrics
)
print(plots_byN[["Misclassification error"]])
print(plots_byN[["Specificity"]])
print(plots_byN[["Precision"]])
print(plots_byN[["AUC pROC"]])

#### MSE and Bias
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

###Plots the MSE
ggplot(
  dl,
  aes(x = size_label, y = MSE,
      color = col_id,
      linetype = MSE_kind,
      group = interaction(col_id, MSE_kind))
) +
  geom_line(linewidth = 0.9) +
  geom_point(size = 2) +
  facet_wrap(~Metric, ncol = 2, scales = "free_y") +
  scale_color_manual(values = c("Average.0.01"="brown2","Pool.0.01"="skyblue","Average.0.1"="indianred","Pool.0.1"="steelblue","Average.0.5"="tomato4","Pool.0.5"="deepskyblue4"), 
                     labels = w,
                     name = "") +
  scale_linetype_manual(values = c(Average = "solid", Pool = "solid"), guide="none",name = NULL) +
  labs(x = "", y = "MSE") +
  theme_bw() +
  theme(legend.position = "bottom",legend.key.width = grid::unit(2.2, "cm") )

#Bias

w<-dl[,c("Metric", "type", "N", "prev", "prev_num", "size_label","Value_Averaged_mean", "Value_Pooled_mean","diff_est_mean","diff_est_q025","diff_est_q975", "est.err_Averaged_mean", "est.err_Pooled_mean")]
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

dl <- melt(
  dl,
  measure.vars  = c("est.err_Averaged_mean", "est.err_Pooled_mean"),
  variable.name = "BIAS_kind",
  value.name    = "bias"
)

dl[, BIAS_kind := fifelse(BIAS_kind == "est.err_Averaged_mean", "Average", "Pool")]
dl[, size_label := factor(size_label, levels = c("Smallest", "Small"  ,  "Medium"  , "Large"  ,  "Largest"))]
dl[, col_id := interaction(BIAS_kind, prev_num, drop = TRUE)]
w<- c("Average.0.01"=paste0("Average, ", "\u03C0","=0.01"),"Pool.0.01"=paste0("Pool, ", "\u03C0","=0.01"), "Average.0.1"=paste0("Average, ", "\u03C0","=0.1"), "Pool.0.1"=paste0("Pool, ", "\u03C0","=0.1") ,"Average.0.5"=paste0("Average, ", "\u03C0","=0.5") ,"Pool.0.5"=paste0("Pool, ", "\u03C0","=0.5"))

###Plots the BIAS
ggplot(
  dl,
  aes(x = size_label, y = bias,
      color = col_id,
      linetype = type,
      group = interaction(col_id, type))
) +
  geom_line(linewidth = 0.9) +
  geom_point(size = 2) +
  facet_wrap(~Metric, ncol = 2, scales = "free_y") +
  scale_color_manual(values = c("Average.0.01"="brown2","Pool.0.01"="skyblue","Average.0.1"="indianred","Pool.0.1"="steelblue","Average.0.5"="tomato4","Pool.0.5"="deepskyblue4"), 
                     labels = w,
                     name = "") +
  guides(
    linetype = guide_legend(ncol = 1, byrow = TRUE,
                            override.aes = list(color = "black"))) +
  scale_linetype_manual(values = c(basic = "solid", stratified = "dashed"),labels = c("basic"="Random", "stratified"="Stratified"), name = NULL) +
  labs(x = "", y = "BIAS") +
  theme_bw() +
  theme(legend.position = "bottom",legend.key.width = grid::unit(2.2, "cm") )

