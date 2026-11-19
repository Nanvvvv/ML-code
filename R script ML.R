# TMJ disc perforation: complete analysis and figures
# Revised manuscript workflow, 19 November 2026.
# Open this script in RStudio and use Source, or run it with Rscript.
# Keep the input and asset folders beside the script.
# The CSV files are the reconciled, deidentified analysis datasets.
# A full run repeats model development, 200 selection bootstraps and all plots.

if (.Platform$OS.type == "windows") invisible(Sys.setlocale("LC_CTYPE", ".UTF-8"))
project_dir <- "C:/Users/lnj82/Desktop/可仿照课题设计/诊断模型/终版投joor/修回/最终上传/终稿整理_20261118/完整R分析_按终稿整理_20261119"
input_dir <- file.path(project_dir, "输入数据")
result_dir <- file.path(project_dir, "结果")
figure_dir <- file.path(project_dir, "图")
dir.create(result_dir, showWarnings = FALSE)
dir.create(figure_dir, showWarnings = FALSE)

seed <- 20261101L
sensitivity_target <- 0.9
models <- c("LR", "RF", "DT", "SVM", "KNN", "MLP", "XGB", "LGBM")
shap_model <- "XGB"
alphas <- c(0.25, 0.5, 0.75, 1)
lambda_grid <- 10^seq(0, -5, length.out = 80)
packages <- c("digest", "glmnet", "randomForest", "rpart", "e1071", "nnet", "xgboost", "lightgbm", "sandwich",
  "ggplot2", "patchwork", "ragg", "fastshap", "png")
missing <- packages[!vapply(packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) stop("Install these packages first: ", paste(missing, collapse = ", "))

# 1. Read and check the analysis datasets ---------------------------------------
predictors <- c("Age", "Duration_since_first_symptom", "Sex", "No_symptom", "TMJ_noise", "Arthralgia",
  "Limited_mouth_opening", "Lesion_laterality", "Morphology_of_disc", "Fluid_location", "Condylar_fossa_contact",
  "Normal", "Erosion", "Sclerosis", "Pseudocysts", "Osteophyte")
continuous <- predictors[1:2]
levels_by_feature <- setNames(rep(list(0:1), 14), predictors[-(1:2)])
levels_by_feature$Morphology_of_disc <- 0:4
levels_by_feature$Fluid_location <- 0:2
EN_NAMES <- c(Age = "Age", Duration_since_first_symptom = "Symptom duration", Sex = "Sex", No_symptom = "Asymptomatic",
  TMJ_noise = "TMJ noise", Arthralgia = "Arthralgia", Limited_mouth_opening = "Limited mouth opening",
  Lesion_laterality = "Lesion laterality", Morphology_of_disc = "Disc morphology", Fluid_location = "Fluid location",
  Condylar_fossa_contact = "Condylar–fossa contact", Normal = "Normal condylar morphology", Erosion = "Erosion",
  Sclerosis = "Sclerosis", Pseudocysts = "Pseudocysts", Osteophyte = "Osteophyte")
LEVEL_NAMES <- list(Sex = c("Male", "Female"), No_symptom = c("No", "Yes"), TMJ_noise = c("No", "Yes"),
  Arthralgia = c("No", "Yes"), Limited_mouth_opening = c("No", "Yes"), Lesion_laterality = c("Unilateral",
    "Bilateral"), Morphology_of_disc = c("Biconcave", "Folded", "Flattened", "Eyeglass-shaped", "Amorphous"),
  Fluid_location = c("None", "Single cavity", "Multiple cavities"), Condylar_fossa_contact = c("No",
    "Yes"), Normal = c("No", "Yes"), Erosion = c("No", "Yes"), Sclerosis = c("No", "Yes"), Pseudocysts = c("No",
    "Yes"), Osteophyte = c("No", "Yes"))
fmt_p <- function(p) ifelse(is.na(p), "NE", ifelse(p < 0.001, "<0.001", sprintf("%.3f", p)))

development <- read.csv(file.path(input_dir, "开发队列.csv"), fileEncoding = "UTF-8")
temporal <- read.csv(file.path(input_dir, "本院时间队列.csv"), fileEncoding = "UTF-8")
geographic <- read.csv(file.path(input_dir, "外院队列.csv"), fileEncoding = "UTF-8")
# patient_id keeps bilateral joints together; row_id identifies the individual joint.
# Original IDs already include cohort/centre namespaces where needed.
for (data in list(development, temporal, geographic)) {
  stopifnot(all(c(predictors, "Perforation", "patient_id", "row_id", "Center") %in% names(data)), !anyNA(data),
    !anyDuplicated(data$row_id), all(data$Perforation %in% 0:1), max(table(data$patient_id)) <= 2L)
  for (variable in names(levels_by_feature)) {
    stopifnot(all(data[[variable]] %in% levels_by_feature[[variable]]))
  }
}
stopifnot(nrow(development) == 306L, length(unique(development$patient_id)) == 216L, sum(development$Perforation) ==
  102L, nrow(temporal) == 141L, length(unique(temporal$patient_id)) == 104L, sum(temporal$Perforation) ==
  24L, nrow(geographic) == 83L, length(unique(geographic$patient_id)) == 76L, sum(geographic$Perforation) ==
  29L)



# 2.5–2.6. Preprocessing, feature selection and fitting operations ----------------
# A few named operations are reused across algorithms and folds.
# No imputation is needed in these complete analysis datasets.
safe_div <- function(a, b) {
  if (b > 0)
    a/b else NA_real_
}

# Keep both joints from each patient in the same outcome-stratified fold.
patient_folds <- function(d, v = 5L, seed = 20261101L) {
  set.seed(seed)
  outcome_by_patient <- tapply(d$Perforation, d$patient_id, max)
  if (length(outcome_by_patient) < v)
    stop("Too few patients for cross-validation.")
  for (attempt in seq_len(100L)) {
    assignment <- setNames(integer(length(outcome_by_patient)), names(outcome_by_patient))
    for (outcome in sort(unique(outcome_by_patient))) {
      ids <- names(outcome_by_patient)[outcome_by_patient == outcome]
      ids <- ids[sample.int(length(ids))]
      assignment[ids] <- rep(sample.int(v), length.out = length(ids))
    }
    folds <- unname(assignment[d$patient_id])
    valid <- logical(v)
    for (k in seq_len(v)) {
      valid[k] <- length(unique(d$Perforation[folds == k])) == 2 && length(unique(d$Perforation[folds !=
        k])) == 2
    }
    if (all(valid))
      break
  }
  if (!all(valid))
    stop("A fold has only one outcome class.")
  for (k in seq_len(v)) {
    overlap <- intersect(d$patient_id[folds == k], d$patient_id[folds != k])
    if (length(overlap))
      stop("A patient appears in both parts of a split.")
  }
  folds
}

# Keep the established random splits reproducible for each training subset.
seed_for <- function(d, salt = 0L) {
  training_ids <- sort(d$row_id)
  fingerprint <- digest::digest(training_ids, algo = "xxhash32")
  number <- strtoi(substr(fingerprint, 1, 7), base = 16L)
  as.integer(number + salt)
}

# Create dummy columns using the recorded clinical categories.
fs_matrix <- function(d, imputation) {
  z <- d[, predictors, drop = FALSE]
  for (v in predictors) z[[v]][is.na(z[[v]])] <- imputation[[v]]
  for (v in names(levels_by_feature)) z[[v]] <- factor(z[[v]], levels = levels_by_feature[[v]])
  tt <- terms(reformulate(predictors))
  mm <- model.matrix(tt, z)
  a <- attr(mm, "assign")[-1]
  mm <- mm[, colnames(mm) != "(Intercept)", drop = FALSE]
  attr(mm, "feature") <- predictors[a]
  mm
}

# Estimate preprocessing in training data; exclude sparse sclerosis and constant columns.
fs_prep <- function(d) {
  imputation <- list()
  for (variable in predictors) {
    values <- na.omit(d[[variable]])
    if (!length(values))
      stop("Completely missing training variable: ", variable)
    if (variable %in% continuous) {
      imputation[[variable]] <- median(values)
    } else {
      counts <- sort(table(values), decreasing = TRUE)
      imputation[[variable]] <- as.numeric(names(counts)[1])
    }
  }
  x <- fs_matrix(d, imputation)
  column_sd <- apply(x, 2, sd)
  feature <- attr(x, "feature")
  keep <- is.finite(column_sd) & column_sd > 0 & feature != "Sclerosis"
  list(imputation = imputation, columns = colnames(x)[keep], feature = setNames(feature[keep], colnames(x)[keep]),
    center = colMeans(x[, keep, drop = FALSE]), scale = column_sd[keep], missing = colSums(is.na(d[,
      predictors])))
}

# Apply a training transformation to the requested rows.
fs_x <- function(d, pp) {
  x <- fs_matrix(d, pp$imputation)[, pp$columns, drop = FALSE]
  centred <- sweep(x, 2, pp$center, "-")
  sweep(centred, 2, pp$scale, "/")
}

# Fit the full elastic-net penalty path, checking convergence.
fs_path <- function(d, pp, alpha, lambda = lambda_grid) {
  ctl <- glmnet::glmnet.control()
  on.exit(do.call(glmnet::glmnet.control, ctl))
  glmnet::glmnet.control(fdev = 0, devmax = 1)
  for (method in c("Newton", "modified.Newton")) {
    f <- suppressWarnings(glmnet::glmnet(fs_x(d, pp), d$Perforation, family = "binomial", alpha = alpha,
      lambda = lambda, standardize = FALSE, thresh = 1e-07, maxit = 10000000L, type.logistic = method))
    if (length(f$lambda) == length(lambda) && f$jerr == 0)
      return(f)
  }
  stop("The full penalty path did not converge.")
}

# Choose alpha by grouped log loss and lambda by the one-standard-error rule.
select_features <- function(d, random_seed) {
  folds <- patient_folds(d, 5L, random_seed)
  preps <- vector("list", 5L)
  for (k in 1:5) preps[[k]] <- fs_prep(d[folds != k, ])
  scores <- list()
  for (a in alphas) {
    loss <- matrix(NA_real_, length(lambda_grid), 5L)
    for (k in 1:5) {
      train <- folds != k
      fit <- fs_path(d[train, ], preps[[k]], a)
      x <- fs_x(d[!train, ], preps[[k]])
      probability <- as.matrix(predict(fit, newx = x, s = lambda_grid, type = "response"))
      y <- d$Perforation[!train]
      log_loss <- -y * log(pmax(probability, 1e-12)) - (1 - y) * log(pmax(1 - probability, 1e-12))
      loss[, k] <- colMeans(log_loss)
    }
    scores[[as.character(a)]] <- data.frame(alpha = a, lambda = lambda_grid, loss = rowMeans(loss),
      se = apply(loss, 1, sd)/sqrt(5L))
  }
  scores <- do.call(rbind, scores)
  best <- scores[order(scores$loss, -scores$alpha, -scores$lambda)[1], ]
  eligible <- scores[scores$alpha == best$alpha & scores$loss <= best$loss + best$se, ]
  one <- eligible[which.max(eligible$lambda), ]
  prep <- fs_prep(d)
  fit <- fs_path(d, prep, one$alpha)
  beta <- as.numeric(coef(fit, s = one$lambda))[-1]
  selected <- intersect(predictors, unname(prep$feature[abs(beta) > 1e-08]))
  list(variables = selected, alpha = one$alpha, lambda = one$lambda, prep = prep)
}

# Apply the selection rule and retain all columns of each selected original variable.
prep_local <- function(d, vars = predictors, strategy = "elastic_net") {
  if (strategy == "elastic_net") {
    s <- if (nrow(d) == 306 && length(unique(d$patient_id)) == 216)
      seed + 500000L else seed_for(d, 117L)
    picked <- select_features(d, s)
    pp <- picked$prep
    vars <- intersect(picked$variables, vars)
    pp$alpha <- picked$alpha
    pp$lambda <- picked$lambda
  } else pp <- fs_prep(d)
  pp$columns <- pp$columns[pp$feature[pp$columns] %in% vars]
  pp$center <- pp$center[pp$columns]
  pp$scale <- pp$scale[pp$columns]
  pp$variables <- vars
  pp
}

# Keep probabilities inside the finite logit range.
clip_prob <- function(p) pmin(pmax(as.numeric(p), 1e-07), 1 - 1e-07)

# Rank AUC with perforation fixed as the positive outcome.
auc_fixed <- function(y, p) {
  if (length(unique(y)) < 2)
    return(NA_real_)
  n1 <- sum(y == 1)
  n0 <- sum(y == 0)
  (sum(rank(p, ties.method = "average")[y == 1]) - n1 * (n1 + 1)/2)/(n1 * n0)
}

# Candidate parameter settings used in Table S2.
model_grids <- function() {
  g <- list(LR = data.frame(method = "MLE"), RF = expand.grid(mtry_frac = c(0.25, 0.5, 1), nodesize = c(3L,
    8L)), DT = expand.grid(cp = c(0.001, 0.01, 0.05), maxdepth = c(2L, 4L, 6L)), SVM = expand.grid(cost = c(0.1,
    1, 10), gamma = c(0.01, 0.1)), KNN = expand.grid(k = c(5L, 9L, 15L, 21L, 31L)), MLP = expand.grid(size = c(2L,
    4L, 6L), decay = c(0.01, 0.1, 1)), XGB = expand.grid(max_depth = c(1L, 2L, 3L), eta = c(0.03,
    0.1), min_child_weight = c(3, 8)), LGBM = expand.grid(num_leaves = c(3L, 7L), learning_rate = c(0.03,
    0.1), min_data_in_leaf = c(10L, 20L)))
  g
}

# Record ordinary-LR convergence warnings and any aliased coefficients.
fit_logistic_mle <- function(x, y) {
  xx <- cbind(`(Intercept)` = 1, as.matrix(x))
  warnings <- character()
  f <- withCallingHandlers(stats::glm.fit(x = xx, y = y, family = stats::binomial(), control = stats::glm.control(maxit = 100L,
    epsilon = 1e-08)), warning = function(w) {
    warnings <<- c(warnings, conditionMessage(w))
    invokeRestart("muffleWarning")
  })
  list(coefficients = f$coefficients, rank = f$rank, converged = f$converged, boundary = f$boundary,
    warnings = unique(warnings), aliased = names(f$coefficients)[is.na(f$coefficients)], method = "unpenalised maximum likelihood")
}

# Predict with estimable ordinary-LR coefficients.
predict_logistic_mle <- function(f, x) {
  b <- f$coefficients
  b[is.na(b)] <- 0
  clip_prob(plogis(as.numeric(cbind(1, as.matrix(x)) %*% b)))
}

# Fit one algorithm at one candidate setting; LR is ordinary maximum likelihood.
fit_engine <- function(model, x, y, h, seed) {
  set.seed(seed)
  yf <- factor(y, levels = 0:1)
  if (ncol(x) < 1)
    return(list(intercept_only = mean(y)))
  if (model == "LR")
    return(fit_logistic_mle(x, y))
  if (model == "RF")
    return(randomForest::randomForest(x = x, y = yf, ntree = 300L, mtry = max(1L, round(ncol(x) *
      h$mtry_frac)), nodesize = h$nodesize))
  if (model == "DT") {
    dd <- data.frame(y = yf, x, check.names = FALSE)
    return(rpart::rpart(y ~ ., dd, method = "class", control = rpart::rpart.control(cp = h$cp, maxdepth = h$maxdepth,
      minsplit = 10, xval = 0)))
  }
  if (model == "SVM")
    return(e1071::svm(x = x, y = yf, kernel = "radial", cost = h$cost, gamma = h$gamma, scale = FALSE,
      probability = TRUE))
  if (model == "KNN")
    return(list(x = x, y = y, k = min(h$k, nrow(x))))
  if (model == "MLP")
    return(nnet::nnet(x = x, y = y, size = h$size, decay = h$decay, entropy = TRUE, maxit = 500,
      MaxNWts = 10000, trace = FALSE))
  if (model == "XGB")
    return(xgboost::xgb.train(params = list(objective = "binary:logistic", eval_metric = "logloss",
      max_depth = h$max_depth, eta = h$eta, min_child_weight = h$min_child_weight, subsample = 0.8,
      colsample_bytree = 0.8, nthread = 1L, seed = seed), data = xgboost::xgb.DMatrix(x, label = y),
      nrounds = 150L, verbose = 0))
  if (model == "LGBM")
    return(lightgbm::lgb.train(params = list(objective = "binary", metric = "binary_logloss", num_leaves = h$num_leaves,
      learning_rate = h$learning_rate, min_data_in_leaf = h$min_data_in_leaf, feature_pre_filter = FALSE,
      verbosity = -1L, num_threads = 1L, seed = seed, deterministic = TRUE, force_col_wise = TRUE),
      data = lightgbm::lgb.Dataset(x, label = y), nrounds = 150L, verbose = -1L))
  stop("Unknown algorithm: ", model)
}

# Return probabilities of perforation for one fitted algorithm.
predict_engine <- function(model, fit, x, h) {
  if (is.list(fit) && !is.null(fit$intercept_only)) {
    return(rep(fit$intercept_only, nrow(x)))
  }
  if (model == "LR")
    p <- predict_logistic_mle(fit, x)
  if (model == "RF")
    p <- predict(fit, x, type = "prob")[, "1"]
  if (model == "DT")
    p <- predict(fit, data.frame(x, check.names = FALSE), type = "prob")[, "1"]
  if (model == "SVM") {
    predicted <- predict(fit, x, probability = TRUE)
    p <- attr(predicted, "probabilities")[, "1"]
  }
  if (model == "KNN") {
    p <- numeric(nrow(x))
    for (i in seq_len(nrow(x))) {
      difference <- sweep(fit$x, 2, x[i, ], "-")
      distance <- rowSums(difference^2)
      boundary <- sort(distance, partial = fit$k)[fit$k]
      p[i] <- mean(fit$y[distance <= boundary])
    }
  }
  if (model == "MLP")
    p <- as.numeric(predict(fit, x, type = "raw"))
  if (model == "XGB")
    p <- as.numeric(predict(fit, xgboost::xgb.DMatrix(x)))
  if (model == "LGBM")
    p <- as.numeric(predict(fit, x))
  if (length(p) != nrow(x) || any(!is.finite(p)))
    stop("Invalid prediction from ", model)
  clip_prob(p)
}

# Store one fitted algorithm with its preprocessing and parameters.
fit_one <- function(d, m, h, pp, seed) {
  x <- fs_x(d, pp)
  fit <- fit_engine(m, x, d$Perforation, h, seed)
  list(model = m, h = h, prep = pp, fit = fit)
}

# Apply stored preprocessing before predicting new rows.
predict_one <- function(f, d) {
  x <- fs_x(d, f$prep)
  predict_engine(f$model, f$fit, x, f$h)
}

# Tune each algorithm by pooled inner-holdout AUC; refit on the supplied training data.
tune_library <- function(d, models_to_fit = models, vars = predictors, strategy = "elastic_net") {
  folds <- patient_folds(d, 5L, seed_for(d, 219L))
  prep <- vector("list", 5L)
  for (k in 1:5) {
    prep[[k]] <- prep_local(d[folds != k, , drop = FALSE], vars, strategy)
  }
  full_prep <- prep_local(d, vars, strategy)
  grids <- model_grids()
  fits <- scores <- list()
  for (model in models_to_fit) {
    grid <- grids[[model]]
    score <- data.frame(grid, AUC = NA_real_, Brier = NA_real_)
    for (setting in seq_len(nrow(grid))) {
      probability <- numeric(nrow(d))
      parameters <- as.list(grid[setting, , drop = FALSE])
      for (k in 1:5) {
        train <- folds != k
        fit <- fit_one(d[train, , drop = FALSE], model, parameters, prep[[k]], seed_for(d, 301L +
          k))
        probability[!train] <- predict_one(fit, d[!train, , drop = FALSE])
      }
      score$AUC[setting] <- auc_fixed(d$Perforation, probability)
      score$Brier[setting] <- mean((d$Perforation - probability)^2)
    }
    best <- order(-score$AUC, score$Brier, seq_len(nrow(grid)))[1]
    parameters <- as.list(grid[best, , drop = FALSE])
    fits[[model]] <- fit_one(d, model, parameters, full_prep, seed_for(d, 401L))
    scores[[model]] <- score
  }
  list(fits = fits, scores = scores, fold = folds, variables = full_prep$variables, prep = full_prep)
}

# Predict all fitted algorithms on the same rows.
predict_library <- function(b, d) {
  probabilities <- matrix(NA_real_, nrow(d), length(b$fits))
  colnames(probabilities) <- names(b$fits)
  for (model in names(b$fits)) probabilities[, model] <- predict_one(b$fits[[model]], d)
  probabilities
}

# Youden threshold for the supplementary operating-point comparison.
choose_cutoff <- function(y, p) {
  cand <- sort(unique(c(0, p, 1)))
  j <- vapply(cand, function(t) mean(p[y == 1] >= t) + mean(p[y == 0] < t) - 1, numeric(1))
  max(cand[j == max(j)])
}

# Highest threshold reaching the requested sensitivity in development predictions.
target_cut <- function(y, p, target) {
  candidates <- sort(unique(c(0, p, 1)))
  sensitivity <- numeric(length(candidates))
  for (i in seq_along(candidates)) sensitivity[i] <- mean(p[y == 1] >= candidates[i])
  max(candidates[sensitivity >= target])
}

# Metrics at a fixed operating threshold.
metric_values <- function(y, p, cutoff = 0.5, calibration = TRUE) {
  p <- clip_prob(p)
  yes <- p >= cutoff
  tp <- sum(yes & y == 1)
  tn <- sum(!yes & y == 0)
  fp <- sum(yes & y == 0)
  fn <- sum(!yes & y == 1)
  o <- order(p, decreasing = TRUE)
  ys <- y[o]
  ps <- p[o]
  j <- which(!duplicated(ps, fromLast = TRUE))
  cum <- cumsum(ys)[j]
  ap <- if (sum(y))
    sum(cum/j * diff(c(0, cum)))/sum(y) else NA_real_
  c1 <- c0 <- cs <- NA_real_
  if (calibration && sum(y == 1) >= 5 && sum(y == 0) >= 5 && sd(p) > 1e-08) {
    lp <- qlogis(p)
    f <- tryCatch(suppressWarnings(glm(y ~ lp, family = binomial(), control = glm.control(maxit = 100))),
      error = function(e) NULL)
    if (!is.null(f) && f$converged && all(abs(coef(f)) < 100)) {
      c0 <- unname(coef(f)[1])
      cs <- unname(coef(f)[2])
    }
    g <- tryCatch(suppressWarnings(glm(y ~ 1 + offset(lp), family = binomial())), error = function(e) NULL)
    if (!is.null(g) && g$converged)
      c1 <- unname(coef(g)[1])
  }
  c(AUC = auc_fixed(y, p), Average_precision = ap, Brier = mean((y - p)^2), Accuracy = mean(yes ==
    y), Sensitivity = safe_div(tp, tp + fn), Specificity = safe_div(tn, tn + fp), PPV = safe_div(tp,
    tp + fp), NPV = safe_div(tn, tn + fn), F1 = safe_div(2 * tp, 2 * tp + fp + fn), Calibration_in_large = c1,
    Calibration_intercept = c0, Calibration_slope = cs)
}

# Score each outer holdout before taking equal-weight means.
fold_scores <- function(d, p, cut, fold) {
  rows <- list()
  keys <- c("AUC", "Accuracy", "Sensitivity", "Specificity", "F1", "Brier")
  for (k in sort(unique(fold))) {
    i <- which(fold == k)
    for (model in colnames(p)) {
      metrics <- metric_values(d$Perforation[i], p[i, model], cut[i, model], FALSE)
      rows[[length(rows) + 1L]] <- data.frame(model = model, fold = k, as.list(metrics[keys]),
        Threshold = unname(cut[i[1], model]), row.names = NULL)
    }
  }
  do.call(rbind, rows)
}

# Sample patients with replacement, retaining all their joint rows.
cluster_draw <- function(d) {
  groups <- split(seq_len(nrow(d)), d$patient_id)
  ids <- sample.int(length(groups), length(groups), replace = TRUE)
  unlist(groups[ids], use.names = FALSE)
}

# Percentile interval from valid bootstrap estimates.
interval <- function(x) unname(quantile(x[is.finite(x)], c(0.025, 0.975)))

# Patient-block confidence intervals conditional on the supplied predictions.
cluster_ci <- function(d, p, cutoff = 0.5, B = 2000L, seed = 20261101L, calibration = TRUE) {
  set.seed(seed)
  point <- metric_values(d$Perforation, p, cutoff, calibration)
  boot <- matrix(NA_real_, B, length(point), dimnames = list(NULL, names(point)))
  for (b in seq_len(B)) {
    i <- cluster_draw(d)
    if (length(unique(d$Perforation[i])) < 2)
      next
    ct <- if (length(cutoff) > 1)
      cutoff[i] else cutoff
    boot[b, ] <- metric_values(d$Perforation[i], p[i], ct, calibration)
  }
  ci <- apply(boot, 2, function(x) if (sum(is.finite(x)) >= 20)
    quantile(x[is.finite(x)], c(0.025, 0.975)) else c(NA, NA))
  data.frame(metric = names(point), estimate = unname(point), lower = ci[1, ], upper = ci[2, ], successful_bootstraps = colSums(is.finite(boot)),
    row.names = NULL)
}

# Matched patient-block draws; preserve development fold membership.
bootstrap_rows <- function(d, fold, B = 2000L) {
  set.seed(seed)
  groups <- split(seq_len(nrow(d)), fold)
  draws <- vector("list", length(groups))
  names(draws) <- names(groups)
  for (k in seq_along(groups)) {
    i <- groups[[k]]
    patients <- split(i, d$patient_id[i])
    draws[[k]] <- vector("list", B)
    for (b in seq_len(B)) {
      sampled <- sample.int(length(patients), replace = TRUE)
      draws[[k]][[b]] <- unlist(patients[sampled], use.names = FALSE)
    }
  }
  draws
}

# Average development fold AUCs within each bootstrap draw.
bootstrap_auc <- function(y, p, draws) {
  total <- numeric(length(draws[[1]]))
  for (fold in draws) {
    values <- numeric(length(fold))
    for (b in seq_along(fold)) {
      i <- fold[[b]]
      values[b] <- auc_fixed(y[i], p[i])
    }
    total <- total + values
  }
  total/length(draws)
}

# Paired Brier differences using the same sampled patients.
bootstrap_brier_delta <- function(y, p1, p2, draws) {
  squared_error_difference <- (p1 - y)^2 - (p2 - y)^2
  total <- numeric(length(draws[[1]]))
  for (fold in draws) {
    values <- numeric(length(fold))
    for (b in seq_along(fold)) values[b] <- mean(squared_error_difference[fold[[b]]])
    total <- total + values
  }
  interval(total/length(draws))
}

# Ordinary LR with the eligible variables from one input group.
fit_modality <- function(d, vars) {
  prep <- prep_local(d, vars, "all_candidates")
  fit <- fit_one(d, "LR", list(method = "MLE"), prep, seed_for(d, 401L))
  list(fits = list(LR = fit))
}

# Patient-clustered logistic Wald test; sparse comparisons are not estimable.
cluster_wald <- function(d, v) {
  dd <- d[, v, drop = FALSE]
  dd$y <- d$Perforation
  if (v %in% names(LEVEL_NAMES))
    dd[[v]] <- factor(dd[[v]], levels = levels_by_feature[[v]])
  if (length(unique(dd[[v]])) < 2 || length(unique(dd$y)) < 2)
    return(NA_real_)
  if (v %in% names(LEVEL_NAMES)) {
    dd[[v]] <- droplevels(dd[[v]])
    tab <- table(dd[[v]], dd$y)
    if (any(tab == 0) || any(rowSums(tab) < 5))
      return(NA_real_)
  }
  tryCatch({
    f <- suppressWarnings(glm(reformulate(v, "y"), dd, family = binomial()))
    b <- coef(f)[-1]
    vc <- sandwich::vcovCL(f, cluster = d$patient_id, type = "HC1")[-1, -1, drop = FALSE]
    if (!f$converged || any(!is.finite(b)) || qr(vc)$rank < length(b))
      return(NA_real_)
    pchisq(as.numeric(t(b) %*% solve(vc, b)), df = length(b), lower.tail = FALSE)
  }, error = function(e) NA_real_)
}

# Median and IQR, or count and percentage, for a table row.
describe_value <- function(x, variable, level = NULL) {
  if (!length(x))
    return("NE")
  if (variable %in% continuous) {
    quantiles <- quantile(x, c(0.25, 0.5, 0.75))
    return(sprintf("%.2f [%.2f, %.2f]", quantiles[2], quantiles[1], quantiles[3]))
  }
  sprintf("%d (%.1f%%)", sum(x == level), 100 * mean(x == level))
}

# Joint-level descriptive summaries, with patient-clustered comparisons.
baseline_table <- function(d) {
  rows <- list(data.frame(variable = "N", level = "", Characteristic = "Joints", Overall = as.character(nrow(d)),
    No_perforation = as.character(sum(d$Perforation == 0)), Perforation = as.character(sum(d$Perforation ==
      1)), P_value = ""))
  for (variable in predictors) {
    p <- fmt_p(cluster_wald(d, variable))
    if (variable %in% continuous) {
      unit <- if (variable == "Age")
        ", years" else ", months"
      rows[[length(rows) + 1L]] <- data.frame(variable = variable, level = "", Characteristic = paste0(EN_NAMES[[variable]],
        unit), Overall = describe_value(d[[variable]], variable), No_perforation = describe_value(d[d$Perforation ==
        0, variable], variable), Perforation = describe_value(d[d$Perforation == 1, variable],
        variable), P_value = p)
    } else {
      rows[[length(rows) + 1L]] <- data.frame(variable = variable, level = "", Characteristic = EN_NAMES[[variable]],
        Overall = "", No_perforation = "", Perforation = "", P_value = p)
      for (level in levels_by_feature[[variable]]) {
        rows[[length(rows) + 1L]] <- data.frame(variable = variable, level = level, Characteristic = paste0("   ",
          LEVEL_NAMES[[variable]][level + 1]), Overall = describe_value(d[[variable]], variable,
          level), No_perforation = describe_value(d[d$Perforation == 0, variable], variable,
          level), Perforation = describe_value(d[d$Perforation == 1, variable], variable, level),
          P_value = "")
      }
    }
  }
  do.call(rbind, rows)
}

# Observed ROC coordinates for the stored probabilities.
roc_data <- function(y, p, model) {
  t <- c(Inf, sort(unique(p), decreasing = TRUE), -Inf)
  data.frame(model = model, FPR = vapply(t, function(a) mean(p[y == 0] >= a), numeric(1)), TPR = vapply(t,
    function(a) mean(p[y == 1] >= a), numeric(1)))
}

# Prediction wrapper used by fastshap with the development background.
predict_shap <- function(object, newdata) {
  for (variable in setdiff(predictors, names(newdata))) {
    newdata[[variable]] <- object$prep$imputation[[variable]]
  }
  predict_one(object, newdata)
}

MODALITIES <- list(clinical = c(continuous, "Sex", "No_symptom", "TMJ_noise", "Arthralgia", "Limited_mouth_opening"),
  MRI = c("Lesion_laterality", "Morphology_of_disc", "Fluid_location", "Condylar_fossa_contact"), CBCT = c("Normal",
    "Erosion", "Sclerosis", "Pseudocysts", "Osteophyte"))
MODALITIES$imaging <- c(MODALITIES$MRI, MODALITIES$CBCT)
MODALITIES$full <- predictors

# 2.6. Patient-grouped nested cross-validation -------------------------------------
# Each outer holdout is predicted by a fit that has never used its patients.
# inner CV: compare hyperparameters using pooled inner-holdout AUC.
# feature screening: rerun elastic-net tuning inside each training partition.
# An additional cross-fitting run inside each outer training set is
# used to determine that fold's operating threshold without holdout labels.
outer <- patient_folds(development, 5L, seed)
history <- vector("list", 5L)
outer_probability <- outer_cutoff <- matrix(NA_real_, nrow(development), length(models), dimnames = list(NULL,
  models))
for (k in 1:5) {
  message("Outer fold ", k, " of 5")
  train <- development[outer != k, , drop = FALSE]
  holdout <- development[outer == k, , drop = FALSE]
  fitted_models <- tune_library(train)
  heldout_probability <- predict_library(fitted_models, holdout)

  # Derive this outer fold's thresholds without its holdout outcomes.
  # Each training cross-fit repeats screening and algorithm tuning.
  training_folds <- patient_folds(train, 5L, seed_for(train, 719L))
  training_probability <- matrix(NA_real_, nrow(train), length(models), dimnames = list(NULL, models))
  for (j in 1:5) {
    training_fit <- tune_library(train[training_folds != j, , drop = FALSE])
    rows <- training_folds == j
    training_probability[rows, ] <- predict_library(training_fit, train[rows, , drop = FALSE])
  }
  rows <- which(outer == k)
  history[[k]] <- list(rows = rows, p = heldout_probability, train_p = training_probability, variables = fitted_models$variables)
  outer_probability[rows, ] <- heldout_probability
  for (model in models) {
    outer_cutoff[rows, model] <- target_cut(train$Perforation, training_probability[, model], sensitivity_target)
  }
}
stopifnot(all(is.finite(outer_probability)), all(is.finite(outer_cutoff)))
fold_results <- fold_scores(development, outer_probability, outer_cutoff, outer)
development_summary <- aggregate(fold_results[c("AUC", "Accuracy", "Sensitivity", "Specificity", "F1",
  "Brier")], list(model = fold_results$model), mean)
# Sorting is for presentation; all eight algorithms continue to assessment.
development_summary <- development_summary[order(-development_summary$AUC), ]
outer_features <- data.frame(fold = 1:5, features = vapply(history, function(x) paste(x$variables, collapse = "; "),
  character(1)))
write.csv(fold_results, file.path(result_dir, "各外层折模型指标.csv"), row.names = FALSE)
write.csv(development_summary, file.path(result_dir, "开发队列五折平均指标.csv"), row.names = FALSE)
write.csv(outer_features, file.path(result_dir, "五个外层训练集入选特征.csv"), row.names = FALSE)

# 2.6. Full-development tuning and final fits -----------
# Use the SAME grid and tuning criterion, but estimate a new setting using all
# development patients. Outer-fold parameter values are not averaged.
# This yields one deployable fitted object per algorithm, not five partial fits.
final_models <- tune_library(development)
final_cutoffs <- apply(outer_probability, 2, target_cut, y = development$Perforation, target = sensitivity_target)
write.csv(data.frame(variable = final_models$variables), file.path(result_dir, "最终入模特征.csv"),
  row.names = FALSE)
write.csv(data.frame(alpha = final_models$prep$alpha, lambda = final_models$prep$lambda, selected_variables = length(final_models$variables)),
  file.path(result_dir, "最终特征选择参数.csv"), row.names = FALSE)
parameter_rows <- list()
for (model in models) {
  grid <- final_models$scores[[model]]
  best <- order(-grid$AUC, grid$Brier, seq_len(nrow(grid)))[1]
  settings <- final_models$fits[[model]]$h
  parameter_rows[[model]] <- data.frame(model = model, parameters = paste(names(settings), unlist(settings),
    sep = "=", collapse = "; "), threshold = final_cutoffs[model], tuning_AUC = grid$AUC[best], tuning_Brier = grid$Brier[best])
  write.csv(grid, file.path(result_dir, paste0("调参记录_", model, ".csv")), row.names = FALSE)
}
write.csv(do.call(rbind, parameter_rows), file.path(result_dir, "模型参数与开发集固定阈值.csv"),
  row.names = FALSE)

# 2.6. Temporal and geographic assessment ----------------------------------
# No assessment outcomes enter fitting, hyperparameter tuning or cutoff choice.
# The 90% development target is not a guarantee of 90% sensitivity in new patients.
z <- list(data = development, outer = outer, history = history, p = outer_probability, cutoffs = outer_cutoff,
  summary = development_summary, final = list(bundle = final_models, cutoffs = final_cutoffs))
ranking <- z$summary
primary <- shap_model
shown <- models
sets <- list(nested = list(d = z$data, p = z$p, cut = z$cutoffs, fold = z$outer))
assessment_cohorts <- list(temporal = temporal, geographic = geographic)
for (cohort in names(assessment_cohorts)) {
  data <- assessment_cohorts[[cohort]]
  probability <- predict_library(z$final$bundle, data)
  cutoff <- matrix(z$final$cutoffs[colnames(probability)], nrow(data), ncol(probability), byrow = TRUE,
    dimnames = dimnames(probability))
  sets[[cohort]] <- list(d = data, p = probability, cut = cutoff, fold = rep(1L, nrow(data)))
}
pairs <- operating <- all_points <- list()
rules <- c("Sensitivity90", "Youden", "Fixed0.5")
for (cohort in names(sets)) {
  q <- sets[[cohort]]
  y <- q$d$Perforation
  draws <- bootstrap_rows(q$d, q$fold)
  ids <- ranking$model
  a <- lapply(setNames(ids, ids), function(m) bootstrap_auc(y, q$p[, m], draws))
  ci <- t(vapply(a, interval, numeric(2)))
  if (cohort == "nested") {
    point <- ranking[, c("model", "AUC", "Accuracy", "Sensitivity", "Specificity", "F1", "Brier")]
    point$Calibration_slope <- NA_real_
  } else {
    value <- do.call(rbind, lapply(ids, function(m) metric_values(y, q$p[, m], q$cut[, m])))
    point <- data.frame(model = ids, value, row.names = NULL)
  }
  point$cohort <- cohort
  point$Lower <- ci[point$model, 1]
  point$Upper <- ci[point$model, 2]
  point$Threshold_min <- apply(q$cut[, point$model, drop = FALSE], 2, min)
  point$Threshold_max <- apply(q$cut[, point$model, drop = FALSE], 2, max)
  all_points[[cohort]] <- point[, c("cohort", "model", "AUC", "Lower", "Upper", "Accuracy", "Sensitivity",
    "Specificity", "F1", "Brier", "Calibration_slope", "Threshold_min", "Threshold_max")]
  for (m in setdiff(models, primary)) {
    value <- point$AUC[match(primary, point$model)] - point$AUC[match(m, point$model)]
    diff_ci <- interval(a[[primary]] - a[[m]])
    brier <- point$Brier[match(primary, point$model)] - point$Brier[match(m, point$model)]
    brier_ci <- bootstrap_brier_delta(y, q$p[, primary], q$p[, m], draws)
    pairs[[paste(cohort, m)]] <- data.frame(cohort = cohort, model1 = primary, model2 = m, difference = value,
      lower = diff_ci[1], upper = diff_ci[2], Brier = brier, Brier_lower = brier_ci[1], Brier_upper = brier_ci[2])
  }
  for (rule in rules) {
    p <- q$p[, primary]
    if (cohort == "nested") {
      cut <- numeric(length(y))
      for (h in z$history) {
        training_y <- y[-h$rows]
        training_p <- h$train_p[, primary]
        threshold <- switch(rule, Sensitivity90 = target_cut(training_y, training_p, sensitivity_target),
          Youden = choose_cutoff(training_y, training_p), Fixed0.5 = 0.5)
        cut[h$rows] <- threshold
      }
    } else {
      threshold <- switch(rule, Sensitivity90 = target_cut(z$data$Perforation, z$p[, primary],
        sensitivity_target), Youden = choose_cutoff(z$data$Perforation, z$p[, primary]), Fixed0.5 = 0.5)
      cut <- rep(threshold, length(y))
    }
    fold_value <- fold_scores(q$d, matrix(p, ncol = 1, dimnames = list(NULL, primary)), matrix(cut,
      ncol = 1, dimnames = list(NULL, primary)), q$fold)
    operating[[paste(cohort, rule)]] <- data.frame(cohort = cohort, rule = rule, threshold_min = min(cut),
      threshold_max = max(cut), Sensitivity = mean(fold_value$Sensitivity), Specificity = mean(fold_value$Specificity),
      F1 = mean(fold_value$F1), detected = sum(y == 1 & p >= cut), perforations = sum(y == 1))
  }
}
evaluation <- list(sets = sets, ranking = ranking, primary = primary, shown = shown, points = do.call(rbind,
  all_points), metrics = do.call(rbind, all_points), paired = do.call(rbind, pairs), thresholds = do.call(rbind,
  operating))
sets <- evaluation$sets
write.csv(evaluation$points, file.path(result_dir, "表2_三队列八模型性能.csv"), row.names = FALSE)
write.csv(evaluation$paired, file.path(result_dir, "配对模型差异与患者bootstrap区间.csv"),
  row.names = FALSE)
write.csv(evaluation$thresholds, file.path(result_dir, "XGBoost不同开发集阈值的表现.csv"),
  row.names = FALSE)
# Bootstrap intervals condition on these predictions; model development is not
# repeated in the 2,000 performance bootstrap samples.
for (cohort in names(sets)) {
  rows <- data.frame(row_id = sets[[cohort]]$d$row_id, patient_id = sets[[cohort]]$d$patient_id, outcome = sets[[cohort]]$d$Perforation,
    sets[[cohort]]$p, check.names = FALSE)
  write.csv(rows, file.path(result_dir, paste0("逐侧预测_", cohort, ".csv")), row.names = FALSE)
}

# 2.7. Feature-selection stability ---------------------------------------------
# Each of 200 patient-block samples repeats penalty tuning and selection.
# All copies of a sampled patient remain together. Inclusion frequency is
# descriptive; it does not replace the full-development elastic-net feature list.
selections <- vector("list", 200L)
for (b in 1:200) {
  if (b%%25 == 0)
    message("Selection bootstrap ", b, " of 200")
  set.seed(seed + 900000L + b)
  sampled_data <- development[cluster_draw(development), ]
  selections[[b]] <- select_features(sampled_data, seed + 910000L + b)$variables
}
selection_count <- integer(length(predictors))
for (i in seq_along(predictors)) {
  for (selected in selections) selection_count[i] <- selection_count[i] + (predictors[i] %in% selected)
}
stability <- data.frame(variable = predictors, selected = selection_count, repetitions = 200L, frequency = selection_count/200)
write.csv(stability, file.path(result_dir, "特征选择bootstrap稳定性.csv"), row.names = FALSE)

# 2.8. Modality ablation using ordinary LR -----------------------------------------
# Hold the algorithm and patient folds fixed; change only the available inputs.
# These fits retain eligible modality variables, rather than the final screened
# subset used by the eight main models. Thus their AUCs are a separate comparison.
ans <- list()
for (nm in names(MODALITIES)) {
  vv <- MODALITIES[[nm]]
  p <- cut <- numeric(nrow(development))
  for (k in 1:5) {
    tr <- outer != k
    train_data <- development[tr, , drop = FALSE]
    f <- patient_folds(train_data, 5L, seed_for(train_data, 819L))
    oo <- numeric(nrow(train_data))
    for (j in 1:5) {
      b <- fit_modality(train_data[f != j, , drop = FALSE], vv)
      oo[f == j] <- predict_library(b, train_data[f == j, , drop = FALSE])[, 1]
    }
    ct <- choose_cutoff(train_data$Perforation, oo)
    b <- fit_modality(train_data, vv)
    p[!tr] <- predict_library(b, development[!tr, , drop = FALSE])[, 1]
    cut[!tr] <- ct
  }
  f <- patient_folds(development, 5L, seed_for(development, 819L))
  oo <- numeric(nrow(development))
  for (j in 1:5) {
    b <- fit_modality(development[f != j, , drop = FALSE], vv)
    oo[f == j] <- predict_library(b, development[f == j, , drop = FALSE])[, 1]
  }
  ans[[nm]] <- list(p = p, cutoffs = cut, final = fit_modality(development, vv), threshold = choose_cutoff(development$Perforation,
    oo))
  message("Completed modality ablation: ", nm)
}
modality <- ans
ablation_rows <- list()
for (cohort in names(sets)) {
  data <- sets[[cohort]]$d
  for (input in names(modality)) {
    fit <- modality[[input]]
    probability <- if (cohort == "nested")
      fit$p else predict_library(fit$final, data)[, 1]
    cutoff <- if (cohort == "nested")
      fit$cutoffs else fit$threshold
    ci <- cluster_ci(data, probability, cutoff, calibration = FALSE)
    auc <- ci[ci$metric == "AUC", ]
    brier <- ci[ci$metric == "Brier", ]
    ablation_rows[[paste(cohort, input)]] <- data.frame(cohort = cohort, input = input, AUC = auc$estimate,
      auc_lower = auc$lower, auc_upper = auc$upper, Brier = brier$estimate, brier_lower = brier$lower,
      brier_upper = brier$upper)
  }
}
ablation <- do.call(rbind, ablation_rows)
write.csv(ablation, file.path(result_dir, "五模态消融_AUC与Brier.csv"), row.names = FALSE)

# 2.9–2.10. Location comparisons and descriptive statistics ---------
cohorts <- list(development = development, temporal = temporal, geographic = geographic)
counts <- do.call(rbind, lapply(cohorts, function(x) data.frame(patients = length(unique(x$patient_id)),
  joints = nrow(x), events = sum(x$Perforation))))
counts$cohort <- rownames(counts)
write.csv(counts, file.path(result_dir, "队列人数与关节数.csv"), row.names = FALSE)
for (cohort in names(cohorts)) {
  write.csv(baseline_table(cohorts[[cohort]]), file.path(result_dir, paste0("分层特征_", cohort,
    ".csv")), row.names = FALSE)
}
for (center in unique(geographic$Center)) {
  write.csv(baseline_table(geographic[geographic$Center == center, ]), file.path(result_dir, paste0("中心分层特征_",
    center, ".csv")), row.names = FALSE)
}
# Location was linked to development patient clusters in the checked input file.
# It is never a predictor or outcome in the main perforation-presence model.
subgroup <- read.csv(file.path(input_dir, "穿孔位置亚组.csv"), fileEncoding = "UTF-8")
stopifnot(nrow(subgroup) == 102L, length(unique(subgroup$patient_id)) == 82L, all(subgroup$Perforation_location %in%
  0:1))
subgroup$Perforation <- subgroup$Perforation_location
location_table <- baseline_table(subgroup)
names(location_table)[names(location_table) == "No_perforation"] <- "Central"
names(location_table)[names(location_table) == "Perforation"] <- "Distal"
location_p <- vapply(predictors, function(v) cluster_wald(subgroup, v), numeric(1))
location_tests <- data.frame(variable = predictors, P = location_p, BH = p.adjust(location_p, "BH"))
write.csv(location_table, file.path(result_dir, "穿孔位置亚组特征.csv"), row.names = FALSE)
write.csv(location_tests, file.path(result_dir, "穿孔位置亚组检验.csv"), row.names = FALSE)

# SHAP interpretation in the three cohorts -----------------------
# XGBoost is the illustration model, not a declaration of a universally best fit.
# All cohorts share the development reference distribution and probability scale.
fitted <- final_models$fits[[shap_model]]
variables <- final_models$variables
reference <- development[, variables, drop = FALSE]
baseline <- mean(predict_shap(fitted, reference))
estimates <- list()
for (cohort in names(sets)) {
  set.seed(seed)
  x <- sets[[cohort]]$d[, variables, drop = FALSE]
  phi <- fastshap::explain(fitted, X = reference, newdata = x, pred_wrapper = predict_shap, nsim = 64,
    adjust = TRUE, baseline = baseline)
  stopifnot(all(is.finite(phi)), max(abs(rowSums(phi) + baseline - predict_shap(fitted, x))) < 1e-07)
  importance <- colMeans(abs(phi))
  estimates[[cohort]] <- list(x = x, phi = phi, summary = data.frame(cohort = cohort, variable = names(importance),
    mean_abs_shap = unname(importance)))
  write.csv(data.frame(row_id = sets[[cohort]]$d$row_id, phi), file.path(result_dir, paste0("SHAP逐侧贡献_",
    cohort, ".csv")), row.names = FALSE)
}
importance <- do.call(rbind, lapply(estimates, "[[", "summary"))
dev_importance <- estimates$nested$summary
ordered <- dev_importance$variable[order(dev_importance$mean_abs_shap)]
importance$feature <- factor(EN_NAMES[importance$variable], levels = EN_NAMES[ordered])
shap_plots <- dependence <- list()
for (cohort in names(sets)) {
  x <- estimates[[cohort]]$x
  phi <- estimates[[cohort]]$phi
  long <- list()
  dependence[[cohort]] <- list()
  for (variable in variables) {
    limits <- range(reference[[variable]])
    value <- if (diff(limits) > 0)
      pmin(1, pmax(0, (x[[variable]] - limits[1])/diff(limits))) else 0
    long[[variable]] <- data.frame(feature = factor(EN_NAMES[variable], levels = EN_NAMES[ordered]),
      SHAP = phi[, variable], value = value)
  }
  for (variable in rev(ordered)) {
    dependence[[cohort]][[variable]] <- data.frame(value = x[[variable]], SHAP = phi[, variable],
      outcome = sets[[cohort]]$d$Perforation)
  }
  shap_plots[[cohort]] <- list(long = do.call(rbind, long), summary = importance[importance$cohort ==
    cohort, ], title = cohort)
}
examples <- list()
dev_probability <- predict_one(fitted, development)
for (outcome in 0:1) {
  indices <- which(development$Perforation == outcome)
  index <- indices[which.min(abs(dev_probability[indices] - median(dev_probability[indices])))]
  examples[[outcome + 1]] <- list(outcome = outcome, prediction = dev_probability[index], phi = estimates$nested$phi[index,
    ])
}
write.csv(importance, file.path(result_dir, "三队列SHAP重要性.csv"), row.names = FALSE)

# Save fitted models and prepare figure data -------------------------
curves <- list()
for (cohort in names(sets)) {
  y <- sets[[cohort]]$d$Perforation
  probability <- sets[[cohort]]$p
  roc <- do.call(rbind, lapply(models, function(m) roc_data(y, probability[, m], m)))
  auc <- vapply(models, function(m) auc_fixed(y, probability[, m]), numeric(1))
  thresholds <- seq(0.05, 0.8, 0.01)
  decision <- list()
  for (model in models) {
    p <- probability[, model]
    net_benefit <- vapply(thresholds, function(t) mean((p >= t) * y) - mean((p >= t) * (1 - y)) *
      t/(1 - t), numeric(1))
    decision[[model]] <- data.frame(model = model, threshold = thresholds, NB = net_benefit)
  }
  dca <- rbind(do.call(rbind, decision), data.frame(model = "All", threshold = thresholds, NB = mean(y) -
    (1 - mean(y)) * thresholds/(1 - thresholds)), data.frame(model = "None", threshold = thresholds,
    NB = 0))
  curves[[cohort]] <- list(roc = roc, auc = auc, dca = dca, prevalence = mean(y))
}
plot_data <- list(summary = counts, points = evaluation$points, n_features = length(variables), feature_names = EN_NAMES,
  levels = LEVEL_NAMES, shap = list(plots = shap_plots, importance = importance, baseline = baseline,
    bound = max(vapply(shap_plots, function(p) max(abs(p$long$SHAP)), numeric(1))), upper = max(importance$mean_abs_shap) *
      1.06), dependence = dependence, examples = examples, curves = curves, stability = stability,
  ablation = ablation)
saveRDS(plot_data, file.path(result_dir, "绘图数据_本次计算.rds"))
saveRDS(list(models = final_models, thresholds = final_cutoffs, outer_folds = outer, outer_probabilities = outer_probability,
  outer_cutoffs = outer_cutoff, outer_history = history, fold_results = fold_results, assessment = evaluation$points,
  stability = stability, modality = modality), file.path(result_dir, "拟合模型与核验结果.rds"))




# Tables in manuscript order ---------------------------------------------------
write.csv(baseline_table(development), file.path(result_dir, "表1_开发队列特征.csv"), row.names = FALSE)
write.csv(baseline_table(temporal), file.path(result_dir, "表S4_时间队列特征.csv"), row.names = FALSE)
write.csv(baseline_table(geographic), file.path(result_dir, "表S5_外院队列特征.csv"), row.names = FALSE)
center_order <- unique(geographic$Center)
for (i in seq_along(center_order)) {
  data <- geographic[geographic$Center == center_order[i], ]
  write.csv(baseline_table(data), file.path(result_dir, paste0("中心特征_", center_order[i], ".csv")),
    row.names = FALSE)
}
all_tuning <- list()
for (model in models) {
  scores <- final_models$scores[[model]]
  parameter_columns <- setdiff(names(scores), c("AUC", "Brier"))
  settings <- character(nrow(scores))
  for (i in seq_len(nrow(scores))) {
    settings[i] <- paste(parameter_columns, as.character(scores[i, parameter_columns]), sep = "=",
      collapse = "; ")
  }
  all_tuning[[model]] <- data.frame(Model = model, Parameters = settings, Tuning_AUC = scores$AUC,
    CV_Brier = scores$Brier)
}
write.csv(do.call(rbind, all_tuning), file.path(result_dir, "表S2_候选参数与调参成绩.csv"),
  row.names = FALSE)
settings <- do.call(rbind, parameter_rows)
write.csv(settings, file.path(result_dir, "表S3_最终参数及阈值.csv"), row.names = FALSE)
selection_table <- stability
selection_table$Candidate <- selection_table$variable != "Sclerosis"
selection_table$Final_inclusion <- selection_table$variable %in% final_models$variables
write.csv(selection_table, file.path(result_dir, "表S7_入选特征与稳定性.csv"), row.names = FALSE)
write.csv(outer_features, file.path(result_dir, "表S8_五个外层训练集入选特征.csv"), row.names = FALSE)
wide_auc <- reshape(fold_results[, c("model", "fold", "AUC")], idvar = "model", timevar = "fold", direction = "wide")
mean_scores <- development_summary
names(mean_scores)[names(mean_scores) == "AUC"] <- "Mean_AUC"
wide_auc <- merge(wide_auc, mean_scores, by = "model", sort = FALSE)
write.csv(wide_auc, file.path(result_dir, "表S9_逐折AUC及平均表现.csv"), row.names = FALSE)
write.csv(evaluation$paired, file.path(result_dir, "表S10_配对模型比较.csv"), row.names = FALSE)
write.csv(evaluation$thresholds, file.path(result_dir, "表S11_XGBoost操作阈值.csv"), row.names = FALSE)
write.csv(ablation, file.path(result_dir, "表S12_模态消融.csv"), row.names = FALSE)
location_table$BH_P <- ""
first_rows <- match(location_tests$variable, location_table$variable)
location_table$BH_P[first_rows] <- fmt_p(location_tests$BH)
write.csv(location_table, file.path(result_dir, "表S13_穿孔位置亚组.csv"), row.names = FALSE)
# Tables S1 and S14 are supplied summaries from the existing supplementary file.
# They are read here for completeness; paired observer ratings are not available
# in the analysis CSVs, so this script does not recalculate agreement coefficients.
file.copy(file.path(input_dir, "原稿观察者一致性.csv"), file.path(result_dir, "表S1_原稿观察者一致性.csv"),
  overwrite = TRUE)
file.copy(file.path(input_dir, "原稿变量定义.csv"), file.path(result_dir, "表S14_变量定义.csv"),
  overwrite = TRUE)

# Figures 1–4 and S1–S7 --------------------------------------------------------
# Common graphical settings. Statistical panels are drawn with ggplot2/grid in R.
shap_blue <- "#009FE3"
shap_red <- "#FF303B"
cohort_titles <- c(nested = "Development", temporal = "Temporal assessment", geographic = "Geographic assessment")
cohort_colors <- c(Development = "#007BB9", `Temporal assessment` = "#78C9ED", `Geographic assessment` = "#EA926F")
PAL <- c(LR = "#0072B2", RF = "#009E73", DT = "#777777", SVM = "#CC79A7", KNN = "#56B4E9", MLP = "#D55E00",
  XGB = "#E69F00", LGBM = "#447744")
feature_names <- plot_data$feature_names
category_labels <- plot_data$levels
outcome_colors <- c(`No perforation` = shap_blue, Perforation = shap_red)

# Use Cairo for offscreen grid capture so Arial metrics match the exported PDF.
capture_device <- function(width, height) {
  grDevices::cairo_pdf(tempfile(fileext = ".pdf"), width = width, height = height)
}

paper_theme <- function(size = 10) {
  ggplot2::theme_bw(base_size = size, base_family = "Arial") + ggplot2::theme(panel.grid.minor = ggplot2::element_blank(),
    legend.position = "bottom", legend.title = ggplot2::element_blank(), plot.title = ggplot2::element_text(face = "bold",
      color = "black"), plot.tag = ggplot2::element_text(face = "bold", color = "black"), strip.background = ggplot2::element_rect(fill = "#E8F0F8"),
    plot.margin = ggplot2::margin(6, 8, 6, 8))
}

save_plot <- function(plot, name, width, height) {
  folder <- file.path(project_dir, "图", name)
  dir.create(folder, recursive = TRUE, showWarnings = FALSE)
  ggplot2::ggsave(file.path(folder, paste0(name, ".pdf")), plot, width = width, height = height, device = grDevices::cairo_pdf,
    bg = "white", limitsize = FALSE)
  ggplot2::ggsave(file.path(folder, paste0(name, ".png")), plot, width = width, height = height, device = ragg::agg_png,
    dpi = 300, bg = "white", limitsize = FALSE)
}

# Display the same feature scales in every cohort. Colours indicate observed outcome.
# Dependence plots describe model contributions, not formal interaction tests.
dependence_panel <- function(cohort, feature) {
  rows <- plot_data$dependence[[cohort]][[feature]]
  rows$Outcome <- factor(rows$outcome, levels = 0:1, labels = names(outcome_colors))
  all_cohorts <- lapply(plot_data$dependence, function(x) x[[feature]])
  y_limits <- range(vapply(all_cohorts, function(x) range(x$SHAP), numeric(2)), 0)
  if (feature %in% c("Age", "Duration_since_first_symptom")) {
    x_limits <- range(vapply(all_cohorts, function(x) range(x$value), numeric(2)))
    panel <- ggplot2::ggplot(rows, ggplot2::aes(value, SHAP)) + ggplot2::geom_hline(yintercept = 0,
      linetype = 2, color = "grey65") + ggplot2::geom_point(ggplot2::aes(color = Outcome), alpha = 0.78,
      size = 1.35) + ggplot2::coord_cartesian(xlim = x_limits, ylim = y_limits) + ggplot2::labs(x = if (feature ==
      "Age")
      "Age (years)" else "Symptom duration (months)")
  } else {
    rows$category <- factor(rows$value, levels = seq_along(category_labels[[feature]]) - 1, labels = category_labels[[feature]])
    stopifnot(!anyNA(rows$category))
    panel <- ggplot2::ggplot(rows, ggplot2::aes(category, SHAP)) + ggplot2::geom_hline(yintercept = 0,
      linetype = 2, color = "grey65") + ggplot2::geom_boxplot(width = 0.3, outlier.shape = NA,
      fill = "white", color = "grey45") + ggplot2::geom_point(ggplot2::aes(color = Outcome), alpha = 0.78,
      size = 1.35, position = ggplot2::position_jitter(width = 0.12, height = 0, seed = seed)) +
      ggplot2::scale_x_discrete(drop = FALSE, guide = ggplot2::guide_axis(n.dodge = 2)) + ggplot2::coord_cartesian(ylim = y_limits) +
      ggplot2::labs(x = unname(feature_names[feature]))
  }
  panel + ggplot2::scale_color_manual(values = outcome_colors, drop = FALSE, guide = ggplot2::guide_legend(override.aes = list(alpha = 1,
    size = 2))) + ggplot2::labs(y = "SHAP value (probability)") + paper_theme(9)
}

dependence_block <- function(cohort, features, columns) {
  panels <- lapply(features, function(feature) dependence_panel(cohort, feature))
  patchwork::wrap_plots(panels, ncol = columns, guides = "collect") & ggplot2::theme(legend.position = "bottom")
}

# A block title uses the same lower-case panel label throughout the figure set.
titled_block <- function(plot, title) {
  patchwork::wrap_elements(full = plot) + patchwork::plot_annotation(title = title, theme = ggplot2::theme(plot.title = ggplot2::element_text(family = "Arial",
    size = 12, face = "bold", color = "black"), plot.margin = ggplot2::margin(6, 6, 6, 6)))
}

waterfall_panel <- function(example) {
  baseline <- plot_data$shap$baseline
  phi <- example$phi[order(abs(example$phi), decreasing = TRUE)]
  finish <- baseline + cumsum(phi)
  start <- c(baseline, head(finish, -1))
  stopifnot(abs(tail(finish, 1) - example$prediction) < 1e-07)
  rows <- data.frame(position = rev(seq_along(phi)), start = start, end = finish, direction = ifelse(phi >=
    0, "Increases probability", "Decreases probability"))
  rows$left <- pmin(rows$start, rows$end)
  rows$right <- pmax(rows$start, rows$end)
  heading <- c("Example without perforation", "Example with perforation")[example$outcome + 1]
  ggplot2::ggplot(rows) + ggplot2::geom_vline(xintercept = baseline, color = "grey55", linetype = 2) +
    ggplot2::geom_vline(xintercept = example$prediction, color = "grey25", linetype = 3) + ggplot2::geom_segment(data = head(rows,
    -1), ggplot2::aes(x = end, xend = end, y = position - 0.34, yend = position - 0.66), color = "grey55") +
    ggplot2::geom_rect(ggplot2::aes(xmin = left, xmax = right, ymin = position - 0.34, ymax = position +
      0.34, fill = direction)) + ggplot2::scale_y_continuous(breaks = rev(seq_along(phi)), labels = unname(feature_names[names(phi)])) +
    ggplot2::scale_fill_manual(values = c(`Increases probability` = shap_red, `Decreases probability` = shap_blue)) +
    ggplot2::labs(title = heading, subtitle = sprintf("Baseline %.3f | Predicted probability %.3f",
      baseline, example$prediction), x = "Predicted probability", y = NULL) + paper_theme(8.5) +
    ggplot2::theme(legend.position = "none")
}

summary_panels <- function(cohort, show_mean = FALSE) {
  p <- plot_data$shap$plots[[match(cohort, names(cohort_titles))]]
  bar <- ggplot2::ggplot(p$summary, ggplot2::aes(feature, mean_abs_shap)) + ggplot2::geom_col(fill = unname(cohort_colors[cohort_titles[cohort]])) +
    ggplot2::coord_flip(ylim = c(0, plot_data$shap$upper)) + ggplot2::labs(x = NULL, y = "Mean absolute SHAP value") +
    paper_theme(10)
  dots <- ggplot2::ggplot(p$long, ggplot2::aes(SHAP, feature, color = value))
  if (show_mean) {
    dots <- dots + ggplot2::geom_col(data = p$summary, ggplot2::aes(mean_abs_shap, feature), inherit.aes = FALSE,
      orientation = "y", width = 0.65, fill = "#CFD9E0")
  }
  dots <- dots + ggplot2::geom_vline(xintercept = 0, color = "grey70") + ggplot2::geom_point(position = ggplot2::position_jitter(height = 0.17,
    width = 0, seed = seed), alpha = 0.78, size = 1) + ggplot2::scale_color_gradient(low = shap_blue,
    high = shap_red, name = "Coded value", limits = c(0, 1)) + ggplot2::coord_cartesian(xlim = c(-plot_data$shap$bound,
    plot_data$shap$bound)) + ggplot2::labs(x = "SHAP value (probability)", y = NULL) + paper_theme(10)
  list(bar = bar, dots = dots)
}

# Figure: manuscript or supplementary panel
# Cohort screening and patient-grouped modelling are shown in one figure.
library(grid)
cohort_map <- png::readPNG(file.path(project_dir, "素材", "图1地图示意.png"))
cohorts <- data.frame(patients = c(216, 104, 76), joints = c(306, 141, 83), perforations = c(102, 24,
  29))
temporal_exclusions <- c(imaging = 13L, psychological = 1L, condylar = 3L)
temporal_screened <- cohorts$patients[2] + sum(temporal_exclusions)
total <- colSums(cohorts)

flow_blue <- c("#3F6392", "#3E9FC9")
flow_warm <- c("#F1D98F", "#DB827F")
flow_mix <- c("#7DAAC2", "#F1D98F", "#DB827F")
flow_font <- "Arial"

flow_box <- function(x, y, w, h, label, colors = flow_blue, size = 13) {
  for (k in 3:1) grid.roundrect(x + 0.0015 * k, y - 0.0025 * k, w, h, r = unit(3, "mm"), gp = gpar(fill = adjustcolor("#444444",
    alpha.f = 0.065), col = NA))
  grid.roundrect(x, y, w, h, r = unit(3, "mm"), gp = gpar(fill = linearGradient(colors, x1 = 0, x2 = 1,
    y1 = 0, y2 = 0), col = "white", lwd = 2.5))
  grid.text(label, x, y, gp = gpar(fontfamily = flow_font, fontsize = size, fontface = "bold", col = "white",
    lineheight = 1.16))
}

flow_details <- function(x, y, label, size = 12) {
  grid.roundrect(x, y, 0.293, 0.163, r = unit(2.5, "mm"), gp = gpar(fill = "#EEF3F9", col = NA))
  grid.text(label, x - 0.129, y, just = "left", gp = gpar(fontfamily = flow_font, fontsize = size,
    col = "#202A35", lineheight = 1.25))
}

flow_arrow <- function(x, y, xx, yy, col = "#4386AF") {
  grid.lines(c(x, xx), c(y, yy), gp = gpar(col = col, fill = col, lwd = 1.8), arrow = arrow(type = "closed",
    length = unit(2, "mm")))
}

draw_cohort <- function() {
  grid.newpage()
  grid.rect(gp = gpar(fill = "white", col = NA))
  flow_box(0.5, 0.915, 0.92, 0.108, sprintf("Multicentre study of surgical ADDwoR patients in China\nIncluded: %d patients · %d joints · %d perforated joints",
    total[1], total[2], total[3]), flow_mix, 16)
  x <- c(0.172, 0.5, 0.828)
  heading <- c("Development · ZJUSS\nOct 2023 – June 2025\nScreened: 243 patients", sprintf("Temporal assessment · ZJUSS\nJune 2025 – June 2026\nScreened: %d patients",
    temporal_screened), "Geographic assessment\nJune 2025 – June 2026\nThree other hospitals")
  middle <- c("Excluded: 27 patients\n• Imaging interval >3 months: 15\n• Inadequate image quality: 8\n• Prior surgery/concurrent disease: 4",
    sprintf("Excluded: %d patients\n• CT–MRI interval >3 months: %d\n• Psychological condition: %d\n• Other condylar disease: %d",
      sum(temporal_exclusions), temporal_exclusions[1], temporal_exclusions[2], temporal_exclusions[3]),
    "HBMU2: 32 patients · 32 joints\nLZU2: 33 patients · 40 joints\nHEH: 11 patients · 11 joints")
  included <- c("Development cohort", "Temporal cohort", "Geographic cohort")
  for (i in seq_along(x)) {
    color <- if (i == 3)
      flow_warm else flow_blue
    arrow_color <- if (i == 3)
      "#B27A5D" else "#4386AF"
    flow_arrow(0.5 + (x[i] - 0.5) * 0.54, 0.856, x[i], 0.821, arrow_color)
    flow_box(x[i], 0.742, 0.298, 0.15, heading[i], color)
    flow_arrow(x[i], 0.662, x[i], 0.623, arrow_color)
    flow_details(x[i], 0.538, middle[i])
    flow_arrow(x[i], 0.452, x[i], 0.413, arrow_color)
    label <- sprintf("%s\n%d patients · %d joints\n%d perforated joints", included[i], cohorts$patients[i],
      cohorts$joints[i], cohorts$perforations[i])
    flow_box(x[i], 0.338, 0.298, 0.142, label, color)
  }
  grid.raster(cohort_map, x = 0.828, y = 0.14, height = 0.237, interpolate = TRUE)
  grid.text(paste("ZJUSS, Stomatology Hospital, Zhejiang University", "HBMU2, Second Hospital of Hebei Medical University",
    "LZU2, Second Hospital and Clinical Medical School, Lanzhou University", "HEH, Hebei Eye Hospital",
    sep = "\n"), 0.043, 0.147, just = "left", gp = gpar(fontfamily = flow_font, fontsize = 10, col = "#38434C",
    lineheight = 1.23))
}

pipe_ink <- "#172D40"
pipe_arrow <- "#536B7D"
pipe_blue <- "#CBE8F8"
pipe_inner <- "#F8E2D6"
pipe_gold <- "#F9EDCF"
pipe_green <- "#DCEBE2"

pipe_box <- function(x, y, w, h, label, fill = "white", bold = FALSE, size = 10) {
  grid.roundrect(x/100, y/100, w/100, h/100, r = unit(1.5, "mm"), gp = gpar(fill = fill, col = "#A0B1BD",
    lwd = 0.75))
  grid.text(label, x/100, y/100, gp = gpar(fontfamily = "sans", fontsize = size, fontface = if (bold)
    "bold" else "plain", col = pipe_ink, lineheight = 1.13))
}
pipe_title <- function(x, y, label, size = 12) grid.text(label, x/100, y/100, gp = gpar(fontfamily = "sans",
  fontsize = size, fontface = "bold", col = pipe_ink))
pipe_link <- function(x, y, end = TRUE, lty = 1) {
  a <- if (end)
    arrow(type = "closed", length = unit(1.8, "mm")) else NULL
  grid.lines(x/100, y/100, arrow = a, gp = gpar(col = pipe_arrow, lwd = 1.25, fill = pipe_arrow, lty = lty))
}
pipe_panel <- function(x, w, fill) grid.roundrect(x/100, 0.535, w/100, 0.79, r = unit(3, "mm"), gp = gpar(fill = fill,
  col = NA))

draw_pipeline <- function() {
  grid.newpage()
  grid.rect(gp = gpar(fill = "white", col = NA))
  pipe_panel(26, 48, "#EDF6FB")
  pipe_panel(74.5, 47, "#F8FAFB")
  pipe_title(50, 96, "Model development and cohort evaluation", 17)
  pipe_title(26, 90, "Patient-grouped nested cross-validation")
  pipe_title(74.5, 90, "Full-development refit and assessment")

  pipe_box(26, 82.3, 38, 8, "Development cohort\n216 patients · 306 joints", pipe_blue, TRUE)
  pipe_box(13.7, 69.8, 21, 8, "Outer training set\n4 of 5 patient folds", pipe_blue, TRUE)
  pipe_box(38.3, 69.8, 20, 8, "Outer holdout\n1 of 5 patient folds", pipe_blue, TRUE)
  pipe_link(c(26, 26, 13.7, 13.7), c(78.3, 76, 76, 73.8))
  pipe_link(c(26, 26, 38.3, 38.3), c(78.3, 76, 76, 73.8))
  pipe_box(14.5, 56.5, 22.6, 13, "Inner 5-fold CV\nTraining-fold preprocessing\nand elastic-net screening\nAlgorithm parameter grid",
    pipe_inner)
  pipe_link(c(13.7, 13.7), c(65.8, 63))
  pipe_box(14.5, 44.8, 22.6, 7.5, "Rank settings by pooled AUC\nLower Brier score breaks ties", pipe_inner)
  pipe_link(c(14.5, 14.5), c(50, 48.55))
  pipe_box(14.5, 33.4, 22.6, 10.2, "Refit each model on the\nouter training set; derive its\nthreshold from cross-fitted\ntraining predictions",
    pipe_blue)
  pipe_link(c(14.5, 14.5), c(41.05, 38.5))
  pipe_box(38.1, 33.4, 20, 10.2, "Apply fitted models and\nthresholds to the\nouter holdout", pipe_blue)
  pipe_link(c(25.8, 28.1), c(33.4, 33.4))
  pipe_link(c(38.3, 38.3), c(65.8, 38.5))
  pipe_box(26, 20.8, 38, 7.8, "Summarize performance across five outer folds\nMean fold AUC, Brier score and threshold metrics",
    pipe_blue, TRUE)
  pipe_link(c(38.1, 38.1, 26, 26), c(28.3, 26.3, 26.3, 24.7))

  pipe_box(65, 82.3, 23.4, 8, "Same development cohort\n216 patients · 306 joints", pipe_blue, TRUE)
  pipe_box(65, 69.5, 23.4, 9, "Repeat screening and\n5-fold tuning within development", pipe_gold)
  pipe_link(c(65, 65), c(78.3, 74))
  pipe_box(65, 57.8, 23.4, 8, "Refit all eight models on\nall development patients", pipe_gold, TRUE)
  pipe_link(c(65, 65), c(65, 61.8))
  pipe_box(65, 45.6, 23.4, 11, "Freeze each fit and threshold\nUse pooled outer holdout predictions\nto target 90% sensitivity",
    pipe_gold)
  pipe_link(c(65, 65), c(53.8, 51.1))
  pipe_box(88.2, 69.5, 17.4, 9, "Same-center temporal\n104 patients · 141 joints", pipe_green, size = 9.3)
  pipe_box(88.2, 57.8, 17.4, 8, "Geographic assessment\n76 patients · 83 joints", pipe_green, size = 9.3)
  pipe_box(76.1, 29.2, 37.7, 10, "Evaluate all eight frozen models\nin each assessment cohort", pipe_green,
    TRUE, 11)
  pipe_link(c(65, 65), c(40.1, 34.2))
  pipe_link(c(96.9, 98.7, 98.7, 95), c(69.5, 69.5, 29.2, 29.2))
  pipe_link(c(88.2, 88.2), c(53.8, 34.2))
  pipe_link(c(48.1, 51.2, 51.2, 53.3), c(33.4, 33.4, 45.6, 45.6), lty = 2)
  grid.text("Patient grouping kept both joints together. Dashed line shows transfer of development holdout probabilities for threshold setting.",
    0.5, 0.08, gp = gpar(fontsize = 9, col = pipe_ink, fontfamily = "sans"))
}

cohort_grob <- grid::grid.grabExpr(draw_cohort(), width = 13, height = 7.7, device = capture_device)
save_plot(cohort_grob, "图1_患者筛查与队列构成", 13, 7.7)
pipeline_grob <- grid::grid.grabExpr(draw_pipeline(), width = 13, height = 8.5, device = capture_device)
save_plot(pipeline_grob, "图2_模型开发与评估流程", 13, 8.5)

# Figure: manuscript or supplementary panel
# Development uses mean outer-fold AUC and its patient-bootstrap interval.
# This avoids presenting pooled ROC AUC as if it were the mean fold AUC.
points <- plot_data$points[plot_data$points$cohort == "nested", ]
points <- points[order(points$AUC), ]
points$label <- factor(paste0(points$model, "  ", sprintf("%.3f", points$AUC)), levels = paste0(points$model,
  "  ", sprintf("%.3f", points$AUC)))
development_auc <- ggplot2::ggplot(points, ggplot2::aes(AUC, label, color = model)) + ggplot2::geom_errorbar(ggplot2::aes(xmin = Lower,
  xmax = Upper), orientation = "y", width = 0.15) + ggplot2::geom_point(size = 2) + ggplot2::scale_color_manual(values = PAL) +
  ggplot2::scale_x_continuous(limits = c(0.5, 1), breaks = seq(0.5, 1, 0.1)) + ggplot2::labs(title = "Development: mean outer-fold AUC",
  x = "AUC (95% CI)", y = NULL) + paper_theme(9) + ggplot2::theme(legend.position = "none")
rocs <- dcas <- list()
for (cohort in names(cohort_titles)) {
  curves <- plot_data$curves[[cohort]]
  model_order <- names(curves$auc)
  labels <- paste0(model_order, " ", sprintf("%.3f", curves$auc))
  rocs[[cohort]] <- ggplot2::ggplot(curves$roc, ggplot2::aes(FPR, TPR, color = model)) + ggplot2::geom_step(linewidth = 0.65) +
    ggplot2::geom_abline(linetype = 2, color = "grey70") + ggplot2::scale_color_manual(values = PAL,
    breaks = model_order, labels = labels) + ggplot2::labs(title = cohort_titles[cohort], x = "1 - Specificity",
    y = "Sensitivity") + paper_theme(9)
  dcas[[cohort]] <- ggplot2::ggplot(curves$dca, ggplot2::aes(threshold, NB, color = model)) + ggplot2::geom_line(linewidth = 0.65) +
    ggplot2::scale_color_manual(values = c(PAL, All = "grey60", None = "black")) + ggplot2::coord_cartesian(ylim = c(-0.05,
    max(0.15, curves$prevalence))) + ggplot2::labs(title = cohort_titles[cohort], x = "Threshold probability",
    y = "Net benefit") + paper_theme(9)
}
figure <- patchwork::wrap_plots(development_auc, rocs$temporal, rocs$geographic, dcas$nested, dcas$temporal,
  dcas$geographic, ncol = 3) + patchwork::plot_annotation(tag_levels = "a", tag_prefix = "(", tag_suffix = ")",
  theme = ggplot2::theme(plot.tag = ggplot2::element_text(family = "Arial", face = "bold", size = 12)))
save_plot(figure, "图3_模型表现", 15, 9)


# Figure: manuscript or supplementary panel
# All 11 final features are displayed; examples use development patients only.
panels <- summary_panels("nested", show_mean = TRUE)
summary <- titled_block(panels$dots, "(a) SHAP summary")
examples <- patchwork::wrap_plots(lapply(rev(plot_data$examples), waterfall_panel), ncol = 1)
individual <- titled_block(examples, "(b) Individual explanations")
continuous_features <- c("Age", "Duration_since_first_symptom")
categorical_features <- setdiff(names(plot_data$dependence$nested), continuous_features)
continuous <- titled_block(dependence_block("nested", continuous_features, 2), "(c) Continuous predictors")
categorical <- titled_block(dependence_block("nested", categorical_features, 3), "(d) Categorical predictors")
left <- patchwork::wrap_plots(patchwork::wrap_elements(full = summary), patchwork::wrap_elements(full = individual),
  ncol = 1, heights = c(1, 1.35))
right <- patchwork::wrap_plots(patchwork::wrap_elements(full = continuous), patchwork::wrap_elements(full = categorical),
  ncol = 1, heights = c(1, 3))
figure <- patchwork::wrap_plots(left, right, ncol = 2, widths = c(1, 2.1))
save_plot(figure, "图4_开发队列SHAP", 20, 14)


# Figure: manuscript or supplementary panel
# Clinical image panels are the original submitted material, not simulated images.
folder <- file.path(project_dir, "图", "补图S1_MRI与CBCT示例")
dir.create(folder, recursive = TRUE, showWarnings = FALSE)
stopifnot(file.copy(file.path(project_dir, "素材", "原稿影像拼图.pdf"), file.path(folder, "补图S1_MRI与CBCT示例.pdf"),
  overwrite = TRUE))


# Figure: manuscript or supplementary panel
# Inclusion frequencies are descriptive, not a second selection threshold.
stability <- plot_data$stability
stability <- stability[stability$variable != "Sclerosis", ]
stability$label <- factor(feature_names[stability$variable], levels = feature_names[stability$variable[order(stability$frequency)]])
figure <- ggplot2::ggplot(stability, ggplot2::aes(label, frequency)) + ggplot2::geom_col(fill = cohort_colors["Development"]) +
  ggplot2::coord_flip() + ggplot2::scale_y_continuous(limits = c(0, 1), labels = function(x) paste0(100 *
  x, "%")) + ggplot2::labs(x = NULL, y = "Selection frequency in 200 patient-block bootstrap samples") +
  paper_theme(11)
save_plot(figure, "补图S2_筛选稳定性", 8, 6)


# Figure: manuscript or supplementary panel
# The background, feature order and SHAP scale are shared across cohorts.
panels <- summary_panels("temporal")
figure <- patchwork::wrap_plots(panels$bar, panels$dots, ncol = 2) + patchwork::plot_annotation(tag_levels = "a",
  tag_prefix = "(", tag_suffix = ")", theme = ggplot2::theme(plot.tag = ggplot2::element_text(family = "Arial",
    face = "bold", size = 12)))
save_plot(figure, "补图S3_时间队列SHAP", 12, 6)


# Figure: manuscript or supplementary panel
# The background, feature order and SHAP scale are shared across cohorts.
panels <- summary_panels("geographic")
figure <- patchwork::wrap_plots(panels$bar, panels$dots, ncol = 2) + patchwork::plot_annotation(tag_levels = "a",
  tag_prefix = "(", tag_suffix = ")", theme = ggplot2::theme(plot.tag = ggplot2::element_text(family = "Arial",
    face = "bold", size = 12)))
save_plot(figure, "补图S4_外院队列SHAP", 12, 6)


# Figure: manuscript or supplementary panel
# Mean absolute contributions use the same fixed XGBoost fit and development background.
importance <- plot_data$shap$importance
importance$Cohort <- factor(cohort_titles[importance$cohort], levels = unname(cohort_titles))
figure <- ggplot2::ggplot(importance, ggplot2::aes(feature, mean_abs_shap, fill = Cohort)) + ggplot2::geom_col(position = "dodge",
  width = 0.74) + ggplot2::coord_flip() + ggplot2::scale_fill_manual(values = cohort_colors) + ggplot2::labs(x = NULL,
  y = "Mean absolute SHAP value (probability)") + paper_theme(11)
save_plot(figure, "补图S5_队列SHAP比较", 9, 6)


# Figure: manuscript or supplementary panel
# Development dependence plots already appear in Figure 4 and are not repeated.
# Panels a and b show continuous features; c and d show the other nine features.
continuous_features <- c("Age", "Duration_since_first_symptom")
categorical_features <- setdiff(names(plot_data$dependence$nested), continuous_features)
stopifnot(length(categorical_features) == 9L)
a <- titled_block(dependence_block("temporal", continuous_features, 2), "(a) Temporal assessment - continuous predictors")
b <- titled_block(dependence_block("geographic", continuous_features, 2), "(b) Geographic assessment - continuous predictors")
c <- titled_block(dependence_block("temporal", categorical_features, 3), "(c) Temporal assessment - categorical predictors")
d <- titled_block(dependence_block("geographic", categorical_features, 3), "(d) Geographic assessment - categorical predictors")
figure <- patchwork::wrap_plots(patchwork::wrap_elements(full = a), patchwork::wrap_elements(full = b),
  patchwork::wrap_elements(full = c), patchwork::wrap_elements(full = d), ncol = 2, heights = c(1,
    3))
save_plot(figure, "补图S6_两评估队列依赖图", 20, 16)


# Figure: manuscript or supplementary panel
# All five input sets use ordinary LR and the same development patient folds.
# Development uses pooled outer-fold predictions for this modality comparison.
labels <- c(clinical = "Clinical only", MRI = "MRI only", CBCT = "CBCT only", imaging = "MRI + CBCT",
  full = "All inputs")
points <- plot_data$ablation
points$Inputs <- factor(labels[points$input], levels = rev(unname(labels)))
points$Cohort <- factor(cohort_titles[points$cohort], levels = unname(cohort_titles))
stopifnot(all(table(points$Cohort, points$Inputs) == 1L))
dodge <- ggplot2::position_dodge(width = 0.62, orientation = "y")
figure <- ggplot2::ggplot(points, ggplot2::aes(AUC, Inputs, color = Cohort, shape = Cohort)) + ggplot2::geom_errorbar(ggplot2::aes(xmin = auc_lower,
  xmax = auc_upper), orientation = "y", width = 0.13, position = dodge, linewidth = 0.7) + ggplot2::geom_point(position = dodge,
  size = 2.7) + ggplot2::scale_color_manual(values = cohort_colors) + ggplot2::scale_shape_manual(values = c(16,
  17, 15)) + ggplot2::scale_x_continuous(breaks = seq(0.5, 1, 0.1)) + ggplot2::labs(x = "AUC (95% CI)",
  y = NULL) + paper_theme(11) + ggplot2::theme(panel.grid.major.y = ggplot2::element_blank())
save_plot(figure, "补图S7_模态消融", 9, 5.5)



