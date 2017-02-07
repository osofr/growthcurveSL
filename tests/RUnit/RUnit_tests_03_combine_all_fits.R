test.combine.all.model.fits <- function() {
  # library("growthcurveSL")
  library("magrittr")
  library("gridisl")
  require("h2o")
  h2o::h2o.init(nthreads = -1)
  options(growthcurveSL.verbose = TRUE)
  data(cpp)
  cpp <- cpp[!is.na(cpp[, "haz"]), ]
  covars <- c("apgar1", "apgar5", "parity", "gagebrth", "mage", "meducyrs", "sexn")

  ## add holdout indicator column:
  cpp_holdout <- gridisl::add_holdout_ind(data = cpp, ID = "subjid", hold_column = "hold", random = TRUE, seed = 12345)
  ## add CV fold column:
  cpp_folds <- gridisl::add_CVfolds_ind(cpp, ID = "subjid", nfolds = 5, seed = 23)

  ## ----------------------------------------------------------------
  ## Hyper parameters for h2o and xgboost GBM:
  ## ----------------------------------------------------------------
  h2o_GBM_hyper <- list(
    ntrees = c(50, 100, 150, 200),
    learn_rate = c(.05, .1),
    max_depth = c(seq(3, 19, 4), 25),
    sample_rate = c(.5, .75, 1),
    col_sample_rate = seq(0.1, 1, 0.05),
    col_sample_rate_per_tree = c(.4, .6, .8, 1),
    col_sample_rate_change_per_level = seq(0.9, 1.1, 0.01)
  )

  xgb_GBM_hyper = list(
    nrounds = c(50, 100, 150, 200),
    learning_rate = c(.05, .1),
    max_depth = c(seq(3, 19, 4), 25),
    subsample = c(.5, .75, 1),
    colsample_bytree = c(.4, .6, .8, 1),
    min_child_weight = c(1, 3, 5, 7),
    gamma = c(.0, .05, seq(.1, .9, by=.2), 1),
    lambda = c(.1, .5, 1, 2, 5), # lambda = c(1,2,5),
    alpha = c(0, .1, .5),
    max_delta_step = c(0, 2, 5, 10)
    )

  ## ----------------------------------------------------------------
  ## Fit brokenstick with random holdout MSE evaluation
  ## ----------------------------------------------------------------
  mfit_BS <- defModel(estimator = "brokenstick__brokenstick",,
                      predict.w.Y = FALSE) %>%
             fit_growth(data = cpp_holdout,
                      method = "holdout",
                      ID = "subjid",
                      t_name = "agedays",
                      x = "agedays",
                      y = "haz",
                      hold_column = "hold")

  all_preds_BS <- predict_all(mfit_BS, cpp_holdout)
  all_preds_BS_hbgd <- all_preds_BS %>%
                           convert_to_hbgd(cpp_holdout, "sex", "brokenstick")

  ## ----------------------------------------------------------------
  ## Define Learners with h2o and xgboost GBM:
  ## ----------------------------------------------------------------
  GRIDparams <- gridisl::defModel(estimator = "h2o__gbm", family = "gaussian",
                                  search_criteria = list(strategy = "RandomDiscrete", max_models = 2, max_runtime_secs = 60*60),
                                  param_grid = h2o_GBM_hyper,
                                  seed = 123456) +

                gridisl::defModel(estimator = "xgboost__gbm", family = "gaussian",
                                  search_criteria = list(strategy = "RandomDiscrete", max_models = 2),
                                  param_grid = xgb_GBM_hyper,
                                  seed = 123456)

  # --------------------------------------------------------------------------------------------
  ## Fit the model with random holdout validation. Use additional special features (summaries) of growth curves.
  mfit_holdSL <- fit_growth(GRIDparams,
                           ID = "subjid",
                           t_name = "agedays",
                           x = c("agedays", covars),
                           y = "haz",
                           data = cpp_holdout,
                           hold_column = "hold",
                           method = "holdout",
                           use_new_features = TRUE)

  ## Returns a single dataset, one row per subject. Combine fits for observed data points, holdout points and grid of points.
  all_preds_holdSL <- predict_all(mfit_holdSL, cpp_holdout)
  all_preds_holdSL_hbgd <- all_preds_holdSL %>%
                           convert_to_hbgd(cpp_holdout, "sex", "brokenstick")

  # --------------------------------------------------------------------------------------------
  ## Fit the model with V fold cross-validation Use additional special features (summaries) of growth curves.
  ## -> Uses internal CV metrics to manually re-score fitted CV models
  mfit_SLcv <- fit_growth(GRIDparams,
                         data = cpp_folds,
                         method = "cv",
                         ID = "subjid",
                         t_name = "agedays",
                         x = c("agedays", covars),
                         y = "haz",
                         fold_column = "fold",
                         use_new_features = TRUE)

  all_preds_cvSL <- predict_all(mfit_SLcv, cpp_holdout)
  all_preds_cvSL_hbgd <- all_preds_cvSL %>%
                         convert_to_hbgd(cpp_holdout, "sex", "brokenstick")


  library("dplyr")
  all_preds_combined <- (all_preds_BS_hbgd %>% rename(brokenstick = fit)) %>%
                         left_join(
                          all_preds_holdSL_hbgd %>% rename(holdoutSL = fit)
                          ) %>%
                         left_join(
                          all_preds_cvSL_hbgd %>% rename(cvSL = fit)
                          )
  save(all_preds_combined, file = "all_preds_combined.Rd")
                     trelliscopejs::trelliscope(name = "test")

}
