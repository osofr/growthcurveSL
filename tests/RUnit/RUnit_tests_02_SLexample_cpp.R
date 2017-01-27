## ------------------------------------------------------------------------------------
## face / brokenstick based on random holdouts
## ------------------------------------------------------------------------------------
test.holdoutfit_FACE_BS_h2o <- function() {
  # require("gridisl")
  # require("data.table")
  require("h2o")
  h2o::h2o.init(nthreads = -1)
  # library("growthcurveSL")
  options(growthcurveSL.verbose = TRUE)
  data(cpp)
  cpp <- cpp[!is.na(cpp[, "haz"]), ]
  # covars <- c("apgar1", "apgar5", "parity", "gagebrth", "mage", "meducyrs", "sexn")
  # define holdout col:
  cpp_holdout <- gridisl::add_holdout_ind(data = cpp, ID = "subjid", hold_column = "hold", random = TRUE, seed = 12345)
  holdout_col <- cpp_holdout[["hold"]]

  ID <- "subjid"
  t_name <- "agedays"
  x <- "agedays"
  y <- "haz"

  run_algo <- function(fit.package, fit.algorithm) {
    # Fit, training on non-holdouts and using holdouts as validation set (for scoring only)
    library("magrittr")
    library("gridisl")
    mfit1 <- defModel(name = "useY",
                      estimator = paste0(fit.package, "__", fit.algorithm),
                      predict.w.Y = TRUE,
                      family = "gaussian") %>%
             fit_growth(data = cpp_holdout,
                        method = "holdout",
                        ID = "subjid",
                        t_name = "agedays",
                        x = "agedays",
                        y = "haz",
                        hold_column = "hold")

    print("Holdout MSE, using the holdout Y for prediction"); print(mfit1$getMSE)
    # FACE MSE: [1] 0.2271609
    # BS MSE: [1] 0.02650036
    # predict for previously used holdout / validation set:
    preds_holdout_1 <- growthcurveSL:::predict_growth(mfit1, holdout = TRUE)
    print(nrow(preds_holdout_1))     # [1] 453
    print(preds_holdout_1[])

    # $m.1.face.useY
    #            preds
    #            <num>
    #   1:  1.13497268
    #   2: -0.94470508
    #   3:  0.47621843
    #   4:  0.66113476
    #  ---
    # 449: -0.38502193
    # 450:  0.89649979
    # 451: -0.07370885
    # 452:  1.01063739
    # 453: -1.15840422

    # m.1.brokenstick.useY
    #           preds
    #           <num>
    #   1:  1.2869520
    #   2: -1.3341364
    #   3:  0.6166304
    #   4:  0.8225809
    #  ---
    # 449: -0.6343973
    # 450:  1.2434923
    # 451:  0.1596915
    # 452:  1.1282024
    # 453: -1.9495589

    ## Obtain predictions for a model trained on all data:
    preds_alldat_train <- predict_growth(mfit1, newdata = cpp_holdout, add_subject_data = TRUE)
    print(preds_alldat_train[])

    mfit2 <-  defModel(name = "correct",
                       estimator = paste0(fit.package, "__", fit.algorithm),
                       predict.w.Y = FALSE,
                       family = "gaussian") %>%
              fit_growth(data = cpp_holdout,
                         method = "holdout",
                         ID = "subjid",
                         t_name = "agedays",
                         x = "agedays",
                         y = "haz",
                         hold_column = "hold")

    print("Holdout MSE, hiding the holdout Y for prediction"); print(mfit2$getMSE)
    # FACE MSE: [1] 1.211989
    # BS MSE: [1] 1.186241
    # speed or reg GLM MSE: [1] 1.813257
    # xgb GLM MSE: [1] 1.619442
    # h2o GLM MSE: [1] 1.619442
    # xgb GBM MSE [1] 1.531563
    # h2o GBM MSE [1] 1.531809
    # xgb DRF MSE [1] 1.682466
    # h2o DRF MSE [1] 1.531855
    # h2o deeplearning MSE [1] 1.643752

    # predict for previously used holdout / validation set:
    preds_holdout_2 <- predict_growth(mfit2, holdout = TRUE)
    print(nrow(preds_holdout_2)) # [1] 453
    print(preds_holdout_2[])

    # $m.1.face.correct:
    #             preds
    #             <num>
    #   1:  0.813006593
    #   2: -0.508011348
    #   3:  0.360798108
    #   4:  0.689897034
    #  ---
    # 449: -0.006218735
    # 450:  0.549329997
    # 451: -0.262319345
    # 452:  0.854839807
    # 453:  0.725667441

    # m.1.brokenstick.correct
    #            preds
    #            <num>
    #   1:  0.55848245
    #   2: -0.74885770
    #   3:  0.62770593
    #   4:  1.03315090
    #  ---
    # 449: -0.02824323
    # 450:  0.94357989
    # 451: -0.14039338
    # 452:  0.98565065
    # 453:  0.70983542

    ## Obtain predictions for a model trained on all data:
    preds_alldat_train <- predict_growth(mfit2, newdata = cpp_holdout, add_subject_data = TRUE)
    print(preds_alldat_train[])

    return(list(mfit1 =  mfit1, mfit2 =  mfit2))
  }

  res_FACE <- run_algo("face", "face")
  res_BS <- run_algo("brokenstick", "brokenstick")
  # res_GLM1 <- run_algo("speedglm", "glm")
  # res_GLM2 <- run_algo("glm", "glm")
  res_xbg_GLM <- run_algo("xgboost", "glm")
  res_GLM3 <- run_algo("h2o", "glm")
  res_xbg_GBM <- run_algo("xgboost", "gbm")
  res_GBM <- run_algo("h2o", "gbm")
  res_xgb_DRF <- run_algo("xgboost", "drf")
  res_DRF <- run_algo("h2o", "randomForest")
  res_DP <- run_algo("h2o", "deeplearning")

  mfits_stack <- gridisl::make_PredictionStack(
                                              res_FACE$mfit1, res_FACE$mfit2,
                                               res_BS$mfit1, res_BS$mfit2,
                                               res_xbg_GLM$mfit2, res_GLM3$mfit2,
                                               res_xbg_GBM$mfit2, res_GBM$mfit2,
                                               res_xgb_DRF$mfit2,  res_DRF$mfit2,
                                               res_DP$mfit2
                                               )

  print(mfits_stack$get_best_MSEs(K = 2))
  print(mfits_stack$get_best_MSE_table(K = 2))

  res_GBM$mfit1$get_modelfits_grid()
  (grids <- mfits_stack$get_modelfits_grid())

  gridisl::make_model_report(mfits_stack, data = cpp_holdout, K = 2,
                  file.name = paste0("BS_ALL_", getOption("growthcurveSL.file.name")),
                  format = "html",
                  # openFile = FALSE)
                  openFile = TRUE)

  tab <- mfits_stack$get_best_MSE_table(K = 10)
  # get the model objects for top K models:
  top_model <- mfits_stack$get_best_models(K = 1)

  train_dat <- get_train_data(res_GBM$mfit2)
  val_dat <- get_validation_data(res_GBM$mfit2)

  preds_tgrid_FACE <- predict_save_tgrid(res_FACE$mfit2, cpp_holdout, ID, t_name, y, tmin = 1, tmax = 500, incr = 50)
  preds_tgrid_FACE[]
  preds_tgrid_BS <- predict_save_tgrid(res_BS$mfit2, cpp_holdout, ID, t_name, y, tmin = 1, tmax = 500, incr = 5)
  preds_tgrid_BS[]

}

test.holdoutSL.GLM <- function() {
  # library("growthcurveSL")
  require("h2o")
  h2o::h2o.init(nthreads = -1)
  options(growthcurveSL.verbose = TRUE)
  data(cpp)
  cpp <- cpp[!is.na(cpp[, "haz"]), ]
  covars <- c("apgar1", "apgar5", "parity", "gagebrth", "mage", "meducyrs", "sexn")

  ## add holdout indicator column
  cpp_holdout <- gridisl::add_holdout_ind(data = cpp, ID = "subjid", hold_column = "hold", random = TRUE, seed = 12345)

  # ----------------------------------------------------------------
  # Perform fitting with regularlized GLMs using h2o.grid
  # ----------------------------------------------------------------
  GRIDparams <- gridisl::defModel(name = "GLM",
                                  estimator = "h2o__glm", family = "gaussian",
                                  search_criteria = list(strategy = "RandomDiscrete", max_models = 5),
                                  param_grid = list(
                                    alpha = c(0,1,seq(0.1,0.9,0.01)),
                                    lambda = c(seq(0.001, 0.009, by = 0.001), seq(0.01, 0.09, by = 0.01))), # , seq(0.1, 0.9, by = 0.02)),
                                  missing_values_handling = c("MeanImputation"),
                                  seed = 23)

  mfit_holdGLM <- fit_growth(GRIDparams, ID = "subjid", t_name = "agedays",
                              y = "haz",
                              # x = c("agedays", covars),
                              x = c("agedays", covars, "nY", "meanY", "medianY", "minY", "maxY"),
                              data = cpp_holdout, method = "holdout",
                              hold_column = "hold") # , use_new_features = TRUE

    print("Holdout MSE, using the holdout Y for prediction"); print(str(mfit_holdGLM$getMSE))
    # List of 1
    #  $ :List of 5
    #   ..$ M.1.h2o.glm.grid.1.GLM: num 1.59
    #   ..$ M.1.h2o.glm.grid.2.GLM: num 1.6
    #   ..$ M.1.h2o.glm.grid.3.GLM: num 1.62
    #   ..$ M.1.h2o.glm.grid.4.GLM: num 1.62
    #   ..$ M.1.h2o.glm.grid.5.GLM: num 1.63
    # NULL

    # holdPredDT <- growthcurveSL:::predict_holdout(mfit_holdGLM, predict_only_bestK_models = 5, add_subject_data = TRUE)
    holdPredDT <- predict_growth(mfit_holdGLM, holdout = TRUE, add_subject_data = TRUE)
    print("GLM holdPredDT"); print(holdPredDT[])

    preds_best_train <- predict_growth(mfit_holdGLM, add_subject_data = TRUE)
    print("GLM preds_best_train"); print(preds_best_train[])

    preds_best_all <- predict_growth(mfit_holdGLM, newdata = cpp_holdout, add_subject_data = TRUE)
    print("GLM preds_best_all"); print(preds_best_all[])

    print("TOP MSE: "); print(min(unlist(mfit_holdGLM$getMSE)), na.rm = TRUE)
    print(mfit_holdGLM$get_best_MSEs(2))
    print(mfit_holdGLM$get_best_MSEs(15))
    BEST_GLM_model <- mfit_holdGLM$get_best_models(K = 1)[[1]]


    grids <- mfit_holdGLM$get_modelfits_grid()
    tab <- mfit_holdGLM$get_best_MSE_table(K = 10)

    gridisl::make_model_report(mfit_holdGLM, data = cpp_holdout,
                    K = 10,
                    file.name = paste0("GLMs_", getOption("growthcurveSL.file.name")),
                    title = paste0("Growth Curve Imputation with GLM"),
                    format = "html",
                    # keep_md = TRUE,
                    openFile = FALSE)
                    # openFile = TRUE)

}

## ------------------------------------------------------------------------------------
## Holdout Growth Curve SL with model scoring based on random holdouts
## ------------------------------------------------------------------------------------
test.holdoutSL.GLM.GBM <- function() {
  # library("growthcurveSL")
  require("h2o")
  h2o::h2o.init(nthreads = -1)
  options(growthcurveSL.verbose = TRUE)
  data(cpp)
  cpp <- cpp[!is.na(cpp[, "haz"]), ]
  covars <- c("apgar1", "apgar5", "parity", "gagebrth", "mage", "meducyrs", "sexn")

  ## add holdout indicator column
  cpp_holdout <- gridisl::add_holdout_ind(data = cpp, ID = "subjid", hold_column = "hold", random = TRUE, seed = 12345)

  ## ----------------------------------------------------------------
  ## Define learners (glm, grid glm and grid gbm)
  gbm_hyper_models <- list(ntrees = c(100, 200, 300, 500),
                           learn_rate = c(0.005, 0.01, 0.03, 0.06),
                           max_depth = c(3, 4, 5, 6, 9),
                           sample_rate = c(0.7, 0.8, 0.9, 1.0),
                           col_sample_rate = c(0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8),
                           balance_classes = c(TRUE, FALSE))

  GRIDparams <- gridisl::defModel(estimator = "h2o__glm", family = "gaussian",
                                  alpha = 0.3,
                                  nlambdas = 50,
                                  lambda_search = TRUE,
                                  seed = 23) +
                gridisl::defModel(name = "GLM",
                                  estimator = "h2o__glm", family = "gaussian",
                                  search_criteria = list(strategy = "RandomDiscrete", max_models = 2),
                                  param_grid = list(
                                    alpha = c(0,1,seq(0.1,0.9,0.1)),
                                    lambda = c(0,1e-7,1e-5,1e-3,1e-1)),
                                  seed = 23,
                                  missing_values_handling = c("MeanImputation")) +
                gridisl::defModel(estimator = "h2o__gbm", family = "gaussian",
                                  search_criteria = list(strategy = "RandomDiscrete", max_models = 1, max_runtime_secs = 60*60),
                                  param_grid = gbm_hyper_models,
                                  seed = 23,
                                  stopping_rounds = 5,
                                  stopping_tolerance = 1e-4,
                                  stopping_metric = "MSE",
                                  score_tree_interval = 10)

  # --------------------------------------------------------------------------------------------
  ## Fit the model based on additional special features (summaries) of the outcomes:
  mfit_hold <- fit_growth(GRIDparams,
                           ID = "subjid",
                           t_name = "agedays",
                           x = c("agedays", covars),
                           y = "haz",
                           data = cpp_holdout,
                           hold_column = "hold",
                           method = "holdout",
                           use_new_features = TRUE)

  ## ------------------------------------------------------------------------------------------------
  ## Combining fits for observed data points, holdout points and grid of points.
  ## Returns a single dataset, one row per subject.
  all_preds <- predict_all(mfit_hold, cpp_holdout, add_grid = TRUE, add_holdout = TRUE)

  ## ------------------------------------------------------------------------------------------------
  ## add a column with subject specific growth trajectories (plots) & plot the whole thing
  library("magrittr")
  all_preds <- all_preds %>%
               add_fit_plots()

  all_preds %>%
    trelliscopejs::trelliscope(name = "test")

  ## ------------------------------------------------------------------------------------------------
  ## Add cognostics and convert into an object that can be accepted by hbgd package
  ## devtools::install_github('hafen/hbgd', ref = "tidy")
  library("magrittr")
  sex_var <- "sex"
  ID <- "subjid"
  add_sex <- cpp_holdout %>%
             dplyr::distinct_(ID, sex_var) %>%
             dplyr::rename_("sex" = sex_var) %>%
             dplyr::left_join(all_preds) %>%
             tibble::as_tibble()

  add_sex <- add_sex %>%
             plyr::mutate(fit = purrr::map2(fit, sex, ~ add_cogs_persubj(fit_dat = .x, sex = .y)))

  add_sex <- add_sex %>%
   hbgd::add_trajectory_plot() %>%
    dplyr::select_("subjid", "panel")

  add_sex %>%
    trelliscopejs::trelliscope(name = "test")


  ## ------------------------------------------------------------------------------------------------
  ## Test MSE evaluation, training / validation data retrieval
  ## ------------------------------------------------------------------------------------------------
  train_dat <- gridisl::get_train_data(mfit_hold)
  print(train_dat)
  checkTrue(nrow(train_dat)==sum(!cpp_holdout[["hold"]]))
  val_dat <- gridisl::get_validation_data(mfit_hold)
  print(val_dat)
  checkTrue(nrow(val_dat)==sum(cpp_holdout[["hold"]]))

  print("Holdout MSE, using the holdout Y for prediction"); print(str(mfit_hold$getMSE))
  # List of 3
  #  $ :List of 1
  #   ..$ M.1.h2o.glm: num 1.46
  #  $ :List of 3
  #   ..$ M.2.h2o.glm.grid.1.GLM: num 1.92
  #   ..$ M.2.h2o.glm.grid.2.GLM: num 1.92
  #   ..$ M.2.h2o.glm.grid.3.GLM: num 1.92
  #  $ :List of 2
  #   ..$ M.3.h2o.gbm.grid.1: num 1.48
  #   ..$ M.3.h2o.gbm.grid.2: num 1.57
  # NULL
  ## OLD:
  # $grid.glm.1
  # [1] 1.797792
  # $grid.glm.2
  # [1] 1.976214
  # $grid.glm.3
  # [1] 1.976399
  # $grid.gbm.4
  # [1] 1.597701
  # $grid.gbm.5
  # [1] 1.700594
  # $h2o.glm.reg03
  # [1] 1.776226

  ## Predictions for new data based on best SL model trained on all data:
  preds_alldat <- predict_growth(mfit_hold, newdata = cpp_holdout, add_subject_data = TRUE)
  preds_alldat[]

  ## Predictions for all holdout data points for all models trained on non-holdout data:
  preds_holdout <- predict_growth(mfit_hold, holdout = TRUE, add_subject_data = TRUE)
  preds_holdout[]

  ## Predictions for training data from models trained on non-holdout data (default):
  ## NOT IMPLEMENTED
  # preds_train <- predict_model(mfit_hold, add_subject_data = TRUE)
  # preds_train[]
  # preds_train <- predict_model(mfit_hold, predict_only_bestK_models = 2, add_subject_data = TRUE)
  # preds_train[]
  # holdPredDT <- growthcurveSL:::predict_holdout(mfit_hold, predict_only_bestK_models = 5, add_subject_data = TRUE)

  ## does not work right now, since it doesn't know how to define the special features!!!!!
  # preds_holdout_alldat <- predict_model(mfit_hold, newdata = cpp_holdout, predict_only_bestK_models = 1, add_subject_data = TRUE)
  ## Instead, have to first manually define features for entire dataset - predictions will be for a model trained on non-holdouts only!
  # cpp_plus <- define_features(cpp_holdout, nodes = mfit_hold$OData_train$nodes, train_set = TRUE, holdout = FALSE)
  cpp_plus <- define_features_drop(cpp_holdout, ID = "subjid", t_name = "agedays", y = "haz", train_set = TRUE)
  preds_holdout_alldat <- gridisl::predict_SL(mfit_hold, newdata = cpp_plus, add_subject_data = TRUE)
  preds_holdout_alldat[]

  print("10 best MSEs among all learners: "); print(mfit_hold$get_best_MSEs(K = 5))
  models <- mfit_hold$get_best_models(K = 10)
  models[[1]]

  print("Top 5 models: "); print(models)
  res_tab <- mfit_hold$get_best_MSE_table(K = 5)
  print("5 best models among all learners: "); print(res_tab)
  gridisl::make_model_report(mfit_hold, data = cpp_holdout, K = 10, format = "html", openFile = FALSE)

  ## ------------------------------------------------------------------------------------------------
  ## Predicting the entire curve (grid)
  ## Check  predict_save_tgrid() works as intended
  ## ------------------------------------------------------------------------------------------------
  cpp_all_train <- define_features_drop(cpp_holdout, ID = "subjid", t_name = "agedays", y = "haz", train_set = TRUE)
  cpp_all_grid <- define_tgrid(cpp_all_train, ID = "subjid", t_name = "agedays", y = "haz", tmin = 1, tmax = 500, incr = 2)
  preds_grid <- predict_growth(mfit_hold, newdata = cpp_all_grid, grid = TRUE, add_subject_data = TRUE)
  preds_grid[]

  ## obtain and save grid preds in one go:
  preds_grid2 <- predict_save_tgrid(mfit_hold,
                                    cpp_holdout,
                                    ID = "subjid",
                                    t_name = "agedays",
                                    y = "haz",
                                    tmin = 1,
                                    tmax = 500,
                                    incr = 2)
  preds_grid2[]

 ## Test reporting
 gridisl::make_model_report(mfit_hold, data = cpp_holdout,
                  K = 10,
                  file.name = paste0("GLMs_", getOption("growthcurveSL.file.name")),
                  title = paste0("Growth Curve Imputation with GLM"),
                  format = "html",
                  # keep_md = TRUE,
                  # openFile = FALSE)
                  openFile = TRUE)


  ## Ask the fit function to determine random holdouts internally:
  mfit_hold3 <- fit_growth(GRIDparams, ID = "subjid", t_name = "agedays", x = c("agedays", covars), y = "haz",
                           data = cpp_holdout,
                           method = "holdout",
                           random_holdout = TRUE,
                           use_new_features = TRUE)
  ## NOT IMPLEMENTED
  # gridisl::save_best_h2o_model(mfit_hold3, file.path = "/Users/olegsofrygin/GoogleDrive/HBGDki/ImputationSL/sofware")


}

## ------------------------------------------------------------------------------------
## Growth Curve SL with model scoring based on full V-FOLD CROSS-VALIDATION
## ------------------------------------------------------------------------------------
test.CV.SL <- function() {
  # library("growthcurveSL")
  require("data.table")
  require("h2o")
  h2o::h2o.init(nthreads = -1)
  options(growthcurveSL.verbose = TRUE)
  options(gridisl.verbose = TRUE)
  data(cpp)
  cpp <- cpp[!is.na(cpp[, "haz"]), ]
  covars <- c("apgar1", "apgar5", "parity", "gagebrth", "mage", "meducyrs", "sexn")

  ID <- "subjid"
  t_name <- "agedays"
  y <- "haz"

  ## define CV folds (respecting that multiple observations per subject must fall within the same fold)
  cpp_folds <- gridisl::add_CVfolds_ind(cpp, ID = "subjid", nfolds = 5, seed = 23)

  ## ------------------------------------------------------------------------------------------------
  ## Define learners (glm, grid glm and grid gbm)
  ## glm grid learner:
  glm_hyper <- list(alpha = c(0,1,seq(0.1,0.9,0.1)),
                    lambda = c(0,1e-7,1e-5,1e-3,1e-1))

  ## gbm grid learner:
  gbm_hyper <- list(
                   ntrees = c(100, 200, 300, 500),
                   learn_rate = c(0.005, 0.01, 0.03, 0.06),
                   max_depth = c(3, 4, 5, 6, 9),
                   sample_rate = c(0.7, 0.8, 0.9, 1.0),
                   col_sample_rate = c(0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8),
                   balance_classes = c(TRUE, FALSE))

  ## this fails on latest build (lambda_search)
  # h2o.glm.reg03 <- function(..., alpha = 0.3, nlambdas = 50, lambda_search = TRUE) gridisl::h2o.glm.wrapper(..., alpha = alpha, nlambdas = nlambdas, lambda_search = lambda_search)
  GRIDparams <-
                # gridisl::defModel(estimator = "h2o__glm", family = "gaussian",
                #            alpha = 0.3,
                #            nlambdas = 50,
                #            lambda_search = TRUE,
                #            seed = 23) +
                gridisl::defModel(estimator = "h2o__glm", family = "gaussian",
                         search_criteria = list(strategy = "RandomDiscrete", max_models = 3),
                         seed = 23,
                         param_grid = glm_hyper,
                         missing_values_handling = c("MeanImputation")) +
                gridisl::defModel(estimator = "h2o__gbm", family = "gaussian",
                         search_criteria = list(strategy = "RandomDiscrete", max_models = 2, max_runtime_secs = 60*60),
                         seed = 23,
                         param_grid = gbm_hyper,
                         stopping_rounds = 5, stopping_tolerance = 1e-4, stopping_metric = "MSE", score_tree_interval = 10)

  ## ------------------------------------------------------------------------------------------------
  ## CV with manually defined fold column, no curve summary features are used as predictors
  ## -> Internal CV metrics must match all manual model scoring results
  mfit_cv <- fit_growth(GRIDparams,
                         data = cpp_folds,
                         method = "cv",
                         ID = "subjid",
                         t_name = "agedays",
                         x = c("agedays", covars),
                         y = "haz",
                         fold_column = "fold",
                         use_new_features = TRUE)

  ## ------------------------------------------------------------------------------------------------
  ## Combining fits for observed data points, holdout points and grid of points.
  ## Returns a single dataset, one row per subject.
  all_preds <- predict_all(mfit_cv, cpp_folds, add_grid = TRUE, add_holdout = TRUE)

  ## ------------------------------------------------------------------------------------------------
  ## add a column with subject specific growth trajectories (plots) & plot the whole thing
  library("magrittr")
  all_preds <- all_preds %>%
               add_fit_plots()

  all_preds %>%
    trelliscopejs::trelliscope(name = "test")

  ## ------------------------------------------------------------------------------------------------
  ## Add cognostics and convert into an object that can be accepted by hbgd package
  ## devtools::install_github('hafen/hbgd', ref = "tidy")
  library("magrittr")
  sex_var <- "sex"
  ID <- "subjid"
  add_sex <- cpp_folds %>%
             dplyr::distinct_(ID, sex_var) %>%
             dplyr::rename_("sex" = sex_var) %>%
             dplyr::left_join(all_preds) %>%
             tibble::as_tibble()

  add_sex <- add_sex %>%
             plyr::mutate(fit = purrr::map2(fit, sex, ~ add_cogs_persubj(fit_dat = .x, sex = .y)))

  add_sex <- add_sex %>%
   hbgd::add_trajectory_plot() %>%
    dplyr::select_("subjid", "panel")

  add_sex %>%
    trelliscopejs::trelliscope(name = "test")


  ## ------------------------------------------------------------------------------------------------
  ## Best (re-trained) model predictions on data used for CV training (default):
  preds_alldat1 <- predict_growth(mfit_cv, add_subject_data = FALSE)
  print(preds_alldat1[])

  ## Best model predictions for new data, must match:
  preds_alldat2 <- predict_growth(mfit_cv, newdata = cpp_folds, add_subject_data = FALSE)
  print(preds_alldat2[])
  checkTrue(all.equal(preds_alldat1, preds_alldat2))

  ## Predictions for best CV model (not re-trained, trained only on non-holdouts), must match:
  ## NOT IMPLEMENTED
  # preds_best_CV <- predict_model(mfit_cv, add_subject_data = TRUE)
  # preds_best_CV[]

  ## return training data used (as data.table)
  train_dat <- gridisl::get_train_data(mfit_cv)
  train_dat[]
  ## return validation data used (as data.table)
  valid_data <- gridisl::get_validation_data(mfit_cv)
  valid_data[]

  ## Obtain out of sample CV predictions for all training data-points (the model used last to obtain these predictions)
  cv_preds <- gridisl::get_out_of_sample_predictions(mfit_cv)
  cv_preds[]
  names(cv_preds) <- "preds"
  ## Validate out-of-sample predictions based on validation data used in CV
  cv_valid_preds_2 <- predict_growth(mfit_cv, holdout = TRUE)
  cv_valid_preds_2[]
  checkTrue(all.equal(cv_preds, cv_valid_preds_2))

  ## ------------------------------------------------------------------------------------------------
  ## Predicting the entire curve (grid)
  ## ------------------------------------------------------------------------------------------------
  cpp_all_train <- define_features_drop(cpp_folds, ID = "subjid", t_name = "agedays", y = "haz", train_set = TRUE)
  cpp_all_grid <- define_tgrid(cpp_all_train, ID = "subjid", t_name = "agedays", y = "haz", tmin = 1, tmax = 500, incr = 2)
  cpp_all_grid[, ("fold") := NULL]
  preds_grid <- predict_growth(mfit_cv, newdata = cpp_all_grid, grid = TRUE, add_subject_data = TRUE)
  preds_grid[]

  ## Re-score the models. By default prev. validation data is used if it was provided:
  (MSE_1 <- unlist(gridisl::eval_MSE(mfit_cv))) ## Use internally h2o-evaluated CV MSE
  # M.1.h2o.glm.grid.1 M.1.h2o.glm.grid.2 M.1.h2o.glm.grid.3 M.2.h2o.gbm.grid.1 M.2.h2o.gbm.grid.2
  #         1.817785           1.837714           1.838126           1.637032           1.516350

  ## MSE evaluated on the training data:
  (MSE_train <- unlist(gridisl::eval_MSE(mfit_cv, gridisl::get_train_data(mfit_cv)))) # re-score using training data
  # M.1.h2o.glm.grid.1 M.1.h2o.glm.grid.2 M.1.h2o.glm.grid.3 M.2.h2o.gbm.grid.1 M.2.h2o.gbm.grid.2
  #          0.5887432          0.5888702          0.5888981          0.3293420          0.4053220

  ## Re-scoring based on validation data used in CV:
  (MSE_2 <- unlist(gridisl::eval_MSE(mfit_cv, gridisl::get_validation_data(mfit_cv)))) ## Rescore on validation data, but use old saved y values
  checkTrue(abs(sum(MSE_1 - MSE_2)) <= 10^-5)


  ## ------------------------------------------------------------------------------------
  ## CHECKING CV IMPLEMENTATION VS. INTERNAL H2O CV
  ## ------------------------------------------------------------------------------------
  train_dat <- get_train_data(mfit_cv)
  valid_data <- get_validation_data(mfit_cv)

  # Re-score on training data (should be equivalent to h2o-based CV metrics):
  preds_cv_check <- gridisl::eval_MSE(mfit_cv, get_train_data(mfit_cv))


  ## Internal CV MSE and CV evaluated on train_data should be the same:
  ## NOT IMPLEMENTED
  # mfit_cv$getmodel_ids
  # for (model_name in names(mfit_cv$getmodel_ids)) {
  #   check_model <- mfit_cv$getmodel_byname(model_name)
  #   # Internal H2O holdout CV predictions by fold:
  #   # cv_preds_fold <- h2o.cross_validation_predictions(check_model[[1]])
  #   # List of individual holdout predictions:
  #   # Holdout (out-of-sample) predictions for all of training data:
  #   cv_preds_all <- as.vector(h2o.cross_validation_holdout_predictions(check_model[[1]])) - preds_cv_check[[model_name]]

  #   # MUST BE TRUE FOR ALL MODELS
  #   print(paste0("CV validity for model: ", model_name))
  #   test1 <- mfit_cv$getMSE[[model_name]] - check_model[[1]]@model$cross_validation_metrics@metrics$MSE < (10^-6)
  #   print(test1); checkTrue(test1)
  #   test2 <- sum(cv_preds_all, na.rm = TRUE) < 10^-8
  #   print(test2); checkTrue(test1)
  # }

  ## Make report, save grid predictions and out of sample predictions
  # fname <- paste0(data.name, "_", "CV_gridSL_")
  gridisl::make_model_report(mfit_cv, K = 10, data = cpp_folds,
                  # file.name = paste0(fname, getOption("growthcurveSL.file.name")),
                  title = paste0("Growth Curve Imputation with cpp Data"),
                  format = "html", keep_md = FALSE,
                  # openFile = TRUE)
                  openFile = FALSE)

}


## ------------------------------------------------------------------------------------
## Holdout Growth Curve SL based on residuals from initial glm regression (model scoring based on random holdouts)
## ------------------------------------------------------------------------------------
test.residual.holdoutSL <- function() {
  # library("growthcurveSL")
  require("h2o")
  h2o::h2o.init(nthreads = -1)
  options(growthcurveSL.verbose = TRUE)
  data(cpp)
  cpp <- cpp[!is.na(cpp[, "haz"]), ]
  covars <- c("apgar1", "apgar5", "parity", "gagebrth", "mage", "meducyrs", "sexn")

  ## ----------------------------------------------------------------
  ## Define learners (glm, grid glm and grid gbm)
  ## ----------------------------------------------------------------
  ## glm grid learner:
  alpha_opt <- c(0,1,seq(0.1,0.9,0.1))
  lambda_opt <- c(0,1e-7,1e-5,1e-3,1e-1)
  glm_hyper_models <- list(
                            alpha = alpha_opt, lambda = lambda_opt,
                            missing_values_handling = c("MeanImputation"))

  ## gbm grid learner:
  gbm_hyper_models <- list(
                           ntrees = c(100, 200, 300, 500),
                           learn_rate = c(0.005, 0.01, 0.03, 0.06),
                           max_depth = c(3, 4, 5, 6, 9),
                           sample_rate = c(0.7, 0.8, 0.9, 1.0),
                           col_sample_rate = c(0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8),
                           balance_classes = c(TRUE, FALSE))

    GRIDparams <-
              gridisl::defModel(estimator = "h2o__glm", family = "gaussian",
                         alpha = 0.3,
                         nlambdas = 50,
                         lambda_search = TRUE,
                         seed = 23) +
                gridisl::defModel(estimator = "h2o__glm", family = "gaussian",
                         search_criteria = list(strategy = "RandomDiscrete", max_models = 3),
                         param_grid = glm_hyper_models, seed = 23) +
                gridisl::defModel(estimator = "h2o__gbm", family = "gaussian",
                         search_criteria = list(strategy = "RandomDiscrete", max_models = 2, max_runtime_secs = 60*60),
                         param_grid = gbm_hyper_models, seed = 23,
                         stopping_rounds = 5, stopping_tolerance = 1e-4, stopping_metric = "MSE", score_tree_interval = 10)


  ## add holdout indicator column
  cpp_holdout <- gridisl::add_holdout_ind(data = cpp, ID = "subjid", hold_column = "hold", random = TRUE, seed = 12345)

  ## fit the model based on additional special features (summaries) of the outcomes:
  # fit_growth

  mfit_resid_hold <- fit_growth(GRIDparams,
                                ID = "subjid",
                                t_name = "agedays",
                                x = c("agedays", covars),
                                y = "haz",
                                data = cpp_holdout,
                                method = 'holdout',
                                hold_column = "hold",
                                use_new_features = TRUE)

  print("Holdout MSE, using the residual holdout Y prediction"); print(mfit_resid_hold$getMSE)
  # [1] "Holdout MSE, using the holdout Y for prediction"
  # $grid.glm.1
  # [1] 1.977161
  # $grid.glm.2
  # [1] 1.977161
  # $grid.glm.3
  # [1] 1.977161
  # $grid.gbm.4
  # [1] 1.435757
  # $grid.gbm.5
  # [1] 1.495895
  # $h2o.glm.reg03
  # [1] 1.776226

  ## Predictions for all holdout data points for all models trained on non-holdout data only:
  preds_holdout_all <- growthcurveSL:::predict_holdout(mfit_resid_hold, add_subject_data = TRUE)
  preds_holdout_all[]
  ## Predictions for new data based on best SL model re-trained on all data:
  preds_alldat <- predict_growth(mfit_resid_hold, newdata = cpp_holdout, add_subject_data = TRUE)
  preds_alldat[]
}
