## ------------------------------------------------------------------------------------
## face / brokenstick without random holdouts
## ------------------------------------------------------------------------------------
test.genericfit_FACE_BS <- function() {
  # library("growthcurveSL")
  options(growthcurveSL.verbose = TRUE)
  data(cpp)
  cpp <- cpp[!is.na(cpp[, "haz"]), ]
  # covars <- c("apgar1", "apgar5", "parity", "gagebrth", "mage", "meducyrs", "sexn")
  # define holdout col:
  cpp_holdout <- add_holdout_ind(data = cpp, ID = "subjid", hold_column = "hold", random = TRUE, seed = 12345)
  holdout_col <- cpp_holdout[["hold"]]

  run_algo <- function(fit.package, fit.algorithm) {
    mfit_useY <- fit_model(ID = "subjid", t_name = "agedays", x = "agedays", y = "haz",
                           train_data = cpp[!holdout_col, ], valid_data = cpp[holdout_col, ],
                           params  = list(fit.package = fit.package, fit.algorithm = fit.algorithm, predict.w.Y = TRUE, name = "useY"))
    fit_preds_1 <- predict_model(mfit_useY, newdata = cpp)
    # print(head(fit_preds_1[]))
    print(mfit_useY$getMSE)

    mfit_cor <- fit_model(ID = "subjid", t_name = "agedays", x = "agedays", y = "haz",
                          train_data = cpp[!holdout_col, ], valid_data = cpp[holdout_col, ],
                          params  = list(fit.package = fit.package, fit.algorithm = fit.algorithm, predict.w.Y = FALSE, name = "correct"))
    print(mfit_cor$getMSE)
    fit_preds_2 <- predict_model(mfit_cor, newdata = cpp)
    # print(head(fit_preds_2[]))
  }

  run_algo("face", "face")
  run_algo("brokenstick", "brokenstick")
  run_algo("speedglm", "glm")
  res_GLM2 <- run_algo("glm", "glm")
  res_GLM3 <- run_algo("h2o", "glm")
  res_GBM <- run_algo("h2o", "gbm")
  res_DRF <- run_algo("h2o", "randomForest")
  res_DP <- run_algo("h2o", "deeplearning")

}

