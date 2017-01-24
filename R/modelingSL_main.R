#' @import R6
#' @import data.table
#' @import ggiraph
#' @import GriDiSL
#' @importFrom graphics axis barplot hist par text  legend plot
#' @importFrom methods is
#' @importFrom stats approx binomial gaussian coef glm.control glm.fit plogis predict qlogis qnorm quantile rnorm terms var predict glm.control
#' @importFrom utils data head str
#' @importFrom stats as.formula glm na.exclude rbinom terms.formula
#' @importFrom assertthat assert_that
#' @importFrom ggplot2 ggplot geom_point geom_errorbar theme_bw coord_flip aes position_dodge alpha
# @importFrom ggiraph geom_point_interactive ggiraph
NULL

#' @rdname fit_growth.ModelStack
#' @export
fit_growth <- function(...) { UseMethod("fit_growth") }

# ---------------------------------------------------------------------------------------
#' Fit Discrete Growth Curve SuperLearner
#'
#' Define and fit discrete SuperLearner for growth curve modeling.
#' Model selection (scoring) can be based on either MSE for a single random (or last)
#' holdout data-point for each subject (method = "holdout") or
#' V-fold cross-validated MSE which uses entire subjects
#' (entire growth curves) for model validation (method = "cv").
#' @param models Parameters specifying the model(s) to fit.
#' This must be a result of calling \code{GriDiSL::defModel(...) + GriDiSL::defModel(...)} functions.
#' See \code{\link{defModel}} for additional information.
#' @param method The type of model selection procedure when fitting several models.
#' Possible options are "none" (no model selection),
#' "cv" (model selection with V-fold cross-validation), and
#' "holdout" (model selection based on validation holdout sample).
#' @param data Input dataset, can be a \code{data.frame} or a \code{data.table}.
#' @param ID A character string name of the column that contains the unique subject identifiers.
#' @param t_name A character string name of the column with integer-valued measurement time-points
#' (in days, weeks, months, etc).
#' @param x A vector containing the names of predictor variables to use for modeling.
#' If x is missing, then all columns except \code{ID}, \code{y} are used.
#' @param y A character string name of the column that represent the response variable in the model.
#' @param nfolds Number of folds to use in cross-validation.
#' @param fold_column The name of the column in the input data that contains the cross-validation
#' fold indicators (must be an ordered factor).
#' @param hold_column The name of the column that contains the holdout observation indicators
#' (TRUE/FALSE) in the input data.
#' This holdout column must be defined and added to the input data prior to calling this function.
#' @param hold_random Logical, specifying if the holdout observations should be selected at random.
#' If FALSE then the last observation for each subject is selected as a holdout.
#' @param seed Random number seed for selecting a random holdout.
#' @param use_new_features Set to \code{TRUE} to use new features (predictors) defined by the growth curve
#' feature-creator function \code{\link{define_features_drop}}.
#' Note that the \code{define_features_drop} function is called automatically, but the features
#' defined inside this function aren't use unless this is set to \code{TRUE}.
#' @param refit Set to \code{TRUE} (default) to refit the best estimator using the entire dataset.
#' When \code{FALSE}, it might be impossible to make predictions from this model fit.
#' @param verbose Set to \code{TRUE} to print messages on status and information to the console.
#' Turn this on by default using \code{options(GriDiSL.verbose=TRUE)}.
#' @param ... Additional arguments that will be passed on to \code{GriDiSL::fit_model} function.
#' @return An R6 object containing the model fit(s).
# @seealso \code{\link{GriDiSL-package}} for the general overview of the package,
# @example tests/examples/1_GriDiSL_example.R
#' @export
fit_growth.ModelStack <- function(models,
                                  method = c("none", "cv", "holdout"),
                                  data,
                                  ID,
                                  t_name,
                                  x,
                                  y,
                                  nfolds = NULL,
                                  fold_column = NULL,
                                  hold_column = NULL,
                                  hold_random = FALSE,
                                  seed = NULL,
                                  use_new_features = FALSE,
                                  refit = TRUE,
                                  verbose = getOption("GriDiSL.verbose"),
                                  ...) {

  method <- method[1L]
  gvars$method <- method
  gvars$verbose <- verbose

  feature_data_cv <- function(ID, t_name, y, train_set) {
    function(data, ...) define_features_drop(data, ID = ID, t_name = t_name, y = y, train_set = train_set, ...)
  }

  feature_data_holdout <- function(ID, t_name, y, train_set, hold_column) {
    function(data, ...) {
      if (train_set) {
        data <- data[!data[[hold_column]], ]
        define_features_drop(data, ID = ID, t_name = t_name, y = y, train_set = train_set, ...)
      } else {
        define_features_drop(data, ID = ID, t_name = t_name, y = y, train_set = train_set, hold_column = hold_column, ...)
      }
    }
  }

  if (!is.ModelStack(models)) stop("argument models must be of class 'ModelStack'")
  if (!(method %in% c("none", "cv", "holdout")))
    stop("argument method must be one of: 'none', 'cv', 'holdout'")
  if (!data.table::is.data.table(data) && !is.DataStorageClass(data))
    stop("argument data must be of class 'data.table, please convert the existing data.frame to data.table by calling 'data.table::as.data.table(...)'")

  if (missing(ID)) ID <- names(data)[1]
  if (missing(t_name)) t_name <- names(data)[2]
  if (missing(y)) y <- names(data)[3]
  if (missing(x)) x <- names(data)[4:ncol(data)]

  nodes <- list(Lnodes = x, Ynode = y, IDnode = ID, tnode = t_name)
  orig_colnames <- colnames(data)

  if (method %in% "none") {

    ## Fit models based on all available data
    data <- define_features_drop(data, ID = ID, t_name = t_name, y = y, train_set = TRUE)
    modelfit <- GriDiSL::fit_model(ID, t_name, x, y, data, models = models, verbose = verbose, ...)

  } else if (method %in% "cv") {

    if (is.null(fold_column)) {
      fold_column <- "fold"
      data <- add_CVfolds_ind(data, ID, nfolds = nfolds, fold_column = fold_column, seed = seed)
    }

    ## Define training data summaries (using all observations):
    fun_train_cv <- feature_data_cv(ID, t_name, y, train_set = TRUE)
    train_data <- fun_train_cv(data)
    ## Define validation data to be used for scoring during CV (each summary row (X_i,Y_i) is created by first dropping this row):
    fun_valid_cv <- feature_data_cv(ID, t_name, y, train_set = FALSE)
    valid_data <- fun_valid_cv(data)

    ## Add new features as predictors?
    if (use_new_features) {
      new_features <- setdiff(colnames(train_data), c(orig_colnames, hold_column))
      x <- c(x, new_features)
    }

    modelfit <- GriDiSL::fit_model(ID, t_name, x, y, train_data = train_data, models = models, fold_column = fold_column, valid_data = valid_data)
    message("...retraining the best model on all data...")
    if (refit) best_fit <- modelfit$refit_best_model(modelfit$OData_train)

  } else if (method %in% "holdout") {

    if (is.null(hold_column)) {
      hold_column <- "hold"
      message("...selecting holdout observations...")
      data <- add_holdout_ind(data, ID, hold_column = hold_column, random = hold_random, seed = seed)
    }

    ## Define training data summaries (excludes holdouts, summaries are created without the holdout observations):
    fun_train <- feature_data_holdout(ID, t_name, y, train_set = TRUE, hold_column = hold_column)
    train_data <- fun_train(data = data)

    ## Define validation data (includes the holdout only, each summary is created without the holdout observation):
    ## by giving the hold_column the non-holdout observations will be automatically dropped (could have also done it manually)
    fun_valid <- feature_data_holdout(ID, t_name, y, train_set = FALSE, hold_column = hold_column)
    valid_data <- fun_valid(data = data)

    # ## Add new features as predictors?
    if (use_new_features) {
      new_features <- setdiff(colnames(train_data), c(orig_colnames, hold_column))
      x <- c(x, new_features)
    }

    ## Perform fitting based on training set (and model scoring based on holdout validation set)
    modelfit <- GriDiSL::fit_model(ID, t_name, x, y, train_data = train_data, valid_data = valid_data, models = models)

    ## Re-fit the best scored model using all available data
    ## Define training data summaries (using all observations):
    if (refit) {
      message("...retraining the best model on all data...")
      data <- define_features_drop(data, ID = ID, t_name = t_name, y = y, train_set = TRUE)
      OData_all <- importData(data = data, ID = ID, t_name = t_name, covars = x, OUTCOME = y) ## Import input data into R6 object, define nodes
      best_fit <- modelfit$refit_best_model(OData_all)
    }

    # modelfit <- fit_holdoutSL(ID, t_name, x, y, data, models = models, hold_column = hold_column, hold_random = hold_random, refit = refit, seed = seed, verbose = verbose, ...)
  }

  return(modelfit)
}

# ---------------------------------------------------------------------------------------
#' Predict from growth-curve SuperLearner fit
#'
#' @param modelfit Model fit object returned by \code{\link{fit_growth}}.
#' @param newdata Subject-specific data for which predictions should be obtained.
#' If missing then the predictions for the training data will be typically returned.
#' See \code{holdout} for discussion of alternative cases.
#' @param add_subject_data Set to \code{TRUE} to add the subject-level data to
#' the resulting predictions (returned as a data.table).
#' When \code{FALSE} (default) only the actual predictions are returned (as a
#' matrix with each column representing predictions from a specific model).
#' @param grid Optional flag, set to \code{TRUE} when making growth
#' curve predictions for an entire grid of equally spaced time-points.
#' @param subset_idx A vector of row indices in \code{newdata} for which the
#' predictions should be obtain.
#' Default is \code{NULL} in which case all observations in \code{newdata} will
#' be used for prediction.
# @param best_refit_only Set to \code{TRUE} to obtained the predictions from the
# best scoring model that was re-trained on all observed data.
#' @param holdout Set to \code{TRUE} for out-of-sample predictions for validation
#' folds (out-of-sample observations) or holdouts.
#'  When \code{newdata} is missing there are two possible types of holdout predictions,
#' depending on the modeling approach.
#'  1. For \code{method = "holdout"} the default holdout predictions will be
#' based on validation data.
#'  2. For \code{method = "cv"} the default is to leave use the previous out-of-sample
#' (holdout) predictions from the training data.
# @param force_data.table Force the output predictions to be a \code{data.table}
#' @param verbose Set to \code{TRUE} to print messages on status and information to
#' the console. Turn this on by default using \code{options(GriDiSL.verbose=TRUE)}.
#' @return A data.table of subject level predictions (subject are rows, columns
#' are different models)
#' or a data.table with subject level covariates added along with model-based predictions.
#' @export
predict_growth <- function(modelfit,
                           newdata,
                           add_subject_data = FALSE,
                           grid = FALSE,
                           subset_idx = NULL,
                           holdout = FALSE,
                           verbose = getOption("growthcurveSL.verbose")) {

  # holdout <- FALSE
  best_only <- TRUE

  if (is.null(modelfit)) stop("must call fit_growth() before obtaining predictions")
  if (is.list(modelfit) && ("modelfit" %in% names(modelfit))) modelfit <- modelfit$modelfit
  assert_that(is.PredictionModel(modelfit) || is.PredictionStack(modelfit))
  gvars$verbose <- verbose
  nodes <- modelfit$OData_train$nodes

  if (missing(newdata) && holdout && !modelfit$runCV) {
    ## For holdout predictions with holdoutSL the default is to use the previous validation data
    # newdata <- modelfit$OData_valid$dat.sVar
    newdata <- modelfit$OData_valid
  } else if (missing(newdata)) {
    ## For all other cases, the default is to use the training data
    # newdata <- modelfit$OData_train$dat.sVar
    newdata <- modelfit$OData_train
  } else if (!grid) {
    newdata <- define_features_drop(newdata, ID = nodes$IDnode, t_name = nodes$tnode, y = nodes$Ynode, train_set = TRUE)
  } else {
    ## in the future might define the grid data internally, right now it is assumed the user correctly defines the grid data
  }

  ## Will use the best model retrained on all data for prediction:
  preds <- predict_generic(modelfit,
                           newdata,
                           # predict_only_bestK_models = 1,
                           add_subject_data = add_subject_data,
                           subset_idx = subset_idx,
                           best_only = best_only,
                           holdout = holdout,
                           force_data.table = TRUE,
                           verbose = verbose)

  # if (best_only) {
  #   best_fit_name <- names(modelfit$getRetrainedfit$model_ids)
  #   data.table::setnames(preds, old = best_fit_name, new = "SL.preds")
  # }

  ## will obtain predictions for all models in the ensemble that were trained on non-holdout observations only:
  # preds2 <- predict_model(modelfit = modelfit, newdata = newdata, add_subject_data = add_subject_data)

  return(preds)
}


# # ---------------------------------------------------------------------------------------
# #' Predict for new dataset
# #'
# #' @param modelfit Model fit object returned by \code{\link{fit_model}} function.
# #' @param newdata Subject-specific data for which predictions should be obtained.
# #' @param predict_only_bestK_models Specify the total number of top-ranked models
# #' (validation or C.V. MSE) for which predictions should be obtained.
# #' Leave missing to obtain predictions for all models that were fit as part of this ensemble.
# #' @param add_subject_data Set to \code{TRUE} to add the subject-level data to the resulting
# #' predictions (returned as a data.table).
# #' When \code{FALSE} (default) only the actual predictions are returned (as a matrix with
# #' each column representing predictions from a specific model).
# #' @param subset_idx ...
# #' @param verbose Set to \code{TRUE} to print messages on status and information to the console.
# #' Turn this on by default using \code{options(growthcurveSL.verbose=TRUE)}.
# #' @return A matrix of subject level predictions (subject are rows, columns are different models)
# #' or a data.table with subject level covariates added along with model-based predictions.
# #' @export
# predict_model <- function(modelfit,
#                           newdata,
#                           predict_only_bestK_models,
#                           add_subject_data = FALSE,
#                           subset_idx = NULL,
#                           verbose = getOption("growthcurveSL.verbose")) {
#   return(predict_generic(modelfit,
#                          newdata,
#                          predict_only_bestK_models,
#                          add_subject_data,
#                          subset_idx,
#                          best_only = FALSE,
#                          holdout = FALSE,
#                          force_data.table = TRUE,
#                          verbose))
# }

# # ---------------------------------------------------------------------------------------
# #' Prediction for holdout (out-of-sample) model
# #'
# #' When \code{newdata} is missing there are two possible types of holdout predictions, depending on the modeling approach.
# #' 1. For \code{fit_holdoutSL} the default holdout predictions will be based on validation data.
# #' 2. For \code{fit_cvSL} the default is to leave use the previous out-of-sample (holdout) predictions from the training data.
# #' @param modelfit Model fit object returned by \code{\link{fit_holdoutSL}} or \code{\link{fit_cvSL}} functions.
# #' @param newdata Subject-specific data for which holdout (out-of-sample) predictions should be obtained.
# #' @param predict_only_bestK_models Specify the total number of top-ranked models (validation or C.V. MSE) for which predictions should be obtained.
# #' Leave missing to obtain predictions for all models that were fit as part of this ensemble.
# #' @param add_subject_data Set to \code{TRUE} to add the subject-level data to the resulting predictions (returned as a data.table).
# #' When \code{FALSE} (default) only the actual predictions are returned (as a matrix with each column representing predictions from a specific model).
# #' @param subset_idx ...
# #' @param verbose Set to \code{TRUE} to print messages on status and information to the console. Turn this on by default using \code{options(growthcurveSL.verbose=TRUE)}.
# #' @return ...
# predict_holdout <- function(modelfit,
#                             newdata,
#                             predict_only_bestK_models,
#                             add_subject_data = FALSE,
#                             subset_idx = NULL,
#                             verbose = getOption("growthcurveSL.verbose")) {
#   if (missing(newdata) && !modelfit$runCV) newdata <- modelfit$OData_valid
#   return(predict_generic(modelfit, newdata, predict_only_bestK_models, add_subject_data, subset_idx,
#                          best_only = FALSE, holdout = TRUE, force_data.table = TRUE, verbose))
# }


# # ---------------------------------------------------------------------------------------
# #' Growth curve SuperLearner with one-out holdout validation
# #'
# #' Define and fit discrete SuperLearner for growth curve modeling.
# #' Model selection (scoring) is based on MSE for a single random (or last) holdout data-point for each subject.
# #' This is in contrast to the model selection with V-fold cross-validated MSE in \code{\link{fit_cvSL}},
# #' which leaves the entire subjects (entire growth curves) outside of the training sample.
# #' @param ID A character string name of the column that contains the unique subject identifiers.
# #' @param t_name A character string name of the column with integer-valued measurement time-points (in days, weeks, months, etc).
# #' @param x A vector containing the names of predictor variables to use for modeling. If x is missing, then all columns except \code{ID}, \code{y} are used.
# #' @param y A character string name of the column that represent the response variable in the model.
# #' @param data Input dataset, can be a \code{data.frame} or a \code{data.table}.
# #' @param models Parameters specifying the type of modeling procedure to be used.
# #' @param hold_column The name of the column that contains the holdout observation indicators (TRUE/FALSE) in the input data.
# #' This holdout column must be defined and added to the input data prior to calling this function.
# #' @param hold_random Logical, specifying if the holdout observations should be selected at random.
# #' If FALSE then the last observation for each subject is selected as a holdout.
# #' @param seed Random number seed for selecting a random holdout.
# # @param expr_to_train Additional logical expression which will further subset observations (rows) for training data.
# # Use this to restrict the model fitting to a specific subsample of the training datset.
# #' @param use_new_features ...
# #' @param verbose Set to \code{TRUE} to print messages on status and information to the console. Turn this on by default using \code{options(growthcurveSL.verbose=TRUE)}.
# #' @return ...
# # @seealso \code{\link{growthcurveSL-package}} for the general overview of the package,
# # @example tests/examples/1_growthcurveSL_example.R
# #' @export
# #expr_to_train = NULL,
# fit_holdoutSL <- function(ID,
#                             t_name,
#                             x,
#                             y,
#                             data,
#                             models,
#                             hold_column = NULL,
#                             hold_random = FALSE,
#                             seed = NULL,
#                             use_new_features = FALSE,
#                             verbose = getOption("growthcurveSL.verbose")) {

#   gvars$verbose <- verbose
#   nodes <- list(Lnodes = x, Ynode = y, IDnode = ID, tnode = t_name)
#   orig_colnames <- colnames(data)

#   if (is.null(hold_column)) {
#     hold_column <- "hold"
#     message("...selecting holdout observations...")
#     data <- add_holdout_ind(data, ID, hold_column = hold_column, random = hold_random, seed = seed)
#   }

#   feature_data_holdout <- function(ID, t_name, y, train_set, hold_column) {
#     function(data, ...) {
#       if (train_set) {
#         data <- data[!data[[hold_column]], ]
#         define_features_drop(data, ID = ID, t_name = t_name, y = y, train_set = train_set, ...)
#       } else {
#         define_features_drop(data, ID = ID, t_name = t_name, y = y, train_set = train_set, hold_column = hold_column, ...)
#       }
#     }
#   }

#   ## Define training data summaries (excludes holdouts, summaries are created without the holdout observations):
#   fun_train <- feature_data_holdout(ID, t_name, y, train_set = TRUE, hold_column = hold_column)
#   ## Define validation data (includes the holdout only, each summary is created without the holdout observation):
#   ## by giving the hold_column the non-holdout observations will be automatically dropped (could have also done it manually)
#   fun_valid <- feature_data_holdout(ID, t_name, y, train_set = FALSE, hold_column = hold_column)
#   train_data <- fun_train(data = data)
#   valid_data <- fun_valid(data = data)

#   ## Add new features as predictors?
#   if (use_new_features) {
#     new_features <- setdiff(colnames(train_data), c(orig_colnames, hold_column))
#     x <- c(x, new_features)
#   }

#   ## Perform fitting based on training set (and model scoring based on holdout validation set)
#   modelfit <- GriDiSL::fit_model(ID, t_name, x, y, train_data = train_data, valid_data = valid_data, models = models)

#   ## Re-fit the best scored model using all available data
#   ## Define training data summaries (using all observations):
#   message("...retraining the best model on all data...")
#   data <- define_features_drop(data, ID = ID, t_name = t_name, y = y, train_set = TRUE)
#   OData_all <- importData(data = data, ID = ID, t_name = t_name, covars = x, OUTCOME = y) ## Import input data into R6 object, define nodes
#   best_fit <- modelfit$refit_best_model(OData_all)
#   return(modelfit)
# }


# # ---------------------------------------------------------------------------------------
# #' Growth curve SuperLearner with V-fold cross-validation.
# #'
# #' Define and fit discrete SuperLearner for growth curve modeling.
# #' Model selection (scoring) is based on V-fold cross-validated MSE that leaves entire subjects outside of the training sample.
# #' This is in contrast to holdout SuperLearner in \code{\link{fit_holdoutSL}} that leaves only a single random (or last) growth measurement  outside of the training sample.
# #' @param ID A character string name of the column that contains the unique subject identifiers.
# #' @param t_name A character string name of the column with integer-valued measurement time-points (in days, weeks, months, etc).
# #' @param x A vector containing the names of predictor variables to use for modeling. If x is missing, then all columns except \code{ID}, \code{y} are used.
# #' @param y A character string name of the column that represent the response variable in the model.
# #' @param data Input dataset, can be a \code{data.frame} or a \code{data.table}.
# #' @param models Parameters specifying the type of modeling procedure to be used.
# #' @param nfolds Number of folds to use in cross-validation.
# #' @param fold_column The name of the column in the input data that contains the cross-validation fold indicators (must be an ordered factor).
# #' @param seed Random number seed for selecting a random holdout.
# # @param expr_to_train Additional logical expression which will further subset observations (rows) for training data.
# # Use this to restrict the model fitting to a specific subsample of the training dataset.
# #' @param use_new_features ...
# #' @param verbose Set to \code{TRUE} to print messages on status and information to the console. Turn this on by default using \code{options(growthcurveSL.verbose=TRUE)}.
# #' @return ...
# # @seealso \code{\link{growthcurveSL-package}} for the general overview of the package,
# # @example tests/examples/1_growthcurveSL_example.R
# #' @export
# # expr_to_train = NULL,
# fit_cvSL <- function(ID,
#                      t_name,
#                      x,
#                      y,
#                      data,
#                      models,
#                      nfolds = 5,
#                      fold_column = NULL,
#                      seed = NULL,
#                      use_new_features = FALSE,
#                      verbose = getOption("growthcurveSL.verbose")) {

#   gvars$verbose <- verbose
#   nodes <- list(Lnodes = x, Ynode = y, IDnode = ID, tnode = t_name)
#   orig_colnames <- colnames(data)

#   if (is.null(fold_column)) {
#     fold_column <- "fold"
#     data <- add_CVfolds_ind(data, ID, nfolds = nfolds, fold_column = fold_column, seed = seed)
#   }
#   # models$fold_column <- fold_column

#   feature_data_cv <- function(ID, t_name, y, train_set) {
#     function(data, ...) define_features_drop(data, ID = ID, t_name = t_name, y = y, train_set = train_set, ...)
#   }

#   ## Define training data summaries (using all observations):
#   fun_train_cv <- feature_data_cv(ID, t_name, y, train_set = TRUE)
#   ## Define validation data to be used for scoring during CV (each summary row (X_i,Y_i) is created by first dropping this row):
#   fun_valid_cv <- feature_data_cv(ID, t_name, y, train_set = FALSE)

#   train_data_2 <- fun_train_cv(data)
#   valid_data_2 <- fun_valid_cv(data)

#   # train_data <- define_features(data, nodes, train_set = TRUE, holdout = FALSE)
#   train_data <- define_features_drop(data, ID = ID, t_name = t_name, y = y, train_set = TRUE)
#   # if (!is.null(expr_to_train)) train_data <- train_data[eval(parse(text=expr_to_train)), ]
#   all.equal(train_data_2, train_data)
#   ## Define validation data to be used for scoring during CV (each summary row (X_i,Y_i) is created by first dropping this row):
#   valid_data <- define_features_drop(data, ID = ID, t_name = t_name, y = y, train_set = FALSE)
#   all.equal(valid_data_2, valid_data)


#   ## ------------------------------------------------------------------------------------------
#   ## Add new features as predictors?
#   ## ------------------------------------------------------------------------------------------
#   if (use_new_features) {
#     new_features <- setdiff(colnames(train_data), c(orig_colnames, fold_column))
#     x <- c(x, new_features)
#   }

#   browser()

#   ## ------------------------------------------------------------------------------------------
#   ## Perform CV fitting (no scoring yet)
#   ## ------------------------------------------------------------------------------------------
#   # modelfit <- GriDiSL::fit_model(ID, t_name, x, y, train_data = train_data, models = models, fold_column = fold_column)

# # Obtaining pre-saved out-of-sample/holdout CV predictions for h2o
# # Obtaining pre-saved out-of-sample/holdout CV predictions for h2o
# # Obtaining pre-saved out-of-sample/holdout CV predictions for h2o
# # [1] "Internally evaluated holdout / CV metrics: "
# #         MSE    MSE.sd     RMSE              model Model_idx order
# # 1: 1.487644 0.1126921 1.219690 M.3.h2o.gbm.grid.1         3     1
# # 2: 1.535146 0.1051345 1.239010        M.1.h2o.glm         1     1
# # 3: 1.535869 0.1049170 1.239302 M.2.h2o.glm.grid.1         2     1
# # 4: 1.536187 0.1127723 1.239430 M.3.h2o.gbm.grid.2         3     2
# # 5: 1.537622 0.1052653 1.240009 M.2.h2o.glm.grid.2         2     2
# # 6: 1.540582 0.1049696 1.241202 M.2.h2o.glm.grid.3         2     3

# # Obtaining pre-saved out-of-sample/holdout CV predictions for h2o
# # Obtaining pre-saved out-of-sample/holdout CV predictions for h2o
# # [1] "Internally evaluated holdout / CV metrics: "
# #         MSE    MSE.sd     RMSE              model Model_idx order
# # 1: 1.458157 0.1067827 1.207542 M.2.h2o.gbm.grid.1         2     1
# # 2: 1.484603 0.1107522 1.218443 M.2.h2o.gbm.grid.2         2     2
# # 3: 1.537324 0.1052537 1.239889 M.1.h2o.glm.grid.1         1     1
# # 4: 1.537618 0.1052652 1.240007 M.1.h2o.glm.grid.2         1     2
# # 5: 1.537622 0.1052653 1.240009 M.1.h2o.glm.grid.3         1     3

#   modelfit <- GriDiSL::fit_model(ID, t_name, x, y, train_data = train_data, models = models, fold_column = fold_column, valid_data = valid_data)

# # [1] "Internally evaluated holdout / CV metrics: "
# #         MSE    MSE.sd     RMSE              model Model_idx order
# # 1: 1.455267 0.1074002 1.206344 M.2.h2o.gbm.grid.1         2     1
# # 2: 1.467470 0.1103306 1.211392 M.2.h2o.gbm.grid.2         2     2
# # 3: 1.534744 0.1049081 1.238848 M.1.h2o.glm.grid.1         1     1
# # 4: 1.546036 0.1050662 1.243397 M.1.h2o.glm.grid.2         1     2
# # 5: 1.549218 0.1049876 1.244676 M.1.h2o.glm.grid.3         1     3

#   ## ------------------------------------------------------------------------------------------
#   ## Score CV models based on validation set
#   ## ------------------------------------------------------------------------------------------
#   # valid_data <- importData(data = valid_data, ID = nodes$IDnode, t_name = nodes$tnode, covars = modelfit$predvars, OUTCOME = modelfit$outvar)
#   # modelfit$OData_valid <- valid_data

#   # # modelfit <- modelfit$score_models(validation_data = valid_data) # returns the modelfit object intself, but does the scoring of each CV model
#   # modelfit$score_models(validation_data = valid_data)
#   # print("CV MSE after manual CV model rescoring: "); print(unlist(modelfit$getMSE))
#   # preds <- predict_CV(modelfit, valid_data) # will return the matrix of predicted cv values (one column per model)
#   # [1] "CV MSE after manual CV model rescoring: "
#   # M.1.h2o.glm.grid.1 M.1.h2o.glm.grid.2 M.1.h2o.glm.grid.3 M.2.h2o.gbm.grid.1 M.2.h2o.gbm.grid.2
#   #           1.534744           1.546036           1.549218           1.455267           1.467470

#   ## ------------------------------------------------------------------------------------------
#   ## Re-fit the best manually-scored model on all data
#   ## Even though we don't need to do this for CV (best model is already trained), we do this to be consistent with holdoutSL (and for additional error checking)
#   ## ------------------------------------------------------------------------------------------
#   best_fit <- modelfit$refit_best_model(modelfit$OData_train)

#   return(modelfit)
# }

# # ---------------------------------------------------------------------------------------
# #' Evaluate MSE for model fits, possibly using new data
# #'
# #' By default this function will extract out-of-sample predictions from original training data (automatically done by h2o) to evaluate the cross-validated MSE.
# #' However, when \code{newdata} is supplied, the predictions for each CV model will be based on this external validation dataset.
# #' These predictions and the outcome stored in \code{newdata} are then used to re-evalute the CV MSE. Note that \code{newdata} must be of the same
# #' dimensionality as the original training data used for fitting the h2o models.
# #' @param modelfit Model fit object returned by \code{\link{fit_model}} function.
# #' @param newdata ...
# #' @param subset_idx ...
# #' @param verbose Set to \code{TRUE} to print messages on status and information to the console. Turn this on by default using \code{options(growthcurveSL.verbose=TRUE)}.
# #' @return ...
# #' @export
# eval_MSE <- function(modelfit, newdata, subset_idx = NULL, verbose = getOption("growthcurveSL.verbose")) {
#   if (is.list(modelfit) && ("modelfit" %in% names(modelfit))) modelfit <- modelfit$modelfit
#   if (is.null(modelfit)) stop("must call get_fit() prior to obtaining predictions")
#   assert_that(is.PredictionModel(modelfit))
#   gvars$verbose <- verbose
#   nodes <- modelfit$OData_train$nodes

#   if (!is.null(subset_idx))
#     if (!is.integer(subset_idx)) stop("subset_idx must be an integer vector, current class: " %+% class(subset_idx))

#   if (missing(newdata)) {
#     ## Use out-of-sample predictions from original training data (automatically done by h2o) to evaluate the CV MSE
#     modelfit <- modelfit$score_models(subset_exprs = subset_idx)
#   } else {
#     newdata <- validate_convert_input_data(newdata, ID = nodes$IDnode, t_name = nodes$tnode,
#                                            x = modelfit$predvars, y = modelfit$outvar,
#                                            useH2Oframe = modelfit$useH2Oframe, dest_frame = "prediction_H2Oframe")
#     ## Get predictions for each CV model based on external validation CV dataset, then use those predictions and outcome in newdata to evalute the CV MSE
#     modelfit <- modelfit$score_models(validation_data = newdata, subset_exprs = subset_idx)
#   }
#   return(modelfit$getMSE)
# }