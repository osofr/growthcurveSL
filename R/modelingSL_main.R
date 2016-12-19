#' @import R6
#' @import data.table
#' @import ggiraph
#' @import longGriDiSL
#' @importFrom graphics axis barplot hist par text  legend plot
#' @importFrom methods is
#' @importFrom stats approx binomial gaussian coef glm.control glm.fit plogis predict qlogis qnorm quantile rnorm terms var predict glm.control
#' @importFrom utils data head str
#' @importFrom stats as.formula glm na.exclude rbinom terms.formula
#' @importFrom assertthat assert_that
#' @importFrom ggplot2 ggplot geom_point geom_errorbar theme_bw coord_flip aes position_dodge alpha
# @importFrom ggiraph geom_point_interactive ggiraph
NULL

# ---------------------------------------------------------------------------------------
#' Growth curve SuperLearner with one-out holdout validation
#'
#' Define and fit discrete SuperLearner for growth curve modeling.
#' Model selection (scoring) is based on MSE for a single random (or last) holdout data-point for each subject.
#' This is in contrast to the model selection with V-fold cross-validated MSE in \code{\link{fit_cvSL}},
#' which leaves the entire subjects (entire growth curves) outside of the training sample.
#' @param ID A character string name of the column that contains the unique subject identifiers.
#' @param t_name A character string name of the column with integer-valued measurement time-points (in days, weeks, months, etc).
#' @param x A vector containing the names of predictor variables to use for modeling. If x is missing, then all columns except \code{ID}, \code{y} are used.
#' @param y A character string name of the column that represent the response variable in the model.
#' @param data Input dataset, can be a \code{data.frame} or a \code{data.table}.
#' @param params Parameters specifying the type of modeling procedure to be used.
#' @param hold_column The name of the column that contains the holdout observation indicators (TRUE/FALSE) in the input data.
#' This holdout column must be defined and added to the input data prior to calling this function.
#' @param random Logical, specifying if the holdout observations should be selected at random.
#' If FALSE then the last observation for each subject is selected as a holdout.
#' @param seed Random number seed for selecting a random holdout.
#' @param expr_to_train Additional logical expression which will further subset observations (rows) for training data.
#' Use this to restrict the model fitting to a specific subsample of the training datset.
#' @param use_new_features ...
#' @param verbose Set to \code{TRUE} to print messages on status and information to the console. Turn this on by default using \code{options(growthcurveSL.verbose=TRUE)}.
#' @return ...
# @seealso \code{\link{growthcurveSL-package}} for the general overview of the package,
# @example tests/examples/1_growthcurveSL_example.R
#' @export
fit_holdoutSL <- function(ID, t_name, x, y, data, params, hold_column = NULL, random = FALSE, seed = NULL, expr_to_train = NULL, use_new_features = FALSE, verbose = getOption("growthcurveSL.verbose")) {
  gvars$verbose <- verbose
  nodes <- list(Lnodes = x, Ynode = y, IDnode = ID, tnode = t_name)
  orig_colnames <- colnames(data)

  if (is.null(hold_column)) {
    hold_column <- "hold"
    message("...selecting holdout observations...")
    data <- add_holdout_ind(data, ID, hold_column = hold_column, random = random, seed = seed)
  }

  # params$hold_column <- hold_column

  train_data <- data[!data[[hold_column]], ]

  ## ------------------------------------------------------------------------------------------
  ## Define training data summaries (excludes holdouts, summaries are created without the holdout observations):
  ## ------------------------------------------------------------------------------------------
  train_data <- define_features_drop(train_data, ID = ID, t_name = t_name, y = y, train_set = TRUE)
  if (!is.null(expr_to_train)) train_data <- train_data[eval(parse(text=expr_to_train)), ]

  ## ------------------------------------------------------------------------------------------
  ## Define validation data (includes the holdout only, each summary is created without the holdout observation):
  ## ------------------------------------------------------------------------------------------
  # by giving the hold_column the non-holdout observations will be automatically dropped (could have also done it manually)
  valid_data <- define_features_drop(data, ID = ID, t_name = t_name, y = y, train_set = FALSE, hold_column = hold_column)

  ## ------------------------------------------------------------------------------------------
  ## Add new features as predictors?
  ## ------------------------------------------------------------------------------------------
  if (use_new_features) {
    new_features <- setdiff(colnames(train_data), c(orig_colnames, hold_column))
    x <- c(x, new_features)
  }

  ## ------------------------------------------------------------------------------------------
  ## Perform fitting based on training set (and model scoring based on holdout validation set)
  ## ------------------------------------------------------------------------------------------
  modelfit <- fit_model(ID, t_name, x, y, train_data = train_data, valid_data = valid_data, params = params)

  ## ------------------------------------------------------------------------------------------
  ## Re-fit the best scored model using all available data
  ## ------------------------------------------------------------------------------------------
  ## Define training data summaries (using all observations):
  data <- define_features_drop(data, ID = ID, t_name = t_name, y = y, train_set = TRUE)
  OData_all <- importData(data = data, ID = ID, t_name = t_name, covars = x, OUTCOME = y) ## Import input data into R6 object, define nodes
  best_fit <- modelfit$refit_best_model(OData_all)
  return(modelfit)
}

# ---------------------------------------------------------------------------------------
#' Growth curve SuperLearner with V-fold cross-validation.
#'
#' Define and fit discrete SuperLearner for growth curve modeling.
#' Model selection (scoring) is based on V-fold cross-validated MSE that leaves entire subjects outside of the training sample.
#' This is in contrast to holdout SuperLearner in \code{\link{fit_holdoutSL}} that leaves only a single random (or last) growth measurement  outside of the training sample.
#' @param ID A character string name of the column that contains the unique subject identifiers.
#' @param t_name A character string name of the column with integer-valued measurement time-points (in days, weeks, months, etc).
#' @param x A vector containing the names of predictor variables to use for modeling. If x is missing, then all columns except \code{ID}, \code{y} are used.
#' @param y A character string name of the column that represent the response variable in the model.
#' @param data Input dataset, can be a \code{data.frame} or a \code{data.table}.
#' @param params Parameters specifying the type of modeling procedure to be used.
#' @param nfolds Number of folds to use in cross-validation.
#' @param fold_column The name of the column in the input data that contains the cross-validation fold indicators (must be an ordered factor).
#' @param seed Random number seed for selecting a random holdout.
#' @param expr_to_train Additional logical expression which will further subset observations (rows) for training data.
#' Use this to restrict the model fitting to a specific subsample of the training dataset.
#' @param use_new_features ...
#' @param verbose Set to \code{TRUE} to print messages on status and information to the console. Turn this on by default using \code{options(growthcurveSL.verbose=TRUE)}.
#' @return ...
# @seealso \code{\link{growthcurveSL-package}} for the general overview of the package,
# @example tests/examples/1_growthcurveSL_example.R
#' @export
fit_cvSL <- function(ID, t_name, x, y, data, params, nfolds = 5, fold_column = NULL, seed = NULL, expr_to_train = NULL, use_new_features = FALSE, verbose = getOption("growthcurveSL.verbose")) {
  gvars$verbose <- verbose
  nodes <- list(Lnodes = x, Ynode = y, IDnode = ID, tnode = t_name)
  orig_colnames <- colnames(data)

  if (is.null(fold_column)) {
    fold_column <- "fold"
    data <- add_CVfolds_ind(data, ID, nfolds = nfolds, fold_column = fold_column, seed = seed)
  }

  # params$fold_column <- fold_column

  ## ------------------------------------------------------------------------------------------
  ## Define training data summaries (using all observations):
  ## ------------------------------------------------------------------------------------------
  # train_data <- define_features(data, nodes, train_set = TRUE, holdout = FALSE)
  train_data <- define_features_drop(data, ID = ID, t_name = t_name, y = y, train_set = TRUE)
  if (!is.null(expr_to_train)) train_data <- train_data[eval(parse(text=expr_to_train)), ]

  ## ------------------------------------------------------------------------------------------
  ## Add new features as predictors?
  ## ------------------------------------------------------------------------------------------
  if (use_new_features) {
    new_features <- setdiff(colnames(train_data), c(orig_colnames, fold_column))
    x <- c(x, new_features)
  }

  ## ------------------------------------------------------------------------------------------
  ## Perform CV fitting (no scoring yet)
  ## ------------------------------------------------------------------------------------------
  modelfit <- fit_model(ID, t_name, x, y, train_data = train_data, params = params, fold_column = fold_column)

  ## ------------------------------------------------------------------------------------------
  ## Score CV models based on validation set
  ## ------------------------------------------------------------------------------------------
  ## Define validation data to be used for scoring during CV (each summary row (X_i,Y_i) is created by first dropping this row):
  valid_data <- define_features_drop(data, ID = ID, t_name = t_name, y = y, train_set = FALSE)
  valid_data <- importData(data = valid_data, ID = nodes$IDnode, t_name = nodes$tnode, covars = modelfit$predvars, OUTCOME = modelfit$outvar)
  modelfit$OData_valid <- valid_data

  # modelfit <- modelfit$score_models(validation_data = valid_data) # returns the modelfit object intself, but does the scoring of each CV model
  modelfit$score_models(validation_data = valid_data)
  print("CV MSE after manual CV model rescoring: "); print(unlist(modelfit$getMSE))
  # preds <- predict_CV(modelfit, valid_data) # will return the matrix of predicted cv values (one column per model)

  ## ------------------------------------------------------------------------------------------
  ## Re-fit the best manually-scored model on all data
  ## Even though we don't need to do this for CV (best model is already trained), we do this to be consistent with holdoutSL (and for additional error checking)
  ## ------------------------------------------------------------------------------------------
  best_fit <- modelfit$refit_best_model(modelfit$OData_train)

  return(modelfit)
}


# ---------------------------------------------------------------------------------------
#' Predict for new dataset
#'
#' @param modelfit Model fit object returned by \code{\link{fit_model}} function.
#' @param newdata Subject-specific data for which predictions should be obtained.
#' @param predict_only_bestK_models Specify the total number of top-ranked models (validation or C.V. MSE) for which predictions should be obtained.
#' Leave missing to obtain predictions for all models that were fit as part of this ensemble.
#' @param add_subject_data Set to \code{TRUE} to add the subject-level data to the resulting predictions (returned as a data.table).
#' When \code{FALSE} (default) only the actual predictions are returned (as a matrix with each column representing predictions from a specific model).
#' @param subset_idx ...
#' @param verbose Set to \code{TRUE} to print messages on status and information to the console.
#' Turn this on by default using \code{options(growthcurveSL.verbose=TRUE)}.
#' @return A matrix of subject level predictions (subject are rows, columns are different models)
#' or a data.table with subject level covariates added along with model-based predictions.
#' @export
predict_model <- function(modelfit, newdata, predict_only_bestK_models,
                          add_subject_data = FALSE,
                          subset_idx = NULL,
                          verbose = getOption("growthcurveSL.verbose")) {
  return(predict_generic(modelfit, newdata, predict_only_bestK_models, add_subject_data, subset_idx,
                         use_best_retrained_model = FALSE, pred_holdout = FALSE, force_data.table = TRUE, verbose))
}

# ---------------------------------------------------------------------------------------
#' Prediction for holdout (out-of-sample) model
#'
#' When \code{newdata} is missing there are two possible types of holdout predictions, depending on the modeling approach.
#' 1. For \code{fit_holdoutSL} the default holdout predictions will be based on validation data.
#' 2. For \code{fit_cvSL} the default is to leave use the previous out-of-sample (holdout) predictions from the training data.
#' @param modelfit Model fit object returned by \code{\link{fit_holdoutSL}} or \code{\link{fit_cvSL}} functions.
#' @param newdata Subject-specific data for which holdout (out-of-sample) predictions should be obtained.
#' @param predict_only_bestK_models Specify the total number of top-ranked models (validation or C.V. MSE) for which predictions should be obtained.
#' Leave missing to obtain predictions for all models that were fit as part of this ensemble.
#' @param add_subject_data Set to \code{TRUE} to add the subject-level data to the resulting predictions (returned as a data.table).
#' When \code{FALSE} (default) only the actual predictions are returned (as a matrix with each column representing predictions from a specific model).
#' @param subset_idx ...
#' @param verbose Set to \code{TRUE} to print messages on status and information to the console. Turn this on by default using \code{options(growthcurveSL.verbose=TRUE)}.
#' @return ...
predict_holdout <- function(modelfit, newdata, predict_only_bestK_models,
                            add_subject_data = FALSE,
                            subset_idx = NULL,
                            verbose = getOption("growthcurveSL.verbose")) {
  if (missing(newdata) && !modelfit$runCV) newdata <- modelfit$OData_valid
  return(predict_generic(modelfit, newdata, predict_only_bestK_models, add_subject_data, subset_idx,
                         use_best_retrained_model = FALSE, pred_holdout = TRUE, force_data.table = TRUE, verbose))
}

# ---------------------------------------------------------------------------------------
#' Predict for holdout or curve SuperLearner fits
#'
#' @param modelfit Model fit object returned by \code{\link{fit_holdoutSL}} or  \code{\link{fit_cvSL}}.
#' @param newdata Subject-specific data for which predictions should be obtained.
#' @param add_subject_data Set to \code{TRUE} to add the subject-level data to the resulting predictions (returned as a data.table).
#' When \code{FALSE} (default) only the actual predictions are returned (as a matrix with each column representing predictions from a specific model).
#' @param grid ...
#' @param subset_idx ...
#' @param verbose Set to \code{TRUE} to print messages on status and information to the console. Turn this on by default using \code{options(growthcurveSL.verbose=TRUE)}.
#' @return ...
#' @export
# predict_SL <- function(modelfit, newdata,
#                       add_subject_data = FALSE,
#                       grid = FALSE,
#                       verbose = getOption("growthcurveSL.verbose")) {
predict_SL <- function(modelfit, newdata,
                       add_subject_data = FALSE,
                       grid = FALSE,
                       subset_idx = NULL,
                       verbose = getOption("growthcurveSL.verbose")) {

  pred_holdout <- FALSE
  use_best_retrained_model <- TRUE

  if (is.null(modelfit)) stop("must call fit_holdoutSL() or fit_cvSL() prior to obtaining predictions")
  if (is.list(modelfit) && ("modelfit" %in% names(modelfit))) modelfit <- modelfit$modelfit
  assert_that(is.PredictionModel(modelfit))
  gvars$verbose <- verbose
  nodes <- modelfit$OData_train$nodes

  if (missing(newdata) && pred_holdout && !modelfit$runCV) {
    ## For holdout predictions with holdoutSL the default is to use the previous validation data
    newdata <- modelfit$OData_valid$dat.sVar
  } else if (missing(newdata)) {
    ## For all other cases, the default is to use the training data
    newdata <- modelfit$OData_train$dat.sVar
  } else if (!grid) {
    newdata <- define_features_drop(newdata, ID = nodes$IDnode, t_name = nodes$tnode, y = nodes$Ynode, train_set = TRUE)
  } else {
    ## in the future might define the grid data internally, right now it is assumed the user correctly defines the grid data
  }

  ## Will use the best model retrained on all data for prediction:
  # modelfit$use_best_retrained_model <- TRUE
  # preds <- predict_model(modelfit = modelfit, newdata = newdata, add_subject_data = add_subject_data)
  preds <- predict_generic(modelfit, newdata,
                           predict_only_bestK_models = 1,
                           add_subject_data = add_subject_data,
                           subset_idx = subset_idx,
                           use_best_retrained_model = use_best_retrained_model,
                           pred_holdout = pred_holdout,
                           force_data.table = TRUE,
                           verbose = verbose)

  if (use_best_retrained_model) {
    best_fit_name <- names(modelfit$getRetrainedfit$model_ids)
    data.table::setnames(preds, old = best_fit_name, new = "SL.preds")
    # modelfit$use_best_retrained_model <- FALSE
  }

  ## will obtain predictions for all models in the ensemble that were trained on non-holdout observations only:
  # preds2 <- predict_model(modelfit = modelfit, newdata = newdata, add_subject_data = add_subject_data)

  return(preds)
}

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