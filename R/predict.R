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
#' the console.
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

  return(preds)
}




# ---------------------------------------------------------------------------------------
#' Obtain all possible predictions with one call
#'
#' This function creates a \code{data.frame} with one row per subject and a column named "fits".
#' The latter column contains the 3 possible subject-level predictions:  fit - training data model predictions,
#' fitgrid - grid-based model predictions, holdout - predictions for holdouts or out-of-sample predictions.
#' Some of these predictions might be empty, e.g. when \code{add_holdout = FALSE}, no holdout predictions are made.
#' @param modelfit Model fit object returned by \code{\link{fit_growth}}.
#' @param newdata Subject-specific data for which predictions should be obtained.
#' @param add_holdout Optional flag, set to \code{TRUE} to add holdout (out-of-sample) predictions
#' used for scoring the current best model fit.
#' @param add_grid Optional flag, set to \code{TRUE} to add a grid of equally spaced
#' predictions over some range of the time variable for the current best model fit.
#' @param tmin Min t value of the grid
#' @param tmax Max t value of the grid
#' @param incr Increment value for the grid of \code{t}'s
#' @param verbose Set to \code{TRUE} to print messages on status and information to
#' the console.
#' @return A data.frame with one row per subject.
#' The column named "fit" contains nested subject-specific observed data and predictions.
#' Each entry (cell) of "fit" is a list of 4 subject-specific data.frames,
#' named "xy", "fit", "fitgrid" and "holdout",
#' containing the training time-outcome values, the model predictions for training data,
#' the grid of equally spaced model
#' prediction and the holdout / out-of-sample predictions, respectively.
#' @export
predict_all <- function(modelfit,
                        newdata,
                        add_holdout = FALSE,
                        add_grid = FALSE,
                        tmin = 1,
                        tmax = 500,
                        incr = 2,
                        verbose = getOption("growthcurveSL.verbose")) {

  nodes <- modelfit$OData_train$nodes
  ID <- nodes$IDnode
  t_name <- nodes$tnode
  y <- nodes$Ynode

  empty_df = data.frame(matrix(vector(), 0, 3,
                        dimnames=list(c(), c(ID, t_name, "preds"))),
                        stringsAsFactors=FALSE) %>%
              data.table::data.table()

  preds_holdout <- preds_grid <- empty_df

  ## generate data with all unique subject IDs that appear in input training data
  unique_subj <- newdata %>%
                 dplyr::distinct_(ID)

  ## Generate data with one row per subj ID
  ## Each row in column "xy" is itself a data.frame
  ## Each "xy" is a data.frame containing all observed (x, y) = (time, outcome) values on one subject
  obs_bysujb <- newdata %>%
                dplyr::rename_("x" = t_name, "y" = y) %>%
                dplyr::select_(ID, "x", "y") %>%
                dplyr::group_by_(ID) %>%
                tidyr::nest(.key = "xy")


  ## Predictions for new data based on best SL model trained on all data:
  preds_alldat <- predict_growth(modelfit, newdata = newdata, add_subject_data = TRUE)

  fit_bysujb <- preds_alldat %>%
                dplyr::rename_("x" = t_name, "y" = "preds") %>%
                dplyr::group_by_(ID) %>%
                tidyr::nest(.key = "fit")


  ## Predictions for all holdout data points for all models trained on non-holdout data:
  if (add_holdout)
    preds_holdout <- predict_growth(modelfit, holdout = TRUE, add_subject_data = TRUE)

  hold_bysujb <-  unique_subj %>%
                  dplyr::left_join(preds_holdout) %>%
                  dplyr::rename_("x" = t_name, "y" = "preds") %>%
                  dplyr::group_by_(ID) %>%
                  tidyr::nest(.key = "holdout")


  ## Predictions for a grid of equally spaced time points
  if (add_grid) {
    train_dat <- define_features_drop(newdata, ID =  ID, t_name = t_name, y = y, train_set = TRUE)
    grid_dat <- define_tgrid(train_dat, ID = ID, t_name = t_name, y = y, tmin = tmin, tmax = tmax, incr = incr)
    preds_grid <- predict_growth(modelfit, newdata = grid_dat, grid = TRUE, add_subject_data = TRUE)
  }

  fitgrid_bysujb <- unique_subj %>%
                    dplyr::left_join(preds_grid) %>%
                    dplyr::rename_("x" = t_name, "y" = "preds") %>%
                    dplyr::group_by_(ID) %>%
                    tidyr::nest(.key = "fitgrid")


  ## combine all 3 datasets:
  fits_all <- observed_bysujb %>%
              dplyr::left_join(fit_bysujb) %>%
              dplyr::left_join(fitgrid_bysujb) %>%
              dplyr::left_join(hold_bysujb)

  fits_all <- fits_all %>%
              dplyr::group_by_(ID) %>%
              tidyr::nest(.key = "fit")

  return(fits_all)
}
