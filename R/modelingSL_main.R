#' @import R6
#' @import data.table
#' @import ggiraph
#' @import GriDiSL
#' @importFrom magrittr %>%
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

