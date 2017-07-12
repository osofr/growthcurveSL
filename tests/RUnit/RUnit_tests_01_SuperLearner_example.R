notest.install.software <- function() {
  ## ------------------------------------------------------------------------
  ## growthcurveSL R package is dependent on a number of open-source R packages.
  ## Some of these packages are not on CRAN
  ## and need to installed directly from github repositories.
  ## ------------------------------------------------------------------------

  ## ------------------------------------------------------------------------
  ## We first install xgboost machine learning toolkit, which we be a part
  ## of our growth curve SuperLearner.
  ## ------------------------------------------------------------------------
  install.packages('xgboost')

  ## ------------------------------------------------------------------------
  ## Next, we install h2o machine learning toolkit, which will also be a part
  ## of our growth curve SuperLearner.
  ## ------------------------------------------------------------------------
  if ("package:h2o" %in% search()) detach("package:h2o", unload=TRUE)
  if ("h2o" %in% rownames(installed.packages())) remove.packages("h2o")
  ## First, we download and install H2O package dependencies:
  pkgs <- c("methods","statmod","stats","graphics","RCurl","jsonlite","tools","utils")
  new.pkgs <- setdiff(pkgs, rownames(installed.packages()))
  if (length(new.pkgs)) install.packages(new.pkgs)
  ## Next, we download and install the h2o R package itself:
  install.packages("h2o", type="source", repos=(c("http://h2o-release.s3.amazonaws.com/h2o/rel-tutte/2/R")))

  ## ------------------------------------------------------------------------
  ## Next, we install brokenstick and face R packages for growth curve modeling:
  ## ------------------------------------------------------------------------
  options(repos = c(
    CRAN = "http://cran.rstudio.com/",
    tessera = "http://packages.tessera.io"))
  install.packages("brokenstick")
  install.packages("face")

  ## ------------------------------------------------------------------------
  ## Installing hbgd R package:
  ## ------------------------------------------------------------------------
  devtools::install_github('hafen/hbgd', ref = "tidy")

  ## ------------------------------------------------------------------------
  ## Installing trelliscopejs for the visualization of the imputed growth trajectories:
  ## ------------------------------------------------------------------------
  devtools::install_github("hafen/trelliscopejs")

  ## ------------------------------------------------------------------------
  ## Finally, below we show how to  install gridisl and growthcurveSL R packages:
  ## ------------------------------------------------------------------------
  devtools::install_github('osofr/gridisl', build_vignettes = FALSE)
  devtools::install_github('osofr/growthcurveSL', build_vignettes = FALSE)

}

test.combine.all.model.fits <- function() {
  library("magrittr")
  library("dplyr")
  library("h2o")
  library("xgboost")
  library("gridisl")
  library("growthcurveSL")
  options(growthcurveSL.verbose = TRUE)
  options(gridisl.verbose = TRUE)

  # options(growthcurveSL.verbose = FALSE)
  # options(gridisl.verbose = FALSE)

  data(cpp)
  cpp <- cpp[!is.na(cpp[, "haz"]), ]
  covars <- c("apgar1", "apgar5", "parity", "gagebrth", "mage", "meducyrs", "sexn")

  ## We start by adding an random indicator for holdout observations.
  ## The holdout SuperLearner will use these holdouts for model comparison.
  # cpp_holdout <- add_holdout_ind(data = cpp, ID = "subjid", hold_column = "hold", random = TRUE, seed = 54321)
  cpp <- add_holdout_ind(data = cpp, ID = "subjid", hold_column = "hold", random = TRUE, seed = 54321)

  ## Similarly, we define the random validation folds (5 folds),
  ## separating each 5th of each subjects into a separate validation fold.
  ## The cross-validated SuperLearner will use each validation fold for model comparison.
  # cpp_folds <- add_CVfolds_ind(cpp, ID = "subjid", nfolds = 5, seed = 23)
  cpp <- add_CVfolds_ind(cpp, ID = "subjid", nfolds = 5, seed = 23)

  ## --------------------------------------------------------------------------------------------
  ## Fitting Growth Trajectories with Random Holdout SuperLearner
  ## --------------------------------------------------------------------------------------------
  ## First we define the grid of hyper parameters for h2o GBM
  h2o_GBM_hyper <- list(
    ntrees = c(20, 50, 100),
    learn_rate = c(.05, .1, .2),
    max_depth = c(3, 6, 10, 15),
    sample_rate = c(.5, .75, .9, 1),
    col_sample_rate = seq(0.5, 1),
    col_sample_rate_per_tree = c(.3, .4, .8, 1)
  )
  ## Similarly, we define the grid of hyper parameters for xgboost
  xgb_GBM_hyper = list(
    nrounds = c(20, 50, 100),
    learning_rate = c(.05, .1, .2),
    max_depth = c(3, 6, 10, 15),
    subsample = c(.5, .75, .9, 1),
    colsample_bytree = c(.3, .4, .8, 1),
    min_child_weight = c(1, 5, 7),
    gamma = c(.0, .05, seq(.1, .9, by=.2), 1),
    lambda = c(.1, .5, 1, 2, 5),
    alpha = c(0, .1, .5, .8, 1)
    )

  ## --------------------------------------------------------------------------------------------
  ## Define an ensemble of learners (SuperLearner) with h2o,  xgboost and brokenstick using the novel
  ## gridisl R package syntax.
  ## Note that any number of learners can be added with "+" syntax.
  ## By setting strategy = "RandomDiscrete", we specify that the models should be selected
  ## at random from the above defined grids of model parameters.
  ## By setting max_models = 10, we specify that at most 10 such models should be considered for h2o grid
  ## (similarly, for xgboost, where we consider 30 randomly drawn model parameters).
  ## Note that for best results max_models should be set substantially higher than 10 or 30 (computational resources permitting).
  ## By setting max_models to higher values we can explore a larger space of tuning parameters,
  ## increasing the chance that we actually find the best performing model in the grid
  ## (i.e., finding the most generalizable model).
  ## As part of our ensemble we also include the brokenstick model (added as the first model).
  ## --------------------------------------------------------------------------------------------
  grid_holdSL <-
    defModel(estimator = "brokenstick__brokenstick", predict.w.Y = FALSE) +

    defModel(estimator = "h2o__gbm", family = "gaussian",
             search_criteria = list(strategy = "RandomDiscrete", max_models = 5),
             param_grid = h2o_GBM_hyper,
             stopping_rounds = 4, stopping_metric = "MSE",
             seed = 123456) +

    defModel(estimator = "xgboost__gbm", family = "gaussian",
             search_criteria = list(strategy = "RandomDiscrete", max_models = 5),
             # search_criteria = list(strategy = "RandomDiscrete", max_models = 30),
             param_grid = xgb_GBM_hyper,
             early_stopping_rounds = 4,
             seed = 123456)

  ## --------------------------------------------------------------------------------------------
  ## Prior to training the model with the SuperLearner, we need to initialize the h2o cluster.
  ## This step is necessarily for modeling with h2o machine learning toolkit.
  ## --------------------------------------------------------------------------------------------
  h2o::h2o.init(nthreads = -1)

  ## --------------------------------------------------------------------------------------------
  ## Fit the discrete holdout SuperLearner using random holdout observations for selecting best model.
  ## The holdout SuperLearner is enabled by specifying the argument method = "holdout".
  ## By setting use_new_features=TRUE, we allow the fitting procedures to use additional
  ## summaries of growth trajectories as predictors.
  ## --------------------------------------------------------------------------------------------
  mfit_holdSL <- fit_growth(grid_holdSL,
                           ID = "subjid",
                           t_name = "agedays",
                           x = c("agedays", covars),
                           y = "haz",
                           # data = cpp_holdout,
                           data = cpp,
                           hold_column = "hold",
                           method = "holdout",
                           use_new_features = TRUE)

  ## --------------------------------------------------------------------------------------------
  ## Obtaining and Visualizing Imputed Growth Trajectories.
  ## --------------------------------------------------------------------------------------------
  ## We now obtain the imputed growth trajectories based on the above model fit.
  ## The resulting dataset consist of a single row per subject.
  ## The column "fit" will contain the subject specific predictions of the growth trajectories.
  ## --------------------------------------------------------------------------------------------
  all_preds_holdSL <- predict_all(mfit_holdSL, cpp) %>%
                      convert_to_hbgd(cpp, "sex", "holdSuperLearner")

  ## --------------------------------------------------------------------------------------------
  ## Create trajectory plots and visualize all growth trajectories with trelliscopejs R package.
  ## These plots will include the predictions for all holdout observations,
  ## which can be used to visually inspect the quality of the model fit.
  ## --------------------------------------------------------------------------------------------
  all_preds_holdSL %>%
    hbgd::add_trajectory_plot() %>%
    dplyr::select_("subjid", "panel") %>%
    trelliscopejs::trelliscope(name = "holdSuperLearner")

  ## --------------------------------------------------------------------------------------------
  ## Fitting Growth Trajectories with Cross-Validation SuperLearner
  ## --------------------------------------------------------------------------------------------

  ## --------------------------------------------------------------------------------------------
  ## Fit the growth trajectories using the novel cross-validated discrete SuperLearner (CV SuperLearner).
  ## The CV SuperLearner is enabled by specifying the argument method = "cv".
  ## Note that in this case the model selection is based on ALL validation growth measurements
  ## (all growth measurements on all subjects).
  ## This is in contrast to the previously used holdout SuperLearner (method = "holdout"),
  ## which selects the best model based on only a single holdout growth measurement for each subject.
  ## Finally, note that it is currently not possible to include brokenstick model as a part this
  ## cross-validated SuperLearner ensemble. This is due to the restriction that each included model
  ## must be able to make predictions for new subjects (i.e., validation subjects that were not used
  ## during the fitting of the model).
  ## --------------------------------------------------------------------------------------------
  grid_cvSL <-
    defModel(estimator = "h2o__gbm", family = "gaussian",
             search_criteria = list(strategy = "RandomDiscrete", max_models = 2),
             param_grid = h2o_GBM_hyper,
             seed = 123456) +

    defModel(estimator = "xgboost__gbm", family = "gaussian",
             search_criteria = list(strategy = "RandomDiscrete", max_models = 5),
             param_grid = xgb_GBM_hyper,
             seed = 123456)

  mfit_SLcv <- fit_growth(grid_cvSL,
                         data = cpp,
                         method = "cv",
                         ID = "subjid",
                         t_name = "agedays",
                         x = c("agedays", covars),
                         y = "haz",
                         fold_column = "fold",
                         use_new_features = TRUE)
  ## --------------------------------------------------------------------------------------------
  ## Next, we create a single data set containing the imputed growth trajectories on each subject
  ## --------------------------------------------------------------------------------------------
  all_preds_cvSL <- predict_all(mfit_SLcv, cpp) %>%
                        convert_to_hbgd(cpp, "sex", "cvSuperLearner")

  ## --------------------------------------------------------------------------------------------
  ## This example shows how the imputed subject trajectories can be visualized with trelliscopejs R package.
  ## These plots will include the predictions for all holdout observations,
  ## which can be used to visually inspect the quality of the model fit.
  ## --------------------------------------------------------------------------------------------
  all_preds_cvSL %>%
    hbgd::add_trajectory_plot() %>%
    dplyr::select_("subjid", "panel") %>%
    trelliscopejs::trelliscope(name = "cvSuperLearner")

  ## --------------------------------------------------------------------------------------------
  ## Fitting Growth Trajectories with random holdout CV-SuperLearner
  ## Perform V-fold CV model fitting, but score the models based on a random holdout observation
  ## for each subject in validation fold
  ## --------------------------------------------------------------------------------------------
  mfit_SLholdcv <- fit_growth(grid_cvSL,
                         data = cpp,
                         method = "holdout_cv",
                         ID = "subjid",
                         t_name = "agedays",
                         x = c("agedays", covars),
                         y = "haz",
                         fold_column = "fold",
                         hold_column = "hold",
                         use_new_features = TRUE)
  all_preds_holdcvSL <- predict_all(mfit_SLholdcv, cpp) %>%
                        convert_to_hbgd(cpp, "sex", "holdcvSuperLearner")

  ## --------------------------------------------------------------------------------------------
  ## Combine all modeling predictions into a single data-base:
  ## --------------------------------------------------------------------------------------------
  all_preds_combined <- (all_preds_BS %>% rename(brokenstick = fit)) %>%
                         left_join(
                          all_preds_holdSL %>% rename(holdoutSL = fit)
                          ) %>%
                         left_join(
                          all_preds_cvSL %>% rename(cvSL = fit)
                          )

  ## --------------------------------------------------------------------------------------------
  ## Assessment of Individual Model Performance
  ## --------------------------------------------------------------------------------------------

  ## --------------------------------------------------------------------------------------------
  ## One can use the function make_model_report to generate an in-depth summary of individual models
  ## used by the  SuperLearner.
  ## This report can be generated in either html, pdf or word formats.
  ## In particular, the report will contain the assessment of the performance of each model used
  ## in the SuperLearner ensemble,
  ## including a plotted comparison of the validation mean-squared-error (CV-MSE) for each model.
  ## For example, for the above random holdout SuperLearner, the plot that describes the performance
  ## of each model is shown below (lower CV-MSE implies a better model fit).
  ## Note that this plot also contains the 95% confidence intervals (CIs) around each estimated CV-MSE.
  ## --------------------------------------------------------------------------------------------
  make_model_report(mfit_holdSL, K = 10, data = cpp,
                     title = paste0("Performance of the holdout SuperLearner for Growth Curve Trajectories with CPP Data"),
                     format = "html", keep_md = FALSE,
                     openFile = TRUE)
  plotMSEs(mfit_holdSL, K = 10, interactive = TRUE)

  ## --------------------------------------------------------------------------------------------
  ## Similarly, the model performance report for CV SuperLearner will contain the plot comparing
  ## the cross-validated mean-squared-error (CV-MSE) for each individual model,
  ## along with the corresponding 95% CIs, as shown below.
  ## --------------------------------------------------------------------------------------------
  make_model_report(mfit_SLcv, K = 10, data = cpp,
                     title = paste0("Performance of CV SuperLearner for Growth Curve Trajectories with CPP Data"),
                     format = "html", keep_md = FALSE,
                     openFile = TRUE)
  plotMSEs(mfit_SLcv, K = 10, interactive = TRUE)

  ## --------------------------------------------------------------------------------------------
  ## Finally, model performance report for holdout-CV SuperLearner
  ## --------------------------------------------------------------------------------------------
  make_model_report(mfit_SLholdcv, K = 10, data = cpp,
                     title = paste0("Performance of CV SuperLearner for Growth Curve Trajectories with CPP Data"),
                     format = "html", keep_md = FALSE,
                     openFile = TRUE)
  plotMSEs(mfit_SLholdcv, K = 10, interactive = TRUE)
}
