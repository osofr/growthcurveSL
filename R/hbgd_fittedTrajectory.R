who_centile2value <- function(x, p = 50, x_var = "agedays", y_var = "htcm",
  sex = "Female", data = NULL) {

  if (!is.null(data)) {
    x <- hbgd:::v_eval(substitute(x), try(x, silent = TRUE), data)
    p <- hbgd:::v_eval(substitute(p), try(p, silent = TRUE), data)
    x_var <- hbgd:::v_eval(substitute(x_var), try(x_var, silent = TRUE), data)
    y_var <- hbgd:::v_eval(substitute(y_var), try(y_var, silent = TRUE), data)
    sex <- hbgd:::v_eval(substitute(sex), try(sex, silent = TRUE), data)
  }

  dat <- data.frame(
    x = x,
    p = p,
    x_var = x_var,
    y_var = y_var,
    sex = sex,
    stringsAsFactors = FALSE
  )

  if (! all(unique(dat$sex) %in% c("Male", "Female")))
    stop("sex must be 'Male' or 'Female'")

  # since coefficients are available only by pair/sex
  # we need to call this for each unique combination
  centile2value_single_pars <- function(x, y, x_var, y_var, sex) {
    pair <- paste(y_var, x_var, sep = "_")
    hbgd:::check_pair(pair)

    coefs <- hbgd::who_coefs[[pair]][[sex]]$data

    # OS added: safety check for missing input
    if (all(is.na(x))) return(NA)

    # subset to neighborhood surrounding input
    idx <- hbgd:::get_coef_idx(x, coefs$x)

    coefs <- coefs[idx, , drop = FALSE] # nolint
    if (nrow(coefs) == 1) {
      coefs <- data.frame(y = y, coefs, row.names = NULL)
    } else {
      coefs <- data.frame(
        x = x,
        y = y,
        l = approx(coefs$x, coefs$l, x)$y,
        m = approx(coefs$x, coefs$m, x)$y,
        s = approx(coefs$x, coefs$s, x)$y)
    }

    with(coefs, m * ((1 + qnorm(y / 100) * l * s)^(1 / l))) # nolint
  }

  dat <- dat %>%
    dplyr::group_by(x_var, y_var, sex) %>%
    dplyr::mutate(res = centile2value_single_pars(x, p, x_var[1], y_var[1], sex[1]))

  dat$res
}

who_zscore2value <- function(x, z = 0, y_var = "htcm", x_var = "agedays",
  sex = "Female", data = NULL) {

  if (!is.null(data)) {
    x <- hbgd:::v_eval(substitute(x), try(x, silent = TRUE), data)
    z <- hbgd:::v_eval(substitute(z), try(z, silent = TRUE), data)
    x_var <- hbgd:::v_eval(substitute(x_var), try(x_var, silent = TRUE), data)
    y_var <- hbgd:::v_eval(substitute(y_var), try(y_var, silent = TRUE), data)
    sex <- hbgd:::v_eval(substitute(sex), try(sex, silent = TRUE), data)
  }

  who_centile2value(p = 100 * pnorm(z), x = x, y_var = y_var, x_var = x_var, sex = sex)
}

who_zscore2htcm <- function(agedays, z = 0, sex = "Female") {
  who_zscore2value(agedays, z, x_var = "agedays", y_var = "htcm", sex = sex)
}
who_zscore2wtkg <- function(agedays, z = 0, sex = "Female") {
  who_zscore2value(agedays, z, x_var = "agedays", y_var = "wtkg", sex = sex)
}
who_zscore2bmi <- function(agedays, z = 0, sex = "Female") {
  who_zscore2value(agedays, z, x_var = "agedays", y_var = "bmi", sex = sex)
}
who_zscore2hcircm <- function(agedays, z = 0, sex = "Female") {
  who_zscore2value(agedays, z, x_var = "agedays", y_var = "hcircm", sex = sex)
}
who_zscore2muaccm <- function(agedays, z = 0, sex = "Female") {
  who_zscore2value(agedays, z, x_var = "agedays", y_var = "muaccm", sex = sex)
}
who_zscore2ssftmm <- function(agedays, z = 0, sex = "Female") {
  who_zscore2value(agedays, z, x_var = "agedays", y_var = "ssftmm", sex = sex)
}
who_zscore2tsftmm <- function(agedays, z = 0, sex = "Female") {
  who_zscore2value(agedays, z, x_var = "agedays", y_var = "tsftmm", sex = sex)
}

who_htcm2zscore <- function(agedays, htcm, sex = "Female") {
  hbgd::who_value2zscore(agedays, htcm, x_var = "agedays", y_var = "htcm", sex = sex)
}
who_wtkg2zscore <- function(agedays, wtkg, sex = "Female") {
  hbgd::who_value2zscore(agedays, wtkg, x_var = "agedays", y_var = "wtkg", sex = sex)
}
who_bmi2zscore <- function(agedays, bmi, sex = "Female") {
  hbgd::who_value2zscore(agedays, bmi, x_var = "agedays", y_var = "bmi", sex = sex)
}
who_hcircm2zscore <- function(agedays, hcircm, sex = "Female") {
  hbgd::who_value2zscore(agedays, hcircm, x_var = "agedays", y_var = "hcircm", sex = sex)
}
who_muaccm2zscore <- function(agedays, muaccm, sex = "Female") {
  hbgd::who_value2zscore(agedays, muaccm, x_var = "agedays", y_var = "muaccm", sex = sex)
}
who_ssftmm2zscore <- function(agedays, ssftmm, sex = "Female") {
  hbgd::who_value2zscore(agedays, ssftmm, x_var = "agedays", y_var = "ssftmm", sex = sex)
}
who_tsftmm2zscore <- function(agedays, tsftmm, sex = "Female") {
  hbgd::who_value2zscore(agedays, tsftmm, x_var = "agedays", y_var = "tsftmm", sex = sex)
}


#' Wrapper to convert fit object into \code{hbgd}-compatible object
#' @param fit_dat_all Dataset produced by calling \code{\link{predict_all}} function.
#' @param data The input data used for training.
#' @param ID The column name in input data containing unique subject identifiers.
#' @param sexvar The column name in the input data containing subject gender coded as "Male" / "Female".
#' @param method The name of the method that produces the growth curve fits.
#' @return ...
#' @export
convert_to_hbgd <- function(fit_dat_all, data, sexvar = "sex", method = "default") {
    nodes <- attr(fit_dat_all, "nodes")
    ID <- nodes$IDnode
    Ynode <- nodes$Ynode
    if (ID != "subjid") stop("The subject identifier must be named 'subjid'")
    # checkpoints = c(1, hbgd::months2days(1:24)),

    if (Ynode %in% "haz") {
      xy_pair_name = c("agedays","htcm")
      fun_y_to_raw = who_zscore2htcm
      fun_y_to_z = function(x, y, ...) return(y)

    } else if (Ynode %in% "htcm") {
      xy_pair_name = c("agedays","htcm")
      fun_y_to_raw = function(x, y, ...) return(y)
      fun_y_to_z = who_htcm2zscore

    } else if (Ynode %in% "waz") {
      xy_pair_name = c("agedays","wtkg")
      fun_y_to_raw = who_zscore2wtkg
      fun_y_to_z = function(x, y, ...) return(y)

    } else if (Ynode %in% "wtkg") {
      xy_pair_name = c("agedays","wtkg")
      fun_y_to_raw = function(x, y, ...) return(y)
      fun_y_to_z = who_wtkg2zscore

    } else if (Ynode %in% "baz") {
      xy_pair_name = c("agedays","bmi")
      fun_y_to_raw = who_zscore2bmi
      fun_y_to_z = function(x, y, ...) return(y)

    } else if (Ynode %in% "bmi") {
      xy_pair_name = c("agedays","bmi")
      fun_y_to_raw = function(x, y, ...) return(y)
      fun_y_to_z = who_bmi2zscore

    } else if (Ynode %in% "hcaz") {
      xy_pair_name = c("agedays","hcircm")
      fun_y_to_raw = who_zscore2hcircm
      fun_y_to_z = function(x, y, ...) return(y)

    } else if (Ynode %in% "hcircm") {
      xy_pair_name = c("agedays","hcircm")
      fun_y_to_raw = function(x, y, ...) return(y)
      fun_y_to_z = who_hcircm2zscore

    } else {
      stop("unrecognized outcome: " %+% Ynode)
    }

    fit_dat_hbgd <- data %>%
      dplyr::distinct_(ID, sexvar) %>%
      dplyr::rename_("sex" = sexvar)

    if (is.integerish(fit_dat_hbgd[["sex"]]))
      fit_dat_hbgd[["sex"]] <- as.integer(fit_dat_hbgd[["sex"]]) %>% dplyr::recode(`0` = "Female", `1` = "Male")



    fit_dat_hbgd <-
      fit_dat_hbgd %>%
      dplyr::left_join(fit_dat_all) %>%
      tibble::as_tibble()

    # MSE <-
    #   fit_dat_hbgd %>%
    #   tidyr::unnest(fit) %>%
    #   dplyr::select_(ID, "MSE") %>%
    #   tidyr::nest(-subjid, .key = MSE)

    xy <-
      fit_dat_hbgd %>%
      tidyr::unnest(fit) %>%
      dplyr::select_(ID, sexvar, "fit") %>%
      tidyr::unnest(fit) %>%
      dplyr::rename_("Y_modeled" = "y") %>%
      dplyr::mutate(y = fun_y_to_raw(x, Y_modeled, sex = sex)) %>%
      dplyr::mutate(z = fun_y_to_z(x, Y_modeled, sex = sex)) %>%
      dplyr::mutate(yfit = fun_y_to_raw(x, preds, sex = sex)) %>%
      dplyr::mutate(zfit = fun_y_to_z(x, preds, sex = sex)) %>%
      dplyr::select_(ID, "x", "y", "z", "yfit", "zfit") %>%
      tidyr::nest(-subjid, .key = xy)

    ## the residuals of the fit
    resid <-
      xy %>%
      tidyr::unnest(xy) %>%
      dplyr::mutate(resid = y - yfit) %>%
      dplyr::select(subjid, resid) %>%
      tidyr::nest(-subjid, .key = resid)

    fitgrid <-
      fit_dat_hbgd %>%
      tidyr::unnest(fit) %>%
      dplyr::select_(ID, sexvar, "fitgrid") %>%
      tidyr::unnest(fitgrid) %>%
      dplyr::mutate(y = fun_y_to_raw(x, preds, sex = sex)) %>%
      dplyr::mutate(z = fun_y_to_z(x, preds, sex = sex)) %>%
      dplyr::select_(ID, "x", "y", "z") %>%
      dplyr::mutate(dy = NA, dz = NA) %>%
      tidyr::nest(-subjid, .key = fitgrid)
    ## add derivative on original and z-score scale
    # if (!is.null(fitgrid$y)) fitgrid$dy <- hbgd::grid_deriv(fitgrid$x, fitgrid$y)
    # if (!is.null(fitgrid$z)) fitgrid$dz <- hbgd::grid_deriv(fitgrid$x, fitgrid$z)
    # res[["fitgrid"]] <- fitgrid

    holdout <-
      fit_dat_hbgd %>%
      tidyr::unnest(fit) %>%
      dplyr::select_(ID, sexvar, "holdout") %>%
      tidyr::unnest(holdout) %>%
      dplyr::mutate(y = fun_y_to_raw(x, preds, sex = sex)) %>%
      dplyr::mutate(z = fun_y_to_z(x, preds, sex = sex)) %>%
      dplyr::select_(ID, "x", "y", "z") %>%
      tidyr::nest(-subjid, .key = holdout)

    checkpoint <-
      fit_dat_hbgd %>%
      tidyr::unnest(fit) %>%
      dplyr::select_(ID, sexvar, "checkpoint") %>%
      tidyr::unnest(checkpoint) %>%
      dplyr::mutate(y = fun_y_to_raw(x, preds, sex = sex)) %>%
      dplyr::mutate(z = fun_y_to_z(x, preds, sex = sex)) %>%
      dplyr::select_(ID, "x", "y", "z") %>%
      dplyr::mutate(zcat = NA) %>%
      tidyr::nest(-subjid, .key = checkpoint)

      # for (i in 1:nrow(checkpoint))
      #   print(all.equal(checkpoint[["checkpoint"]][[i]], fit_dat_hbgd_slow[["fit"]][[i]][["checkpoint"]]))
      # fit_dat_hbgd_slow[["fit"]][[1]][["xy"]]

    # browser()

    fit_dat_hbgd_fast <-
      fit_dat_hbgd %>%
      dplyr::select(subjid, sex) %>%
      dplyr::left_join(fit_dat_all %>% tidyr::unnest(fit) %>% dplyr::select_(ID, "MSE")) %>%
      dplyr::left_join(xy) %>%
      dplyr::left_join(resid) %>%
      dplyr::left_join(fitgrid) %>%
      dplyr::left_join(holdout) %>%
      dplyr::left_join(checkpoint) %>%
      dplyr::mutate(zcat = NA, pars = NA, x_var = xy_pair_name[1], y_var = xy_pair_name[2], method = method) %>%
      tidyr::nest(-subjid, .key= fit) %>%
      dplyr::mutate(fit = purrr::map(fit, function(res) {
        # browser();
        res <- list(
          xy = res[["xy"]][[1]],
          resid = res[["resid"]][[1]],
          fitgrid = res[["fitgrid"]][[1]],
          holdout = res[["holdout"]][[1]],
          checkpoint = res[["checkpoint"]][[1]],
          zcat = res[["zcat"]],
          MSE = res[["MSE"]],
          pars = res[["pars"]], x_var = res[["x_var"]], y_var = res[["y_var"]], method = res[["method"]], sex = res[["sex"]]
          )
        class(res) <- "fittedTrajectory";
        res}))

    # fit_dat_hbgd_fast[2, "fit"][["fit"]]

    fit_dat_hbgd_fast

    # fit_dat_hbgd[1:10, ][["fit"]]
    # fit_dat_hbgd_fast[["fit"]]
    # fit_dat_hbgd_fast[1, ][["fit"]]

    # fit_dat_hbgd_slow <-
    #   fit_dat_hbgd %>%
    #   plyr::mutate(fit = purrr::map2(fit, sex,
    #       ~ add_cogs_persubj(fit_dat = .x, sex = .y, method = method)))
    # fit_dat_hbgd_slow
    # fit_dat_hbgd_slow[1, ][["fit"]]

  # browser()
  # checkpoint[1, "checkpoint"][[1]][[1]]
  # unlist(fit_dat_hbgd_fast[1, "holdout"][[1]])
  # fit_dat_hbgd_fast[1, "fit"][[1]]
  # fit_dat_hbgd_slow[1, "fit"][[1]]
}

#' Create a fittedTrajectory object with fits for a single subject
#'
#' Creates objects of class \code{fittedTrajectory} with growth curve predictions, one object for subject ID in the data.
#' This can be used for further growth curve analysis with hbgd R package.
#' Requires the the dataset containing the predicted growth curves for all subjects \code{grid_fits_dat}.
#' @param fit_dat Subject specific data to be added to the final object.
#' This must also contain the observed outcome on each subect under a nested dataset "fit", in column named "y".
#' Similarly, the model predictions must be under the nested dataset "fit", in column named "preds".
#' Additionally, the holdout and grid predictions can be stored in nested datasets "holdout" and "fitgrid", respectively.
#' @param sex The gender (as character "Male" / "Female") for this subject.
#' @param method A name of the modeling approach for future model comparison.
#' @param xy_pair_name ... Same as in hbgd package...
#' @param fun_y_to_raw The function to convert the scale of the observed outcome and predictions into the raw scale (height / weight).
#' If the fits are already on the raw scale then specify this argument as \code{function(x, y, ...) return(y)}.
#' @param fun_y_to_z The function to convert the scale of the observed outcome and predictions into the z-scale.
#' If the fits are already on the z-scale then specify this argument as \code{function(x, y, ...) return(y)}.
#' @return ...
#' @export
add_cogs_persubj <- function(fit_dat, sex,
                             method = "default",
                             xy_pair_name = c("agedays","htcm"),
                             # checkpoints = c(1, hbgd::months2days(1:24)),
                             fun_y_to_raw = who_zscore2htcm,
                             # fun_y_to_raw = hbgd::who_zscore2htcm,
                             fun_y_to_z = function(x, y, ...) return(y)) {

    res <- list()
    res[["xy"]] <- fit_dat %>%
                tidyr::unnest(fit) %>%
                dplyr::rename_("Y_modeled" = "y") %>%
                dplyr::mutate(y = fun_y_to_raw(x, Y_modeled, sex = sex)) %>%
                dplyr::mutate(z = fun_y_to_z(x, Y_modeled, sex = sex)) %>%
                dplyr::mutate(yfit = fun_y_to_raw(x, preds, sex = sex)) %>%
                dplyr::mutate(zfit = fun_y_to_z(x, preds, sex = sex)) %>%
                dplyr::select_("x", "y", "z", "yfit", "zfit")

    ## the residuals of the fit
    res[["resid"]] <- res[["xy"]] %>% dplyr::transmute(resid = y - yfit)

    ## a data frame with the model fits applied across an equally-spaced grid of “x” points with columns (x,y,z), z is optional z-scores of y
    ## transform to raw-scale if fitted on Z scale, otherwise keep as is
    fitgrid <- fit_dat %>%
                tidyr::unnest(fitgrid) %>%
                dplyr::mutate(y = fun_y_to_raw(x, preds, sex = sex)) %>%
                dplyr::mutate(z = fun_y_to_z(x, preds, sex = sex)) %>%
                dplyr::select_("x", "y", "z")
    ## add derivative on original and z-score scale
    if (!is.null(fitgrid$y)) fitgrid$dy <- hbgd::grid_deriv(fitgrid$x, fitgrid$y)
    if (!is.null(fitgrid$z)) fitgrid$dz <- hbgd::grid_deriv(fitgrid$x, fitgrid$z)
    res[["fitgrid"]] <- fitgrid

    res[["holdout"]] <- fit_dat %>%
                tidyr::unnest(holdout) %>%
                dplyr::mutate(y = fun_y_to_raw(x, preds, sex = sex)) %>%
                dplyr::mutate(z = fun_y_to_z(x, preds, sex = sex)) %>%
                dplyr::select_("x", "y", "z")

    res[["checkpoint"]] <- fit_dat %>%
                tidyr::unnest(checkpoint) %>%
                dplyr::mutate(y = fun_y_to_raw(x, preds, sex = sex)) %>%
                dplyr::mutate(z = fun_y_to_z(x, preds, sex = sex)) %>%
                dplyr::select_("x", "y", "z") %>%
                dplyr::mutate(zcat = NA)
    # res[["checkpoint"]] <- data.frame(x = checkpoints, y = NA, z = NA, zcat = NA)
    # res[["checkpoint"]] <- checkpoint

    # a z-score categorization (e.g. <-2 or >-2) for where the subject’s growth falls in the z-score scale at each checkpoint
    res[["zcat"]] <- NA
    # a named list of fitting parameters supplied to the particular fitting method
    res[["pars"]] <- NA
    # all columns of the data argument for this subject (preserved for later analyses)

    # the subject’s sex
    res[["sex"]] <- sex

    # # the names of the variables used for “x” and “y” in the trajectory fitting
    res[["x_var"]] <- xy_pair_name[1]
    res[["y_var"]] <- xy_pair_name[2]

    # # the name of the fitting method used
    res[["method"]] <- method

    class(res) <- "fittedTrajectory"
    res
}
