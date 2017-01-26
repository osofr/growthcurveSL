
## -----------------------------------------------------------------------------
## Class Membership Tests
## -----------------------------------------------------------------------------
is.DataStorageClass <- function(DataStorageClass) "DataStorageClass"%in%class(DataStorageClass)
is.PredictionModel <- function(PredictionModel) "PredictionModel"%in%class(PredictionModel)
is.PredictionStack <- function(PredictionStack) "PredictionStack"%in%class(PredictionStack)
is.ModelStack <- function(obj) "ModelStack" %in% class(obj)


# ---------------------------------------------------------------------------------------
# Import data, define nodes (columns), define dummies for factor columns and define input data R6 object
#
# @param data Input dataset, can be a \code{data.frame} or a \code{data.table}.
# @param ID A character string name of the column that contains the unique subject identifiers.
# @param t_name A character string name of the column with integer-valued measurement time-points (in days, weeks, months, etc).
# @param covars Names of predictors (covariates) in the data.
# @param OUTCOME Character name of the column containing outcomes.
# @param verbose Set to \code{TRUE} to print messages on status and information to the console. Turn this on by default using \code{options(growthcurveSL.verbose=TRUE)}.
# @return An R6 object that contains the input data. This can be passed as an argument to \code{get_fit} function.
# @example tests/examples/1_growthcurveSL_example.R
# @export
importData <- function(data, ID = "Subject_ID", t_name = "time_period", covars, OUTCOME = "Y", verbose = getOption("growthcurveSL.verbose")) {
  gvars$verbose <- verbose

  if (missing(covars)) { # define time-varing covars (L) as everything else in data besides these vars
    covars <- setdiff(colnames(data), c(ID, OUTCOME))
  }

  nodes <- list(Lnodes = covars, Ynode = OUTCOME, IDnode = ID, tnode = t_name)
  OData <- DataStorageClass$new(Odata = data, nodes = nodes)

  for (Lnode in nodes$Lnodes) CheckVarNameExists(OData$dat.sVar, Lnode)
  return(OData)
}

## -----------------------------------------------------------------------------
## General utilities / Global Vars
## -----------------------------------------------------------------------------
`%+%` <- function(a, b) paste0(a, b)
is.integerish <- function (x) is.integer(x) || (is.numeric(x) && all(x == as.integer(x)))

# Return the left hand side variable of formula f as a character
LhsVars <- function(f) {
  f <- as.formula(f)
  return(as.character(f[[2]]))
}
# Return the right hand side variables of formula f as a character vector
RhsVars <- function(f) {
  f <- as.formula(f)
  return(all.vars(f[[3]]))
}

checkpkgs <- function(pkgs) {
  for (pkg in pkgs) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      stop(pkg %+% " package needed for this function to work. Please install it.", call. = FALSE)
    }
  }
}

#if warning is in ignoreWarningList, ignore it; otherwise post it as usual
SuppressGivenWarnings <- function(expr, warningsToIgnore) {
  h <- function (w) {
    if (w$message %in% warningsToIgnore) invokeRestart( "muffleWarning" )
  }
  withCallingHandlers(expr, warning = h )
}

GetWarningsToSuppress <- function(update.step=FALSE) {
  warnings.to.suppress <- c("glm.fit: fitted probabilities numerically 0 or 1 occurred",
                            "prediction from a rank-deficient fit may be misleading",
                            "non-integer #successes in a binomial glm!",
                            "the matrix is either rank-deficient or indefinite")
  if (update.step) {
    warnings.to.suppress <- c(warnings.to.suppress, "glm.fit: algorithm did not converge")
  }
  return(warnings.to.suppress)
}

# returns NULL if no factors exist, otherwise return the name of the factor variable(s)
CheckExistFactors <- function(data) {
  testvec <- unlist(lapply(data, is.factor))
  if (any(testvec)) {
    return(names(data)[which(testvec)])
  } else {
    return(NULL)
  }
}

# throw exception if 1) varname doesn't exist; 2) more than one varname is matched
CheckVarNameExists <- function(data, varname) {
  idvar <- names(data) %in% varname
  if (sum(idvar) < 1) stop("variable name " %+% varname %+% " not found in data input")
  if (sum(idvar) > 1) stop("more than one column in the input data has been matched to name "
                            %+% varname %+% ". Consider renaming some of the columns: " %+%
                            paste0(names(data)[idvar], collapse=","))
  return(invisible(NULL))
}

