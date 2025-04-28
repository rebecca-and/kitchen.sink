# Used in process_rds_data -> function from RDS Analyst
wtd.contingency.tables <- function (row.vars, col.vars, stratum.var, weights = NULL, data = NULL,
                                    missing.include = FALSE, subset = NULL, norm.weights = FALSE)
{
  subset <- eval(substitute(subset), data, parent.frame())
  if (!is.null(subset))
    stop("subset not yet implemented")
  arguments <- as.list(match.call()[-1])
  if (missing(row.vars) || missing(col.vars))
    stop("Please specify the row variables (row.vars), and column variables (col.vars)")
  row.vars <- eval(substitute(row.vars), data, parent.frame())
  col.vars <- eval(substitute(col.vars), data, parent.frame())
  if (length(dim(row.vars)) < 1.5) {
    row.vars <- as.data.frame(row.vars)
    names(row.vars) <- as.character(arguments$row.vars)
  }
  if (length(dim(col.vars)) < 1.5) {
    col.vars <- as.data.frame(col.vars)
    names(col.vars) <- as.character(arguments$col.vars)
  }
  if (!missing(stratum.var))
    stratum.var <- eval(substitute(stratum.var), data, parent.frame())
  else stratum.var <- NULL
  vector.x <- FALSE
  num.row.vars <- dim(row.vars)[2]
  num.col.vars <- dim(col.vars)[2]
  if (!is.null(weights)) {
    weights <- eval(substitute(weights), data, parent.frame())
  }
  else {
    weights <- rep(1, nrow(data))
  }
  single.table <- function(dat, dnn, weights) {
    x <- dat[[1]]
    y <- dat[[2]]
    if (is.null(stratum.var))
      stratum.var <- rep("No Strata", length(x))
    if (length(x) != length(y))
      stop("all row.vars and col.vars must have the same length")
    if (missing.include) {
      x <- factor(x, exclude = c())
      y <- factor(y, exclude = c())
      strata <- factor(stratum.var, exclude = c())
    }
    else {
      x <- factor(x)
      y <- factor(y)
      strata <- factor(stratum.var)
    }
    lev <- levels(strata)
    table.list <- list()
    for (level in lev) {
      temp.x <- x[strata == level]
      temp.y <- y[strata == level]
      temp.w <- weights[strata == level]
      if (norm.weights) {
        temp.n <- sum(!is.na(temp.x) & !is.na(temp.y))
        temp.w <- temp.n * temp.w/sum(temp.w)
      }
      t <- xtabs(temp.w ~ temp.x + temp.y)
      if (sum(t) == 0) {
        stop(paste("No non-missing values for", paste(dnn,
                                                      collapse = " and ")))
      }
      names(dimnames(t)) <- dnn
      CPR <- prop.table(t, 1)
      CPC <- prop.table(t, 2)
      CPT <- prop.table(t)
      RS <- round(rowSums(t), 2)
      CS <- round(colSums(t), 2)
      t <- t[RS > 0, CS > 0, drop = FALSE]
      RS <- round(rowSums(t), 2)
      CS <- round(colSums(t), 2)
      CST <- try(suppressWarnings(chisq.test(t, correct = FALSE)),
                 silent = TRUE)
      if (class(CST) != "htest")
        CST <- list(expected = t * NA)
      GT <- sum(t)
      if (length(dim(x) == 2))
        TotalN <- GT
      else TotalN <- length(temp.x)
      table.list[[level]] <- list(table = round(t, 1),
                                  row.sums = RS, col.sums = CS, total = GT, row.prop = CPR,
                                  col.prop = CPC, total.prop = CPT, expected = CST$expected)
      class(table.list[[level]]) <- "single.table"
    }
    class(table.list) <- "contin.table"
    table.list
  }
  result <- list()
  count <- 1
  for (i in 1:num.row.vars) {
    for (j in 1:num.col.vars) {
      result[[paste(names(row.vars)[i], "by", names(col.vars)[j])]] <- single.table(data.frame(row.vars[,
                                                                                                        i], col.vars[, j]), dnn = c(names(row.vars)[i],
                                                                                                                                    names(col.vars)[j]), weights = weights)
      count <- count + 1
    }
  }
  class(result) <- "contingency.tables"
  attr(result, "rowNames") = paste(as.character(arguments$row.vars),
                                   collapse = ", ")
  attr(result, "colNames") = paste(as.character(arguments$col.vars),
                                   collapse = ", ")
  if (!missing(stratum.var))
    attr(result, "strata.name") = as.character(arguments$stratum.var)
  result
}


# Used in process_rds_data -> function from RDS Analyst
d <- function (..., row.names = NULL, check.rows = FALSE, check.names = FALSE,
               stringsAsFactors = FALSE)
{
  data.frame(..., row.names = row.names, check.rows = check.rows,
             check.names = check.names, stringsAsFactors = stringsAsFactors)
}


# Obtain RDS weighted estimates of age distributions + HIV prev by age
process_rds_data <- function(df, n_bootstrap) {
  library(RDS)

  df <- df[, !(grepl("^coupon", names(df)) & sapply(df, function(x) all(is.na(x))))]

  pop_size <- unique(df$pop_size)
  pse_method <- unique(df$pse_meth)
  survey_id <- unique(df$survey_id)
  survey_city <- unique(df$survey_city)
  coupons <- colnames(df)[grep("^coupon", colnames(df))]
  num_coupons <- length(coupons)

  df$recruiter.id <- rid.from.coupons(df,
                                      subject.id = 'rowid',
                                      subject.coupon = 'own_coupon',
                                      coupon.variables = coupons)

  df <- as.rds.data.frame(df, id = 'rowid',
                          recruiter.id = 'recruiter.id',
                          network.size = 'network_size',
                          population.size = pop_size,
                          max.coupons = num_coupons,
                          notes = paste0("PSE Method:", pse_method),
                          time = NULL)

  df$seed <- get.seed.id(df)
  df$wave <- get.wave(df)


  unstratified_weighted_estimates <- RDS.bootstrap.intervals(df,
                                                             outcome.variable = c("age", "hiv"),
                                                             weight.type = "Gile's SS",
                                                             uncertainty = "Gile",
                                                             confidence.level = 0.95,
                                                             number.of.bootstrap.samples = n_bootstrap,
                                                             to.factor = TRUE,
                                                             N = pop_size,
                                                             number.ss.samples.per.iteration = 1000)


  stratified_weighted_estimates <- wtd.contingency.tables(
    row.vars = d(age),
    col.vars = d(hiv),
    data = df,
    weights = compute.weights(df, 'RDS-II'),
    norm.weights = TRUE
  )

  stratified_wtd_counts <- data.frame(
    stratified_weighted_estimates[["age by hiv"]][["No Strata"]][["table"]]
  ) %>%
    dplyr::mutate(survey_id = survey_id, survey_city = survey_city)

  prev_by_age <- data.frame(
    stratified_weighted_estimates[["age by hiv"]][["No Strata"]][["row.prop"]]
  ) %>%
    dplyr::mutate(survey_id = survey_id, survey_city = survey_city, pse_method = pse_method)

  return(list(
    processed_df = df,
    unstratified_estimates = unstratified_weighted_estimates,
    stratified_counts = stratified_wtd_counts,
    prevalence_by_age = prev_by_age
  ))
}
