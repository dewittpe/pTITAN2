#' Permuted TITAN
#'
#' A wrapper function to generate permustated results from
#' \code{\link[TITAN2]{titan}}.
#'
#' @param taxa a list of \code{data.frame}s with the taxa.
#' @param envs a list of \code{data.frame}s with the environmental gradients.
#' @param sid a character vector of length one with the name of the column
#' identifying the station id.
#' @param n the minimum number of occurrences, see \code{\link{occurrences}}.
#' @param permutations number of permutations to generate.
#' @param titan_args a list of named arguments to pass to \code{\link[TITAN2]{titan}}.
#'
#' @seealso \code{\link{permute}}, \code{\link{occurrences}},
#' \code{\link[TITAN2]{titian}}
#'
#' @examples
#'
#' eg <-
#'   ptitan(taxa = list(CD_06_Mall_wID, CN_06_Mall_wID, CN_06_Mall_wID, CN_06_Mall_wID),
#'          envs = list(C_IC_D_06_wID, C_IC_N_06_wID, C_IC_N_06_wID, C_IC_N_06_wID),
#'          sid  = "StationID",
#'          permutations = 20,
#'          titan_args = list(ncpus = 8, numPerm = 50, nBoot = 50, messaging = FALSE)
#'          )
#'
#' eg
#'
#' summary(eg)
#'
#' ggplot2::ggplot(data = eg[eg$metric == "sumz+", ]) +
#'   ggplot2::aes(x = `0.50`, color = treatment, fill = treatment) +
#'   ggplot2::geom_density(alpha = 0.5) +
#'   ggplot2::geom_vline(mapping = ggplot2::aes(xintercept = `0.50`), data = eg[eg$metric == "sumz+" & eg$permutation == 0, ])
#'
#' @export
ptitan <- function(taxa, envs, sid, n = 6L, permutations = 100, titan_args = list(numPerm = 50, nBoot = 50)) {

  stopifnot(length(taxa) == length(envs))

  # run TITAN2::titan on the observed data before permutations
  obs <- vector("list", length = length(taxa))
  names(obs) <- paste0("Treatment", seq_along(obs))

  for (i in seq_along(obs)) {
    this_env <- envs[[i]][which(names(envs[[i]]) != sid)]
    this_txa <- subset(taxa[[i]], select = occurrences(taxa[[i]][which(names(taxa[[i]]) != sid)], n = n)$taxon)
    obs[[i]] <-
      do.call(TITAN2::titan, c(list(  env = this_env, txa = this_txa), titan_args))
  }

  obs_sumz.cp <- lapply(obs, getElement, "sumz.cp")
  obs_sumz.cp <- lapply(obs_sumz.cp, as.data.frame)
  obs_sumz.cp <- do.call(rbind, obs_sumz.cp)
  rownames(obs_sumz.cp) <- paste0("Observed.", rownames(obs_sumz.cp))

  # permute the data and run TITAN2::titan on the permuted data
  perms <- vector("list", length = permutations)
  for (i in seq_along(perms)) {
    p <- permute(taxa = taxa, envs = envs, sid = sid)
    trt_codes <-
      lapply(p, function(x) {subset(x$taxa, select = occurrences(x$taxa, n = n)$taxon)})
    titan_runs <- vector("list", length = length(p))
      for (j in seq_along(titan_runs)) {
        titan_runs[[j]] <-
          do.call(TITAN2::titan, c(list(env = p[[j]][["env"]], txa = trt_codes[[j]]), titan_args))
      }
    perms[[i]] <- lapply(titan_runs, getElement, "sumz.cp")
    perms[[i]] <- lapply(perms[[i]], as.data.frame)
    names(perms[[i]]) <- paste0("Treatment", seq_along(titan_runs))
    perms[[i]] <- do.call(rbind, perms[[i]])
  }

  names(perms) <- paste0("Permutation_", formatC(seq_along(perms), format = "d", flag = "0", width = nchar(permutations)))

  rtn <- rbind(obs_sumz.cp, do.call(rbind, perms))

  idcols <- do.call(rbind, strsplit(rownames(rtn), "\\.") )
  idcols <- as.data.frame(idcols)
  idcols[["permutation"]] <-
            sapply(strsplit(idcols[, 1], "_"),
                   function(x) {
                     x <- x[length(x)]
                     x[x == "Observed"] <- "0"
                     as.integer(x)
                   })

  names(idcols) <- c("set", "treatment", "metric", "permutation")
  rownames(rtn) <- NULL
  rtn <- cbind(rtn, idcols[, -1])
  class(rtn) <- c("pTITAN2_ptitan", class(rtn))
  attr(rtn, "n_treatments") <- length(taxa)
  attr(rtn, "call") <- match.call()
  rtn
}


#' @export
print.pTITAN2_ptitan <- function(x, ...) {
  NextMethod()
}

#' @export
summary.pTITAN2_ptitan <- function(object, ...) {

  d <-
    stats::reshape(
                   object[c("metric", "treatment", "0.50", "permutation")]
                   , direction = "wide"
                   , idvar = c("metric", "permutation")
                   , timevar = "treatment")

  d <- split(d, d$metric)

  d <-
    lapply(d,
           function(x) {
             rtn <- c()
             for (i in seq(3, ncol(x) - 1)) {
               for (j in seq(i + 1, ncol(x))) {
                 a <- abs(x[, i] - x[, j])
                 p <- mean(a[-1] > a[1])
                 names(p) <- paste("treatment", i - 2, "vs", j - 2)
                 rtn <- c(rtn, p)
               }
             }
             rtn
           })

  message("p-values for median change point differing between treatments:")
  do.call(rbind, d)

}



