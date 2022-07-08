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
#' @param nBoot number of bootstraps, see \code{\link[TITAN2]{titan}}.
#' @param permutations number of permutations to generate.
#'
#' @seealso \code{\link{permute}}, \code{\link{occurrences}},
#' \code{\link[TITAN2]{titian}}
#'
#' @export
ptitan <- function(taxa, envs, sid, n = 6L, nBoot = 5L, permutations) {
  p <- permute(taxa = taxa, envs = envs, sid = "StationID")
  trt_codes <-
    lapply(p, function(x) {subset(x$taxa, select = occurrences(x$taxa, n = n)$taxon)})

  i <- 1

  out <-
    # foreach::`%dopar%`(
    # foreach::foreach(i = seq_along(p))
    # ,
    # {
          TITAN2::titan(env = p[[i]][["env"]][[1]],
                        txa = trt_codes[[i]],
                        boot = FALSE,
                        nBoot = 100,
                        ncpus = 4)
  # }
  # )
  out
}


#library(doParallel)
#workers = makeCluster(4, type = "SOCK")
#registerDoParallel(workers)

out <- 
  ptitan(taxa = list(CD_06_Mall_wID, CN_06_Mall_wID, CD_06_Mall_wID),
         envs = list(C_IC_D_06_wID, C_IC_N_06_wID, C_IC_D_06_wID),
         sid = "StationID"
         )

str(out, max.level = 1)
out$sumz.cp
names(out)

out$metricArray[,,1]


head(out$ivz)

getivz(spp = out$sppmax)


getivz(out[[1]])

out[[1]]$sumz.cp
nrow(out[[1]]$env)

cpsumz_neg <- 


str(out[[1]], max.level = 2)





# library(pTITAN2)
# p <- permute(taxa = list(CD_06_Mall_wID, CN_06_Mall_wID, CD_06_Mall_wID),
#              envs = list(C_IC_D_06_wID, C_IC_N_06_wID, C_IC_D_06_wID),
#              sid = "StationID")
# 
# 
# 
# 
# #?foreach
# 
# #library(doMC)
# #registerDoMC(cores = 4)
# 
# 
# 
# 
# 
# Map(f = function(env, taxa, nBoot) {
#       TITAN2::titan(env = env, txa = taxa, boot = TRUE, nBoot = nBoot)
#              }, 
#       env = lapply(p, function(x) {x$env[[1]]}),
#       taxa = trt_codes,
#       nBoot = 5)
# 
#        
# 
# 
#        p[[1]][["env"]][[1]]
#        trt_codes[[1]]
# 
