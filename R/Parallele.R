
autoStopCluster <- function(cl) {
  try(detectCores())
  stopifnot(inherits(cl, "cluster"))
  env <- new.env()
  env$cluster <- cl
  attr(cl, "gcMe") <- env
  reg.finalizer(env, function(e) {
    #message("Finalizing cluster ...")
    #message(capture.output(print(e$cluster)))
    try(parallel::stopCluster(e$cluster), silent = FALSE)
    
  })
  message("Finalizing cluster ... done")
  cl
}


automakeCluster <- function(performance=1) {
  
  ncores=try(detectCores())
  
  cl=parallel::makeCluster(round(ncores*performance))
  message("Finalizing cluster ... done")
  env <- new.env()
  env$cluster <- cl
  attr(cl, "gcMe") <- env


  cl
}


do_Parallel <- function(ncores = parallel::detectCores() - 1,pop, fun,  ...) {
  require(foreach)
  cl <- parallel::makeCluster(ncores)
  on.exit(parallel::stopCluster(cl), add = TRUE)
  doParallel::registerDoParallel(cl)
  foreach(i = pop) %dopar% fun(i, ...)
}
