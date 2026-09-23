library(SuperSurv)
data("metabric", package = "SuperSurv")
dat <- metabric[1:70, ]
X <- dat[, grep("^x", names(dat))[1:3], drop = FALSE]
learners <- list(surv.coxph = list(), surv.weibull = list(), surv.km = list())
if (requireNamespace("randomForestSRC", quietly = TRUE)) {
  learners$surv.rfsrc <- list(ntree = 10, nodesize = 5)
}
if (requireNamespace("ranger", quietly = TRUE)) {
  learners$surv.ranger <- list(num.trees = 10, min.node.size = 5)
}
for (learner in names(learners)) {
  for (n in c(1L, 3L)) {
    for (times in list(40, c(20, 40, 60))) {
      set.seed(6384)
      fit <- do.call(learner, c(list(
        time = dat$duration, event = dat$event, X = X,
        newdata = X[seq_len(n), , drop = FALSE], new.times = times,
        obsWeights = rep(1, nrow(X)), id = NULL
      ), learners[[learner]]))
      prediction <- predict(fit$fit, newdata = X[seq_len(n), , drop = FALSE],
                            new.times = times)
      stopifnot(identical(dim(fit$pred), c(n, length(times))),
                identical(dim(prediction), c(n, length(times))),
                max(abs(prediction - fit$pred)) < 1e-10,
                all(prediction >= 0 & prediction <= 1))
    }
  }
}
