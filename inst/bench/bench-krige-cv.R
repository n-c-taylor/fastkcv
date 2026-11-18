# Benchmark: kcv() vs gstat::krige.cv(), global kriging.
# Run from the package root with: Rscript inst/bench/bench-krige-cv.R
# Not run on CRAN (inst/bench is in .Rbuildignore).

devtools::load_all(quiet = TRUE)
library(gstat)

set.seed(2026)

# meuse, n = 155 -------------------------------------------------------------
data(meuse, package = "sp")
v_meuse <- vgm(0.5, "Sph", 900, 0.1)
folds10_meuse <- sample(rep(1:10, length.out = nrow(meuse)))

# simulated, n = 1000 --------------------------------------------------------
n_sim <- 1000
v_sim <- vgm(1, "Exp", 200, 0.1)
sim <- data.frame(x = runif(n_sim, 0, 1000), y = runif(n_sim, 0, 1000))
C_sim <- fastkcv:::covariance_matrix(v_sim, as.matrix(sim))
sim$z <- drop(crossprod(chol(C_sim), rnorm(n_sim))) + 5
folds10_sim <- sample(rep(1:10, length.out = n_sim))

run <- function(label, n, expr_kcv, expr_gstat, iter_gstat) {
  b1 <- bench::mark(kcv = eval(expr_kcv), iterations = 10, check = FALSE,
    filter_gc = FALSE)
  b2 <- bench::mark(krige.cv = eval(expr_gstat), iterations = iter_gstat,
    check = FALSE, filter_gc = FALSE)
  t1 <- as.numeric(b1$median)
  t2 <- as.numeric(b2$median)
  data.frame(
    case = label, n = n,
    kcv_median_s = signif(t1, 3), krige.cv_median_s = signif(t2, 3),
    speedup = round(t2 / t1, 1),
    kcv_mem = format(b1$mem_alloc), krige.cv_mem = format(b2$mem_alloc)
  )
}

res <- rbind(
  run("meuse LOO", 155,
    quote(kcv(log(zinc) ~ 1, ~ x + y, meuse, model = v_meuse)),
    quote(krige.cv(log(zinc) ~ 1, ~ x + y, meuse, model = v_meuse, verbose = FALSE)),
    10),
  run("meuse 10-fold", 155,
    quote(kcv(log(zinc) ~ 1, ~ x + y, meuse, model = v_meuse, nfold = folds10_meuse)),
    quote(krige.cv(log(zinc) ~ 1, ~ x + y, meuse, model = v_meuse,
      nfold = folds10_meuse, verbose = FALSE)),
    10),
  run("sim LOO", n_sim,
    quote(kcv(z ~ 1, ~ x + y, sim, model = v_sim)),
    quote(krige.cv(z ~ 1, ~ x + y, sim, model = v_sim, verbose = FALSE)),
    2),
  run("sim 10-fold", n_sim,
    quote(kcv(z ~ 1, ~ x + y, sim, model = v_sim, nfold = folds10_sim)),
    quote(krige.cv(z ~ 1, ~ x + y, sim, model = v_sim, nfold = folds10_sim,
      verbose = FALSE)),
    3)
)
print(res, row.names = FALSE)
cat("\nR", as.character(getRversion()), "| gstat", as.character(packageVersion("gstat")),
  "| BLAS:", extSoftVersion()[["BLAS"]], "\n")
