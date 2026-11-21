# fastkcv itself is single-threaded, but a multithreaded BLAS (OpenBLAS by
# default starts one thread per core) makes the tests use several cores,
# which CRAN's check machines flag (CPU time > 2.5 times elapsed time). On
# CRAN (NOT_CRAN unset), limit the BLAS to one thread if RhpcBLASctl is
# available. Local and CI runs (NOT_CRAN = "true") keep the BLAS default.
if (!identical(Sys.getenv("NOT_CRAN"), "true") &&
  requireNamespace("RhpcBLASctl", quietly = TRUE)) {
  old_blas_threads <- RhpcBLASctl::blas_get_num_procs()
  RhpcBLASctl::blas_set_num_threads(1)
  withr::defer(
    RhpcBLASctl::blas_set_num_threads(old_blas_threads),
    teardown_env()
  )
}
