#' @keywords internal
"_PACKAGE"

## usethis namespace: start
#' @useDynLib fastkcv, .registration = TRUE
#' @importFrom Rcpp sourceCpp
## usethis namespace: end

# Models come from gstat::vgm(), and the documentation and examples rely on
# gstat, so it stays in Imports although only the tests call its functions.
#' @importFrom gstat vgm
NULL

## mockable bindings: start
## mockable bindings: end
NULL
