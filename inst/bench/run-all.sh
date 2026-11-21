#!/bin/sh
# Runs every benchmark under OpenBLAS and then the reference BLAS, from a
# clean optimized install, and prints the summary tables. From the package
# root:  sh inst/bench/run-all.sh
# OpenBLAS runs first: it creates the cached data sets and the n = 20,000
# prepared object that the reference-BLAS run reuses. Logs go to
# inst/bench/cache/. The whole run takes several hours.
# Optional arguments run only the given steps, e.g.
#   sh inst/bench/run-all.sh openblas/rss openblas/buffer-large refblas/folds
LIB=/usr/lib/x86_64-linux-gnu
OPENBLAS=$LIB/openblas-pthread:/usr/lib/R/lib:$LIB
REFBLAS=$LIB/blas:$LIB/lapack:/usr/lib/R/lib:$LIB

rm -f src/*.o src/*.so
R CMD INSTALL --preclean . > inst/bench/cache/install.log 2>&1 || exit 1
steps="$*"
if [ -z "$steps" ]; then
  for blas in openblas refblas; do
    for s in folds dice auto-rule rss buffer-large multi; do
      steps="$steps $blas/$s"
    done
  done
fi
for step in $steps; do
  blas=${step%%/*}
  s=${step#*/}
  if [ $blas = openblas ]; then path=$OPENBLAS; else path=$REFBLAS; fi
  echo "$(date '+%H:%M:%S') $blas $s"
  R_LD_LIBRARY_PATH=$path Rscript inst/bench/bench-$s.R \
    > inst/bench/cache/$s-$blas.log 2>&1 || echo "  failed, see the log"
done
echo "$(date '+%H:%M:%S') done"
Rscript inst/bench/summarize.R
