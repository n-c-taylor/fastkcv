# fastkcv: Fast Exact Cross-Validation for Global Kriging

Fast, exact leave-one-out, k-fold and buffered leave-one-out
cross-validation for global (unique neighborhood) simple, ordinary and
universal kriging, computed from a single factorization of the
covariance matrix (Dubrule (1983)
[doi:10.1007/BF01033232](https://doi.org/10.1007/BF01033232) ;
Ginsbourger and Schärer (2024)
[doi:10.1080/10618600.2024.2353633](https://doi.org/10.1080/10618600.2024.2353633)
). Takes the same inputs as 'gstat' krige.cv() and reproduces its output
to numerical precision, with nugget effects, nested and anisotropic
variogram models. Includes random, spatial block and k-means fold
builders and summary statistics for comparing fold designs.

## See also

Useful links:

- <https://n-c-taylor.github.io/fastkcv/>

- <https://github.com/n-c-taylor/fastkcv>

- Report bugs at <https://github.com/n-c-taylor/fastkcv/issues>

## Author

**Maintainer**: Nicholas Taylor <nicholasclevelandtaylor@gmail.com>

Authors:

- Nicholas Taylor <nicholasclevelandtaylor@gmail.com>
