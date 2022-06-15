# ChopBLAS - A Matlab library for mixed-precision and stochastically rounded linear algebra

[![CI](https://github.com/imciner2/ChopBLAS/actions/workflows/ci.yml/badge.svg)](https://github.com/imciner2/ChopBLAS/actions/workflows/ci.yml)
![License](https://img.shields.io/github/license/imciner2/ChopBLAS)

ChopBLAS is a MATLAB library that implements the functions contained inside the BLAS (Basic Linear Algebra Subsystem) specification
with each operation rounded using the [chop](https://github.com/higham/chop) MATLAB library.
This allows for the easy simulation of linear algebra operations/algorithms in low precision floating-point and with stochastic rounding.

## BLAS Functions

### Level 1

| Function | Description                                                           |
|----------|:----------------------------------------------------------------------|
| chscal   | Scale all entries of the vector x by alpha                            |
| chaxpy   | Add the scaled vector x to the vector y                               |
| chdot    | Compute the chopped dot product between x and y                       |
| chnrm2   | Compute the 2-norm of the vector x                                    |
| chasum   | Compute the sum of the absolute value of the elements of the vector x |


## License

This library is licensed under the 2-Clause BSD license, provided inside the `LICENSE` file in this repository.
