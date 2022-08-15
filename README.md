# ChopBLAS - A Matlab library for mixed-precision and stochastically rounded linear algebra

[![CI](https://github.com/imciner2/ChopBLAS/actions/workflows/ci.yml/badge.svg)](https://github.com/imciner2/ChopBLAS/actions/workflows/ci.yml)
![License](https://img.shields.io/github/license/imciner2/ChopBLAS)

ChopBLAS is a MATLAB library that implements the functions contained inside the BLAS (Basic Linear Algebra Subsystem) specification
with each operation rounded using the [chop](https://github.com/higham/chop) MATLAB library.
This allows for the easy simulation of linear algebra operations/algorithms in low precision floating-point and with stochastic rounding.

## BLAS Functions

### Level 1

| Function | Operation                                | Description                                                           |
|----------|:-----------------------------------------|:----------------------------------------------------------------------|
| chscal   | $x_{out} = \alpha x$                     | Scale all entries of the vector x by alpha                            |
| chaxpy   | $x_{out} = \alpha x + y$                 | Add the scaled vector x to the vector y                               |
| chdot    | $x_{out} = x'y$                          | Compute the dot product between x and y                               |
| chnrm2   | $x_{out} = \lVert x \rVert_{2}$          | Compute the 2-norm of the vector x                                    |
| chasum   | $x_{out} = \sum_{i} \lvert x_{i} \rvert$ | Compute the sum of the absolute value of the elements of the vector x |

### Level 2

| Function | Operation                                                             | Description                                                           |
|----------|:----------------------------------------------------------------------|:----------------------------------------------------------------------|
| chgemv   | $x_{out} = \alpha A x + \beta y$ or $x_{out} = \alpha A' x + \beta y$ | Compute the matrix-vector product $Ax + y$                             |
| chtrmv   | $x_{out} = A x$ or $x_{out} = A' x$                                   | Compute the matrix-vector product $Ax$ when $A$ is a triangular matrix      |
| chtrsv   | $A x = b$ or $A' x = b$                                               | Compute the solution to the triangular system of equations given by $A$ when $A$ is a triangular matrix      |

## License

This library is licensed under the BSD 2-Clause license, provided inside the `LICENSE` file in this repository.
