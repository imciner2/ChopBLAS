function [s] = chopblas_recursive_sum_mat( x, roundfunc, opts )
%CHOPBLAS_RECURSIVE_SUM_MAT Reduce the rows of a matrix to a vector with operation-level rounding
%
% Reduce a matrix to a vector by adding all the entries in a row in a recursive fasion, with
% rounding after each operation. This sum starts with the first column and iterates
% in order to the last column.

% Created by: Ian McInerney
% Created on: August 15, 2022
% SPDX-License-Identifier: BSD-2-Clause

s = x(:,1);
for i=2:size(x,2)
    s = roundfunc( s + x(:,i), opts );
end

end
