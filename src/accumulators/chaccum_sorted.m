function [s] = chaccum_sorted( x, roundfunc, opts, dir )
%CHOPBLAS_SORTED_SUM_MAT Accumulate into a column vector by iterating through the sorted columns of the matrix
%
% Reduce a matrix to a column vector by adding all the entries in a row recursive fasion, with
% rounding after each operation. First the vector x is sorted in the specified order
% then the elements are recursively added.
%
% 'dir' can either be 'ascend' or 'descend' to sort by ascending or descending values, respectively.

% Created by: Ian McInerney
% Created on: August 15, 2022
% SPDX-License-Identifier: BSD-2-Clause

% Sorting is done along either the vector length or the rows of the matrix
if iscolumn( x )
    sx  = sort(x, dir, 'ComparisonMethod', 'abs');
else
    sx  = sort(x, 2, dir, 'ComparisonMethod', 'abs');
end

% This really is just recursive accumulation after the sorting
s = chaccum_recursive( sx, roundfunc, opts );

end
