function [s] = chopblas_sorted_sum_mat( x, dir, roundfunc, opts )
%CHOPBLAS_SORTED_SUM_MAT Reduce a matrix to a column vector with operation-level rounding
%
% Reduce a matrix to a column vector by adding all the entries in a row recursive fasion, with
% rounding after each operation. First the vector x is sorted in the specified order
% then the elements are recursively added.
%
% If dir==1, sort in increasing direction, otherwise sort in decreasing direction.

% Created by: Ian McInerney
% Created on: August 15, 2022
% SPDX-License-Identifier: BSD-2-Clause

if dir == 1
    s = sort(x, 2, 'ascend', 'ComparisonMethod', 'abs');
else
    s = sort(x, 2, 'descend', 'ComparisonMethod', 'abs');
end

while size(s,2) > 1
    t = roundfunc( s(:,1) + s(:,2), opts );
    s = s(:,2:end);
    s(:,1) = t;
end

end
