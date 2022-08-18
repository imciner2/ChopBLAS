function [s] = chopblas_insertion_sum_mat( x, dir, roundfunc, opts )
%CHOPBLAS_INSERTION_SUM_MAT Reduce a matrix to a vector with operation-level rounding
%
% Reduce a matrix to a column vector by adding all the entries in a row together in
% a recursive fasion, with rounding after each operation. First the rows of the matrix
% x are sorted in the specified order, then the first two elements are added together
% and the result is placed in the proper sorted position, and the process repeats.
%
% If dir==1, sort in increasing direction, otherwise sort in decreasing direction.

% Created by: Ian McInerney
% Created on: August 18, 2022
% SPDX-License-Identifier: BSD-2-Clause

if dir == 1
    sfunc = @(x) sort(x, 2, 'ascend', 'ComparisonMethod', 'abs');
else
    sfunc = @(x) sort(x, 2, 'descend', 'ComparisonMethod', 'abs');
end

s = sfunc(x);

while size(s,2) > 1
    t = roundfunc( s(:,1) + s(:,2), opts );
    s = s(:,2:end);
    s(:,1) = t;
    s = sfunc(s);
end

end
