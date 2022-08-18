function [s] = chopblas_sorted_sum_vec( x, dir, roundfunc, opts )
%CHOPBLAS_SORTED_SUM_VEC Reduce a vector to a scalar with operation-level rounding
%
% Reduce a vector to a scalar by adding all the entries in a recursive fasion, with
% rounding after each operation. First the vector x is sorted in the specified order
% then the elements are recursively added.
%
% If dir==1, sort in increasing direction, otherwise sort in decreasing direction.

% Created by: Ian McInerney
% Created on: August 15, 2022
% SPDX-License-Identifier: BSD-2-Clause

if dir == 1
    s = sort(x, 'ascend', 'ComparisonMethod', 'abs');
else
    s = sort(x, 'descend', 'ComparisonMethod', 'abs');
end

while length(s) > 1
    t = roundfunc( s(1) + s(2), opts );
    s = s(2:end);
    s(1) = t;
end

end
