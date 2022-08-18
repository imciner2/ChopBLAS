function [s] = chopblas_insertion_sum_vec( x, dir, roundfunc, opts )
%CHOPBLAS_INSERTION_SUM_VEC Reduce a vector to a scalar with operation-level rounding
%
% Reduce a vector to a scalar by adding all the entries in a recursive fasion, with
% rounding after each operation. First the vector x is sorted in the specified order
% then the first two elements are added together and its result is placed in the
% proper sorted position, and the process repeats.
%
% If dir==1, sort in increasing direction, otherwise sort in decreasing direction.

% Created by: Ian McInerney
% Created on: August 18, 2022
% SPDX-License-Identifier: BSD-2-Clause

if dir == 1
    sfunc = @(x) sort(x, 'ascend', 'ComparisonMethod', 'abs');
else
    sfunc = @(x) sort(x, 'descend', 'ComparisonMethod', 'abs');
end

s = sfunc(x);

while length(s) > 1
    t = roundfunc( s(1) + s(2), opts );
    s = s(2:end);
    s(1) = t;
    s = sfunc(s);
end

end
