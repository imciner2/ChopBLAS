function [s] = chopblas_recursive_sum_vec( x, roundfunc, opts )
%CHOPBLAS_RECURSIVE_SUM_VEC Reduce a vector to a scalar with operation-level rounding
%
% Reduce a vector to a scalar by adding all the entries in a recursive fasion, with
% rounding after each operation. This sum starts with the first element and iterates
% in order to the last element.

% Created by: Ian McInerney
% Created on: August 8, 2022
% SPDX-License-Identifier: BSD-2-Clause

s = x(1);
for i=2:length(x)
    s = roundfunc( s + x(i), opts );
end

end
