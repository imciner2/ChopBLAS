function [s] = chopblas_pairwise_sum( x, roundfunc, opts )
%CHOPBLAS_PAIRWISE_SUM Reduce a vector to a scalar with operation-level rounding
%
% Reduce a vector to a scalar by adding all the entries in a pairwise fasion, with
% rounding after each operation.

% Created by: Ian McInerney
% Created on: August 8, 2022
% SPDX-License-Identifier: BSD-2-Clause

lx = length(x);
isodd = mod( length(x), 2 );

% Extract only an even number of elements, leaving one element behind if there are
% an odd number
s = x(1:end-isodd);

while length(s) > 1
    s = roundfunc( s(1:2:end-1) + s(2:2:end), opts );
end

% The pairwise summation always saves the last element until the very end
if isodd
    s = roundfunc( s + x(end), opts );
end

end
