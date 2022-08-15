function [s] = chopblas_pairwise_sum_vec( x, roundfunc, opts )
%CHOPBLAS_PAIRWISE_SUM_VEC Reduce a vector to a scalar with operation-level rounding
%
% Reduce a vector to a scalar by adding all the entries in a pairwise fasion, with
% rounding after each operation.

% Created by: Ian McInerney
% Created on: August 8, 2022
% SPDX-License-Identifier: BSD-2-Clause

s = x;

while length(s) > 1
    % Extract only an even number of elements to use in the sum this iteration,
    % leaving one element behind if there are an odd number. This element
    % is carried forward to the next iteration unchanged.
    ls = length(s);
    isodd = mod( ls, 2 );

    t = zeros( ceil(ls/2), 1 );

    if isodd
        t(end) = s(end);
    end

    t(1:1:end-isodd) = roundfunc( s(1:2:(ls-1-isodd)) + s(2:2:(ls-isodd)), opts );
    s = t;
end

end
