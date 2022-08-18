function [s] = chopblas_pairwise_sum_mat( x, roundfunc, opts )
%CHOPBLAS_PAIRWISE_SUM_MAT Reduce the rows of a matrix to a vector with operation-level rounding
%
% Reduce a matrix to a vector by adding all the entries in a row together
% in a pairwise fasion, with rounding after each operation.

% Created by: Ian McInerney
% Created on: August 15, 2022
% SPDX-License-Identifier: BSD-2-Clause

s = x;
nr = size(s, 1);

while size(s, 2) > 1
    % Extract only an even number of elements to use in the sum this iteration,
    % leaving one element behind if there are an odd number. This element
    % is carried forward to the next iteration unchanged.
    ls = size(s, 2);
    isodd = mod( ls, 2 );

    t = zeros( nr, ceil(ls/2) );

    if isodd
        t(:,end) = s(:,end);
    end

    t(:,1:1:end-isodd) = roundfunc( s(:,1:2:(ls-1-isodd)) + s(:,2:2:(ls-isodd)), opts );
    s = t;
end

end
