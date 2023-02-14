function [s] = chaccum_compensated( x, roundfunc, opts )
%CHACCUM_COMPENSATED Perform a compensated summation of the vector using the given rounding mode
%
% Reduce a matrix to a vector by adding all the entries in a row in a recursive fashion
% using the Kahan compensated summation algorithm to incorporate the computational errors
% in the partial sums.
%
% Created by: Ian McInerney
% Created on: February 14, 2023
% SPDX-License-Identifier: BSD-2-Clause

% Setup the various parts that must differentiate between a column vector and row vector/matrix
if iscolumn( x )
    el  = @(i) x(i);
    len = length(x);
else
    el  = @(i) x(:,i);
    len = size(x, 2);
end

s = 0;
e = 0;

% Do the actual summation
for i=1:len
    temp = s;
    y = roundfunc( el(i) + e, opts );
    s = roundfunc( temp + y, opts );
    e = roundfunc( roundfunc( temp - s, opts ) + y, opts );
end

end
