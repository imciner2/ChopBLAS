function [s] = chaccum_recursive( x, roundfunc, opts )
%CHACCUM_RECURSIVE Accumulate into a column vector by just iterating through the columns of the matrix
%
% Reduce a matrix to a vector by adding all the entries in a row in a recursive fasion, with
% rounding after each operation. This sum starts with the first column and iterates
% in order to the last column.

% Created by: Ian McInerney
% Created on: August 15, 2022
% SPDX-License-Identifier: BSD-2-Clause

% Setup the various parts that must differentiate between a column vector and row vector/matrix
if iscolumn( x )
    el  = @(i) x(i);
    len = length(x);
else
    el  = @(i) x(:,i);
    len = size(x, 2);
end

s = el(1);

% Do the actual summation
for i=2:len
    s = roundfunc( s + el(i), opts );
end

end
