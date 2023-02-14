function [s] = chaccum_doubly_compensated( x, roundfunc, opts )
%CHACCUM_DOUBLY_COMPENSATED Perform a doubly compensated summation of the vector using the given rounding mode
%
% Reduce a matrix to a vector by adding all the entries in a row in a recursive fashion
% using the Priest doubly compensated summation algorithm to incorporate the computational errors
% in the partial sums.
%
% Created by: Ian McInerney
% Created on: February 14, 2023
% SPDX-License-Identifier: BSD-2-Clause

% Sort the input array such that the larger number is at the front.
% For QOL also turn column vectors into row vectors so we can only operate
% in a row-wise fashion.
if iscolumn( x )
    sx = sort(x', 'descend', 'ComparisonMethod', 'abs');
else
    sx = sort(x, 2, 'descend', 'ComparisonMethod', 'abs');
end

s  = sx(:, 1);
sl = sx(:, 1);
c  = 0;
cl = 0;

% Do the actual summation
for i=2:size(sx, 2)
    % Perform the compensated summation
    y = roundfunc( cl + sx(:, i), opts );
    u = roundfunc( sx(:, i) - roundfunc( y - cl, opts ), opts );
    t = roundfunc( y + sl, opts );
    v = roundfunc( y - roundfunc( t - sl, opts ), opts );
    z = roundfunc( u + v, opts );
    s = roundfunc( t + z, opts );
    c = roundfunc( z - roundfunc( s - t, opts ), opts );

    % Store values for next iteration
    sl = s;
    cl = c;
end

end
