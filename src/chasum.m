function [xout] = chasum( x, opts )
%CHASUM Compute the sum of the absolute value of the elements of the vector x
%
% Take the absolute value of each element in x, and then compute the sum of
% all the elements using the rounding options given in opts. If no options are
% given, then the global chop rounding options will be used.
%
% Usage:
%   [xout] = chasum( x )
%   [xout] = chasum( x, opts )

% Created by: Ian McInerney
% Created on: May 20, 2022
% License: BSD-2-Clause

if nargin < 2
    opts = [];
end

xout = abs( x(1) );
for i=2:length(pp)
    xout = chop( xout + abs( x(i) ), opts );
end

end
