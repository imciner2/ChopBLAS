function [xout] = chscal( alpha, x, opts )
%CHSCAL Scale all entries of the vector x by alpha
%
% Scale each element of the vector x by the scale factor alpha, then
% use chop to round each entry.
%
% If opts is empty/not provided, the rounding will be done using
% the global options for chop. Otherwise, the rounding will use opts.
%
% Usage:
%   [xout] = chscal( alpha, x )
%   [xout] = chscal( alpha, x, opts )
%
% Created by: Ian McInerney
% Created on: May 16, 2022
% License: BSD-2-Clause

if nargin < 3
    opts = [];
end

xout = chop( alpha.*x, opts );

end
