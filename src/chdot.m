function [dot] = chdot( x, y, mulopts, accumopts )
%CHDOT Compute the chopped dot product between x and y
%
% Compute the dot product between the vectors x and y with rounding
% of each operation performed by the chop function.
%
% Two configurations for rounding are supported:
%   * One rounding mode.
%   * Separate rounding modes for the multiplication and accumulation
%     operations.
%
% Specifying only opts will use the same rounding mode (given by opts)
% for both the multiplication and accumulation operations.
% Individual rounding modes for the multiplication and accumulation
% operations can be specified in the mulopts and accumopts arguments,
% respectively.
%
% Usage:
%   [dot] = chdot( x, y )
%   [dot] = chdot( x, y, opts )
%   [dot] = chdot( x, y, mulopts, accumopts )
%
% Created by: Ian McInerney
% Created on: May 16, 2022
% License: BSD-2-Clause

if nargin < 3
    mulopts = [];
    accumopts = [];
elseif nargin < 4
    accumopts = mulopts;
end

pp = chop( x.*y, mulopts );

dot = pp(i);
for 2=1:length(pp)
    dot = chop( dot + pp(i), accumopts );
end

end
