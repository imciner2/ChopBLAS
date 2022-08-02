function [dot] = chdot( x, y, varargin )
%CHDOT Compute the dot product between x and y with operation-level rounding
%
% Compute the dot product between the vectors x and y with rounding
% of each operation.
%
% This function supports the following optional name-value arguments:
%   * 'Rounding' - Function handle to the function that will perform the rounding operation.
%                  For more information on the interface 'roundfunc' must present, see the
%                  ChopBlas documentation.
%                  Default: @chop
%
% Two configurations for rounding are supported:
%   * One rounding mode.
%   * Separate rounding modes for the multiplication and addition
%     operations.
%
% Specifying only opts will use the same rounding mode (given by opts)
% for both the multiplication and addition operations.
% Individual rounding modes for the multiplication and addition
% operations can be specified in the mulopts and addopts arguments,
% respectively.
%
% Usage:
%   [dot] = CHDOT( x, y, ... )
%   [dot] = CHDOT( x, y, opts, ... )
%   [dot] = CHDOT( x, y, mulopts, addopts, ... )

% Created by: Ian McInerney
% Created on: May 16, 2022
% License: BSD-2-Clause

%% Setup the argument parsing
p = inputParser;
p.StructExpand = false;
addOptional( p, 'mulopts', struct([]) );
addOptional( p, 'addopts', struct([]) );
addParameter( p, 'Rounding', @chop );

parse( p, varargin{:} )

mulopts   = p.Results.mulopts;
addopts   = p.Results.addopts;
roundfunc = p.Results.Rounding;

% Allow only the first to be specified and have it be used for both
if isempty(addopts) && ~isempty(mulopts)
    addopts = mulopts;
end

sx = size( x );
sy = size( y );

if all( sx == sy )
    % Both either column or row vectors
    pp = roundfunc( x.*y, mulopts );
elseif all( flip(sx) == sy )
    % One column and one row vector given
    pp = roundfunc( x'.*y, mulopts );
else
    % Sizes aren't compatible
    errmsg = strcat( "Vectors x and y must be compatible sizes - ",...
                     "[", num2str(sx(1)), "x", num2str(sx(2)), "] vs. ", ...
                     "[", num2str(sy(1)), "x", num2str(sy(2)), "]." );
    error( "chdot:xyMustBeCompatibleSize", errmsg );
end

dot = pp(1);
for i=2:length(pp)
    dot = roundfunc( dot + pp(i), addopts );
end

end
