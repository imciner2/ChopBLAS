function [nrm] = chnrm2( x, varargin )
%CHNRM2 Compute the 2-norm of the vector x with operation-level rounding
%
% Compute the 2-norm of the vector x with rounding after each operation.
%
% This function supports the following optional name-value arguments:
%   * 'Rounding' - Function handle to the function that will perform the rounding operation.
%                  For more information on the interface 'roundfunc' must present, see the
%                  ChopBlas documentation.
%                  Default: @chop
%
% Two configurations for rounding are supported:
%   * One rounding mode.
%   * Separate rounding modes for the multiplication, accumulation
%     and square-root operations.
%
% Specifying only opts will use the same rounding mode (given by opts)
% for all three operations.
% Individual rounding modes for the multiplication, accumulation and
% square-root operations can be specified in the mulopts, accumopts,
% and sqrtopts arguments, respectively.
%
% Usage:
%   [nrm] = CHNRM2( x, ... )
%   [nrm] = CHNRM2( x, opts, ... )
%   [nrm] = CHNRM2( x, mulopts, accumopts, sqrtopts, ... )

% Created by: Ian McInerney
% Created on: May 20, 2022
% License: BSD-2-Clause


%% Setup the argument parsing
p = inputParser;
p.StructExpand = false;
addOptional( p, 'mulopts', struct([]) );
addOptional( p, 'addopts', struct([]) );
addOptional( p, 'sqrtopts', struct([]) );
addParameter( p, 'Rounding', @chop );

parse( p, varargin{:} )

mulopts   = p.Results.mulopts;
addopts   = p.Results.addopts;
roundfunc = p.Results.Rounding;

% Allow only the first to be specified and have it be used for both
if ( isempty(addopts) || isempty(sqrtopts) ) && ~isempty(mulopts)
    addopts  = mulopts;
    sqrtopts = mulopts;
end

pp = roundfunc( x.*x, mulopts );

dot = pp(1);
for i=2:length(pp)
    dot = roundfunc( dot + pp(i), accumopts );
end

nrm = roundfunc( sqrt( dot ), sqrtopts );

end
