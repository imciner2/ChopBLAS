function [nrm] = chnrm2( x, varargin )
%CHNRM2 Compute the 2-norm of the vector x with operation-level rounding
%
% Compute the 2-norm of the vector x with rounding after each operation.
%
% This function supports the following optional name-value arguments:
%   * 'Rounding'    - Function handle to the function that will perform the rounding operation.
%                     For more information on the interface this function must have, see the
%                     ChopBlas documentation.
%                     Default: @chop
%   * 'Accumulator' - The algorithm to use when performing the additions.
%                     Default: @chaccum_recursive
%
% Two configurations for rounding are supported:
%   * One rounding mode.
%   * Separate rounding modes for the multiplication, accumulation
%     and square-root operations.
%
% Specifying only opts will use the same rounding mode (given by opts)
% for all three operations.
% Individual rounding modes for the multiplication, accumulation and
% square-root operations can be specified in the mulopts, addopts,
% and sqrtopts arguments, respectively.
%
% Usage:
%   [nrm] = CHNRM2( x, ... )
%   [nrm] = CHNRM2( x, opts, ... )
%   [nrm] = CHNRM2( x, mulopts, addopts, sqrtopts, ... )

% Created by: Ian McInerney
% Created on: May 20, 2022
% SPDX-License-Identifier: BSD-2-Clause


%% Setup the argument parsing
p = inputParser;
p.StructExpand = false;
addOptional( p, 'mulopts', struct([]) );
addOptional( p, 'addopts', struct([]) );
addOptional( p, 'sqrtopts', struct([]) );
addParameter( p, 'Rounding', @chop );
addParameter( p, 'Accumulator', @chaccum_recursive );

parse( p, varargin{:} )

accum     = p.Results.Accumulator;
mulopts   = p.Results.mulopts;
addopts   = p.Results.addopts;
sqrtopts  = p.Results.sqrtopts;
roundfunc = p.Results.Rounding;

% Allow only the first to be specified and have it be used for both
if ( isempty(addopts) || isempty(sqrtopts) ) && ~isempty(mulopts)
    addopts  = mulopts;
    sqrtopts = mulopts;
end

pp = roundfunc( x.*x, mulopts );

dot = accum( pp, roundfunc, addopts );

nrm = roundfunc( sqrt( dot ), sqrtopts );

end
