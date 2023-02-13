function [xout] = chasum( x, varargin )
%CHASUM Compute the sum of the absolute value of the elements of the vector x with operation-level rounding
%
% Take the absolute value of each element in x, and then compute the sum of
% all the elements with rounding after each operation using the rounding
% options given in opts. If no options are given, then the global rounding
% options will be used.
%
% This function supports the following optional name-value arguments:
%   * 'Rounding'    - Function handle to the function that will perform the rounding operation.
%                     For more information on the interface this function must have, see the
%                     ChopBlas documentation.
%                     Default: @chop
%   * 'Accumulator' - The algorithm to use when performing the additions.
%                     Default: @chaccum_recursive
%
% Usage:
%   [xout] = CHASUM( x, ... )
%   [xout] = CHASUM( x, opts, ... )

% Created by: Ian McInerney
% Created on: May 20, 2022
% SPDX-License-Identifier: BSD-2-Clause

%% Setup the argument parsing
p = inputParser;
p.StructExpand = false;
addOptional( p, 'opts', struct([]) );
addParameter( p, 'Rounding', @chop );
addParameter( p, 'Accumulator', @chaccum_recursive );

parse( p, varargin{:} )

opts      = p.Results.opts;
accum     = p.Results.Accumulator;
roundfunc = p.Results.Rounding;

x = abs( x );

xout = accum( x, roundfunc, opts );

end
