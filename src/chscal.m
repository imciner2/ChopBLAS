function [xout] = chscal( alpha, x, varargin )
%CHSCAL Scale all entries of the vector x by alpha with operation-level rounding
%
% Scale each element of the vector x by the scale factor alpha, then
% round each entry.
%
% This function supports the following optional name-value arguments:
%   * 'Rounding' - Function handle to the function that will perform the rounding operation.
%                  For more information on the interface 'roundfunc' must present, see the
%                  ChopBlas documentation.
%                  Default: @chop
%
% If opts is empty/not provided, the rounding will be done using
% the global options for the rounding function. Otherwise, the rounding will use opts.
%
% Usage:
%   [xout] = CHSCAL( alpha, x, ... )
%   [xout] = CHSCAL( alpha, x, opts, ... )

% Created by: Ian McInerney
% Created on: May 16, 2022
% SPDX-License-Identifier: BSD-2-Clause

%% Setup the argument parsing
p = inputParser;
p.StructExpand = false;
addOptional( p, 'opts', struct([]) );
addParameter( p, 'Rounding', @chop );

parse( p, varargin{:} )

opts   = p.Results.opts;
roundfunc = p.Results.Rounding;

if ~isscalar( alpha )
    error( "chscal:AlphaMustBeScalar", "alpha must be a scalar." );
end

xout = roundfunc( alpha.*x, opts );

end
