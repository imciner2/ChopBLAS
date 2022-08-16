function [X] = chger( alpha, x, y, A, varargin )
%CHGER Perform the rank-1 update of A with operation-level rounding
%
% Perform the rank-1 update of A using the vectors x and y such that
%     X = alpha*x*y' + A
% with x an m element column vector, y an n element column vector, alpha a
% scalar and A an m by n matrix.
%
% This function supports the following optional name-value arguments:
%   * 'Rounding' - Function handle to the function that will perform the rounding operation.
%                  For more information on the interface this function must have, see the
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
%   [X] = CHGER( alpha, x, y, A, ... )
%   [X] = CHGER( alpha, x, y, A, opts, ... )
%   [X] = CHGER( alpha, x, y, A, mulopts, addopts, ... )

% Created by: Ian McInerney
% Created on: August 16, 2022
% SPDX-License-Identifier: BSD-2-Clause

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

% Verify arguments
sx = size( x );
sy = size( y );
sA = size( A );

if ~isscalar( alpha )
    error( "chger:AlphaMustBeScalar", "alpha must be a scalar." );
end
if sx(1) ~= sA(1)
    errmsg = strcat( "Length of x must be same as the number of rows of A - ",...
                     num2str(sx(1)), " versus ", num2str(sA(1)) );
    error( "chger:xAMustHaveCompatibleSize", errmsg );
end
if sy(1) ~= sA(2)
    errmsg = strcat( "Length of y must be same as the number of columns of A - ",...
                     num2str(sx(1)), " versus ", num2str(sA(2)) );
    error( "chger:yAMustHaveCompatibleSize", errmsg );
end

if alpha == 0
    % If we get an alpha of zero, then just return A unrounded
    X = A;
else
    X = roundfunc( x*y', mulopts );

    if alpha ~= 1
        X = roundfunc( alpha.*X, mulopts );
    end

    X = roundfunc( X + A, addopts );
end

end
