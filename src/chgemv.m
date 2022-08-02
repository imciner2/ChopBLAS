function [xout] = chgemv( alpha, A, x, beta, y, varargin )
%CHGEMV Perform a matrix-vector product and an addition with operation-level rounding
%
% Perform the matrix-vector computation
%   xout = alpha*A*x + beta*y   - when 'Transpose' is false, or
%   xout = alpha*A'*x + beta*y  - when 'Transpose' is true
% where A is a matrix, alpha and beta are scalars, x and y are vectors, and
% 'Transpose' is an optional argument.
% If beta is zero or y is empty, then the addition with y is not performed.
%
% This function supports the following optional name-value arguments
%   * 'Transpose'       - If true, the computation is A'*x instead of A*x
%                         Default: false
%   * 'Rounding'        - Function handle to the function that will perform the rounding operation.
%                         For more information on the interface 'roundfunc' must present, see the
%                         ChopBlas documentation.
%                         Default: @chop
%
% The order of operations for this function are as follows:
%   1) Populate xout with the y vector
%      If beta==1, then no multiplication/rounding is done
%      if beta==0 or y is empty, xout is loaded with 0.
%   2) Scale x by alpha (if alpha==1, no multiplication/rounding is done)
%   3) Compute A*x or A'*x with the scaled x, accumulating onto xout.
%
% Two configurations for rounding are supported:
%   * One rounding mode.
%   * Separate rounding modes for the multiplication and addition
%     operations.
%
% Specifying only roundopts will use the same rounding mode (given by
% roundopts) for both the multiplication and addition operations.
% Individual rounding modes for the multiplication and addition
% operations can be specified in the mulopts and addopts arguments,
% respectively.

% Usage:
%   [xout] = CHGEMV( alpha, A, x, beta, y, ... )
%   [xout] = CHGEMV( alpha, A, x, beta, y, roundopts, ... )
%   [xout] = CHGEMV( alpha, A, x, beta, y, mulopts, addopts, ... )

% Created by: Ian McInerney
% Created on: June 16, 2022
% License: BSD-2-Clause


%% Setup the argument parsing
isboolean = @(x) islogical(x) && isscalar(x);
p = inputParser;
p.StructExpand = false;
addOptional( p, 'mulopts', struct([]) );
addOptional( p, 'addopts', struct([]) );
addParameter( p, 'Transpose', false, isboolean );
addParameter( p, 'Rounding', @chop );

parse( p, varargin{:} )

mulopts   = p.Results.mulopts;
addopts   = p.Results.addopts;
trans     = p.Results.Transpose;
roundfunc = p.Results.Rounding;

% Allow only the first to be specified and have it be used for both
if isempty(addopts) && ~isempty(mulopts)
    addopts = mulopts;
end

% Initialize output using scaled y vector
if ( beta == 1 ) && ( ~isempty(y) )
    xout = y;
elseif ( beta ~= 0 ) && ( ~isempty(y) )
    xout = roundfunc( beta.*y, mulopts );
else
    xout = zeros(length(x), 1);
end

if alpha == 0
    % Short circuit return
    return;
elseif alpha ~= 1
    % Apply the scaling on the matrix-vector product
    x = roundfunc( alpha.*x, mulopts );
end

if trans
    % Matrix indexing needed to compute the transposed product
    matind = @(i) A(i,:)';
else
    % Matrix indexing needed to compute the non-transposed product
    matind = @(i) A(:,i);
end

for i=1:1:length(x)
    % matind() will return a column vector of the matrix elements in the proper position
    % for the final add (based on the transpose option)
    t = roundfunc( matind(i).*x(i), mulopts );
    xout = roundfunc( xout + t, addopts );
end

end
