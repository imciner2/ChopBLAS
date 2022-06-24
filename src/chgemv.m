function [xout] = chgemv( alpha, A, x, beta, y, varargin )
%CHGEMV Perform a matrix-vector product and an addition
%
% Perform the matrix-vector computation
%   xout = alpha*A*x + beta*y   - when 'Transpose' is false, or
%   xout = alpha*A'*x + beta*y  - when 'Transpose' is true
% where A is a matrix, alpha and beta are scalars, x and y are vectors, and
% 'Transpose' is an optional argument.
% If beta is zero or y is empty, then the addition with y is not performed.
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
%
% Usage:
%   [xout] = chgemv( alpha, A, x, beta, y, ... )
%   [xout] = chgemv( alpha, A, x, beta, y, roundopts, ... )
%   [xout] = chgemv( alpha, A, x, beta, y, mulopts, addopts, ... )

% Created by: Ian McInerney
% Created on: June 16, 2022
% License: BSD-2-Clause


%% Setup the argument parsing
isboolean = @(x) islogical(x) && isscalar(x);
p = inputParser;
p.StructExpand = false;
addOptional( p, 'mulopts', [] );
addOptional( p, 'addopts', [] );
addParameter( p, 'Transpose', false, isboolean );

parse( p, varargin{:} )

mulopts  = p.Results.mulopts;
addopts  = p.Results.addopts;
trans    = p.Results.Transpose;

% Allow only the first to be specified and have it be used for both
if isempty(addopts) && ~isempty(mulopts)
    addopts = mulopts;
end

% Initialize output using scaled y vector
if ( beta == 1 ) && ( ~isempty(y) )
    xout = y;
elseif ( beta ~= 0 ) && ( ~isempty(y) )
    xout = chop( beta.*y, mulopts );
else
    xout = zeros(length(x), 1);
end

if alpha == 0
    % Short circuit return
    return;
elseif alpha ~= 1
    % Apply the scaling on the matrix-vector product
    x = chop( alpha.*x, mulopts );
end

if trans
    % Compute the transposed product
    for i=1:1:size(A,2)
        A(:,i) = chop( A(:,i).*x, mulopts );
    end

    for i=1:1:size(A,1)
        xout = chop( xout + A(i,:)', addopts );
    end
else
    % Compute the non-transposed product
    for i=1:1:size(A,1)
        A(i,:) = chop( A(i,:).*x', mulopts );
    end

    for i=1:1:size(A,2)
        xout = chop( xout + A(:,i), addopts );
    end
end

end
