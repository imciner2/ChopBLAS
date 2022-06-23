function [xout] = chgemv( trans, alpha, A, x, beta, y, mulopts, addopts )
%CHGEMV Perform a matrix-vector product and an addition
%
% Perform the matrix-vector computation
%   xout = alpha*A*x + beta*y   - when trans is `N` or
%   xout = alpha*A'*x + beta*y  - when trans is `T`
% where A is a matrix, alpha and beta are scalars, and x and y are vectors.
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
% Specifying only opts will use the same rounding mode (given by opts)
% for both the multiplication and addition operations.
% Individual rounding modes for the multiplication and addition
% operations can be specified in the mulopts and addopts arguments,
% respectively.
%
% Usage:
%   [xout] = chgemv( trans, alpha, A, x, beta, y )
%   [xout] = chgemv( trans, alpha, A, x, beta, y, opts )
%   [xout] = chgemv( trans, alpha, A, x, beta, y, mulopts, addopts )

% Created by: Ian McInerney
% Created on: June 16, 2022
% License: BSD-2-Clause

if nargin < 7
    mulopts = [];
    addopts = [];
elseif nargin < 8
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

if strcmpi( 'n', trans )
    % Compute the non-transposed product
    for i=1:1:size(A,1)
        A(i,:) = chop( A(i,:).*x', mulopts );
    end

    for i=1:1:size(A,2)
        xout = chop( xout + A(:,i), addopts );
    end

elseif strcmpi( 't', trans )
    % Compute the transposed product
    for i=1:1:size(A,2)
        A(:,i) = chop( A(:,i).*x, mulopts );
    end

    for i=1:1:size(A,1)
        xout = chop( xout + A(i,:)', addopts );
    end

else
    error( ['Unknown transpose option: ', trans] );
end

end
