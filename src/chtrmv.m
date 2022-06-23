function [xout] = chtrmv( uplo, trans, unitdiag, A, x, mulopts, addopts )
%CHTRMV Perform a matrix-vector product for a triangular matrix
%
% Perform the matrix-vector computation
%   xout = A*x   - when trans is `N` or
%   xout = A'*x  - when trans is `T`
% where A is a triangular matrix and x is a vector. If uplo is `U` (`L`) then
% A is an upper(lower)-triangular matrix. If unitdiag is `U` then A is assumed
% to be unit triangular (1 on the main diagonal), if unitdiag is `N`, then A is
% not assumed to be unit triangular.
%
% The order of operations for this function are as follows:
%   1) Populate xout with the diagonal multiplication
%      If unitdiag is `U`, then xout = x with no rounding.
%      If unitdiag is `N`, then x is multiplied by the main diagonal and stored in xout.
%   2) Compute the remaining terms of A*x or A'*x, accumulating onto xout.
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
%   [xout] = chtrmv( uplo, trans, unitdiag, A, x )
%   [xout] = chtrmv( uplo, trans, unitdiag, A, x, opts )
%   [xout] = chtrmv( uplo, trans, unitdiag, A, x, mulopts, addopts )

% Created by: Ian McInerney
% Created on: June 23, 2022
% License: BSD-2-Clause

if nargin < 6
    mulopts = [];
    addopts = [];
elseif nargin < 7
    addopts = mulopts;
end

% Preload the output vector with the main diagonal result
if strcmpi( 'n', unitdiag )
    xout = chop( diag( A, 0 ).*x, mulopts );
elseif strcmpi( 'u', unitdiag )
    xout = x;
else
    error( ['Unknown unitdiag option: ', unitdiag] );
end

[nr, nc] = size(A);

if strcmpi( 'n', trans )
    if strcmpi( 'u', uplo )
        % Compute the non-transposed upper-triangular product
        for i=2:1:nc
            % The diagonal elements are already in xout, so skip them (and the first column)
            r = 1:1:(i-1);
            t = chop( A(r, i)*x(i), mulopts );

            xout(r) = chop( xout(r) + t, addopts );
        end

%        for i=1:1:nr
%            r = i+1:1:nc;
%            t = chop( A(i, r).*x(r)', mulopts );
%
%            for j=1:1:length(t)
%                xout(i) = chop( xout(i) + t(j), addopts );
%            end
%        end

    elseif strcmpi( 'l', uplo )
        % Compute the non-transposed lower-triangular product
        for i=1:1:(nc-1)
            % The diagonal elements are already in xout, so skip them (and the last column)
            r = (i+1):1:nr;
            t = chop( A(r,i)*x(i), mulopts );

            xout(r) = chop( xout(r) + t, addopts );
        end

    else
        error( ['Unknown uplo option: ', uplo] );
    end

elseif strcmpi( 't', trans )
    if strcmpi( 'u', uplo )
        % Compute the transposed upper-triangular product
        for i=1:1:(nr-1)
            % The diagonal elements are already in xout, so skip them (and the first column)
            r = (i+1):1:nc;
            t = chop( A(i, r)*x(i), mulopts );

            xout(r) = chop( xout(r) + t', addopts );
        end

    elseif strcmpi( 'l', uplo )
        % Compute the transposed lower-triangular product
        for i=2:1:nr
            % The diagonal elements are already in xout, so skip them (and the last column)
            r = 1:1:(i-1);
            t = chop( A(i, r)*x(i), mulopts );

            xout(r) = chop( xout(r) + t', addopts );
        end

    else
        error( ['Unknown uplo option: ', uplo] );
    end

else
    error( ['Unknown transpose option: ', trans] );
end

end
