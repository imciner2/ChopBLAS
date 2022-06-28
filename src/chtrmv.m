function [xout] = chtrmv( A, x, varargin )
%CHTRMV Perform a matrix-vector product for a triangular matrix
%
% Perform the matrix-vector computation
%   xout = A*x   - when 'Transpose' is false, or
%   xout = A'*x  - when 'Transpose' is true
% where A is a triangular matrix and x is a vector.
%
% Four optional boolean name-value pairs can be provided:
%   * 'LowerTriangular' - If true the matrix A is lower triangular, false means upper triangular
%                         Default: false
%   * 'Transpose'       - If true, the computation is A'*x instead of A*x
%                         Default: false
%   * 'UnitTriangular'  - If true, the matrix A is assumed to be unit triangular
%                         (1 on the main diagonal) and the main diagonal isn't used.
%                         Default: false
%   * 'Rounding'        - Function handle to the function that will perform the rounding operation.
%                         For more information on the interface 'roundfunc' must present, see the
%                         ChopBlas documentation.
%                         Default: @chop
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
% Specifying only roundopts will use the same rounding mode (given by
% roundopts) for both the multiplication and addition operations.
% Individual rounding modes for the multiplication and addition
% operations can be specified in the mulopts and addopts arguments,
% respectively.
%
% Usage:
%   [xout] = chtrmv( A, x, ... )
%   [xout] = chtrmv( A, x, roundopts, ... )
%   [xout] = chtrmv( A, x, mulopts, addopts, ... )

% Created by: Ian McInerney
% Created on: June 23, 2022
% License: BSD-2-Clause

%% Setup the argument parsing
isboolean = @(x) islogical(x) && isscalar(x);
p = inputParser;
p.StructExpand = false;
addOptional( p, 'mulopts', struct([]) );
addOptional( p, 'addopts', struct([]) );
addParameter( p, 'LowerTriangular', false, isboolean );
addParameter( p, 'Transpose', false, isboolean );
addParameter( p, 'UnitTriangular', false, isboolean );
addParameter( p, 'Rounding', @chop );

parse( p, varargin{:} )

mulopts   = p.Results.mulopts;
addopts   = p.Results.addopts;
trans     = p.Results.Transpose;
lowertri  = p.Results.LowerTriangular;
unitdiag  = p.Results.UnitTriangular;
roundfunc = p.Results.Rounding;

% Allow only the first to be specified and have it be used for both
if isempty(addopts) && ~isempty(mulopts)
    addopts = mulopts;
end

%% Preload the output vector with the main diagonal result
if unitdiag
    xout = x;
else
    xout = roundfunc( diag( A, 0 ).*x, mulopts );
end

[nr, nc] = size(A);

%% Actually compute the result
if trans
    if lowertri
        % Compute the transposed lower-triangular product
        for i=2:1:nr
            % The diagonal elements are already in xout, so skip them (and the last column)
            r = 1:1:(i-1);
            t = roundfunc( A(i, r)*x(i), mulopts );

            xout(r) = roundfunc( xout(r) + t', addopts );
        end
    else
        % Compute the transposed upper-triangular product
        for i=1:1:(nr-1)
            % The diagonal elements are already in xout, so skip them (and the first column)
            r = (i+1):1:nc;
            t = roundfunc( A(i, r)*x(i), mulopts );

            xout(r) = roundfunc( xout(r) + t', addopts );
        end
    end
else
    if lowertri
        % Compute the non-transposed lower-triangular product
        for i=1:1:(nc-1)
            % The diagonal elements are already in xout, so skip them (and the last column)
            r = (i+1):1:nr;
            t = roundfunc( A(r,i)*x(i), mulopts );

            xout(r) = roundfunc( xout(r) + t, addopts );
        end
    else
        % Compute the non-transposed upper-triangular product
        for i=2:1:nc
            % The diagonal elements are already in xout, so skip them (and the first column)
            r = 1:1:(i-1);
            t = roundfunc( A(r, i)*x(i), mulopts );

            xout(r) = roundfunc( xout(r) + t, addopts );
        end

%        for i=1:1:nr
%            r = i+1:1:nc;
%            t = chop( A(i, r).*x(r)', mulopts );
%
%            for j=1:1:length(t)
%                xout(i) = chop( xout(i) + t(j), addopts );
%            end
%        end
    end
end

end
