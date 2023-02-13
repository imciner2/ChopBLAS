function [x] = chtrsv( A, b, varargin )
%CHTRSV Solves a triangular system of linear equations
%
% Solves the systems of equations
%   b = A*x   - when 'Transpose' is false, or
%   b = A'*x  - when 'Transpose' is true
% where A is a triangular matrix, b is a vector and x is unknown.
%
% This function supports the following optional name-value arguments
%   * 'LowerTriangular' - If true the matrix A is lower triangular, false means upper triangular
%                         Default: false
%   * 'Transpose'       - If true, the system is A'*x instead of A*x
%                         Default: false
%   * 'UnitTriangular'  - If true, the matrix A is assumed to be unit triangular
%                         (1 on the main diagonal) and the given main diagonal isn't used.
%                         Default: false
%   * 'RoundAll'        - If true, ensure all elements of b have been rounded using the divison rounding
%                         mode before being substituted into the system. This only has an effect when
%                         'UnitTriangular' is true.
%                         Default: true
%   * 'Rounding'        - Function handle to the function that will perform the rounding operation.
%                         For more information on the interface this function must have, see the
%                         ChopBlas documentation.
%                         Default: @chop
%
% Two configurations for rounding are supported:
%   * One rounding mode.
%   * Separate rounding modes for the multiplication, addition and division
%     operations.
%
% Specifying only roundopts will use the same rounding mode (given by
% roundopts) for all the operations. Individual rounding modes for the
% multiplication, addition and division  operations can be specified in
% the mulopts, addopts and divopts arguments, respectively.
%
% Usage:
%   [xout] = CHTRSV( A, x, ... )
%   [xout] = CHTRSV( A, x, roundopts, ... )
%   [xout] = CHTRSV( A, x, mulopts, addopts, divopts, ... )

% Created by: Ian McInerney
% Created on: August 2, 2022
% SPDX-License-Identifier: BSD-2-Clause

%% Setup the argument parsing
isboolean = @(x) islogical(x) && isscalar(x);
p = inputParser;
p.StructExpand = false;
addOptional( p, 'mulopts', struct([]) );
addOptional( p, 'addopts', struct([]) );
addOptional( p, 'divopts', struct([]) );
addParameter( p, 'LowerTriangular', false, isboolean );
addParameter( p, 'Transpose', false, isboolean );
addParameter( p, 'UnitTriangular', false, isboolean );
addParameter( p, 'Rounding', @chop );
addParameter( p, 'RoundAll', true, isboolean );

parse( p, varargin{:} )

mulopts   = p.Results.mulopts;
addopts   = p.Results.addopts;
divopts   = p.Results.divopts;
trans     = p.Results.Transpose;
lowertri  = p.Results.LowerTriangular;
unitdiag  = p.Results.UnitTriangular;
roundfunc = p.Results.Rounding;
roundall  = p.Results.RoundAll;

% Allow only the first to be specified and have it be used for both
if ( isempty(addopts) || isempty(divopts) ) && ~isempty(mulopts)
    addopts = mulopts;
    divopts = mulopts;
end

% Verify arguments and their sizes
sb = size(b);
sA = size(A);

if sb(2) ~= 1
    error( "chtrsv:bMustBeColumnVector", "The b vector must be a column vector." );
end
if ( ~trans && ( sb(1) ~= sA(2) ) ) || ( trans && ( sb(1) ~= sA(1) ) )
    errmsg = strcat( "A and b must be compatible sizes - ",...
                     "[", num2str(sA(1)), "x", num2str(sA(2)), "] vs. ", ...
                     "[", num2str(sb(1)), "x", num2str(sb(2)), "]." );
    error( "chtrsv:AbMustBeCompatibleSizes", errmsg );
end

% When there is a unit diagonal, we don't need to do any actual division.
% This means we can just leave b alone, but that would mean these entries would
% be used in future computations in the precision given by the user instead
% of divopts (which would be different from all future values of b).
% So provide an option to round them and ensure they are in the divopts precision.
if unitdiag && roundall
    b = roundfunc( b, divopts );
end

[nr, nc] = size(A);

%% Actually compute the result
if trans
    if lowertri
        % Compute the transposed lower-triangular product
        for i=nr:-1:2
            r = 1:1:(i-1);

            if unitdiag && roundall
                b(i) = roundfunc( b(i), divopts );
            elseif ~unitdiag
                b(i) = roundfunc( b(i)/A(i,i), divopts );
            end

            t = roundfunc( b(i).*A(i, r), mulopts );

            b(r) = roundfunc( b(r) - t', addopts );
        end

        fe = 1;
    else
        % Compute the transposed upper-triangular product
        for i=1:1:(nr-1)
            r = (i+1):1:nc;

            if unitdiag && roundall
                b(i) = roundfunc( b(i), divopts );
            elseif ~unitdiag
                b(i) = roundfunc( b(i)/A(i,i), divopts );
            end

            t = roundfunc( b(i).*A(i, r), mulopts );

            b(r) = roundfunc( b(r) - t', addopts );
        end

        fe = nr;
    end
else
    if lowertri
        % Compute the non-transposed lower-triangular product
        for i=1:1:(nc-1)
            r = (i+1):1:nr;

            if unitdiag && roundall
                b(i) = roundfunc( b(i), divopts );
            elseif ~unitdiag
                b(i) = roundfunc( b(i)/A(i,i), divopts );
            end

            t = roundfunc( b(i).*A(r,i), mulopts );

            b(r) = roundfunc( b(r) - t, addopts );
        end

        fe = nr;
    else
        % Compute the non-transposed upper-triangular product
        for i=nc:-1:2
            r = 1:1:(i-1);

            if unitdiag && roundall
                b(i) = roundfunc( b(i), divopts );
            elseif ~unitdiag
                b(i) = roundfunc( b(i)/A(i,i), divopts );
            end

            t = roundfunc( b(i).*A(r, i), mulopts );

            b(r) = roundfunc( b(r) - t, addopts );
        end

        fe = 1;
    end
end

% Final element
if unitdiag && roundall
    b(fe) = roundfunc( b(fe), divopts );
elseif ~unitdiag
    b(fe) = roundfunc( b(fe)/A(fe,fe), divopts );
end

x = b;

end
