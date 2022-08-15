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
%   * 'Transpose' - If true, the computation is A'*x instead of A*x
%                   Default: false
%   * 'Rounding'  - Function handle to the function that will perform the rounding operation.
%                   For more information on the interface this function must have, see the
%                   ChopBlas documentation.
%                   Default: @chop
%   * 'Summation' - The algorithm to use when performing the additions.
%                   Supported algorithms: 'recursive', 'pairwise', 'increasing', 'decreasing'
%                   Default: 'recursive'
%   * 'BlockSize' - The number of rows of A to process in each internal iteration.
%                   Default: All rows of A.
%
% The order of operations for this function are as follows:
%   1) Populate xout with the y vector
%      If beta==1, then no multiplication/rounding is done
%      if beta==0 or y is empty, xout is loaded with 0.
%   2) Scale x by alpha (if alpha==1, no multiplication/rounding is done)
%   3) Compute A*x or A'*x with the scaled x by row, summing into xout using the
%      summation method selected.
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
% The BlockSize option defines a compute vs memory tradeoff for the function,
% with more rows requiring more working memory but fewer calls into the
% Rounding function, and fewer rows requiring more calls to the Rounding
% function but less working memory. For instance, processing all rows makes
% another copy of A for temporary use and only calls the Rounding function at most
% n+1 times (for a recursive summation method), while processing one row at a time
% only creates a new vector for temporary use but calls the rounding function upto
% n^2 + 1 times (for a recursive summation method).
%
% Usage:
%   [xout] = CHGEMV( alpha, A, x, beta, y, ... )
%   [xout] = CHGEMV( alpha, A, x, beta, y, roundopts, ... )
%   [xout] = CHGEMV( alpha, A, x, beta, y, mulopts, addopts, ... )

% Created by: Ian McInerney
% Created on: June 16, 2022
% SPDX-License-Identifier: BSD-2-Clause


%% Setup the argument parsing
isboolean = @(x) islogical(x) && isscalar(x);
p = inputParser;
p.StructExpand = false;
addOptional( p, 'mulopts', struct([]) );
addOptional( p, 'addopts', struct([]) );
addParameter( p, 'Transpose', false, isboolean );
addParameter( p, 'Rounding', @chop );
addParameter( p, 'Summation', 'recursive' );
addParameter( p, 'BlockSize', size(A,1) );

parse( p, varargin{:} )

mulopts   = p.Results.mulopts;
addopts   = p.Results.addopts;
trans     = p.Results.Transpose;
roundfunc = p.Results.Rounding;
algorithm = p.Results.Summation;
blocksize = p.Results.BlockSize;

% Allow only the first to be specified and have it be used for both
if isempty(addopts) && ~isempty(mulopts)
    addopts = mulopts;
end

% Verify arguments and their sizes
sx = size(x);
sy = size(y);
sA = size(A);

if ~isscalar( alpha )
    error( "chgemv:AlphaMustBeScalar", "alpha must be a scalar." );
end
if ~isscalar( beta )
    error( "chgemv:BetaMustBeScalar", "beta must be a scalar." );
end
if ( ( sx(2) ~= 1 ) && ~isempty(x) ) || ( ( sy(2) ~= 1 ) && ~isempty(y) )
    error( "chgemv:xyMustBeColumnVectors", "The x and y vectors must be column vectors." );
end
if ~( all( sx == sy ) ) && ~isempty(x) && ~isempty(y)
    error( "chgemv:xyMustBeSameSize", "The x and y vectors must be the same size." );
end
if ( ~trans && ( sx(1) ~= sA(2) ) ) || ( trans && ( sx(1) ~= sA(1) ) )
    errmsg = strcat( "A and x must be compatible sizes - ", ...
                     "[", num2str(sA(1)), "x", num2str(sA(2)), "] vs. ", ...
                     "[", num2str(sx(1)), "x", num2str(sx(2)), "]." );
    error( "chgemv:AxMustBeCompatibleSizes", errmsg );
end
if blocksize > sA(1)
    errmsg = strcat( "Block size too large, resetting to ", num2str(sA(1)) );
    warning( "chgemv:BlockSizeTooLarge", errmsg );
    blocksize = sA(1);
end
if blocksize < 1
    warning( "chgemv:BlockSizeTooSmall", "Block size too small, resetting to 1" );
    blocksize = 1;
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
    matind = @(i) A(:,i)';
else
    % Matrix indexing needed to compute the non-transposed product
    matind = @(i) A(i,:);
end

lx = length(x);
for i=1:blocksize:lx
    inds = i:1:min(i+blocksize-1, lx);

    % matind() will return a column vector of the matrix elements in the proper position
    % for the final add (based on the transpose option)
    t = [xout(inds), roundfunc( matind(inds).*x', mulopts )];

    if strcmpi( algorithm, 'recursive' )
        xout(inds) = chopblas_recursive_sum_mat( t, roundfunc, addopts );
    elseif strcmpi( algorithm, 'pairwise' )
        xout(inds) = chopblas_pairwise_sum_mat( t, roundfunc, addopts );
    elseif strcmpi( algorithm, 'increasing' )
        xout = chopblas_sorted_sum_mat( t, 1, roundfunc, addopts );
    elseif strcmpi( algorithm, 'decreasing' )
        xout = chopblas_sorted_sum_mat( t, 0, roundfunc, addopts );
    else
        errmsg = strcat( "Unknown summation algorithm: ", algorithm );
        error( "chgemv:unknownSummationAlgorithm", errmsg );
    end
end

end
