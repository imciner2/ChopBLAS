function xout = chopblas_sum_mat( x, algorithm, roundfunc, opts )
%CHOPBLAS_SUM_MAT Sum all elements of a matrix into a vector
%
% This function will reduce a matrix to a vector by summing all values
% in a row  together using the specified algorithm with operation level
% rounding given by the roundfunc and opts arguments.

% Created by: Ian McInerney
% Created on: August 18, 2022
% SPDX-License-Identifier: BSD-2-Clause

if strcmpi( algorithm, 'recursive' )
    xout = chopblas_recursive_sum_mat( x, roundfunc, opts );
elseif strcmpi( algorithm, 'pairwise' )
    xout = chopblas_pairwise_sum_mat( x, roundfunc, opts );
elseif strcmpi( algorithm, 'increasing' )
    xout = chopblas_sorted_sum_mat( x, 1, roundfunc, opts );
elseif strcmpi( algorithm, 'decreasing' )
    xout = chopblas_sorted_sum_mat( x, 0, roundfunc, opts );
elseif strcmpi( algorithm, 'insertion' )
    xout = chopblas_insertion_sum_mat( x, 1, roundfunc, opts );
else
    errmsg = strcat( "Unknown summation algorithm: ", algorithm );
    error( "chopblas:unknownSummationAlgorithm", errmsg );
end

end
