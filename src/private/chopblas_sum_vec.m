function xout = chopblas_sum_vec( x, algorithm, roundfunc, opts )
%CHOPBLAS_SUM_VEC Sum all elements of a vector into a scalar value
%
% This function will reduce a vector to a scalar by summing all values
% together using the specified algorithm with operation level rounding
% given by the roundfunc and opts arguments.

% Created by: Ian McInerney
% Created on: August 15, 2022
% SPDX-License-Identifier: BSD-2-Clause

if strcmpi( algorithm, 'recursive' )
    xout = chopblas_recursive_sum_vec( x, roundfunc, opts );
elseif strcmpi( algorithm, 'pairwise' )
    xout = chopblas_pairwise_sum_vec( x, roundfunc, opts );
elseif strcmpi( algorithm, 'increasing' )
    xout = chopblas_sorted_sum_vec( x, 1, roundfunc, opts );
elseif strcmpi( algorithm, 'decreasing' )
    xout = chopblas_sorted_sum_vec( x, 0, roundfunc, opts );
else
    errmsg = strcat( "Unknown summation algorithm: ", algorithm );
    error( "chopblas:unknownSummationAlgorithm", errmsg );
end

end
