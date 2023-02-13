function [s] = chaccum_insertion_sort( x, roundfunc, opts, dir )
%CHACCUM_INSERTION_SORT Accumulate into a column vector using the insertion sort algorithm
%
% Reduce a matrix to a column vector by adding all the entries in a row together in
% a recursive fasion, with rounding after each operation. First the rows of the matrix
% x are sorted in the specified order, then the first two elements are added together
% and the result is placed in the proper sorted position, and the process repeats.
%
% 'dir' can either be 'ascend' or 'descend' to sort by ascending or descending values, respectively.

% Created by: Ian McInerney
% Created on: August 18, 2022
% SPDX-License-Identifier: BSD-2-Clause


sfunc = @(x) sort(x, 2, dir, 'ComparisonMethod', 'abs');

% It is much easier to work with a row vector (since we work with row matrices), and since
% we already need a copy of the vector, just transpose a column vector to a row vector for
% QOL.
if iscolumn( x )
    s = sfunc( x' );
else
    s = sfunc( x );
end

while size(s,2) > 1
    t = roundfunc( s(:,1) + s(:,2), opts );
    s = s(:,2:end);
    s(:,1) = t;
    s = sfunc(s);
end

end
