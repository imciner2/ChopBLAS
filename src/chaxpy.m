function [xout] = chaxpy( alpha, x, y, mulopts, addopts )
%CHAXPY Add the scaled vector x to the vector y
%
% Scale each element of the vector x by the scale factor alpha then add
% the result to the elements in vector y with rounding using chop after
% each operations.
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
%   [xout] = chaxpy( alpha, x, y )
%   [xout] = chaxpy( alpha, x, opts )
%   [xout] = chaxpy( alpha, x, mulopts, addopts )
%
% Created by: Ian McInerney
% Created on: June 15, 2022
% License: BSD-2-Clause

if nargin < 4
    mulopts = [];
    addopts = [];
elseif nargin < 5
    addopts = mulopts;
end

xout = chop( alpha.*x, mulopts );
xout = chop( xout+y, addopts );

end
