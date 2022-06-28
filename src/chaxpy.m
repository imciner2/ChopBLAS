function [xout] = chaxpy( alpha, x, y, varargin )
%CHAXPY Add the scaled vector x to the vector y with operation-level rounding
%
% Scale each element of the vector x by the scale factor alpha then add
% the result to the elements in vector y with rounding after each operation.
%
% This function supports the following optional name-value arguments:
%   * 'Rounding' - Function handle to the function that will perform the rounding operation.
%                  For more information on the interface 'roundfunc' must present, see the
%                  ChopBlas documentation.
%                  Default: @chop
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
%   [xout] = CHAXPY( alpha, x, y, ... )
%   [xout] = CHAXPY( alpha, x, opts, ... )
%   [xout] = CHAXPY( alpha, x, mulopts, addopts, ... )

% Created by: Ian McInerney
% Created on: June 15, 2022
% License: BSD-2-Clause

%% Setup the argument parsing
p = inputParser;
p.StructExpand = false;
addOptional( p, 'mulopts', struct([]) );
addOptional( p, 'addopts', struct([]) );
addParameter( p, 'Rounding', @chop );

parse( p, varargin{:} )

mulopts   = p.Results.mulopts;
addopts   = p.Results.addopts;
roundfunc = p.Results.Rounding;

% Allow only the first to be specified and have it be used for both
if isempty(addopts) && ~isempty(mulopts)
    addopts = mulopts;
end

xout = chop( alpha.*x, mulopts );
xout = chop( xout+y, addopts );

end
