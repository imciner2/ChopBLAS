function [nrm] = chnrm2( x, mulopts, accumopts, sqrtopts )
%CHNRM2 Compute the 2-norm of the vector x
%
% Compute the 2-norm of the vector x.
%
% Two configurations for rounding are supported:
%   * One rounding mode.
%   * Separate rounding modes for the multiplication, accumulation
%     and square-root operations.
%
% Specifying only opts will use the same rounding mode (given by opts)
% for all three operations.
% Individual rounding modes for the multiplication, accumulation and
% square-root operations can be specified in the mulopts, accumopts,
% and sqrtopts arguments, respectively.
%
% Usage:
%   [nrm] = chnrm2( x )
%   [nrm] = chnrm2( x, opts )
%   [nrm] = chnrm2( x, mulopts, accumopts, sqrtopts )
%
% Created by: Ian McInerney
% Created on: May 20, 2022
% License: BSD-2-Clause

if nargin < 2
    mulopts = [];
    accumopts = [];
    sqrtopts = [];
elseif nargin < 4
    accumopts = mulopts;
    sqrtopts = mulopts;
end

pp = chop( x.*x, mulopts );

dot = pp(1);
for i=2:length(pp)
    dot = chop( dot + pp(i), accumopts );
end

nrm = chop( sqrt( dot ), sqrtopts );

end
