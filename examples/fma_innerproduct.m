%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% fma_innerproduct.m
%
% Simulate a 16-bit inner product computation done using a fused-multiply-add
% (FMA) computational unit.
%
% Created by: Ian Mcinerney
% Created on: May 20, 2022
% License: BSD-2-Clause
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Create chop options for a 16-bit rounding mode
hopts.format = 'h';


%% Create chop options for the 64-bit "no-rounding" mode
dopts.format = 'd';


%% Generate random 16-bit data
n = 10;

x = chop( rand(n, 1), hopts );
y = chop( rand(n, 1), hopts );



%% Call the inner product
% To simulate FMA, the multiply rounding options are given as double precision
% while the addition rounding options are given as half precision. This simulates
% the absence of rounding in the half-precision inner product computation's
% multiplication and only rounds to the 16-bit format after the addition.

dot = chdot( x, y, dopts, hopts );
