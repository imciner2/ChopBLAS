%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% armv7_mvm.m
%
% Simulate a the ARMv7 NEON vector unit floating-point operations.
% This unit only supported single precision values, and flushed
% subnormals to 0 before and after the operations.
%
% Created by: Ian Mcinerney
% Created on: August 19, 2022
% License: BSD-2-Clause
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Seed the RNG for reproducibility
rng(15678);

% Generate some subnormal data
ms = double( realmin( 'single' ) );

N = 10;
A = randn( N, N );
x = logspace( -15, -5, N)';

As = A.*ms


%% Compute using the ARMv7 floating-point options
xarm = armv7_mvm( As, x )


%% Compute with subnormals allowed
roundopts.format = 'fp32';
roundopts.subnormal = 1;

xieee = chgemv( 1.0, As, x, 0.0, [], roundopts )

xarm - xieee
