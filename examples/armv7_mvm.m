function [xout] = armv7_mvm(A, x)
    roundopts.format = 'fp32';     % Use single precision
    roundopts.subnormal = 0;       % Flush subnormals to zero

    % Flush subnormals to 0 in the inputs
    A = chop( A, roundopts );
    x = chop( x, roundopts );

    xout = chgemv( 1.0, A, x, 0.0, [], roundopts );
end
