function [dot] = fp16_fma_dot(x, y)
    hopts.format = 'fp16'; % Addition rounding
    dopts.format = 'fp64'; % Multiplication "no-rounding"
    dot = chdot( x, y, dopts, hopts );
end
