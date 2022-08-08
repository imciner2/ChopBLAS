classdef kernel_recursive_sum < matlab.unittest.TestCase
%KERNEL_RECURSIVE_SUM Test for the recursive summation compute kernel

    properties
        sopts
        hopts
        dopts

        rf

        xodd
        xeven

        tol
    end

    methods(TestMethodSetup)
        function setup_test(testCase)
            testCase.sopts.format = 's';
            testCase.hopts.format = 'h';
            testCase.dopts.format = 'd';

            % Choose the rounding function to use
            if strcmpi( getenv('CHOPBLAS_ROUND_FUNC'), 'cpfloat' )
                testCase.rf = @cpfloat;
            else
                % Default to chop
                testCase.rf = @chop;
            end

            % These have been designed to give different results for the different summation algorithms
            testCase.xodd = [1.0005; 2.3458954; 3.242456; 4; 5.25678];
            testCase.xeven = [1.0005; 2.3458954; 3.242456; 4];

            % Set a tolerance for all the tests
            testCase.tol = 1e-4;
        end
    end

    methods(Test)
        % Test with a vector with an even number of elements
        function even_length(testCase)
            % Test with global rounding options
            % In double precision, this is the same value as the normal sum
            testCase.rf( [], testCase.dopts );

            res = sum( testCase.xeven );
            z   = chopblas_recursive_sum( testCase.xeven, testCase.rf, struct([]) );
            testCase.verifyEqual( z, res );

            % Test with a half-precision rounding mode
            res = testCase.xeven(1);
            for i=2:1:length(testCase.xeven)
                res = double( half( res + testCase.xeven(i) ) );
            end

            z = chopblas_recursive_sum( testCase.xeven, testCase.rf, testCase.hopts );
            testCase.verifyEqual( z, res );
        end

        % Test with a vector with an odd number of elements
        function odd_length(testCase)
            % Test with global rounding options
            % In double precision, this is the same value as the normal sum
            testCase.rf( [], testCase.dopts );

            res = sum( testCase.xodd );
            z   = chopblas_recursive_sum( testCase.xodd, testCase.rf, struct([]) );
            testCase.verifyEqual( z, res );

            % Test with a half-precision rounding mode
            res = testCase.xodd(1);
            for i=2:1:length(testCase.xodd)
                res = double( half( res + testCase.xodd(i) ) );
            end

            z = chopblas_recursive_sum( testCase.xodd, testCase.rf, testCase.hopts );
            testCase.verifyEqual( z, res );
        end
    end
end
