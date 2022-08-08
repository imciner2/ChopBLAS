classdef kernel_pairwise_sum < matlab.unittest.TestCase
%KERNEL_PAIRWISE_SUM Test for the pairwise summation compute kernel

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

            testCase.xodd = [1.0005; 2.34; 3.24; 4; 5.25678];
            testCase.xeven = [1.0005; 2.34; 3.24; 4];

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
            z   = chopblas_pairwise_sum( testCase.xeven, testCase.rf, struct([]) );
            testCase.verifyEqual( z, res );

            % Test with a half-precision rounding mode
            x = testCase.xeven;
            i1  = double( half( x(1) + x(3) ) );
            i2  = double( half( x(2) + x(4) ) );
            res = double( half( i1 + i2 ) );

            z = chopblas_recursive_sum( testCase.xeven, testCase.rf, testCase.hopts );
            testCase.verifyEqual( z, res );
        end

        % Test with a vector with an odd number of elements
        function odd_length(testCase)
            % Test with global rounding options
            % In double precision, this is the same value as the normal sum
            testCase.rf( [], testCase.dopts );

            res = sum( testCase.xodd );
            z   = chopblas_pairwise_sum( testCase.xodd, testCase.rf, struct([]) );
            testCase.verifyEqual( z, res );

            % Test with a half-precision rounding mode
            x = testCase.xodd;

            i1  = double( half( x(1) + x(3) ) );
            i2  = double( half( x(2) + x(4) ) );
            i3  = double( half( i1 + i2 ) );
            res = double( half( i3 + x(5) ) );

            z = chopblas_pairwise_sum( testCase.xodd, testCase.rf, testCase.hopts );
            testCase.verifyEqual( z, res );
        end
    end
end
