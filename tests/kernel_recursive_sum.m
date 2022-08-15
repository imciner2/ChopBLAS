classdef kernel_recursive_sum < matlab.unittest.TestCase
%KERNEL_RECURSIVE_SUM Test for the recursive summation compute kernel

    properties
        sopts
        hopts
        dopts

        rf

        xodd
        xeven

        Xodd
        Xeven

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

            % Matrix examples
            testCase.Xodd  = repmat( testCase.xodd', 5, 1 );
            testCase.Xeven = repmat( testCase.xeven', 4, 1 );

            % Set a tolerance for all the tests
            testCase.tol = 1e-4;
        end
    end

    methods(Test)
        % Test with a vector with an even number of elements
        function even_length_vec(testCase)
            % Test with global rounding options
            % In double precision, this is the same value as the normal sum
            testCase.rf( [], testCase.dopts );

            res = sum( testCase.xeven );
            z   = chopblas_recursive_sum_vec( testCase.xeven, testCase.rf, struct([]) );
            testCase.verifyEqual( z, res );

            % Test with a half-precision rounding mode
            res = testCase.xeven(1);
            for i=2:1:length(testCase.xeven)
                res = double( half( res + testCase.xeven(i) ) );
            end

            z = chopblas_recursive_sum_vec( testCase.xeven, testCase.rf, testCase.hopts );
            testCase.verifyEqual( z, res );
        end

        % Test with a matrix with an even number of elements in each row
        function even_length_mat(testCase)
            % Test with global rounding options
            % In double precision, this is the same value as the normal sum
            testCase.rf( [], testCase.dopts );

            res = sum( testCase.Xeven, 2 );
            z   = chopblas_recursive_sum_mat( testCase.Xeven, testCase.rf, struct([]) );
            testCase.verifyEqual( z, res );

            % Test with a half-precision rounding mode
            res = testCase.Xeven(:,1);
            for i=2:1:size(testCase.Xeven,2)
                res = double( half( res + testCase.Xeven(:,i) ) );
            end

            z = chopblas_recursive_sum_mat( testCase.Xeven, testCase.rf, testCase.hopts );
            testCase.verifyEqual( z, res );
        end

        % Test with a vector with an odd number of elements
        function odd_length_vec(testCase)
            % Test with global rounding options
            % In double precision, this is the same value as the normal sum
            testCase.rf( [], testCase.dopts );

            res = sum( testCase.xodd );
            z   = chopblas_recursive_sum_vec( testCase.xodd, testCase.rf, struct([]) );
            testCase.verifyEqual( z, res );

            % Test with a half-precision rounding mode
            res = testCase.xodd(1);
            for i=2:1:length(testCase.xodd)
                res = double( half( res + testCase.xodd(i) ) );
            end

            z = chopblas_recursive_sum_vec( testCase.xodd, testCase.rf, testCase.hopts );
            testCase.verifyEqual( z, res );
        end

        % Test with a matrix with an odd number of elements in each row
        function odd_length_mat(testCase)
            % Test with global rounding options
            % In double precision, this is the same value as the normal sum
            testCase.rf( [], testCase.dopts );

            res = sum( testCase.Xodd, 2 );
            z   = chopblas_recursive_sum_mat( testCase.Xodd, testCase.rf, struct([]) );
            testCase.verifyEqual( z, res );

            % Test with a half-precision rounding mode
            res = testCase.Xodd(:,1);
            for i=2:1:size(testCase.Xodd, 2)
                res = double( half( res + testCase.Xodd(:,i) ) );
            end

            z = chopblas_recursive_sum_mat( testCase.Xodd, testCase.rf, testCase.hopts );
            testCase.verifyEqual( z, res );
        end
    end
end
