classdef kernel_pairwise_sum < matlab.unittest.TestCase
%KERNEL_PAIRWISE_SUM Test for the pairwise summation compute kernel

    properties
        sopts
        hopts
        dopts

        rf

        xodd
        xeven
        xevenodd
        xevenodd14

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

            testCase.xodd  = [1.0005; 2.34; 3.24; 4; 5.25678];
            testCase.xeven = [1.0005; 2.34; 3.24; 4];

            % Examples that go from an even to an odd length during summation
            testCase.xevenodd   = [0.05; 1.0005; 2.34; 3.24; 4; 5.25678];
            testCase.xevenodd14 = [testCase.xevenodd; testCase.xevenodd; 2.2568; 1.2348];

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
            z   = chopblas_pairwise_sum_vec( testCase.xeven, testCase.rf, struct([]) );
            testCase.verifyEqual( z, res );

            % Test with a half-precision rounding mode
            x = testCase.xeven;
            i1  = double( half( x(1) + x(3) ) );
            i2  = double( half( x(2) + x(4) ) );
            res = double( half( i1 + i2 ) );

            z = chopblas_recursive_sum_vec( testCase.xeven, testCase.rf, testCase.hopts );
            testCase.verifyEqual( z, res );
        end

        % Test with a vector with an odd number of elements
        function odd_length(testCase)
            % Test with global rounding options
            % In double precision, this is the same value as the normal sum
            testCase.rf( [], testCase.dopts );

            res = sum( testCase.xodd );
            z   = chopblas_pairwise_sum_vec( testCase.xodd, testCase.rf, struct([]) );
            testCase.verifyEqual( z, res );

            % Test with a half-precision rounding mode
            x = testCase.xodd;

            i1  = double( half( x(1) + x(3) ) );
            i2  = double( half( x(2) + x(4) ) );
            i3  = double( half( i1 + i2 ) );
            res = double( half( i3 + x(5) ) );

            z = chopblas_pairwise_sum_vec( testCase.xodd, testCase.rf, testCase.hopts );
            testCase.verifyEqual( z, res );
        end

        % Test with a vector with an even number of elements to start, but in the process
        % gets an odd number of elements
        function evenodd_length(testCase)
            % Test with global rounding options
            % In double precision, this is the same value as the normal sum (within tolerance)
            testCase.rf( [], testCase.dopts );

            res = sum( testCase.xevenodd );
            z   = chopblas_pairwise_sum_vec( testCase.xevenodd, testCase.rf, struct([]) );
            testCase.verifyEqual( z, res, 'AbsTol', testCase.tol );

            res = sum( testCase.xevenodd14 );
            z   = chopblas_pairwise_sum_vec( testCase.xevenodd14, testCase.rf, struct([]) );
            testCase.verifyEqual( z, res, 'AbsTol', testCase.tol );

            %% Test with a half-precision rounding mode
            x = testCase.xevenodd;

            % First batch goes from 6 values to 3
            i1  = double( half( x(1) + x(2) ) );
            i2  = double( half( x(3) + x(4) ) );
            i3  = double( half( x(5) + x(6) ) );

            ii1  = double( half( i1 + i2 ) );
            res = double( half( ii1 + i3 ) );

            z = chopblas_pairwise_sum_vec( testCase.xevenodd, testCase.rf, testCase.hopts );
            testCase.verifyEqual( z, res );

            x = testCase.xevenodd14;

            % First batch goes from 14 values to 7
            i1  = double( half( x(1)  + x(2) ) );
            i2  = double( half( x(3)  + x(4) ) );
            i3  = double( half( x(5)  + x(6) ) );
            i4  = double( half( x(7)  + x(8) ) );
            i5  = double( half( x(9)  + x(10) ) );
            i6  = double( half( x(11) + x(12) ) );
            i7  = double( half( x(13) + x(14) ) );

            % Next goes from 7 to 4
            ii1 = double( half( i1 + i2 ) );
            ii2 = double( half( i3 + i4 ) );
            ii3 = double( half( i5 + i6 ) );
            ii4 = i7;

            % Next goes from 4 to 2
            iii1 = double( half( ii1 + ii2 ) );
            iii2 = double( half( ii3 + ii4 ) );

            % Final result
            res = double( half( iii1 + iii2 ) );

            z = chopblas_pairwise_sum_vec( testCase.xevenodd14, testCase.rf, testCase.hopts );
            testCase.verifyEqual( z, res );
        end
    end
end
