classdef accumulate_insertion_sort < matlab.unittest.TestCase
%accumulate_insertion_sort Test for the insertion sort accumulation compute kernel

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
            testCase.xodd  = [0.5; 1.0; 0.5; 0.75; 0.5];
            testCase.xeven = [0.5; 1.0; 0.5; 0.75; 0.5; 0.5];

            testCase.Xodd  = repmat( testCase.xodd', 5, 1 );
            testCase.Xeven = repmat( testCase.xeven', 6, 1 );

            % Randomly permute the order of the elements in each row of the matrices
            % Each row should be ordered differently so that we can ensure they get
            % sorted properly
            [nr, nc] = size( testCase.Xodd );
            for i = 1:1:nr
                testCase.Xodd(i,:)  = testCase.Xodd( i, randperm( nc ) );
            end

            [nr, nc] = size( testCase.Xeven );
            for i = 1:1:nr
                testCase.Xeven(i,:)  = testCase.Xeven( i, randperm( nc ) );
            end

            % Set a tolerance for all the tests
            testCase.tol = 1e-4;
        end
    end

    methods(Test)
        % Test with a column vector with an even number of elements
        function even_length_column_vec(testCase)
            % Test with global rounding options
            % In double precision, this is the same value as the normal sum
            testCase.rf( [], testCase.dopts );

            res = sum( testCase.xeven );
            z   = chaccum_insertion_sort( testCase.xeven, testCase.rf, struct([]), 'ascend' );
            testCase.verifyEqual( z, res );

            % Test with round to nearest (tie to up)
            z = chaccum_insertion_sort( testCase.xeven, @(x, s) round(x), struct([]), 'ascend' );
            testCase.verifyEqual( z, 4.0 );
        end

        % Test with a row vector with an even number of elements
        function even_length_row_vec(testCase)
            % Test with global rounding options
            % In double precision, this is the same value as the normal sum
            testCase.rf( [], testCase.dopts );

            res = sum( testCase.xeven );
            z   = chaccum_insertion_sort( testCase.xeven', testCase.rf, struct([]), 'ascend' );
            testCase.verifyEqual( z, res );

            % Test with round to nearest (tie to up)
            z = chaccum_insertion_sort( testCase.xeven', @(x, s) round(x), struct([]), 'ascend' );
            testCase.verifyEqual( z, 4.0 );
        end

        % Test with a matrix with an even number of elements
        function even_length_mat(testCase)
            % Test with global rounding options
            % In double precision, this is the same value as the normal sum
            testCase.rf( [], testCase.dopts );

            res = sum( testCase.Xeven, 2 );
            z   = chaccum_insertion_sort( testCase.Xeven, testCase.rf, struct([]), 'ascend' );
            testCase.verifyEqual( z, res );

            % Test with round to nearest (tie to up)
            z = chaccum_insertion_sort( testCase.Xeven, @(x, s) round(x), struct([]), 'ascend' );
            testCase.verifyEqual( z, 4.0*ones( size( testCase.Xeven, 1 ), 1 ) );
        end

        % Test with a column vector with an odd number of elements
        function odd_length_column_vec(testCase)
            % Test with global rounding options
            % In double precision, this is the same value as the normal sum
            testCase.rf( [], testCase.dopts );

            res = sum( testCase.xodd );
            z   = chaccum_insertion_sort( testCase.xodd, testCase.rf, struct([]), 'ascend' );
            testCase.verifyEqual( z, res );

            % Test with round to nearest (tie to up)
            z = chaccum_insertion_sort( testCase.xodd, @(x, s) round(x), struct([]), 'ascend' );
            testCase.verifyEqual( z, 3.0 );
        end

        % Test with a row vector with an odd number of elements
        function odd_length_row_vec(testCase)
            % Test with global rounding options
            % In double precision, this is the same value as the normal sum
            testCase.rf( [], testCase.dopts );

            res = sum( testCase.xodd );
            z   = chaccum_insertion_sort( testCase.xodd', testCase.rf, struct([]), 'ascend' );
            testCase.verifyEqual( z, res );

            % Test with round to nearest (tie to up)
            z = chaccum_insertion_sort( testCase.xodd', @(x, s) round(x), struct([]), 'ascend' );
            testCase.verifyEqual( z, 3.0 );
        end

        % Test with a matrix with an odd number of elements
        function odd_length_mat(testCase)
            % Test with global rounding options
            % In double precision, this is the same value as the normal sum
            testCase.rf( [], testCase.dopts );

            res = sum( testCase.Xodd, 2 );
            z   = chaccum_insertion_sort( testCase.Xodd, testCase.rf, struct([]), 'ascend' );
            testCase.verifyEqual( z, res );

            % Test with round to nearest (tie to up)
            z = chaccum_insertion_sort( testCase.Xodd, @(x, s) round(x), struct([]), 'ascend' );
            testCase.verifyEqual( z, 3.0*ones( size( testCase.Xodd, 1 ), 1 ) );
        end
    end
end
