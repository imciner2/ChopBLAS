classdef accumulate_sorted < matlab.unittest.TestCase
%ACCUMULATE_SORTED Test for the recursive sorted accumulation compute kernel

    properties
        sopts
        hopts
        dopts

        rf

        xoddsorted
        xevensorted

        xodd
        xeven

        Xoddsorted
        Xevensorted

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
            testCase.xoddsorted  = [1.0005; 2.3458954; 3.242456; 4; 5.25678];
            testCase.xevensorted = [1.0005; 2.3458954; 3.242456; 4];

            testCase.xodd  = testCase.xoddsorted( randperm( length(testCase.xoddsorted) ) );
            testCase.xeven = testCase.xevensorted( randperm( length(testCase.xevensorted) ) );

            % Matrix examples
            testCase.Xoddsorted  = repmat( testCase.xoddsorted', 5, 1 );
            testCase.Xevensorted = repmat( testCase.xevensorted', 4, 1 );

            % Randomly permute the order of the elements in each row of the matrices
            % Each row should be ordered differently so that we can ensure they get
            % sorted properly
            [nr, nc] = size( testCase.Xoddsorted );
            for i = 1:1:nr
                testCase.Xodd(i,:) = testCase.Xoddsorted( i, randperm( nc ) );
            end

            [nr, nc] = size( testCase.Xevensorted );
            for i = 1:1:nr
                testCase.Xeven(i,:) = testCase.Xevensorted( i, randperm( nc ) );
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

            % Increasing order
            z = chaccum_sorted( testCase.xeven, testCase.rf, struct([]), 'ascend' );
            testCase.verifyEqual( z, res, 'AbsTol', testCase.tol );

            % Decreasing order
            z = chaccum_sorted( testCase.xeven, testCase.rf, struct([]), 'descend' );
            testCase.verifyEqual( z, res, 'AbsTol', testCase.tol );

            % Test with a half-precision rounding mode
            res = testCase.xevensorted(1);
            for i=2:1:length(testCase.xevensorted)
                res = double( half( res + testCase.xevensorted(i) ) );
            end

            z = chaccum_sorted( testCase.xeven, testCase.rf, testCase.hopts, 'ascend' );
            testCase.verifyEqual( z, res );
        end

        % Test with a row vector with an even number of elements
        function even_length_row_vec(testCase)
            % Test with global rounding options
            % In double precision, this is the same value as the normal sum
            testCase.rf( [], testCase.dopts );

            res = sum( testCase.xeven );

            % Increasing order
            z = chaccum_sorted( testCase.xeven', testCase.rf, struct([]), 'ascend' );
            testCase.verifyEqual( z, res, 'AbsTol', testCase.tol );

            % Decreasing order
            z = chaccum_sorted( testCase.xeven', testCase.rf, struct([]), 'descend' );
            testCase.verifyEqual( z, res, 'AbsTol', testCase.tol );

            % Test with a half-precision rounding mode
            res = testCase.xevensorted(1);
            for i=2:1:length(testCase.xevensorted)
                res = double( half( res + testCase.xevensorted(i) ) );
            end

            z = chaccum_sorted( testCase.xeven', testCase.rf, testCase.hopts, 'ascend' );
            testCase.verifyEqual( z, res );
        end

        % Test with a matrix with an even number of elements in each row
        function even_length_mat(testCase)
            % Test with global rounding options
            % In double precision, this is the same value as the normal sum
            testCase.rf( [], testCase.dopts );

            res = sum( testCase.Xeven, 2 );

            % Increasing order
            z = chaccum_sorted( testCase.Xeven, testCase.rf, struct([]), 'ascend' );
            testCase.verifyEqual( z, res, 'AbsTol', testCase.tol );

            % Decreasing order
            z = chaccum_sorted( testCase.Xeven, testCase.rf, struct([]), 'descend' );
            testCase.verifyEqual( z, res, 'AbsTol', testCase.tol );

            % Test with a half-precision rounding mode
            res = testCase.Xevensorted(:,1);
            for i=2:1:size(testCase.Xevensorted,2)
                res = double( half( res + testCase.Xevensorted(:,i) ) );
            end

            z = chaccum_sorted( testCase.Xeven, testCase.rf, testCase.hopts, 'ascend' );
            testCase.verifyEqual( z, res );
        end

        % Test with a column vector with an odd number of elements
        function odd_length_column_vec(testCase)
            % Test with global rounding options
            % In double precision, this is the same value as the normal sum
            testCase.rf( [], testCase.dopts );

            res = sum( testCase.xodd );

            % Increasing order
            z = chaccum_sorted( testCase.xodd, testCase.rf, struct([]), 'ascend' );
            testCase.verifyEqual( z, res, 'AbsTol', testCase.tol );

            % Decreasing order
            z = chaccum_sorted( testCase.xodd, testCase.rf, struct([]), 'descend' );
            testCase.verifyEqual( z, res, 'AbsTol', testCase.tol );

            % Test with a half-precision rounding mode
            res = testCase.xoddsorted(1);
            for i=2:1:length(testCase.xoddsorted)
                res = double( half( res + testCase.xoddsorted(i) ) );
            end

            z = chaccum_sorted( testCase.xodd, testCase.rf, testCase.hopts, 'ascend' );
            testCase.verifyEqual( z, res );
        end

        % Test with a row vector with an odd number of elements
        function odd_length_row_vec(testCase)
            % Test with global rounding options
            % In double precision, this is the same value as the normal sum
            testCase.rf( [], testCase.dopts );

            res = sum( testCase.xodd );

            % Increasing order
            z = chaccum_sorted( testCase.xodd', testCase.rf, struct([]), 'ascend' );
            testCase.verifyEqual( z, res, 'AbsTol', testCase.tol );

            % Decreasing order
            z = chaccum_sorted( testCase.xodd', testCase.rf, struct([]), 'descend' );
            testCase.verifyEqual( z, res, 'AbsTol', testCase.tol );

            % Test with a half-precision rounding mode
            res = testCase.xoddsorted(1);
            for i=2:1:length(testCase.xoddsorted)
                res = double( half( res + testCase.xoddsorted(i) ) );
            end

            z = chaccum_sorted( testCase.xodd', testCase.rf, testCase.hopts, 'ascend' );
            testCase.verifyEqual( z, res );
        end

        % Test with a matrix with an odd number of elements in each row
        function odd_length_mat(testCase)
            % Test with global rounding options
            % In double precision, this is the same value as the normal sum
            testCase.rf( [], testCase.dopts );

            res = sum( testCase.Xodd, 2 );

            % Increasing order
            z = chaccum_sorted( testCase.Xodd, testCase.rf, struct([]), 'ascend' );
            testCase.verifyEqual( z, res, 'AbsTol', testCase.tol );

            % Decreasing order
            z = chaccum_sorted( testCase.Xodd, testCase.rf, struct([]), 'descend' );
            testCase.verifyEqual( z, res, 'AbsTol', testCase.tol );

            % Test with a half-precision rounding mode
            res = testCase.Xoddsorted(:,1);
            for i=2:1:size(testCase.Xoddsorted, 2)
                res = double( half( res + testCase.Xoddsorted(:,i) ) );
            end

            z = chaccum_sorted( testCase.Xodd, testCase.rf, testCase.hopts, 'ascend' );
            testCase.verifyEqual( z, res );
        end
    end
end
