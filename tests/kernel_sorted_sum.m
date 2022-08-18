classdef kernel_sorted_sum < matlab.unittest.TestCase
%KERNEL_SORTED_SUM Test for the recursive sorted summation compute kernel

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
            % Ensure the kernels are on the path (they are private to the source code)
            % We can't use the PathFixture here, because it is private so MATLAB prevents
            % adding it to the path. Instead change the folder to the private folder with
            % the kernels in it.
            import matlab.unittest.fixtures.CurrentFolderFixture
            testCase.applyFixture( CurrentFolderFixture( ["../src/private"] ) );

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
        % Test with a vector with an even number of elements
        function even_length_vec(testCase)
            % Test with global rounding options
            % In double precision, this is the same value as the normal sum
            testCase.rf( [], testCase.dopts );

            res = sum( testCase.xeven );

            % Increasing order
            z = chopblas_sorted_sum_vec( testCase.xeven, 1, testCase.rf, struct([]) );
            testCase.verifyEqual( z, res, 'AbsTol', testCase.tol );

            % Decreasing order
            z = chopblas_sorted_sum_vec( testCase.xeven, 0, testCase.rf, struct([]) );
            testCase.verifyEqual( z, res, 'AbsTol', testCase.tol );

            % Test with a half-precision rounding mode
            res = testCase.xevensorted(1);
            for i=2:1:length(testCase.xevensorted)
                res = double( half( res + testCase.xevensorted(i) ) );
            end

            z = chopblas_sorted_sum_vec( testCase.xeven, 1, testCase.rf, testCase.hopts );
            testCase.verifyEqual( z, res );
        end

        % Test with a matrix with an even number of elements in each row
        function even_length_mat(testCase)
            % Test with global rounding options
            % In double precision, this is the same value as the normal sum
            testCase.rf( [], testCase.dopts );

            res = sum( testCase.Xeven, 2 );

            % Increasing order
            z = chopblas_sorted_sum_mat( testCase.Xeven, 1, testCase.rf, struct([]) );
            testCase.verifyEqual( z, res, 'AbsTol', testCase.tol );

            % Decreasing order
            z = chopblas_sorted_sum_mat( testCase.Xeven, 0, testCase.rf, struct([]) );
            testCase.verifyEqual( z, res, 'AbsTol', testCase.tol );

            % Test with a half-precision rounding mode
            res = testCase.Xevensorted(:,1);
            for i=2:1:size(testCase.Xevensorted,2)
                res = double( half( res + testCase.Xevensorted(:,i) ) );
            end

            z = chopblas_sorted_sum_mat( testCase.Xeven, 1, testCase.rf, testCase.hopts );
            testCase.verifyEqual( z, res );
        end

        % Test with a vector with an odd number of elements
        function odd_length_vec(testCase)
            % Test with global rounding options
            % In double precision, this is the same value as the normal sum
            testCase.rf( [], testCase.dopts );

            res = sum( testCase.xodd );

            % Increasing order
            z = chopblas_sorted_sum_vec( testCase.xodd, 1, testCase.rf, struct([]) );
            testCase.verifyEqual( z, res, 'AbsTol', testCase.tol );

            % Decreasing order
            z = chopblas_sorted_sum_vec( testCase.xodd, 0, testCase.rf, struct([]) );
            testCase.verifyEqual( z, res, 'AbsTol', testCase.tol );

            % Test with a half-precision rounding mode
            res = testCase.xoddsorted(1);
            for i=2:1:length(testCase.xoddsorted)
                res = double( half( res + testCase.xoddsorted(i) ) );
            end

            z = chopblas_sorted_sum_vec( testCase.xodd, 1, testCase.rf, testCase.hopts );
            testCase.verifyEqual( z, res );
        end

        % Test with a matrix with an odd number of elements in each row
        function odd_length_mat(testCase)
            % Test with global rounding options
            % In double precision, this is the same value as the normal sum
            testCase.rf( [], testCase.dopts );

            res = sum( testCase.Xodd, 2 );

            % Increasing order
            z = chopblas_sorted_sum_mat( testCase.Xodd, 1, testCase.rf, struct([]) );
            testCase.verifyEqual( z, res, 'AbsTol', testCase.tol );

            % Decreasing order
            z = chopblas_sorted_sum_mat( testCase.Xodd, 0, testCase.rf, struct([]) );
            testCase.verifyEqual( z, res, 'AbsTol', testCase.tol );

            % Test with a half-precision rounding mode
            res = testCase.Xoddsorted(:,1);
            for i=2:1:size(testCase.Xoddsorted, 2)
                res = double( half( res + testCase.Xoddsorted(:,i) ) );
            end

            z = chopblas_sorted_sum_mat( testCase.Xodd, 1, testCase.rf, testCase.hopts );
            testCase.verifyEqual( z, res );
        end
    end
end
