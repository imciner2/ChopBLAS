classdef kernel_insertion_sum < matlab.unittest.TestCase
%KERNEL_INSERTION_SUM Test for the insertion summation compute kernel

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
        % Test with a vector with an even number of elements
        function even_length_vec(testCase)
            % Test with global rounding options
            % In double precision, this is the same value as the normal sum
            testCase.rf( [], testCase.dopts );

            res = sum( testCase.xeven );
            z   = chopblas_insertion_sum_vec( testCase.xeven, 1, testCase.rf, struct([]) );
            testCase.verifyEqual( z, res );

            % Test with round to nearest (tie to up)
            z = chopblas_insertion_sum_vec( testCase.xeven, 1, @(x, s) round(x), struct([]) );
            testCase.verifyEqual( z, 4.0 );
        end

        % Test with a matrix with an even number of elements
        function even_length_mat(testCase)
            % Test with global rounding options
            % In double precision, this is the same value as the normal sum
            testCase.rf( [], testCase.dopts );

            res = sum( testCase.Xeven, 2 );
            z   = chopblas_insertion_sum_mat( testCase.Xeven, 1, testCase.rf, struct([]) );
            testCase.verifyEqual( z, res );

            % Test with round to nearest (tie to up)
            z = chopblas_insertion_sum_mat( testCase.Xeven, 1, @(x, s) round(x), struct([]) );
            testCase.verifyEqual( z, 4.0*ones( size( testCase.Xeven, 1 ), 1 ) );
        end

        % Test with a vector with an odd number of elements
        function odd_length_vec(testCase)
            % Test with global rounding options
            % In double precision, this is the same value as the normal sum
            testCase.rf( [], testCase.dopts );

            res = sum( testCase.xodd );
            z   = chopblas_insertion_sum_vec( testCase.xodd, 1, testCase.rf, struct([]) );
            testCase.verifyEqual( z, res );

            % Test with round to nearest (tie to up)
            z = chopblas_insertion_sum_vec( testCase.xodd, 1, @(x, s) round(x), struct([]) );
            testCase.verifyEqual( z, 3.0 );
        end

        % Test with a matrix with an odd number of elements
        function odd_length_mat(testCase)
            % Test with global rounding options
            % In double precision, this is the same value as the normal sum
            testCase.rf( [], testCase.dopts );

            res = sum( testCase.Xodd, 2 );
            z   = chopblas_insertion_sum_mat( testCase.Xodd, 1, testCase.rf, struct([]) );
            testCase.verifyEqual( z, res );

            % Test with round to nearest (tie to up)
            z = chopblas_insertion_sum_mat( testCase.Xodd, 1, @(x, s) round(x), struct([]) );
            testCase.verifyEqual( z, 3.0*ones( size( testCase.Xodd, 1 ), 1 ) );
        end
    end
end
