classdef gemv < matlab.unittest.TestCase
%GEMV Unit tests for the xGEMV class of functions

    properties
        sopts
        hopts
        dopts
        tol

        Aint
        Adec
        Aeye

        xint
        yint

        xdec
        ydec
        alpha
        beta

        rf
    end

    methods(TestMethodSetup)
        function setup_test(testCase)
            testCase.sopts.format = 's';
            testCase.hopts.format = 'h';
            testCase.dopts.format = 'd';

            testCase.xint = [1; 2; 3; 4; 5; -3; -2];
            testCase.yint = [2; 3; 4; 5; 6;  7;  8];

            testCase.xdec = [1.0005; 2.34; 3.24; 4];
            testCase.ydec = [2.00125; 3.00225; 4.00014; 5.0025];

            yint = testCase.yint;
            ydec = testCase.ydec;

            testCase.Aeye = eye(4);
            testCase.Adec = [ydec, ydec, ydec, ydec];
            testCase.Aint = [yint, yint, yint, yint, yint, yint, yint];

            testCase.alpha = 3;
            testCase.beta  = 4;

            % Choose the rounding function to use
            if strcmpi( getenv('CHOPBLAS_ROUND_FUNC'), 'cpfloat' )
                testCase.rf = @cpfloat;
            else
                % Default to chop
                testCase.rf = @chop;
            end

            % Set a tolerance for all the tests
            testCase.tol = 1e-4;
        end
    end

    methods(Test)
        function argument_checks(testCase)
            alpha = 1.0;
            beta = 2.0;
            A = [2, 3, 4;
                 5, 6, 7];
            x = [2; 3; 4];
            y = [4; 5; 6];

            testCase.verifyError( @() chgemv( [2; 2], A, x, beta, y ), "chgemv:AlphaMustBeScalar" );
            testCase.verifyError( @() chgemv( alpha, A, x, [2; 2], y ), "chgemv:BetaMustBeScalar" );

            % Wrong vector lengths
            testCase.verifyError( @() chgemv( alpha, A, x, beta, y(1:2) ), "chgemv:xyMustBeSameSize" );
            testCase.verifyError( @() chgemv( alpha, A, x(1:2), beta, y ), "chgemv:xyMustBeSameSize" );

            % Wrong vector orientations
            testCase.verifyError( @() chgemv( alpha, A, x, beta, y' ), "chgemv:xyMustBeColumnVectors" );
            testCase.verifyError( @() chgemv( alpha, A, x', beta, y ), "chgemv:xyMustBeColumnVectors" );
            testCase.verifyError( @() chgemv( alpha, A, x', beta, y' ), "chgemv:xyMustBeColumnVectors" );

            % Wrong sizes for A and x
            testCase.verifyError( @() chgemv( alpha, A(:,1:2), x, beta, y ), "chgemv:AxMustBeCompatibleSizes" );
            testCase.verifyError( @() chgemv( alpha, A, x, beta, y, 'Transpose', true ), "chgemv:AxMustBeCompatibleSizes" );
        end

        function normal_operation(testCase)
            %% Test the full function with full precision
            testCase.rf( [], testCase.dopts );

            z = chgemv( testCase.alpha, testCase.Aint, testCase.xint, testCase.beta, testCase.yint, 'Rounding', testCase.rf );
            testCase.verifyEqual( z, (testCase.alpha.*(testCase.Aint*testCase.xint) ) + testCase.beta.*testCase.yint, 'AbsTol', testCase.tol );

            %% Test the rounding in the matrix-vector product with only one precision
            z = chgemv( 1, testCase.Adec, testCase.xdec, 0, testCase.ydec, testCase.hopts );

            c = zeros(length(testCase.xdec), 1);
            for i=1:length(c)
                for j=1:1:size(testCase.Adec,2)
                    c(i) = testCase.rf( c(i) + testCase.rf( testCase.Adec(i,j) * testCase.xdec(j), testCase.hopts ), testCase.hopts );
                end
            end

            testCase.verifyEqual( z, c );

            %% Test the rounding in the matrix-vector product with two different precisions
            z = chgemv( 1, testCase.Adec, testCase.xdec, 0, testCase.ydec, testCase.sopts, testCase.hopts, 'Rounding', testCase.rf );

            c = zeros(length(testCase.xdec), 1);
            for i=1:length(c)
                for j=1:1:size(testCase.Adec,2)
                    c(i) = testCase.rf( c(i) + testCase.rf( testCase.Adec(i,j) * testCase.xdec(j), testCase.sopts ), testCase.hopts );
                end
            end
            c = testCase.rf( c, testCase.sopts );

            testCase.verifyEqual( z, c );
        end

        % No options specified, use global default of the double mode
        function rounding_function(testCase)
            %% Test with default rounding function
            chop( [], testCase.dopts );

            z = chgemv( testCase.alpha, testCase.Adec, testCase.xdec, 0, testCase.ydec );
            testCase.verifyEqual( z, (testCase.alpha.*(testCase.Adec*testCase.xdec) ), 'AbsTol', testCase.tol );

            % Test with trivial rounding function
            z = chgemv( testCase.alpha, testCase.Adec, testCase.xdec, 2, [], 'Rounding', @(x, y) zeros(length(x), 1) );
            testCase.verifyEqual( z, zeros(length(testCase.xdec), 1) );
        end

        function transposed_operation(testCase)
            %% Test the transposed function with full precision
            testCase.rf( [], testCase.dopts );

            z = chgemv( testCase.alpha, testCase.Aint, testCase.xint, testCase.beta, testCase.yint, 'Transpose', true, 'Rounding', testCase.rf );
            testCase.verifyEqual( z, (testCase.alpha.*(testCase.Aint'*testCase.xint) ) + testCase.beta.*testCase.yint, 'AbsTol', testCase.tol );

            %% Test the rounding in the matrix-vector product with only one precision
            z = chgemv( 1, testCase.Adec, testCase.xdec, 0, testCase.ydec, testCase.hopts, 'Transpose', true, 'Rounding', testCase.rf );

            c = zeros(length(testCase.xdec), 1);
            for i=1:length(c)
                for j=1:1:size(testCase.Adec,2)
                    c(i) = testCase.rf( c(i) + testCase.rf( testCase.Adec(j,i) * testCase.xdec(j), testCase.hopts ), testCase.hopts );
                end
            end

            testCase.verifyEqual( z, c );

            %% Test the rounding in the matrix-vector product with two different precisions
            z = chgemv( 1, testCase.Adec, testCase.xdec, 0, testCase.ydec, testCase.sopts, testCase.hopts, 'Transpose', true, 'Rounding', testCase.rf );

            c = zeros(length(testCase.xdec), 1);
            for i=1:length(c)
                for j=1:1:size(testCase.Adec,2)
                    c(i) = testCase.rf( c(i) + testCase.rf( testCase.Adec(j,i) * testCase.xdec(j), testCase.sopts ), testCase.hopts );
                end
            end

            testCase.verifyEqual( z, c );
        end

        function accumulate_rounding(testCase)
            %% Only one precision used
            % No rounding for beta=1
            z = chgemv( 1, testCase.Aeye, testCase.xdec, 1, testCase.ydec, testCase.hopts, 'Rounding', testCase.rf );

            % The specification is that y is preloaded into the output and then accumulated onto
            c = testCase.ydec;
            for i=1:length(c)
                for j=1:1:size(testCase.Aeye, 2)
                    c(i) = testCase.rf( c(i) + testCase.rf( testCase.Aeye(i,j) * testCase.xdec(j), testCase.hopts ), testCase.hopts );
                end
            end

            testCase.verifyEqual( z, c );

            % Rounding for non-unit beta
            z = chgemv( 1, testCase.Aeye, testCase.xdec, 3, testCase.ydec, testCase.hopts, 'Rounding', testCase.rf );

            % The specification is that y is preloaded into the output and then accumulated onto
            c = testCase.rf( 3.*testCase.ydec, testCase.hopts );
            for i=1:length(c)
                for j=1:1:size(testCase.Aeye, 2)
                    c(i) = testCase.rf( c(i) + testCase.rf( testCase.Aeye(i,j) * testCase.xdec(j), testCase.hopts ), testCase.hopts );
                end
            end

            testCase.verifyEqual( z, c );

            %% Two precisions used
            % No rounding for beta=1
            z = chgemv( 1, testCase.Aeye, testCase.xdec, 1, testCase.ydec, testCase.sopts, testCase.hopts, 'Rounding', testCase.rf );

            % The specification is that y is preloaded into the output and then accumulated onto
            c = testCase.ydec;
            for i=1:length(c)
                for j=1:1:size(testCase.Aeye, 2)
                    c(i) = testCase.rf( c(i) + testCase.rf( testCase.Aeye(i,j) * testCase.xdec(j), testCase.sopts ), testCase.hopts );
                end
            end

            testCase.verifyEqual( z, c );

            % Rounding for non-unit beta
            z = chgemv( 1, testCase.Aeye, testCase.xdec, 3, testCase.ydec, testCase.sopts, testCase.hopts, 'Rounding', testCase.rf );

            % The specification is that y is preloaded into the output with rounding done using mulopts, and then accumulated onto
            c = testCase.rf( 3.*testCase.ydec, testCase.sopts );
            for i=1:length(c)
                for j=1:1:size(testCase.Aeye, 2)
                    c(i) = testCase.rf( c(i) + testCase.rf( testCase.Aeye(i,j) * testCase.xdec(j), testCase.sopts ), testCase.hopts );
                end
            end

            testCase.verifyEqual( z, c );
        end

        % No options specified, use global default of the double mode
        function mvm_rounding(testCase)
            %% Test with no accumulate vector with full precision
            testCase.rf( [], testCase.dopts );

            z = chgemv( testCase.alpha, testCase.Adec, testCase.xdec, 0, testCase.ydec, 'Rounding', testCase.rf );
            testCase.verifyEqual( z, (testCase.alpha.*(testCase.Adec*testCase.xdec) ), 'AbsTol', testCase.tol );

            z = chgemv( testCase.alpha, testCase.Adec, testCase.xdec, 2, [], 'Rounding', testCase.rf );
            testCase.verifyEqual( z, (testCase.alpha.*(testCase.Adec*testCase.xdec) ), 'AbsTol', testCase.tol );

            z = chgemv( testCase.alpha, testCase.Adec, testCase.xdec, 0, testCase.ydec, 'Transpose', true, 'Rounding', testCase.rf );
            testCase.verifyEqual( z, (testCase.alpha.*(testCase.Adec'*testCase.xdec) ), 'AbsTol', testCase.tol );

            z = chgemv( testCase.alpha, testCase.Adec, testCase.xdec, 2, [], 'Transpose', true, 'Rounding', testCase.rf );
            testCase.verifyEqual( z, (testCase.alpha.*(testCase.Adec'*testCase.xdec) ), 'AbsTol', testCase.tol );

            %% Test the rounding in the scaling of the matrix-vector result with one precision
            % Normal matrix
            z = chgemv( testCase.alpha, testCase.Adec, testCase.xdec, 0, testCase.ydec, testCase.hopts, 'Rounding', testCase.rf );

            c  = zeros(length(testCase.xdec), 1);
            tx = testCase.rf( testCase.alpha.*testCase.xdec, testCase.hopts );
            for i=1:length(c)
                for j=1:1:size(testCase.Adec, 2)
                    c(i) = testCase.rf( c(i) + testCase.rf( testCase.Adec(i,j) * tx(j), testCase.hopts ), testCase.hopts );
                end
            end

            testCase.verifyEqual( z, c );

            % Transposed matrix
            z = chgemv( testCase.alpha, testCase.Adec, testCase.xdec, 0, testCase.ydec, testCase.hopts, 'Transpose', true, 'Rounding', testCase.rf );

            c  = zeros(length(testCase.xdec), 1);
            tx = testCase.rf( testCase.alpha.*testCase.xdec, testCase.hopts );
            for i=1:length(c)
                for j=1:1:size(testCase.Adec,2)
                    c(i) = testCase.rf( c(i) + testCase.rf( testCase.Adec(j,i) * tx(j), testCase.hopts ), testCase.hopts );
                end
            end

            testCase.verifyEqual( z, c );

            %% Test the rounding in the scaling of the matrix-vector result with two precisions
            % Normal matrix
            z = chgemv( testCase.alpha, testCase.Adec, testCase.xdec, 0, testCase.ydec, testCase.sopts, testCase.hopts, 'Rounding', testCase.rf );

            c = zeros(length(testCase.xdec), 1);
            tx = testCase.rf( testCase.alpha.*testCase.xdec, testCase.sopts );
            for i=1:length(c)
                for j=1:1:size(testCase.Adec, 2)
                    c(i) = testCase.rf( c(i) + testCase.rf( testCase.Adec(i,j) * tx(j), testCase.sopts ), testCase.hopts );
                end
            end

            testCase.verifyEqual( z, c );

            % Transposed matrix
            z = chgemv( testCase.alpha, testCase.Adec, testCase.xdec, 0, testCase.ydec, testCase.sopts, testCase.hopts, 'Transpose', true, 'Rounding', testCase.rf );

            c  = zeros(length(testCase.xdec), 1);
            tx = testCase.rf( testCase.alpha.*testCase.xdec, testCase.hopts );
            for i=1:length(c)
                for j=1:1:size(testCase.Adec,2)
                    c(i) = testCase.rf( c(i) + testCase.rf( testCase.Adec(j,i) * tx(j), testCase.sopts ), testCase.hopts );
                end
            end

            testCase.verifyEqual( z, c );
        end
    end
end
