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

            % Set the default format for chop to double precision (mimic no rounding)
            chop( [], testCase.dopts );

            % Set a tolerance for all the tests
            testCase.tol = 1e-4;
        end
    end

    methods(Test)
        function normal_operation(testCase)
            alpha = 3;
            beta = 4;

            %% Test the full function with full precision
            chop( [], testCase.dopts );

            z = chgemv( alpha, testCase.Aint, testCase.xint, beta, testCase.yint );
            testCase.verifyEqual( z, (alpha.*(testCase.Aint*testCase.xint) ) + beta.*testCase.yint, 'AbsTol', testCase.tol );

            %% Test the rounding in the matrix-vector product with only one precision
            z = chgemv( 1, testCase.Adec, testCase.xdec, 0, testCase.ydec, testCase.hopts );

            c = zeros(length(testCase.xdec), 1);
            for i=1:length(c)
                for j=1:1:size(testCase.Adec,2)
                    c(i) = chop( c(i) + chop( testCase.Adec(i,j) * testCase.xdec(j), testCase.hopts ), testCase.hopts );
                end
            end

            testCase.verifyEqual( z, c );

            %% Test the rounding in the matrix-vector product with two different precisions
            z = chgemv( 1, testCase.Adec, testCase.xdec, 0, testCase.ydec, testCase.sopts, testCase.hopts );

            c = zeros(length(testCase.xdec), 1);
            for i=1:length(c)
                for j=1:1:size(testCase.Adec,2)
                    c(i) = chop( c(i) + chop( testCase.Adec(i,j) * testCase.xdec(j), testCase.sopts ), testCase.hopts );
                end
            end
            c = chop( c, testCase.sopts );

            testCase.verifyEqual( z, c );
        end

        function transposed_operation(testCase)
            alpha = 3;
            beta = 4;

            %% Test the transposed function with full precision
            chop( [], testCase.dopts );

            z = chgemv( alpha, testCase.Aint, testCase.xint, beta, testCase.yint, 'Transpose', true );
            testCase.verifyEqual( z, (alpha.*(testCase.Aint'*testCase.xint) ) + beta.*testCase.yint, 'AbsTol', testCase.tol );

            %% Test the rounding in the matrix-vector product with only one precision
            z = chgemv( 1, testCase.Adec, testCase.xdec, 0, testCase.ydec, testCase.hopts, 'Transpose', true );

            c = zeros(length(testCase.xdec), 1);
            for i=1:length(c)
                for j=1:1:size(testCase.Adec,2)
                    c(i) = chop( c(i) + chop( testCase.Adec(j,i) * testCase.xdec(j), testCase.hopts ), testCase.hopts );
                end
            end

            testCase.verifyEqual( z, c );

            %% Test the rounding in the matrix-vector product with two different precisions
            z = chgemv( 1, testCase.Adec, testCase.xdec, 0, testCase.ydec, testCase.sopts, testCase.hopts, 'Transpose', true );

            c = zeros(length(testCase.xdec), 1);
            for i=1:length(c)
                for j=1:1:size(testCase.Adec,2)
                    c(i) = chop( c(i) + chop( testCase.Adec(j,i) * testCase.xdec(j), testCase.sopts ), testCase.hopts );
                end
            end

            testCase.verifyEqual( z, c );
        end

        function accumulate_rounding(testCase)
            alpha = 3;
            beta = 4;

            %% Only one precision used
            % No rounding for beta=1
            z = chgemv( 1, testCase.Aeye, testCase.xdec, 1, testCase.ydec, testCase.hopts );

            % The specification is that y is preloaded into the output and then accumulated onto
            c = testCase.ydec;
            for i=1:length(c)
                for j=1:1:size(testCase.Aeye, 2)
                    c(i) = chop( c(i) + chop( testCase.Aeye(i,j) * testCase.xdec(j), testCase.hopts ), testCase.hopts );
                end
            end

            testCase.verifyEqual( z, c );

            % Rounding for non-unit beta
            z = chgemv( 1, testCase.Aeye, testCase.xdec, 3, testCase.ydec, testCase.hopts );

            % The specification is that y is preloaded into the output and then accumulated onto
            c = chop( 3.*testCase.ydec, testCase.hopts );
            for i=1:length(c)
                for j=1:1:size(testCase.Aeye, 2)
                    c(i) = chop( c(i) + chop( testCase.Aeye(i,j) * testCase.xdec(j), testCase.hopts ), testCase.hopts );
                end
            end

            testCase.verifyEqual( z, c );

            %% Two precisions used
            % No rounding for beta=1
            z = chgemv( 1, testCase.Aeye, testCase.xdec, 1, testCase.ydec, testCase.sopts, testCase.hopts );

            % The specification is that y is preloaded into the output and then accumulated onto
            c = testCase.ydec;
            for i=1:length(c)
                for j=1:1:size(testCase.Aeye, 2)
                    c(i) = chop( c(i) + chop( testCase.Aeye(i,j) * testCase.xdec(j), testCase.sopts ), testCase.hopts );
                end
            end

            testCase.verifyEqual( z, c );

            % Rounding for non-unit beta
            z = chgemv( 1, testCase.Aeye, testCase.xdec, 3, testCase.ydec, testCase.sopts, testCase.hopts );

            % The specification is that y is preloaded into the output with rounding done using mulopts, and then accumulated onto
            c = chop( 3.*testCase.ydec, testCase.sopts );
            for i=1:length(c)
                for j=1:1:size(testCase.Aeye, 2)
                    c(i) = chop( c(i) + chop( testCase.Aeye(i,j) * testCase.xdec(j), testCase.sopts ), testCase.hopts );
                end
            end

            testCase.verifyEqual( z, c );
        end

        % No options specified, use global default of the double mode
        function mvm_rounding(testCase)
            alpha = 3;
            beta = 4;

            %% Test with no accumulate vector with full precision
            chop( [], testCase.dopts );

            z = chgemv( alpha, testCase.Adec, testCase.xdec, 0, testCase.yint );
            testCase.verifyEqual( z, (alpha.*(testCase.Adec*testCase.xdec) ), 'AbsTol', testCase.tol );

            z = chgemv( alpha, testCase.Adec, testCase.xdec, 2, [] );
            testCase.verifyEqual( z, (alpha.*(testCase.Adec*testCase.xdec) ), 'AbsTol', testCase.tol );

            z = chgemv( alpha, testCase.Adec, testCase.xdec, 0, testCase.yint, 'Transpose', true );
            testCase.verifyEqual( z, (alpha.*(testCase.Adec'*testCase.xdec) ), 'AbsTol', testCase.tol );

            z = chgemv( alpha, testCase.Adec, testCase.xdec, 2, [], 'Transpose', true );
            testCase.verifyEqual( z, (alpha.*(testCase.Adec'*testCase.xdec) ), 'AbsTol', testCase.tol );

            %% Test the rounding in the scaling of the matrix-vector result with one precision
            % Normal matrix
            z = chgemv( alpha, testCase.Adec, testCase.xdec, 0, testCase.ydec, testCase.hopts );

            c  = zeros(length(testCase.xdec), 1);
            tx = chop( alpha.*testCase.xdec, testCase.hopts );
            for i=1:length(c)
                for j=1:1:size(testCase.Adec, 2)
                    c(i) = chop( c(i) + chop( testCase.Adec(i,j) * tx(j), testCase.hopts ), testCase.hopts );
                end
            end

            testCase.verifyEqual( z, c );

            % Transposed matrix
            z = chgemv( alpha, testCase.Adec, testCase.xdec, 0, testCase.ydec, testCase.hopts, 'Transpose', true );

            c  = zeros(length(testCase.xdec), 1);
            tx = chop( alpha.*testCase.xdec, testCase.hopts );
            for i=1:length(c)
                for j=1:1:size(testCase.Adec,2)
                    c(i) = chop( c(i) + chop( testCase.Adec(j,i) * tx(j), testCase.hopts ), testCase.hopts );
                end
            end

            testCase.verifyEqual( z, c );

            %% Test the rounding in the scaling of the matrix-vector result with two precisions
            % Normal matrix
            z = chgemv( alpha, testCase.Adec, testCase.xdec, 0, testCase.ydec, testCase.sopts, testCase.hopts );

            c = zeros(length(testCase.xdec), 1);
            tx = chop( alpha.*testCase.xdec, testCase.sopts );
            for i=1:length(c)
                for j=1:1:size(testCase.Adec, 2)
                    c(i) = chop( c(i) + chop( testCase.Adec(i,j) * tx(j), testCase.sopts ), testCase.hopts );
                end
            end

            testCase.verifyEqual( z, c );

            % Transposed matrix
            z = chgemv( alpha, testCase.Adec, testCase.xdec, 0, testCase.ydec, testCase.sopts, testCase.hopts, 'Transpose', true );

            c  = zeros(length(testCase.xdec), 1);
            tx = chop( alpha.*testCase.xdec, testCase.hopts );
            for i=1:length(c)
                for j=1:1:size(testCase.Adec,2)
                    c(i) = chop( c(i) + chop( testCase.Adec(j,i) * tx(j), testCase.sopts ), testCase.hopts );
                end
            end

            testCase.verifyEqual( z, c );
        end
    end
end