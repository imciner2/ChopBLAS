classdef trmv < matlab.unittest.TestCase
%TRMV Unit tests for the xTRMV class of functions

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
        function normal_operation(testCase)
            %% Test the full function with full precision
            testCase.rf( [], testCase.dopts );

            %% Test the different triangular options
            % Upper-triangular, normal, non-unit triangular matrix
            z = chtrmv( testCase.Aint, testCase.xint, ...
                        'Rounding', testCase.rf );
            testCase.verifyEqual( z, triu( testCase.Aint )*testCase.xint, 'AbsTol', testCase.tol );

            z = chtrmv( testCase.Aint, testCase.xint, ...
                        'LowerTriangular', false, ...
                        'Rounding', testCase.rf );
            testCase.verifyEqual( z, triu( testCase.Aint )*testCase.xint, 'AbsTol', testCase.tol );

            % Lower-triangular, normal, non-unit triangular matrix
            z = chtrmv( testCase.Aint, testCase.xint, ...
                        'LowerTriangular', true, ...
                        'Rounding', testCase.rf );
            testCase.verifyEqual( z, tril( testCase.Aint )*testCase.xint, 'AbsTol', testCase.tol );

            %% Test the unit-traingular options
            Au = testCase.Aint;
            Au(eye(length(Au)) == 1) = 1;

            % Upper-triangular, normal, unit triangular matrix
            z = chtrmv( testCase.Aint, testCase.xint, ...
                        'UnitTriangular', true, ...
                        'Rounding', testCase.rf );
            testCase.verifyEqual( z, triu( Au )*testCase.xint, 'AbsTol', testCase.tol );

            % Lower-triangular, normal, unit triangular matrix
            z = chtrmv( testCase.Aint, testCase.xint, ...
                        'LowerTriangular', true, ...
                        'UnitTriangular', true, ...
                        'Rounding', testCase.rf );
            testCase.verifyEqual( z, tril( Au )*testCase.xint, 'AbsTol', testCase.tol );
        end

        function normal_operation_rounding_function(testCase)
            rf = @(x, y) zeros(length(x), 1);
            z0 = zeros(length(testCase.xint), 1);
            z1 = z0;
            ze = z0;

            z1(1) = testCase.xint(1);
            ze(end) = testCase.xint(end);

            % Upper-triangular, normal, non-unit triangular matrix
            z = chtrmv( testCase.Aint, testCase.xint, 'Rounding', rf );
            testCase.verifyEqual( z, z0, 'AbsTol', testCase.tol );

            z = chtrmv( testCase.Aint, testCase.xint, ...
                        'LowerTriangular', false, ...
                        'Rounding', rf );
            testCase.verifyEqual( z, z0, 'AbsTol', testCase.tol );

            % Lower-triangular, normal, non-unit triangular matrix
            z = chtrmv( testCase.Aint, testCase.xint, ...
                        'LowerTriangular', true, ...
                        'Rounding', rf );
            testCase.verifyEqual( z, z0, 'AbsTol', testCase.tol );

            %% Test the unit-triangular options
            % Upper-triangular, normal, unit triangular matrix, don't round last xout element
            z = chtrmv( testCase.Aint, testCase.xint, ...
                        'UnitTriangular', true, ...
                        'Rounding', rf );
            testCase.verifyEqual( z, ze, 'AbsTol', testCase.tol );

            % Upper-triangular, normal, unit triangular matrix, round last xout element
            z = chtrmv( testCase.Aint, testCase.xint, ...
                        'UnitTriangular', true, ...
                        'Rounding', rf, ...
                        'RoundAll', true );
            testCase.verifyEqual( z, z0, 'AbsTol', testCase.tol );

            % Lower-triangular, normal, unit triangular matrix, don't round first xout element
            z = chtrmv( testCase.Aint, testCase.xint, ...
                        'LowerTriangular', true, ...
                        'UnitTriangular', true, ...
                        'Rounding', rf );
            testCase.verifyEqual( z, z1, 'AbsTol', testCase.tol );

            % Lower-triangular, normal, unit triangular matrix, round first xout element
            z = chtrmv( testCase.Aint, testCase.xint, ...
                        'LowerTriangular', true, ...
                        'UnitTriangular', true, ...
                        'Rounding', rf, ...
                        'RoundAll', true );
            testCase.verifyEqual( z, z0, 'AbsTol', testCase.tol );
        end

        %% Test the rounding in the matrix-vector product with only one precision
        function normal_operation_single_precision(testCase)
            % Upper-triangular, not unit-triangular
            z = chtrmv( testCase.Adec, testCase.xdec, testCase.hopts, ...
                        'Rounding', testCase.rf );

            c = testCase.rf( diag( testCase.Adec, 0 ) .* testCase.xdec, testCase.hopts );
            for i=1:length(c)
                for j=i:1:size(testCase.Adec,2)
                    if i~=j
                        c(i) = testCase.rf( c(i) + testCase.rf( testCase.Adec(i,j) * testCase.xdec(j), testCase.hopts ), testCase.hopts );
                    end
                end
            end

            testCase.verifyEqual( z, c );

            % Upper-triangular, unit-triangular
            z = chtrmv( testCase.Adec, testCase.xdec, testCase.hopts, ...
                        'UnitTriangular', true, ...
                        'Rounding', testCase.rf );

            c = testCase.xdec;
            for i=1:length(c)
                for j=i:1:size(testCase.Adec,2)
                    if i~=j
                        c(i) = testCase.rf( c(i) + testCase.rf( testCase.Adec(i,j) * testCase.xdec(j), testCase.hopts ), testCase.hopts );
                    end
                end
            end

            testCase.verifyEqual( z, c );

            % Lower-triangular, not unit-triangular
            z = chtrmv( testCase.Adec, testCase.xdec, testCase.hopts, ...
                        'LowerTriangular', true, ...
                        'Rounding', testCase.rf );

            c = testCase.rf( diag( testCase.Adec, 0 ) .* testCase.xdec, testCase.hopts );
            for i=1:length(c)
                for j=1:1:i
                    if i~=j
                        c(i) = testCase.rf( c(i) + testCase.rf( testCase.Adec(i,j) * testCase.xdec(j), testCase.hopts ), testCase.hopts );
                    end
                end
            end

            testCase.verifyEqual( z, c );

            % Lower-triangular, unit-triangular
            z = chtrmv( testCase.Adec, testCase.xdec, testCase.hopts, ...
                        'LowerTriangular', true, ...
                        'UnitTriangular', true, ...
                        'Rounding', testCase.rf );

            c = testCase.xdec;
            for i=1:length(c)
                for j=1:1:i
                    if i~=j
                        c(i) = testCase.rf( c(i) + testCase.rf( testCase.Adec(i,j) * testCase.xdec(j), testCase.hopts ), testCase.hopts );
                    end
                end
            end
        end

        %% Test the rounding in the matrix-vector product with two precisions
        function normal_operation_two_precision(testCase)
            % Upper-triangular, not unit-triangular
            z = chtrmv( testCase.Adec, testCase.xdec, testCase.sopts, testCase.hopts, ...
                        'Rounding', testCase.rf );

            c = testCase.rf( diag( testCase.Adec, 0 ) .* testCase.xdec, testCase.sopts );
            for i=1:length(c)
                for j=i:1:size(testCase.Adec,2)
                    if i~=j
                        c(i) = testCase.rf( c(i) + testCase.rf( testCase.Adec(i,j) * testCase.xdec(j), testCase.sopts ), testCase.hopts );
                    end
                end
            end

            testCase.verifyEqual( z, c );

            % Upper-triangular, unit-triangular
            z = chtrmv( testCase.Adec, testCase.xdec, testCase.sopts, testCase.hopts, ...
                        'UnitTriangular', true, ...
                        'Rounding', testCase.rf );

            c = testCase.xdec;
            for i=1:length(c)
                for j=i:1:size(testCase.Adec,2)
                    if i~=j
                        c(i) = testCase.rf( c(i) + testCase.rf( testCase.Adec(i,j) * testCase.xdec(j), testCase.sopts ), testCase.hopts );
                    end
                end
            end

            testCase.verifyEqual( z, c );

            % Lower-triangular, not unit-triangular
            z = chtrmv( testCase.Adec, testCase.xdec, testCase.sopts, testCase.hopts, ...
                        'LowerTriangular', true, ...
                        'Rounding', testCase.rf );

            c = testCase.rf( diag( testCase.Adec, 0 ) .* testCase.xdec, testCase.sopts );
            for i=1:length(c)
                for j=1:1:i
                    if i~=j
                        c(i) = testCase.rf( c(i) + testCase.rf( testCase.Adec(i,j) * testCase.xdec(j), testCase.sopts ), testCase.hopts );
                    end
                end
            end

            testCase.verifyEqual( z, c );

            % Lower-triangular, unit-triangular
            z = chtrmv( testCase.Adec, testCase.xdec, testCase.sopts, testCase.hopts, ...
                        'LowerTriangular', true, ...
                        'UnitTriangular', true, ...
                        'Rounding', testCase.rf );

            c = testCase.xdec;
            for i=1:length(c)
                for j=1:1:i
                    if i~=j
                        c(i) = testCase.rf( c(i) + testCase.rf( testCase.Adec(i,j) * testCase.xdec(j), testCase.sopts ), testCase.hopts );
                    end
                end
            end
        end

        function transposed_operation(testCase)
            alpha = 3;
            beta = 4;

            testCase.rf( [], testCase.dopts );

            % Upper-triangular, transposed, non-unit triangular matrix
            z = chtrmv( testCase.Aint, testCase.xint, ...
                        'Transpose', true, ...
                        'Rounding', testCase.rf );
            testCase.verifyEqual( z, triu( testCase.Aint )'*testCase.xint, 'AbsTol', testCase.tol );

            % Lower-triangular, transposed, non-unit triangular matrix
            z = chtrmv( testCase.Aint, testCase.xint, ...
                        'Transpose', true, ...
                        'LowerTriangular', true, ...
                        'Rounding', testCase.rf );
            testCase.verifyEqual( z, tril( testCase.Aint )'*testCase.xint, 'AbsTol', testCase.tol );

            %% Test the unit-triangular options
            Au = testCase.Aint;
            Au(eye(length(Au)) == 1) = 1;

            % Upper-triangular, transposed, unit triangular matrix
            z = chtrmv( testCase.Aint, testCase.xint, ...
                        'Transpose', true, ...
                        'UnitTriangular', true, ...
                        'Rounding', testCase.rf );
            testCase.verifyEqual( z, triu( Au )'*testCase.xint, 'AbsTol', testCase.tol );

            % Lower-triangular, transposed, unit triangular matrix
            z = chtrmv( testCase.Aint, testCase.xint, ...
                        'LowerTriangular', true, ...
                        'Transpose', true, ...
                        'UnitTriangular', true, ...
                        'Rounding', testCase.rf );
            testCase.verifyEqual( z, tril( Au )'*testCase.xint, 'AbsTol', testCase.tol );
        end

        function transposed_operation_rounding_function(testCase)
            rf = @(x, y) zeros(length(x), 1);
            z0 = zeros(length(testCase.xint), 1);
            z1 = z0;
            ze = z0;

            z1(1) = testCase.xint(1);
            ze(end) = testCase.xint(end);

            % Upper-triangular, transposed, non-unit triangular matrix
            z = chtrmv( testCase.Aint, testCase.xint, ...
                        'Transpose', true, ...
                        'Rounding', rf );
            testCase.verifyEqual( z, z0, 'AbsTol', testCase.tol );

            % Lower-triangular, transposed, non-unit triangular matrix
            z = chtrmv( testCase.Aint, testCase.xint, ...
                        'Transpose', true, ...
                        'LowerTriangular', true, ...
                        'Rounding', rf );
            testCase.verifyEqual( z, z0, 'AbsTol', testCase.tol );

            %% Test the unit-triangular options
            % Upper-triangular, normal, unit triangular matrix, don't round first xout element
            z = chtrmv( testCase.Aint, testCase.xint, ...
                        'Transpose', true, ...
                        'UnitTriangular', true, ...
                        'Rounding', rf );
            testCase.verifyEqual( z, z1, 'AbsTol', testCase.tol );

            % Upper-triangular, normal, unit triangular matrix, round first xout element
            z = chtrmv( testCase.Aint, testCase.xint, ...
                        'Transpose', true, ...
                        'UnitTriangular', true, ...
                        'Rounding', rf, ...
                        'RoundAll', true );
            testCase.verifyEqual( z, z0, 'AbsTol', testCase.tol );

            % Lower-triangular, normal, unit triangular matrix, don't round last xout element
            z = chtrmv( testCase.Aint, testCase.xint, ...
                        'Transpose', true, ...
                        'LowerTriangular', true, ...
                        'UnitTriangular', true, ...
                        'Rounding', rf );
            testCase.verifyEqual( z, ze, 'AbsTol', testCase.tol );

            % Lower-triangular, normal, unit triangular matrix, round last xout element
            z = chtrmv( testCase.Aint, testCase.xint, ...
                        'Transpose', true, ...
                        'LowerTriangular', true, ...
                        'UnitTriangular', true, ...
                        'Rounding', rf, ...
                        'RoundAll', true );
            testCase.verifyEqual( z, z0, 'AbsTol', testCase.tol );
        end

        %% Test the rounding in the matrix-vector product with only one precision
        function transposed_operation_single_precision(testCase)
            % Upper-triangular, not unit-triangular
            z = chtrmv( testCase.Adec, testCase.xdec, testCase.hopts, ...
                        'Transpose', true, ...
                        'Rounding', testCase.rf );

            c = testCase.rf( diag( testCase.Adec, 0 ) .* testCase.xdec, testCase.hopts );
            for i=1:length(c)
                for j=1:1:i
                    if i~=j
                        c(i) = testCase.rf( c(i) + testCase.rf( testCase.Adec(j,i) * testCase.xdec(j), testCase.hopts ), testCase.hopts );
                    end
                end
            end

            testCase.verifyEqual( z, c );

            % Upper-triangular, unit-triangular
            z = chtrmv( testCase.Adec, testCase.xdec, testCase.hopts, ...
                        'Transpose', true, ...
                        'UnitTriangular', true, ...
                        'Rounding', testCase.rf );

            c = testCase.xdec;
            for i=1:length(c)
                for j=1:1:i
                    if i~=j
                        c(i) = testCase.rf( c(i) + testCase.rf( testCase.Adec(j,i) * testCase.xdec(j), testCase.hopts ), testCase.hopts );
                    end
                end
            end

            testCase.verifyEqual( z, c );

            % Lower-triangular, not unit-triangular
            z = chtrmv( testCase.Adec, testCase.xdec, testCase.hopts, ...
                        'Transpose', true, ...
                        'LowerTriangular', true, ...
                        'Rounding', testCase.rf );

            c = testCase.rf( diag( testCase.Adec, 0 ) .* testCase.xdec, testCase.hopts );
            for i=1:length(c)
                for j=i:1:size(testCase.Adec,2)
                    if i~=j
                        c(i) = testCase.rf( c(i) + testCase.rf( testCase.Adec(j,i) * testCase.xdec(j), testCase.hopts ), testCase.hopts );
                    end
                end
            end

            testCase.verifyEqual( z, c );

            % Lower-triangular, unit-triangular
            z = chtrmv( testCase.Adec, testCase.xdec, testCase.hopts, ...
                        'Transpose', true, ...
                        'LowerTriangular', true, ...
                        'UnitTriangular', true, ...
                        'Rounding', testCase.rf );

            c = testCase.xdec;
            for i=1:length(c)
                for j=i:1:size(testCase.Adec,2)
                    if i~=j
                        c(i) = testCase.rf( c(i) + testCase.rf( testCase.Adec(j,i) * testCase.xdec(j), testCase.hopts ), testCase.hopts );
                    end
                end
            end
        end

        %% Test the rounding in the matrix-vector product with two precisions
        function transposed_operation_two_precision(testCase)
            % Upper-triangular, not unit-triangular
            z = chtrmv( testCase.Adec, testCase.xdec, testCase.sopts, testCase.hopts, ...
                        'Transpose', true, ...
                        'Rounding', testCase.rf );

            c = testCase.rf( diag( testCase.Adec, 0 ) .* testCase.xdec, testCase.sopts );
            for i=1:length(c)
                for j=1:1:i
                    if i~=j
                        c(i) = testCase.rf( c(i) + testCase.rf( testCase.Adec(j,i) * testCase.xdec(j), testCase.sopts ), testCase.hopts );
                    end
                end
            end

            testCase.verifyEqual( z, c );

            % Upper-triangular, unit-triangular
            z = chtrmv( testCase.Adec, testCase.xdec, testCase.sopts, testCase.hopts, ...
                        'Transpose', true, ...
                        'UnitTriangular', true, ...
                        'Rounding', testCase.rf );

            c = testCase.xdec;
            for i=1:length(c)
                for j=1:1:i
                    if i~=j
                        c(i) = testCase.rf( c(i) + testCase.rf( testCase.Adec(j,i) * testCase.xdec(j), testCase.sopts ), testCase.hopts );
                    end
                end
            end

            testCase.verifyEqual( z, c );

            % Lower-triangular, not unit-triangular
            z = chtrmv( testCase.Adec, testCase.xdec, testCase.sopts, testCase.hopts, ...
                        'Transpose', true, ...
                        'LowerTriangular', true, ...
                        'Rounding', testCase.rf );

            c = testCase.rf( diag( testCase.Adec, 0 ) .* testCase.xdec, testCase.sopts );
            for i=1:length(c)
                for j=i:1:size(testCase.Adec,2)
                    if i~=j
                        c(i) = testCase.rf( c(i) + testCase.rf( testCase.Adec(j,i) * testCase.xdec(j), testCase.sopts ), testCase.hopts );
                    end
                end
            end

            testCase.verifyEqual( z, c );

            % Lower-triangular, unit-triangular
            z = chtrmv( testCase.Adec, testCase.xdec, testCase.sopts, testCase.hopts, ...
                        'Transpose', true, ...
                        'LowerTriangular', true, ...
                        'UnitTriangular', true, ...
                        'Rounding', testCase.rf );

            c = testCase.xdec;
            for i=1:length(c)
                for j=i:1:size(testCase.Adec,2)
                    if i~=j
                        c(i) = testCase.rf( c(i) + testCase.rf( testCase.Adec(j,i) * testCase.xdec(j), testCase.sopts ), testCase.hopts );
                    end
                end
            end
        end
    end
end
