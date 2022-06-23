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
            %% Test the full function with full precision
            chop( [], testCase.dopts );

            %% Test the different triangular options
            % Upper-triangular, normal, non-unit triangular matrix
            z = chtrmv( 'u', 'n', 'n', testCase.Aint, testCase.xint );
            testCase.verifyEqual( z, triu( testCase.Aint )*testCase.xint, 'AbsTol', testCase.tol );

            z = chtrmv( 'U', 'n', 'n', testCase.Aint, testCase.xint );
            testCase.verifyEqual( z, triu( testCase.Aint )*testCase.xint, 'AbsTol', testCase.tol );

            % Lower-triangular, normal, non-unit triangular matrix
            z = chtrmv( 'l', 'n', 'n', testCase.Aint, testCase.xint );
            testCase.verifyEqual( z, tril( testCase.Aint )*testCase.xint, 'AbsTol', testCase.tol );

            z = chtrmv( 'L', 'n', 'n', testCase.Aint, testCase.xint );
            testCase.verifyEqual( z, tril( testCase.Aint )*testCase.xint, 'AbsTol', testCase.tol );

            %% Test the unit-traingular options
            Au = testCase.Aint;
            Au(eye(length(Au)) == 1) = 1;

            % Upper-triangular, normal, unit triangular matrix
            z = chtrmv( 'u', 'n', 'u', testCase.Aint, testCase.xint );
            testCase.verifyEqual( z, triu( Au )*testCase.xint, 'AbsTol', testCase.tol );

            z = chtrmv( 'u', 'n', 'U', testCase.Aint, testCase.xint );
            testCase.verifyEqual( z, triu( Au )*testCase.xint, 'AbsTol', testCase.tol );

            % Lower-triangular, normal, unit triangular matrix
            z = chtrmv( 'l', 'n', 'u', testCase.Aint, testCase.xint );
            testCase.verifyEqual( z, tril( Au )*testCase.xint, 'AbsTol', testCase.tol );

            z = chtrmv( 'l', 'n', 'U', testCase.Aint, testCase.xint );
            testCase.verifyEqual( z, tril( Au )*testCase.xint, 'AbsTol', testCase.tol );
        end

        %% Test the rounding in the matrix-vector product with only one precision
        function normal_operation_single_precision(testCase)
            % Upper-triangular, not unit-triangular
            z = chtrmv( 'u', 'n', 'n', testCase.Adec, testCase.xdec, testCase.hopts );

            c = chop( diag( testCase.Adec, 0 ) .* testCase.xdec, testCase.hopts );
            for i=1:length(c)
                for j=i:1:size(testCase.Adec,2)
                    if i~=j
                        c(i) = chop( c(i) + chop( testCase.Adec(i,j) * testCase.xdec(j), testCase.hopts ), testCase.hopts );
                    end
                end
            end

            testCase.verifyEqual( z, c );

            % Upper-triangular, unit-triangular
            z = chtrmv( 'u', 'n', 'u', testCase.Adec, testCase.xdec, testCase.hopts );

            c = testCase.xdec;
            for i=1:length(c)
                for j=i:1:size(testCase.Adec,2)
                    if i~=j
                        c(i) = chop( c(i) + chop( testCase.Adec(i,j) * testCase.xdec(j), testCase.hopts ), testCase.hopts );
                    end
                end
            end

            testCase.verifyEqual( z, c );

            % Lower-triangular, not unit-triangular
            z = chtrmv( 'l', 'n', 'n', testCase.Adec, testCase.xdec, testCase.hopts );

            c = chop( diag( testCase.Adec, 0 ) .* testCase.xdec, testCase.hopts );
            for i=1:length(c)
                for j=1:1:i
                    if i~=j
                        c(i) = chop( c(i) + chop( testCase.Adec(i,j) * testCase.xdec(j), testCase.hopts ), testCase.hopts );
                    end
                end
            end

            testCase.verifyEqual( z, c );

            % Lower-triangular, unit-triangular
            z = chtrmv( 'u', 'n', 'u', testCase.Adec, testCase.xdec, testCase.hopts );

            c = testCase.xdec;
            for i=1:length(c)
                for j=1:1:i
                    if i~=j
                        c(i) = chop( c(i) + chop( testCase.Adec(i,j) * testCase.xdec(j), testCase.hopts ), testCase.hopts );
                    end
                end
            end
        end

        %% Test the rounding in the matrix-vector product with two precisions
        function normal_operation_two_precision(testCase)
            % Upper-triangular, not unit-triangular
            z = chtrmv( 'u', 'n', 'n', testCase.Adec, testCase.xdec, testCase.sopts, testCase.hopts );

            c = chop( diag( testCase.Adec, 0 ) .* testCase.xdec, testCase.sopts );
            for i=1:length(c)
                for j=i:1:size(testCase.Adec,2)
                    if i~=j
                        c(i) = chop( c(i) + chop( testCase.Adec(i,j) * testCase.xdec(j), testCase.sopts ), testCase.hopts );
                    end
                end
            end

            testCase.verifyEqual( z, c );

            % Upper-triangular, unit-triangular
            z = chtrmv( 'u', 'n', 'u', testCase.Adec, testCase.xdec, testCase.sopts, testCase.hopts );

            c = testCase.xdec;
            for i=1:length(c)
                for j=i:1:size(testCase.Adec,2)
                    if i~=j
                        c(i) = chop( c(i) + chop( testCase.Adec(i,j) * testCase.xdec(j), testCase.sopts ), testCase.hopts );
                    end
                end
            end

            testCase.verifyEqual( z, c );

            % Lower-triangular, not unit-triangular
            z = chtrmv( 'l', 'n', 'n', testCase.Adec, testCase.xdec, testCase.sopts, testCase.hopts );

            c = chop( diag( testCase.Adec, 0 ) .* testCase.xdec, testCase.sopts );
            for i=1:length(c)
                for j=1:1:i
                    if i~=j
                        c(i) = chop( c(i) + chop( testCase.Adec(i,j) * testCase.xdec(j), testCase.sopts ), testCase.hopts );
                    end
                end
            end

            testCase.verifyEqual( z, c );

            % Lower-triangular, unit-triangular
            z = chtrmv( 'u', 'n', 'u', testCase.Adec, testCase.xdec, testCase.sopts, testCase.hopts );

            c = testCase.xdec;
            for i=1:length(c)
                for j=1:1:i
                    if i~=j
                        c(i) = chop( c(i) + chop( testCase.Adec(i,j) * testCase.xdec(j), testCase.sopts ), testCase.hopts );
                    end
                end
            end
        end

        function transposed_operation(testCase)
            alpha = 3;
            beta = 4;

            chop( [], testCase.dopts );

            % Upper-triangular, transposed, non-unit triangular matrix
            z = chtrmv( 'u', 't', 'n', testCase.Aint, testCase.xint );
            testCase.verifyEqual( z, triu( testCase.Aint )'*testCase.xint, 'AbsTol', testCase.tol );

            z = chtrmv( 'u', 'T', 'n', testCase.Aint, testCase.xint );
            testCase.verifyEqual( z, triu( testCase.Aint )'*testCase.xint, 'AbsTol', testCase.tol );

            % Lower-triangular, transposed, non-unit triangular matrix
            z = chtrmv( 'l', 't', 'n', testCase.Aint, testCase.xint );
            testCase.verifyEqual( z, tril( testCase.Aint )'*testCase.xint, 'AbsTol', testCase.tol );

            z = chtrmv( 'l', 'T', 'n', testCase.Aint, testCase.xint );
            testCase.verifyEqual( z, tril( testCase.Aint )'*testCase.xint, 'AbsTol', testCase.tol );

            %% Test the unit-traingular options
            Au = testCase.Aint;
            Au(eye(length(Au)) == 1) = 1;

            % Upper-triangular, transposed, unit triangular matrix
            z = chtrmv( 'u', 't', 'u', testCase.Aint, testCase.xint );
            testCase.verifyEqual( z, triu( Au )'*testCase.xint, 'AbsTol', testCase.tol );

            z = chtrmv( 'u', 't', 'U', testCase.Aint, testCase.xint );
            testCase.verifyEqual( z, triu( Au )'*testCase.xint, 'AbsTol', testCase.tol );

            % Lower-triangular, transposed, unit triangular matrix
            z = chtrmv( 'l', 't', 'u', testCase.Aint, testCase.xint );
            testCase.verifyEqual( z, tril( Au )'*testCase.xint, 'AbsTol', testCase.tol );

            z = chtrmv( 'l', 't', 'U', testCase.Aint, testCase.xint );
            testCase.verifyEqual( z, tril( Au )'*testCase.xint, 'AbsTol', testCase.tol );
        end

        %% Test the rounding in the matrix-vector product with only one precision
        function transposed_operation_single_precision(testCase)
            % Upper-triangular, not unit-triangular
            z = chtrmv( 'u', 't', 'n', testCase.Adec, testCase.xdec, testCase.hopts );

            c = chop( diag( testCase.Adec, 0 ) .* testCase.xdec, testCase.hopts );
            for i=1:length(c)
                for j=1:1:i
                    if i~=j
                        c(i) = chop( c(i) + chop( testCase.Adec(j,i) * testCase.xdec(j), testCase.hopts ), testCase.hopts );
                    end
                end
            end

            testCase.verifyEqual( z, c );

            % Upper-triangular, unit-triangular
            z = chtrmv( 'u', 't', 'u', testCase.Adec, testCase.xdec, testCase.hopts );

            c = testCase.xdec;
            for i=1:length(c)
                for j=1:1:i
                    if i~=j
                        c(i) = chop( c(i) + chop( testCase.Adec(j,i) * testCase.xdec(j), testCase.hopts ), testCase.hopts );
                    end
                end
            end

            testCase.verifyEqual( z, c );

            % Lower-triangular, not unit-triangular
            z = chtrmv( 'l', 't', 'n', testCase.Adec, testCase.xdec, testCase.hopts );

            c = chop( diag( testCase.Adec, 0 ) .* testCase.xdec, testCase.hopts );
            for i=1:length(c)
                for j=i:1:size(testCase.Adec,2)
                    if i~=j
                        c(i) = chop( c(i) + chop( testCase.Adec(j,i) * testCase.xdec(j), testCase.hopts ), testCase.hopts );
                    end
                end
            end

            testCase.verifyEqual( z, c );

            % Lower-triangular, unit-triangular
            z = chtrmv( 'u', 't', 'u', testCase.Adec, testCase.xdec, testCase.hopts );

            c = testCase.xdec;
            for i=1:length(c)
                for j=i:1:size(testCase.Adec,2)
                    if i~=j
                        c(i) = chop( c(i) + chop( testCase.Adec(j,i) * testCase.xdec(j), testCase.hopts ), testCase.hopts );
                    end
                end
            end
        end

        %% Test the rounding in the matrix-vector product with two precisions
        function transposed_operation_two_precision(testCase)
            % Upper-triangular, not unit-triangular
            z = chtrmv( 'u', 't', 'n', testCase.Adec, testCase.xdec, testCase.sopts, testCase.hopts );

            c = chop( diag( testCase.Adec, 0 ) .* testCase.xdec, testCase.sopts );
            for i=1:length(c)
                for j=1:1:i
                    if i~=j
                        c(i) = chop( c(i) + chop( testCase.Adec(j,i) * testCase.xdec(j), testCase.sopts ), testCase.hopts );
                    end
                end
            end

            testCase.verifyEqual( z, c );

            % Upper-triangular, unit-triangular
            z = chtrmv( 'u', 't', 'u', testCase.Adec, testCase.xdec, testCase.sopts, testCase.hopts );

            c = testCase.xdec;
            for i=1:length(c)
                for j=1:1:i
                    if i~=j
                        c(i) = chop( c(i) + chop( testCase.Adec(j,i) * testCase.xdec(j), testCase.sopts ), testCase.hopts );
                    end
                end
            end

            testCase.verifyEqual( z, c );

            % Lower-triangular, not unit-triangular
            z = chtrmv( 'l', 't', 'n', testCase.Adec, testCase.xdec, testCase.sopts, testCase.hopts );

            c = chop( diag( testCase.Adec, 0 ) .* testCase.xdec, testCase.sopts );
            for i=1:length(c)
                for j=i:1:size(testCase.Adec,2)
                    if i~=j
                        c(i) = chop( c(i) + chop( testCase.Adec(j,i) * testCase.xdec(j), testCase.sopts ), testCase.hopts );
                    end
                end
            end

            testCase.verifyEqual( z, c );

            % Lower-triangular, unit-triangular
            z = chtrmv( 'u', 't', 'u', testCase.Adec, testCase.xdec, testCase.sopts, testCase.hopts );

            c = testCase.xdec;
            for i=1:length(c)
                for j=i:1:size(testCase.Adec,2)
                    if i~=j
                        c(i) = chop( c(i) + chop( testCase.Adec(j,i) * testCase.xdec(j), testCase.sopts ), testCase.hopts );
                    end
                end
            end
        end
    end
end
