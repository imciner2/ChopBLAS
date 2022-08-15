classdef trsv < matlab.unittest.TestCase
%TRSV Unit tests for the xTRSV class of functions

    properties
        sopts
        hopts
        dopts
        tol

        Adec
        Aeye

        xdec
        bdec

        rf
    end

    methods(TestMethodSetup)
        function setup_test(testCase)
            testCase.sopts.format = 's';
            testCase.hopts.format = 'h';
            testCase.dopts.format = 'd';

            testCase.bdec = [2.00125; 3.00225; 4.00014; 5.0025];

            bdec = testCase.bdec;

            testCase.Aeye = eye(4);
            testCase.Adec = [bdec, bdec, bdec, bdec];

            % Ensure there is a solution to the system Adec*x=bdec by making it positive definite
            testCase.Adec = testCase.Adec'*testCase.Adec;

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
            A = [2, 3, 4;
                 0, 6, 7];
            b = [2; 3; 4];

            % Wrong vector orientation
            testCase.verifyError( @() chtrsv( A, b' ), "chtrsv:bMustBeColumnVector" );

            % Wrong sizes for A and x
            testCase.verifyError( @() chtrsv( A(:,1:2), b ), "chtrsv:AbMustBeCompatibleSizes" );
            testCase.verifyError( @() chtrsv( A, b, 'Transpose', true ), "chtrsv:AbMustBeCompatibleSizes" );
        end

        function normal_operation(testCase)
            %% Test the full function with full precision
            testCase.rf( [], testCase.dopts );

            %% Test the different triangular options
            % Upper-triangular, normal, non-unit triangular matrix
            z = chtrsv( testCase.Adec, testCase.bdec, ...
                        'Rounding', testCase.rf );
            testCase.verifyEqual( z, triu( testCase.Adec )\testCase.bdec, 'AbsTol', testCase.tol );

            z = chtrsv( testCase.Adec, testCase.bdec, ...
                        'LowerTriangular', false, ...
                        'Rounding', testCase.rf );
            testCase.verifyEqual( z, triu( testCase.Adec )\testCase.bdec, 'AbsTol', testCase.tol );

            % Lower-triangular, normal, non-unit triangular matrix
            z = chtrsv( testCase.Adec, testCase.bdec, ...
                        'LowerTriangular', true, ...
                        'Rounding', testCase.rf );
            testCase.verifyEqual( z, tril( testCase.Adec )\testCase.bdec, 'AbsTol', testCase.tol );

            %% Test the unit-traingular options
            Au = testCase.Adec;
            Au(eye(length(Au)) == 1) = 1;

            % Upper-triangular, normal, unit triangular matrix
            z = chtrsv( testCase.Adec, testCase.bdec, ...
                        'UnitTriangular', true, ...
                        'Rounding', testCase.rf );
            testCase.verifyEqual( z, triu( Au )\testCase.bdec, 'AbsTol', testCase.tol );

            % Lower-triangular, normal, unit triangular matrix
            z = chtrsv( testCase.Adec, testCase.bdec, ...
                        'LowerTriangular', true, ...
                        'UnitTriangular', true, ...
                        'Rounding', testCase.rf );
            testCase.verifyEqual( z, tril( Au )\testCase.bdec, 'AbsTol', testCase.tol );
        end

        function normal_operation_rounding_function(testCase)
            rf = @(x, y) zeros(length(x), 1);
            z0 = zeros(length(testCase.bdec), 1);
            z1 = z0;
            ze = z0;

            z1(1) = testCase.bdec(1);
            ze(end) = testCase.bdec(end);

            % Upper-triangular, normal, non-unit triangular matrix
            z = chtrsv( testCase.Adec, testCase.bdec, 'Rounding', rf );
            testCase.verifyEqual( z, z0, 'AbsTol', testCase.tol );

            z = chtrsv( testCase.Adec, testCase.bdec, ...
                        'LowerTriangular', false, ...
                        'Rounding', rf );
            testCase.verifyEqual( z, z0, 'AbsTol', testCase.tol );

            % Lower-triangular, normal, non-unit triangular matrix
            z = chtrsv( testCase.Adec, testCase.bdec, ...
                        'LowerTriangular', true, ...
                        'Rounding', rf );
            testCase.verifyEqual( z, z0, 'AbsTol', testCase.tol );

            %% Test the unit-triangular options
            % Upper-triangular, normal, unit triangular matrix, round last x element
            z = chtrsv( testCase.Adec, testCase.bdec, ...
                        'UnitTriangular', true, ...
                        'Rounding', rf );
            testCase.verifyEqual( z, z0, 'AbsTol', testCase.tol );

            % Upper-triangular, normal, unit triangular matrix, Don't round last x element
            z = chtrsv( testCase.Adec, testCase.bdec, ...
                        'UnitTriangular', true, ...
                        'Rounding', rf, ...
                        'RoundAll', false );
            testCase.verifyEqual( z, ze, 'AbsTol', testCase.tol );

            % Lower-triangular, normal, unit triangular matrix, round first xout element
            z = chtrsv( testCase.Adec, testCase.bdec, ...
                        'LowerTriangular', true, ...
                        'UnitTriangular', true, ...
                        'Rounding', rf );
            testCase.verifyEqual( z, z0, 'AbsTol', testCase.tol );

            % Lower-triangular, normal, unit triangular matrix, don't round first xout element
            z = chtrsv( testCase.Adec, testCase.bdec, ...
                        'LowerTriangular', true, ...
                        'UnitTriangular', true, ...
                        'Rounding', rf, ...
                        'RoundAll', false );
            testCase.verifyEqual( z, z1, 'AbsTol', testCase.tol );
        end

        %% Test the rounding in the matrix-vector product with only one precision
        function normal_operation_single_precision(testCase)
            A = [2, 3, 0, 0;
                 5, 2, 3, 0;
                 0, 5, 2, 3;
                 0, 0, 5, 2];
            b = [4; 5; 6; 7];

            rf = @(x, s) floor(x);

            %% Upper-triangular, not unit-triangular
            z = chtrsv( A, b, testCase.hopts, ...
                        'Rounding', rf );

            % Verify only certain entries
            testCase.verifyEqual( z(end), 3.0 );        % 7/2 rounded down to 3
            testCase.verifyEqual( z(end-1), -2.0 );     % ( 6 - 3*3 ) / 2 rounded down to -2

            %% Upper-triangular, unit-triangular
            z = chtrsv( A, b, testCase.hopts, ...
                        'UnitTriangular', true, ...
                        'Rounding', rf );

            % Verify only certain entries
            testCase.verifyEqual( z(end), 7.0 );        % 7/1 rounded down to 7
            testCase.verifyEqual( z(end-1), -15.0 );    % ( 6 - 7*3 ) / 1 rounded down to -15

            % Lower-triangular, not unit-triangular
            z = chtrsv( A, b, testCase.hopts, ...
                        'LowerTriangular', true, ...
                        'Rounding', rf );

            % Verify only certain entries
            testCase.verifyEqual( z(1), 2.0 );      % 4/2 rounded down to 2
            testCase.verifyEqual( z(2), -3.0 );     % ( 5 - 5*2 ) / 2 rounded down to -3

            % Lower-triangular, unit-triangular
            z = chtrsv( A, b, testCase.hopts, ...
                        'LowerTriangular', true, ...
                        'UnitTriangular', true, ...
                        'Rounding', rf );

            testCase.verifyEqual( z(1), 4.0 );      % 4/1 rounded down to 4
            testCase.verifyEqual( z(2), -15.0 );    % ( 5 - 5*4 ) / 1 rounded down to -15
        end

        %% Test the rounding in the matrix-vector product with three precisions
        function normal_operation_three_precision(testCase)
            A = [2, 3, 0, 0;
                 5, 2, 3, 0;
                 0, 5, 2, 3;
                 0, 0, 5, 2];
            b = [4; 5; 6; 7];

            rf = @(x, s) round(x, s.num);
            addopts.num = 0;
            mulopts.num = 0;
            divopts.num = 1;

            %% Upper-triangular, not unit-triangular
            z = chtrsv( A, b, mulopts, addopts, divopts, ...
                        'Rounding', rf );

            % Verify only certain entries
            testCase.verifyEqual( z(end), 3.5 );        % 7/2 rounded to 3.5
            testCase.verifyEqual( z(end-1), -2.5 );     % ( 6 - 3.5*3 ) / 2 rounded to -2.5 (3.5*3 is rounded to 11)

            %% Upper-triangular, unit-triangular
            z = chtrsv( A, b, mulopts, addopts, divopts, ...
                        'UnitTriangular', true, ...
                        'Rounding', rf );

            % Verify only certain entries
            testCase.verifyEqual( z(end), 7.0 );        % 7/1 rounded to 7
            testCase.verifyEqual( z(end-1), -15.0 );    % ( 6 - 7*3 ) / 1 rounded to -15

            % Lower-triangular, not unit-triangular
            z = chtrsv( A, b, mulopts, addopts, divopts, ...
                        'LowerTriangular', true, ...
                        'Rounding', rf );

            % Verify only certain entries
            testCase.verifyEqual( z(1), 2.0 );      % 4/2 rounded to 2
            testCase.verifyEqual( z(2), -2.5 );     % ( 5 - 5*2 ) / 2 rounded to -2.5

            % Lower-triangular, unit-triangular
            z = chtrsv( A, b, mulopts, addopts, divopts, ...
                        'LowerTriangular', true, ...
                        'UnitTriangular', true, ...
                        'Rounding', rf );

            testCase.verifyEqual( z(1), 4.0 );      % 4/1 rounded to 4
            testCase.verifyEqual( z(2), -15.0 );    % ( 5 - 5*4 ) / 1 rounded to -15
        end

        function transposed_operation(testCase)
            alpha = 3;
            beta = 4;

            testCase.rf( [], testCase.dopts );

            % Upper-triangular, transposed, non-unit triangular matrix
            z = chtrsv( testCase.Adec, testCase.bdec, ...
                        'Transpose', true, ...
                        'Rounding', testCase.rf );
            testCase.verifyEqual( z, triu( testCase.Adec )'\testCase.bdec, 'AbsTol', testCase.tol );

            % Lower-triangular, transposed, non-unit triangular matrix
            z = chtrsv( testCase.Adec, testCase.bdec, ...
                        'Transpose', true, ...
                        'LowerTriangular', true, ...
                        'Rounding', testCase.rf );
            testCase.verifyEqual( z, tril( testCase.Adec )'\testCase.bdec, 'AbsTol', testCase.tol );

            %% Test the unit-triangular options
            Au = testCase.Adec;
            Au(eye(length(Au)) == 1) = 1;

            % Upper-triangular, transposed, unit triangular matrix
            z = chtrsv( testCase.Adec, testCase.bdec, ...
                        'Transpose', true, ...
                        'UnitTriangular', true, ...
                        'Rounding', testCase.rf );
            testCase.verifyEqual( z, triu( Au )'\testCase.bdec, 'AbsTol', testCase.tol );

            % Lower-triangular, transposed, unit triangular matrix
            z = chtrsv( testCase.Adec, testCase.bdec, ...
                        'LowerTriangular', true, ...
                        'Transpose', true, ...
                        'UnitTriangular', true, ...
                        'Rounding', testCase.rf );
            testCase.verifyEqual( z, tril( Au )'\testCase.bdec, 'AbsTol', testCase.tol );
        end

        function transposed_operation_rounding_function(testCase)
            rf = @(x, y) zeros(length(x), 1);
            z0 = zeros(length(testCase.bdec), 1);
            z1 = z0;
            ze = z0;

            z1(1) = testCase.bdec(1);
            ze(end) = testCase.bdec(end);

            % Upper-triangular, transposed, non-unit triangular matrix
            z = chtrsv( testCase.Adec, testCase.bdec, ...
                        'Transpose', true, ...
                        'Rounding', rf );
            testCase.verifyEqual( z, z0, 'AbsTol', testCase.tol );

            % Lower-triangular, transposed, non-unit triangular matrix
            z = chtrsv( testCase.Adec, testCase.bdec, ...
                        'Transpose', true, ...
                        'LowerTriangular', true, ...
                        'Rounding', rf );
            testCase.verifyEqual( z, z0, 'AbsTol', testCase.tol );

            %% Test the unit-triangular options
            % Upper-triangular, normal, unit triangular matrix, round first xout element
            z = chtrsv( testCase.Adec, testCase.bdec, ...
                        'Transpose', true, ...
                        'UnitTriangular', true, ...
                        'Rounding', rf );
            testCase.verifyEqual( z, z0, 'AbsTol', testCase.tol );

            % Upper-triangular, normal, unit triangular matrix, don't round first xout element
            z = chtrsv( testCase.Adec, testCase.bdec, ...
                        'Transpose', true, ...
                        'UnitTriangular', true, ...
                        'Rounding', rf, ...
                        'RoundAll', false );
            testCase.verifyEqual( z, z1, 'AbsTol', testCase.tol );

            % Lower-triangular, normal, unit triangular matrix, round last xout element
            z = chtrsv( testCase.Adec, testCase.bdec, ...
                        'Transpose', true, ...
                        'LowerTriangular', true, ...
                        'UnitTriangular', true, ...
                        'Rounding', rf );
            testCase.verifyEqual( z, z0, 'AbsTol', testCase.tol );

            % Lower-triangular, normal, unit triangular matrix, don't round last xout element
            z = chtrsv( testCase.Adec, testCase.bdec, ...
                        'Transpose', true, ...
                        'LowerTriangular', true, ...
                        'UnitTriangular', true, ...
                        'Rounding', rf, ...
                        'RoundAll', false );
            testCase.verifyEqual( z, ze, 'AbsTol', testCase.tol );
        end

        %% Test the rounding in the matrix-vector product with only one precision
        function transposed_operation_single_precision(testCase)
            A = [2, 3, 0, 0;
                 5, 2, 3, 0;
                 0, 5, 2, 3;
                 0, 0, 5, 2];
            b = [4; 5; 6; 7];

            rf = @(x, s) floor(x);

            %% Upper-triangular, not unit-triangular
            z = chtrsv( A, b, testCase.hopts, ...
                        'Transpose', true, ...
                        'Rounding', rf );

            % Verify only certain entries
            testCase.verifyEqual( z(1), 2.0 );        % 4/2 rounded down to 2
            testCase.verifyEqual( z(2), -1.0 );       % ( 5 - 2*3 ) / 2 rounded down to -1

            %% Upper-triangular, unit-triangular
            z = chtrsv( A, b, testCase.hopts, ...
                        'Transpose', true, ...
                        'UnitTriangular', true, ...
                        'Rounding', rf );

            % Verify only certain entries
            testCase.verifyEqual( z(1), 4.0 );        % 4/1 rounded down to 4
            testCase.verifyEqual( z(2), -7.0 );       % ( 5 - 4*3 ) / 1 rounded down to -7

            % Lower-triangular, not unit-triangular
            z = chtrsv( A, b, testCase.hopts, ...
                        'Transpose', true, ...
                        'LowerTriangular', true, ...
                        'Rounding', rf );

            % Verify only certain entries
            testCase.verifyEqual( z(end), 3.0 );      % 7/2 rounded down to 3
            testCase.verifyEqual( z(end-1), -5.0 );   % ( 6 - 5*3 ) / 2 rounded down to -5

            % Lower-triangular, unit-triangular
            z = chtrsv( A, b, testCase.hopts, ...
                        'Transpose', true, ...
                        'LowerTriangular', true, ...
                        'UnitTriangular', true, ...
                        'Rounding', rf );

            testCase.verifyEqual( z(end), 7.0 );      % 7/1 rounded down to 7
            testCase.verifyEqual( z(end-1), -29.0 );  % ( 6 - 7*5 ) / 1 rounded down to -29
        end

        %% Test the rounding in the matrix-vector product with three precisions
        function transposed_operation_three_precision(testCase)
            A = [2, 3, 0, 0;
                 5, 2, 3, 0;
                 0, 5, 2, 3;
                 0, 0, 5, 2];
            b = [4; 5; 6; 7];

            rf = @(x, s) round(x, s.num);
            addopts.num = 0;
            mulopts.num = 0;
            divopts.num = 1;

            %% Upper-triangular, not unit-triangular
            z = chtrsv( A, b, mulopts, addopts, divopts, ...
                        'Transpose', true, ...
                        'Rounding', rf );

            % Verify only certain entries
            testCase.verifyEqual( z(1), 2.0 );        % 4/2 rounded to 2
            testCase.verifyEqual( z(2), -0.5 );       % ( 5 - 2*3 ) / 2 rounded to -0.5

            %% Upper-triangular, unit-triangular
            z = chtrsv( A, b, mulopts, addopts, divopts, ...
                        'Transpose', true, ...
                        'UnitTriangular', true, ...
                        'Rounding', rf );

            % Verify only certain entries
            testCase.verifyEqual( z(1), 4.0 );        % 4/1 rounded to 4
            testCase.verifyEqual( z(2), -7.0 );       % ( 5 - 4*3 ) / 1 rounded to -7

            % Lower-triangular, not unit-triangular
            z = chtrsv( A, b, mulopts, addopts, divopts, ...
                        'Transpose', true, ...
                        'LowerTriangular', true, ...
                        'Rounding', rf );

            % Verify only certain entries
            testCase.verifyEqual( z(end), 3.5 );      % 7/2 rounded  to 3.5
            testCase.verifyEqual( z(end-1), -6.0 );   % ( 6 - 5*3.5 ) / 2 rounded to -6 (5*3.5 is rounded to 18)

            % Lower-triangular, unit-triangular
            z = chtrsv( A, b, mulopts, addopts, divopts, ...
                        'Transpose', true, ...
                        'LowerTriangular', true, ...
                        'UnitTriangular', true, ...
                        'Rounding', rf );

            testCase.verifyEqual( z(end), 7.0 );      % 7/1 rounded to 7
            testCase.verifyEqual( z(end-1), -29.0 );  % ( 6 - 7*5 ) / 1 rounded to -29
        end
    end
end
