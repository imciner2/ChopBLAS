classdef ger < matlab.unittest.TestCase
%GER Unit tests for the xGER class of functions

    properties
        sopts
        hopts
        dopts

        rf

        xsq
        ysq
        Asq

        xrect
        yrect
        Arect

        alpha

        tol
        dtol
    end

    methods(TestMethodSetup)
        function setup_test(testCase)
            testCase.sopts.format = 's';
            testCase.hopts.format = 'h';
            testCase.dopts.format = 'd';

            if strcmpi( getenv('CHOPBLAS_ROUND_FUNC'), 'cpfloat' )
                testCase.rf = @cpfloat;
            else
                % Default to chop
                testCase.rf = @chop;
            end

            testCase.xsq = [1; 2; 3; 4; 5; -3; -2];
            testCase.ysq = [2; 3; 4; 5; 6;  7;  8];
            testCase.Asq = repmat( testCase.xsq', length(testCase.xsq), 1 );

            testCase.xrect = [1.0005; 2.34; 3.24; 4];
            testCase.yrect = [2.00125; 3.00225];
            testCase.Arect = repmat( testCase.xrect', 2, 1 )';

            testCase.alpha = 3;

            % Set a tolerance for all the tests
            testCase.tol  = 1e-4;
            testCase.dtol = 1e-9;
        end
    end

    methods(Test)
        % Test argument validation
        function argument_checks(testCase)
            % Alpha not a scalar
            testCase.verifyError( @() chger( [2; 2], [1; 2], [3; 4], eye(2) ), "chger:AlphaMustBeScalar" );

            % Wrong vector lengths
            testCase.verifyError( @() chger( [2], [1; 2; 3], [3; 4], eye(2) ), "chger:xAMustHaveCompatibleSize" );
            testCase.verifyError( @() chger( [2], [1; 2], [3; 4; 5], eye(2) ), "chger:yAMustHaveCompatibleSize" );
        end

        % No options specified, use global default of the double mode
        function chop_no_opts(testCase)
            testCase.rf( [], testCase.dopts );

            % Square matrix
            res = (testCase.alpha.*testCase.xsq*testCase.ysq') + testCase.Asq;
            z   = chger( testCase.alpha, testCase.xsq, testCase.ysq, testCase.Asq, ...
                         'Rounding', testCase.rf );
            testCase.verifyEqual( z, res, 'AbsTol', testCase.dtol );

            % Square matrix, alpha of 0
            z   = chger( 0, testCase.xsq, testCase.ysq, testCase.Asq, ...
                         'Rounding', testCase.rf );
            testCase.verifyEqual( z, testCase.Asq, 'AbsTol', testCase.dtol );

            % Rectangular matrix
            res = (testCase.alpha.*testCase.xrect*testCase.yrect') + testCase.Arect;
            z   = chger( testCase.alpha, testCase.xrect, testCase.yrect, testCase.Arect, ...
                         'Rounding', testCase.rf );
            testCase.verifyEqual( z, res, 'AbsTol', testCase.dtol );

            % Rectangular matrix, alpha of 0
            z   = chger( 0, testCase.xrect, testCase.yrect, testCase.Arect, ...
                         'Rounding', testCase.rf );
            testCase.verifyEqual( z, testCase.Arect, 'AbsTol', testCase.dtol );
        end

        % Specify the rounding function
        function chop_round_func(testCase)
            % Use the default rounding function (chop)
            chop( [], testCase.dopts );

            res = (testCase.alpha.*testCase.xsq*testCase.ysq') + testCase.Asq;
            z   = chger( testCase.alpha, testCase.xsq, testCase.ysq, testCase.Asq );
            testCase.verifyEqual( z, res );

            % Test a trivial rounding function to ensure it overrides it
            z = chger(  testCase.alpha, testCase.xsq, testCase.ysq, testCase.Asq, ...
                        'Rounding', @(x, y) zeros(size(x)) );
            testCase.verifyEqual( z, zeros(size(testCase.Asq)) );
        end

        % Only one mode specified, it is used in both operations
        function chop_one_mode(testCase)
            r.num = 0;
            rf = @(x, s) round(x, s.num);

            % Rounding function for alpha computation
            z = chger( 0.5, [1; 2], [3; 4], zeros(2,2), r, 'Rounding', rf );
            testCase.verifyEqual( z(1,1), 2.0 );    % (1*3) * 0.5 = 1.5, rounded up to 2.0
            testCase.verifyEqual( z(2,2), 4.0 );    % (2*4) * 0.5 = 4 rounded down to 4.0

            % Rounding function for outer product computation
            z = chger( 1.0, [1.5; 2], [3; 4], zeros(2,2), r, 'Rounding', rf );
            testCase.verifyEqual( z(1,1), 5.0 );    % (1.5*3) = 4.5, rounded up to 2.0
            testCase.verifyEqual( z(1,2), 6.0 );    % (1.5*4) = 6 rounded down to 6.0

            % Rounding function for addition with A
            z = chger( 1.0, [1; 1], [1; 1], 0.25*ones(2,2), r, 'Rounding', rf );
            testCase.verifyEqual( z(1,1), 1.0 );    % 1.0 + 0.25 = 1.25, rounded down to 1.0
            testCase.verifyEqual( z(1,2), 1.0 );    % 1.0 + 0.25 = 1.25, rounded down to 1.0
        end

        % Two modes specified, each operation has a different mode
        function chop_two_modes(testCase)
            rm.num = 0;     % Rounding for mulopts
            ra.num = 1;     % Rounding for addopts
            rf = @(x, s) round(x, s.num);

            % Rounding function for alpha computation
            z = chger( 0.5, [1; 2], [3; 4], zeros(2,2), rm, ra, 'Rounding', rf );
            testCase.verifyEqual( z(1,1), 2.0 );    % (1*3) * 0.5 = 1.5, rounded up to 2.0
            testCase.verifyEqual( z(2,2), 4.0 );    % (2*4) * 0.5 = 4 rounded down to 4.0

            % Rounding function for outer product computation
            z = chger( 1.0, [1.5; 2], [3; 4], zeros(2,2), rm, ra, 'Rounding', rf );
            testCase.verifyEqual( z(1,1), 5.0 );    % 1.5*3 = 4.5, rounded up to 2.0
            testCase.verifyEqual( z(1,2), 6.0 );    % 1.5*4 = 6 rounded down to 6.0

            % Rounding function for addition with A
            z = chger( 1.0, [1; 1], [1; 1], 0.25*ones(2,2), rm, ra, 'Rounding', rf );
            testCase.verifyEqual( z(1,1), 1.3 );    % 1.0 + 0.25 = 1.3, rounded up to 1.3
            testCase.verifyEqual( z(1,2), 1.3 );    % 1.0 + 0.25 = 1.3, rounded up to 1.3
        end
    end
end
