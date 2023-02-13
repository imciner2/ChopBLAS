classdef nrm2 < matlab.unittest.TestCase
%NRM2 Unit tests for the xNRM2 class of functions

    properties
        sopts
        hopts
        dopts

        rf

        xint
        xdec

        xodd
        xoddunsorted

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

            testCase.xint = [1; 2; 3; 4; 5; -3; -2];
            testCase.xdec = [1.0005; 2.34; 3.24; 4];

            testCase.xodd = [1.0005; 2.34; 3.24; 4; 5.25678];

            testCase.xoddunsorted = testCase.xodd( randperm( length(testCase.xodd) ) );

            % Set a tolerance for all the tests
            testCase.tol = 1e-7;
        end
    end

    methods(Test)
        % No options specified, use global default of the double mode
        function chop_no_opts(testCase)
            testCase.rf( [], testCase.dopts );;

            % Test with global rounding options
            z = chnrm2( testCase.xint, ...
                        'Rounding', testCase.rf );
            testCase.verifyEqual( z, sqrt( testCase.xint'*testCase.xint ) );
        end

        % Use custom rounding function
        function chop_round_func(testCase)
            % Test with default of chop
            chop( [], testCase.dopts );

            z = chnrm2( testCase.xint );
            testCase.verifyEqual( z, sqrt( testCase.xint'*testCase.xint ) );

            % Test with trivial rounding function
            z = chnrm2( testCase.xint, ...
                        'Rounding', @(x, y) zeros(length(x), 1) );
            testCase.verifyEqual( z, 0 );
        end

        % Change the addition algorithm
        function chop_add_algorithm(testCase)
            % Test with default of chop
            chop( [], testCase.dopts );

            z = chnrm2( testCase.xodd, ...
                        'Rounding', testCase.rf, ...
                        'Accumulator', @chaccum_recursive );
            testCase.verifyEqual( z, sqrt( testCase.xodd'*testCase.xodd ), 'AbsTol', testCase.tol );

            % Test with recursive algorithm
            x = half( testCase.xodd.*testCase.xodd );
            res = x(1);
            for i=2:1:length(x)
                res = double( half( res + x(i) ) );
            end
            res = double( half( sqrt( res ) ) );

            z = chnrm2( testCase.xodd, testCase.hopts, ...
                        'Rounding', testCase.rf, ...
                        'Accumulator', @chaccum_recursive );
            testCase.verifyEqual( z, res );

            % Test with pairwise algorithm
            x = half( testCase.xodd.*testCase.xodd );
            i1  = double( half( x(1) + x(3) ) );
            i2  = double( half( x(2) + x(4) ) );
            i3  = double( half( i1 + i2 ) );
            res = double( half( i3 + x(5) ) );
            res = double( half( sqrt( res ) ) );

            z = chnrm2( testCase.xodd, testCase.hopts, ...
                        'Rounding', testCase.rf, ...
                        'Accumulator', @chaccum_pairwise );
            testCase.verifyEqual( z, res );

            % Test with recursive algorithm
            x = half( testCase.xodd.*testCase.xodd );
            res = x(1);
            for i=2:1:length(x)
                res = double( half( res + x(i) ) );
            end
            res = double( half( sqrt( res ) ) );

            z = chnrm2( testCase.xoddunsorted, testCase.hopts, ...
                        'Rounding', testCase.rf, ...
                        'Accumulator', @(x, func, opts) chaccum_sorted(x, func, opts, 'ascend' ) );
            testCase.verifyEqual( z, res );

            % Test with recursive algorithm
            x = half( testCase.xodd.*testCase.xodd );
            res = x(end);
            for i=length(x)-1:-1:1
                res = double( half( res + x(i) ) );
            end
            res = double( half( sqrt( res ) ) );

            z = chnrm2( testCase.xoddunsorted, testCase.hopts, ...
                        'Rounding', testCase.rf, ...
                        'Accumulator', @(x, func, opts) chaccum_sorted(x, func, opts, 'descend' ) );
            testCase.verifyEqual( z, res );
        end

        % Only one mode specified, it is used in all operations
        function chop_one_mode(testCase)
            z = chnrm2( testCase.xodd, testCase.hopts, ...
                        'Rounding', testCase.rf );

            c = 0;
            for i=1:length(testCase.xodd)
                i1 = double( half( testCase.xodd(i) * testCase.xodd(i) ) );
                c  = double( half( c + i1 ) );
            end
            c = double( half( sqrt( c ) ) );

            testCase.verifyEqual( z, c );
        end

        % Two modes specified, each operation has a different mode
        function chop_three_modes(testCase)
            c = 0;
            for i=1:length(testCase.xodd)
                c = half( c + double( single( testCase.xodd(i) * testCase.xodd(i) ) ) );
            end
            c = sqrt( double( c ) );

            zall = chnrm2( testCase.xodd, testCase.sopts, testCase.hopts, testCase.dopts, ...
                           'Rounding', testCase.rf );
            testCase.verifyEqual( zall, c );

            c = 0;
            for i=1:length(testCase.xodd)
                i1 = double( half( testCase.xodd(i) * testCase.xodd(i) ) );
                c  = double( half( c + i1 ) );
            end
            c = sqrt( c );

            z = chnrm2( testCase.xodd, testCase.hopts, testCase.hopts, testCase.dopts, ...
                        'Rounding', testCase.rf );
            testCase.verifyEqual( z, c );

            % The two test cases above should not produce the same result due to the half-precision rounding
            testCase.verifyNotEqual( z, zall );
        end
    end
end
