%
% matlab.unittest automatic test code for bicas.proc.L1L2.demuxer.
%
% Very basic tests at this stage. Could be improved but unsure how much is
% meaningful.
%
%
% Author: Erik P G Johansson, IRF, Uppsala, Sweden
% First created 2021-09-08, using older test code.
%
classdef demuxer___UTEST < matlab.unittest.TestCase



    %##############
    %##############
    % TEST METHODS
    %##############
    %##############
    methods(Test)

        
        
        % Test two function in combination.
        %
        % IMPLEMENTATION NOTE: The design is for historical reasons before the
        % two functions were split up.
        function test_get_routings_calibrated_BLTSs_to_ASRs(testCase)
            
            function test(demuxMode, dlrUsing12, bltsSamplesAVolt, ExpSrcArray, ExpAsrSamplesAVolt)
                ActRoutingArray    = bicas.proc.L1L2.demuxer.get_routings(demuxMode, dlrUsing12);
                ActAsrSamplesAVolt = bicas.proc.L1L2.demuxer.calibrated_BLTSs_to_ASRs(ActRoutingArray, bltsSamplesAVolt);
                testCase.verifyEqual([ActRoutingArray.src], ExpSrcArray)
                testCase.verifyEqual(ActAsrSamplesAVolt,    ExpAsrSamplesAVolt)
            end

            % =========
            % Test data
            % =========
            V1   = 10;
            V2   = 11;
            V3   = 12;
            V12  = V1-V2;
            V13  = V1-V3;
            V23  = V2-V3;
            V12a = 45-56;
            V13a = 45-69;
            V23a = 56-69;

            % Create ASRs using constants. Arguments determine which ASRs should
            % be NaN instead of constants.
            function AsrSamplesVolt = ASR_samples(varargin)
                import bicas.proc.L1L2.demuxer___UTEST.as
                
                % varargin{i} == 0 or 1. Determines whether constant or NaN will
                % be used.
                assert(nargin == 9)
                AsrSamplesVolt = struct(...
                    'dcV1',  as(varargin{1}, V1), ...
                    'dcV2',  as(varargin{2}, V2), ...
                    'dcV3',  as(varargin{3}, V3), ...
                    'dcV12', as(varargin{4}, V12), ...
                    'dcV13', as(varargin{5}, V13), ...
                    'dcV23', as(varargin{6}, V23), ...
                    'acV12', as(varargin{7}, V12a), ...
                    'acV13', as(varargin{8}, V13a), ...
                    'acV23', as(varargin{9}, V23a));
            end
            %===================================================================
            import bicas.proc.L1L2.demuxer___UTEST.BLTS_src_array
            %===================================================================

            % ===========================
            % mux = 0, drUsing12 = [0, 1]
            % ===========================
            test(...
                0, false, {V1, V13, V23, V13a, V23a}, ...
                BLTS_src_array(...
                    {'DC single', 'DC diff', 'DC diff', 'AC diff', 'AC diff'}, ...
                    {[1], [1 3], [2 3], [1 3], [2 3]}), ...
                ASR_samples(1,1,1, 1,1,1, 1,1,1)...
            )
            test(...
                0, true, {V1, V12, V23, V12a, V23a}, ...
                BLTS_src_array(...
                    {'DC single', 'DC diff', 'DC diff', 'AC diff', 'AC diff'}, ...
                    {[1], [1 2], [2 3], [1 2], [2 3]}), ...
                ASR_samples(1,1,1, 1,1,1, 1,1,1)...
            )
        
            % ======================
            % mux = 1, drUsing12 = 0
            % ======================
            test(...
                1, false, {V2, V3, V23, V13a, V23a}, ...
                BLTS_src_array(...
                    {'DC single', 'DC single', 'DC diff', 'AC diff', 'AC diff'}, ...
                    {[2], [3], [2 3], [1 3], [2 3]} ...
                ), ...
                ASR_samples(0,1,1, 0,0,1, 1,1,1) ...
            )
        
            % ====================================
            % mux = 4 (calibration), drUsing12 = 1
            % ====================================
            test(...
                4, true, {V1, V2, V3, V12a, V23a}, ...
                BLTS_src_array(...
                    {'DC single', 'DC single', 'DC single', 'AC diff', 'AC diff'}, ...
                    {[1], [2], [3], [1 2], [2 3]} ...
                ), ...
                ASR_samples(1,1,1, 1,1,1, 1,1,1) ...
            )
        end
        
        
        
        function test_complement_ASR(testCase)
            
            % Local utility function.
            function assert_relation(A, B, C)
                % NOTE: Uses testCase. ==> Do not make static function.
                b = ~isnan(A) & ~isnan(B) & ~isnan(C);

                testCase.verifyEqual( A(b), B(b) + C(b) )
            end

        
        
            % NOTE: Only verifies the correct relationships between the return
            % results.
            function test(inputFieldsCa)
                A = bicas.proc.L1L2.demuxer.complement_ASR( struct(inputFieldsCa{:}) );

                % Test all possible relationsships.
                %
                % NOTE: Implicitly asserts that all fields are present.
                % NOTE: Must account for that some fields may be NaN, and
                %       therefore can not be checked against relations.
                assert_relation(A.dcV1,  A.dcV12, A.dcV2 )
                assert_relation(A.dcV1,  A.dcV13, A.dcV3 )
                assert_relation(A.dcV2,  A.dcV23, A.dcV3 )
                assert_relation(A.dcV13, A.dcV12, A.dcV23)    % DC. All diffs
                %
                assert_relation(A.acV13, A.acV12, A.acV23)    % AC. All diffs
            end
            %===================================================================

            % TODO: dlrUsing12

            test({'dcV1', 19, 'dcV12', 27, 'dcV23', 33,    'acV12', 54, 'acV23', 75});    % mux=0, dlrUsing12=1
            test({'dcV1', 19, 'dcV13', 27, 'dcV23', 33,    'acV13', 54, 'acV23', 75});    % mux=0, dlrUsing12=0
            test({'dcV2', 19, 'dcV3',  27, 'dcV23', 19-27, 'acV12', 54, 'acV23', 75});    % mux=1
            test({'dcV1', 2   'dcV2',  7,  'dcV3',  32,    'acV12', 74, 'acV23', 85});    % mux=4

        end
        
        
        
    end    % methods(Test)
        
        
    
    %########################
    %########################
    % PRIVATE STATIC METHODS
    %########################
    %########################
    methods(Static, Access=private)
        
        
        
        % Local utility function.
        % as = assign. Effectively implements ~ternary operator + constant (NaN).
        function V = as(b,V)
            assert(isscalar(b) && ismember(b, [0,1]))
            if b; V = V;
            else  V = NaN;
            end
        end
        
        
        
        function BltsSrcArray = BLTS_src_array(categoryArray, antennasArray)
            assert( numel(categoryArray) == numel(antennasArray) )

            for i =1:numel(categoryArray)
                BltsSrcArray(i) = bicas.proc.L1L2.PhysicalSignalSrcDest(...
                    categoryArray{i}, ...
                    antennasArray{i});
            end
        end
        
        
        
    end    % methods(Static, Access=private)

    
    
end
