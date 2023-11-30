%
% matlab.unittest automatic test code for bicas.proc.L1L2.demuxer.
%
% Could be improved but unsure how much is meaningful. Seems to complicated.
%
%
% Author: Erik P G Johansson, IRF, Uppsala, Sweden
% First created 2021-09-08, using older test code.
%
classdef demuxer___UTEST < matlab.unittest.TestCase
    
    
    
    methods(Static, Access=private)
        
        function [channelIdStringsCa, testSamplesCa] = create_channel_test_data(sampleSize)
            % (iChannel, 1) = ASID ID string.
            % (iChannel, 2) = Sample value (that is consistent with other
            %                 channels).
            A = bicas.proc.L1L2.AntennaSignalId.C;
            TEST_DATA_CA = { ...
                A.DC_V1.s,  10; ...
                A.DC_V2.s,  11; ...
                A.DC_V3.s,  13; ...
                A.DC_V12.s, 10-11; ...
                A.DC_V13.s, 10-13; ...
                A.DC_V23.s, 11-13; ...
                A.AC_V12.s, 45-56; ...
                A.AC_V13.s, 45-69; ...
                A.AC_V23.s, 56-69 ...
            };
            
            % Multiply the sample values with matrix to test multiple records
            % with "snapshots" (SPR>1).
            TEST_DATA_CA(:, 2) = cellfun(@(x) (x * ones(sampleSize)), TEST_DATA_CA(:, 2), 'UniformOutput', false);
            
            channelIdStringsCa = TEST_DATA_CA(:, 1);
            testSamplesCa      = TEST_DATA_CA(:, 2);
        
        end

    end
    


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

            A = bicas.proc.L1L2.AntennaSignalId.C;
            R = bicas.proc.L1L2.Routing.C;

            % =========
            % Test data
            % =========
            TEST_DATA_UNKNOWN = [999];   % Data from unknown source.
            SAMPLES_SIZE      = [3,2];
            
            nRows = SAMPLES_SIZE(1);
            [channelIdStringsCa, testSamplesCa] = bicas.proc.L1L2.demuxer___UTEST.create_channel_test_data(SAMPLES_SIZE);
            
            AsidTestSamplesSrm = containers.Map(channelIdStringsCa, testSamplesCa);
            
            
            
            % Test any BDM.
            function test(bdmFloatNan, dlrFloatNan, bltsSamplesAVolt, ExpRoutingArray, ExpAsrSamplesAVoltSrm)
                assert(numel(ExpRoutingArray) == 5)
                
                dlrFpa = bicas.utils.FPArray.floatNan2logical(dlrFloatNan);
                bdmFpa = bicas.utils.FPArray.floatNan2int(bdmFloatNan, 'uint8');
                
                % CALL FUNCTIONS
                ActRoutingArray       = bicas.proc.L1L2.demuxer.get_routings(...
                    bdmFpa, dlrFpa);
                ActAsrSamplesAVoltSrm = bicas.proc.L1L2.demuxer.calibrated_BLTSs_to_ASRs(...
                    [ActRoutingArray.sdid], bltsSamplesAVolt);
                
                % ASSERTIONS
                testCase.assertEqual(ActRoutingArray, ExpRoutingArray)
                testCase.assertTrue(ActAsrSamplesAVoltSrm == ExpAsrSamplesAVoltSrm)
            end

            

            % Function for testing BDM 0-4. All those all those map (route) ASR
            % to ASR (no GNS, no 2.5V REF, no "unknown", no "nowhere").
            function test_BDM01234(bdmFloatNan, dlrFloatNan, ExpRoutingArray, ExpAsrSamplesAVoltSrm)
                assert(isa(ExpAsrSamplesAVoltSrm, 'bicas.utils.SameRowsMap'))
                assert(ismember(bdmFloatNan, [0:4]))
                
                % Autogenerate bltsSamplesCa (test argument) using
                % ExpRoutingArray (only possible for BDM 0-4.
                tempBltsSamplesAVolt = gen_BLTS_samples(ExpRoutingArray);
                
                test(bdmFloatNan, dlrFloatNan, tempBltsSamplesAVolt, ExpRoutingArray, ExpAsrSamplesAVoltSrm)
            end
            
            
            
            function tempBltsSamplesAVolt = gen_BLTS_samples(RoutingArray)
                tempBltsSamplesAVolt = zeros(SAMPLES_SIZE);
                
                for i = 1:numel(RoutingArray)
                    routing = RoutingArray(i);
                    if routing.ssid.is_ASR()
                        tempBltsSamplesAVolt(:, :, i) = AsidTestSamplesSrm(routing.ssid.asid.s);
                    else
                        tempBltsSamplesAVolt(:, :, i) = TEST_DATA_UNKNOWN;
                    end
                end
            end



            % Create samples per ASID using constants. Arguments determine for
            % which ASIDs samples should be NaN instead of data.
            %
            % varargin{i} == 0 or 1. Determines whether constant or NaN will
            % be used.
            function AsSrm = get_ASR_samples(varargin)
                assert(nargin == 9)
                
                % Define which varargin{i} corresponds to which ASID.
                ASID_CA = {...
                    A.DC_V1,  A.DC_V2,  A.DC_V3,  ...
                    A.DC_V12, A.DC_V13, A.DC_V23, ...
                    A.AC_V12, A.AC_V13, A.AC_V23 ...
                };
                AsSrm = bicas.utils.SameRowsMap('char', nRows, 'EMPTY');
                
                for iAsid = 1:9
                    asidName = ASID_CA{iAsid}.s;

                    samplesAVolt = AsidTestSamplesSrm(asidName);
                    if ~varargin{iAsid}
                        samplesAVolt = samplesAVolt * NaN;
                    end
                    AsSrm.add(asidName, samplesAVolt);
                end
            end



            % ==========================
            % bdm = 0, dlr = [0, 1, NaN]
            % ==========================
            test_BDM01234(...
                0, 0, ...
                [R.DC_V1, R.DC_V12, R.DC_V23, R.AC_V12, R.AC_V23], ...
                get_ASR_samples(1,1,1, 1,1,1, 1,1,1)...
            )
            test_BDM01234(...
                0, 1, ...
                [R.DC_V1, R.DC_V13, R.DC_V23, R.AC_V13, R.AC_V23], ...
                get_ASR_samples(1,1,1, 1,1,1, 1,1,1)...
            )
            test_BDM01234(...
                0, NaN, ...
                [R.DC_V1, R.UNKNOWN_TO_NOWHERE, R.DC_V23, R.UNKNOWN_TO_NOWHERE, R.AC_V23], ...
                get_ASR_samples(1,0,0, 0,0,1, 0,0,1)...
            )

            % =======================
            % bdm = 1, dlr = [0, NaN]
            % =======================
            test_BDM01234(...
                1, 1, ...
                [R.DC_V2, R.DC_V3, R.DC_V23, R.AC_V13, R.AC_V23], ...
                get_ASR_samples(0,1,1, 0,0,1, 1,1,1) ...
            )
            test_BDM01234(...
                1, NaN, ...
                [R.DC_V2, R.DC_V3, R.DC_V23, R.UNKNOWN_TO_NOWHERE, R.AC_V23], ...
                get_ASR_samples(0,1,1, 0,0,1, 0,0,1) ...
            )

            % ==============================
            % bdm = 4 (calibration), dlr = 1
            % ==============================
            test_BDM01234(...
                4, 0, ...
                [R.DC_V1, R.DC_V2, R.DC_V3, R.AC_V12, R.AC_V23], ...
                get_ASR_samples(1,1,1, 1,1,1, 1,1,1) ...
            )


            % ==============
            % BDM = 5, DLR 1
            % ==============
            bltsSamplesAVolt(:, :, 1) = AsidTestSamplesSrm(A.DC_V1.s);
            bltsSamplesAVolt(:, :, 2) = AsidTestSamplesSrm(A.DC_V2.s);
            bltsSamplesAVolt(:, :, 3) = AsidTestSamplesSrm(A.DC_V3.s);
            bltsSamplesAVolt(:, :, 4) = AsidTestSamplesSrm(A.AC_V13.s);
            bltsSamplesAVolt(:, :, 5) = AsidTestSamplesSrm(A.AC_V23.s);
            test(5, 1, ...
                bltsSamplesAVolt, ...
                [R.REF25V_TO_DC_V1, R.REF25V_TO_DC_V2, R.REF25V_TO_DC_V3, R.AC_V13, R.AC_V23], ...
                get_ASR_samples(1,1,1, 1,1,1, 1,1,1))
            
            % ==============
            % BDM = 6, DLR 0
            % ==============
            bltsSamplesAVolt(:, :, 1) = AsidTestSamplesSrm(A.DC_V1.s);
            bltsSamplesAVolt(:, :, 2) = AsidTestSamplesSrm(A.DC_V2.s);
            bltsSamplesAVolt(:, :, 3) = AsidTestSamplesSrm(A.DC_V3.s);
            bltsSamplesAVolt(:, :, 4) = AsidTestSamplesSrm(A.AC_V12.s);
            bltsSamplesAVolt(:, :, 5) = AsidTestSamplesSrm(A.AC_V23.s);
            test(6, 0, ...
                bltsSamplesAVolt, ...
                [R.GND_TO_DC_V1, R.GND_TO_DC_V2, R.GND_TO_DC_V3, R.AC_V12, R.AC_V23], ...
                get_ASR_samples(1,1,1, 1,1,1, 1,1,1))
            
            % ==============
            % BDM = Unknwon, DLR Unknown
            % ==============
            bltsSamplesAVolt(:, :, 1) = AsidTestSamplesSrm(A.DC_V1.s);
            bltsSamplesAVolt(:, :, 2) = AsidTestSamplesSrm(A.DC_V2.s);
            bltsSamplesAVolt(:, :, 3) = AsidTestSamplesSrm(A.DC_V3.s);
            bltsSamplesAVolt(:, :, 4) = AsidTestSamplesSrm(A.AC_V12.s);
            bltsSamplesAVolt(:, :, 5) = AsidTestSamplesSrm(A.AC_V23.s);
            test(NaN, NaN, ...
                bltsSamplesAVolt, ...
                [R.UNKNOWN_TO_NOWHERE, R.UNKNOWN_TO_NOWHERE, R.UNKNOWN_TO_NOWHERE, R.UNKNOWN_TO_NOWHERE, R.AC_V23], ...
                get_ASR_samples(0,0,0, 0,0,0, 0,0,1))
        end
        
        
        
        function test_complement_ASR(testCase)
            
            C = bicas.proc.L1L2.AntennaSignalId.C;


            
            % Local utility function.
            function assert_relation(A, B, C)
                % NOTE: Uses testCase. ==> Do not make static function.
                b = ~isnan(A) & ~isnan(B) & ~isnan(C);

                testCase.verifyEqual( A(b), B(b) + C(b) )
            end


        
            % NOTE: Only verifies the correct relationships between the return
            % results. Does not verify entire return results.
            % BUG/NOTE: Will fail if function returns NaN when it should not!
            function test(inputFieldsCa)
                nRows = size(inputFieldsCa{2}, 1);
                AsSrm = bicas.utils.SameRowsMap('char', nRows, 'EMPTY');
                for i = 1:(numel(inputFieldsCa)/2)
                    asid  = inputFieldsCa{2*i-1};
                    value = inputFieldsCa{2*i  };
                    %AsMap(key) = value;
                    AsSrm.add(asid.s, value)
                end
                
                % RUN FUNCTION TO BE TESTED
                bicas.proc.L1L2.demuxer.complement_ASR(AsSrm);
                ActAsSrm = AsSrm;

                % Test all possible relationsships.
                %
                % NOTE: Implicitly asserts that all fields are present.
                % NOTE: Must account for that some fields may be NaN, and
                %       therefore can not be checked against relations.
                assert_relation(ActAsSrm(C.DC_V1.s),  ActAsSrm(C.DC_V12.s), ActAsSrm(C.DC_V2.s))
                assert_relation(ActAsSrm(C.DC_V2.s),  ActAsSrm(C.DC_V23.s), ActAsSrm(C.DC_V3.s))
                assert_relation(ActAsSrm(C.DC_V1.s),  ActAsSrm(C.DC_V13.s), ActAsSrm(C.DC_V3.s))
                
                % DC. All diffs
                assert_relation(ActAsSrm(C.DC_V13.s), ActAsSrm(C.DC_V12.s), ActAsSrm(C.DC_V23.s))
                
                % AC. All diffs
                assert_relation(ActAsSrm(C.AC_V13.s), ActAsSrm(C.AC_V12.s), ActAsSrm(C.AC_V23.s))
            end



            test({C.DC_V1, 19, C.DC_V12, 27, C.DC_V23, 33,    C.AC_V12, 54, C.AC_V23, 75});    % bdm=0, dlr=0
            test({C.DC_V1, 19, C.DC_V13, 27, C.DC_V23, 33,    C.AC_V13, 54, C.AC_V23, 75});    % bdm=0, dlr=1
            test({C.DC_V2, 19, C.DC_V3,  27, C.DC_V23, 19-27, C.AC_V12, 54, C.AC_V23, 75});    % bdm=1
            test({C.DC_V1,  2, C.DC_V2,   7, C.DC_V3,  32,    C.AC_V12, 74, C.AC_V23, 85});    % bdm=4

        end
        
        
        
    end    % methods(Test)



end
