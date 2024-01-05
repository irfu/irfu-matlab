%
% Author: Erik P G Johansson, IRF, Uppsala, Sweden
%
classdef RctTypeTdsRswf < bicas.proc.L1L2.cal.rct.RctType



    %#####################
    %#####################
    % INSTANCE PROPERTIES
    %#####################
    %#####################
    properties(Constant, GetAccess=public)
        filenameRegexpSettingKey = 'PROCESSING.RCT_REGEXP.TDS-LFM-RSWF';
    end



    %#######################
    %#######################
    % PUBLIC STATIC METHODS
    %#######################
    %#######################
    methods(Static)



        % NOTE: The TDS RSWF RCT contains ITFs, not FTFs.
        %
        function RctData = read_RCT(filePath)

            Do = dataobj(filePath);

            try
                % ASSUMPTION: Exactly 1 CDF record.
                % IMPLEMENTATION NOTE: Does not want to rely one dataobj special
                % behaviour for 1 record case.
                % ==> Remove leading singleton dimensions (just to be
                % format-tolerant), many assertions.
                freqsHz  = shiftdim(Do.data.CALIBRATION_FREQUENCY.data);   % 1x512 --> 512x1
                amplIvpt = shiftdim(Do.data.CALIBRATION_AMPLITUDE.data);   % 3x512 --> 3x512
                phaseDeg = shiftdim(Do.data.CALIBRATION_PHASE.data);       % 3x512 --> 3x512

                % ASSERTIONS: Check CDF array sizes, no change in format.
                nFreqs = irf.assert.sizes(...
                    freqsHz,  [-1,  1], ...
                    amplIvpt, [ 3, -1], ...
                    phaseDeg, [ 3, -1]);
                assert(nFreqs >= bicas.proc.L1L2.cal.rct.RctType.TF_TABLE_MIN_LENGTH)

                for iBlts = 1:3
                    % NOTE: RCT contains ITF, not FTF.
                    ItfIvpt = irf.utils.tabulated_transform(...
                        freqsHz * 2*pi, ...
                        amplIvpt(        iBlts, :), ...
                        deg2rad(phaseDeg(iBlts, :)));

                    % ASSERTION: INVERTED TF
                    assert(~ItfIvpt.toward_zero_at_high_freq(), ...
                        ['TDS RSWF transfer function appears to go toward', ...
                        ' zero at high frequencies. Has it not been', ...
                        ' inverted/made backward in time, i.e. does it not', ...
                        ' describe physical output-to-physical input?'])

                    ItfIvptList{iBlts} = ItfIvpt;
                end

                RctData = [];
                RctData.ItfIvptList = ItfIvptList;

            catch Exc1
                Exc2 = MException(...
                    'BICAS:FailedToReadInterpretRCT', ...
                    'Error when interpreting calibration file (TDS team''s', ...
                        ' LFM RSWF RCT for BIAS/BICAS) "%s"', ...
                    filePath);
                Exc2 = Exc2.addCause(Exc1);
                throw(Exc2);
            end
        end



        function RctData2 = modify_RCT_data(RctData1)
            RctData2 = [];

            % Modify tabulated TDS-RSWF TFs.
            for iBlts = 1:numel(RctData1.ItfIvptList)
                % NOTE: Overwriting.

                ItfRctIvpt = RctData1.ItfIvptList{iBlts};

                % Store tabulated ITF EXACTLY AS THEY ARE in the RCT (before
                % modification).
                % NOTE: Struct field does not need to be
                % pre-initialized/pre-allocated.
                RctData2.ItfRctIvptCa{iBlts} = ItfRctIvpt;

                % MODIFY __tabulated__ ITF
                % (Does NOT wrap function handle in function handle.)
                ItfModifIvpt = bicas.proc.L1L2.cal.utils.extrapolate_tabulated_TF_to_zero_Hz(ItfRctIvpt);

                % MODIFY tabulated ITF --> function ITF
                %
                % NOTE: Use zero outside of tabulated frequencies (beyond
                % already made extrapolation). TDS-RSWF data requires
                % extrapolation.
                VALUE_OUTSIDE_TABLE = 0;
                itfModifIvpt = @(omegaRps) (bicas.proc.L1L2.cal.utils.eval_tabulated_TF(...
                    ItfModifIvpt, omegaRps, VALUE_OUTSIDE_TABLE));

                RctData2.itfModifIvptCa{iBlts} = itfModifIvpt;

            end

        end



        function log_RCT(RctData, L)
            % TODO: Log tabulated TFs.

            FREQ_HZ = 0;

            for iBlts = 1:3
                itfNamePrefix = sprintf('TDS RSWF, BLTS/BIAS_%i, ITF', iBlts);

                bicas.proc.L1L2.cal.utils.log_TF_tabulated(...
                    bicas.proc.L1L2.cal.rct.RctType.RCT_DATA_LL, ...
                    sprintf('%s (as in RCT)', itfNamePrefix), ...
                    RctData.ItfRctIvptCa{iBlts}, ...
                    L);

                bicas.proc.L1L2.cal.utils.log_TF_function_handle(...
                    bicas.proc.L1L2.cal.rct.RctType.RCT_DATA_LL, ...
                    sprintf('%s (modif., interp.)', itfNamePrefix), ...
                    'ivolt/TM unit', FREQ_HZ, RctData.itfModifIvptCa{iBlts}, L)
            end
        end



    end    % methods(Static)



end