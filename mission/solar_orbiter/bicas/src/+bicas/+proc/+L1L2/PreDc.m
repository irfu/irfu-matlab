%
% Store of all input information for processing (DC = Demuxing+Calibration) that
% is common for all L1/L1R-->L2 LFR+TDS processing.
%
%
% Author: Erik P G Johansson, IRF, Uppsala, Sweden
%
classdef PreDc
    % PROPOSAL: Proper abbreviation for class PreDc and PostDc.
    %   PRDC, PODC
    %   before/after
    %   calibration
    %   demultiplexing
    %   data
    %
    % IMPLEMENTATION NOTE: bltsSamplesTm
    % ==================================
    % PreDc always represents (has variables/elements for) all five BLTS's,
    % despite that only three are used at any given time. The channels not used
    % are set to NaN. Which ones are actually used can be switched at any given
    % time due to lfrRx changing.
    %
    % This handling does in principle differ from the handling of other "data
    % parameters" (BW, LSF, isTdsCwf, freqHz, hasSnapshotFormat, etc.). lfrRx
    % determines which channels contain actual data. The BLTS index into the
    % array is used to determine whether one has DC diff or AC diff data.
    %
    % PROPOSAL: Refactor code to only represent 3 BLTS's simultaneously, which
    %           all (nominally) have data. Store iBlts for the diff channels so
    %           that one can use the demultiplexer routings.
    %   CON: Has to split these channels/arrays to set output dataset ZVs which
    %        are greater in number.
    %       CON: Should be easy to implement.
    %   CON: This increases the assumptions which the ~demultiplexer code, PreDc
    %        make.
    %       CON: No, it does not. Which?
    %           CON: Always 3 channels. BLTS 2-3/4-5 are different from BLTS 1
    %                and can change meaning due to Rx.
    %   PROPOSAL: Split BLTSs: One DC single (BLTS 1) and two diffs which are
    %             either DC or AC (BLTS 2-3 or 4-5).
    %       CON: Can not iterate over all 5 BLTSs.
    %       CON: BLTS 2-3 ARE NOT ALWAYS DIFFS!!!
    
    
    
    %#####################
    %#####################
    % INSTANCE PROPERTIES
    %#####################
    %#####################
    properties(SetAccess=immutable)
        Zv
        Ga
        hasSnapshotFormat
        isLfr
        isTdsCwf
    end



    %#########################
    %#########################
    % PUBLIC INSTANCE METHODS
    %#########################
    %#########################
    methods(Access=public)

        function obj = PreDc(Zv, Ga, hasSnapshotFormat, isLfr, isTdsCwf)

            irf.assert.struct(Zv, ...
                {'Epoch', 'bltsSamplesTm', 'freqHz', 'nValidSamplesPerRecord', ...
                'bdmFpa', 'biasHighGainFpa', 'dlrFpa', ...
                'iLsf', ...
                'QUALITY_BITMASK', 'QUALITY_FLAG', 'SYNCHRO_FLAG', ...
                'DELTA_PLUS_MINUS', 'CALIBRATION_TABLE_INDEX', ...
                'ufv', 'lfrRx'}, ...
                {'BW'});
            bicas.proc.utils.assert_struct_num_fields_have_same_N_rows(Zv);
            assert(size(Zv.bltsSamplesTm, 3) == 5)
            assert(isa(Zv.freqHz,         'double' ))
            assert(isa(hasSnapshotFormat, 'logical'))
            assert(isa(isLfr,             'logical'))
            assert(isa(isTdsCwf,          'logical'))

            obj.Zv                = Zv;
            obj.Ga                = Ga;
            obj.hasSnapshotFormat = hasSnapshotFormat;
            obj.isLfr             = isLfr;
            obj.isTdsCwf          = isTdsCwf;
        end

    end    % methods(Access=public)

end
