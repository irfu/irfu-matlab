%
% Store of all input information for processing (DC = Demuxing+Calibration) that
% is common for all L1/L1R-->L2 LFR+TDS processing.
%
%
% Author: Erik P G Johansson, IRF, Uppsala, Sweden
%
classdef PreDc

    %#####################
    %#####################
    % INSTANCE PROPERTIES
    %#####################
    %#####################
    properties(SetAccess=immutable)
        Zv
        ZvFpa
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

        function obj = PreDc(Zv, ZvFpa, Ga, hasSnapshotFormat, isLfr, isTdsCwf)

            irf.assert.struct(Zv, ...
                {'Epoch', 'bltsSamplesTm', 'freqHz', 'nValidSamplesPerRecord', ...
                'iLsf', 'DIFF_GAIN', ...
                'MUX_SET', 'QUALITY_BITMASK', 'QUALITY_FLAG', 'SYNCHRO_FLAG', ...
                'DELTA_PLUS_MINUS', 'CALIBRATION_TABLE_INDEX', ...
                'ufv', 'lfrRx'}, ...
                {'BW'});
            irf.assert.struct(ZvFpa, ...
                {'HK_BIA_MODE_DIFF_PROBE'}, {})
            bicas.proc.utils.assert_struct_num_fields_have_same_N_rows(Zv);
            assert(isa(Zv.freqHz,         'double' ))
            assert(isa(hasSnapshotFormat, 'logical'))
            assert(isa(isLfr,             'logical'))
            assert(isa(isTdsCwf,          'logical'))

            obj.Zv                = Zv;
            obj.ZvFpa             = ZvFpa;
            obj.Ga                = Ga;
            obj.hasSnapshotFormat = hasSnapshotFormat;
            obj.isLfr             = isLfr;
            obj.isTdsCwf          = isTdsCwf;
        end

    end    % methods(Access=public)

end
