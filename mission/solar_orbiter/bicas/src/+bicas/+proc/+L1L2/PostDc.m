%
% Store of all output information from processing (DC = Demuxing+Calibration)
% that is common for all L1/L1R-->L2 LFR+TDS processing.
%
% Class is admittedly small, but the corresponding data structure has
% historically been more larger(?), is partly an analogue to
% bicas.proc.L1L2.PreDc, and does indeed clarify the code.
%
%
% Author: Erik P G Johansson, IRF, Uppsala, Sweden
%
classdef PostDc
    
    %#####################
    %#####################
    % INSTANCE PROPERTIES
    %#####################
    %#####################
    % IMPLEMENTATION NOTE: Can not be immutable(!) since
    % bicas.proc.L1L2.process_quality_filter_L2() modifies it!
    properties   % (SetAccess=immutable)
        Zv
    end



    %#########################
    %#########################
    % PUBLIC INSTANCE METHODS
    %#########################
    %#########################
    methods(Access=public)

        function obj = PostDc(Zv)
            irf.assert.struct(Zv, ...
                {'DemuxerOutput', 'currentAAmpere'}, {'L2_QUALITY_BITMASK'});
            
            obj.Zv = Zv;

            bicas.proc.utils.assert_struct_num_fields_have_same_N_rows(obj.Zv);
        end

    end    % methods(Access=public)

end
