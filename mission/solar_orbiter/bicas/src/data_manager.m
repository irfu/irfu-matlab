% Author: Erik P G Johansson, IRF-U, Uppsala, Sweden
% First created 2016-06-10
%
% "Data manager". Does the actual processing of datasets.
%
classdef data_manager
%###################################################################################################
% PROPOSAL: Use other class name that implices processing. "processing_manager"?
%
% QUESTION: Should this class at all handle reading and writing cdf _FILES_?
%    Should that not be done in some separate layer?!
%    PROPOSAL: Only work with input/output data types that correspond tightly to cdf files.
%       CON: Would work against the idea of combining code with the S/W descriptor info.
%
% QUESTION: Should this code have some kind of facility for handling multiple dataset versions in
% principle, eventhough not formally required? Could be useful.
%    NOTE: Could have multiple internal datatypes, remnants from earlier implementations.
%    
%
% PROPOSAL: Use a containers.Map object for the actual data instead of separate properties.
%    PRO: Can use more arbitrary string identifiers, use e.g. dash.
%       PRO: Can recycle arbitrary string identifiers from elsewhere, e.g. dataset IDs
%          Want to use the same string identifiers as the execute_sw_mode parameter uses.
%       NOTE: Intermediate data types would not have dataset IDs.
%    PRO: Can probably recycle code because of it.
%       CON: Which?
%    CON: Want to use shorter strings anyway.
%    CON: Less clear which variables exists.
%
% PROPOSAL: Use properties to store actual data. Field name = data type identifier string.
%    CON: Unwise to share same namespace with other properties.
%       CON: Can distinguish them with prefix.
%           CON: That's an ugly hack...
%
% PROPOSAL: private class instance variable storing global ERROR_CODES.
%
% PROPOSAL: Different name for internal data variable. data, datas, data_types, internal_data,
% data_network, data_chain, chain_data, process_data, processing_data. A one-letter name?
%
% PROPOSAL: Implement master cdfs as process data types?!!
%    NOTE: Should fit with S/W description.
%
% NOTE: Both BIAS HK and LFR SURV CWF contain MUX data (only LFR has one timestamp per snapshot). True also for other input datasets?
%
% PROPOSAL: Move out the reading of input cdf files?! Of reading master cdfs?!!
%###################################################################################################


    properties(Access=private)
        % Convention: properties is empty = property has not been set yet.
        
        
        % NOTE: Only set keys here. This makes sure that no disallowed keys are ever used.
        % For data corresponding to cdfs: Use dataset IDs. Possibly with version numbers.
        process_data = containers.Map('KeyType', 'char', 'ValueType', 'any');
        
        test_property = [123];
    end
    
    %###############################################################################################

    methods (Access=public)
        
        % Constructor
        function obj = data_manager()
            INPUT_CDFS  = {'ROC-SGSE_L2R_RPW-LFR-SURV-CWF',   'ROC-SGSE_L2R_RPW-LFR-SURV-SWF', 'ROC-SGSE_HK_RPW-BIA'};
            OUTPUT_CDFS = {'ROC-SGSE_L2S_RPW-LFR-SURV-CWF-E', 'ROC-SGSE_L2S_RPW-LFR-SURV-SWF-E'};
            ALL_DATA_TYPES = {INPUT_CDFS{:}, OUTPUT_CDFS{:}};
            
            for data_type = ALL_DATA_TYPES
                obj.process_data(data_type{1}) = [];
            end
        end

        
        
        % NOTE: Not preventing from setting data types/properties which are not cdf files.
        function set_input_cdf(obj, data_type, file_path)
            irf.log('n', sprintf('data_type=%s: file_path=%s', data_type, file_path))    % NOTE: irf.log adds the method name.
            global ERROR_CODES
            
            data = dataobj(file_path);
            obj.set_process_data(data_type, data);
        end
        
        
        
        % Like get_process_data but tries to fill in values recursively using other process data.
        function data = get_data(obj, data_type)
            global ERROR_CODES
            
            data = get_process_data(obj, data_type);   % Implicitly checks if the data type exists.
            if ~isempty(data)
                return
            end
            
            C = bicas_constants.get_constants();
            output = bicas_constants.get_cdf_output_constants({data_type});
            switch(data_type)
                case 'ROC-SGSE_L2S_RPW-LFR-SURV-CWF-E'
                    % TODO: Read master cdf!
                    input_BIAS_HK = get_data(obj, 'ROC-SGSE_HK_RPW-BIA');
                    input_LFR_CWF = get_data(obj, 'ROC-SGSE_L2R_RPW-LFR-SURV-CWF');
                    
                    BIAS_HK_ACQUISITION_TIME = input_BIAS_HK.data.ACQUISITION_TIME.data;
                    BIAS_HK_mux              = input_BIAS_HK.data.HK_BIA_MODE_MUX_SET.data;
                    LFR_ACQUISITION_TIME     = input_LFR_CWF.data.ACQUISITION_TIME.data;
                    LFR_BIAS_mux             = input_LFR_CWF.data.BIAS_MODE_MUX_SET.data;
                    
                    
                    
                    data = [];
                    
                otherwise
                    errorp(E, 'Can not produce this data type, "%s".', data_type)
            end
        end
        
    end   % methods
    
    %###############################################################################################

    methods(Access=private)
        
        % Return what is in the internal data structure, without trying to fill it with data.
        function data = get_process_data(obj, data_type)
            global ERROR_CODES
            
            if ~obj.process_data.isKey(data_type)
                errorp(ERROR_CODES.SW_MODE_PROCESSING_ERROR, 'There is no such process data type, "%s".', data_type);
            end
            data = obj.process_data(data_type)
        end
        
        % Set the value associated with a key in the process_data.
        % Will make sure that an unused key is not set.
        function set_process_data(obj, data_type, data)
            global ERROR_CODES
            
            if ~obj.process_data.isKey(data_type)
                errorp(ERROR_CODES.SW_MODE_PROCESSING_ERROR, 'There is no such process data type, "%s".', data_type);
            elseif ~isempty(obj.process_data(data_type))
                errorp(ERROR_CODES.SW_MODE_PROCESSING_ERROR, 'There is aldready process data for that data type, "%s".', data_type);
            end
            obj.process_data(data_type) = data;
        end
        
    end   % methods
    
end
