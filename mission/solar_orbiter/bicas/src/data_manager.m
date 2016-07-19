% Author: Erik P G Johansson, IRF-U, Uppsala, Sweden
% First created 2016-06-10
%
% "Data manager". Does the actual processing of datasets.
%
classdef data_manager
%###################################################################################################
% PROPOSAL: Use a containers.Map object for the actual data instead of separate properties.
%    PRO: Can use more arbitrary string identifiers, use e.g. dash.
%       PRO: Can recycle arbitrary string identifiers from elsewhere, e.g. dataset IDs
%          Want to use the same string identifiers as the execute_sw_mode parameter uses.
%       NOTE: Intermediate data types would not have dataset IDs.
%    CON: Want to use shorter strings anyway.
%    CON: Less clear which variables exists.
%    PRO: Can probably recycle code because of it.
%       CON: Which?
%
% PROPOSAL: Use properties to store actual data. Field name = data type identifier string.
%    CON: Unwise to share same namespace with other properties.
%       CON: Can distinguish them with prefix.
%           CON: That's an ugly hack...
%
% PROPOSAL: Use other class name that implices processing. "processing_manager"?
%
% QUESTION: Should this class at all handle reading and writing cdf _FILES_?
%    Should that not be done in some separate layer?!
%    PROPOSAL: Only work with input/output data types that correspond tightly to cdf files.
%       CON: Would work against the idea of combining code with the S/W descriptor info.
%
% TODO: Set access restrictions on properties, methods.
%
%###################################################################################################


    properties   %(Access=private)
        % Convention: properties is empty = property has not been set yet.
        
        input_cdfs = containers.Map({'cdf_BIAS_HK', 'cdf_LFR_CWF'}, {[], []});       
        output_cdf_LFR_CWF = [];
        
        %input_cdf_LFR_SWF
        %output_cdf_LFR_SWF
        
        test_property = [123];
    end
    
    %###############################################################################################

    methods (Access=public)
        
        % NOTE/PROBLEM: Can not prevent from setting data types/properties which are not "input_*".
        function set_input_cdf(obj, data_type, file_path)
            irf.log('n', sprintf('data_type=%s: file_path=%s', data_type, file_path))    % NOTE: irf.log adds the method name.
            
            data = dataobj(file_path);            
            
            if ~obj.input_cdfs.isKey(data_type)
                errorp(ERROR_CODES.SW_MODE_PROCESSING_ERROR, 'Can not set non-existent data type data_type "%s".', data_type);
            end
            obj.input_cdfs(data_type) = data;
        end
        
        % IMPLEMENTATION NTOE: Replacement for using get.* property access methods for all input
        % data types. If we implemented separate methods, they would be bacially identical anyway.
        function data = get_input_cdf_data(obj, data_type)
            global ERROR_CODES

            data = obj.input_cdfs(data_type);   % NOTE: Will give error for non-existant data types.
            if ~isempty(data)
                return
            end
            errorp(ERROR_CODES.SW_MODE_PROCESSING_ERROR, 'Can not find data type "%s".', data_type)
        end
        
        function write_output_LFR_CWF_cdf(obj, file_path)
            data = obj.output_cdf_LFR_CWF;
            
        end
        
        %==================================================
        
        % TEST
        %function test_method(obj, data_type, x)
        %    obj.(data_type) = x;
        %end
    end
    
    %###############################################################################################

    methods
        function data = get.output_cdf_LFR_CWF(obj)
            if ~isempty(output_cdf_LFR_CWF)
                data = output_cdf_LFR_CWF;
                return
            end
            
            input_cdf_BIAS_HK = get_input_cdf_data('input_cdf_BIAS_HK');
            input_cdf_LFR_CWF = get_input_cdf_data('input_cdf_LFR_CWF');
        end
        
        % TEST
        %function r = get.test_property(obj)
        %    r = obj.test_property;
        %end
    end
    
end
