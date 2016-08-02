% Author: Erik P G Johansson, IRF-U, Uppsala, Sweden
% First created 2016-07-06
%
% Standalone "validation" tool for BIAS master CDFs, i.e. master (template) CDF files for any of the datasets that BIAS
% produces.
%
% The code tries to verify that multiple master cdfs satsify certain selected criteria instead of manually checking
% them. The code is primarily intended as a standalone tool, separate from the pipeline. It is neither intended to be
% "complete", nor to ever be constant and "finished".
% 
% This is useful when
% (1) manually creating master CDFs based on other CDFs (e.g. from other teams, other levels), or
% (2) the formal specification for master CDFs is modified/reinterpreted.
%
% ARGUMENTS:
% dir_path        : Directory path
% file_name_regex : Regular expression which specifies which files in the directory to include. ^ (beginning of string)
% and $ (end of string) are added automatically.
%
% RETURN VALUE:
% Cell array of dataobj for the validated CDF files.
%
function do_list = validate_BIAS_master_CDFs(dir_path, file_name_regex)
    
    irf('check_path');    
    irf.log('critical')      % Set log level.
    
    file_name_regex = ['^', file_name_regex, '$'];
    
    file_info_list = dir(dir_path);
    do_list = {};
    N_validated_files = 0;
    for i = 1:length(file_info_list)
        file_name = file_info_list(i).name;
        if ~isempty(regexp(file_name, file_name_regex))
            file_path = [dir_path, filesep, file_name];
            do_list{end+1} = validate_one_BIAS_master_CDF(file_path);
            N_validated_files = N_validated_files + 1;
        end
    end
    
    % Warning if no files were validated.
    % Useful if the caller types a regex with no matches.
    % Could otherwise be mistaken for no validation warnings.
    if N_validated_files == 0
        warning('No files to validate.')
    end
end



function do = validate_one_BIAS_master_CDF(file_path)
%
% PROPOSAL: Check for pad values? How?
% PROPOSAL: Some naming convention to exactly mimic CDF attributes?
% PROPOSAL: Print relevant values. Optional?
% PROPOSAL: Be able to check the BIAS HK master?
%
% Naming convention: CDF_* : Refers to something "concrete" in the CDF file.
%

    % NOTE: All initials capitalized.
    receiver_texts     = containers.Map({'LFR', 'TDS'}, {'Low Frequency Receiver', 'Time Domain Sampler'});
    mode_texts         = containers.Map({'SBM1', 'SBM2', 'SURV', 'LFM'}, {'Selective Burst Mode 1', 'Selective Burst Mode 2', 'Survey Mode', 'Low Frequency Mode'});
    data_product_texts = containers.Map({'CWF', 'SWF', 'RSWF'}, {'Continuous Waveform', 'Snapshot Waveform', 'Regular Snapshot Waveform'});
    
    fprintf('Validating %s\n', file_path)

    [parent_path, filename_base, filename_suffix] = fileparts(file_path);
    filename = [filename_base, filename_suffix];
    

            
    do = dataobj(file_path); 
    
    ga = do.GlobalAttributes;
    CDF_Descriptor       = ga.Descriptor{1};       % E.g.              "RPW-LFR-SURV-CWF-E>RPW Low Frequency Receiver Continuous Waveforms in survey mode. Electric component."
    CDF_DATASET_ID       = ga.DATASET_ID{1};       % E.g. "ROC-SGSE_L2S_RPW-LFR-SURV-CWF-E"
    CDF_SKELETON_PARENT  = ga.SKELETON_PARENT{1};  % E.g. "ROC-SGSE_L2S_RPW-LFR-SURV-CWF-E_V01.xlsx"
    CDF_Skeleton_version = ga.Skeleton_version{1}; % E.g.                                  "01"
    CDF_Level            = ga.Level{1};            % E.g.          "L2S>Level 2S"
    CDF_TEXT             = ga.TEXT{1};             % E.g. "This file contains RPW LFR level 2S continous waveform of electric data in survey mode."
    CDF_Logical_source             = ga.Logical_source{1};
    CDF_Logical_source_description = ga.Logical_source_description{1};  % E.g. "Solar Orbiter Radio/Plasma Wave, LFR L2S electric parameters"

    CDF_DATASET_ID_descriptor         = splitstr(CDF_DATASET_ID,            '_', 3, [3], 'Can not split DATASET_ID.');
    [receiver, mode, data_product]    = splitstr(CDF_DATASET_ID_descriptor, '-', 5, [2 3 4], 'Can not split DATASET_ID descriptor.');
    derived_CDF_DATASET_ID_descriptor = splitstr(CDF_Descriptor,            '>', 2, [1], 'Can not split Descriptor.');
    
    % receiver       % LFR, TDS.
    % mode           % SBM1, SBM2, SURV, LFM
    % data_product   % CWF, SWF, RSWF        
    derived_CDF_Descriptor = sprintf('%s>RPW %s %s in %s. Electric component.', ...
        derived_CDF_DATASET_ID_descriptor, ...
        get_map_value(receiver_texts, receiver), ...
        get_map_value(data_product_texts, data_product), ...
        lower(get_map_value(mode_texts, mode))   );    
    
    derived_CDF_TEXT = sprintf('This file contains RPW %s level 2S %s of electric data in %s.', ...
        receiver, ...
        lower(get_map_value(data_product_texts, data_product)), ...
        lower(get_map_value(mode_texts, mode)));
    
    derived_Logical_source_description = sprintf('Solar Orbiter Radio/Plasma Wave, %s L2S electric parameters', ...
        receiver);
        
        
    
    %=============================
    % Check filename
    %=============================
    % NOTE: Does not compare with skeleton parent, but does indirectly. Is that appropriate?
    derived_filename = sprintf('%s_V%s.cdf', CDF_DATASET_ID, CDF_Skeleton_version);  % NOTE: Forces lower case filename suffix.
    if ~strcmp(filename, derived_filename)
        validation_warning('Filename does not match global attributes.', filename, derived_filename)
    end
    
    
    
    %=============================
    % Check CDF global attributes
    %=============================
    validate_value(CDF_DATASET_ID_descriptor, derived_CDF_DATASET_ID_descriptor, 'DATASET_ID descriptor');
    validate_value(CDF_Descriptor, derived_CDF_Descriptor, 'Descriptor');
    
    % NOTE: Forces lower-case file suffix.
    derived_CDF_SKELETON_PARENT = sprintf('%s_V%s.xlsx', CDF_DATASET_ID, CDF_Skeleton_version);
    validate_value(CDF_SKELETON_PARENT, derived_CDF_SKELETON_PARENT, 'SKELETON_PARENT');    
    
    validate_value(CDF_Level, 'L2S>Level 2S', 'Level');    
    validate_value(CDF_TEXT, derived_CDF_TEXT, 'TEXT')
    
    validate_value(CDF_Logical_source, CDF_DATASET_ID, 'Logical_source')
    validate_value(CDF_Logical_source_description, derived_Logical_source_description, 'Logical_source_description');

    
    %======================================
    % Check absence/prescence of variables
    %======================================
    validate_zVariable_presence(do, 'ACQUISITION_TIME', 0)
    validate_zVariable_presence(do, 'ACQUISITION_TIME_UNITS', 1)
    validate_zVariable_presence(do, 'ACQUISITION_TIME_LABEL', 1)
    
    % Should ideally look at do.Variables(:,1) ?
    if any(~cellfun(@isempty, regexp(fieldnames(ga), 'PACKET_.*')))
        validation_warning('Found at least one global attribute PACKET_*.')
    end    
    
    % Check for misspelled variable attributes (i.e. variable NAMES).
    if isfield(do.VariableAttributes, 'CATEDESC')
        validation_warning('Found at least one instance of misspelled variable attribute CATEDESC.')
    end
    
    
    
    %=================================
    % Check "attribute" record values
    %=================================
    validate_value(strtrim(do.data.ACQUISITION_TIME_UNITS.data(2,:)), 's / 65536', 'ACQUISITION_TIME_UNITS(1,2) (trimmed)')
    
end



function validate_zVariable_presence(do, zvar_name, validate_nonempty)
    if ~any(strcmp(do.Variables(:,1), zvar_name))
        validation_warning(sprintf('zVariable %s is missing.', zvar_name));
        return
    end

    if validate_nonempty && isempty(do.data.(zvar_name).data)
        validation_warning(sprintf('Can not find any expected data/values for zVariable %s.', zvar_name));        
    end
    
    % Checks one variable attribute as a proxy for all variable attributes.
    if ~any(strcmp(do.VariableAttributes.FIELDNAM(:,1), zvar_name))
        validation_warning(sprintf('Can not find any FIELDNAM variable attribute for %s.', zvar_name));
    end
    
end



function validate_value(value, comparison_value, var_name)
    if ~strcmp(value, comparison_value)
        validation_warning(sprintf('%s does not match expected value.', var_name), value, comparison_value)
    end
end



% Get key in Map with default value rather than error if key does not exist.
function value = get_map_value(map, key)
    if map.isKey(key)
        value = map(key);
    else
        % Best to return something so that something is inserted into strings. Absence can be
        % ambiguous, e.g. when deriving Descriptor or TEXT.
        % NOTE: The caller might change the case.
        value = '<Unknown>';   
    end
end



% Split string using delimiter into substrings and return a specified subset of these. Verifies number of substrings.
% Useful for analyzing e.g. dataset IDs.
%
% str       : String to be split.
% delimiter : Delimiter used for splitting (duh!).
% N_parts   : The expected number of substrings (error if wrong).
% i_parts   : The substrings to be returned (numeric array).
% 
function varargout = splitstr(str, delimiter, N_parts, i_parts, val_warning_msg)
    if nargout ~= length(i_parts)
        error('Number of return values does not match input. This indicates a pure bug.')
    end

    parts = strsplit(str, delimiter);   % Returns ROW vector.
    if length(parts) ~= N_parts
        validation_warning(val_warning_msg, str)
        parts = [parts, cell(1, N_parts-length(parts))];   % Extend ROW vector.
    end
    
    for i_out = 1:length(i_parts)
        varargout{i_out} = parts{i_parts(i_out)};
    end
    
end



% Print standardized warning.
%
% ARGUMENTS: Number of arguments = 1, 2, or 3
% (msg)
% (msg, CDF_value_found)
% (msg, CDF_value_found, comparison_value)
function validation_warning(msg, varargin)
    LF = sprintf('\n');
    
    if length(varargin) >= 1
        msg = sprintf('%s\n           CDF value = "%s"', msg, varargin{1});
    end
    if length(varargin) == 2
        msg = sprintf('%s\n    Comparison value = "%s"', msg, varargin{2});
    end
        
    msg = strrep(msg, LF, [LF, '   ']);
    disp([' * ', msg])    
end

