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
% file_name_regex : Regular expression which match the entire filenames of the files to validate in the directory. The
%                   function adds ^ (beginning of string) and $ (end of string) are added automatically to the regular expression.
%
% RETURN VALUE:
% varargout : Optional single return value. Cell array of "dataobj" objects for the validated CDF files. This is useful
%             for manually inspecting the CDF files to investigate what is wrong.
%
function [varargout] = validate_BIAS_master_CDFs(dir_path, file_name_regex)

    %test_code ; return
    
    
    irf('check_path');    
    
    if ~ischar(file_name_regex)
        error('"file_name_regex" is not a string.')
    end
    file_name_regex = ['^', file_name_regex, '$'];
    
    
    
    file_info_list = dir(dir_path);
    do_list = {};
    N_validated_files = 0;
    for i = 1:length(file_info_list)
        file_name = file_info_list(i).name;
        if ~isempty(regexp(file_name, file_name_regex, 'once'))
            file_path = fullfile(dir_path, file_name);
            do_list{end+1} = validate_one_BIAS_master_CDF(file_path);    % do = dataobj
            N_validated_files = N_validated_files + 1;
        end
    end
    
    % Warning if no files were validated.
    % Useful if the caller types a regex with no matches.
    % This could otherwise be mistaken for no validation warnings.
    if N_validated_files == 0
        warning('No files to validate.')
    end
    
    
    
    if nargout == 0
        % Do nothing
    elseif nargout == 1
        varargout = {do_list};
    else
        error('Illegal number of return values.')
    end
end



% Informal test code
function test_code
    validate_value(123, 234, 'Xxxx', 0)
end



function do = validate_one_BIAS_master_CDF(file_path)
%
% PROPOSAL: Some naming convention to exactly mimic CDF attributes?
%
% PROPOSAL: Split validation into part applicable to all datasets, and part specifically for BICAS.
%    PRO: Can compare the compliance of master CDF files from other groups.
%    NOTE: Complicated(?) for semiconditional zvariable attributes for which choices have been made for BICAS datasets.
%
% PROPOSAL: Check that zVariables are correct w.r.t. sample/snapshot per record.
%    NOTE: Need list of zVariables which change with sample/snapshot per record: V, E, EAC, BIAS1/2/3.
%    PROPOSAL: Check that zVariables use the right dimension for snapshots.
%
% PROPOSAL: Check that contains the right set of zVariables (no more, no less).
% PROPOSAL: Check that the mandatory (by spec) zVariables are present.
% PROPOSAL: Check that every zVariable contains the mandatory (by spec) zVariable attributes.
% PROPOSAL: Check that all mandatory global attributes are present.
%
% PROPOSAL: Check that IBIAS1-3 are analogous/almost identical (within a dataset).
% PROPOSAL: Have a master CDF for the master CDF files?!!
% PROPOSAL: Some means of checking that large subsets of different CDF files are identical.
%    NOTE: Breaks present model of comparing CDF files independently.
%    QUESTION: Of a set of CDF files, compare which datasets with which?
%        PROPOSAL: Let one be special that serves as "master" for the other master files?
%    QUESTION: How define subsets that are be expected to be equal?
%        PROPOSAL: All attributes for a given zVariable.
%        PROPOSAL: Specific subset of global attributes.
%
% Naming convention: CDF_* : Refers to something "concrete" in the CDF file.
%
    PRINT_LATEST_MODS = 1;

    % NOTE: All initials capitalized.
    receiver_texts     = containers.Map({'LFR', 'TDS'}, {'Low Frequency Receiver', 'Time Domain Sampler'});
    mode_texts         = containers.Map({'SBM1', 'SBM2', 'SURV', 'LFM'}, {'Selective Burst Mode 1', 'Selective Burst Mode 2', 'Survey Mode', 'Low Frequency Mode'});
    data_product_texts = containers.Map({'CWF', 'SWF', 'RSWF'}, {'Continuous Waveform', 'Snapshot Waveform', 'Regular Snapshot Waveform'});
    
    fprintf('Validating %s\n', file_path)

    [parent_path, filename_base, filename_suffix] = fileparts(file_path);
    filename = [filename_base, filename_suffix];
    

            
    do = dataobj(file_path);    % do = dataobj   
    zVar_names = do.Variables(:,1);    
    ga = do.GlobalAttributes;
    
    if PRINT_LATEST_MODS
        % Print latest MOD record.
        % Useful for checking that all MODs have been updated.
        fprintf('   Last MODS = "%s"\n', do.GlobalAttributes.MODS{end});
    end
    
    
    % Seems values are always strings in cdfs, although ROC-TST-GSE-NTT-00017-LES thinks it should be a number.
    % Resolution, Bin_location should probably be numbers.
    % Probably a deficiency in xlsx2skt.
    %EPOCH_zVariable_name = 'EPOCH';
    EPOCH_zVariable_name = 'Epoch';   % Prescribed by the ROC-TST-GSE-NTT-00017-LES, iss2rev0.
    zVar_EPOCH_attributes = {...
        'FIELDNAM',     'Epoch'; ...
        'CATDESC',      'Default time'; ...
        'FILLVAL',      '9999-12-31T23:59:59.999999999'; ...
        'LABLAXIS',     'Epoch'; ...
        'UNITS',        'ns'; ...
        'VALIDMIN',     '2000-01-01T00:00:00.000000000'; ...
        'VALIDMAX',     '2050-12-31T23:59:59.999000000'; ...
        'SCALEMIN',     '1990-01-01T00:00:00.000000000'; ...
        'SCALEMAX',     '2050-12-31T23:59:59.999000000'; ...
        'VAR_TYPE',     'support_data'; ...
        'SCALETYP',     'linear'; ...      % Not misspelled according to ISTP/IACG Variable Attributes
        'TIME_BASE',    'Spacecraft clock'; ...
        'TIME_SCALE',   'Spacecraft clock'; ...
        'REFERENCE_POSITION', 'MEB GSE'; ...
        'Resolution',   '15258'; ...
        'Bin_location', '0.5'; ...
        % 'VAR_NOTES',    'Primary time used a reference in the file.'   % Prescribed in docs, but misspelled.
        'VAR_NOTES',    'Primary time used as a reference in the file.'
    };
    for i = 1:size(zVar_EPOCH_attributes,1)
        validate_zVariable_attribute_value(do, EPOCH_zVariable_name, zVar_EPOCH_attributes{i,1}, zVar_EPOCH_attributes{i,2})
    end
    
    
    
    %=====================================
    % Validate global attributes presence
    %=====================================
    % Docs specifies "Acknowledgment", but ISTP specifies "Acknowledgement" (different spelling).
    % NOTE: Incomplete list of mandatory global attributes (too long).
    MANDATORY_L2S_GLOBAL_ATTRIBUTES = {'ACCESS_FORMAT', 'ACCESS_URL', 'Acknowledgement', 'Calibration_version', 'Dataset_ID', 'Data_type', 'Data_version', 'Descriptor', 'Discipline', 'File_naming_convention', 'FILE_UUID', 'Generated_by', 'HTTP_LINK', 'Instrument_type', 'Level', 'LINK_TEXT'};
    for i = 1:length(MANDATORY_L2S_GLOBAL_ATTRIBUTES)
        validate_global_attribute_presence(do, MANDATORY_L2S_GLOBAL_ATTRIBUTES{i})
    end
    
    
    
    CDF_Descriptor       = ga.Descriptor{1};       % E.g.              "RPW-LFR-SURV-CWF-E>RPW Low Frequency Receiver Continuous Waveforms in survey mode. Electric component."
    CDF_DATASET_ID       = ga.Dataset_ID{1};       % E.g. "ROC-SGSE_L2S_RPW-LFR-SURV-CWF-E"
    CDF_SKELETON_PARENT  = ga.SKELETON_PARENT{1};  % E.g. "ROC-SGSE_L2S_RPW-LFR-SURV-CWF-E_V01.xlsx"
    CDF_Skeleton_version = ga.Skeleton_version{1}; % E.g.                                  "01"
    CDF_Level            = ga.Level{1};            % E.g.          "L2S>Level 2S"
    CDF_TEXT             = ga.TEXT{1};             % E.g. "This file contains RPW LFR level 2S continous waveform of electric data in survey mode."
    CDF_Logical_source             = ga.Logical_source{1};
    CDF_Logical_source_description = ga.Logical_source_description{1};  % E.g. "Solar Orbiter Radio/Plasma Wave, LFR L2S electric parameters"

    CDF_DATASET_ID_descriptor         = splitstr(CDF_DATASET_ID,            '_', 3, [3], 'Can not split DATASET_ID.');
    [receiver, mode, data_product]    = splitstr(CDF_DATASET_ID_descriptor, '-', 5, [2 3 4], 'Can not split DATASET_ID descriptor.');
    derived_CDF_DATASET_ID_descriptor = splitstr(CDF_Descriptor,            '>', 2, [1], 'Can not split Descriptor.');
    
    % Receivers      % LFR, TDS.
    % Modes          % SBM1, SBM2, SURV, LFM
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
    % NOTE: Forces lower case filename suffix.    
    % NOTE: Does not compare with skeleton parent, but does indirectly. Is that appropriate?
    derived_filename = sprintf('%s_V%s.cdf', CDF_DATASET_ID, CDF_Skeleton_version);
    if ~strcmp(filename, derived_filename)
        validation_value('Filename', filename, derived_filename)
    end
    

    
    %=======================================
    % Validate CDF global attributes values
    %=======================================
    validate_value(CDF_DATASET_ID_descriptor, derived_CDF_DATASET_ID_descriptor, 'DATASET_ID descriptor');
    validate_value(CDF_Descriptor,            derived_CDF_Descriptor,            'Descriptor');
    
    % NOTE: Forces lower-case file suffix.
    derived_CDF_SKELETON_PARENT = sprintf('%s_V%s.xlsx', CDF_DATASET_ID, CDF_Skeleton_version);
    validate_value(CDF_SKELETON_PARENT, derived_CDF_SKELETON_PARENT, 'SKELETON_PARENT');    
    
    validate_value(CDF_Level, 'L2S>Level 2S',   'Level');    
    validate_value(CDF_TEXT,  derived_CDF_TEXT, 'TEXT')
    
    validate_value(CDF_Logical_source,             CDF_DATASET_ID,                     'Logical_source')
    validate_value(CDF_Logical_source_description, derived_Logical_source_description, 'Logical_source_description');

    % Should ideally look at do.Variables(:,1) ?
    if any(~cellfun(@isempty, regexp(fieldnames(ga), 'PACKET_.*')))
        validation_warning('Found at least one global attribute PACKET_*.')
    end
    
    
    
    switch(CDF_DATASET_ID)
        case {'ROC-SGSE_L2S_RPW-LFR-SURV-SWF-E', 'ROC-SGSE_L2S_RPW-TDS-LFM-RSWF-E'}
            snapshot = 1;
        case {'ROC-SGSE_L2S_RPW-LFR-SBM1-CWF-E', 'ROC-SGSE_L2S_RPW-LFR-SBM2-CWF-E', 'ROC-SGSE_L2S_RPW-LFR-SURV-CWF-E', 'ROC-SGSE_L2S_RPW-TDS-LFM-CWF-E'}
            snapshot = 0;
        otherwise
            validation_warning('Does not recognize DATASET_ID="%s"', CDF_DATASET_ID)
            return
    end
    
    
    %=========================================
    % Validate absence/presence of zVariables
    %=========================================
    
    % Disabled since it appears that TIME_SYNCHRO_FLAG is not mandatory for L2S after all (misunderstood documentation?).
    %validate_zVariable_presence(do, 'TIME_SYNCHRO_FLAG', 0)
    validate_zVariable_presence(do, 'ACQUISITION_TIME', 0)
    validate_zVariable_presence(do, 'ACQUISITION_TIME_UNITS', 1)
    validate_zVariable_presence(do, 'ACQUISITION_TIME_LABEL', 1)
    validate_zVariable_presence(do, 'QUALITY_BITMASK', 0)
    validate_zVariable_presence(do, 'QUALITY_FLAG', 0)
    validate_zVariable_presence(do, 'IBIAS1', 0)
    validate_zVariable_presence(do, 'IBIAS2', 0)
    validate_zVariable_presence(do, 'IBIAS3', 0)
    if snapshot
        validate_zVariable_presence(do, 'F_SAMPLE', 0)
    end
    
    %==================================
    % Validate zVariable record values
    %==================================
    validate_value(strtrim(do.data.ACQUISITION_TIME_UNITS.data(2,:)), 's / 65536', 'ACQUISITION_TIME_UNITS(1,2) (trimmed)')

    %===============================
    % Validate zVariable attributes
    %===============================
    % Check for presence of misspelled zvariable attribute.
    if isfield(do.VariableAttributes, 'CATEDESC')
        validation_warning('Found at least one instance of misspelled variable attribute CATEDESC.')
    end
    if isfield(do.VariableAttributes, 'SCALETYPE')
        validation_warning('Found at least one instance of misspelled variable attribute SCALETYPE.')
    end


    validate_zVariable_attribute_value(do, 'DELTA_PLUS_MINUS', 'CATDESC', 'Time between sample timestamp and beginning/end of integration. Total integration time is twice this value.')

    % ROC-TST-GSESPC-00017-LES, iss01, rev 03 requires a number of zVar attributes, many of them conditionally.
    % SRDB_PARAM_ID, SRDB_ENUM_ID - mandatory.
    % LABLAXIS - conditional
    % DISPLAY_TYPE - conditional
    % FORMAT or FORM_PTR
    % 
    % Disabled since there are so many errors (misunderstood documentation?). Awaiting validation at ROC instead.    
    if 0
        for i = 1:length(zVar_names)
            zVar_name = zVar_names{i};
            validate_zVariable_attribute_presence(do, zVar_name, 'FIELDNAM')
            validate_zVariable_attribute_presence(do, zVar_name, 'CATDESC')
            validate_zVariable_attribute_presence(do, zVar_name, 'DEPEND_0')   % Not for zvar Epoch.
            validate_zVariable_attribute_presence(do, zVar_name, 'DISPLAY_TYPE')
            validate_zVariable_attribute_presence(do, zVar_name, 'FILLVAL')
            validate_zVariable_attribute_presence(do, zVar_name, 'FORMAT')
            validate_zVariable_attribute_presence(do, zVar_name, 'VAR_NOTES')
            validate_zVariable_attribute_presence(do, zVar_name, 'LABLAXIS')
        end
    end



    %==================
    % Check pad values
    %==================
    % NOTE: F_SAMPLE only present in snapshot datasets.
    zVAR_PAD_VALUES_CHECK = {'V', 'E', 'EAC', 'IBIAS1', 'IBIAS2', 'IBIAS3'};
    for i = 1:length(zVAR_PAD_VALUES_CHECK)
        zVar_name = zVAR_PAD_VALUES_CHECK{i};
        
        validate_zVariable_presence(do, zVar_name, 0)
        
        i_zVar = find(strcmp(do.Variables(:,1), zVar_name));
        if numel(i_zVar) == 1
            pad_value = do.Variables{i_zVar, 9};
            validate_value(pad_value, -1e30, sprintf('pad value for zVariable "%s"', zVar_name), 0.0);
        end
    end
end



function validate_global_attribute_presence(do, attribute_name)
    if ~isfield(do.GlobalAttributes, attribute_name)
        validation_warning('Can not find global attribute "%s"', attribute_name)
    end
end



function validate_zVariable_attribute_presence(do, zVar_name, attribute_name)
    if isfield(do.VariableAttributes, attribute_name)
        i = find(strcmp(zVar_name, do.VariableAttributes.(attribute_name){:, 1}));
        if ~isempty(i)
            return
        end
    end
    validation_warning('Can not find zVariable attribute %s:%s', zVar_name, attribute_name)
end



function value = get_zVariable_attr(do, zVar_name, attribute_name)
% QUESTION: How handle non-presence of attribute? Error? Warning? (Could be second time for same warning.)
%    If warning (attribute not found), then how should validation proceed?

    if isfield(do.VariableAttributes, attribute_name)
        i = find(strcmp(zVar_name, do.VariableAttributes.(attribute_name)(:, 1)));
        if ~isempty(i)
            value = do.VariableAttributes.(attribute_name){i, 2};
            return
        end        
    end
    validation_warning('Can not find zVariable attribute %s:%s', zVar_name, attribute_name)
    value = [];
end



function validate_zVariable_attribute_value(do, zVar_name, attribute_name, comparison_value)
    % QUESTION: How handle comparing different data types (MATLAB classes)?
    value = get_zVariable_attr(do, zVar_name, attribute_name);
    if ~isempty(value)
        validate_value(value, comparison_value, [zVar_name, ':', attribute_name])
    end
end



% do : dataobj object
function validate_zVariable_presence(do, zvar_name, validate_nonempty)
    % PROPOSAL: Separate validate empty.
    %   PRO: Iterate over Zvars that should be present.
    if ~any(strcmp(do.Variables(:,1), zvar_name))
        validation_warning('Can not find zVariable %s.', zvar_name);
        return
    end

    if validate_nonempty && isempty(do.data.(zvar_name).data)
        validation_warning('Can not find any expected data/values for zVariable %s.', zvar_name)        
    end
    
    % Checks one variable attribute as a proxy for all variable attributes.
    if ~any(strcmp(do.VariableAttributes.FIELDNAM(:,1), zvar_name))
        validation_warning('Can not find any FIELDNAM variable attribute for %s.', zvar_name);
    end
    
end



%======================================================================================
% Validate a specific value by comparing it to another value that should be identical.
% 
% NOTE: Special case for comparing zero.
% NOTE: Approximate comparison of numbers which can permit error even for permitted error zero (example: compare pad
% value -1e30 with -1e30). This is an unintended effect of a lack of precision.
%
% varargin : For strings - No argument
%            For numbers - Natural logarithm of permitted error.
%======================================================================================
function validate_value(value, comparison_value, var_name, varargin)
    msg = sprintf('%s does not match expected value.', var_name);

    if ischar(value) && ischar(comparison_value)
        % CASE: Comparing strings
        
        if numel(varargin) ~= 0
            error('Wrong number of arguments')
        end
        if ~strcmp(value, comparison_value)
            validate_value___warning(msg, sprintf('"%s"', value), sprintf('"%s"', comparison_value))
        end
        
    elseif isnumeric(value) && isnumeric(comparison_value)
        % CASE: Comparing numbers
        
        if numel(varargin) ~= 1
            error('Wrong number of arguments')
        end
        permitted_error = varargin{1};
        
        if value == 0 && comparison_value == 0
            return
        elseif value ~= 0 && comparison_value ~= 0
            if abs(log(value/comparison_value)) > permitted_error   % Can not handle comparing zero.
                validate_value___warning(msg, sprintf('%f', value), sprintf('%f', comparison_value))
            end
        else
            validate_value___warning(msg, sprintf('%f', value), sprintf('%f', comparison_value))
        end
        
    else
        error('Can not compare these types of variables.')
    end
    
    %----------------------------------------------------------------------------------------------
    % Print standardized warning when comparing two values.
    % Printed line 1-N : Message string
    % Printed line N+1 : Value
    % Printed line N+2 : Comparison value
    %
    % NOTE: Does not automatically quote values so that the caller can chose quoted/unquoted (e.g. strings/numbers).
    function validate_value___warning(msg, value, comparison_value)
        % NOTE: The value may be the filename, i.e. not something INSIDE the file. The phrases should be chosen
        % accordingly.
        % NOTE: Extra indentation in addition to what validation_warning adds.
        % NOTE: Example of long var_name: "ACQUISITION_TIME_UNITS(1,2) (trimmed)" (37 characters)
        msg = [msg, sprintf('\n   %-37s = %s',            var_name, value)];
        msg = [msg, sprintf('\n   Comparison value                      = %s', comparison_value)];
        validation_warning(msg)
    end
    %----------------------------------------------------------------------------------------------
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



%======================================================================================
% Split string using delimiter into substrings and return a specified subset of these. Verifies number of substrings.
% Useful for analyzing e.g. dataset IDs.
%
% str       : String to be split.
% delimiter : Delimiter used for splitting.
% N_parts   : The expected number of substrings (error if wrong).
% i_parts   : The substrings to be returned (numeric array).
%======================================================================================
function varargout = splitstr(str, delimiter, N_parts, i_parts, val_warning_msg)
    if nargout ~= length(i_parts)
        error('Number of return values does not match input. This indicates a pure bug.')
    end

    parts = strsplit(str, delimiter);   % Returns ROW vector.
    if length(parts) ~= N_parts
        validation_warning([val_warning_msg, '. str="', str, '"'])
        parts = [parts, cell(1, N_parts-length(parts))];   % Extend ROW vector.
    end
    
    for i_out = 1:length(i_parts)
        varargout{i_out} = parts{i_parts(i_out)};
    end
    
end



%======================================================================================
% Print warning on standardized format.
%
% Adds its own indentation.
% Can handle multiline messages.
%======================================================================================
function validation_warning(msg, varargin)
    LF = sprintf('\n');
    msg = sprintf(msg, varargin{:});
    msg = strrep(msg, LF, [LF, '   ']);   % Necessary for indentation.
    disp([' * ', msg])    
end



