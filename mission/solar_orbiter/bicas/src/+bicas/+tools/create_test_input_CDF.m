% Author: Erik P G Johansson, IRF-U, Uppsala, Sweden
% First created 2016-08-05
%
% Utility for modifying a genuine dataset CDF file to contain varied (random) data useful for manually comparing BICAS'
% input and output CDF files. This is useful for testing.
%
% NOTE: EM2_CAL_BIAS_SWEEP_LFR_CONF1_1M_2016-04-15_Run1__e1d0a9a__CNES/ROC-SGSE_HK_RPW-BIA_e1d0a9a_CNE_V01.cdf has
% DATASET_ID = ' ', i.e. can not determine the dataset ID.
%
function create_test_input_CDF(input_file_path, output_file_path)
%
% NOTE: Certain values must be constrained. Ex: MUX.
% NOTE: Certain zVariables should maybe never be set to the fill value. Ex: MUX.
% NOTE: BICAS takes advantage of continuous sequences with constant settings.
%       ==> Wants sequences (of random length) of constant values.
%     QUESTION: How handle multiple such variables which should be constant together?
% NOTE: MUX can be set both in SCI CDF file and HK CDF file. If randomizing them, then they will be inconsistent if
%       BICAS ever checks for that.
%
% PROPOSAL: Set to pad values?
% PROPOSAL: Take dataset version into account.
% PROPOSAL: Remove from BICAS/irfu-matlab.
%
% QUESTION: How handle files without DATASET_ID? 1) How determine the DATASET_ID. Set the DATASET_ID in the output
% files?



% Randomizer settings.
BULK_SETTINGS = struct('fraction_fill_values', 0.1, 'min', 0, 'max', 9, 'record_transition_probability', 1);
MUX_SETTINGS  = struct('fraction_fill_values', 0.0, 'min', 0, 'max', 4, 'record_transition_probability', 0.4);



[out, info] = spdfcdfread(input_file_path, 'Structure', 1, 'KeepEpochAsIs', 1);



%================================================
% Determine which variables to randomize and how
%================================================
datasetId = info.GlobalAttributes.DATASET_ID{1};
if isempty(strtrim(datasetId))                      
    % Ugly temporary "fix" for incomplete input files.
    fprintf('No DATASET_ID for %s.\n   Replacing with .Logical_source .\n', input_file_path);
    datasetId = info.GlobalAttributes.Logical_source{1};
end

zVars_to_randomize = struct('zVar_name', {}, 'randomizer_settings', {});
switch(datasetId)
    case 'ROC-SGSE_HK_RPW-BIA'
        zVars_to_randomize(end+1) = struct('zVar_name', 'HK_BIA_MODE_MUX_SET',  'randomizer_settings', MUX_SETTINGS);
    case 'ROC-SGSE_L2R_RPW-LFR-SURV-CWF'
        zVars_to_randomize(end+1) = struct('zVar_name', 'POTENTIAL',  'randomizer_settings', BULK_SETTINGS);
        zVars_to_randomize(end+1) = struct('zVar_name', 'ELECTRICAL', 'randomizer_settings', BULK_SETTINGS);        
    case 'ROC-SGSE_L2R_RPW-LFR-SURV-SWF'
        zVars_to_randomize(end+1) = struct('zVar_name', 'POTENTIAL',  'randomizer_settings', BULK_SETTINGS);
        zVars_to_randomize(end+1) = struct('zVar_name', 'ELECTRICAL', 'randomizer_settings', BULK_SETTINGS);
    otherwise
        error('Can not handle this file. datasetId = "%s".', datasetId)
end



%================================================
% Randomize variables
%================================================
for i = 1:length(zVars_to_randomize)
    zVar_name           = zVars_to_randomize(i).zVar_name;
    randomizer_settings = zVars_to_randomize(i).randomizer_settings;
    
    j = find(strcmp(info.VariableAttributes.FIELDNAM(:,1), zVar_name));
    
    fill_val = info.VariableAttributes.FILLVAL{...
        find(...
            strcmp(info.VariableAttributes.FILLVAL(:,1), ...
            zVar_name)), ...
        2};
    
    out(j).Data = replace_with_random_data(out(j).Data, fill_val, randomizer_settings);
end



bicas.utils.write_CDF(output_file_path, out, info, 'fill_empty')

end



% Create random data variable of the same size and the same MATLAB class.
%
% settings:
%   .min, .max                     : inclusive range
%   .record_transition_probability : Probability that consecutive values (over records) are separately random values.
%                                    (Note that they can still be identical, if generated separately)
function random_data = replace_with_random_data(actual_data, fill_value, settings)

% Set to random numerical values.
random_data = floor(settings.min + (settings.max - settings.min+1) * rand(size(actual_data)));

% Make consecutive records similar or not.
% PROPOSAL: Combine with generating random numbers.
i_value = 1;
for i = 1:size(random_data, 1)
    if rand < settings.record_transition_probability;
        i_value = i;
    end
    random_data(i,:) = random_data(i_value,:);
end

% Replace random components with fill value.
random_data(rand(size(actual_data)) < settings.fraction_fill_values) = fill_value;

% Set MATLAB class to the original one. This is important when writing to the CDF file.
random_data = cast(random_data, class(actual_data));

end



