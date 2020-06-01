%
% Utility for modifying a genuine dataset CDF file to contain varied (random) data useful for manually comparing BICAS'
% input and output CDF files. This is useful for testing.
%
% NOTE: Unmaintained. Likely needs updating.
%
% NOTE: EM2_CAL_BIAS_SWEEP_LFR_CONF1_1M_2016-04-15_Run1__e1d0a9a__CNES/ROC-SGSE_HK_RPW-BIA_e1d0a9a_CNE_V01.cdf has
% DATASET_ID = ' ', i.e. can not determine the dataset ID.
%
% Author: Erik P G Johansson, IRF, Uppsala, Sweden
% First created 2016-08-05
%
function create_test_input_CDF(inputFilePath, outputFilePath)
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
    BULK_SETTINGS = struct('fractionFillValues', 0.1, 'min', 0, 'max', 9, 'recordTransitionProbability', 1);
    MUX_SETTINGS  = struct('fractionFillValues', 0.0, 'min', 0, 'max', 4, 'recordTransitionProbability', 0.4);



    [out, info] = spdfcdfread(inputFilePath, 'Structure', 1, 'KeepEpochAsIs', 1);



    %================================================
    % Determine which variables to randomize and how
    %================================================
    datasetId = info.GlobalAttributes.Dataset_ID{1};
    if isempty(strtrim(datasetId))
        % Ugly temporary "fix" for incomplete input files.
        fprintf('No global attribute "Dataset_ID" for %s.\n   Replacing with "Logical_source" .\n', inputFilePath);
        datasetId = info.GlobalAttributes.Logical_source{1};
    end
    
    zVarsToRandomize = struct('zvName', {}, 'RandomizerSettings', {});
    switch(datasetId)
        case 'ROC-SGSE_HK_RPW-BIA'
            zVarsToRandomize(end+1) = struct('zvName', 'HK_BIA_MODE_MUX_SET',  'RandomizerSettings', MUX_SETTINGS);
        case 'ROC-SGSE_L2R_RPW-LFR-SURV-CWF'
            zVarsToRandomize(end+1) = struct('zvName', 'POTENTIAL',  'RandomizerSettings', BULK_SETTINGS);
            zVarsToRandomize(end+1) = struct('zvName', 'ELECTRICAL', 'RandomizerSettings', BULK_SETTINGS);
        case 'ROC-SGSE_L2R_RPW-LFR-SURV-SWF'
            zVarsToRandomize(end+1) = struct('zvName', 'POTENTIAL',  'RandomizerSettings', BULK_SETTINGS);
            zVarsToRandomize(end+1) = struct('zvName', 'ELECTRICAL', 'RandomizerSettings', BULK_SETTINGS);
        otherwise
            error('Can not handle this file. datasetId = "%s".', datasetId)
    end
    
    
    
    %================================================
    % Randomize variables
    %================================================
    for i = 1:length(zVarsToRandomize)
        zvName             = zVarsToRandomize(i).zvName;
        RandomizerSettings = zVarsToRandomize(i).RandomizerSettings;
        
        j = find(strcmp(info.VariableAttributes.FIELDNAM(:,1), zvName));
        
        fillValue = info.VariableAttributes.FILLVAL{...
            find(...
            strcmp(info.VariableAttributes.FILLVAL(:,1), ...
            zvName)), ...
            2};
        
        out(j).Data = replace_with_random_data(out(j).Data, fillValue, RandomizerSettings);
    end
    
    
    
    bicas.utils.write_CDF(outputFilePath, out, info, 'fill_empty')
    
end



% Create random data variable of the same size and the same MATLAB class.
%
% settings:
%   .min, .max                   : Inclusive range
%   .recordTransitionProbability : Probability that consecutive values (over records) are separately random values.
%                                  (Note that they can still be identical, if generated separately)
function randomData = replace_with_random_data(actualData, fillValue, Settings)
    
    % Set to random numerical values.
    randomData = floor(Settings.min + (Settings.max - Settings.min+1) * rand(size(actualData)));
    
    % Make consecutive records similar or not.
    % PROPOSAL: Combine with generating random numbers.
    iValue = 1;
    for i = 1:size(randomData, 1)
        if rand < Settings.recordTransitionProbability
            iValue = i;
        end
        randomData(i,:) = randomData(iValue,:);
    end
    
    % Replace random components with fill value.
    randomData(rand(size(actualData)) < Settings.fractionFillValues) = fillValue;
    
    % Set MATLAB class to the original one. This is important when writing to the CDF file.
    randomData = cast(randomData, class(actualData));
    
end
