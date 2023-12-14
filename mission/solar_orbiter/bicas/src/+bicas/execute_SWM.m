%
% Execute a "S/W mode" as (indirectly) specified by the CLI arguments. This
% function should be agnostic of CLI syntax.
%
%
% ARGUMENTS AND RETURN VALUES
% ===========================
% SwmInfo
% InputFilePathMap  : containers.Map with
%    key   = prodFuncInputKey
%    value = Path to input file
% OutputFilePathMap : containers.Map with
%    key   = prodFuncOutputKey
%    value = Path to output file
%
%
% "BUGS"
% ======
% - Sets GlobalAttributes.Generation_date in local time (no fixed time zone,
%   e.g. UTC+0).
%
%
% Author: Erik P G Johansson, IRF, Uppsala, Sweden
% First created 2016-06-09
%
function execute_SWM(...
        SwmInfo, InputFilePathMap, OutputFilePathMap, ...
        masterCdfDir, rctDir, NsoTable, Bso, L)

    % TODO-NI: How verify dataset ID and dataset version against constants?
    %    NOTE: Need to read CDF first.
    %    NOTE: Need S/W mode.
    %
    % PROPOSAL: Verify output zVariables against master CDF zVariable dimensions
    %           (accessible with dataobj, even for zero records).
    %   PROPOSAL: function matchesMaster(DataObj, MasterDataobj)
    %       PRO: Want to use dataobj to avoid reading file (input dataset) twice.
    %
    % NOTE: Things that need to be done when writing PDV-->CDF
    %       Read master CDF file.
    %       Compare PDV variables with master CDF variables (only write a subset).
    %       Check variable types, sizes against master CDF.
    %       Write GlobalAttributes:
    %           Parents, Parent_version,
    %           Generation_date, Logical_file_id,
    %           Software_version, SPECTRAL_RANGE_MIN/-MAX (optional?), TIME_MIN/-MAX
    %       Write VariableAttributes: pad value? (if master CDF does not contain
    %           a correct value), SCALE_MIN/-MAX
    %
    % PROPOSAL: Print variable statistics also for zVariables which are created with fill values.
    %   NOTE: These do not use NaN, but fill values.



    % ASSERTION: Check that all input & output dataset paths (strings) are
    % unique.
    % NOTE: Manually entering CLI argument, or copy-pasting BICAS call, can
    % easily lead to reusing the same path by mistake, and e.g. overwriting an
    % input file.
    datasetFileList = [InputFilePathMap.values(), OutputFilePathMap.values()];
    assert(numel(unique(datasetFileList)) == numel(datasetFileList), ...
        'BICAS:CLISyntax', ...
        ['Input and output dataset paths are not all unique.', ...
        ' This hints of a manual mistake', ...
        ' in the CLI arguments in call to BICAS.'])



    %=============================
    % READ CDFs
    % ---------
    % Iterate over all INPUT CDFs
    %=============================
    InputDatasetsMap = containers.Map();
    for i = 1:length(SwmInfo.inputsList)
        prodFuncInputKey = SwmInfo.inputsList(i).prodFuncInputKey;
        inputFilePath    = InputFilePathMap(prodFuncInputKey);

        %=======================
        % Read dataset CDF file
        %=======================
        InputDataset = bicas.read_dataset_CDF(inputFilePath, Bso, L);
        
        InputDatasetsMap(prodFuncInputKey) = InputDataset;



        %============================================
        % ASSERTIONS: Check global attributes values
        %============================================
        % NOTE: Can not use
        % bicas.proc.utils.assert_struct_num_fields_have_same_N_rows(Zv) since
        % not all zVariables have same number of records.
        % Ex: Metadata such as ACQUISITION_TIME_UNITS.
        if ~isfield(InputDataset.Ga, 'Dataset_ID')
            error('BICAS:Assertion:DatasetFormat', ...
                ['Input dataset does not contain (any accepted variation', ...
                ' of) the global attribute Dataset_ID.\n    File: "%s"'], ...
                inputFilePath)
        end
        cdfDsi = InputDataset.Ga.Dataset_ID{1};

        if ~strcmp(cdfDsi, SwmInfo.inputsList(i).dsi)
            [settingValue, settingKey] = Bso.get_fv(...
                'INPUT_CDF.GA_DSI_MISMATCH_POLICY');
            anomalyDescrMsg = sprintf(...
                ['The input CDF dataset''s stated DSI does', ...
                ' not match value expected from the S/W mode.\n', ...
                '    File: %s\n', ...
                '    Global attribute InputDataset.Ga.Dataset_ID{1} : "%s"\n', ...
                '    Expected value:                                : "%s"\n'], ...
                inputFilePath, cdfDsi, SwmInfo.inputsList(i).dsi);
            bicas.default_anomaly_handling(L, ...
                settingValue, settingKey, ...
                'E+W+illegal', anomalyDescrMsg, 'BICAS:DatasetFormat')
        end
    end



    %==============
    % PROCESS DATA
    %==============
    % NPEF = No Processing, Empty File
    [settingNpefValue, settingNpefKey] = Bso.get_fv(...
        'OUTPUT_CDF.NO_PROCESSING_EMPTY_FILE');
    if ~settingNpefValue
        %==========================
        % CALL PRODUCTION FUNCTION
        %==========================
        OutputDatasetsMap = SwmInfo.Swmp.production_function(InputDatasetsMap, rctDir, NsoTable, Bso, L);
    else
        L.logf('warning', ...
            'Disabled processing due to setting %s.', settingNpefKey)
        
        % IMPLEMENTATION NOTE: Needed for passing assertion. Maybe to be
        % considered a hack?!
        OutputDatasetsMap = struct(...
            'keys', {{SwmInfo.outputsList.prodFuncOutputKey}});
    end



    %==============================
    % WRITE CDFs
    % ----------
    % Iterate over all OUTPUT CDFs
    %==============================
    % ASSERTION: The output datasets generated by processing equal the datasets
    % required by the s/w mode.
    irf.assert.castring_sets_equal(...
        OutputDatasetsMap.keys, ...
        {SwmInfo.outputsList.prodFuncOutputKey});
    %
    for iOutputCdf = 1:length(SwmInfo.outputsList)
        OutputInfo = SwmInfo.outputsList(iOutputCdf);

        prodFuncOutputKey = OutputInfo.prodFuncOutputKey;
        outputFilePath    = OutputFilePathMap(prodFuncOutputKey);



        %========================
        % Write dataset CDF file
        %========================
        masterCdfPath = fullfile(...
            masterCdfDir, ...
            bicas.get_master_CDF_filename(...
                OutputInfo.dsi, ...
                OutputInfo.skeletonVersion));

        if ~settingNpefValue
            % CASE: Nominal
            OutputDataset = OutputDatasetsMap(OutputInfo.prodFuncOutputKey);

            ZvsSubset = OutputDataset.Zv;

            GaSubset = bicas.derive_output_dataset_GAs(...
                InputDatasetsMap, OutputDataset, ...
                irf.fs.get_name(outputFilePath), OutputInfo.dsi, ...
                Bso, L);
        else
            % CASE: No processing.
            ZvsSubset = [];
            GaSubset  = struct();
        end
        bicas.write_dataset_CDF( ...
            ZvsSubset, GaSubset, outputFilePath, masterCdfPath, ...
            Bso, L );
    end



end   % execute_SWM
