%
% Function for dynamically deriving global attributes for a specific output
% dataset given the global attributes of multiple input datasets.
%
% NOTE: Some of the global attribute values determined here are
%   (1) unique for this particular output dataset,
%   (2) common for all output datasets for the current s/w mode,
%   (3) common for alla output datasets.
%
%
% SOOP_TYPE, Datetime, OBS_ID
% ===========================
% XB on RCS telecon 2020-09-17: SOOP_TYPE, Datetime, OBS_ID should be taken from
% L1 (not HK, unless implicit that it should).
% --
% Global attributes Datetime, OBS_ID, SOOP_TYPE appear to be present in BICAS
% input L1R datasets, CURRENT datasets, and BIAS HK datasets. Not true for old
% SBM1 datasets (at least).
% Exception: OBS_ID is not in BIAS HK. /2020-09-17
%
%
% ARGUMENTS
% =========
% InputDatasetsMap
%       NOTE: This function does not really need all of InputDatasetsMap as
%       input (contains zVars) but the function uses that input argument since
%       it is easily accessible where this function is called.
% OutputDataset
%       Struct from processing with fields
%           .Ga.(globAttrName)
%               Subset of global attribute values that should be used.
%           .Zv.(zvName)
%               zVariables.
% outputFilename
%       Output dataset filename. Could potentially be used for deriving
%       Glob.attrs. Datetime (time interval string), Data_version,
%       (DATASET_ID).
%
%
% RETURN VALUE
% ============
% OutGaSubset
%       Struct where each field name corresponds to a CDF global atttribute.
%       NOTE: Deviates from the usual variable naming conventions.
%       GlobalAttributesSubset field names have
%               the exact names of CDF global attributes.
%
%
% Author: Erik P G Johansson, IRF, Uppsala, Sweden
% First created 2021-03-31, based on two functions broken out of
% execute_sw_mode().
%
function OutGaSubset = derive_output_dataset_GlobalAttributes(...
        InputDatasetsMap, OutputDataset, outputFilename, outputDatasetId, ...
        SETTINGS, L)

    % PROPOSAL: Automatic test code.
    
    % ASSERTIONS
    irf.assert.struct(OutputDataset.Ga, ...
        {'OBS_ID', 'SOOP_TYPE'}, {'Misc_calibration_versions'})

    
    
    OutGaSubset = OutputDataset.Ga;

    OutGaSubset.Parent_version = {};
    OutGaSubset.Parents        = {};
    OutGaSubset.Provider       = {};
    
    %=============================
    % Iterate over INPUT datasets
    %=============================
    keysCa = InputDatasetsMap.keys;
    for i = 1:numel(keysCa)
        
        InputDatasetInfo = InputDatasetsMap(keysCa{i});
        InputGa          = InputDatasetInfo.Ga;

        % ASSERTION
        % NOTE: ROC DFMD is not completely clear on which version number should
        % be used.
        % NOTE: Stores all values to be safe.
        assert(isscalar(InputGa.Data_version), ...
            'BICAS:DatasetFormat', ...
            ['Global attribute "Data_version" for input dataset', ...
            ' with key=%s is not a MATLAB scalar (i.e. the global attribute is', ...
            ' not exactly ONE string).'], ...
            keysCa{i})

        % 2020-12-16, EJ: Has found input datasets to have global
        % attribute "Data_version" values which are either NUMERIC or STRINGS
        % (e.g. "02"). Varies.
        %
        %-----------------------------------------------------------------------
        % Ex: solo_L2_rpw-lfr-surv-bp1-cdag_20200625_V10.cdf:
        % NOTE: Seems to take last two digits from basename as Parent_version
        % values, even for RCTs (i.e. a BUG)!!!
        % Parent_version (9 entries):
        %       0 (CDF_CHAR/2):         "09"
        %       1 (CDF_CHAR/2):         "08"
        %       2 (CDF_CHAR/2):         "08"
        %       3 (CDF_CHAR/2):         "20"
        %       4 (CDF_CHAR/2):         "20"
        %       5 (CDF_CHAR/2):         "43"
        %       6 (CDF_CHAR/2):         "00"
        %       7 (CDF_CHAR/2):         "07"
        %       8 (CDF_CHAR/2):         "00"
        % Parents (9 entries):
        %       0 (CDF_CHAR/46):        "solo_L1_rpw-lfr-surv-bp1-cdag_20200625_V09.cdf"
        %       1 (CDF_CHAR/32):        "solo_HK_rpw-lfr_20200625_V08.cdf"
        %       2 (CDF_CHAR/32):        "solo_HK_rpw-bia_20200625_V08.cdf"
        %       3 (CDF_CHAR/40):        "SOLO_CAL_RCT-LFR-SCM_V20190123171020.cdf"
        %       4 (CDF_CHAR/41):        "SOLO_CAL_RCT-LFR-BIAS_V20190123171020.cdf"
        %       5 (CDF_CHAR/40):        "SOLO_CAL_RCT-LFR-VHF_V20200720165743.cdf"
        %       6 (CDF_CHAR/55):        "SOLO_CAL_RCT-SCM_RPW_SCM-FM-MEB-PFM_V20190519120000.cdf"
        %       7 (CDF_CHAR/35):        "SOLO_CAL_RPW_BIAS_V202003101607.cdf"
        %       8 (CDF_CHAR/42):        "SOLO_CAL_RPW-HF-PREAMP_V20200624000000.cdf"
        %-----------------------------------------------------------------------
        % Ex: solo_L1_rpw-lfr-surv-bp2-cdag_20200625_V09.cdf
        % Parent_version (1 entry):
        %   0 (CDF_INT8/1):         8
        % Parents (1 entry):
        %   0 (CDF_CHAR/33):        "CDF>solo_L0_rpw-cdag_20200625_V08"
        %-----------------------------------------------------------------------
        % Ex: solo_L1_rpw-lfr-surv-bp2-cdag_20201225_V01.cdf, at internal
        % reprocessing (2021-01), not part of regular versioning.
        % Parent_version (1 entry):
        %      0 (CDF_CHAR/2):         "05"
        % Parents (1 entry):
        %      0 (CDF_CHAR/33):        "CDF>solo_L0_rpw-cdag_20201225_V05"
        %-----------------------------------------------------------------------
        % NOTE: Skeletons imply that Parent_version should be strings
        % (CDF_CHAR), though that setting should be inherited from some other
        % dataset, which might be a good or bad source.
        % Ex: SOLO_L2_RPW-LFR-SBM1-CWF-E_V11.skt:
        %   "Parent_version"      1:    CDF_CHAR     { " " }
        %-----------------------------------------------------------------------
        
        % NOTE: Using Data_version to set Parent_version.
        %OutGaSubset.Parent_version{end+1} = InputGa.Data_version{1};   % Number, not string. Correct?!
        %OutGaSubset.Parents       {end+1} = ['CDF>', InputGa.Logical_file_id{1}];
        if isfield(InputGa, 'Provider')
            OutGaSubset.Provider = union(OutGaSubset.Provider, InputGa.Provider);
        else
            % IMPLEMENTATION NOTE: MAG datasets have been observed to not have
            % glob.attr. "Provider". VHT datasets have MAG datasets as parents.
            % /2021-05-05
            % Ex: solo_L2_mag-srf-normal_20200701_V02.cdf
            L.logf('warning', ...
                'Input dataset "%s"\ndoes not have CDF global attribute "Provider".\n', ...
                InputDatasetInfo.filePath)
        end
        
        % NOTE: Parsing INPUT dataset filename to set some global attributes.
        [logicalFileId, ~, dataVersionStr, ~] = parse_dataset_filename(...
            irf.fs.get_name(InputDatasetInfo.filePath));
        % Sets string, not number. Correct?
        OutGaSubset.Parent_version{end+1} = dataVersionStr;
        OutGaSubset.Parents       {end+1} = ['CDF>', logicalFileId];

    end



    OutGaSubset.Software_name       = bicas.constants.SWD_METADATA('SWD.identification.name');
    OutGaSubset.Software_version    = bicas.constants.SWD_METADATA('SWD.release.version');
    
    % BUG? Assigns local time, not UTC!!! ROC DFMD does not mention time zone.
    OutGaSubset.Generation_date     = char(datetime("now","Format","uuuu-MM-dd'T'HH:mm:ss"));
    
    % NOTE: Parsing OUTPUT dataset filename to set some global attributes.
    [logicalFileId, logicalSource, dataVersionStr, timeIntervalStr] = parse_dataset_filename(outputFilename);

    OutGaSubset.Logical_file_id     = logicalFileId;
    % Logical_source: Overwrites skeleton value. Can otherwise not handle -cdag.
    OutGaSubset.Logical_source      = logicalSource;    
    OutGaSubset.Data_version        = dataVersionStr;
    OutGaSubset.Datetime            = timeIntervalStr;
    
    %DataObj.GlobalAttributes.SPECTRAL_RANGE_MIN
    %DataObj.GlobalAttributes.SPECTRAL_RANGE_MAX

    %---------------------------------------------------------------------------
    % "Metadata Definition for Solar Orbiter Science Data", SOL-SGS-TN-0009:
    %   "TIME_MIN   The date and time of the beginning of the first acquisition
    %               for the data contained in the file"
    %   "TIME_MAX   The date and time of the end of the last acquisition for the
    %               data contained in the file"
    %   States that TIME_MIN, TIME_MAX should be "Julian day" (not "modified
    %   Julian day", which e.g. OVT uses internally).
    %
    % NOTE: Implementation does not consider the integration time of each
    % sample.
    % NOTE: juliandate() is consistent with Julian date converter at
    % https://www.onlineconversion.com/julian_date.htm
    % NOTE: ZvsSubset.Epoch already asserted to be monotonically increasing.
    %
    % NOTE: Exact format unclear from documentation, autochecks.
    % NOTE: Issue for autochecks on L3:
    %       https://gitlab.obspm.fr/ROC/DataPool/-/issues/16
    % check_cdf_istp.solo_L3_rpw-bia.txt:
    %   """"
    % 	Global attribute TIME_MAX is of type CDF_DOUBLE.
    % 	    Datatypes other than CDF_CHAR may be problematic.
    % 	Global attribute TIME_MIN is of type CDF_DOUBLE.
    % 	    Datatypes other than CDF_CHAR may be problematic.""""
    % NOTE: ROC data reprocessed ~2021-01-25,
    % solo_L1_rpw-bia-current-cdag_20201201-20201231_V01.cdf (version number
    % probably not part of official versioning) uses
    %     TIME_MIN (1 entry):
    %         0 (CDF_CHAR/17):        "2459184.982450046"
    %     TIME_MAX (1 entry):
    %         0 (CDF_CHAR/17):        "2459215.007218565"
    % Note the number of decimals. No exponent. Other files with ten decimals.
    %
    % PROPOSAL: Copy values from the corresponding values from the relevant input dataset.
    %   CON: Does not work for downsampled.
    %   CON: There has historically been problems with copying bad values from
    %        not-up-to-date input datasets.
    %---------------------------------------------------------------------------
    % NOTE: Choosing 10 decimals (instead of 9) so that time resolution is
    % higher than highest LFR sampling frequency (not sure of highest for
    % TDS-LFM).
    TIME_MINMAX_FORMAT = '%.10f';
    gaTimeMinNbr = juliandate(irf.cdf.TT2000_to_datevec(OutputDataset.Zv.Epoch(1  )));
    gaTimeMaxNbr = juliandate(irf.cdf.TT2000_to_datevec(OutputDataset.Zv.Epoch(end)));
    OutGaSubset.TIME_MIN = sprintf(TIME_MINMAX_FORMAT, gaTimeMinNbr);
    OutGaSubset.TIME_MAX = sprintf(TIME_MINMAX_FORMAT, gaTimeMaxNbr);
    


    enableMods = SETTINGS.get_fv('OUTPUT_CDF.GA_MODS_ENABLED');
    
    MODS = bicas.constants.GA_MODS_DB.get_MODS_strings_CA(outputDatasetId);
    if ~isempty(MODS) && enableMods
        OutGaSubset.MODS = MODS;
    end
    
    
    
    % ROC DFMD hints that value should not be set dynamically. (See meaning of
    % non-italic black text for global attribute name in table.)
    %DataObj.GlobalAttribute.CAVEATS = ?!!

    
    
    % ~ASSERTION
    if ~isscalar(OutGaSubset.Parents)
        [settingValue, settingKey] = SETTINGS.get_fv(...
            'INPUT_CDF.GA_PROVIDER_MISMATCH_POLICY');
        bicas.default_anomaly_handling(...
            L, settingValue, settingKey, 'E+W+illegal', ...
            ['The value of the input CDF files'' global attribute "Provider"', ...
            ' differ (and they should not, or?).'], ...
            'BICAS:DatasetFormat')
        % NOTE: Maybe wrong choice of error ID "DatasetFormat".
    end

    % ASSERTION: Required subset for every dataset.
    irf.assert.struct(OutGaSubset, ...
        {'Parents', 'Parent_version', 'Provider', ...
        'Datetime', 'OBS_ID', 'SOOP_TYPE'}, 'all')
end



% NOTE: Only works correctly for files that follow the official filenaming
% scheme. logicalFileId does not work for e.g. IRFU-internal filenaming
% extension, as e.g. for test files that might be sent to ROC as part of
% official RCS test package.
%
% NOTE: Does not change case.
%
function [logicalFileId, logicalSource, dataVersionStr, timeIntervalStr] ...
        = parse_dataset_filename(filename)
    [~, basename, ~] = fileparts(filename);

    % NOTE: Will include IRFU-internal filenaming extension.
    logicalFileId = basename;

    R = solo.adm.parse_dataset_filename(filename);
    assert(~isempty(R), 'BICAS:Assertion', ...
        ['Can not parse dataset filename "%s" and therefore not', ...
        ' derive values for global attributes', ...
        ' (Logical_source, Data_version, Datetime). Filename does not appear', ...
        ' to follow filenaming conventions.'], filename)

    logicalSource   = R.fnDatasetIdCdag;
    dataVersionStr  = R.versionStr;
    timeIntervalStr = R.timeIntervalStr;
end
