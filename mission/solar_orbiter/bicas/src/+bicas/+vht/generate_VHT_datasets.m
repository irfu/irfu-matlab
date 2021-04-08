%
% NOTE: THIS CODE IS NEVER CALLED BY BICAS PROPER. It does however use shared
% BICAS functionality for constructing datasets.
%
%
% Code for generating L3 VHT datasets from .mat file, ONLY. .mat file with data
% is produced by Konrad Steinvall & Yuri Khotyaintsev (2021-03-31).
%
% IMPLEMENTATION NOTE
% ===================
% The production of VHT datasets does not fit into BICAS's and ROC's model for
% production of datasets and can therefore not be performed by BICAS as of
% 2021-03-31.
% Reasons:
% (1) The RCS interface (2021-03-31) can not handle this case in principle:
%     Multiple input datasets (one for every day of the month) for one output
%     dataset (one per month). Can not construct such s/w mode.
%     Footnote: Erik Johansson's BICAS batch processing utilities (not part of
%     BICAS) also do not cover this case. Can not find combinations of datasets
%     that determine what to process.
% (2) output datasets can not be generated from arbitrary (versions of) input
%     datasets. Exact input dataset versions are linked to the content of
%     .mat file.
%
% NOTE: Future versions may try to read the original (input) datasets that were
% used to produce the .mat file, to complement the output datasets with metadata
% and quality variables.
%
%
% ARGUMENTS
% =========
% emptyDatasetPolicy
%       How to handle months without data.
%
%
% RETURN VALUES
% =============
% (none)
%
%
% Author: Erik P G Johansson, IRF Uppsala, Sweden
% First created 2021-03-26.
%
function generate_VHT_datasets(...
        matFilePath, masterCdfDir, yearMonth, outputDir, outputVerNbr, ...
        emptyDatasetPolicy)
%
% TODO-DEC: Behaviour fro time interval without data?
%   PROPOSAL: Error
%   PROPOSAL: Produce empty dataset
%   PROPOSAL: Do nothing.
%
% PROPOSAL: Write code so that it can be transplanted/moved to BICAS proper.
%   CON: Can not be done if requires multiple input dataset of same
%        DATASET_ID.
%
% TODO-DEC: How specify time interval?
%   PROPOSAL: Year-month
%   PROPOSAL: Timestamp object. Any timestamp in month interval suffices.
%       PROPOSAL: EpochTT
%       PROPOSAL: datetime
%       CON: Easy to extract year+month from timestamp formats.
%
% PROPOSAL: Separate into multiple scripts.
%   PROPOSAL: Collect in separate package.
%       PROPOSAL: +vht/
%       TODO-DEC: Location?
%           PROPOSAL: root-level.
%               CON: Collisions.
%           PROPOSAL: bicas.vht
%
%   PROPOSAL: Batch script: Select output dataset versions automatically.
%   PROPOSAL: Process one file in-->one file out.
%       NOTE: Will read entire .mat file every time. ~Inefficient.
%
% TODO-DEC: How set QUALITY_FLAG, QUALITY_BITMASK, L2_QUALITY_BITMASK ?
%   TODO-NI: Required? Does YK want them?
%   PROPOSAL: Set L2_QUALITY_BITMASK from NSOPS.
%   PROPOSAL: Set using relevant L2 input file behind data?
%       NOTE: Input L2 files use different time resolution.
%
% TODO: Derive parents.
%
% TEST CALL:
% bicas.vht.generate_VHT_datasets('/home/erjo/temp/L3/V_RPW.mat', '/nonhome_data/work_files/SOLAR_ORBITER/DataPool/SOLO/RPW/CDF/Master', [2020,07], '/home/erjo/temp/L3', 2, 'ignore empty')
    
    DATASET_ID             = 'SOLO_L3_RPW-BIA-VHT';
    MASTER_CDF_VERSION_STR = '01';
    DELTA_PLUS_MINUS_NS    = int64(3*3600*1e9);
    
    % ASSERTIONS
    assert(ischar(matFilePath))
    assert((length(yearMonth) == 2) && isnumeric(yearMonth))
    EJ_library.assert.dir_exists(masterCdfDir)
    
    
    
    BICAS_SETTINGS = bicas.create_default_SETTINGS();
    BICAS_SETTINGS.make_read_only();
    BICAS_L = bicas.logger('human-readable', false);
    
    
    
    %================
    % READ .mat FILE
    %================
    load(matFilePath, 'V_RPW');
    % ASSERTION: .mat file
    mostCommonTimeDiffSec = mode(diff(V_RPW.time.ttns));
    assert(...
        mostCommonTimeDiffSec == 2*DELTA_PLUS_MINUS_NS, ...
        ['Timestamps in %s (mostCommonTimeDiffSec=%i) do not seem consistent', ...
        ' with the expected DELTA_PLUS_MINUS = %i'], ...
        mostCommonTimeDiffSec, DELTA_PLUS_MINUS_NS)

    %==============================================
    % Only keep data for the specified time period
    %==============================================
    % Beginning & end of month.
    dv1 = datevec(datetime([yearMonth(1), yearMonth(2),   1]));
    dv2 = datevec(datetime([yearMonth(1), yearMonth(2)+1, 1]));
    % IMPLEMENTATION NOTE: Slight hack using intermediate UTC string, but there
    % is no (?) smooth way of converting date/time vector-->EpochTT.
    timeIntStr = sprintf(...
        '%04i-%02i-01T00:00:00/%04i-%02i-01T00:00:00', ...
        dv1(1:2), dv2(1:2));
    V_RPW = V_RPW.tlim(irf.tint(timeIntStr));
    
    
    
    %========================================
    % Handle datasets/months with empty data
    %========================================
    if isempty(V_RPW)
        switch(emptyDatasetPolicy)
            case 'assert non-empty'
                error(['Trying to create empty dataset.', ...
                    ' There is no data for yearMonth=[%d, %d]'], ...
                    yearMonth(:))
                
            case 'ignore empty'
                fprintf('There is no data for yearMonth=[%d, %d]. Ignoring.\n', yearMonth(:))
                return
                
            otherwise
                error('Illegal argument="%s".', emptyDatasetPolicy)
        end
    end
    
    
    
    %==========================
    % CONSTRUCT DATASET STRUCT
    %==========================
    Zv = [];
    Zv.Epoch            = V_RPW.time.ttns;
    Zv.Vx_SRF           = V_RPW.data;
    Zv.DELTA_PLUS_MINUS = DELTA_PLUS_MINUS_NS + Zv.Epoch*0;

    Ga = [];
    Ga.OBS_ID    = ' ';
    Ga.SOOP_TYPE = ' ';
    
    OutputDataset = [];
    OutputDataset.Zv = Zv;
    OutputDataset.Ga = Ga;



    %=====================
    % Create dataset file
    %=====================
    %GaSubset = struct();
    
    outputFile = fullfile(...
        outputDir, ...
        bicas.vht.derive_dataset_filename(yearMonth, outputVerNbr));
    masterCdfFileName = bicas.get_master_CDF_filename(...
                DATASET_ID, ...
                MASTER_CDF_VERSION_STR);
    masterCdfPath = fullfile(masterCdfDir, masterCdfFileName);
    
    InputDatasetsMap = containers.Map();    % NO PARENTS! -- TEMP    
    %---------------------------------------------------------------------------
    % IMPORTANT NOTE: BICAS uses
    % execute_sw_mode:derive_output_dataset_GlobalAttributes() to derive many
    % global attributes.
    %   NOTE: OutGaSubset = derive_output_dataset_GlobalAttributes(...
    %       InputDatasetsMap, OutputDataset, outputFilename, SETTINGS, L)
    %   Ex: Generation_date, Parents, Software_name (BICAS), Datetime (time
    %   interval string from filename)    
    %---------------------------------------------------------------------------
    GaSubset = bicas.derive_output_dataset_GlobalAttributes(...
        InputDatasetsMap, OutputDataset, ...
        EJ_library.fs.get_name(outputFile), BICAS_SETTINGS, BICAS_L);
    
    bicas.write_dataset_CDF(...
        Zv, GaSubset, outputFile, masterCdfPath, BICAS_SETTINGS, BICAS_L)

end