%
% NOTE: THIS CODE IS NEVER CALLED BY BICAS PROPER. It does however use shared
% BICAS functionality for constructing datasets.
%
%
% Code for generating one L3 VHT dataset from .mat file, ONLY. The .mat file
% with data is produced by Konrad Steinvall & Yuri Khotyaintsev (2021-03-31).
%
%
% IMPLEMENTATION NOTE
% ===================
% The production of VHT datasets does not fit into BICAS's and ROC's model for
% production of datasets and can therefore not be performed by BICAS proper as
% of 2021-03-31.
% Reasons:
% (1) The RCS interface (2021-03-31) can not handle this case in principle:
%     Multiple input datasets (one for every day of the month) for one output
%     dataset (one per month). Can not construct such s/w mode.
%     Footnote: Erik P G Johansson's BICAS batch processing utilities (not part of
%     BICAS) also do not cover this case. Can not find combinations of datasets
%     that determine what to process.
% (2) Output datasets can not be generated from arbitrary (versions of) input
%     datasets. Exact input dataset versions are linked to the content of
%     .mat file.
%
% NOTE: Future versions may try to read the original (input) datasets that were
% used to produce the .mat file, to complement the output datasets with metadata
% and quality variables.
% VHT files therefore has no ZVs QUALITY_FLAG, QUALITY_BITMASK,
% L2_QUALITY_BITMASK.
%
%
% ARGUMENTS
% =========
% yearMonth
%       1D vector: [yearNbr, monthNbr]
% emptyDatasetPolicy
%       How to handle months without data.
%
%
% RETURN VALUES
% =============
% (none)
%
%
% Author: Erik P G Johansson, IRF, Uppsala, Sweden
% First created 2021-03-26.
%
function generate_VHT_dataset(...
        masterCdfPath, yearMonth, matFilePath, InputDatasetsMap, outputFile, ...
        emptyDatasetPolicy)
%
% PROPOSAL: Write code so that it can be transplanted/moved to BICAS proper.
%   CON: Can not be done since it requires multiple input datasets of same
%        DSI.
%
% TODO-DEC: How specify month?
%   PROPOSAL: Year+month array, explicitly
%   PROPOSAL: Some timestamp object. Any timestamp in month interval suffices.
%       PROPOSAL: EpochTT
%       PROPOSAL: datetime
%       CON: Easy to extract year+month from timestamp formats.
%
% TODO-DEC: How set QUALITY_FLAG, QUALITY_BITMASK, L2_QUALITY_BITMASK ?
%   TODO-NI: Required? Does YK want them?
%   PROPOSAL: Set L2_QUALITY_BITMASK from NSOPS.
%   PROPOSAL: Set using relevant L2 input file behind data?
%       NOTE: Input L2 files use different time resolution.
%
% TEST CALL:
% bicas.vht.generate_VHT_dataset('/home/erjo/temp/L3/V_RPW.mat', '/nonhome_data/work_files/SOLAR_ORBITER/DataPool/SOLO/RPW/CDF/Master', [2020,07], '/home/erjo/temp/L3', 2, 'ignore empty')

%     DSI                  = 'SOLO_L3_RPW-BIA-VHT';
%     MASTER_CDF_VERSION_STR      = '01';
    EXPECTED_SAMPLE_INTERVAL_NS = int64(10*60*1e9);    % For assertion.
    DELTA_PLUS_MINUS_NS         = int64(1800*1e9);

    % Used for assertion on data.
    % NOTE: Velocity is negative due to coordinate system.
    VX_SRF_MIN_KMPS = -1500;
    VX_SRF_MAX_KMPS =  0;



    % ASSERTIONS
    assert(ischar(matFilePath))
    assert((length(yearMonth) == 2) && isnumeric(yearMonth))
    assert(isa(InputDatasetsMap, 'containers.Map'))



    Bso     = bicas.create_default_BSO();
    Bso.make_read_only();
    BICAS_L = bicas.Logger('human-readable', false);



    %================
    % READ .mat FILE
    %================
    load(matFilePath, 'V_RPW');
    % ASSERTION: .mat file
    mostCommonTimeDiffSec = mode(diff(V_RPW.time.ttns));
    assert(...
        mostCommonTimeDiffSec == EXPECTED_SAMPLE_INTERVAL_NS, ...
        ['Timestamps in %s (mostCommonTimeDiffSec=%i) do not seem consistent', ...
        ' with the expected time intervals between samples,', ...
        ' EXPECTED_SAMPLE_INTERVAL_NS = %i'], ...
        matFilePath, mostCommonTimeDiffSec, EXPECTED_SAMPLE_INTERVAL_NS)
    % NOTE: bicas.write_dataset_CDF() should replace NaN-->Fill value, but given
    % that VHT .mat file contains data gaps as absence of timestamps, it would
    % be surprising if it contained NaN.
    assert(all(~isnan(V_RPW.data)), 'Found NaN in V_RPW.data.')
    assert(all((VX_SRF_MIN_KMPS <= V_RPW.data) & (V_RPW.data <= VX_SRF_MAX_KMPS)))


    %==============================================
    % Only keep data for the specified time period
    %==============================================
    % Beginning & end of calendar month.
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
                fprintf(...
                    'There is no data for yearMonth=[%d, %d]. Ignoring.\n', ...
                    yearMonth(:))
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
    Zv.VX_SRF           = V_RPW.data;
    Zv.DELTA_PLUS_MINUS = DELTA_PLUS_MINUS_NS + Zv.Epoch*0;

    Ga = [];
    Ga.OBS_ID    = ' ';
    Ga.SOOP_TYPE = ' ';

    OutputDataset    = [];
    OutputDataset.Zv = Zv;
    OutputDataset.Ga = Ga;



    %=====================
    % Create dataset file
    %=====================

    %InputDatasetsMap = containers.Map();    % NO PARENT DATASETS! -- TEMP

    %---------------------------------------------------------------------------
    % IMPORTANT NOTE: BICAS uses
    % execute_SWM:derive_output_dataset_GAs() to derive many
    % global attributes.
    %   NOTE: OutGaSubset = derive_output_dataset_GAs(...
    %       InputDatasetsMap, OutputDataset, outputFilename, Bso, L)
    %   Ex: Generation_date, Parents, Software_name (BICAS), Datetime (time
    %   interval string from filename)
    %---------------------------------------------------------------------------
    GaSubset = bicas.derive_output_dataset_GAs(...
        InputDatasetsMap, OutputDataset, ...
        irf.fs.get_name(outputFile), Bso, BICAS_L);

    bicas.write_dataset_CDF(...
        Zv, GaSubset, outputFile, masterCdfPath, ...
        Bso, BICAS_L)

end
