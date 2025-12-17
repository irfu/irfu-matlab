%
% NOTE: THIS CODE IS NEVER CALLED BY BICAS PROPER. It does however use shared
% BICAS functionality for constructing datasets.
%
%
% Code for generating *one* L3 VHT dataset from .mat file, ONLY. The .mat file
% with data is produced by Konrad Steinvall & Yuri Khotyaintsev (2021-03-31).
%
%
% IMPLEMENTATION NOTE
% ===================
% The production of VHT datasets does not fit into BICAS's and ROC's model for
% production of datasets and can therefore not be performed by BICAS proper as
% of 2021-03-31.
% Reasons:
% (1) The RCS interface (as of 2025-01-21) can not handle this case even in
%     principle:
%     Multiple input datasets (one for every day of the month) for one output
%     dataset (one per month). One can not construct such s/w mode.
%     Footnote: BICAS batch processing utilities (bicas.tools.batch) also does
%     not cover this case. Can not find combinations of datasets that determine
%     what to process.
% (2) Output datasets can not be generated from arbitrary (versions of) input
%     datasets. Exact input dataset versions are (should be) the ones used to
%     create the content of .mat file.
%
% NOTE: Future versions may try to read the original (input) datasets that were
% used to produce the .mat file, to complement the output datasets with metadata
% and quality variables.
% VHT files therefore have no ZVs QUALITY_FLAG, QUALITY_BITMASK, or
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
% PROPOSAL: Move to bicas.tools?
%   PRO: Is standalone code, separate from BICAS proper.
%   PRO: Easier to find code.
%   CON: Code is not a "tool", but rather the proper way of generating
%        official datasets.
%
% PROPOSAL: Write code so that it can be transplanted/moved to BICAS proper.
%   CON: Can not be done since it requires multiple input datasets of same
%        DSID.
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
% PROPOSAL: Convert emptyDatasetPolicy string constants to SCREAMING_SNAKE_CASE.
%
% PROBLEM: Can not find .mat file to use!!! /2025-01-21
%   NOTE: VHT datasets (June-Dec 2020, V01) have Epoch values for
%     the time interval:
%       2020-06-01T02:10:00.000000000 to
%       2020-12-31T11:40:00.000000000
%     The files do contain data (not fill values) for the beginning and end
%     timestamps (zVariable "VX_SRF").
%   NOTE: Code assertions on the .mat file:
%     The most common timestamp interval is 10 min.
%     File contains variable "V_RPW".
%   --
%   brain:/data/solo/data_yuri/V_RPW.mat
%     Contains V_RPW but most common time interval is 6 h.
%       mode(diff(b.V_RPW.time.ttns)) / (1e9*60*60*6) == 1
%     Covers the wrong time interval:
%       2020-06-01T02:10:00.000000000Z
%       2020-06-01T08:10:00.000000000Z
%       ... skipped 306 records ...
%       2020-12-05T12:00:00.000000000Z
%       2020-12-05T18:00:00.000000000Z
%   brain:/data/solo/data_yuri/V_RPW_1h.mat
%     Does not contain variable "V_RPW" but does contain "V_RPW_1h".
%     V_RPW_1h has the most common time interval 1h.
%       mode(diff(a.V_RPW_1h.time.ttns)) / (1e9*60*60) == 1
%     Covers the wrong time interval:
%       2020-06-01T02:10:00.000000000Z
%       2020-06-01T03:10:00.000000000Z
%       ... skipped 942 records ...
%       2020-12-05T21:00:00.000000000Z
%       2020-12-05T22:00:00.000000000Z
%   brain:/data/solo/data_yuri/Vp_SWA_PAS.mat
%     Does not contain variable "V_RPW" but does contain "V_pas".
%     V_pas has the most common time interval 4 s.
%       mode(diff(a.V_pas.time.ttns)) / (1e9*4) == 1
%     Covers the wrong time interval:
%       2020-04-14T12:35:09.665000000Z
%       2020-04-14T12:35:13.665000000Z
%       ... skipped 1172549 records ...
%       2020-08-30T23:53:01.575000000Z
%       2020-08-30T23:53:05.575000000Z
%
% TEST CALL:
% bicas.vht.generate_VHT_dataset('/home/erjo/temp/L3/V_RPW.mat', '/nonhome_data/work_files/SOLAR_ORBITER/DataPool/SOLO/RPW/CDF/Master', [2020,07], '/home/erjo/temp/L3', 2, 'ignore empty') -- OBSOLETE /2025-01-21
% bicas.vht.generate_VHT_dataset('/nonhome_data/work_files/SOLAR_ORBITER/DataPool/SOLO/RPW/CDF/Master/SOLO_L3_RPW-BIA-VHT_V03.cdf', [2020,07], "/data/solo/data_yuri/V_RPW.mat", containers.Map, '/home/erjo/temp/vht/solo_L3_rpw-bia-vht-cdag_20200701-20200731_V01.cdf', 'ignore empty') - NEW AND INCOMPLETE /2025-01-21

%     DSID                  = 'SOLO_L3_RPW-BIA-VHT';
%     MASTER_CDF_VERSION_STR      = '01';
EXPECTED_SAMPLE_INTERVAL_NS = int64(10*60*1e9);    % 10 minutes in ns. For assertion on correct .mat file.
DELTA_PLUS_MINUS_NS         = int64(1800*1e9);     % 30 minutes in ns. Contradicts EXPECTED_SAMPLE_INTERVAL_NS?!!

% Used for assertion on data.
% NOTE: Velocity is negative due to coordinate system.
VX_SRF_MIN_KMPS = -1500;
VX_SRF_MAX_KMPS =  0;



% ASSERTIONS
assert(isstring(matFilePath), 'matFilePath is not a string.')
assert((length(yearMonth) == 2) && isnumeric(yearMonth))
assert(isa(InputDatasetsMap, 'containers.Map'))



Bso     = bicas.create_default_BSO();
Bso.make_read_only();
BICAS_L = bicas.Logger('HUMAN_READABLE', false);



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
timeIntervalStr = sprintf(...
  '%04i-%02i-01T00:00:00/%04i-%02i-01T00:00:00', ...
  dv1(1:2), dv2(1:2));
V_RPW = V_RPW.tlim(irf.tint(timeIntervalStr));



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

%---------------------------------------------------------------------------
% IMPORTANT NOTE: BICAS uses
% execute_SWM:get_output_dataset_GAs() to derive many
% global attributes.
%   NOTE: OutGaSubset = get_output_dataset_GAs(...
%       InputDatasetsMap, OutputDataset, outputFilename, Bso, L)
%   Ex: Generation_date, Parents, Software_name (BICAS), Datetime (time
%   interval string from filename)
%---------------------------------------------------------------------------
GaSubset = bicas.ga.get_output_dataset_GAs(...
  InputDatasetsMap, OutputDataset, ...
  irf.fs.get_name(outputFile), Bso, BICAS_L);

bicas.write_dataset_CDF(...
  Zv, GaSubset, outputFile, masterCdfPath, ...
  Bso, BICAS_L)

end
