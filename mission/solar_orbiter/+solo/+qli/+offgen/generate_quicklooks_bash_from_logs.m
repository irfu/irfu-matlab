%
% Wrapper around solo.qli.offgen.generate_quicklooks() for generating quicklooks
% for dates which are automatically derived from implicitly specified log files,
% which should indicate updated datasets.
%
%
% NOTE: This script is NOT intended to be called from MATLAB by the average
%       user. See irfu-matlab/mission/solar_orbiter/+solo/+qli/README.TXT.
%
% NOTE: This function is designed to potentially be called from bash/the OS.
%
% NOTE: The dataset IDs referred to in the implementation must be consistent
%       with the use of datasets in solo.qli.generate_quicklooks_batch().
%
%
% ARGUMENTS
% =========
% outputDir
%       Path to output directory.
% varargin
%       String IDs representing different dataset source directories. The
%       function will use the logs for the specified directories to derive dates
%       for which QLI should be generated.
%
%
% RETURN VALUES
% =============
% (None)
%
%
% Author: Erik P G Johansson, IRF, Uppsala, Sweden
%
function generate_quicklooks_bash_from_logs(outputDir, varargin)

LESIA_LOG_FILE_DIR_PATTERN = '/home/erjo/logs/pull.so.data.cron.brain.*.log';
SOAR_LOG_FILE_DIR_PATTERN  = '/home/erjo/logs/so_soar_irfu_mirror_sync.*.log';



AllDaysDtArray = datetime.empty(0, 1);
AllDaysDtArray.TimeZone = 'UTCLeapSeconds';

for i = 1:numel(varargin)
  datasetsSource = varargin{i};
  assert(ischar(datasetsSource), 'datasetsSource is not a string.')

  switch(datasetsSource)

    %     case 'IRFU'
    %       DaysDtArray = solo.qli.extract_dataset_dates_from_logs(...
    %         '/home/erjo/logs/so_bicas_batch_cron.*.log', ...
    %         {...
    %           'solo_L3_rpw-bia-efield-10-seconds', ...
    %           'solo_L3_rpw-bia-density-10-seconds'...
    %         });

    case 'LESIA'
      %       DaysDtArray  = solo.qli.extract_dataset_dates_from_logs(...
      %         '/home/erjo/logs/pull.so.data.cron.brain.*.log', ...
      %         {'solo_L2_rpw-tnr-surv-cdag'});
      DaysDtArray  = solo.qli.extract_dataset_dates_from_logs(...
        LESIA_LOG_FILE_DIR_PATTERN, ...
        {
        'solo_L2_rpw-tnr-surv', ...
        'solo_L3_rpw-bia-efield-10-seconds', ...
        'solo_L3_rpw-bia-density-10-seconds'...
        });

    case 'SOAR'
      DaysDtArray  = solo.qli.extract_dataset_dates_from_logs(...
        SOAR_LOG_FILE_DIR_PATTERN, ...
        {
        'solo_L2_mag-rtn-normal', ...
        'solo_L2_mag-rtn-normal-1-minute', ...
        'solo_L2_swa-pas-eflux', ...
        'solo_L2_swa-pas-grnd-mom', ...
        });

    otherwise
      error('Illegal datasetsSource="%s" (argument %i)', datasetsSource, i)
  end

  AllDaysDtArray = [AllDaysDtArray; DaysDtArray];
end

AllDaysDtArray = unique(AllDaysDtArray);

solo.qli.offgen.generate_quicklooks(outputDir, true, true, AllDaysDtArray)

end
