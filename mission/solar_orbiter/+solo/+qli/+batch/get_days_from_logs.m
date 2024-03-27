%
% Generate array of dates derived from implicitly specified log files, which
% should indicate updated datasets.
%
%
% NOTE: Several arguments are designed to partly handle arguments deriving from
%       the bash/the OS and are therefore on a string format.
%
% NOTE: The dataset IDs referred to in the implementation must be consistent
%       with the use of datasets in
%       solo.qli.generate_quicklooks_24h_6h_2h_using_DB_SPICE() and
%       solo.qli.generate_quicklook_7days().
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
% DaysDtArray
%
%
% Author: Erik P G Johansson, IRF, Uppsala, Sweden
%
function DaysDtArray = get_days_from_logs(Settings, varargin)
% PROPOSAL: Move dataset IDs to solo.qli.batch.const (new)

DaysDtArray = solo.qli.const.EMPTY_DT_ARRAY;

for i = 1:numel(varargin)
  datasetsSourceId = varargin{i};
  assert(ischar(datasetsSourceId), 'datasetsSource is not a string.')

  switch(datasetsSourceId)

    %     case 'IRFU'
    %       DaysDtArray = solo.qli.batch.extract_dataset_dates_from_logs(...
    %         '/home/erjo/logs/so_bicas_batch_cron.*.log', ...
    %         {...
    %           'solo_L3_rpw-bia-efield-10-seconds', ...
    %           'solo_L3_rpw-bia-density-10-seconds'...
    %         });

    case 'LESIA'
      SourceDaysDtArray = solo.qli.batch.extract_dataset_dates_from_logs(...
        Settings.lesiaLogFileDirPattern, ...
        {
        'solo_L2_rpw-tnr-surv', ...
        'solo_L3_rpw-bia-efield-10-seconds', ...
        'solo_L3_rpw-bia-density-10-seconds'...
        });

    case 'SOAR'
      SourceDaysDtArray = solo.qli.batch.extract_dataset_dates_from_logs(...
        Settings.soarLogFileDirPattern, ...
        {
        'solo_L2_mag-rtn-normal', ...
        'solo_L2_mag-rtn-normal-1-minute', ...
        'solo_L2_swa-pas-eflux', ...
        'solo_L2_swa-pas-grnd-mom', ...
        });

    otherwise
      error('Illegal datasetsSourceId="%s"', datasetsSourceId)
  end

  DaysDtArray = [DaysDtArray; SourceDaysDtArray];
end

DaysDtArray = unique(DaysDtArray);

end
