%
% Demonstration code for how to generate "IRFU quicklooks" ("QLI") locally.
%
%
% NOTE: The configuration is system-dependent! The function is configured for
% being run on brain/spis at IRFU, or any system where /data/solo/ has been
% mounted to the same location.
%
% NOTE: GENERATING QUICKLOOKS IS TIME-CONSUMING, in particular when MAG data is
% available (80% of the time). Depending on the machine, it may take on the
% order of 20 minutes per day of data with MAG data.
%
% NOTE: The exact locations and sizes of text and panels on quicklooks generated
% locally may be different or outright wrong (for text) when compared to
% officially generated quicklooks.
%
%
% ARGUMENTS
% =========
% outputDir
%       Path to output directory where quicklooks will be generated.
%
%
% Author: Erik P G Johansson, IRF, Uppsala, Sweden
%
function generate_quicklooks_demo(outputDir)

VHT_DIR                       = '/data/solo/data_yuri/';
GENERATE_NONWEEKLY_QUICKLOOKS = true;   % Whether to generate 24h, 6h, and 2h quicklooks.
GENERATE_WEEKLY_QUICKLOOKS    = true;

% datetime column array. Array of UTC midnights representing the beginning of
% days for which to generate quicklooks. Specifies explicit individual days (not
% a time range).
DAYS_DATETIME_COLUMN_ARRAY    = datetime({...
  '2022-11-30T00:00:00Z'; ...
  '2022-12-03T00:00:00Z'}, 'TimeZone', 'UTCLeapSeconds');



%===============================================================
% Configure Solar Orbiter database from which data will be used
%===============================================================
% NOTE: Only uses /data/solo/remote/data/ and /data/solo/soar/ (and not
%       /data/solo/data_irfu/) since data_irfu/ has less reliable data
%       ("bleeding edge").
%
% IMPORTANT NOTE
% --------------
% The quicklook generation relies on the path specified submitted to
% solo.db_init(), but it also relies on a hardcoded path to
% /data/solo/remote/data/L2/thr/ in solo.read_TNR(). There is currently no way
% of changing this without changing that hardcoded path and therefore currently
% no way of generating quicklooks without this path available. This restriction
% may be eliminated eventually.

irf()    % Needed to make "SolO DB" work.
solo.db_init('local_file_db', '/data/solo/');
% Setup cache
solo.db_init('db_cache_size_max', 4096)   % Unit: MiB.
solo.db_cache('on', 'save')



%=====================
% Generate quicklooks
%=====================
solo.qli.generate_quicklooks_all_types(...
  [], VHT_DIR, outputDir, ...
  GENERATE_NONWEEKLY_QUICKLOOKS, GENERATE_WEEKLY_QUICKLOOKS, ...
  DAYS_DATETIME_COLUMN_ARRAY)
end