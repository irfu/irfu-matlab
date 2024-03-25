%
% Demonstration code for how to generate "IRFU quicklooks" ("QLI") locally. The
% user can copy and modify this to generate quicklooks her-/himself.
%
%
% NOTE: The example implementation below is system-dependent! The function is
% configured for being run on brain/spis at IRFU, or any system where
% /data/solo/ has been mounted to the same location.
%
% NOTE: GENERATING QUICKLOOKS IS TIME-CONSUMING, in particular 24h, 6h, 2h
% quicklooks when MAG data is available (80% of the mission time). Depending on
% the machine, it may take on the order of 20-30 minutes per day of data when
% there is MAG data!
%
% NOTE: The exact locations and sizes of text and panels on quicklooks as seen
% in windows (figures) in the OS GUI may be different or outright wrong when
% compared to file versions of quicklooks (for text, e.g. bad location).
%
%
% ARGUMENTS
% =========
% outputDir
%       Path to output directory where quicklooks will be generated.
%       Different types of quicklooks will be generated under separate
%       subdirectories under this directory.
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
% DAYS_DATETIME_COLUMN_ARRAY    = datetime({...
%   '2022-12-01T00:00:00Z'; ...
%   '2022-12-03T00:00:00Z'; ...
%   '2022-12-05T00:00:00Z'}, 'TimeZone', 'UTCLeapSeconds');
DAYS_DATETIME_COLUMN_ARRAY    = datetime({...
  '2022-12-05T00:00:00Z'}, 'TimeZone', 'UTCLeapSeconds');



%===============================================================
% Configure Solar Orbiter database from which data will be used
%===============================================================
% NOTE: Only uses /data/solo/remote/data/ and /data/solo/soar/ (and not
%       /data/solo/data_irfu/) since data_irfu/ has less reliable data
%       ("bleeding edge").

irf()    % Needed to make "SolO DB" work.
solo.db_init('local_file_db', '/data/solo/');
% Setup cache
solo.db_init('db_cache_size_max', 4096)   % Unit: MiB.
solo.db_cache('on', 'save')



% Enable more log messages. Not necessary, but is useful for debugging..
irf.log('notice')

%=====================
% Generate quicklooks
%=====================
solo.qli.batch.generate_quicklooks(...
  [], VHT_DIR, outputDir, ...
  GENERATE_NONWEEKLY_QUICKLOOKS, GENERATE_WEEKLY_QUICKLOOKS, ...
  DAYS_DATETIME_COLUMN_ARRAY)
end
