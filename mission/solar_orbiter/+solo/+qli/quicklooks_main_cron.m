%
% Wrapper around solo.qli.quicklooks_main() intended for being used by being
% called from system scripts (e.g. bash) for the purpose of cron jobs on
% brain/spis. The arguments have also been designed for this purpose and
% therefore all strings.
%
% NOTE: This script is NOT intended to be called from MATLAB by the average
%       user. See solo.qli.quicklooks_main() instead.
%
%
% ARGUMENTS
% =========
% logoPath
%       Path to IRF logo image.
%       Normally located in irfu-matlab:
%       irfu-matlab/mission/solar_orbiter/+solo/irf_logo.png
%       Empty ==> Do not plot any logo.
% vhtDataDir
%       Path to directory containing VHT (velocity) .mat files.
% outputDir
%       Plots will be placed in subdirectories under this directory.
%       NOTE: Will create subdirectories if not pre-existing.
% runNonweeklyPlots, runWeeklyPlots
%       NOTE: STRINGS.
%       Whether to run ("1") or not run ("0") the resp. groups of plots.
% utcBegin, utcEnd : Strings.
%       Defines time interval for which quicklooks should be generated.
%
%
%
% Author: Erik P G Johansson, IRF, Uppsala, Sweden
% First created 2022-08-30.
%
function quicklooks_main_cron(...
        logoPath, vhtDataDir, outputDir, ...
        runNonweeklyPlots, runWeeklyPlots, utcBegin, utcEnd)

    runNonweeklyPlots = interpret_argument_flag(runNonweeklyPlots);
    runWeeklyPlots    = interpret_argument_flag(runWeeklyPlots);

    % IMPLEMENTATION NOTE: Needed to make "DB" work. Necessary when calling from
    % bash.
    irf

    %=================================
    % Configure Solar Orbiter database
    %=================================
    % NOTE: System-dependent configuration!
    solo.db_init('local_file_db', '/data/solo/');
    solo.db_init('local_file_db', '/data/solo/data_irfu');
    % Setup cache
    solo.db_init('db_cache_size_max', 4096)
    solo.db_cache('on', 'save')

    solo.qli.quicklooks_main(...
        logoPath, vhtDataDir, outputDir, ...
        runNonweeklyPlots, runWeeklyPlots, utcBegin, utcEnd)

end



% Interpret argument for main function interface. Intended accept and normalize
% arguments which are either
% (1) MATLAB-friendly (numeric/logical), or
% (2) bash script-friendly (strings).
%
function value = interpret_argument_flag(arg)
    assert(isscalar(arg), 'Flag argument is not scalar.')

    if ischar(arg) && arg=='0'
        value = false;
    elseif ischar(arg) && arg=='1'
        value = true;
    else
        error('Can not interpret argument flag. Illegal format.')
    end
end
