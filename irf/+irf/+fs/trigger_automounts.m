%
% Try to trigger any automounting of disks (if there is any automount) in order
% to avoid file-reading problems for long batch runs. The function tries to
% trigger automount(s) by simply accessing specified paths without actually
% doing anything else. Since the code has no effect on non-automounts, it can
% potentially be used on all paths which are used by a code.
%
% This function is meant to be called for representative paths (e.g. "root
% directories") before every period of "continous" disk access by code which
% does long unattended runs.
%
%
% CONTEXT
% =======
% The automounting which should make IRFU's /data/*/ directories available on
% brain/spis (IRFU) does not always work as intended. It appears to not always
% work fast enough for reading files on rare occassions, leading to apparent
% can-not-read-file errors for existing files when making long runs, in
% particular during out-of-office hours(?). Other code which has faced this
% problem has seemingly resolved it by using Linux commands to list the
% directory contents and ignore the result before using /data/solo/ (such as
% calling "ls /data/solo/ >> /dev/null" in bash scripts).
%
%
% ARGUMENTS
% =========
% pathsCa
%       Column cell array of paths to files and directories.
%       NOTE: Exception is raised if any path is not valid.
%
%
% Author: Erik P G Johansson, IRF, Uppsala, Sweden
%
function trigger_automounts(pathsCa)
assert(iscell(pathsCa) & iscolumn(pathsCa))

for i = 1:numel(pathsCa)
  path = pathsCa{i};

  % irf.log('n', sprintf(...
  %   'Trying to trigger automounting, if not already mounted: %s', path ...
  %   ))
  junk = dir(path);
end

end
