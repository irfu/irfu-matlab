function load_mkernel(mission, varargin)
% LOAD_MKERNEL load a locally adapted metakernel for SPICE computations
%
%	LOAD_MKERNEL(mission, [flown_or_predicted]) load locally adopted SPICE
%   metakernel for "mission", in order to do orbit computations.
%
%   'mission' - mission to load (must be one of the following):
%               'bepicolombo'
%               'juice'
%               'parkersolarprobe'
%               'rosetta'
%               'solarorbiter'
%               'comet-interceptor'
%   ['flown_or_predicted'] - optional indicator of the metakernel to load:
%          'flown'     = only the reconstructed actual flown orbit kernels.
%          'predicted' = predicted orbit kernels files, (for past dates
%                        reconstructed files are loaded as well).
%           If omitted will default to 'predicted'.
%
% First time loading a new mission (on any computer) load_mkernel will ask
% user for local (root) paths to the SPICE kernel repository.
% These paths will be remebered for future Matlab sessions, however it will
% be reset for "mission" if it is no longer present as a path on the system
% during next session it is trying to load metakernels for "mission".
%
%	Example:
%		load_mkernel('solarorbiter', 'predicted');
%   % Then do some calculations (position of 'solo'):
%   Tint = irf.tint('2020-07-31T00:00:00Z/2020-08-01T15:00:00Z');
%   et = Tint.start.tts:3600:Tint.stop.tts;
%   pos = cspice_spkpos('solo', et, 'ECLIPJ2000', 'LT+s', 'Sun');
%
% Internally use a global variable "LoadedSpiceKernel" to keep track of
% what has been loaded already.
%
% Important Note: The order of loaded kernels is important for the result!
% For computation which use body a found in multiple kernels it is up to
% YOU to ensure it is loaded correctly if YOU chose to load multiple
% missions without clearing the loaded kernels in between.
% Each missions own metakernel should be selfconsistent, but no guarantee
% is given for inter-mission loading of SPICE metakernels.
%
% See "SPICE Kernel Required Reading", available from NASA NAIF website.

p = inputParser;
addRequired(p, 'mission', ...
  @(x) any(validatestring(x, ...
  {'bepicolombo', 'juice', 'parkersolarprobe', 'rosetta', 'solarorbiter', 'comet-interceptor'})));
addOptional(p, 'flown_or_predicted', 'predicted', ...
  @(x) any(validatestring(x, {'flown', 'predicted'})));
parse(p, mission, varargin{:});

global LoadedSpiceKernels % keep track of loaded kernels

if ispc
  errStr = 'Wnidows system not yet supported.';
  irf.log('critical', errStr);
  error(errStr);
end

% Get total number of loaded metakernels
nMeta = cspice_ktotal('meta');
if nMeta == 0, LoadedSpiceKernels = []; end
for iMeta = 1:nMeta
  % Get info of each loaded metakernel
  file = cspice_kdata(iMeta, 'meta');
  if exist(file, 'file')
    % Everything is OK, file still exists on our system
  else
    % We have some inconsistency, a loaded MK no longer exists on system.
    % Unload all and start fresh.
    cspice_kclear();
    LoadedSpiceKernels = [];
    break; % Break out of for loop
  end
end

if ~isempty(LoadedSpiceKernels)
  if isfield(LoadedSpiceKernels, p.Results.mission) && ~isfield(LoadedSpiceKernels.(p.Results.mission), p.Results.flown_or_predicted)
    % We have been requested a different flown_or_predicted orbit for
    % already loaded mission! This is not OKEY!
    irf.log('warning', 'We had already loaded a kernel for mission with different flown_or_predicted status. NOT OKEY! Clearing kernels.');
    % Unload all and start fresh.
    cspice_kclear();
    LoadedSpiceKernels = [];
  elseif isfield(LoadedSpiceKernels, p.Results.mission) && ~exist(LoadedSpiceKernels.(p.Results.mission).(p.Results.flown_or_predicted), 'file')
    % We have been requested a previously loaded mission and flown status,
    % but its local file no longer exist.
    irf.log('warning', 'We had previously loaded a kernel for mission with this flown_or_predicted status, but that file is no longer present. NOT OKEY! Clearing kernels.');
    % Unload all and start fresh.
    cspice_kclear();
    LoadedSpiceKernels = [];
  end
end

if isempty(LoadedSpiceKernels) || ~isfield(LoadedSpiceKernels, p.Results.mission)
  % New loading (entirly new, or new mission)
  rootKernelPath = get_spice_root(p.Results.mission);

  mkPath = [rootKernelPath, filesep, 'kernels', filesep, 'mk'];
  if ~exist(mkPath, 'dir')
    errStr = ['Did not find the expected metakernel path: ', mkPath];
    irf.log('critical', errStr);
    datastore('spice_paths', p.Results.mission, []); % Clear out stored path
    error(errStr);
  end

  % Default names for "flown" and "predicted" kernels, (each mission have
  % their own naming standard).
  flown.rosetta          = 'ROS_OPS_*.TM'; % Rosetta flown mk name standard
  flown.juice            = 'juice_crema_5_0b23_1.tm'; % FIXME: Some other orbit senario?/Update when JUICE has launched.
  flown.bepicolombo      = 'bc_ops_*.tm'; % BepiColombo flown
  flown.parkersolarprobe = 'test*.tm'; % FIXME: UPDATE WHEN PSP sync script is tested
  flown.solarorbiter     = 'solo_ANC_soc-flown-mk_v*.tm'; % SolO flown
  flown.cometinterceptor = 'interceptor_study_v02.tm'; % FIXME: Some other orbit senario?/Update when CI has launched.

  pred.rosetta           = flown.rosetta; % Rosetta EOL, no more predicted
  pred.juice             = 'juice_crema_5_0b23_1.tm'; % FIXME: Some other orbit senario?/Update when JUICE has launched.
  pred.bepicolombo       = 'bc_plan_*.tm'; % BepiColombo predicted
  pred.parkersolarprobe  = 'test*.tm'; % FIXME: UPDATE WHEN PSP sync script is tested
  pred.solarorbiter      = 'solo_ANC_soc-pred-mk_v*.tm'; % SolO predicted
  pred.cometinterceptor  = 'interceptor_study_v02.tm'; % FIXME: Some other orbit senario?/Update when CI has launched.

  switch p.Results.flown_or_predicted
    case 'predicted'
      srcMKfile = dir([mkPath, filesep, pred.(p.Results.mission)]);
    case 'flown'
      srcMKfile = dir([mkPath, filesep, flown.(p.Results.mission)]);
  end
  if size(srcMKfile, 1) > 1
    % Multiple kernels could be found if executing this script at the same
    % time as syncing new kernel files
    error('Found multiple metakernels, please check your SPICE folder.');
  elseif isempty(srcMKfile)
    irf.log('warning', 'Did not find any SPICE metakernels.');
    return
  end
  metakernel = adoptMetakernel(srcMKfile.folder, srcMKfile.name);

  cspice_furnsh(metakernel);
  % Keep its named stored (along with fields of what we have stored)
  LoadedSpiceKernels.(p.Results.mission).(p.Results.flown_or_predicted) = metakernel;
else
  irf.log('debug', 'It appears we have already loaded a kernel for requested "mission" and "flown_or_predicted" status.');
end


%% Help function
  function rootKernelPath = get_spice_root(mission)
    narginchk(1,1);
    % Load previously stored paths for mission in question
    rootKernelPath = datastore('spice_paths', mission);
    % Check if empty (ie. not defined previously) or if it no longer exists on
    % system.
    if isempty(rootKernelPath) || ~exist(rootKernelPath, 'dir')
      % It was not defined, ask user for path
      prompt = ['Please specify root path to where SPICE kernels for ', ...
        mission, ' is found. [i.e. "/share/SPICE/rosetta/"]: '];
      inputPath = input(prompt, 's');
      % Check it is a path and store it
      if exist(inputPath, 'dir')
        datastore('spice_paths', mission, inputPath);
        rootKernelPath = inputPath;
      else
        irf.log('critical', 'The path you specified did not exits.');
      end
    end
  end

  function localMKfile = adoptMetakernel(srcPath, srcFile)
    % Adopt a local metakernel to local mounted paths
    srcKernel = fileread([srcPath, filesep, srcFile]); % Read remote kernel file
    % Local temporary file (placed in "/tmp/" on Unix systems), in which we
    % can correct the local paths to the various SPICE kernel files.
    localMKfile = tempname;
    % Full local absolute path (metakernel are in subfolder "mk", so up one level)
    localSpicePath = what([srcPath, filesep, '..']);

    % Replace relative or remote paths with absolute path on locally mounted
    % system
    expression = 'PATH_VALUES\s{0,}=\s{0,}(\s{0,}''[a-zA-Z_0-9\.\/]*''\s{0,}';
    localKernel = regexprep(srcKernel, ...
      expression, ...
      ['PATH_VALUES = ( ''', localSpicePath.path, ''' ']);
    if ispc
      % FIXME:
      % Windows system use "\" as filesep, Unix use "/", check OS and if PC then
      % change all of the lines like:
      % KERNELS_TO_LOAD   = (
      %           '$KERNELS/ck/solo_ANC_soc-sc-iboom-ck_20180930-21000101_V01.bc'
      % to use "\" instead.
      % FIXME (second point):
      % Windows use EOL "CRLF" while Unix systems and most official mission repos have "LF" only.
      % SPICE require native OS EOL char for non-binary (i.e. text) kernel files.
      irf.log('critical', 'NOTE: Windows File separations in the metakernel is not yet implemented');
    end
    % Write the local metakernel to file ("cspice_furnsh" do not like it as a
    % string nor a file with relative paths)
    fileID = fopen(localMKfile, 'w');
    fwrite(fileID, localKernel);
    fclose(fileID);
  end

end
