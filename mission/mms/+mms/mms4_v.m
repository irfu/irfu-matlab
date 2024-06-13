function v = mms4_v(t, coord_sys)
%MMS4_V   Calculate velocity from timing between 4 spacecraft
% NOT VALID FOR Cluster! (use c_v instead)
%
% v  = MMS4_V(t, [coord_sys]);
% dt = MMS4_V(v, [coord_sys]);
%
% Calculate velocity from timing between 4 spacecraft
% or calculate timing from the velocity
%
%   Inputs:
% t  = [t1, t2, t3, t4]; in EpochTT
% v  = TSeries of [vx, vy, vz]
%   ie: v = irf.ts_vec_xyz(t, [vx, vy, vz]); t in EpochTT and v in GSE ref frame.
% dt = [0, t2-t1, t3-t1, t4-t1];
%
% coord_sys = 'gsm'; % when calculate in GSM reference frame
%
%   See also c_v

% min #args=1, max #args=2
narginchk(1,2);

% R and V will retain their value between function calls
persistent R V

% set R and V if they don't exist or are empty
% struct() will create an array where R.time corresponds to the next
% element. i.e. s = struct(field1,value1,...,fieldN,valueN)
if ~exist('R','var') || isempty(R)
  R = struct('time',[], 'gseR1',[], 'gseR2',[], 'gseR3',[], 'gseR4',[], ...
    'gsmR1',[], 'gsmR2',[], 'gsmR3',[], 'gsmR4',[]);
end
if ~exist('V','var') || isempty(V)
  V = struct('time',[], 'gseV1',[], 'gseV2',[], 'gseV3',[], 'gseV4',[], ...
    'gsmV1',[], 'gsmV2',[], 'gsmV3',[], 'gsmV4',[]);
end

% i.e. GSE are default coords
if nargin==1, coord_sys = 'gse'; end

if isa(t, 'TSeries')
  flag = 'dt_from_v';
  tint = irf.tint([t.time.tts-61, t.time.tts+61]);
  vOrig = [t.time.epochUnix, t.data];
  if strcmpi(coord_sys, 'gsm') && (isempty(t.coordinateSystem) || ...
      strcmpi(t.coordinateSystem, 'gse'))
    v = irf_gse2gsm(vOrig, -1);
  else
    v = vOrig;
  end
  t = v(1);
elseif isa(t, 'GenericTimeArray')
  flag = 'v_from_t';
  tint = irf.tint([min(t.tts)-61, max(t.tts)+61]);
  t = t.epochUnix;
else
  errStr='Unexpected input';
  irf.log('critical', [errStr,'s in t: ', t]);
  error(errStr);
end

% Begin by looking for locally saved mat files of positon and velocity to load (quickest)
if ~is_R_ok && exist(['.' filesep, 'mmsR.mat'],'file') && ...
    exist(['.' filesep, 'mmsV.mat'],'file')
  % Local files found
  load(['.' filesep, 'mmsR.mat']);
  load(['.' filesep, 'mmsv.mat']);
  irf.log('warning',['MMS position loaded from current folder (.',filesep, 'mmsR.dat)']);
end

if ~is_R_ok && exist([filesep, 'data', filesep, 'mms', filesep, 'irfu', filesep, 'mmsR.mat'], 'file') && ...
    exist([filesep, 'data', filesep, 'mms', filesep, 'irfu', filesep, 'mmsV.mat'], 'file')
  % Local files found
  load([filesep, 'data', filesep, 'mms', filesep, 'irfu', filesep, 'mmsR.mat']);
  load([filesep, 'data', filesep, 'mms', filesep, 'irfu', filesep, 'mmsV.mat']);
  irf.log('warning',['MMS position loaded from local repository folder (', ...
    filesep, 'data', filesep, 'mms', filesep, 'irfu', filesep, 'mmsR.mat']);
end

if ~is_R_ok && exist('mmsR.mat', 'file') && exist('mmsV.mat', 'file')
  % Local files found
  load('mmsR.mat');
  load('mmsV.mat');
  irf.log('warning','MMS position loaded from somewhere in the path (mmsR.mat)')
end

% In case no locally saved mat files found, try to read ancillary data
if ~is_R_ok
  R_gse = mms.get_data('R_gse', tint);
  V_gse = mms.get_data('V_gse', tint);
  R_gsm = mms.get_data('R_gsm', tint); %#ok<NASGU>
  V_gsm = mms.get_data('V_gsm', tint); %#ok<NASGU>
  V.time = V_gse.time;
  c_eval('V.gseV? = V_gse.gseV?;');
  c_eval('V.gsmV? = V_gsm.gsmV?;');
  R.time = R_gse.time;
  c_eval('R.gseR? = R_gse.gseR?;');
  c_eval('R.gsmR? = R_gsm.gsmR?;');
end

if ~is_R_ok
  irf.log('warning','!!! Could not obtain position data !!!');
  return
end

%remove NaN
ind = isnan([R.gseR1(:,1), R.gseR2(:,1), R.gseR3(:,1), R.gseR4(:,1)]);
[~, col] = max(sum(ind, 1)); % Most restrictive
c_eval('R.gseR?(ind(:,col), :) = [];');
c_eval('R.gsmR?(ind(:,col), :) = [];');
R.time = R.time(~ind(:,col));

ind = isnan([V.gseV1(:,1), V.gseV2(:,1), V.gseV3(:,1), V.gseV4(:,1)]);
[~, col] = max(sum(ind, 1)); % Most restrictive
c_eval('V.gseV?(ind(:,col), :) = [];');
c_eval('V.gsmV?(ind(:,col), :) = [];');
V.time = V.time(~ind(:,col));

t_center = 0.5*t(1)+0.5*t;

if strcmp(flag, 'v_from_t')
  for ic='1234'
    i = ic-'0';
    R.(['vsc', ic]) = irf_resamp([V.time.epochUnix, ...
      V.([coord_sys, 'V', ic])], ...
      t_center', 'spline');
    R.(['drsc',ic]) = irf_resamp([R.time.epochUnix, ...
      R.([coord_sys, 'R', ic]) - R.([coord_sys, 'R1'])], ...
      t(i), 'spline');
    R.(['dr', ic])  = R.(['drsc', ic]) + [0, (t(i)-t(1))*R.(['vsc', ic])(1, 2:4)];
    R.dt(i)         = t(i) - t(1);
    R.(['sdt', ic]) = num2str(R.dt(i), 3);
  end
  D = [R.dr2(2:4); R.dr3(2:4); R.dr4(2:4)];
  T = [R.dt(2), R.dt(3), R.dt(4)]';
  m = D\T;
  clear v
  v = m/norm(m)/norm(m);
  v = v';	% velocity vector of the boundary
  disp(EpochUnix(t(1)).toUtc)
  vn = irf_norm(v);
  fprintf('dt=[%5.2f, %5.2f, %5.2f, %5.2f] s. dt=[t1-t1, t2-t1, ...]\n', R.dt);
  fprintf('V=%3.2f [%5.2f %5.2f %5.2f] km/s %s\n', irf_abs(v,1), vn(end-2:end), coord_sys);

elseif strcmp(flag, 'dt_from_v')
  for ic='1234'
    i = ic-'0';
    R.(['vsc', ic]) = irf_resamp([V.time.epochUnix, ...
      V.([coord_sys, 'V', ic])], ...
      t_center', 'spline');
    R.(['v', ic])  = v(2:4) - dot(R.(['vsc', ic])(2:4), v(2:4)) .* v(2:4) ./ norm(v(2:4))^2;
    R.(['dr', ic]) = irf_resamp([R.time.epochUnix, ...
      R.([coord_sys, 'R', ic]) - R.([coord_sys, 'R1'])], ...
      t, 'spline');
    R.dt(i) = irf_dot(R.(['dr', ic]), R.(['v', ic]), 1) ./ norm(R.(['v', ic]))^2;
  end
  % print result
  disp(EpochUnix(t(1)).toUtc)
  vn = irf_norm(vOrig);
  fprintf('V=%3.2f*[ %5.2f %5.2f %5.2f] km/s %s\n', irf_abs(vOrig,1), vn(end-2:end), coord_sys);
  fprintf('dt=[%5.2f, %5.2f, %5.2f, %5.2f] s. dt=[t1-t1, t2-t1, ...]\n', R.dt);
  v = R.dt; % output variable is v
end

% LOCAL FUNCTION
  function answer = is_R_ok(sc)
    % check if position data are ok for spacecraft number 'sc'
    % if input argument not given check if ok for all spacecraft that needs
    % to be plotted.
    if nargin == 0, scList = 1:4; else, scList = sc; end
    for iSc=scList
      if numel(R.([coord_sys, 'R', num2str(iSc)])) < 8 % less than 2 time points
        answer = false; return
      else
        if(R.time(1).epochUnix>min(t)) || (R.time(end).epochUnix<max(t))
          answer = false; return
        end
      end
    end
    answer = true;
  end
end