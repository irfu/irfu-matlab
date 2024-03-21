function caa_sh_3h_int(st,et)
%CAA_SH_3H_INT  create common 3h intervals
%
% caa_sh_3h_int(st,et)

if ischar(st), st = iso2epoch(st); end
if ischar(et), et = iso2epoch(et); end
if et-1 < st, error('interval too short'), end

tt = epoch2iso(st);
year = str2double(tt(1:4));

% Sanity check
tt = epoch2iso(et-1);
if year ~=str2double(tt(1:4)), error('ST and ET must be from the same year'), end

ORB =[]; MP1 = []; MP2 = []; MP3 = []; MP4 = [];

if exist('./mPlan.mat','file'), load ./mPlan.mat
else, error('No MPlan.mat found')
end

for cl_id = 1:4
  v_s = sprintf('MP%dY%d',cl_id,year);
  if ~exist(v_s,'var'), error(['Cannot load ' v_s]), end
  c_eval([ 'MP?=' v_s ';'],cl_id);
end

v_s = sprintf('ORB%dY%d',1,year);
if ~exist(v_s,'var'), error(['Cannot load ' v_s]), end
eval([ 'ORB=' v_s ';'])

irf_log('proc',['creating 3h intervals for ' epoch2iso(st,1) ' -- ' ...
  epoch2iso(et,1)])

ORB = ORB((ORB(:,1)+ORB(:,2))>=st & ORB(:,1)<et, :);
if isempty(ORB), error('no ORB for the requested interval'), end

MP = [];
for o=1:size(ORB,1)
  st = ORB(o,1);
  et = ORB(o,1)+ORB(o,2);
  t_out = min(MP1(MP1(:,1)>st & MP1(:,1)<et,1));
  t_out2 = min(MP2(MP2(:,1)>st & MP2(:,1)<et,1));
  t_out3 = min(MP3(MP3(:,1)>st & MP3(:,1)<et,1));
  t_out4 = min(MP4(MP4(:,1)>st & MP4(:,1)<et,1));
  if(isempty(t_out) || (~isempty(t_out2) && t_out2<t_out)), t_out=t_out2; end
  if(isempty(t_out) || (~isempty(t_out3) && t_out3<t_out)), t_out=t_out3; end
  if(isempty(t_out) || (~isempty(t_out4) && t_out4<t_out)), t_out=t_out4; end
  t_in   = max(MP1(MP1(:,2)>st & MP1(:,2)<et,2));
  t_in2  = max(MP2(MP2(:,2)>st & MP2(:,2)<et,2));
  t_in3  = max(MP3(MP3(:,2)>st & MP3(:,2)<et,2));
  t_in4  = max(MP4(MP4(:,2)>st & MP4(:,2)<et,2));
  if(isempty(t_in) || (~isempty(t_in2) && t_in2>t_in)), t_in=t_in2; end
  if(isempty(t_in) || (~isempty(t_in3) && t_in3>t_in)), t_in=t_in3; end
  if(isempty(t_in) || (~isempty(t_in4) && t_in4>t_in)), t_in=t_in4; end

  if ~isempty(t_out) && ~isempty(t_in)
    tt = fromepoch(t_out);
    t_out = toepoch([tt(1:3) fix(tt(4)/3)*3 0 0]);
    tt = fromepoch(t_in);
    t_in = toepoch([tt(1:3) 0 0 0]) + 3600*ceil(tt(4)/3)*3;
    MP = [MP; [t_out, t_in]];
  else
    irf_log('proc',['no MP for ' epoch2iso(st,1) ' -- ' ...
      epoch2iso(et,1)])
  end
end

if isempty(MP), irf_log('proc','no MP :('), end

v_s = sprintf('MPauseY%d',year);
eval([v_s '=MP; save ./mPlan.mat ' v_s ' -append'])
disp([v_s ' -> ./mPlan.mat'])
