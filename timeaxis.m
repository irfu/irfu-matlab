function res = timeaxis(limit)
% res = timeaxis(limit)
%   returns tick mark locations and labels for a time axis
%   in a cell array. limit is a 2-element vector in seconds.

dtime = limit(2)-limit(1);

% get the time difference between two ticks and the number of unlabeled ticks
if dtime>3600*24*500
  dticv = 3600*24*100;
  mtics = 2;
elseif dtime>3600*24*200
  dticv = 3600*24*50;
  mtics = 5;
elseif dtime>3600*24*100
  dticv = 3600*24*20;
  mtics = 4;
elseif dtime>3600*24*50
  dticv = 3600*24*10;
  mtics = 2;
elseif dtime>3600*24*20
  dticv = 3600*24*5;
  mtics = 5;
elseif dtime>3600*24*10
  dticv = 3600*24*2;
  mtics = 2;
elseif dtime>3600*24*5
  dticv = 3600*24;
  mtics = 4;
elseif dtime>3600*24*2
  dticv = 3600*6;
  mtics = 3;
elseif dtime>3600*24
  dticv = 3600*4;
  mtics = 4;
elseif dtime>3600*12
  dticv = 3600*2;
  mtics = 2;
elseif dtime>3600*5
  dticv = 3600;
  mtics = 6;
elseif dtime>=3600*2
  dticv = 1800;
  mtics = 3;
elseif dtime>3600
  dticv = 1200;
  mtics = 2;
elseif dtime>60*30
  dticv = 600;
  mtics = 5;
elseif dtime>60*10
  dticv = 300;
  mtics = 5;
elseif dtime>60*3
  dticv = 60;
  mtics = 6;
elseif dtime>60*2
  dticv = 30;
  mtics = 3;
elseif dtime>60
  dticv = 20;
  mtics = 4;
elseif dtime>30
  dticv = 15;
  mtics = 3;
elseif dtime>20
  dticv = 10;
  mtics = 5;
elseif dtime>10
  dticv = 5;
  mtics = 5;
elseif dtime>6
  dticv = 2;
  mtics = 4;
elseif dtime>3
  dticv = 1.;
  mtics = 5.;
elseif dtime>1
  dticv = 0.5;
  mtics = 5;
elseif dtime>.6
  dticv = .2;
  mtics = 4;
elseif dtime>.3
  dticv = .1;
  mtics = 5;
elseif dtime>.1
  dticv = .05;
  mtics = 5;
elseif dtime>.06
  dticv = .03;
  mtics = 6;
elseif dtime>.04
  dticv = .02;
  mtics = 4;
elseif dtime>.02
  dticv = .01;
  mtics = 5;
elseif dtime>.01
  dticv = .005;
  mtics = 5;
elseif dtime>.006
  dticv = .003;
  mtics = 6;
elseif dtime>.004
  dticv = .002;
  mtics = 4;
elseif dtime>.002
  dticv = .001;
  mtics = 5;
elseif dtime>.001
  dticv = .0005;
  mtics = 5;
else
  dticv = 0.0001;
  mtics = 10;
end
%disp(['dticv=' num2str(dticv) ', mtics=' num2str(mtics)]);
% calculate the time value of the first major tick
% tbeg = dticv*ceil(limit(1)/dticv);

% calculate the time value of the first minor tick
dmort = dticv/mtics;
tbeg = dmort*ceil(limit(1)/dmort);

% calculate the number of ticks
ttic = tbeg;
ntics = 0;
while ttic<=limit(2)
%  ttic = ttic+dticv;
  ttic = ttic+dmort;
  ntics = ntics+1;
  if ntics>100
    warning, 'too many ticks in timeaxis'
    break;
  end
end

% generate array with the time values of the major ticks
% tictv = tbeg + dticv.*[0:ntics-1];
% generate array with the time values of the ticks
tictv = tbeg + dmort*[0:ntics-1];

% generate the time strings for the labels,
ticstr = cell(1, ntics);
ticval = mod(tictv, 86400);
% use the long format hh:mm:ss, if more than one label within one second,
%     else use hh:mm
%n = find(tictv-floor(tictv)<2e-7);
n=1:ntics;
hour = floor(ticval(n)/3600);
minute = floor(mod(ticval(n), 3600)/60);
if dticv>=3600*24
  hhmmss = 'datestr(datenum(fromepoch(tictv(j))),1)';
elseif dticv>=60
  format = '%02d:%02d';
  hhmmss = [hour; minute];
elseif dticv<.001
  format = '%02d:%02d:%07.4f';
  hhmmss = [hour; minute; mod(ticval(n), 60)];
elseif dticv<.01
  format = '%02d:%02d:%06.3f';
  hhmmss = [hour; minute; mod(ticval(n), 60)];
elseif dticv<.1
  format = '%02d:%02d:%05.2f';
  hhmmss = [hour; minute; mod(ticval(n), 60)];
elseif dticv<1
  format = '%02d:%02d:%04.1f';
  hhmmss = [hour; minute; mod(ticval(n), 60)];
else
  format = '%02d:%02d:%02.0f';
  hhmmss = [hour; minute; mod(ticval(n), 60)];
end

ind_labels=find(abs(mod(tictv,dticv))<1e-6); % NOTE does not work for tick distance below a few microseconds
for j=n,ticstr{j} = ' ';end
if dticv>=3600*24,
  for j=ind_labels, ticstr{j} = eval(hhmmss);end
else,
  for j=ind_labels, ticstr{j} = sprintf(format, hhmmss(:,j));end
end

if dticv>=.1,
  ind_ms_labels=find(abs(mod(tictv(ind_labels),1))>1e-6);
  if length(ind_ms_labels) < length(ind_labels),
  for j=ind_labels(ind_ms_labels), ticstr{j} = sprintf('.%01.0f', mod(tictv(j),1)*10);end
  end
elseif dticv>=.01,
  ind_ms_labels=find(abs(mod(tictv(ind_labels),.1))>1e-6);
  if length(ind_ms_labels) < length(ind_labels),
  for j=ind_labels(ind_ms_labels), ticstr{j} = sprintf('.%02.0f', mod(tictv(j),1)*100);end
  end
end
res{1} = tictv;
res{2} = ticstr;
