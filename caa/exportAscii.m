function exportAscii(var,vs)
%exportAscii export variable into ascii file with comments
% exportAscii(var,var_name_s)
%
% Information written into comment is based on var_name_s string
%
% Example:
%	exportAscii(diE3p1234,'diE3p1234')
%
% $Revision$  $Date$
%

% Copyright 2004 Yuri Khotyaintsev


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% raw E p12 and p34
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if regexp(vs,'wE[1-4]p(12|34)')==1

	cl_id = vs(3);
	inst = 'EFW';
	sig = 'E';
	sen = vs(4:6);
	frame = 'SC';
	var_labels = {['E' sen]};

	% remove averages
	cp = ClusterProc('./');
	var = corrADCOffset(cp,var);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% spin fits E p12 and p34
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elseif regexp(vs,'diEs[1-4]p(12|34)')==1

	cl_id = vs(5);
	inst = 'EFW';
	sig = 'E';
	sen = ['spin fits ' vs(6:8)];
	frame = 'DSI,  Ez==0 : not measured';
	var_labels = {'Ex','Ey','Ez'};

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% despun full resolution E
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elseif regexp(vs,'diE[1-4]p1234')==1

	cl_id = vs(4);
	inst = 'EFW';
	sig = 'E';
	sen = vs(5:9);
	frame = 'DSI,  Ez==0 : not measured';
	var_labels = {'Ex','Ey','Ez'};

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% phase
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elseif regexp(vs,'A[1-4]')

	cl_id = vs(2);
	inst = 'ephemeris';
	sig = 'phase';
	sen = '';
	frame = 'SC';
	var_labels = {'phase'};

else
	warning('Wariable name not recognized, will do nothing.')
	return
end

t0_s = datestr(datenum(fromepoch(var(1,1))));
var(:,1) = var(:,1) - var(1,1); 

sz = size(var);
n_data = sz(2) - 1; % number of data columns - time
var_s = 'time [sec] ';
mask = '%10.4f ';

if		strcmp(sig,'E'),		units = 'mV/m';
elseif	strcmp(sig,'phase'),	units = 'deg';
else,	units = 'undef';
end

for i=1:n_data
	var_s = [var_s var_labels{i} ' [' units '] '];
	mask = [mask '%8.3f '];
end

fid = fopen([vs '.dat'],'w');
fprintf(fid,['%% file created on ' date ' \n%%\n']);
fprintf(fid,['%% SC:        Cluster ' cl_id ' \n']);
fprintf(fid,['%% Intrument: ' inst ' \n']);
fprintf(fid,['%% Signal:    ' sig ' \n']);
fprintf(fid,['%% Sensor:    ' sen ' \n']);
fprintf(fid,['%% Coord Sys: ' frame ' \n%%\n']);
fprintf(fid,['%% Time from: ' t0_s ' \n%%\n']);
fprintf(fid,['%% ' var_s ' \n']);
fprintf(fid,[mask '\n'],var');
fclose(fid);

