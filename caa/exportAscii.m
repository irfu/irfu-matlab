function exportAscii(var,vs)
%exportAscii(var,var_name_s)

if strcmp(vs,'wE1p12') | strcmp(vs,'wE1p12')
	cl_id = vs(3);
	inst = 'EFW';
	sig = 'E';
	sen = vs(4:6);
	frame = 'SC';
	var_s = ['time [sec]  E' sig];
end

t0_s = datestr(datenum(fromepoch(var(1,1))));
var(:,1) = var(:,1) - var(1,1); 

fid = fopen([vs '.dat'],'w');
fprintf(fid,['%% SC:        Cluster ' cl_id' \n']);
fprintf(fid,['%% Intrument: ' inst ' \n']);
fprintf(fid,['%% Signal:    ' sig ' \n']);
fprintf(fid,['%% Sensor:    ' sen ' \n']);
fprintf(fid,['%% Coord Sys: ' frame ' \n']);
fprintf(fid,['%% Time from: ' t0_s ' \n']);
fprintf(fid,['%% ' var_s ' \n']);

