function data=c_wbd_read(start_time,dt,cl_id)
%C_WBD_READ   Read WBD data using R_WBD_WF
%
% http://www-pw.physics.uiowa.edu/plasma-wave/istp/cluster/dvd/SOFTWARE/
%
% data = c_wbd_read(start_time,dt,cl_id)
%
% Input:
%	start_epoch - start time (isdat epoch)
%	dt - time interval in sec
%	cl_id - SC#
%
% Output:
%	data - Cluster IRFU format [t_epoch data()]
%
% see also C_GET
%

% Copyright 2004 Yuri Khotyaintsev

narginchk(3,3)

cl_id_s = num2str(cl_id);

[s,m] = unix('uname');

DP_S = c_ctl(0,'data_path');
R_WBD_WF = [DP_S '/WBD/bin/' m(2:end-1) '/r_wbd_wf'];
tmp_file = tempname;

ts_s = epoch2iso(start_time,1);
te_s = epoch2iso(start_time + dt,1);

ux_cmd = [R_WBD_WF ' ' cl_id_s ' 0 0 ' ts_s ' ' te_s '>' tmp_file];
[s,m] = unix(ux_cmd);
if s~=0, error('Error running R_WBD_WF'), end
irf_log('dsrc',m)

fid=fopen(tmp_file);
if fid==-1, error(['cannot read temporary file ' tmp_file]), end

data = fread(fid,[2 inf],'float32');
fclose(fid);
unix(['rm ' tmp_file]);

data = data';
data(:,1) = data(:,1) + start_time;

irf_log('dsrc',['c_wbd_read: ' num2str(size(data,1)) ' data points'])
