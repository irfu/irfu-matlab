function data=readWBD(start_time,dt,cl_id)
%readWBD read WBD data using R_WBD_WF
% http://www-pw.physics.uiowa.edu/plasma-wave/istp/cluster/dvd/SOFTWARE/
%
% data = readWBD(start_time,dt,cl_id)
%
% Input:
%	start_epoch - start time (isdat epoch)
%	dt - time interval in sec
%	cl_id - SC#
%
% Output:
%	data - Cluster AV format [t_epoch data()]
%
% $Id$
%
% see also C_GET

warning('caa:cleanup',...
'Function %s is deprecated and will be removed on May 1, 2004.\nUse %s instead',
...
mfilename,'c_wbd_read')

% Copyright 2004 Yuri Khotyaintsev

error(nargchk(3,3,nargin))

cl_id_s = num2str(cl_id);

[s,m] = unix('uname');

R_WBD_WF = ['/data/cluster/WBD/bin/' m(2:end-1) '/r_wbd_wf'];
tmp_file = tempname;

ts_s = c_epoch2str(start_time);
te_s = c_epoch2str(start_time + dt);

ux_cmd = [R_WBD_WF ' ' cl_id_s ' 0 0 ' ts_s ' ' te_s '>' tmp_file];
[s,m] = unix(ux_cmd);
if s~=0, error('Error running R_WBD_WF'), end
c_log('dsrc',m)

fid=fopen(tmp_file);
if fid==-1, error(['cannot read temporary file ' tmp_file]), end
	
data = fread(fid,[2 inf],'float32');
fclose(fid);
unix(['rm ' tmp_file]);

data = data';
data(:,1) = data(:,1) + start_time;

c_log('dsrc',['readWBD: ' num2str(size(data,1)) ' data points'])
