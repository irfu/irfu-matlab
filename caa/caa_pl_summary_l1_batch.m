function caa_pl_summary_l1_batch(iso_st,iso_et,dir)
%CAA_PL_SUMMARY_L1_BATCH: generate CAA summary plots for a given interval.
%
% caa_pl_summary_l1_batch(iso_st,iso_et,[dir])
% Inputs:
%   iso_st,iso_et: start and end times for the interval (ISO format)
%   dir: directory in which to put the plots. Default: /data/caa/sp/YYYY
%
% Example:
%  caa_pl_summary_l1_batch('2007-01-01T00:00:00Z', '2007-02-01T00:00:00Z')
%
% Plots always span 3 hours, although the interval may be as long as you
% wish.
error(nargchk(2,3,nargin))
old_pwd=pwd;
if nargin<3
	default_dir=1;
else
	cd(dir);
	defualt_dir=0;
end
BASE_DIR='/data/caa/l1';

t0=iso2epoch(iso_st);
t1=iso2epoch(iso_et);
t0=t0-mod(t0,3*3600);
t1=t1+mod(t0,3*3600);
for t=t0:3*3600:t1
	y=fromepoch(t);
	l1dir = [BASE_DIR '/' num2str(y(1)) '/' irf_fname(t)];
	if default_dir, cd(['/data/caa/sp/' num2str(y(1))]), end
	caa_pl_summary_l1(-1,-1,l1dir,'savepdf')
end

cd(old_pwd);