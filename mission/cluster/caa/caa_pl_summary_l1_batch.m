function caa_pl_summary_l1_batch(iso_st,iso_et,sp_dir)
%CAA_PL_SUMMARY_L1_BATCH: generate CAA summary plots for a given interval.
%
% caa_pl_summary_l1_batch(iso_st,iso_et,[sp_dir])
% Inputs:
%   iso_st,iso_et: start and end times for the interval (ISO format)
%   sp_dir: directory in which to put the plots. Default: /data/caa/sp/YYYY
%
% Example:
%  caa_pl_summary_l1_batch('2007-01-01T00:00:00Z', '2007-02-01T00:00:00Z')
%
% Plots always span 3 hours, although the interval may be as long as you
% wish.
narginchk(2,3)
old_pwd=pwd;
if nargin<3
	default_dir=1;
else
	cd(sp_dir);
	default_dir=0;
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
	if isempty(dir([l1dir '/C*']))
		figure(78)
		clf
		text(0.2,0.5,['No data for '  epoch2iso(t)]);
		fn = sprintf('EFW_SPLOT_L1ESPEC__%s',irf_fname(t));
		fne = sprintf('EFW_SPLOT_L1ERSPEC__%s',irf_fname(t));
		fnq = sprintf('EFW_SPLOT_L1QUAL__%s',irf_fname(t));
		fone = sprintf('EFW_SPLOT_L1__%s',irf_fname(t));
		irf_log('save','saving  blank pdf plots.')
		print( 78, '-dpdf', fn), print( 78, '-dpdf', fne), print( 78, '-dpdf', fnq)
		if exist('/usr/local/bin/pdfjoin','file')
			irf_log('save',['joining to ' fone '.pdf'])
			s = unix(['LD_LIBRARY_PATH="" /usr/local/bin/pdfjoin ' fn '.pdf ' fne '.pdf ' fnq '.pdf --outfile ' fone '.pdf']);
			if s~=0, irf_log('save','problem with pdfjoin'), end
		else
			irf_log('proc',...
				'cannot join PDFs: /usr/local/bin/pdfjoin does not exist')
		end
    else
%l1dir %nl
        di=l1dir(length(l1dir)-12:end);
        it=[di(1:4) '-' di(5:6) '-' di(7:8) 'T' di(10:11) ':' di(12:13) ':00Z'];
		caa_pl_summary_l1(it,10800,l1dir,'savepdf')
%		caa_pl_summary_l1(-1,-1,l1dir,'savepdf')
	end
end

cd(old_pwd);