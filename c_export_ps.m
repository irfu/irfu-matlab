function c_export_ps(st,sc_list,varargin)
%C_EXPORT_PS export figures 1:4 into PS and PDF
%
% function c_export_ps(st,[sc_list],[option,value])
% Input:
%	st - isdat epoch, defines filename YYYYMMDD_HHMM
%	sc_list - list of sc [optional], default 1:4
%	option [optional]:
%	sp - storage path
%	suf - suffix to add // file name becomes YY...HMM_SUF
% 
% Example:
%	c_export_ps(toepoch([2002 03 04 10 00 00]),[2 4],'suf','zoom1')
%
% See also IRF_FNAME
%
% $Id$

% Copyright 2004 Yuri Khotyaintsev (yuri@irfu.se)

if nargin<2, sc_list = 1:4; end
if nargin>2, have_options = 1; args = varargin;
else, have_options = 0;
end

sp = '.';
suf = '';

while have_options
	l = 2;
	if length(args)>1
		switch(args{1})
		case 'sp'
			if ischar(args{2}), sp = args{2};
			else, irf_log('fcal','wrong ArgType : sp must be string')
			end
		case 'suf'
			if ischar(args{2}), suf = ['_' args{2}];
			else, irf_log('fcal','wrong ArgType : suf must be string')
			end
		otherwise
        	irf_log('fcal',['Option ''' args{i} '''not recognized'])
    	end
		if length(args) > l, args = args(l+1:end);
		else break
		end
	else
		error('caa:wrongArgType','use c_export_ps(..,''option'',''value'')')
	end
end

old_pwd=pwd;
cd(sp)

for cl_id=sc_list
	figure(cl_id)
	fn = sprintf('EFW_C%d_%s%s',cl_id,irf_fname(st),suf);
	irf_log('save',['saving ' fn])
	print( gcf, '-dpsc2', fn) 
	unix(['/usr/local/bin/ps2pdf12 ' fn '.ps']);
end

cd(old_pwd)
