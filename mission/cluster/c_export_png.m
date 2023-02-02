function c_export_png(st,sc_list,varargin)
%C_EXPORT_PNG export figures 1:4 into PNG
%
% function c_export_png(st,[sc_list],[option,value])
% Input:
%	st - isdat epoch, defines filename YYYYMMDD_HHMM
%	if st is omitted, we try to guess it from userdata of figure 1
%	sc_list - list of sc [optional], default 1:4
%	option [optional]:
%	sp - storage path
%	suf - suffix to add // file name becomes YY...HMM_SUF
%
% Example:
%	c_export_png(toepoch([2002 03 04 10 00 00]),[2 4],'suf','zoom1')
%
% See also C_EXPORT_PS, IRF_FNAME
%

% Copyright 2004,2005 Yuri Khotyaintsev (yuri@irfu.se)

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
        irf_log('fcal',['Option ''' args{1} '''not recognized'])
    end
    if length(args) > l, args = args(l+1:end);
    else, break
    end
  else
    error('caa:wrongArgType','use c_export_ps(..,''option'',''value'')')
  end
end

old_pwd=pwd;
cd(sp)
if nargin<1
  ud=get(1,'userdata');
  if isfield(ud,'t_start_epoch'), st=ud.t_start_epoch;
  else, error('please specify ts')
  end
end

for cl_id=sc_list
  figure(cl_id)
  fn = sprintf('EFW_C%d_%s%s',cl_id,irf_fname(st),suf);
  irf_log('save',['saving ' fn])
  print( gcf, '-depsc2', fn)
  [s,w] = unix(['/usr/local/bin/eps2png ' fn '.eps; rm -f ' fn '.eps']);
  if s~=0, irf_log('save','problem with eps2png'), end
end

cd(old_pwd)
