function caa_sp_batch(fname)
%CAA_SP_BATCH  produce summary plots
%
% caa_sp_batch(fname)
%
% $Id$

% Copyright 2007 Yuri Khotyaintsev

old_pwd = pwd;
BASE_DIR = '/data/caa/l1';
dirs = textread(fname,'%s');

if isempty(dirs), disp('NO DIRS'), cd(old_pwd), return, end

for d=1:length(dirs)
	s = dirs{d};
	caa_pl_summary_l1(...
		sprintf('%s-%s-%sT%s:00:00Z',s(1:4),s(10:11),s(12:13),s(15:16)),...
		3*3600,	[BASE_DIR '/' s],'saveps')
end
