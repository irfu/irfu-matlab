function caa_export_batch_l1(cl_id,out_path)
%CAA_EXPORT_BATCH_L1 run CAA_EXPORT in a batch script
%
% caa_export_batch_l1(cl_id,[out_path])
%
% $Id$

% Copyright 2005 Yuri Khotyaintsev

if nargin<2, out_path='.'; end

QUAL=3;

l1p = {'P1','P2','P3','P4'};
l1e = {'P12','P32','P34'};

if ~(exist('./mPR.mat','file') | exist('./mER.mat','file'))
	c_log('load','no data'), return
end

v = caa_read_version;
if isempty(v), v = 0; end
v = v + 1;
v_s = num2str(v);
if v<10, v_s = ['0'  v_s]; end
[s,w] = unix(...
	['echo "' v_s ' ' epoch2iso(date2epoch(now),1) '">>.version'],'-echo');
if s~=0, irf_log('save','problem updating version'),cd(old_pwd),return, end
disp(['new version is : ' v_s])

sp = pwd;
cd(out_path)

% L1
if exist([sp '/mPR.mat'],'file')
	for k=1:length(l1p)
		disp(['Level 1 : ' l1p{k}])
		caa_export(1,l1p{k},cl_id,QUAL,v_s,sp) 
	end
end
if exist([sp '/mER.mat'],'file')
	for k=1:length(l1e)
		disp(['Level 1 : ' l1e{k}])
		caa_export(1,l1e{k},cl_id,QUAL,v_s,sp) 
	end
end

% EF
if exist([sp '/mEDSIf.mat'],'file')
	disp('Level 2 : EF')
	caa_export(2,'EF',cl_id,QUAL,v_s,sp)
end

% P
if exist([sp '/mP.mat'],'file')
	for lev=2:3
		disp(['Level ' num2str(lev) ' : P'])
		caa_export(lev,'P',cl_id,QUAL,v_s,sp)
	end
end

cd(sp)
