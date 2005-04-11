function caa_export_batch(sc_list)
%CAA_EXPORT_BATCH run CAA_EXPORT in a batch script
%
% caa_export_batch([sc_list])
%
% $Id$

% Copyright 2004,2005 Yuri Khotyaintsev

if nargin<1, sc_list=1:4; end

if exist('./mInfo.mat','file')
	load -mat mInfo caa_q
end
if ~exist('caa_q','var')
	disp('cannot load quality information from mInfo.mat')
	disp('please run_caa_quality')
	return
end

l1 = {'P1','P2','P3','P4','P12','P32','P34'};
l2 = {'P','E'};

for j=sc_list
	disp(['Cluster ' num2str(j)])
	for k=1:length(l1)
		disp(['Level 1 : ' l1{k}])
		if k<=4, q = caa_q.p(j,1);
		else, q = caa_q.e(j,1);
		end
		caa_export(1,l1{k},j,q) 
	end
	for lev=2:3, for k=1:length(l2)
		disp(['Level ' num2str(lev) ' : ' l2{k}])
		if k==1, q = caa_q.p(j,lev);
		else, q = caa_q.e(j,lev);
		end
		caa_export(lev,l2{k},j,q)
	end, end
end
