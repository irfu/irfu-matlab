function caa_export_batch(sc_list)
%CAA_EXPORT_BATCH run CAA_EXPORT in a batch script
%
% caa_export_batch([sc_list])
%
% $Id$

% Copyright 2004,2005 Yuri Khotyaintsev

if nargin<1, sc_list=1:4; end

l1 = {'P1','P2','P3','P4','P12','P32','P34'};
l2 = {'P','E'};

for j=sc_list
	for k=1:length(l1), disp(['Level 1 : ' l1{k}]), caa_export(1,l1{k},j), end
	for lev=2:3, for k=1:length(l2)
		disp(['Level ' num2str(lev) ' : ' l2{k}])
		caa_export(lev,l2{k},j)
	end, end
end
