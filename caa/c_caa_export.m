function c_caa_export(sc_list)

warning('caa:cleanup',...
'Function %s is deprecated and will be removed on May 1, 2004.\nUse %s instead',...
mfilename,'caa_export_batch')

l1 = {'P1','P2','P3','P4','P12','P32','P34'};
l2 = {'P','E'};

for j=sc_list
	for k=1:length(l1), disp(['Level 1 : ' l1{k}]), export2caa(1,l1{k},j), end
	for lev=2:3, for k=1:length(l2)
		disp(['Level ' num2str(lev) ' : ' l2{k}])
		export2caa(lev,l2{k},j)
	end, end
end
