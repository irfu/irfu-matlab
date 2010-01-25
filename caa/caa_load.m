%Script to load data downloaded from the CAA in CDF format.
%Downloaded zip file must be unpacked, and script must be run 
%from a CAA_Download_YYYYMMDD_hhmm directory.
%
% $Id$

% ----------------------------------------------------------------------------
% "THE BEER-WARE LICENSE" (Revision 42):
% <yuri@irfu.se> wrote this file.  As long as you retain this notice you
% can do whatever you want with this stuff. If we meet some day, and you think
% this stuff is worth it, you can buy me a beer in return.   Yuri Khotyaintsev
% ----------------------------------------------------------------------------

nloaded = 0;
dirs = dir;
old_pwd = pwd;
for j = 1:numel(dirs)
	if regexp(dirs(j).name,'^C[1-4]_(C|P)P_')
        var_name = dirs(j).name;
        d3 = findstr(var_name,'-');
        if d3, var_name( d3 ) = '_'; end
		disp(['loading ' var_name]);
		try
			cd(dirs(j).name)
			eval([var_name '=dataobj(''*.cdf'');'])
			nloaded = nloaded + 1;
		catch
			disp(['error loading ' var_name]);
		end
		cd(old_pwd)
	end
end
if nloaded, disp(sprintf('CAA_LOAD : loaded %d variables',nloaded));
else disp('CAA_LOAD : nothing to load')
end
clear nloaded dirs old_pwd