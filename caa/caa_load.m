function c=caa_load(varargin)
%Script to load data downloaded from the CAA in CDF format.
%Downloaded zip file must be unpacked, and script must be run 
%from a CAA_Download_YYYYMMDD_hhmm directory.
% Examples:
% caa_load
% caa_load CIS 
% caa_load C1
% $Id$

% ----------------------------------------------------------------------------
% "THE BEER-WARE LICENSE" (Revision 42):
% <yuri@irfu.se> wrote this file.  As long as you retain this notice you
% can do whatever you want with this stuff. If we meet some day, and you think
% this stuff is worth it, you can buy me a beer in return.   Yuri Khotyaintsev
% ----------------------------------------------------------------------------

flag_filter = 0;
if nargin > 0, % filter which variables to load 
    i=1;clear variable_filter;
    for j=1:length(varargin),
        if ischar(varargin{j}),
            variable_filter{i}=varargin{j};
            i=i+1;flag_filter=1;
        end
    end
end

nloaded = 0;
dirs = dir;
for j = 1:numel(dirs)
	if regexp(dirs(j).name,'^C[1-4]_(C|P)P_')
        var_name = dirs(j).name;
        d3 = findstr(var_name,'-');
        if d3, var_name( d3 ) = '_'; end
        flag_load_variable=1;
        if flag_filter==1, % if there is name filtering required check if to load variable
            for jj=1:length(variable_filter),
                if isempty(findstr(variable_filter{jj},var_name)),
                    flag_load_variable=0;
                end
            end
        end
        if flag_load_variable,
            try
                disp(['loading ' dirs(j).name filesep '*.cdf']);
                evalin('caller',[var_name '=dataobj(''' dirs(j).name filesep '*.cdf'');']);
                nloaded = nloaded + 1;
            catch
                disp(['error loading ' var_name]);
            end
        end
	end
end
if nloaded, disp(sprintf('CAA_LOAD : loaded %d variables',nloaded));
else disp('CAA_LOAD : nothing to load')
end