function caa_load(varargin)
%Script to load data downloaded from the CAA in CDF format.
%Downloaded zip file must be unpacked, and script must be run 
%from a CAA_Download_YYYYMMDD_hhmm directory.
%  Examples:
%   caa_load
%   caa_load CIS 
%   caa_load C1 CIS
% $Id$

% ----------------------------------------------------------------------------
% "THE BEER-WARE LICENSE" (Revision 42):
% <yuri@irfu.se> wrote this file.  As long as you retain this notice you
% can do whatever you want with this stuff. If we meet some day, and you think
% this stuff is worth it, you can buy me a beer in return.   Yuri Khotyaintsev
% ----------------------------------------------------------------------------

flag_filter = 0;
if nargin > 0, % filter which variables to load 
    i=1;
    variable_filter=cell(nargin,1);
    for j=1:length(varargin),
        if ischar(varargin{j}),
            variable_filter{i}=varargin{j};
            i=i+1;
            flag_filter=1;
        end
    end
end

if isdir('CAA'), % check if is CAA folder, then assume data are there
    dirs = dir('CAA');
else % otherwise assume one is in the CAA data folder
    dirs = dir;
end

nloaded = 0;
for j = 1:numel(dirs)
	if regexp(dirs(j).name,'^C[1-4]_(C|P)P_')
        var_name = dirs(j).name;
        var_name(strfind(var_name,'-'))='_'; % substitute '-' to '_'
        flag_load_variable=1;
        if flag_filter==1, % if there is name filtering required check if to load variable
            for jj=1:length(variable_filter),
                if isempty(strfind(var_name,variable_filter{jj})),
                    flag_load_variable=0;
                end
            end
        end
        if flag_load_variable,
            try
                disp(['caa_load ' var_name]);
                if evalin('caller',['exist(''' var_name ''',''var'')']),
                    disp('Variable exist in memory. NOT LOADING FROM FILE!')
                else
                    evalin('caller',[var_name '=dataobj(''' dirs(j).name filesep '*.cdf'');']);
                end
                nloaded = nloaded + 1;
            catch
                disp(['error loading ' var_name]);
            end
        end
	end
end
if nloaded, fprintf('\nCAA_LOAD : loaded %d variables\n',nloaded);
else disp('CAA_LOAD : nothing to load')
end