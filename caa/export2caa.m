function export2caa(sc_list,varargin)
%export2caa export data into caa CEF files.
% function export2caa(sc_list)
% Input:
%	sc_list - list of sc [optional], default 1:4
%
% We need to export the following data:
% Level 1: wEp[cl_id][12,p34] (mER), P10Hz[cl_id]p[1-4] (mP)
% Level 2: diE[cl_id]p1234 (mEDSI), P[cl_id] (mP)
% Level 3: diEs[cl_id] (mEdB), diVExBs[cl_id] (mEdB), <P[cl_id]> NOT DONE
%
% $Id$
%
% See also exportAscii

% Copyright 2004 Yuri Khotyaintsev (yuri@irfu.se)

if nargin<1, sc_list = 1:4; end

f_name = {'mER', 'mP', 'mEDSI', 'mEdB'};
vars{1} = {'wE?p12', 'wE?p34'};
vars{2} = {'P10Hz?p1', 'P10Hz?p2', 'P10Hz?p3', 'P10Hz?p4', 'P?'};
vars{3} = {'diE?p1234'}; 
vars{4} = {'diEs?', 'diVExBs?'};

for cl_id=sc_list
	for j = 1:length(f_name)
		if exist([f_name{j} '.mat'],'file')
			for k = 1:length(vars{j})
				c_log('load',av_ssub([vars{j}{k} ' <- ' f_name{j}],...
					cl_id))
					c_eval(['load ' f_name{j} ' ' vars{j}{k}],cl_id)
				if exist(av_ssub(vars{j}{k},cl_id),'var')
					c_eval(['exportAscii(' vars{j}{k} ',''mode'',''caa'')'],cl_id)
				else
					c_log('load',av_ssub(['No ' vars{j}{k} ' in ' f_name{j}],...
						cl_id))
				end
			end
		end
	end
end
