function doSPlot(sp,sc_list)
%doSPlot make EFW summary plots
% doSPlot(sp,[sc_list])
% Input:
% sp - storage path
% sc_list - list of SC [optional]
%
% It is necessary to obtain the Sunward offset and amplitude factor
% before running this routine (use ClusterProc/corrSOffsetM)
%
% Example:
% sc_list = 2:4;
% sp = '/home/yuri/caa-data/20020304';
% %load the necessary data
% getRData(sp,toepoch([2002 03 04 10 00 00]),30*60,sc_list);
% %correct the S-offset and amplitude
% for i=sc_list, corrSOffsetM(ClusterProc(sp),i), end
% %do plots
% doSPlot(sp,sc_list)
%
% $Revision$  $Date$
%
% See also ClusterProc/summaryPlot, ClusterProc/getData

% Copyright 2004 Yuri Khotyaintsev

error(nargchk(1,2,nargin))

if nargin<2, sc_list=1:4; end

cp=ClusterProc(sp);

vars = {'dies','die'};
for cl_id=1:4
	getData(cp,cl_id,'edbs','ang_blank','ang_limit',15);
	getData(cp,cl_id,'edb','ang_blank','ang_limit',15);
	getData(cp,cl_id,'vedbs');

	figure(cl_id)
	summaryPlot(cp,cl_id,'dsi')
end
