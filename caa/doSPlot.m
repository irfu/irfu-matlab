function doSPlot(sp,sc_list,options)
%doSPlot make EFW summary plots
% doSPlot([sp],[sc_list],[options])
% Input:
%   sp - storage path. [default '.']
%   sc_list - list of SC [default 1:4]
%   options :
%      'noproc' - do not recalculate anything
%
% It is necessary to obtain the Sunward offset and amplitude factor
% before running this routine (use ClusterProc/corrSOffsetM)
%
% Examples:
%   doSPlot
%   % make plots using data in the current directory,
%   % for SC 1 to 4, with reprocessing of data.
%
%   doSPlot('noproc')
%   % same as above, but without reprocessing
%
%   doSPlot([2 4],'noproc')
%   % make plots for SC 2 and 4 without reprocessing
%
%   sc_list = 2:4;
%   sp = '/home/yuri/caa-data/20020304';
%   %load the necessary data
%   getRData(sp,toepoch([2002 03 04 10 00 00]),30*60,sc_list);
%   %correct the S-offset and amplitude
%   for i=sc_list, corrSOffsetM(ClusterProc(sp),i), end
%   %do plots
%   doSPlot(sp,sc_list)
%
% $Revision$  $Date$
%
% See also ClusterProc/summaryPlot, ClusterProc/getData

% Copyright 2004 Yuri Khotyaintsev

do_proc = 1;
error(nargchk(0,2,nargin))
switch nargin
case 0
	sp = '.';
	sc_list = 1:4;
case 1
	if isstr(sp), 
		if strcmp(sp,'noproc')
			sp = '.';
			do_proc = 0; 
		end
		sc_list = 1:4;
	elseif isnumeric(sp)
		sc_list = sp; 
		sp = '.';
	else, error('Input must be eather a storage path or sc_list (1:4)')
	end
case 2
	if isnumeric(sp)
		if strcmp(sc_list,'noproc'), do_proc = 0; end
		sc_list = sp; 
		sp = '.';
	end
case 3
	if strcmp(options,'noproc'), do_proc = 0; end
end

cp=ClusterProc(sp);

vars = {'dies','die'};
for cl_id=sc_list
	if do_proc
		getData(cp,cl_id,'edbs','ang_blank','ang_limit',15);
		getData(cp,cl_id,'edb','ang_blank','ang_limit',15);
		getData(cp,cl_id,'vedbs');
	end

	figure(cl_id)
	summaryPlot(cp,cl_id,'dsi')
end
