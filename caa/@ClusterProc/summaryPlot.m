function out=summaryPlot(cp,cl_id,cs,st,dt)
% summaryPlot make EFW summary plot
% h = summaryPlot(cp,cl_id,[cs],[st,dt])
% Input:
% cp - ClusterProc object
% cl_id - SC#
% cs is a coordinate system : 'dsi' [default] of 'gse'
% st, dt - start time and interval length [optional]
% 
% Output:
% h - axes handles // can be omitted
%
% Example:
% summaryPlot(ClusterProc('/home/yuri/caa-data/20020304'),1,'gse')
%
% $Id$

% Copyright 2004 Yuri Khotyaintsev
error(nargchk(2,5,nargin))

if nargin<3, cs = 'dsi'; %DSI
end

if ~strcmp(cs,'dsi') & ~strcmp(cs,'gse')
   disp('unknown CS. defaulting to DSI')
	cs= 'dsi';
end

% load data
if strcmp(cs,'dsi') 
	q_list = {'P?','diBPP?','diE?','diEs?','diVExBs?'};
	l_list = {'SC pot [-V]','B PP DSI [nT]','E DSI [mV/m]','E DSI [mV/m]','V=ExB DSI [km/s]'};
else
	q_list = {'P?','BPP?','E?','Es?','VExBs?'};
	l_list = {'SC pot [-V]','B PP GSE [nT]','E GSE [mV/m]','E GSE [mV/m]','V=ExB GSE [km/s]'};
end
f_list = {'mP','mBPP','mEdB','mEdB','mEdB'};

old_pwd = pwd;
cd(cp.sp) %enter the storage directory

n_plots = 0;
data = {};
labels = {};
for k=1:length(q_list)
	if exist(['./' f_list{k} '.mat'],'file')
		eval(av_ssub(['load ' f_list{k} ' ' q_list{k}],cl_id))
		if exist(av_ssub(q_list{k},cl_id))
			n_plots = n_plots + 1;
			if k==2 % B-field
				c_eval(['data{n_plots}=av_abs(' q_list{k} '(:,1:4));'],cl_id)
				labels{n_plots} = l_list{k};
			elseif k==3 % E-field
				eval(av_ssub(['data{n_plots}=' q_list{k} '(:,1:4);'],cl_id)) 
				labels{n_plots} = l_list{k};
				n_plots = n_plots + 1;
				eval(av_ssub(['data{n_plots}=' q_list{k} '(:,[1 5]);'],cl_id)) 
				labels{n_plots} = '\theta (B,spin) [deg]';
			else
				eval(av_ssub(['d_t=' q_list{k} ';'],cl_id))
				labels{n_plots} = l_list{k};
				if min(size(d_t))> 4
					data{n_plots} = d_t(:,1:4);
				else	
					data{n_plots} = d_t;
				end
				clear d_t
			end
		end
	end
end

cd(old_pwd)

if n_plots==0, return, end %nothing to plot

% define time limits
if nargin<4,
	t_st = 1e32;
	t_end = 0;
	for k=1:n_plots
		t_st = min(t_st,data{k}(1,1));
		t_end = max(t_end,data{k}(end,1));
	end
else
	t_st = st;
	t_end = st + dt;
end

%Plotting
clf
orient tall

dummy = 'data{1}';
for k=2:n_plots, dummy = [dummy ',data{1}'];, end
eval(['h = av_tplot({' dummy '});']) 
clear dummy

for k=1:n_plots
	axes(h(k))
	av_tplot(data{k});
	av_zoom([t_st t_end],'x',h(k))
	set(gca,'YLim',get(gca,'YLim')*.99)
	ylabel(labels{k})
	if k==1, title(['EFW, Cluster ' num2str(cl_id,'%1d')]), end
	if k<n_plots, xlabel(''),set(gca,'XTickLabel',[]), end		
end

addPlotInfo

for k=n_plots:-1:1
	if min(size(data{k}))>2, legend(h(k),'X','Y','Z','Location','NorthEastOutside'), end
end

if nargout>0, out=h;,end
