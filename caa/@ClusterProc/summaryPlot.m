function summaryPlot(cp,cl_id,cs)
% summaryPlot make a summary plot
% summaryPlot(cp,cl_id,cs) make a summary plot
% cs is a coordinate system : 'dsi' [default] of 'gse'
%
% $Revision$  $Date$

% Copyright 2004 Yuri Khotyaintsev
error(nargchk(2,3,nargin))

if nargin < 3, cs = 'dsi'; %DSI
end

if ~strcmp(cs,'dsi') & ~strcmp(cs,'gse')
   disp('unknown CS. defaulting to DSI')
	cs= 'dsi';
end

% load data
if strcmp(cs,'dsi') 
	q_list = {'P?','diE?','diVs?'};
	l_list = {'SC pot [-V]','E DSI [mV/m]','V DSI [km/s]'};
else
	q_list = {'P?','E?','Vs?'};
	l_list = {'SC pot [-V]','E GSE [mV/m]','V GSE [km/s]'};
end
f_list = {'mP','mEdB','mEdB'};

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
			if k==2 % E-field
				eval(av_ssub(['data{n_plots}=' q_list{k} '(:,1:4);'],cl_id)) 
				labels{n_plots} = l_list{k};
				n_plots = n_plots + 1;
				eval(av_ssub(['data{n_plots}=' q_list{k} '(:,[1 5]);'],cl_id)) 
				labels{n_plots} = '\alpha(B,spin) [deg]';
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

%Plotting
clf

for k=1:n_plots
	h{k} = subplot(n_plots,1,k);
	av_tplot(data{k});
	ylabel(labels{k})
	if k==1, title(['Cluster ' num2str(cl_id,'%1d')]), end
end

addPlotInfo

for k=n_plots:-1:1
	axes(h{k})
	if min(size(data{k}))>2, legend('X','Y','Z'), end
end
