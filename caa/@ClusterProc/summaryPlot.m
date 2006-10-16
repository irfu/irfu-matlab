function out=summaryPlot(cp,cl_id,varargin)
% summaryPlot make EFW summary plot
%
% h = summaryPlot(cp,cl_id,[options])
%
% Input:
%   cp - ClusterProc object
%   cl_id - SC#
%   Options: go in pair 'option', value
%    'cs' - coordinate system : 'dsi' [default] of 'gse'
%    'st', 'dt' - start time (ISDAT epoch) and interval length (sec)
%    'fullb' - use full resolution B FGM
%    'leavewhip' - plot time intervals with Whisper pulses
% 
% Output:
%   h - axes handles // can be omitted
%
% Example:
%   summaryPlot(ClusterProc('/home/yuri/caa-data/20020304'),1,'cs','gse')
%   summaryPlot(ClusterProc('.'),2,'st',toepoch([2004 1 4 12 47 0]),'dt',60,'fullb')
%
% $Id$

% Copyright 2004 Yuri Khotyaintsev
error(nargchk(2,9,nargin))

if nargin>2, have_options = 1; args = varargin;
else have_options = 0;
end

% Default values
cs = 'dsi';
st = 0;
dt = 0;
have_tint = 0;
use_fullb = 'rs';
flag_rmwhip = 1; 

while have_options
	l = 2;
	if length(args)>=1
		switch(args{1})
		case 'cs'
			if ischar(args{2}), 
				cs = args{2};
				if ~strcmp(cs,'dsi') && ~strcmp(cs,'gse')
					irf_log('fcal','unknown CS. defaulting to DSI')
					cs = 'dsi';
				end
            else irf_log('fcal','wrongArgType : CS must be a string')
			end
		case 'st'
			if isnumeric(args{2}), st = args{2};
            else irf_log('fcal','wrongArgType : ST must be numeric')
			end
		case 'dt'
			if isnumeric(args{2}), dt = args{2};
            else irf_log('fcal','wrongArgType : DT must be numeric')
			end
		case 'fullb'
			use_fullb = '';	l = 1;
		case 'leavewhip'
			flag_rmwhip = 0;
		otherwise
        	irf_log('fcal',['Option ''' args{1} '''not recognized'])
    	end
		if length(args) > l, args = args(l+1:end);
		else break
		end
	else
		error('caa:wrongArgType','use summaryPlot(..,''option'',''value'')')
	end
end

if st && dt, have_tint = 1; end

% Define variables we want to plot
if strcmp(cs,'dsi') 
	q_list = {'P?',['diB' use_fullb '?'],'diE?p1234','diEs?','diVExBs?'};
	l_list = {'SC pot [-V]','B DSI [nT]','E DSI [mV/m]','E DSI [mV/m]','V=ExB DSI [km/s]'};
else
	q_list = {'P?',['B' use_fullb '?'],'E?','Es?','VExBs?'};
	l_list = {'SC pot [-V]','B GSE [nT]','E GSE [mV/m]','E GSE [mV/m]','V=ExB GSE [km/s]'};
end

old_pwd = pwd;
cd(cp.sp)

n_plots = 0;
data = {};
labels = {};

if flag_rmwhip
	if exist('./mFDM.mat','file')
		c_eval('load mFDM WHIP?',cl_id)
	end
	
end

% Load data
for k=1:length(q_list)
	if c_load(q_list{k},cl_id)
		if have_tint
			c_eval([q_list{k} '=irf_tlim(' q_list{k} ',st+[0 dt]);'],cl_id)
		end
		
		if flag_rmwhip && (k==1 || k==3)
			if exist(irf_ssub('WHIP?',cl_id),'var')
				irf_log('proc','not using times with Whisper pulses')
				c_eval([q_list{k} '=caa_rm_blankt(' q_list{k}  ',WHIP? );'],cl_id)
			end
		end
		
		n_plots = n_plots + 1;
		if k==2 % B-field
			c_eval(['data{n_plots}=irf_abs(' q_list{k} '(:,1:4));'],cl_id)
			labels{n_plots} = l_list{k};
		elseif k==3 % E-field
			if strcmp(cs,'dsi') 
				% correct DSI offsets
				dsiof = c_ctl(cl_id,'dsiof');
				if isempty(dsiof), dsiof = [1+0i 1]; end
				[ok,Dxy] = c_load('Ddsi?',cl_id);
				if ~ok, Dxy = dsiof(1); end
				[ok,Da] = c_load('Damp?',cl_id);
				if ~ok, Da = dsiof(2); end
				c_eval([q_list{k} '=caa_corof_dsi(' q_list{k} ',Dxy,Da);'],cl_id)
				clear dsiof Dxy
			end

			c_eval(['data{n_plots}=' q_list{k} '(:,1:3);'],cl_id)
			labels{n_plots} = l_list{k};
		elseif k==4 % Es-field
			c_eval(['d_t=' q_list{k} ';'],cl_id)
			data{n_plots} = d_t(:,[1 5]); 
			labels{n_plots} = '\theta (B,spin) [deg]';
			n_plots = n_plots + 1;
			data{n_plots} = d_t(:,1:4);
			labels{n_plots} = l_list{k};
			clear d_t
		else
			c_eval(['d_t=' q_list{k} ';'],cl_id)
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

cd(old_pwd)

if n_plots==0, return, end % Nothing to plot

% Define time limits
if have_tint
	t_st = st;
	t_end = st + dt;
else
	t_st = 1e32;
	t_end = 0;
	for k=1:n_plots
		t_st = min(t_st,data{k}(1,1));
		t_end = max(t_end,data{k}(end,1));
	end
end

% Plotting
clf
orient tall

dummy = 'data{1}';
for k=2:n_plots, dummy = [dummy ',data{1}']; end
eval(['h = irf_plot({' dummy '});']) 
clear dummy

for k=1:n_plots
	axes(h(k))
	irf_plot(data{k});
	irf_zoom([t_st t_end],'x',h(k))
	set(gca,'YLim',get(gca,'YLim')*.99)
	ylabel(labels{k})
	if k==1, title(['EFW, Cluster ' num2str(cl_id,'%1d')]), end
	if k<n_plots, xlabel(''),set(gca,'XTickLabel',[])
    else
        add_timeaxis;
		% This magic is needed for correct location of panels on printouts
		ttt = get(gca,'XTickLabel'); 
		ttt(end) = {' '}; 
		set(gca,'XTickLabel',ttt)
	end
	if k==n_plots
		YLIM = 990;
		yl = get(gca,'YLim');
		if yl(1)<-YLIM, yl(1) = -YLIM; end
		if yl(2)>YLIM, yl(2) = YLIM; end
		set(gca,'YLim',yl);
		clear yl YLIM
	end
end

irf_pl_add_info

lyy = 0;
for k=n_plots:-1:1
	if min(size(data{k}))>2
        legend(h(k),'X','Y','Z','Location','NorthEastOutside')
        if lyy==0, pos = get(h(k),'Position'); lyy = pos(3); clear pos, end
    end
end

if lyy
    for k=n_plots:-1:1
        pos = get(h(k),'Position'); 
        set(h(k),'Position', [pos(1) pos(2) lyy pos(4)])
    end
end

if nargout>0, out=h; end
