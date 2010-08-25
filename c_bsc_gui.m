function c_bsc_gui(varargin)
%C_BSC_GUI STAFF B GUI
%
% This program will load STAFF-SC B data
%
% $Id$

% ----------------------------------------------------------------------------
% "THE BEER-WARE LICENSE" (Revision 42):
% <yuri@irfu.se> wrote this file.  As long as you retain this notice you
% can do whatever you want with this stuff. If we meet some day, and you think
% this stuff is worth it, you can buy me a beer in return.   Yuri Khotyaintsev
% ----------------------------------------------------------------------------

persistent h0;
%persistent inprog_mtx; % this is a kinda mutex to hold "in progress" status

main_fig_id = 28;

if nargin, action = varargin{1};
else action = 'init';
end

if regexp(action,'^click_[X-Z]axes$')
    curr_ax = action(7);
	action = 'click_axes';
end

switch action
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% init
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    case 'init'
        %inprog_mtx = 0;
        
        % Create figure
        h0 = figure(main_fig_id);
        clf
        set(main_fig_id,'Name', 'CLUSTER B-SA GUI')
        
        hnd = guihandles(h0);
        
        hnd.BSCData = [];
        hnd.BSCDataAppend = [];
        hnd.ts_marker = [];
        hnd.cl_id = [];
        
        %Load data
        for cl_id = 1:4
            [ok, hnd.BSCData] = c_load('wBSC?',cl_id);
            if ok && ~isempty(hnd.BSCData), break, end
        end
        
        if isempty(hnd.BSCData), error('no BSC data to load'), end
        hnd.cl_id = cl_id;
        
        [iso_st,dt] = caa_read_interval;
        if isempty(iso_st)
            hnd.tint = [hnd.BSCData(1,1) hnd.BSCData(1,1)];
        else
            hnd.tint = iso2epoch(iso_st) + [0 dt];
        end
        
        % Create Data Axes
        hnd.Xaxes = irf_subplot(3,1,-1);
        hnd.Yaxes = irf_subplot(3,1,-2);
        hnd.Zaxes = irf_subplot(3,1,-3);
        
        guidata(h0,hnd);
        
        c_bsc_gui('replot_all')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% replot_all
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    case 'replot_all'
        % Replot all data
        hnd = guidata(h0);
        h = [hnd.Xaxes hnd.Yaxes hnd.Zaxes];
        
        leg='xyz';
        for comp=1:3
            plot(h(comp),hnd.BSCData(:,1),hnd.BSCData(:,comp+1));
            if ~isempty(hnd.BSCDataAppend)
                hold on
                plot(h(comp),...
                    hnd.BSCDataAppend(:,1),hnd.BSCDataAppend(:,comp+1),'r');
                hold off
            end
            ylabel(h(comp),['B_' leg(comp) ' [nT]'])
        end
        
        % Time span
        irf_zoom(hnd.tint,'x',h);
        
        set(hnd.Xaxes,'Tag','Xaxes',...
            'ButtonDownFcn','c_bsc_gui(''click_Xaxes'')');
        set(hnd.Yaxes,'Tag','Yaxes',...
            'ButtonDownFcn','c_bsc_gui(''click_Yaxes'')');
        set(hnd.Zaxes,'Tag','Zaxes',...
            'ButtonDownFcn','c_bsc_gui(''click_Zaxes'')');
        
        guidata(h0,hnd);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% click_axes
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    case 'click_axes'
        hnd = guidata(h0);
        
        t = get(eval(['hnd.' curr_ax 'axes']),'CurrentPoint');
        t = t(1);
        
        hnd.ts_marker.t = t;
        hnd.ts_marker = replot_t_marker(hnd,hnd.ts_marker);
        
        guidata(h0,hnd);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function hide_t_marker
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function newmarker = hide_t_marker(hnd,marker) %#ok<INUSL>
newmarker.t = [];
if ~isempty(marker)
	newmarker.t = marker.t;
	% Hide the marker if marker.h exists
	if isfield(marker,'h')
		lasterr('')
		try
			for j=3:-1:1, delete(marker.h(j)),end
		catch
			disp('missed marker')
		end
	end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function replot_t_marker
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function newmarker = replot_t_marker(hnd,marker)
h = [hnd.Xaxes hnd.Yaxes hnd.Zaxes];

% Try to hide the marker
newmarker = hide_t_marker(hnd,marker);

for j=1:length(h)
    yy = get(h(j),'YLim');
    hold(h(j),'on')
    newmarker.h(j) = plot(h(j),[newmarker.t newmarker.t],yy,...
        'LineWidth',1,'Color','red');
    hold(h(j),'off')
end