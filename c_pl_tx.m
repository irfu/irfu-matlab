function out=c_pl_tx(varargin)
%C_PL_TX   Plot data from all four Cluster spacecraft in the same plot
%
% c=c_pl_tx('x?',column,[linestyle],[dt1 dt2 dt3 dt4])
% c=c_pl_tx(x1,x2,x3,x4,[column],[linestyle],[dt1 dt2 dt3 dt4])
%   plot all 4 cluster values, time is in the first column
%
%   column - gives which column to plot. All columns will be plotted 
%            if set to empty string or ommited.
%   linestyle - string or cell (size 4) in format accepted by plot.
%            Usefull to set line style and marker (but not color).
%   dt1 dt2 dt3 dt4 - timeshifts array
%
%   Example:
%      c_pl_tx('irf_abs(B?)') 
%      % plot 3 components + magnitude of B1:4.
%      c_pl_tx('B?',3:4)
%      % plot 3th and 4th components of B1:4.
%      c_pl_tx('B?',3:4,[0 2 3 .5])
%      % plot 3th and 4th components of B1:4 with timeshifts
%      c_pl_tx('B?','',[0 2 3 .5])
%      % plot all components of B1:4 with timeshifts
%      c_pl_tx('B?','.-')
%      % plot all components of B1:4 using '.-' (line with dot markers)
%      c_pl_tx('B?','',[0 2 3 .5],{'.-','*','+','-'})
%      % plot all components of B1:4 using timeshifts and individual 
%      % linestyles for each sc
%
% See also IRF_PLOT, PLOT
%
% $Id$

error(nargchk(1,8,nargin))

args = varargin;

if isstr(args{1})
% We have 4 arguments
	for cl_id=1:4
		ttt = evalin('caller',irf_ssub(args{1},cl_id),'[]'); 
		eval(irf_ssub('x? =ttt;',cl_id)); clear ttt
	end
	if length(args) > 1, args = args(2:end); 
	else, args = ''; end
else
	if length(args)<4, error('use c_pl_tx(x1,x2,x3,x4) or c_pl_tx(''x?'')'), end
	% We have x1,x2..x4
	c_eval('x? = args{?};');
	if length(args) > 4, args = args(5:end); 
	else, args = ''; end
end

% Check for deprecated 1
if length(args)>1
	if length(args{2})==1 & args{2}==1
		irf_log('fcal','this usage of c_pl_tx is deprecated. see help c_pl_tx')
		args(2) = [];
	end
end

column = [];
if length(args)>0
	if isnumeric(args{1})
		column = args{1};
		args = args(2:end);
	elseif isstr(args{1})
		% empty string means default matrix size
		if isempty(args{1}), args = args(2:end); end
	end
end
if isempty(column)
	% try to guess the size of the matrix
	column = size(x1,2);
	if column > 2, column = 2:column; end
end

delta_t = [];
l_style = {};

if length(args)>0, have_args = 1;
else, have_args = 0;
end

while have_args
	if isstr(args{1})
		%linestyle
		if isempty(l_style), c_eval('l_style(?)={args{1}};')
		else, irf_log('fcal','L_STYLE is already set')
		end
		args = args(2:end);
	elseif iscell(args{1}) & length(args{1})==4
		%individual linestyles for each sc
		if isempty(l_style), l_style = args{1};
		else, irf_log('fcal','L_STYLE is already set')
		end
		args = args(2:end);
	elseif isnumeric(args{1}) & length(args{1})==4
		%dt
		if isempty(delta_t), delta_t = args{1};
		else, irf_log('fcal','DELTA_T is already set')
		end
		args = args(2:end);
	else
		irf_log('fcal','ignoring input argument')
		args = args(2:end);
	end
	
	if length(args)>0, have_args = 1;
	else, have_args = 0;
	end
end

if isempty(delta_t), delta_t = [0 0 0 0]; end

if isempty(l_style), l_style= {'k','r','g','b'};
else
	cls='krgb';
	c_eval('l_style(?)={[cls(?) l_style{?}]};')
	clear cls
end

% t_start_epoch is saved in figures user_data variable
% check first if it exist otherwise assume zero
ud=get(gcf,'userdata');
if isfield(ud,'t_start_epoch'), 
	t_start_epoch=ud.t_start_epoch;
elseif x1(1,1)> 1e8 | x2(1,1)> 1e8 | x3(1,1)> 1e8 | x4(1,1)> 1e8, 
	% set start_epoch if time is in isdat epoch, warn about changing t_start_epoch
	t_start_epoch=min([x1(1,1) x2(1,1) x3(1,1) x4(1,1)]);
	ud.t_start_epoch=t_start_epoch;set(gcf,'userdata',ud);
	irf_log('proc',['user_data.t_start_epoch is set to ' epoch2iso(t_start_epoch)]);
else
	t_start_epoch=0;
end

c_eval('ts?=t_start_epoch+delta_t(?);')

if length(column) == 1,
	pl = '';
	for jj=1:4
		if eval(irf_ssub('~isempty(x?)',jj))
			c_eval(['s_s=''(x?(:,1)-ts?),x?(:,column),''''' l_style{jj} ''''''';'],jj);
			if isempty(pl), pl = s_s; else, pl = [pl ',' s_s]; end
			clear s_s
		end
	end 
	eval(['h=plot(' pl ');'])
	c=get(h(1),'parent');
	grid on
else
	clf;
	for j=1:length(column),
		c(j)=irf_subplot(length(column),1,-j);
		pl = '';
		for jj=1:4
			if eval(irf_ssub('~isempty(x?)',jj))
				c_eval(['s_s=''(x?(:,1)-ts?),x?(:,column(j)),''''' l_style{jj} ''''''';'],jj);
				if isempty(pl), pl = s_s; else, pl = [pl ',' s_s]; end
				clear s_s
			end
		end 
		eval(['plot(' pl ')'])
		grid on
        ud=get(gcf,'userdata');ud.subplot_handles=c;set(gcf,'userdata',ud);
	end
end

add_timeaxis(c);
irf_figmenu;

if nargout > 0, out = c; end

