function out=irf_pl_tx(varargin)
%IRF_PL_TX   Plot data from all four Cluster spacecraft in the same plot
%
% IRF_PL_TX(x1,x2,x3,x4,[column],[linestyle],[dt1 dt2 dt3 dt4])
% IRF_PL_TX('x?',[column],[linestyle],[dt1 dt2 dt3 dt4])
%	plot variables x1,x2,x3,x4 with time shift dt1...dt4
%	time is 1st column, default plot 2nd column
% IRF_PL_TX(...,'sc_list',sc_list) specify list of spacecraft to plot
% H=IRF_PL_TX(..) return handle to plot
% IRF_PL_TX(AX,...) plot in the specified axis
% IRF_PL_TX(X) where X is structure plots fields X.C1, X.C2, X.C3, X.C4
%
%   column - gives which column to plot. All columns will be plotted
%            in separate panels if set to empty string or ommited.
%   linestyle - string or cell (size 4) in format accepted by plot.
%            Usefull to set line style and marker (but not color).
%   dt1 dt2 dt3 dt4 - timeshifts array
%
%   Example:
%      IRF_PL_TX('irf_abs(B?)')
%      % plot 3 components + magnitude of B1:4.
%      IRF_PL_TX('B?',3:4)
%      % plot 3th and 4th components of B1:4.
%      IRF_PL_TX('B?',3:4,[0 2 3 .5])
%      % plot 3th and 4th components of B1:4 with timeshifts
%      IRF_PL_TX('B?','',[0 2 3 .5])
%      % plot all components of B1:4 with timeshifts
%      IRF_PL_TX('B?','.-')
%      % plot all components of B1:4 using '.-' (line with dot markers)
%      IRF_PL_TX('B?','',[0 2 3 .5],{'.-','*','+','-'})
%      % plot all components of B1:4 using timeshifts and individual
%      % linestyles for each sc
%
% See also IRF_PLOT, PLOT

[ax,args,nargs] = axescheck(varargin{:});
if isempty(ax) % if empty axes
    ax=gca;
end
%hcf=get(ax,'parent'); % get figure handle

if nargs == 0 % show help if no input parameters
    help IRF_PL_TX;
    return
end

% Defaults
X = struct('x1',[],'x2',[],'x3',[],'x4',[],'x5',[]);
sc_list=1:5;	% default plot all s/c data
column = [];	% column to plot 
delta_t = [];	% time shift
line_style = {};% line styles
flagCluster = false;
flagMMS = false;

% Check which are input variables
if ischar(args{1})
    % Variables defined in form 'B?'
    getVariablesFromCaller = true;
    variableNameInCaller   = args{1};
    args(1) = [];
elseif isstruct(args{1}) % format vector.C1, vector.C2,...
  for cc = '1':'4'
    if ~flagMMS % CLuster
      f = ['C',cc];
      if assignField(), flagCluster = true;
      else
        f = ['c',cc];
        if assignField(), flagCluster = true; end
      end
      if flagCluster, continue; end
    end
    %MMS
    f = ['MMS',cc];
    if assignField(), flagMMS = true;
    else
      f = ['mms',cc];
      if assignField(), flagMMS = true; end
    end
  end
  if ~flagCluster && ~flagMMS % THEMIS
    for cc = '1':'5'
        f = ['C',char(cc+15)];
        if ~assignField(), f = lower(f); assignField(), end
    end
  end
  args(1) = [];
  getVariablesFromCaller = false;
else
    % Variables given as 4 input paramters
    if length(args)<4, error('use IRF_PL_TX(x1,x2,x3,x4) or IRF_PL_TX(''x?'')'), end
    X.x1 = args{1}; X.x2 = args{2}; X.x3 = args{3}; X.x4 = args{4}; % assign x1,x2..x4
    args = args(5:end);
    getVariablesFromCaller = false;
end

% Check if column to plot is specified as input
if ~isempty(args)
    if isnumeric(args{1})
        column = args{1};
        args(1) = [];
    elseif isempty(args{1}) % empty string means default matrix size
        args(1) = [];
    end
end

% check for other input parameters
while ~isempty(args)
    if ischar(args{1})
        if strcmp(args{1},'sc_list')
            args(1) = [];
            sc_list=args{1};
			if isempty(sc_list)
				irf_log('fcal','sc_list empty');
				return;
			end
        else
            % assume that argument defines Linestyle
            if isempty(line_style), c_eval('line_style(?)={args{1}};')
            else, irf_log('fcal','L_STYLE is already set')
            end
        end
    elseif iscell(args{1}) && length(args{1})==4
        % Individual linestyles for each sc
        if isempty(line_style), line_style = args{1};
        else, irf_log('fcal','L_STYLE is already set')
        end
    elseif iscell(args{1})
        % Individual linestyles for each sc
        irf_log('fcal','L_STYLE must be a cell with 4 elements')
    elseif isnumeric(args{1}) && length(args{1})==4
        % dt1..dt4
        if isempty(delta_t), delta_t = args{1};
        else, irf_log('fcal','DELTA_T is already set')
        end
    else
        irf_log('fcal',['ignoring input argument: ' args{1}])
    end
	args(1) = [];
end
if isempty(delta_t), delta_t = [0 0 0 0 0]; end

% Get variable values from caller if needed
if getVariablesFromCaller
    for cl_id=sc_list
        ttt = evalin('caller',irf_ssub(variableNameInCaller,cl_id),'[]');
        X.(sprintf('x%d',cl_id)) = ttt; clear ttt
    end
end

% If column empty, check which columns to plot
if isempty(column)
	for cl_id=sc_list
    if isa(X.(fId),'TSeries'), nCol = size(X.(fId).data,2);
    else, nCol = size(X.(fId),2) - 1;
    end
		if ~isempty(nCol) && nCol > 0, column = 1:nCol; break, end
	end
end
if isempty(column)
    irf_log('fcal','all inputs are empty')
    return
end

% check which spacecraft data are available
sc_list_with_data=[];
for cl_id=1:5
  if ~isempty(X.(fId)), sc_list_with_data=[sc_list_with_data cl_id]; end %#ok<AGROW>
end

% if more than one column reset figure 
if length(column) > 1 && numel(ax) ~= numel(column)
	ax = irf_plot(length(column),'reset'); 
end

% define Cluster colors
cluster_colors={[0 0 0]; [1 0 0];[0 0.5 0];[0 0 1];[0 1 1]};

for j=1:length(column)
  for cl_id=sc_list_with_data
    if isempty(X.(fId)), continue, end
    if isa(X.(fId),'TSeries')
      dataTmp = [X.(fId).time.epochUnix-delta_t(cl_id) ...
        double(X.(fId).data(:,column(j)))];
    else, dataTmp = [X.(fId)(:,1)-delta_t(cl_id), X.(fId)(:,column(j)+1)];
    end
    if isempty(line_style)
      hl = irf_plot(ax(j), dataTmp, 'color', cluster_colors{cl_id});
    else
      hl = irf_plot(ax(j), dataTmp, line_style{cl_id}, ...
        'color', cluster_colors{cl_id});
    end
    hold(ax(j),'on');
    set(hl,'Tag','C?');
  end
  hold(ax(j),'off');
  irf_zoom(ax(j),'y'); % optimize Y zoom to skip labels at top and bottom
  grid(ax(j),'on');
end
irf_timeaxis(ax);
irf_figmenu;


if nargout > 0, out = ax; end

  function res = assignField()
    res = false;
    if isfield(args{1},f)
      X.(['x' cc']) =  args{1}.(f);
      res = true; return
    end
  end
  
  function res = fId
   res = sprintf('x%d',cl_id);
  end
end
