function res = plot(dobj,varargin)
%PLOT([H], dobj, var, [comp], [options])  plot a variable 'var' in dataobj 'dobj'
%   dobj - data object, see also CAA_LOAD
%   var  - can be a string of variable name in data object
%        - or can be a a variable itself (structure, see caa form in C_CAA_VAR_GET)
%   comp - components to plot
%
% OPTIONS - one of the following:
%	'AX'         - axis handles to use
%   'COMP'       - components to plot
%   'SUM_DIM1'   - average over first dimension (frequency, azimuthal angle, energy)
%   'SUM_DIM1PITCH'   - average over first dimension (pitch angle)
%   'COMP_DIM1'  - form subplots from that component
%					similar 'SUM_DIM2','SUM_DIM3','SUM_DIM2PITCH',....
%   'nolabels'   - only plot data, do not add any labels or text
%   'NoColorbar' - do not plot colorbar
%   'ColorbarLabel' - specify colorbar label
%   'FitColorbarLabel' - fit the text of colorbar label to colobar height
%   'FillSpectrogramGaps' - fill data gaps with previous value (makes spectrograms look nice on screen)
%   'ClusterColors' - use Cluster colors C1-black, C2-red, C3-green, C4-blue
%
% See also
%		CAA_META
%
% for common cluster variables see: https://bit.ly/pKWVKh

% ----------------------------------------------------------------------------
% SPDX-License-Identifier: Beerware
% "THE BEER-WARE LICENSE" (Revision 42):
% <yuri@irfu.se> wrote this file.  As long as you retain this notice you
% can do whatever you want with this stuff. If we meet some day, and you think
% this stuff is worth it, you can buy me a beer in return.   Yuri Khotyaintsev
% ----------------------------------------------------------------------------

%narginchk(2,14)

[ax,args,~] = axescheck(varargin{:});
if isempty(ax)
  ax=gca;
  create_axes = 1;
else
  create_axes = 0;
end
var_s=args{1};
args=args(2:end);

LCOMP = 3;

if isempty(dobj.data)
  irf_log('fcal','Dataobject empty, nothing to plot');
  return
end
if ischar(var_s)             % input is the name of variable
  data = getv(dobj,var_s);
else                         % input is variable itself
  data=var_s;
  var_s=data.name;
end

if isempty(data), error('VAR_S not found'), end
if isempty(data.data), error('Nothing to plot (Empty dataset)'), end
dim = length(data.variance(3:end));
dep = getdep(dobj,var_s);
units = corr_latex(getunits(dobj,var_s));
fieldnam = findva(dobj,'FIELDNAM',var_s);
ii = regexp(fieldnam,'_'); fieldnam(ii) = ' '; % Get rid of underscores
lablaxis = getlablaxis(dobj,var_s);
ii = regexp(lablaxis,'_'); lablaxis(ii) = ' '; % Get rid of underscores
cs = getcs(dobj,var_s);
if ~isempty(cs), cs = [' ' cs]; end
fillv = getfillval(dobj,var_s);
data.data(data.data==fillv) = NaN;


%% INPUT ARGUMENTS

% Default values that can be override by options
sum_dim = 0;  % along which dimension to sum
comp_dim = 0; % component dimension used to separate subplots
use_comp = 0; % pick up separate component values in separate plots
comp = [];    % index of component vector values to pick up
ydim = 0;     % default dimension of data used for y axis (0- data value itself);
plot_properties=cell(0);
flag_lineplot = 0;
flag_spectrogram = 0;
flag_labels_is_on=1;
flag_colorbar_is_on=1;
flag_colorbar_label_is_manually_specified=0;
flag_colorbar_label_fit_to_colorbar_height_is_on=0;
flag_fill_spectrogram_gaps=0;
line_color=''; % default line color; can be changed with flags, e.g. clustercolors
flag_use_cluster_colors=0;
flag_log = 1;
flagPitchAngleAverage = false; % default averages are simple means, pitch angle average
%			has to take into account the size of spherical angle for each pitch angle

arg_pos = 0;
while ~isempty(args)
  arg_pos = arg_pos + 1;
  l = 1;
  if arg_pos==1 && isnumeric(args{1})
    use_comp = 1;
    comp = args{1};
  else
    switch(lower(args{1}))
      case 'ax'
        l = 2;
        if all(ishandle(args{2}))
          ax = args{2};
          create_axes = 0;
        else, disp('invalid value for AX')
        end
      case 'comp'
        l = 2;
        if isnumeric(args{2})
          use_comp = 1;
          comp = args{2};
        else
          disp('invalid value for COMP')
        end
      case 'clustercolors'
        flag_use_cluster_colors = 1;
        if strfind(var_s,'C1'), line_color='k'; %#ok<STRIFCND>
        elseif strfind(var_s,'C2'), line_color='r'; %#ok<STRIFCND>
        elseif strfind(var_s,'C3'), line_color='g'; %#ok<STRIFCND>
        elseif strfind(var_s,'C4'), line_color='b'; %#ok<STRIFCND>
        else, flag_use_cluster_colors = 0;
        end
      case 'nolabels'
        flag_labels_is_on = 0;
      case 'nocolorbar'
        flag_colorbar_is_on = 0;
      case 'colorbarlabel'
        l=2;
        if ischar(args{2}) || iscell(args{2})
          flag_colorbar_label_is_manually_specified=1;
          colorbar_label = args{2};
        else
          disp('invalid value for ColorbarLabel in PLOT')
        end
      case 'fitcolorbarlabel'
        flag_colorbar_label_fit_to_colorbar_height_is_on=1;
      case 'fillspectrogramgaps'
        flag_fill_spectrogram_gaps=1;
      case 'lin'
        flag_log = 0; % for spectrograms
      case 'sum_dim1'
        sum_dim = 1;
      case 'sum_dim2'
        sum_dim = 2;
      case 'sum_dim3'
        sum_dim = 3;
      case 'sum_dim1pitch'
        sum_dim = 1;flagPitchAngleAverage = true;
      case 'sum_dim2pitch'
        sum_dim = 2;flagPitchAngleAverage = true;
      case 'sum_dim3pitch'
        sum_dim = 3;flagPitchAngleAverage = true;
      case 'comp_dim1'
        comp_dim = 1;
      case 'comp_dim2'
        comp_dim = 2;
      case 'comp_dim3'
        comp_dim = 3;
      otherwise
        disp('unknown argument')
        disp('the rest or arguments are passed to plot routines');
        plot_properties=args;
        break
    end
  end
  args = args(l+1:end);
end

if comp_dim ~=0 && comp_dim == sum_dim
  error('SUM_DIM and COMP_DIM must be different')
end

%% DATA PROCESSING


% define summing dimension and component to plot when component not defined
if comp_dim == 0
  if dim == 1, comp_dim = 1;
  elseif dim == 2 && sum_dim == 0,    comp_dim = 2; ydim = 1;
  elseif dim == 2 && sum_dim == 1,    comp_dim = 2; ydim = 0;
  elseif dim == 2 && sum_dim == 2,    comp_dim = 1;
  elseif dim == 3 && sum_dim == 0,    comp_dim = 3; sum_dim = 1;
  elseif dim == 3 && sum_dim == 1,    comp_dim = 3;
  elseif dim == 3 && sum_dim == 2,    comp_dim = 3;
  elseif dim == 3 && sum_dim == 3,    comp_dim = 2; ydim = 1;
  end
end

if dim ==3
  switch comp_dim
    case 1
      if sum_dim == 2, ydim = 3; else, ydim = 2; end
    case 2
      if sum_dim == 1, ydim = 3; else, ydim = 1; end
    case 3
      if sum_dim == 1, ydim = 2; else, ydim = 1; end
  end
end


if dim == 0
  plot_data = {double(data.data)};
  flag_lineplot = 1;
  
elseif dim == 1
  if use_comp
    plot_data = cell(size(comp));
    for i=1:length(comp)
      plot_data{i} = double(data.data(:,comp(i)));
    end
    flag_lineplot = 1;
  else
    plot_data = {double(data.data)};
    if dim == 1
      if isfield(data,'TENSOR_ORDER')
        if isnumeric(data.TENSOR_ORDER)
          tensorOrder = data.TENSOR_ORDER;
        elseif ischar(data.TENSOR_ORDER)
          tensorOrder = str2double(data.TENSOR_ORDER); % TODO TENSOR_ORDER should be defined numeric
        else
          irf.log('critical','data.TENSOR_ORDER of unknown type!');
          error('dataobj/plot: data.TENSOR_ORDER of unknown type!');
        end
        if data.dim(tensorOrder+1) > 1
          flag_spectrogram = 1;
          ydim = 1;
        else
          flag_lineplot = 1;
        end
      elseif isfield(data,'DEPEND_1')
        flag_spectrogram = 1;
        ydim = 1;
      else
        flag_lineplot = 1;
      end
    end
  end
  
elseif dim == 2
  if sum_dim > 0
    if flagPitchAngleAverage
      irf.log('notice','Using pitch angle average!');
      pitchAngles=getfield(get(dobj,dep.DEPEND_X{sum_dim}),'data');
      data.data = irf.pitch_angle_average(double(data.data),...
        pitchAngles,[],[],sum_dim+1);
    else
      data.data = irf.nanmean(double(data.data),sum_dim+1);
    end
    if sum_dim == 1, comp_dim = 2;
    else, comp_dim = 1;
    end
  end
  
  if sum_dim > 0 && isfield(data,'DEPEND_1') && ~use_comp
    plot_data = {squeeze(data.data)};
    ydim = comp_dim;
  else
    ndim = data.dim(comp_dim);
    if ~use_comp, comp=1:ndim; end
    plot_data = cell(size(comp));
    for i=1:length(comp)
      switch comp_dim
        case 1
          plot_data{i} = squeeze(data.data(:,comp(i),:));
        case 2
          plot_data{i} = squeeze(data.data(:,:,comp(i)));
        otherwise
          error('smth wrong')
      end
    end
    if comp_dim == 2, ydim =1; else, ydim = 2; end
  end
  
  plot_f = 2;
  if ~isfield(data,'DEPEND_1'), plot_f = plot_f-1; end
  if sum_dim > 0, plot_f = plot_f-1; end
  if use_comp, plot_f = plot_f-1; end
  if plot_f > 0
    flag_spectrogram = 1;
  else
    flag_lineplot = 1;
  end
  
elseif dim == 3
  if use_comp == 0 && sum_dim == 0
    if comp_dim==2
      sum_dim = 1;
    else
      sum_dim = 2;
    end
  end
  
  % ignore NaNs when averaging
  data.data = irf.nanmean(double(data.data),sum_dim+1);
  
  ndim = data.dim(comp_dim);
  if ~use_comp, comp=1:ndim; end
  if ndim == 1
    plot_data = {data.data};
  else
    plot_data = cell(size(comp));
    for i=1:length(comp)
      switch comp_dim
        case 1
          plot_data{i} = squeeze(data.data(:,comp(i),:,:));
        case 2
          plot_data{i} = squeeze(data.data(:,:,comp(i),:));
        case 3
          plot_data{i} = squeeze(data.data(:,:,:,comp(i)));
        otherwise
          error('smth wrong')
      end
    end
  end
  flag_spectrogram = 1;
else
  error('plotting not implememnted')
end


if flag_lineplot
  %% PLOTTING -- LINE PLOT
  if isfield(dep,'DEPEND_O')
    if strcmpi(dep.DEPEND_O.type,'tt2000')
      timeLine = EpochTT(dep.DEPEND_O.data).epochUnix;
    else, timeLine = dep.DEPEND_O.data;
    end
    h = irf_plot(ax,[timeLine plot_data{:}],line_color,plot_properties{:});
  else
    h = plot(ax,data.data,line_color,plot_properties{:});
  end
  flab = getlablaxis(dobj,var_s);
  lab_1 = '';
  if ~isempty(dep.DEPEND_X) %
    dep_x_s = dep.DEPEND_X{comp_dim,1};
    dep_x = getv(dobj,dep_x_s);
    if ~isempty(dep_x)
      if strcmp(dep_x.type,'char')
        if use_comp % pick up components, data should be char
          lab_1 = ['(' dep_x.data(1,:,comp) ')'];
        else % data are values, label under LABLAXIS
          legend(ax,num2cell(dep_x.data(1,:,:),2), 'Location','NorthWest')
          legend(ax,'boxoff')
        end
      else
        if ~isempty(comp) && isfield(dep_x,'UNITS')
          lab_1 = ['(' num2str(dep_x.data(1,comp),'%6.2f') dep_x.UNITS ')'];
        else
          lab_1 = ['(' num2str(dep_x.data(1,comp),'%6.2f') ')'];
        end
      end
    end
  end
  ylabel(ax,sprintf('%s%s [%s]', flab, lab_1, units))
  
  if isfield(dobj.GlobalAttributes,'OBSERVATORY')
    text_s = [dobj.GlobalAttributes.OBSERVATORY{1} ' > '];
  elseif isfield(dobj.GlobalAttributes,'Source_name')
    text_s = [dobj.GlobalAttributes.Source_name{1} ' > '];
  else, text_s = '';
  end
  if isfield(dobj.GlobalAttributes,'INSTRUMENT_NAME')
    text_s = [text_s ...
      dobj.GlobalAttributes.INSTRUMENT_NAME{1} ' > '];
  elseif isfield(dobj.GlobalAttributes,'Data_type')
    text_s = [text_s ...
      dobj.GlobalAttributes.Data_type{1} ' > '];
  end
  text_s = [text_s fieldnam];
  if ~isempty(cs), text_s = [text_s ' [' shorten_cs(cs) ']']; end
  if flag_labels_is_on
    add_text(ax,text_s);
  end
  
elseif flag_spectrogram
  %% PLOT -- SPECTROGRAM
  
  dep_x=cell(size(dep.DEPEND_X,1));
  for d = 1:length(dep_x)
    dep_x{d} = getv(dobj,dep.DEPEND_X{d,1});
    dep_x{d}.s = dep.DEPEND_X{d,1};
    try
      dep_x{d}.fillv = getfillval(dobj,dep_x{d}.s);
    catch
      dep_x{d}.fillv = '';
    end
    if ~strcmp(dep_x{d}.type,'char')
      dep_x{d}.data(dep_x{d}.data==dep_x{d}.fillv) = NaN;
    end
    dep_x{d}.units = getunits(dobj,dep_x{d}.s);
    dep_x{d}.lab = getlablaxis(dobj,dep_x{d}.s);
    % check if DELTA_PLUS and  DELTA_MINUS are given
    if isfield(dep_x{d},'DELTA_PLUS') && isfield(dep_x{d},'DELTA_MINUS')
      dep_x{d}.df=struct('plus',dep_x{d}.DELTA_PLUS,'minus',dep_x{d}.DELTA_MINUS);
      if ischar(dep_x{d}.DELTA_PLUS)
        deltaplus= getv(dobj,dep_x{d}.DELTA_PLUS);
        deltaminus= getv(dobj,dep_x{d}.DELTA_MINUS);
        dep_x{d}.df.plus=deltaplus.data;
        dep_x{d}.df.minus=deltaminus.data;
      end
    else, dep_x{d}.df=[];
    end
  end
  
  % Obtain time DELTA_PLUS and  DELTA_MINUS if given
  % Also do necessary tome conversion if needed
  if strcmpi(dep.DEPEND_O.type,'tt2000')
    timeLine = EpochTT(dep.DEPEND_O.data).epochUnix; factor = 1e9;
  else, timeLine = dep.DEPEND_O.data; factor = 1;
  end
  timevar=getv(dobj,dobj.VariableAttributes.DEPEND_0{1,2});
  if isfield(timevar,'DELTA_PLUS') && isfield(timevar,'DELTA_MINUS')
    dep.dt=struct('plus',timevar.DELTA_PLUS,'minus',timevar.DELTA_MINUS);
    if ischar(timevar.DELTA_PLUS)
      deltaplus = getv(dobj,timevar.DELTA_PLUS);
      dep.dt.plus = double(deltaplus.data(1,:))/factor;
    elseif isnumeric(timevar.DELTA_PLUS)
      dep.dt.plus = double(timevar.DELTA_PLUS)/factor;
    end
    if ischar(timevar.DELTA_MINUS)
      deltaminus = getv(dobj,timevar.DELTA_MINUS);
      dep.dt.minus = double(deltaminus.data(1,:))/factor;
    elseif isnumeric(timevar.DELTA_MINUS)
      dep.dt.minus = double(timevar.DELTA_MINUS)/factor;
    end
  end
  if flag_fill_spectrogram_gaps==1 && isfield(dep,'dt') % fill gaps, disregard delta_plus and delta_minus for each data point
    dep=rmfield(dep,'dt');
  end
  if sum_dim > 0
    irf.log('notice',sprintf('Summing over dimension %d (%s)\n', ...
      sum_dim, dep_x{sum_dim}.lab));
  end
  if flag_log, plot_type='log'; else, plot_type='lin'; end
  specrec = struct('t',timeLine,'f',dep_x{1}.data,'f_unit',...
    dep_x{1}.units,'p',[],'df',dep_x{1}.df,'plot_type',plot_type);
  if isfield(dep,'dt')
    specrec.dt=dep.dt;
  end
  lab_2 ='';
  if length(dep_x)>1 && ~isempty(dep_x{comp_dim})
    if strcmp(dep_x{comp_dim}.type,'char') && strcmp(dep_x{comp_dim}.variance,'F/T')...
        && strfind(dep_x{comp_dim}.s,'LABEL_2')
      %            reclen = size(dep_x{comp_dim}.data,2)/length(dep.DEPEND_O);
      lab_2 = shiftdim(dep_x{comp_dim}.data(:,:,:),1)';
      %            lab_2 = dep_x{comp_dim}.data(:,1:reclen);
    elseif strcmp(dep_x{comp_dim}.type,'single') && ...
        (strcmp(dep_x{comp_dim}.variance,'F/T') || ...
        strcmp(dep_x{comp_dim}.variance,'T/T'))
      lab_2 = num2str(dep_x{comp_dim}.data(comp,1)',['%.2f ' dep_x{comp_dim}.units '\n']);
    else
      error('BAD type for DEPEND_X')
    end
  end
  
  if isfield(dobj.GlobalAttributes,'OBSERVATORY')
    text_s = [dobj.GlobalAttributes.OBSERVATORY{1} ' > '];
  elseif isfield(dobj.GlobalAttributes,'Source_name')
    text_s = [dobj.GlobalAttributes.Source_name{1} ' > '];
  else, text_s = '';
  end
  if isfield(dobj.GlobalAttributes,'INSTRUMENT_NAME')
    text_s = [text_s ...
      dobj.GlobalAttributes.INSTRUMENT_NAME{1} ' > '];
  elseif isfield(dobj.GlobalAttributes,'Data_type')
    text_s = [text_s ...
      dobj.GlobalAttributes.Data_type{1} ' > '];
  end
  text_s = [text_s fieldnam];
  if ~isempty(cs), text_s = [text_s ' [' shorten_cs(cs) ']']; end
  
  if isempty(comp), comp = 1; end
  ncomp = length(comp);
  h = gobjects(1,ncomp);
  if create_axes, ax = gobjects(1, ncomp); end
  
  
  if ydim > 1
    specrec.f = dep_x{ydim}.data;
    specrec.f_unit = dep_x{ydim}.units;
    specrec.df = dep_x{ydim}.df;
  end
  
  % special case for degrees
  ytick = [];
  if strcmpi(dep_x{ydim}.units,'degrees') || strcmpi(dep_x{ydim}.units,'deg')
    frange = abs(max(specrec.f(:))-min(specrec.f(:)));
    if frange > 80 && frange <=150, da = 15;
    elseif frange > 150 && frange <=200, da = 45;
    elseif frange > 200 && frange <=380, da = 90;
    else, da = [];
    end
    if ~isempty(da)
      ytick = round(min(specrec.f(:))/da):round(max(specrec.f(:))/da);
      ytick = ytick*da;
    end
  end
  
  for i=1:ncomp
    specrec.p = plot_data(i);
    if create_axes, ax(i) = irf_subplot(length(comp),1,-i); end
    h(i) = irf_spectrogram(ax(i),specrec);
    if ~isempty(ytick), set(ax(i) ,'YTick',ytick), end
    %if ~isempty(lab_2), lab_2s = [' (' lab_2(i,:) ')'];
    %else lab_2s = '';
    %end
    if flag_labels_is_on
      if ncomp<=LCOMP % Small number of components
        ylabel(h(i),sprintf('%s [%s]', dep_x{ydim}.lab, dep_x{ydim}.units))
        if ~isempty(lab_2), lab_2s = [text_s ' > ' lab_2(i,:)];
        else, lab_2s = text_s;
        end
        add_text(h(i),lab_2s);
      else % Large number of components
        if i==1, title(h(i),text_s), end
        if i==fix(ncomp/2)+1, ylabel(h(i),sprintf('%s [%s]', dep_x{ydim}.lab, dep_x{ydim}.units))
        else, ylabel('')
        end
        add_text(h(i),lab_2(i,:));
      end
    end
  end
  % Add colorbar
  if flag_colorbar_is_on
    i=fix(ncomp/2)+1;
    if isa(h(i),'handle'), hcb = colorbar(h(i)); % HG2
    else, hcb = colorbar('peer',h(i));
    end
    posCb = get(hcb,'Position');
    posAx = get(ax(i),'Position');
    dy = posAx(3);
    if ncomp>1
      set(hcb,'Position',...
        [posCb(1) posCb(2)-posCb(4)*(ncomp-fix(ncomp/2)-1) ...
        posCb(3) posCb(4)*ncomp]);
    else
      set(hcb,'TickDir','out','Position',...
        [posCb(1) posCb(2)+posCb(4)*0.05 posCb(3)*.75 posCb(4)*0.9])
    end
    set(ax(i),'Position',posAx)
    if flag_labels_is_on || flag_colorbar_label_is_manually_specified
      if ~flag_colorbar_label_is_manually_specified
        colorbar_label=[lablaxis ' [' units ']' ];
        if flag_log, colorbar_label = ['Log ' colorbar_label]; end
      end
      set(hcb,'YTickLabel','0.0'); % the distance to colorlabel defined by width of 0.0 (stupid workaround to nonfunctioning automatic distance)
      ylabel(hcb,colorbar_label);
      if flag_colorbar_label_fit_to_colorbar_height_is_on
        irf_colorbar_fit_label_height(hcb);
      end
      set(hcb,'yticklabelmode','auto');
    else
      ylabel(hcb,'');
    end
    % Resize all panels after addition of the colorbar
    if ~isempty(dy)
      for i=1:ncomp
        tt = get(ax(i),'Position');
        set(ax(i),'Position',[tt(1) tt(2) dy tt(4)])
      end
    end
  end
  set(ax(1:ncomp-1),'XTickLabel',[]);
  for i=1:1:ncomp-1, xlabel(ax(i),'');end
  
end

if nargout > 0, res = h; end

function add_text(h,txt)
text(0.99, 0.97, [' ' txt],'HorizontalAlignment','right','VerticalAlignment','top',...
  'units','normalized','fontsize',6,'parent',h)

function cs = shorten_cs(cs)

if isempty(cs), return, end

% Remove leading spaces
while cs(1) == ' ', cs(1) = []; end

if strcmpi(cs(1:3),'GSE'), cs = 'GSE'; end

% Try to correct latex
function s = corr_latex(s)
expr = {'\^-[1-3]','\^[2-3]'};
exprl = [2 1];
for i=1:length(expr)
  while 1
    ii = regexp(s,expr{i});
    if isempty(ii), break, end
    ii = ii(1);
    l = length(s);
    s_tmp = [s(1:ii) '{' s(ii+1:ii+exprl(i)) '}'];
    if l > ii+2, s = [s_tmp s(ii+exprl(i)+1:end)];
    else, s = s_tmp;
    end
  end
end

