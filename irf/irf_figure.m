function c=irf_figure(varargin)
%IRF_FIGURE  Create a new figure or reset existing figure
%  
%  H = IRF_FIGURE(NPANELS)
%  Create new figure with NPANELS new panels
%  Will use current figure, if empty.
%
%  H = IRF_FIGURE(FIGURE_NUMBER,NPANELS)
%  Create new figure with FIGURE_NUMBER and NPANELS new panels 
%
%  H = IRF_FIGURE(HF,NPANELS)
%  Initialize existing figure with handle HF with NPANELS new panels 
%  
%  H = IRF_FIGURE(HF,NPANELS,'reset')
%  Initialize figure HF with NPANELS new panels
%  and reset it's size and other properties to default
%
% See also: irf_panel, irf_plot

narginchk(1, 3)

args = varargin;
nargs = nargin;
hcf=[]; flag = '';
switch nargs
  case 1, nSubplots = args{1};
  case 2, hcf = args{1}; nSubplots = args{2};
  case 3, hcf = args{1}; nSubplots = args{2}; flag = args{3};
  otherwise, error('invalid number of input arguments')
end
    
if ~isnumeric(nSubplots), error('nSubplots must be a positive number'), end
if length(nSubplots)>1 || nSubplots<1 || nSubplots>20
  error('nSubplots must be >=1 && <=20')
end

flagReset = false;
if ~ischar(flag), error('FLAG must be a string'), end
switch lower(flag)
  case ''
  case 'reset', flagReset = true;
  otherwise, error('unrecognized FLAG')
end

if isempty(hcf) 
  if isempty(get(0,'CurrentFigure')) || ~isempty(get(gcf,'children')) 
    hcf = figure; flagReset = true;
  else
    irf.log('notice','Using current empty figure (no reset)')
    hcf = gcf;
  end
else
  if numel(hcf)>1 || ~isnumeric(hcf)
    error('HF must be a valid figure handle or figure number')
  end
  if ~ishghandle(hcf,'figure'), flagReset = true; end
  hcf = figure(hcf);
end

if nSubplots>=1 && nSubplots<=20
  nSubplots = floor(nSubplots);
  c = gobjects(1,nSubplots);
  if flagReset
		xSize = 11;
		ySize = 5+5*sqrt(nSubplots);
		xLeft = (21-xSize)/2; yTop = (30-ySize)/2;
		set(hcf,'PaperPosition',[xLeft yTop xSize ySize])
		un = get(0,'units');
		set(0,'units','pixels');
		sz = get(0,'screensize');
		xx = min(min(700,sz(3))/xSize,min(900,sz(4))/ySize); % figure at least 600 wide or 900 height but not outside screen
		set(hcf,'Position',[10 10 xSize*xx ySize*xx])
		set(0,'units',un);
		clear xSize sLeft ySize yTop
    % XXX not sure this belongs here
    set(hcf,'color','white');
    set(hcf,'defaultAxesColorOrder',[0 0 0;0 0 1;1 0 0;0.3 0.3 0.3;0 1 1 ;1 0 1; 1 1 0])
	end
  clf;
  all_axis_position=[0.17 0.1 0.9 0.95]; % xmin ymin xmax ymax
  subplot_width=all_axis_position(3)-all_axis_position(1);
  subplot_height=(all_axis_position(4)-all_axis_position(2))/nSubplots;
  for j=1:nSubplots
    c(j)=axes('position',[all_axis_position(1) ...
      all_axis_position(4)-j*subplot_height ...
      subplot_width subplot_height]); % [x y dx dy]
    cla(c(j));
    set(c(j),'box','on','tag','');
  end
  user_data = get(gcf,'userdata');
  user_data.subplot_handles = c;
  user_data.current_panel = 0;
  set(hcf,'userdata',user_data);
  figure(hcf); % bring figure to front
end
end