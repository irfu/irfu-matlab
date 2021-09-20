function [outhandle, colr]=irf_pl_mark(varargin)
%IRF_PL_MARK   Mark time intervals or instants
% Marks time intervals with transparent rectangle or time instants with line
%
%   IRF_PL_MARK(tlim) mark time interval
%
%   IRF_PL_MARK(AX,tlim) mark in the specified axis (default gca)
%
%   IRF_PL_MARK(AX,tlim,color) mark with the specified color
%
%   IRF_PL_MARK(AX,tlim,color,'Property1',PropertyValue1,..) specify properties
%
%   [H,COLOR]=IRF_PL_MARK(...) returns handle to patch or line and the color table
%
% tlim - time interval or array of intervals to mark
%        or time intervals as string (ISO 8601 time interval 1. format)
%        or column vector of time instants to draw lines
%        or a GenericTime interval (as created by irf.tint)
% color - string, rgb 1x3, nx3, or 1xnx3 specifying color(s);
%         if omitted colors are chosen randomly.
%
%
% WARNING!!! IRF_PL_MARK has changed (2011-03-30) the order of input parameters
%  now more compliant to MATLAB


[ax,args,nargs] = irf.axescheck(varargin{:});

if nargs == 0 % show only help
  help irf_pl_mark;
  return
end

if isempty(ax)
  if any(ishandle(args{1})) % first argument is axis handles
    ax=args{1};
    args=args(2:end);
    nargs=nargs-1;
  else
    % call irf_pl_mark recursively with GCA as handle
    [H,COLOR]= irf_pl_mark(gca,varargin{:});
    if nargout > 0, outhandle = H; end
    if nargout > 1, colr = COLOR; end
    return;
  end
end

tlim=args{1};
% if mark time instants instead of time intervals (only one time given)
if ischar(tlim) % time interval specified in iso format
  tlim=irf_time(tlim, 'utc>tint');
end
if(isa(tlim,'GenericTimeArray'))
  % Convert to epochUnix tint as is used by this function
  tlim = [tlim.start.epochUnix, tlim.stop.epochUnix];
end
if size(tlim,2) == 1, tlim(:,2)=tlim(:,1); end

if nargs == 1 % if only time intervals given, specify color
  if size(tlim,1) == 1
    color='yellow';
  else % choose random colors
    color=rand(size(tlim,1),3);
  end
end

if nargs >= 2 % color i specified
  color = args{2};
end

if nargs > 2 && (rem(nargs,2) ~= 0)
  error('IRFU_MATLAB:irf_pl_mark:InvalidNumberOfInputs','Incorrect number of input arguments')
end

% properties specified
pvpairs = args(3:end);


% create 1 x n x 3 color matrix
if ischar(color)
  color = repmat(color, size(tlim,1), 1);
elseif size(color,1)<size(tlim,1)
  color=repmat(color(1,:),size(tlim,1),1);
end


fig = get(ax(1),'Parent');
ud = get(fig,'userdata');
if isfield(ud,'t_start_epoch'),  tlim=tlim-ud.t_start_epoch;end


tpoints = [tlim(:,1) tlim(:,2) tlim(:,2) tlim(:,1)];

%tlim = reshape( tlim, 1, prod(size(tlim)) );

h = reshape( ax, 1, numel(ax) );
hp=gobjects(length(h),size(tlim,1)); % predefine patch handles
for j=1:length(h)
  ylim=get(h(j),'ylim');
  ypoints=zeros(size(tpoints));
  ypoints(:,1:2) = ylim(1);
  ypoints(:,3:4) = ylim(2);
  zpoints = zeros(size(ypoints,1),4); % to put patches under all plots
  for jj=1:size(tpoints,1)
    if tlim(jj,1)==tlim(jj,2) % draw line instead of patch (in this case draw line above everything, therefore "+2" in the next line)
      hp(j,jj)=line(tpoints(jj,1:2), ypoints(jj,[1 3]), zpoints(jj,[1 3]),'color',color(jj,:),'parent',h(j),pvpairs{:});
      set(hp(j,jj),'Tag','irf_pl_mark')
    else % make patch
      %
      % transparency yet does not work on spectrograms (work only in
      % opengl renderer mode in which other things does not work).
      % therefore we put interval marking only in figures where are no
      % spectrograms, patches or surface plots (except irf_pl_mark marking
      % itself)
      %
      if ~isempty(findobj(h(j),'tag','irf_pl_mark')) || ~any(~isempty([findobj(h(j),'Type','surface') findobj(h(j),'Type','patch')])) % put mark under everything
        hold(h(j),'on');
        hp(j,jj)=patch(tpoints(jj,:)', ypoints(jj,:)', zpoints(jj,:)', color(jj,:),'edgecolor','none','parent',h(j),'facealpha',1,'tag','irf_pl_mark',pvpairs{:});
        set(h(j),'children',circshift(get(h(j),'children'),-1)); % move patch to be the first children (below other plots)
        fc=get(hp(j,jj),'facecolor');
        fc=[1 1 1]-([1 1 1]-fc)*0.5; % make facecolor lighter
        set(hp(j,jj),'facecolor',fc);
      end
    end
  end
  set(h(j),'layer','top');
end

% do not display the mark as a legend entry
% hp might not have been set if used on a spectrogram
if isgraphics(hp)
  hMarks = get(hp,'Annotation');
  if numel(hMarks) > 1
    for iMark = 1:numel(h)
      set(get(hMarks{iMark},'LegendInformation'),'IconDisplayStyle','off')
    end
  else
    set(get(hMarks,'LegendInformation'),'IconDisplayStyle','off')
  end
end



if nargout > 0
  outhandle = hp;
end

if nargout > 1
  colr = color;
end



