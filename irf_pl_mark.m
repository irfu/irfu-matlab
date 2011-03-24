function [outhandle colr]=irf_pl_mark(tlim,inhandle,color,varargin)
%IRF_PL_MARK   Mark that time interval(s) with color backround
%
% outhandle=irf_pl_mark(tlim,inhandle,color);
% mark that time interval(s) with color backround
% output is the handle of the patch
%
% [outhandle color]=irf_pl_mark(tlim,inhandle);
% returns also the color table
%
% tlim - time interval or array of intervals to mark
% inhandle - is the handle of axes
% color - string, rgb 1x3, nx3, or 1xnx3 specifying color(s); 
%         if omitted colors are chosen randomly.
%
% $Id$

if nargin<1, help irf_pl_mark;return;end
if nargin == 1, inhandle=gca;end
if nargin>=2, if isempty(inhandle), inhandle=gca;end;end

% transpose tlim if needed
%if size(tlim,1) == 2 && size(tlim,2) ~= 2, tlim=tlim';end; 

% if mark time instants instead of time intervals (only one time given)
if size(tlim,2) == 1, tlim(:,2)=tlim(:,1); end

if nargin < 3
  if size(tlim,1) == 1
    color='yellow';
  else % choose random colors
    color=rand(size(tlim,1),3);
  end
end   

% create 1 x n x 3 color matrix
if ischar(color),
    color = repmat(color, size(tlim,1), 1);
end


ud=get(gcf,'userdata');
if isfield(ud,'t_start_epoch'),  tlim=tlim-ud.t_start_epoch;end


tpoints = [tlim(:,1) tlim(:,2) tlim(:,2) tlim(:,1)];

%tlim = reshape( tlim, 1, prod(size(tlim)) );

h = reshape( inhandle, 1, numel(inhandle) );
hp=zeros(length(h),size(tlim,1)); % predefine patch handles
for j=1:length(h)
  ylim=get(h(j),'ylim');
  ypoints=zeros(size(tpoints));
  ypoints(:,1:2) = ylim(1);
  ypoints(:,3:4) = ylim(2);
  zpoints = -1 * ones(size(ypoints,1),4); % to put patches under all plots
  for jj=1:size(tpoints,1),
      if tlim(jj,1)==tlim(jj,2) % draw line instead of patch
          hp(j,jj)=line(tpoints(jj,1:2), ypoints(jj,[1 3]), zpoints(jj,[1 3]),'color',color(jj,:),'parent',h(j),varargin{:});
      else % make patch
          hp(j,jj)=patch(tpoints(jj,:)', ypoints(jj,:)', zpoints(jj,:)', color(jj,:),'edgecolor','none','parent',h(j),varargin{:});
          fc=get(hp(j,jj),'facecolor');
          fc=[1 1 1]-([1 1 1]-fc)/3; % make facecolor lighter
          set(hp(j,jj),'facecolor',fc);
      end
  end
  set(h(j),'layer','top');
end

if nargout > 0
  outhandle = hp;
end

if nargout > 1
  colr = color;
end
  


