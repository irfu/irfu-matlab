function [outhandle colr]=irf_pl_mark(tlim,inhandle,color);
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

% transpose tlim if needed
if size(tlim,1) == 2 && size(tlim,2) ~= 2, tlim=tlim';end; 

if nargin < 3
  if size(tlim,1) == 1
    color='yellow';
  else % choose random colors
    color=rand(size(tlim,1),3);
  end
end   

% create 1 x n x 3 color matrix
if ndims(color) == 2
  color = reshape(color, 1, size(color,1), size(color,2));
end


ud=get(gcf,'userdata');
if isfield(ud,'t_start_epoch'),  tlim=tlim-ud.t_start_epoch;end


tpoints = [tlim(:,1) tlim(:,2) tlim(:,2) tlim(:,1)];

%tlim = reshape( tlim, 1, prod(size(tlim)) );

h = reshape( inhandle, 1, prod(size(inhandle)) );
for j=1:length(h)
  axes( h(j) );
  ylim=get(h(j),'ylim');
  ypoints=zeros(size(tpoints));
  ypoints(:,1:2) = ylim(1);
  ypoints(:,3:4) = ylim(2);
  zpoints = -1 * ones(size(ypoints,1),4); % to put patches under all plots
  hp(j)=patch(tpoints', ypoints', zpoints', color,'edgecolor','none');
  set(h(j),'layer','top');
end

if nargout > 0
  outhandle = hp;
end

if nargout > 1
  colr = color;
end
  


