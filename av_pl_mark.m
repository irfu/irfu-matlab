function outhandle=av_pl_mark(tlim,inhandle,color);
%function outhandle=av_pl_mark(tlim,inhandle,color);
% mark that time interval with color backround
% output is the handle of the patch
% tlim - time interval to mark
% inhandle - is the handle of axes
% color - string specifying color

warning('caa:cleanup',...
'Function %s is deprecated and will be removed on May 1, 2004.\nUse %s instead',...
mfilename,'irf_pl_mark')

if nargin<1, help av_pl_mark;end
if nargin == 1, inhandle=gca;end
if nargin < 2, color='yellow';end

tlim = reshape( tlim, 1, prod(size(tlim)) );
h = reshape( inhandle, 1, prod(size(inhandle)) );
for j=1:length(h)
  axes( h(j) );
  ylim=get(h(j),'ylim');
  hp=patch([tlim tlim(2) tlim(1)],[ylim(1) ylim(1) ylim(2) ylim(2)],[-1 -1 -1 -1],color,'facecolor',color,'edgecolor','none');
  set(h(j),'layer','top');
end


