function outhandle=irf_pl_mark(tlim,inhandle,color);
%IRF_PL_MARK   Mark that time interval with color backround
%
% outhandle=irf_pl_mark(tlim,inhandle,color);
% mark that time interval with color backround
% output is the handle of the patch
% tlim - time interval to mark
% inhandle - is the handle of axes
% color - string specifying color
%
% $Id$

if nargin<1, help irf_pl_mark;return;end
if nargin == 1, inhandle=gca;end
if nargin < 2, color='yellow';end

ud=get(gcf,'userdata');
if isfield(ud,'t_start_epoch'),  tlim=tlim-ud.t_start_epoch;end
  

tlim = reshape( tlim, 1, prod(size(tlim)) );
h = reshape( inhandle, 1, prod(size(inhandle)) );
for j=1:length(h)
  axes( h(j) );
  ylim=get(h(j),'ylim');
  hp=patch([tlim tlim(2) tlim(1)],[ylim(1) ylim(1) ylim(2) ylim(2)],[-1 -1 -1 -1],color,'facecolor',color,'edgecolor','none');
  set(h(j),'layer','top');
end


