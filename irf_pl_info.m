function outhandle=irf_pl_info(string,inhandle,position);
%IRF_PL_INFO   Add information to the plot
%
% outhandle=irf_pl_info(string,inhandle,position);
% outhandle=irf_pl_info(string,inhandle);
% add information in string to the plot given by handle inhandle
% output is the handle of the text
% position - [x y], where [0,1] is the left upper corner
%
% Example:
%   ht=irf_pl_info(['c_pl_sc_orient() ' datestr(now)],gca,[0,1 ]); 
%   set(ht,'interpreter','none','FontSize', 10);
%
% $Id$

if nargin == 1, inhandle=gca;position=[.02 1];end
if nargin == 2, position=[.02 1];end
if inhandle ==0, % add to the whole figure if inhandle = 0
    h0 = gca;
    h1 = axes('Units','normalized', 'Position',[0 0 1 1], 'Visible','off', ...
        'Tag','BackgroundAxes', 'HitTest','off');
    outhandle = text('Units','normalized','FontSize',6,'HorizontalAlignment','left', ...
        'Position',[0.01 0.97], 'String', string, 'Tag','CreatedText');
    axes(h0)
else% add to the given or current axis
    outhandle=text(0,0,string);
    axes(inhandle);
    set(outhandle,'Units','normalized','Position',position, ...
        'HorizontalAlignment', 'left','VerticalAlignment', 'bottom','FontSize', 5);
end
