function outhandle=av_pl_info(string,inhandle,position);
%function outhandle=av_pl_info(string,inhandle,position);
%function outhandle=av_pl_info(string,inhandle);
% add information in string to the plot given by handle inhandle
% output is the handle of the text
% position - [x y], where [0,1] is the left upper corner
% varargin - pairs of property and proeprty values for the text
% Example:
%   ht=av_pl_info(['c_pl_sc_orientation() ' datestr(now)],gca,[0,1 ]); set(ht,'interpreter','none','FontSize', 10);

if nargin == 1, inhandle=gca;position=[.02 1];end
if nargin == 2, position=[.02 1];end
outhandle=text(0,0,string);
axes(inhandle);
set(outhandle,'Units','normalized','Position',position,'HorizontalAlignment', 'left','VerticalAlignment', 'bottom','FontSize', 5);

