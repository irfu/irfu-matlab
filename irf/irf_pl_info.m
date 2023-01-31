function outhandle=irf_pl_info(string,inhandle,position)
%IRF_PL_INFO   Add information to the plot
%
% Will be removed!!! Use IRF_LEGEND instead!
%
% outhandle=irf_pl_info(string,inhandle,position);
% outhandle=irf_pl_info(string,inhandle);
% add information in string to the plot given by handle inhandle
% output is the handle of the text
% position - [x y], where [0,1] is the left upper corner
% if inhandle is 0 use the whole figure area as coordinates
%
% Example:
%   ht=irf_pl_info(['c_pl_sc_orient() ' char(datetime("now","Format","dd-MMM-uuuu HH:mm:ss"))],gca,[0,1 ]);
%   set(ht,'interpreter','none','FontSize', 10);
%

disp('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!');
disp('IRF_PL_INFO will be removed! Use IRF_LEGEND instead.');
disp('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!');

if nargin == 1, inhandle=gca;position=[.02 1];end
if nargin == 2, position=[.02 1];end
if inhandle ==0 % add to the whole figure if inhandle = 0
  h00 = gca;
  h1 = axes('Units','normalized', 'Position',[0 0 1 1], 'Visible','off', ...
    'Tag','BackgroundAxes', 'HitTest','off');
  outhandle = text('Units','normalized','FontSize',6,'HorizontalAlignment','left', ...
    'Position',[0.01 0.97], 'String', string, 'Tag','CreatedText');
  axes(h00)
else% add to the given or current axis
  axes(inhandle);
  outhandle=text(0,0,string);
  set(outhandle,'Units','normalized','Position',position, ...
    'HorizontalAlignment', 'left','VerticalAlignment', 'bottom','FontSize', 5);
end
