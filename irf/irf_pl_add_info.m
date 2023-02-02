function irf_pl_add_info
%IRF_PL_ADD_INFO add a string "created dd-mm-yy by \n user@host"
%

% ----------------------------------------------------------------------------
% SPDX-License-Identifier: Beerware
% "THE BEER-WARE LICENSE" (Revision 42):
% <yuri@irfu.se> wrote this file.  As long as you retain this notice you
% can do whatever you want with this stuff. If we meet some day, and you think
% this stuff is worth it, you can buy me a beer in return.   Yuri Khotyaintsev
% ----------------------------------------------------------------------------


h00 = gca;
axes('Units','normalized', ...
  'Position',[0 0 1 1], ...
  'Visible','off', ...
  'Tag','BackgroundAxes', ...
  'HitTest','off');
created_string = [ 'Created ' date ];
if isunix
  [~,u] = unix('whoami');
  [~,h] = unix('hostname');
  u = clean_unix(u);
  h = clean_unix(h);
  created_string = [ created_string ' by ' u '@' h];
end
text('Units','normalized', ...
  'FontSize',6, ...
  'HorizontalAlignment','left', ...
  'Position',[0.01 0.97], ...
  'String',created_string, ...
  'Tag','CreatedText');

axes(h00)
end

function s = clean_unix(s)
% Help function to remove unwanted symbols before and after the string
if isempty(s), return, end
if ~isletter(s(1)), s = s(2:end); end
if isempty(s), return, end
if ~isletter(s(end)), s = s(1:end-1); end
end