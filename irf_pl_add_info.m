function irf_pl_add_info
%IRF_PL_ADD_INFO add a string "created dd-mm-yy by \n user@host"
%
% $Id$

% Copyright 2004-2007 Yuri Khotyaintsev

h00 = gca;
axes('Units','normalized', ...
	'Position',[0 0 1 1], ...
	'Visible','off', ...
	'Tag','BackgroundAxes', ...
	'HitTest','off');
	created_string = [ 'Created ' date ];
if isunix
	[s,u] = unix('whoami');
	[s,h] = unix('hostname');
	u = clean_unix(u);
	h = clean_unix(h);
	created_string = [ created_string ' by ' u '@' h];
end
h1 = text('Units','normalized', ...
	'FontSize',6, ...
	'HorizontalAlignment','left', ...
	'Position',[0.01 0.97], ...
	'String',created_string, ...
	'Tag','CreatedText');

% Make some other axes active
ch = get(gcf,'children');
for j=1:length(ch)
	if ch(j)==h1, continue, end
	lasterr('')
	try
		axes(ch(j))
	catch
		% Do nothing
	end
end 
axes(h00)
end

function s = clean_unix(s)
% Help function to remove unwanted symbols before and after the string 
if isempty(s), return, end
if ~isletter(s(1)), s = s(2:end); end
if isempty(s), return, end
if ~isletter(s(end)), s = s(1:end-1); end
end