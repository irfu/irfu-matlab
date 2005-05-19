function irf_pl_add_info
%IRF_PL_ADD_INFO add a string "created dd-mm-yy by \n user@host"
%
% $Revision$  $Date$

% Copyright 2004 Yuri Khotyaintsev

h00 = gca;
h1 = axes('Units','normalized', ...
	'Position',[0 0 1 1], ...
	'Visible','off', ...
	'Tag','BackgroundAxes', ...
	'HitTest','off');
	created_string = [ 'Created ' date ];
if isunix
	[s,u] = unix('whoami');
	[s,h] = unix('hostname');
	u=u(1:end-1);h=h(2:end-1);
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
	end
end 
axes(h00)
