function res = plot(dobj,var_s,comp,ax)
%PLOT(dobj, var_s)  plot a variable
%
% $Id$

% ----------------------------------------------------------------------------
% "THE BEER-WARE LICENSE" (Revision 42):
% <yuri@irfu.se> wrote this file.  As long as you retain this notice you
% can do whatever you want with this stuff. If we meet some day, and you think
% this stuff is worth it, you can buy me a beer in return.   Yuri Khotyaintsev
% ----------------------------------------------------------------------------

error(nargchk(2,4,nargin))

if ~ischar(var_s), error('VAR_S must be a stirng'), end

data = getv(dobj,var_s);
dim = length(data.variance(3:end));
dep = getdep(dobj,var_s);
units = getunits(dobj,var_s);
fieldnam = findva(dobj,'FIELDNAM',var_s);
%lablaxis = getlablaxis(dobj,var_s);
cs = getcs(dobj,var_s);
if ~isempty(cs), cs = [' ' cs]; end
fillv = getfillval(dobj,var_s);
data.data(data.data==fillv) = NaN;


if nargin == 2, comp = []; ax = [];
elseif nargin == 3, ax = []; 
end

if isempty(comp), use_comp = 0;
else use_comp = 1;
end

if dim == 1
	if ~isempty(dep.DEPEND_X)
		if strcmp(dep.DEPEND_X{1,2},'DEPEND_1')
			dim = 2;
		end
	end
end

if dim == 0 || dim == 1 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%% LINEAR PLOT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		
		plot_data = double(data.data)';
		if use_comp, plot_data = plot_data(:,comp); end

		if isfield(dep,'DEPEND_O')
			if ishandle(ax)
				cax = gca;
				axes(ax)
				h = irf_plot([dep.DEPEND_O plot_data]);
				axes(cax)
			else
				if ~isempty(ax), disp('AX is not an axis handle'), end
				h = irf_plot([dep.DEPEND_O plot_data]);
			end
		else
			if ishandle(ax), h = plot(ax,data.data);
			else
				if ~isempty(ax), disp('AX is not an axis handle'), end
				h = plot(data.data);
			end
		end
		
		lab_1 = '';
		if ~isempty(dep.DEPEND_X)
			dep_x_s = dep.DEPEND_X{1,1};
			dep_x = getv(dobj,dep_x_s);
			if ~isempty(dep_x)
				if strcmp(dep_x.type,'char') && strcmp(dep_x.variance,'F/T')...
						&& strcmp(dep.DEPEND_X{1,2},'LABEL_1')
					reclen=size(dep_x.data,2)/length(dep.DEPEND_O);
					if use_comp, lab_1 = ['(' dep_x.data(comp,1:reclen) ')'];
					else
						legend({dep_x.data(:,1:reclen)}, 'location','NorthWest')
						legend('boxoff')
					end
				else
					error('BAD type for DEPEND_X')
				end
			end
		end
		ylabel(sprintf('%s%s [%s]', fieldnam, lab_1, units))
		
		text_s = [dobj.GlobalAttributes.OBSERVATORY{1} ' > ' ...
			dobj.GlobalAttributes.INSTRUMENT_NAME{1} ' > ' fieldnam];
		if ~isempty(cs), text_s = [text_s ' [' shorten_cs(cs) ']']; end
		add_text(h,text_s);
		
else %%%%%%%%%%%%%%%%%%%%%%%%%%%%%% SPECTROGRAM %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	if dim == 2
		plot_data = double(data.data)';
	else
		error('plotting not implememnted')
	end
	f = getv(dobj,dep.DEPEND_X{1,1});
	fillv = getfillval(dobj,dep.DEPEND_X{1,1});
	f.data(f.data==fillv) = NaN;
	funits = getunits(dobj,dep.DEPEND_X{1,1});
	flab = getlablaxis(dobj,dep.DEPEND_X{1,1});
	if size(f.data,2) == length(dep.DEPEND_O)
		if isempty(ax), ax = gca; end
		h = caa_spectrogram(ax, dep.DEPEND_O, plot_data, f.data(:,1));
	end
	ylabel(sprintf('%s [%s]', flab, funits))
end
		
if nargout > 0, res = h; end

function add_text(h,txt)

ylim = get(h,'YLim');
xlim = get(h,'XLim');
text(xlim(2) - my_range(xlim)*.05, ylim(2) - my_range(ylim)*.05, [' ' txt],...
	'HorizontalAlignment','right')

function cs = shorten_cs(cs)

if isempty(cs), return, end

% Remove leading spaces
while cs(1) == ' ', cs(1) = []; end
	
if strcmpi(cs(1:3),'GSE'), cs = 'GSE'; end

% Help function
function r=my_range(x)
r=max(x)-min(x);

