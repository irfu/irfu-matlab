function res = plot(dobj,var_s,comp,ax)
% plot(dobj, var_s) Plot a variable
%
% $Revision$  $Date$

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
lablaxis = getlablaxis(dobj,var_s);
cs = getcs(dobj,var_s);
if ~isempty(cs), cs = [' ' cs]; end
fillv = getfillval(dobj,var_s);
data.data(data.data==fillv) = NaN;

use_comp = 1;
if nargin == 2, use_comp = 0; ax = [];
elseif nargin == 3, ax = []; 
end

if dim == 1
	if ~isempty(dep.DEPEND_X)
		if strcmp(dep.DEPEND_X{1,2},'DEPEND_1')
			dim = 2;
		end
	end
end

if dim == 0
		error('cannot be plotted')
elseif dim == 1 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%% LINEAR PLOT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		
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
						legend({dep_x.data(:,1:reclen)})
					end
				else
					error('BAD type for DEPEND_X')
				end
			end
		end
		
		ylabel(sprintf('%s%s%s [%s]', lablaxis, lab_1, cs, units))
		
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
		h = caa_spectrogram(ax, dep.DEPEND_O, plot_data, f.data(:,1));
	end
	ylabel(sprintf('%s [%s]', flab, funits))
end
		
if nargout > 0, res = h; end


