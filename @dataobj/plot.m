function res = plot(dobj,var_s,comp,ax)
%PLOT(dobj, var_s, [comp], [ax])  plot a variable
%
% comp - which component of vector to plot
% ax - axis handle where to plot 
% $Id$

% ----------------------------------------------------------------------------
% "THE BEER-WARE LICENSE" (Revision 42):
% <yuri@irfu.se> wrote this file.  As long as you retain this notice you
% can do whatever you want with this stuff. If we meet some day, and you think
% this stuff is worth it, you can buy me a beer in return.   Yuri Khotyaintsev
% ----------------------------------------------------------------------------

error(nargchk(2,4,nargin))

LCOMP = 5;

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
		flab = getlablaxis(dobj,var_s);
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
		ylabel(sprintf('%s%s [%s]', flab, lab_1, units))
		
		text_s = [dobj.GlobalAttributes.OBSERVATORY{1} ' > ' ...
			dobj.GlobalAttributes.INSTRUMENT_NAME{1} ' > ' fieldnam];
		if ~isempty(cs), text_s = [text_s ' [' shorten_cs(cs) ']']; end
		add_text(h,text_s);
		
else %%%%%%%%%%%%%%%%%%%%%%%%%%%%%% SPECTROGRAM %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	if dim == 2
		if data.dim(2)>1
			
			% This is a hack for STAFF B
			if any(any(data.data<0))
				disp('Data contains negative values!!!')
				data.data(data.data<0) = NaN;
			end
			
			ndim = data.dim(2);
			ii = 1:length(dep.DEPEND_O);
			ii = (ii-1)*ndim;
			if ~use_comp, comp=1:ndim; end
            plot_data = cell(size(comp));
			for i=comp
				plot_data{i} = data.data(:,ii+i)';
			end
		else
			plot_data = double(data.data)';
		end
	else
		error('plotting not implememnted')
	end
	f = getv(dobj,dep.DEPEND_X{1,1});
	fillv = getfillval(dobj,dep.DEPEND_X{1,1});
	f.data(f.data==fillv) = NaN;
	funits = getunits(dobj,dep.DEPEND_X{1,1});
	flab = getlablaxis(dobj,dep.DEPEND_X{1,1});
	if size(f.data,2) == length(dep.DEPEND_O)
        if isempty(ax), create_axes = 1;
        else create_axes = 0;
        end
        
 		lab_2 ='';
		if ~isempty(dep.DEPEND_X) && size(dep.DEPEND_X,1)>1
			dep_x_s = dep.DEPEND_X{2,1};
			dep_x = getv(dobj,dep_x_s);
			if ~isempty(dep_x)
				if strcmp(dep_x.type,'char') && strcmp(dep_x.variance,'F/T')...
						&& strcmp(dep.DEPEND_X{2,2},'LABEL_2')
					reclen=size(dep_x.data,2)/length(dep.DEPEND_O);
					lab_2 = dep_x.data(:,1:reclen);
                elseif strcmp(dep_x.type,'single') && strcmp(dep_x.variance,'F/T')
                    dep_x_units = getunits(dobj,dep.DEPEND_X{2,1});
                    lab_2 = num2str(dep_x.data(:,1),['%.2f ' dep_x_units]);
				else
					error('BAD type for DEPEND_X')
				end
			end
		end
		
		text_s = [dobj.GlobalAttributes.OBSERVATORY{1} ' > ' ...
			dobj.GlobalAttributes.INSTRUMENT_NAME{1} ' > ' fieldnam];
		if ~isempty(cs), text_s = [text_s ' [' shorten_cs(cs) ']']; end
		
		specrec = struct('t',dep.DEPEND_O,'f',f.data(:,1),'f_unit',funits,'p',[]);
		if isempty(comp)
			specrec.p = {plot_data};
            if create_axes, ax = gca; end
			h = caa_spectrogram(ax(1),specrec);
			if ~isempty(lab_2), lab_2 = [' (' lab_2 ')']; end
			ylabel(sprintf('%s [%s]%s', flab, funits,lab_2))
			add_text(h,text_s);
        else
            ncomp = length(comp);
			h = zeros(1,ncomp);
            % special case for degrees
            ytick = [];
            if strcmp(funits,'degrees')
                frange = abs(min(specrec.f) - max(specrec.f));
                if frange > 80 && frange <=150, da = 15;
                elseif frange > 150 && frange <=200, da = 45;
                elseif frange > 200 && frange <=380, da = 90;
                else da = [];
                end
                if ~isempty(da)
                    ytick = round(min(specrec.f)/da):round(max(specrec.f)/da);
                    ytick = ytick*da;
                end
            end
            
			for i=1:ncomp
				specrec.p = plot_data(i);
                if create_axes, ax(i) = irf_subplot(length(comp),1,-i); end
				h(i) = caa_spectrogram(ax(i),specrec);
                if ~isempty(ytick), set(ax(i) ,'YTick',ytick), end
                if ~isempty(lab_2), lab_2s = [' (' lab_2(i,:) ')'];
                else lab_2s = '';
                end
                if ncomp<=LCOMP % small number of components
                    ylabel(sprintf('%s [%s]%s', flab, funits, lab_2s))
                    if ~isempty(lab_2), lab_2s = [text_s ' > ' lab_2(i,:)];
                    else lab_2s = text_s;
                    end
                    add_text(h(i),lab_2s);
                else %large number of components
                    if i==1, title(text_s), end
                    if i==fix(ncomp/2), ylabel(sprintf('%s [%s]', flab, funits))
                    else ylabel('')
                    end
                    add_text(h(i),lab_2(i,:));
                end
			end
			set(ax(1:length(comp)-1),'XTickLabel',[])
		end
	end
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

