function res = plot(dobj,var_s,varargin)
%PLOT(dobj, var_s, [comp], [options])  plot a variable
%
% OPTIONS - one of the following:
%	'AX'       - axis handles to use
%   'COMP'     - components to plot
%   'SUM_DIM1' - sum over first dimension (frequency, pitch angle)
%
% $Id$

% ----------------------------------------------------------------------------
% "THE BEER-WARE LICENSE" (Revision 42):
% <yuri@irfu.se> wrote this file.  As long as you retain this notice you
% can do whatever you want with this stuff. If we meet some day, and you think
% this stuff is worth it, you can buy me a beer in return.   Yuri Khotyaintsev
% ----------------------------------------------------------------------------

error(nargchk(2,14,nargin))

LCOMP = 3;

if ~ischar(var_s), error('VAR_S must be a stirng'), end

data = getv(dobj,var_s);
if isempty(data), error('VAR_S not found'), end
dim = length(data.variance(3:end));
dep = getdep(dobj,var_s);
units = corr_latex(getunits(dobj,var_s));
fieldnam = findva(dobj,'FIELDNAM',var_s);
ii = regexp(fieldnam,'_'); fieldnam(ii) = ' '; % Get rid of underscores
lablaxis = getlablaxis(dobj,var_s);
ii = regexp(lablaxis,'_'); lablaxis(ii) = ' '; % Get rid of underscores
cs = getcs(dobj,var_s);
if ~isempty(cs), cs = [' ' cs]; end
fillv = getfillval(dobj,var_s);
data.data(data.data==fillv) = NaN;


%% INPUT ARGUMENTS

have_options = 0;
arg_pos = 0;
args = varargin; 
if nargin > 2, have_options = 1; end

% Default values that can be override by options
sum_dim = 0;
use_comp = 0;
comp = [];
ax = [];
create_axes = 1;

while have_options
    arg_pos = arg_pos + 1;
    l = 1;
    if arg_pos==1 && isnumeric(args{1})
        use_comp = 1; 
        comp = args{1};
    else
        switch(lower(args{1}))
            case 'ax'
                l = 2;
                if all(ishandle(args{2}))
                    ax = args{2};
                    create_axes = 0;
                else disp('invalid value for AX')
                end
            case 'comp'
                l = 2;
                if isnumeric(args{2})
                    use_comp = 1; 
                    comp = args{2};
                else
                    disp('invalid value for COMP')
                end
            case 'sum_dim1'
                sum_dim = 1;
            otherwise
                disp('unknown argument')
                break
        end
    end
	args = args(l+1:end);
	if isempty(args), break, end
end

if dim == 1
	if ~isempty(dep.DEPEND_X)
		if strcmp(dep.DEPEND_X{1,2},'DEPEND_1')
			dim = 2;
		end
	end
end

%% PLOT
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
	
    dep_x_1 = getv(dobj,dep.DEPEND_X{1,1});
    dep_x_1.s = dep.DEPEND_X{1,1};
	dep_x_1.fillv = getfillval(dobj,dep_x_1.s);
    dep_x_1.data(dep_x_1.data==dep_x_1.fillv) = NaN;
	dep_x_1.units = getunits(dobj,dep_x_1.s);
	dep_x_1.lab = getlablaxis(dobj,dep_x_1.s);
    specrec = struct('t',dep.DEPEND_O,'f',dep_x_1.data(:,1),'f_unit',dep_x_1.units,'p',[]);
    
    if size(dep.DEPEND_X,1)>1
        dep_x_2 = getv(dobj,dep.DEPEND_X{2,1});
        dep_x_2.s = dep.DEPEND_X{2,1};
        dep_x_2.fillv = getfillval(dobj,dep_x_2.s);
        dep_x_2.data(dep_x_2.data==dep_x_2.fillv) = NaN;
        dep_x_2.units = getunits(dobj,dep_x_2.s);
        dep_x_2.lab = getlablaxis(dobj,dep_x_2.s);
    else
        dep_x_2 = [];
    end
    
    if dim == 2
		if data.dim(2)>1
			
			% This is a hack for STAFF B
			if any(any(data.data<=0))
				disp('Data contains negative & zero values!!!')
				data.data(data.data<0) = NaN;
			end
			
			ndim = data.dim(2);
			ii = 1:length(dep.DEPEND_O);
			ii = (ii-1)*ndim;
			if ~use_comp, comp=1:ndim; end
            plot_data = cell(size(comp));
            for i=1:length(comp)
                plot_data{i} = data.data(:,ii+comp(i))';
            end
            if sum_dim==1 % sum over frequency, pitch angle
                pd_tmp = zeros(length(dep.DEPEND_O),length(comp));
                for i=1:length(comp)
                    plot_data{i}(isnan(plot_data{i})) = 0;
                    pd_tmp(:,i) = sum(plot_data{i},2);
                end
                clear plot_data
                plot_data{1} = pd_tmp; clear pd_tmp
                dep_x_1 = dep_x_2; dep_x_2 = [];
                specrec.f = dep_x_1.data(comp,1);
                specrec.f_unit = dep_x_1.units;
                comp = [];
            end
		else
			plot_data{1} = double(data.data)';
		end
	else
		error('plotting not implememnted')
    end
        
    if size(dep_x_1.data,2) ~= length(dep.DEPEND_O)
        error('bad size for DEPEND_X_1')
    end
     
    lab_2 ='';
    if ~isempty(dep_x_2)
        if strcmp(dep_x_2.type,'char') && strcmp(dep_x_2.variance,'F/T')...
                && strcmp(dep_x_2.s,'LABEL_2')
            reclen = size(dep_x_2.data,2)/length(dep.DEPEND_O);
            lab_2 = dep_x_2.data(:,1:reclen);
        elseif strcmp(dep_x_2.type,'single') && ...
                (strcmp(dep_x_2.variance,'F/T') || strcmp(dep_x_2.variance,'T/T'))
            lab_2 = num2str(dep_x_2.data(comp,1),['%.2f ' dep_x_2.units]);
        else
            error('BAD type for DEPEND_X')
        end
    end
    
    text_s = [dobj.GlobalAttributes.OBSERVATORY{1} ' > ' ...
        dobj.GlobalAttributes.INSTRUMENT_NAME{1} ' > ' fieldnam];
    if ~isempty(cs), text_s = [text_s ' [' shorten_cs(cs) ']']; end
    
    if isempty(comp), comp = 1; end
    ncomp = length(comp);
    h = zeros(1,ncomp);
    if create_axes, ax = zeros(1, ncomp); end
    
    % special case for degrees
    ytick = [];
    if strcmpi(dep_x_1.units,'degrees')
        frange = abs(my_range(specrec.f));
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
        if create_axes, ax(i) = irf_subplot(length(comp),1,-i); end %#ok<AGROW>
        h(i) = caa_spectrogram(ax(i),specrec);
        if ~isempty(ytick), set(ax(i) ,'YTick',ytick), end
        %if ~isempty(lab_2), lab_2s = [' (' lab_2(i,:) ')'];
        %else lab_2s = '';
        %end
        if ncomp<=LCOMP % Small number of components
            ylabel(sprintf('%s [%s]', dep_x_1.lab, dep_x_1.units))
            if ~isempty(lab_2), lab_2s = [text_s ' > ' lab_2(i,:)];
            else lab_2s = text_s;
            end
            add_text(h(i),lab_2s);
        else % Large number of components
            if i==1, title(text_s), end
            if i==fix(ncomp/2), ylabel(sprintf('%s [%s]', dep_x_1.lab, dep_x_1.units))
            else ylabel('')
            end
            add_text(h(i),lab_2(i,:));
        end
        % Add colorbar
        if i==fix(ncomp/2)+1
            hcb = colorbar;
            dy = get(ax(i),'Position'); dy = dy(3);
            pcb = get(hcb,'Position');
            ylabel(hcb,['Log ' lablaxis ' [' units ']' ])
            set(hcb,'Position',[pcb(1) pcb(2)-pcb(4)*(ncomp-fix(ncomp/2)-1) pcb(3) pcb(4)*ncomp])
        end
    end
    % Resize all panels aftre addition of a colorbar
    if ~isempty(dy)
        for i=1:ncomp
            tt = get(ax(i),'Position');
            set(ax(i),'Position',[tt(1) tt(2) dy tt(4)])
        end
    end
    set(ax(1:length(comp)-1),'XTickLabel',[])

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

% Try to correct latex
function s = corr_latex(s)

expr = {'\^-[1-3]','\^[2-3]'};
exprl = [2 1];
for i=1:length(expr)
    while 1
        ii = regexp(s,expr{i});
        if isempty(ii), break, end
        ii = ii(1);
        l = length(s);
        s_tmp = [s(1:ii) '{' s(ii+1:ii+exprl(i)) '}'];
        if l > ii+2, s = [s_tmp s(ii+exprl(i)+1:end)];
        else s = s_tmp;
        end
    end
end

