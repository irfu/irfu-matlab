function cmap1=irf_colormap(varargin)
% IRF_COLORMAP return colormap by name or apply and freeze the colormap
%
% CMAP = IRF_COLORMAP(colormap_name)
%  Colormap_names:
%       'standard'  - (default), same as 'space','cmap' (commonly used showing space data)
%       'poynting'  - white in center and blue/green for negative and red/black for positive values
%       'poynting_gray'  - gray in center and blue/green for negative and red/black for positive values
%
% IRF_COLORMAP(AX,colormap_name) - apply colormap to axis AX
%

[ax,args,nargs] = axescheck(varargin{:});

if nargs == 0 % show only help
    help irf_colormap;
    return
end

% check which axis to apply
if isempty(ax)
    axes(gca);
else
    axes(ax(1));
end

colormap_name=args{1};

load caa/cmap.mat % default map
if nargs > 0
    switch lower(colormap_name)
        case 'poynting'
            it=0:.02:1;it=it(:);
            cmap=[ [0*it flipud(it) it];[it it it*0+1];[it*0+1 flipud(it) flipud(it)]; [flipud(it) 0*it 0*it]]; clear it;
        case {'poynting_grey','poynting_gray'}
            it=0:.02:1;it=it(:);
            cmap=[ [0*it flipud(it) it];...
				[it*.8 it*.8 it*0+1];...
				[it*0+1 flipud(it*.8) flipud(it*.8)];...
				[flipud(it) 0*it 0*it]]; 
			clear it;
        case 'solo'
            it=0:.02:1;it=it(:);
            cmap=[ [it it it*0+1];[it*0+1 flipud(it) flipud(it)]; [flipud(it) 0*it 0*it]]; clear it;
			case {'parula','jet','hsv','hot','cool','spring','summer','autumn',...
					'winter','gray','bone','copper','pink','lines','colorcube','prism','flag','white'}
            cmap = colormap(colormap_name);
      case {'bluered'}
            rr = interp1([1 64 128 192 256],[0.0  0.5 0.75 1.0 0.75],1:256);
            gg = interp1([1 64 128 192 256],[0.0  0.5 0.75 0.5 0.00],1:256);
            bb = interp1([1 64 128 192 256],[0.75 1.0 0.75 0.5 0.00],1:256);
            cmap = [rr' gg' bb'];
        case 'waterfall' % fancy-schmancy
            c = [55,137,187;...
                106,193,165;...
                172,220,166;...
                230,244,157;...
                255,254,194;...
                253,223,144;...
                251,173,104;...
                242,109,074;...
                211,064,082]/255;
            cmap = interp1(linspace(1,64,size(c,1)),c,1:64);
            
    end
end

if nargout == 0 % apply the colormap and freeze
    colormap(cmap);
    freezeColors;
    hcb = cbhandle;
    if hcb % workaround cbfreeze bug that cbfreeze removes cblabel
        hy=get(hcb,'ylabel');
        ylabel_string=get(hy,'string');
        ylabel_fontsize=get(hy,'fontsize');
        new_hcb = cbfreeze(hcb);
        new_hy=get(new_hcb,'ylabel');
        set(new_hy,'string',ylabel_string,'fontsize',ylabel_fontsize);
    end
    %    cbfreeze;
elseif nargout == 1 % only return colormap
    cmap1=cmap;
end

