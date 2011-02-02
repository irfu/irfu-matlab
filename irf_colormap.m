function cmap1=irf_colormap(colormap_name)
% IRF_COLORMAP return colormap by name or apply and freeze the colormap
% 
%  Colormap_names:
%       standard - in space physics (default)
%       cmap     - same as standard
%       pynting  - white in center and blue/green for negative and red/black for positive values

load caa/cmap.mat % default map
if nargin == 1, 
    switch lower(colormap_name)
        case 'poynting'
            it=0:.02:1;it=it(:);
            cmap=[ [0*it flipud(it) it];[it it it*0+1];[it*0+1 flipud(it) flipud(it)]; [flipud(it) 0*it 0*it]]; clear it;
    end
end

if nargout == 0, % apply the colormap and freeze
    colormap(cmap);
    freezeColors;
    cbfreeze;
elseif nargout == 1, % only return colormap
    cmap1=cmap;
end

