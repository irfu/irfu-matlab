
set(0,'defaultAxesFontName', 'Arial')
set(0,'defaultLineLineWidth', 1.5)

scrsz = get(0,'ScreenSize');
set(0,'DefaultFigurePosition', [5 scrsz(4)*.1 scrsz(3)*.4 scrsz(4)*.8]);
clear scrsz;

irf_units % defines globa variable Units

load caa/cmap.mat
irf_colormap.standard=cmap;
irf_colormap.cmap=cmap;
irf_colormap.default=cmap;

it=0:.02:1;it=it(:); 
irf_colormap.poynting=[ [0*it flipud(it) it];[it it it*0+1];[it*0+1 flipud(it) flipud(it)]; [flipud(it) 0*it 0*it]]; clear it;
clear it

colormap(irf_colormap.standard);

