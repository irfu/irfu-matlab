
set(0,'defaultAxesFontName', 'Arial')
set(0,'defaultLineLineWidth', 1.5)

scrsz = get(0,'ScreenSize');
set(0,'DefaultFigurePosition', [5 scrsz(4)*.1 scrsz(3)*.4 scrsz(4)*.8]);
clear scrsz;

Units=irf_units; % defines globa variable Units

%irf_colormap('standard'); % defines default colormap

