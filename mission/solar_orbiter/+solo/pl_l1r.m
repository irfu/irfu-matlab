function hout = pl_l1r(fName)
%SOLO.PL_L2  Fast plot L2 data
%
% H = PL_L2(fName)

%fName = 'solo_L2_rpw-lfr-surv-cwf-e-cdag_20200312_V01.cdf';
d =dataobj(fName);

muxSet = get_ts(d,'BIAS_MODE_MUX_SET');
El1 = get_ts(d,'E');
Vl1 = get_ts(d,'V');

Tint = irf.tint(El1.time);

%%

h = irf_plot(3,'reset');

irf_plot(h(1),El1,'.');
ylabel(h(1),'E')
%legend(h(1),'V12','V13','V23')
title(h(1),fName,'Interpreter','none')

irf_plot(h(2),Vl1,'.');
ylabel(h(2),'V')
%legend(h(2),'V1','V2','V3')

irf_plot(h(3),muxSet,'.');
ylabel(h(3),'mux')

irf_zoom(h,'x',Tint)

if nargout, hout = h; end