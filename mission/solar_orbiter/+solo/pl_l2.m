function hout = pl_l2(fName)
%SOLO.PL_L2  Fast plot L2 data
%
% H = PL_L2(fName)

%fName = 'solo_L2_rpw-lfr-surv-cwf-e-cdag_20200312_V01.cdf';
%% load data
d =dataobj(fName);

El2 = get_ts(d,'EDC');
Vl2 = get_ts(d,'VDC');

%% plot data

h = irf_plot(2,'reset');

irf_plot(h(1),El2);
legend(h(1),'V12','V13','V23')
[~,fName_s,~] = fileparts(fName);
title(h(1),fName_s,'Interpreter','none')

irf_plot(h(2),Vl2);
legend(h(2),'V1','V2','V3')

%%
if nargout, hout = h; end 