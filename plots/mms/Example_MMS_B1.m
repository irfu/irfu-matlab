%% Load Data
dfgFile = dataobj('/data/mms/mms1/dfg/srvy/l2pre/2015/09/mms1_dfg_srvy_l2pre_20150911_v2.0.0.cdf');
gsmB1 = get_ts(dfgFile,'mms1_dfg_srvy_l2pre_gsm');
gsmR1 = get_ts(dfgFile,'mms1_pos_gsm');

%% Plot
tintStr = '2015-09-11T07:30:00Z/2015-09-11T08:30:00Z';
h = irf_plot(gsmB1,'comp');
irf_zoom(h,'x',irf.tint(tintStr));
irf_zoom(h,'y')
add_position(h(end),gsmR1), xlabel(h(end),'')
title(h(1),tintStr)
ylabel(h(1),'Bx GSM [nT]')
ylabel(h(2),'By GSM [nT]')
ylabel(h(3),'Bz GSM [nT]')
