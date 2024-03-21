%% Pitch angle distributions
mms_id = 1;
tint_zoom = irf.tint('2017-07-06T13:54:05.50Z/2017-07-06T13:54:05.65Z') + [-1 1];
tint_pa = irf.tint('2017-07-06T13:54:05.50Z/2017-07-06T13:54:05.65Z');

c_eval('ePitch?_! = ePDist?.tlim(tint_pa).pitchangles(dmpaB?,180 -11.25*[!:-1:0]);',mms_id,1:16)
c_eval('ePitch?_! = ePDist?.tlim(tint_pa).pitchangles(dmpaB?,0 + 11.25*[0:1:!]);',mms_id,1:16)

%c_eval('ePitch?_!bins = ePDist?.tlim(tint_zoom).pitchangles(dmpaB?,180 + [-11.25*! 0]);',1:4,1:4)
%c_eval('eFlux?_4 = ePitch?_4.flux;',1:4)
%c_eval('ePitch?_16 = ePDist?.tlim(tint_zoom).pitchangles(dmpaB?,0:11.25:180);',1:4,1:4)
%c_eval('eFlux?_fov_!bins = ePDist?.tlim(tint_zoom).flux.pitchangles(dmpaB?,180 + [-11.25*! 0]);',1:4,1:4)
%c_eval('eFlux?_4_Eedi = irf.ts_scalar(eFlux?_4.time,squeeze(eFlux?_4.elim(E_edi).data));',1:4)

%%
h = setup_subplots(4,4,'horizontal');
isub = 0;
for ipitch = 1:16  
  variable = eval(sprintf('ePitch%.0f_%.0f',mms_id,ipitch));
  isub = isub + 1; hca = h(isub);
  hpitch = variable.plot_pad_polar(hca,'tint',tint_pa,'scpot',scPot1.resample(variable));  
  %isub = isub + 1; hca = h(isub);
  %hpitch = variable.elim([100 600]).plot_pad_polar(hca,'tint',tint_pa,'scpot',scPot1.resample(variable));  
end
hl_CLim = linkprop(h,'CLim');
hl_XLim = linkprop(h,'XLim');
hl_YLim = linkprop(h,'YLim');

h(1).XLim = [-5 0];
h(1).YLim = [-5 5];
colormap(irf_colormap('waterfall'))
%colormap(pic_colors('candy'))