% A routine that identifies strong currents from brst data from 1 Nov 2015
% to present and then plots E, B, J, EdotJ for the specific intervals.
% Written by Elin Eriksson
%
% Loads data and calculates currents (curlometer method).
%Plots B, J (curlometer method), E, JxB electric field, and J.E as an average
%of all spacecraft values and a seperate plot with the same variables where only the E and
%B are from a specific spacecraft given by the number of ic.

%% Set time interval to look in, the specific spacecraft of interest and set limits for Null method.
% taken from
ic=1; %Gives number of spacecraft where density is taken for Hall field calculations.
currentLim=500E-9;
intervalStart=irf_time([2015 11 01 00 00 00]);

currentIntervals=mms.strong_current_search_brst(ic,currentLim,intervalStart);
%%
for i=1:length(currentIntervals(:,1))
  Tint  = irf.tint(EpochUnix(currentIntervals(i,1)),EpochUnix(currentIntervals(i,2)));
  % Load magnetic field and spacecraft positional data
  % Magnetic Field
  disp('Loading Magnetic fields');
  c_eval('B?=mms.db_get_ts(''mms?_dfg_brst_l2pre'',''mms?_dfg_brst_l2pre_gse'',Tint);',1:4);
  if isempty(B1) || isempty(B2) || isempty(B3) || isempty(B4)
    continue
  else
    %c_eval('B?=removerepeatpnts(B?)',[1:4]); %Removes data points with too small a time difference between them which buggs the resample functions
    timeStart=max([B1.time.start.epochUnix B2.time.start.epochUnix B3.time.start.epochUnix B4.time.start.epochUnix],[],2);
    timeStart=EpochUnix(timeStart);
    timeStop=min([B1.time.stop.epochUnix B2.time.stop.epochUnix B3.time.stop.epochUnix B4.time.stop.epochUnix],[],2);
    timeStop=EpochUnix(timeStop);
    timeLog=B1.time >= timeStart & B1.time <=timeStop;
    newTime=B1.time(timeLog,:);
    c_eval('B? = B?.resample(newTime);',1:4);
  end
  % Spacecraft Position
  disp('Loading Spacecraft Position');
  R  = mms.get_data('R_gse',Tint);%Cailbrated position
  if length(R.gseR1(1,:))==4 || length(R.gseR2(1,:))==4 || length(R.gseR3(1,:))==4 || length(R.gseR4(1,:))==4
    c_eval('R? =irf.ts_vec_xyz(R.time, R.gseR?(:,1:3));',1:4);
  else
    c_eval('R? =irf.ts_vec_xyz(R.time, R.gseR?);',1:4);
  end
  clear R
  % Checks if there is any position data missing in that case go to predicted
  % spacecraft position
  if isempty(R1) || isempty(R2) || isempty(R3) || isempty(R4)
    disp('Loading predicted spacecraft position');
    c_eval('R?=mms.db_get_ts(''mms?_dfg_srvy_ql'',''mms?_ql_pos_gse'',Tint);',1:4) %Predicted s/c position
  end
  %If the spacecraft data is still missing for a spacecraft give error
  if isempty(R1) || isempty(R2) || isempty(R3) || isempty(R4)
    error('Missing spacecraft position from at least one spacecraft');
  else
    c_eval('R? = R?.resample(B1);',1:4);
  end

  %% Average calculations. Loads electric fields, density and calculates J and JxB with curlometer method


  disp('Loads Electric fields')
  % Electric field
  c_eval('Efield=mms.db_get_ts(''mms?_edp_brst_ql_dce2d'',''mms?_edp_dce_xyz_dsl'',Tint);',ic);
  %Removes 1.5 mV/m offset from Ex
  Efield.data(:,1) = Efield.data(:,1)-1.5;
  c_eval('Bfield=B?;',ic);

  disp('Load electron fpi data')
  c_eval('eEnSp_pX=mms.db_get_variable(''mms?_fpi_fast_sitl'',''mms?_fpi_eEnergySpectr_pX'',Tint);',ic);%positiveX
  c_eval('eEnSp_pY=mms.db_get_variable(''mms?_fpi_fast_sitl'',''mms?_fpi_eEnergySpectr_pY'',Tint);',ic);
  c_eval('eEnSp_pZ=mms.db_get_variable(''mms?_fpi_fast_sitl'',''mms?_fpi_eEnergySpectr_pZ'',Tint);',ic);
  c_eval('eEnSp_mX=mms.db_get_variable(''mms?_fpi_fast_sitl'',''mms?_fpi_eEnergySpectr_mX'',Tint);',ic);%negativeX
  c_eval('eEnSp_mY=mms.db_get_variable(''mms?_fpi_fast_sitl'',''mms?_fpi_eEnergySpectr_mY'',Tint);',ic);
  c_eval('eEnSp_mZ=mms.db_get_variable(''mms?_fpi_fast_sitl'',''mms?_fpi_eEnergySpectr_mZ'',Tint);',ic);
  disp('Load ion fpi data')
  c_eval('iEnSp_pX=mms.db_get_variable(''mms?_fpi_fast_sitl'',''mms?_fpi_iEnergySpectr_pX'',Tint);',ic);%positiveX
  c_eval('iEnSp_pY=mms.db_get_variable(''mms?_fpi_fast_sitl'',''mms?_fpi_iEnergySpectr_pY'',Tint);',ic);
  c_eval('iEnSp_pZ=mms.db_get_variable(''mms?_fpi_fast_sitl'',''mms?_fpi_iEnergySpectr_pZ'',Tint);',ic);
  c_eval('iEnSp_mX=mms.db_get_variable(''mms?_fpi_fast_sitl'',''mms?_fpi_iEnergySpectr_mX'',Tint);',ic);%negativeX
  c_eval('iEnSp_mY=mms.db_get_variable(''mms?_fpi_fast_sitl'',''mms?_fpi_iEnergySpectr_mY'',Tint);',ic);
  c_eval('iEnSp_mZ=mms.db_get_variable(''mms?_fpi_fast_sitl'',''mms?_fpi_iEnergySpectr_mZ'',Tint);',ic);
  disp('Makes ion energy omni directional')
  if  isempty(eEnSp_pX) || isempty(eEnSp_pY) || isempty(eEnSp_pZ) || isempty(eEnSp_mX) || isempty(eEnSp_mY) || isempty(eEnSp_mZ) || isempty(iEnSp_pX) || isempty(iEnSp_pY) || isempty(iEnSp_pZ) || isempty(iEnSp_mX) || isempty(iEnSp_mY) || isempty(iEnSp_mZ)
    disp('No fpi data available')
    continue
  else
    tmp1 = (iEnSp_pX.data + iEnSp_pY.data + iEnSp_pZ.data + iEnSp_mX.data + iEnSp_mY.data + iEnSp_mZ.data); %/6;
    [~,energy] = hist([log10(10),log10(30e3)],32);
    energy = 10.^energy;
    speciEnSp = struct('t', irf_time(iEnSp_pX.DEPEND_0.data, 'ttns>epoch'));
    speciEnSp.p = double(tmp1);
    speciEnSp.p_label={'log'};
    speciEnSp.f_label = {'Energy (eV)'};
    speciEnSp.f = single(energy);

    disp('Makes electron energy omni directional')
    tmp2 = (eEnSp_pX.data + eEnSp_pY.data + eEnSp_pZ.data + eEnSp_mX.data + eEnSp_mY.data + eEnSp_mZ.data); %/6;
    [~,energy] = hist([log10(10),log10(30e3)],32);
    energy = 10.^energy;
    speceEnSp = struct('t', irf_time(eEnSp_pX.DEPEND_0.data, 'ttns>epoch'));
    speceEnSp.p = double(tmp2);
    speceEnSp.p_label={'log'};
    speceEnSp.f_label = {'Energy (ev)'};
    speceEnSp.f = single(energy);
  end
  disp('Calculates current')
  % Assuming GSE and DMPA are the same coordinate system calculates j and jXB.
  c_eval('R.C? = R?;',1:4);
  c_eval('B.C? = B?;',1:4);
  curlB = c_4_grad(R,B,'curl');
  j     = curlB.data./1.0e3.*1e-9./(4*pi*1e-7); %A/m^2 if B in nT and R in km

  j=[B1.time.epochUnix j];
  j(:,2:4) = j(:,2:4).*1e9; %nA/m^2


  %% Plot data should only be used if nulls are found. Plots are only focused around the intervals of nulls found because of above segment.
  disp('Plotting figure')
  h = irf_plot(5,'newfigure');

  hca = irf_panel('BMMSav');
  irf_plot(hca,Bfield);
  ylabel(hca,{'B_{DMPA}','(nT)'},'Interpreter','tex');
  irf_legend(hca,{'B_{x}','B_{y}','B_{z}'},[0.88 0.10])
  %irf_legend(hca,'(a)',[0.99 0.98],'color','k')
  grid(hca,'off');
  set(hca,'xticklabel',[]);

  hca = irf_panel('J');
  irf_plot(hca,j);
  ylabel(hca,{'J_{DMPA}','(nA m^{-2})'},'Interpreter','tex');
  irf_legend(hca,{'J_{x}','J_{y}','J_{z}'},[0.88 0.10])
  %irf_legend(hca,'(c)',[0.99 0.98],'color','k')
  grid(hca,'off');
  set(hca,'xticklabel',[]);

  hca = irf_panel('EMMS');
  irf_plot(hca,Efield);
  ylabel(hca,{'E_{DSL}','(mV m^{-1})'},'Interpreter','tex');
  irf_legend(hca,{'E_{x}','E_{y}','E_{z}'},[0.88 0.10])
  %irf_legend(hca,'(b)',[0.99 0.98],'color','k')
  grid(hca,'off');
  set(hca,'xticklabel',[]);

  hca=irf_panel('iEnSp');
  irf_spectrogram(hca, speciEnSp, 'log', 'donotfitcolorbarlabel');
  colormap(hca,jet)
  set(hca,'yscale','log');
  set(hca,'ytick',[1e1 1e2 1e3 1e4]);
  %scaxis(hca,[-1, 2])
  ylabel(hca,{'E_{i} (eV)'},'Interpreter','tex');

  hca=irf_panel('eEnSp');
  irf_spectrogram(hca, speceEnSp, 'log', 'donotfitcolorbarlabel');
  colormap(hca,jet)
  set(hca,'yscale','log');
  set(hca,'ytick',[1e1 1e2 1e3 1e4]);
  %scaxis(hca,[-1, 2])
  ylabel(hca,{'E_{e} (eV)'},'Interpreter','tex');

  title(h(1),strcat(['MMS', num2str(ic)]));
  irf_plot_axis_align(h);
  irf_pl_number_subplots(h,[0.99, 0.95]);
  irf_zoom(h,'x',Tint);


  %export figure as eps picture with name filename.
  filename=['MMS',num2str(ic),'_CurrentOverview_', num2str(i)];
  irf_print_fig(h,filename,'png')
end
% 2) from terminal convert to eps file without white margins
% > epstool --copy --bbox Jan11pic10.eps Jan11pic10_crop.eps
% 3) convert eps file to pdf, result is in Current_crop.pdf
% > ps2pdf -dEPSFitPage -dEPSCrop -dAutoRotatePages=/None Jan11pic10_crop.eps

