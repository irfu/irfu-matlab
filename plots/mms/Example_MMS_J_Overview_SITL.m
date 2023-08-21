% A routine that identifies strong currents from SITL data and then
%plots E, B, J, EdotJ for the specific intervals.
% Written by Elin Eriksson
%
% Loads data and calculates currents (curlometer method).
%Plots B, J (curlometer method), E, JxB electric field, and J.E as an average
%of all spacecraft values and a seperate plot with the same variables where only the E and
%B are from a specific spacecraft given by the number of ic.

%TODO: Fix resampling when data gaps exist when looking at larger
%intervals.

%% Set time interval to look in, the specific spacecraft of interest and set limits for Null method.
% taken from
ic=1; %Gives number of spacecraft where density is taken for Hall field calculations.
currentLim=500E-9;
Tint  = irf.tint('2016-01-26T20:00:00Z/2016-01-27T15:00:00Z');

disp('Loading Magnetic fields');
c_eval('B?=mms.db_get_ts(''mms?_dfg_srvy_ql'',''mms?_dfg_srvy_dmpa'',Tint);',1:4);
%Resamples according to the time line where all spacecraft has data
timeStart=max([B1.time.start.epochUnix B2.time.start.epochUnix B3.time.start.epochUnix B4.time.start.epochUnix],[],2);
timeStart=EpochUnix(timeStart);
timeStop=min([B1.time.stop.epochUnix B2.time.stop.epochUnix B3.time.stop.epochUnix B4.time.stop.epochUnix],[],2);
timeStop=EpochUnix(timeStop);
timeLog=B1.time >= timeStart & B1.time <=timeStop;
newTime=B1.time(timeLog,:);
c_eval('B? = B?.resample(newTime);',1:4);
% Spacecraft Position
disp('Loading Spacecraft Position');
R  = mms.get_data('R_gse',Tint); %Cailbrated position
if   length(R.gseR1(1,:))==4 && length(R.gseR2(1,:))==4 && length(R.gseR3(1,:))==4 && length(R.gseR4(1,:))==4
  % Assume the first column is time
  c_eval('R? =irf.ts_vec_xyz(R.time, R.gseR?(:,1:3));',1:4);
  clear R
  c_eval('R? = R?.resample(B1);',1:4);
else
  c_eval('R? =irf.ts_vec_xyz(R.time, R.gseR?);',1:4);
  clear R
  c_eval('R? = R?.resample(B1);',1:4);
end

disp('Looking for strong currents in brst data');
curlB = c_4_grad('R?','B?','curl');
j     = curlB.data./1.0e3.*1e-9./(4*pi*1e-7); %A/m^2 if B in nT and R in km
jabs=irf_abs([curlB.time.epochUnix j],1);
strongCurrents=jabs>currentLim;
if sum(strongCurrents) >= 1
  IntervalStrongCurrent=[curlB.time(strongCurrents).epochUnix jabs(strongCurrents,1)];
else
  error('No strong currents in this interval')
end


%Sorts the currents found from highest to lowest size
[~,IndexSize] = sort(IntervalStrongCurrent(:,2),'descend');
IntervalStrongCurrent=IntervalStrongCurrent(IndexSize,:); %[currentTime currentAbs]

%TODO - still some overlapping between the intervals
jTemp = [IntervalStrongCurrent ones(size(IntervalStrongCurrent,1),1)];
for ii = 1:size(jTemp,1)
  if jTemp(ii,3)
    jTemp(jTemp(:,1)<jTemp(ii,1)+30 & jTemp(:,1)>jTemp(ii,1)-30,3)=0;
    jTemp(ii,3)=1;
  end
end

IntervalStrongCurrent=IntervalStrongCurrent(logical(jTemp(:,3)),:);

%%

disp('Loads Electric fields')
% Electric field
c_eval('Efield=mms.db_get_ts(''mms?_edp_fast_ql_dce2d'',''mms?_edp_dce_xyz_dsl'',Tint);',ic);
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
  error('No fpi data available')
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
j=[curlB.time.epochUnix j];
j(:,2:4) = j(:,2:4).*1e9; %nA/m^2


%% Plot data should only be used if nulls are found. Plots are only focused around the intervals of nulls found because of above segment.
disp('Plotting figure')
h = irf_plot(5,'newfigure');

hca = irf_panel('BMMSav');
irf_plot(hca,Bfield);
ylabel(hca,{'B_{DMPA}','(nT)'},'Interpreter','tex');
irf_legend(hca,{'B_{x}','B_{y}','B_{z}'},[0.88 0.10],'Interpreter','tex')
%irf_legend(hca,'(a)',[0.99 0.98],'color','k')
grid(hca,'off');
set(hca,'xticklabel',[]);

hca = irf_panel('J');
irf_plot(hca,j);
ylabel(hca,{'J_{DMPA}','(nA m^{-2})'},'Interpreter','tex');
irf_legend(hca,{'J_{x}','J_{y}','J_{z}'},[0.88 0.10],'Interpreter','tex')
%irf_legend(hca,'(c)',[0.99 0.98],'color','k')
grid(hca,'off');
set(hca,'xticklabel',[]);

hca = irf_panel('EMMS');
irf_plot(hca,Efield);
ylabel(hca,{'E_{DSL}','(mV m^{-1})'},'Interpreter','tex');
irf_legend(hca,{'E_{x}','E_{y}','E_{z}'},[0.88 0.10],'Interpreter','tex')
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
tmarks=IntervalStrongCurrent(:,1); %Makes tmarks for all null data points.
irf_pl_mark(h,tmarks,[0.8 0.8 0.8]) %Setting the color lines to grey
irf_pl_number_subplots(h,[0.99, 0.95]);
tint_zoom=[Bfield.time.start.epochUnix Bfield.time.stop.epochUnix];
irf_zoom(h,'x',tint_zoom);

% 2) from terminal convert to eps file without white margins
% > epstool --copy --bbox Jan11pic10.eps Jan11pic10_crop.eps
% 3) convert eps file to pdf, result is in Current_crop.pdf
% > ps2pdf -dEPSFitPage -dEPSCrop -dAutoRotatePages=/None Jan11pic10_crop.eps

