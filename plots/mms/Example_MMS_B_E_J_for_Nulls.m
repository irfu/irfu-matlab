% Loads data and calculates nulls that has a BoxLim of 70km and a Current
% Lim of 500E-9 [nA/m^2 if B is in nT and R is in km].
%Plots of B, J (curlometer method), E, JxB electric field, and J.E. 
% Written by Elin Eriksson

%% Set time interval to look in and set the spacecraft that density will be
% taken from
ic=1; %Gives number of spacecraft where density is taken for Hall field calculations.
Tint  = irf.tint('2015-12-10T23:00:00Z/2015-12-11T23:59:59Z');
boxLim=70;
currentLim=500E-9;
%% Load magnetic field and spacecraft positional data
mms.db_init('local_file_db','/data/mms'); %Point towards the folder with the mms data

% Magnetic Field
disp('Loading Magnetic fields');
c_eval('B?=mms.db_get_ts(''mms?_dfg_srvy_ql'',''mms?_dfg_srvy_dmpa'',Tint);',1:4);
%c_eval('B?=removerepeatpnts(B?)',[1:4]); %Removes data points with too small a time difference between them which buggs the resample functions
c_eval('B? = B?.resample(B1);',2:4);

% Spacecraft Position
disp('Loading Spacecraft Position');
R  = mms.get_data('R_gse',Tint);%Cailbrated position 
c_eval('R? =irf.ts_vec_xyz(R.time, R.gseR?);',1:4);
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
%% Looks for Nulls and if they are found makes average calculations and Loads electric fields, density and calculates J and JxB with curlometer method
disp('Looking for Nulls');
% Assuming GSE and DMPA are the same coordinate system.
Nulls=c_4_null(R1,R2,R3,R4,B1,B2,B3,B4,'boxLim',boxLim,'strong',currentLim);

if isempty(Nulls.t)
    disp('No Nulls are found')
else
% Possibility to save Null data
if 0
    filename='YearMonthDay_Nulls';
    save(filename,'Nulls');
end
%% Average calculations. Loads electric fields, density and calculates J and JxB with curlometer method
disp('Calculates average Bfield')
Bav = (B1.data+B2.data+B3.data+B4.data)/4;
Bav=[B1.time.epochUnix double(Bav)];
disp('Loads Electric fields')
% Electric field
c_eval('E?_orig=mms.db_get_ts(''mms?_edp_fast_ql_dce2d'',''mms?_edp_dce_xyz_dsl'',Tint);');
c_eval('E? = E?_orig.resample(E1_orig);',1:4);
c_eval('E? = E?.resample(B1);',1:4); %Resampling to B-field due to the index searching later to make smaller filesizes.
% Calibrating E-field (TO-DO)
%c_eval('E?.data(find(abs(E?.data) > 100)) = NaN;',1:4); %Remove some questionable fields
%c_eval('[E?,~]=irf_edb(E?,B?,15,''E.B=0'');',1:4); % Removes some wake fields
% Average E-field for EdotJ
disp('Calculates average E-fields')
Eav = (E1.data+E2.data+E3.data+E4.data)/4;
Eav=[E1.time.epochUnix double(Eav)];
disp(['Loads density from spacecraft',num2str(ic)])
%Density from spacecraft ic for Hall field calculation
c_eval('ni=mms.db_get_ts(''mms?_fpi_fast_sitl'',''mms?_fpi_DISnumberDensity'',Tint);',ic);
ni = TSeries(ni.time,ni.data,'to',1);
ni = ni.resample(B1);
ni=[ni.time.epochUnix double(ni.data)];

disp('Calculates current')
% Assuming GSE and DMPA are the same coordinate system calculates j and jXB.
[j,divB,B,jxB,divTshear,divPb] = c_4_j('R?','B?');
divovercurl = divB;
divovercurl.data = abs(divovercurl.data)./j.abs.data;
divovercurl=[B1.time.epochUnix double(divovercurl.data)];

j=[j.time.epochUnix double(j.data)];
j(:,2:4) = j(:,2:4).*1e9;
jxB=[jxB.time.epochUnix double(jxB.data)];
jxB(:,2:4) = jxB(:,2:4)./[ni(:,2) ni(:,2) ni(:,2)]; 
jxB(:,2:4) = jxB(:,2:4)./1.6e-19./1000; %Convert to (mV/m)
jxB((irf_abs(jxB(:,2:4),1) > 100),1:4) = NaN; % Remove some questionable fields
disp('Calculates EdotJ')
% Calculates EdotJ
EdotJ = dot(Eav(:,2:4),j(:,2:4),2)./1000; %J (nA/m^2), E (mV/m), E.J (nW/m^3)
EdotJ = [Eav(:,1) EdotJ];
end
%% Only picks out a smaller time interval around the nulls (to keep filesize small)  
%Possibility to pick out interesting time interval to keep the file size small
if true
if isempty(Nulls.t)
    error('No Nulls are found')
else
index=find(Bav(:,1)>(Nulls.t(1,1)-1) & Bav(:,1)<Nulls.t(end,1)+1);
Bav=Bav(index,:);
Eav=Eav(index,:);
ni=ni(index,:);
j=j(index,:);
jxB=jxB(index,:);
EdotJ=EdotJ(index,:);
divovercurl=divovercurl(index,:);
end
end
%% Plot data should only be used if nulls are found. Plots are only focused around the intervals of nulls found because of above segment.
if isempty(Nulls.t)
    error('Could not plot because no nulls have been found') 
else
h = irf_plot(6,'newfigure');

hca = irf_panel('BMMSav');
irf_plot(hca,Bav);
ylabel(hca,{'B_{DMPA}','(nT)'},'Interpreter','tex');
irf_legend(hca,{'B_{x}','B_{y}','B_{z}'},[0.88 0.10])
%irf_legend(hca,'(a)',[0.99 0.98],'color','k')
grid(h(1),'off');

hca = irf_panel('J');
irf_plot(hca,j);
ylabel(hca,{'J_{DMPA}','(nA m^{-2})'},'Interpreter','tex');
irf_legend(hca,{'J_{x}','J_{y}','J_{z}'},[0.88 0.10])
%irf_legend(hca,'(c)',[0.99 0.98],'color','k')
grid(hca,'off');

hca = irf_panel('divovercurl');
irf_plot(hca,divovercurl);
ylabel(hca,{'|\nabla . B| ','|\nabla \times B|'},'Interpreter','tex');
%irf_legend(hca,'(e)',[0.99 0.98],'color','k')
grid(hca,'off');

hca = irf_panel('EMMSav');
irf_plot(hca,Eav);
ylabel(hca,{'E_{DSL}','(mV m^{-1})'},'Interpreter','tex');
irf_legend(hca,{'E_{x}','E_{y}','E_{z}'},[0.88 0.10])
%irf_legend(hca,'(b)',[0.99 0.98],'color','k')
grid(hca,'off');

hca = irf_panel('jxB');
irf_plot(hca,jxB);
ylabel(hca,{'J \times B/n_{e} q_{e}','(mV m^{-1})'},'Interpreter','tex');
%irf_legend(hca,'(f)',[0.99 0.98],'color','k')
grid(hca,'off');

hca = irf_panel('jdotE');
irf_plot(hca,EdotJ);
ylabel(hca,{'E . J','(nW m^{-3})'},'Interpreter','tex');
%irf_legend(hca,'(g)',[0.99 0.98],'color','k')
grid(hca,'off');

title(h(1),strcat('MMS averaged - Current density and fields'));
tmarks=Nulls.t; %Makes tmarks for all null data points.
irf_pl_mark(h,tmarks,[0.8 0.8 0.8]) %Setting the color lines to grey
irf_pl_number_subplots(h,[0.05, 0.9]); 
irf_plot_axis_align(1,h)
Tint=[Bav(1,1) Bav(end,1)];
irf_zoom(h,'x',Tint);
end

% 2) from terminal convert to eps file without white margins
% > epstool --copy --bbox 15DecNullBoring.eps 15DecNullBoring_crop.eps
% 3) convert eps file to pdf, result is in Current_crop.pdf
% > ps2pdf -dEPSFitPage -dEPSCrop -dAutoRotatePages=/None 15DecNullBoring_crop.eps

