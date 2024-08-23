%% Turbulence model spectra and noise       // SolO noise turb
% 20100622 version (Milans antenna parameters)
Units=irf_units;
if 1 % parameters
  f_range=[1e-3 3.9e5];
  f_noise_range=10.^(log10(f_range(1)*5):.1:log10(f_range(2)));
  f=f_noise_range';
  PE_range=[1e-19 1e-6]; % (V^/m^)/Hz
  title_text='';
  %Rsolo=0.28;UV=1;R_plasma_nobias=3e6;R_plasma_bias=0.1e6;plasma=getfield(solar_orbiter('plasma'),'perihelion');
  Rsolo=1.0;UV=1;
  T_plasma_eV=30; % plasma temperature in eV
  n = 20*1e6;      % plasma density m^-3
  title_text=[title_text 'T_p=' num2str(T_plasma_eV) ' eV, '];
  % HFA
  HFA_antenna_eff_length=2.5/2; % efficient distance between antennas
  HFA_R_plasma_nobias=14e6;HFA_R_plasma_bias=5e6;
  HFA_probe=lp.default_lprobe('THOR_HFA');
  HFA_C_antenna=HFA_probe.capacitance; % antenna capacitance in F
  title_text=[title_text 'HFAC_{ant}=' num2str(HFA_probe.capacitance/1e-12) ' pF, '];
  A_antenna=HFA_probe.Area.total; % antenna area in m2
  title_text=[title_text 'HFAA_{ant}=' num2str(HFA_probe.Area.total,'%5.2f') ' m^2, '];

  % HFA-D
  HFD_antenna_eff_length=2/2; % efficient distance between antennas % 2m-long dipole
  HFD_R_plasma_nobias=30e6;HFD_R_plasma_bias=10e6;
  HFD_probe=lp.default_lprobe('THOR_HFD');
  HFD_C_antenna=HFD_probe.capacitance; % antenna capacitance in F
  title_text=[title_text 'HFDC_{ant}=' num2str(HFD_probe.capacitance/1e-12) ' pF, '];
  A_antenna=HFD_probe.Area.total; % antenna area in m2
  title_text=[title_text 'HFDA_{ant}=' num2str(HFD_probe.Area.total,'%5.2f') ' m^2, '];

  % SDP
  SDP_antenna_eff_length=50; % efficient distance between antennas
  SDP_R_plasma_nobias=916e6;SDP_R_plasma_bias=25e6; SDP_R_plasma_bias_low=10e6;
  SDP_probe=lp.default_lprobe('THOR_SDP');
  SDP_C_antenna=SDP_probe.capacitance; % antenna capacitance in F
  title_text=[title_text 'SDPC_{ant}=' num2str(SDP_probe.capacitance/1e-12) ' pF, '];
  A_antenna=HFA_probe.Area.total; % antenna area in m2
  title_text=[title_text 'SDPA_{ant}=' num2str(SDP_probe.Area.total,'%5.2f') ' m^2, '];

end
if 1 % read Alexandrova et al. 2009 spectra
  %file reading
  [xx,yy]=textread('alexandrova2009_fig2.dat','%s%s','headerlines',3);
  for j=1:size(xx,1)
    alexandrova2009_Bspectra(j,1)=10^str2num(xx{j});
    alexandrova2009_Bspectra(j,2)=10^str2num(yy{j})*1e-6; % (V/m)^2/Hz
  end
  clear xx yy;
end
if 1 % read Sahraoui et al. 2009 spectra
  %file reading
  [xx,yy]=textread('sahraoui2009_fig4.dat','%s%s','headerlines',3);
  for j=1:size(xx,1)
    sahraoui2009_Espectra(j,1)=10^str2num(xx{j});
    sahraoui2009_Espectra(j,2)=10^str2num(yy{j})*1e-6; % (V/m)^2/Hz
  end
  clear xx yy;
end
if 1 % read Yuris spectra for Sahraoui et al. 2009 event
  %file reading
  [xx,yy]=textread('Cluster_yuri.dat','%s%s','headerlines',3);
  for j=1:size(xx,1)
    Cluster_Espectra(j,1)=10^str2num(xx{j});
    Cluster_Espectra(j,2)=10^str2num(yy{j})*1e-6; % (V/m)^2/Hz
  end
  clear xx yy;
end
if 1 % read Helios2 noise spectra for Sahraoui et al. 2009 event
  %file reading
  [xx,yy]=textread('Helios2_noise.dat','%s%s','headerlines',3);
  for j=1:size(xx,1)
    Helios2_Espectra(j,1)=10^str2num(xx{j});
    Helios2_Espectra(j,2)=10^str2num(yy{j}); % (V/m)^2/Hz
  end
  clear xx yy;
end
if 0 % read example solar wind spectra
  %file reading
  [xx,yy]=textread('THEMIS_solar_wind_example_spectra.txt','%s%s','headerlines',4);
  for j=1:size(xx,1)
    SW_example_Espectra_THEMIS(j,1)=str2num(xx{j}); % Hz
    SW_example_Espectra_THEMIS(j,2)=str2num(yy{j})*1e-6; % (V/m)^2/Hz
  end
  clear xx yy;
end
if 1 % read example solar wind spectra
  %file reading
  [xx,yy]=textread('Cluster4_Ey_20070130_0000_0120.dat','%s%s','headerlines',4);
  for j=1:size(xx,1)
    SW_example_Espectra_Cluster4(j,1)=str2num(xx{j}); % Hz
    SW_example_Espectra_Cluster4(j,2)=str2num(yy{j})*1e-6; % (V/m)^2/Hz
  end
  clear xx yy;
end
if 1 % read example solar wind spectra
  %file reading
  [xx,yy]=textread('Cluster4_EySW_20070130_0000_0120.dat','%s%s','headerlines',4);
  for j=1:size(xx,1)
    SW_example_Espectra_Cluster4SW(j,1)=str2num(xx{j}); % Hz
    SW_example_Espectra_Cluster4SW(j,2)=str2num(yy{j})*1e-6; % (V/m)^2/Hz
  end
  clear xx yy;
end

if 1 % read example solar wind spectra
  %file reading
  [xx,yy]=textread('MMS4_brst_EySW_20151228_0518_0519.dat','%s%s','headerlines',4);
  for j=1:size(xx,1)
    SW_example_Espectra_MMS4brst(j,1)=str2num(xx{j}); % Hz
    SW_example_Espectra_MMS4brst(j,2)=str2num(yy{j})*1e-6; % (V/m)^2/Hz
  end
  clear xx yy;
end
if 1 % read example solar wind spectra
  %file reading
  [xx,yy]=textread('MMS4_fast_EySW_20151228_0511_0526.dat','%s%s','headerlines',4);
  for j=1:size(xx,1)
    SW_example_Espectra_MMS4fast(j,1)=str2num(xx{j}); % Hz
    SW_example_Espectra_MMS4fast(j,2)=str2num(yy{j})*1e-6; % (V/m)^2/Hz
  end
  clear xx yy;
end
if 1 % read example solar wind spectra
  %file reading
  [xx,yy]=textread('MMS4_fast_EySC_20151228_0511_0526.dat','%s%s','headerlines',4);
  for j=1:size(xx,1)
    SWsc_example_Espectra_MMS4fast(j,1)=str2num(xx{j}); % Hz
    SWsc_example_Espectra_MMS4fast(j,2)=str2num(yy{j})*1e-6; % (V/m)^2/Hz
  end
  clear xx yy;
end
if 1 % read example magnetopause whistler spectra from LeContell 2016
  %file readind
  xx = load('whistlerspec.mat');
  MP_example_Espectra_MMS1brstMPwhi(:,1) = xx.whistlerspec.freq; % Hz
  MP_example_Espectra_MMS1brstMPwhi(:,2) = xx.whistlerspec.power*1e-6; % (V/m)^2/Hz
  clear xx ;
end
if 1 % read example magnetopause IAW spectra from Steinvall 2021
  %file readind
  xx = load('acousticspec.mat');
  MP_example_Espectra_MMS1brstMPaiw(:,1) = xx.acousticspec.freq; % Hz
  MP_example_Espectra_MMS1brstMPaiw(:,2) = xx.acousticspec.power*1e-6; % (V/m)^2/Hz
  clear xx ;
end
if 1 % read example shock spectra
  %file reading
  [xx,yy]=textread('MMS4_brst_EyShock_20151222_071324.dat','%s%s','headerlines',4);
  for j=1:size(xx,1)
    SW_example_Espectra_MMS4brstShock(j,1)=str2num(xx{j}); % Hz
    SW_example_Espectra_MMS4brstShock(j,2)=str2num(yy{j})*1e-6; % (V/m)^2/Hz
  end
  clear xx yy;
end
if 1 % read example solar wind spectra
  %file reading
  [xx,yy]=textread('MMS4_hmfe_EySW_20151228_051852.dat','%s%s','headerlines',4);
  for j=1:size(xx,1)
    SW_example_Espectra_MMS4hmfe(j,1)=str2num(xx{j}); % Hz
    SW_example_Espectra_MMS4hmfe(j,2)=str2num(yy{j})*1e-6; % (V/m)^2/Hz
  end
  clear xx yy;
end
if 1 % instrument noise calculations HFA
  if 1 % preamp parameters
    HFA_preamp_noise=150e-9; % preamplifier noise 4nV/Hz1/2
    HFA_preamp_noise_level=(HFA_preamp_noise/HFA_antenna_eff_length)^2;
    f_break=300; % transition frequency at which 1/f noise is starting
  end
  if 1 % preamplifier noise
    HFA_instr_noise=[f(1) f_break f(end)]';
    HFA_instr_noise(:,2)=2*HFA_preamp_noise_level*[(f_break/f(1))^2; 1;  1];
  end
  if 1 % photoelectron thermal noise  S=4kTZ
    T_eV=1; % photoelectron temperature
    HFA_thermal_noise_bias=[f 4*Units.e*T_eV*sqrt(HFA_R_plasma_bias^2./(1+(2*pi*f).^2*HFA_R_plasma_bias^2*HFA_C_antenna^2))/HFA_antenna_eff_length^2];
    HFA_thermal_noise_nobias=[f 4*Units.e*T_eV*sqrt(HFA_R_plasma_nobias^2./(1+(2*pi*f).^2*HFA_R_plasma_nobias^2*HFA_C_antenna^2))/HFA_antenna_eff_length^2];
  end
  if 1 % shot noise plasma and photoelectrons
    % SOFI_shot_noise_plasma XXX WE IGNORE THIS
    %nu=n/2*sqrt(8*Units.e*T_plasma_eV/pi/Units.me)*A_antenna;
    %SOFI_shot_noise_bias=[f 2*Units.e^2*nu*(R_plasma_bias^2./(1+(2*pi*f).^2*R_plasma_bias^2*C_antenna^2))/antenna_eff_length^2];
    %SOFI_shot_noise_nobias=[f 2*Units.e^2*nu*(R_plasma_nobias^2./(1+(2*pi*f).^2*R_plasma_nobias^2*C_antenna^2))/antenna_eff_length^2];
    % SOFI_shot_noise_photoelectron
    HFA_I=lp.current(HFA_probe,-1,Rsolo,UV,[]);

    HFA_shot_noise_photoelectron_bias=[f 2*Units.e*abs(HFA_I.photo)*(HFA_R_plasma_bias^2./(1+(2*pi*f).^2*HFA_R_plasma_bias^2*HFA_C_antenna^2))/HFA_antenna_eff_length^2];
    HFA_shot_noise_photoelectron_nobias=[f 2*Units.e*abs(HFA_I.photo)*(HFA_R_plasma_nobias^2./(1+(2*pi*f).^2*HFA_R_plasma_nobias^2*HFA_C_antenna^2))/HFA_antenna_eff_length^2];
  end
  if 1 % total noise = shot noise photo + thermal
    HFA_total_noise_bias=irf_add(1,HFA_shot_noise_photoelectron_bias,1,HFA_thermal_noise_bias);
    HFA_total_noise_nobias=irf_add(1,HFA_shot_noise_photoelectron_nobias,1,HFA_thermal_noise_nobias);
  end
  if 1 % bit noise
    f_sampling=25e3; % DC up to kHz
    tmunit=15.2e-6; % in V/m for gain=1
    tmrange=0.5;    % in V/m
    SOFI_bit_noise_gain1=[f_range(1) f_range(1)*10]';
    SOFI_bit_noise_gain1(:,2)=(tmunit/2)^2/(f_sampling/2)/3;
  end
end

if 1 % instrument noise calculations HFA-D
  if 1 % preamp parameters
    HFD_preamp_noise=150e-9; % preamplifier noise 4nV/Hz1/2
    HFD_preamp_noise_level=(HFD_preamp_noise/HFD_antenna_eff_length)^2;
    f_break=300; % transition frequency at which 1/f noise is starting
  end
  if 1 % preamplifier noise
    HFD_instr_noise=[f(1) f_break f(end)]';
    HFD_instr_noise(:,2)=2*HFD_preamp_noise_level*[(f_break/f(1))^2; 1;  1];
  end
  if 1 % photoelectron thermal noise  S=4kTZ
    T_eV=1; % photoelectron temperature
    HFD_thermal_noise_bias=[f 4*Units.e*T_eV*sqrt(HFD_R_plasma_bias^2./(1+(2*pi*f).^2*HFD_R_plasma_bias^2*HFD_C_antenna^2))/HFD_antenna_eff_length^2];
    HFD_thermal_noise_nobias=[f 4*Units.e*T_eV*sqrt(HFD_R_plasma_nobias^2./(1+(2*pi*f).^2*HFD_R_plasma_nobias^2*HFD_C_antenna^2))/HFD_antenna_eff_length^2];
  end
  if 1 % shot noise plasma and photoelectrons
    % SOFI_shot_noise_plasma XXX WE IGNORE THIS
    %nu=n/2*sqrt(8*Units.e*T_plasma_eV/pi/Units.me)*A_antenna;
    %SOFI_shot_noise_bias=[f 2*Units.e^2*nu*(R_plasma_bias^2./(1+(2*pi*f).^2*R_plasma_bias^2*C_antenna^2))/antenna_eff_length^2];
    %SOFI_shot_noise_nobias=[f 2*Units.e^2*nu*(R_plasma_nobias^2./(1+(2*pi*f).^2*R_plasma_nobias^2*C_antenna^2))/antenna_eff_length^2];
    % SOFI_shot_noise_photoelectron
    HFD_I=lp.current(HFD_probe,-1,Rsolo,UV,[]);

    HFD_shot_noise_photoelectron_bias=[f 2*Units.e*abs(HFD_I.photo)*(HFD_R_plasma_bias^2./(1+(2*pi*f).^2*HFD_R_plasma_bias^2*HFD_C_antenna^2))/HFD_antenna_eff_length^2];
    HFD_shot_noise_photoelectron_nobias=[f 2*Units.e*abs(HFD_I.photo)*(HFD_R_plasma_nobias^2./(1+(2*pi*f).^2*HFD_R_plasma_nobias^2*HFD_C_antenna^2))/HFD_antenna_eff_length^2];
  end
  if 1 % total noise = shot noise photo + thermal
    HFD_total_noise_bias=irf_add(1,HFD_shot_noise_photoelectron_bias,1,HFD_thermal_noise_bias);
    HFD_total_noise_nobias=irf_add(1,HFD_shot_noise_photoelectron_nobias,1,HFD_thermal_noise_nobias);
  end
  if 1 % bit noise
    f_sampling=25e3; % DC up to kHz
    tmunit=15.2e-6; % in V/m for gain=1
    tmrange=0.5;    % in V/m
    SOFI_bit_noise_gain1=[f_range(1) f_range(1)*10]';
    SOFI_bit_noise_gain1(:,2)=(tmunit/2)^2/(f_sampling/2)/3;
  end
end

if 1 % instrument noise calculations SDP
  if 1 % preamp parameters
    SDP_preamp_noise=150e-9; % preamplifier noise 70nV/Hz1/2
    SDP_preamp_noise_level=(SDP_preamp_noise/SDP_antenna_eff_length)^2;
    f_break=300; % transition frequency at which 1/f noise is starting
  end
  if 1 % preamplifier noise
    SDP_instr_noise=[f(1) f_break f(end)]';
    SDP_instr_noise(:,2)=2*SDP_preamp_noise_level*[(f_break/f(1))^2; 1;  1];
  end
  if 1 % photoelectron thermal noise  S=4kTZ
    T_eV=1; % photoelectron temperature
    SDP_thermal_noise_bias=[f 4*Units.e*T_eV*sqrt(SDP_R_plasma_bias^2./(1+(2*pi*f).^2*SDP_R_plasma_bias^2*SDP_C_antenna^2))/SDP_antenna_eff_length^2];
    SDP_thermal_noise_bias_low=[f 4*Units.e*T_eV*sqrt(SDP_R_plasma_bias_low^2./(1+(2*pi*f).^2*SDP_R_plasma_bias_low^2*SDP_C_antenna^2))/SDP_antenna_eff_length^2];
  end
  if 1 % shot noise plasma and photoelectrons
    % SOFI_shot_noise_plasma XXX WE IGNORE THIS
    %nu=n/2*sqrt(8*Units.e*T_plasma_eV/pi/Units.me)*A_antenna;
    %SOFI_shot_noise_bias=[f 2*Units.e^2*nu*(R_plasma_bias^2./(1+(2*pi*f).^2*R_plasma_bias^2*C_antenna^2))/antenna_eff_length^2];
    %SOFI_shot_noise_nobias=[f 2*Units.e^2*nu*(R_plasma_nobias^2./(1+(2*pi*f).^2*R_plasma_nobias^2*C_antenna^2))/antenna_eff_length^2];
    % SOFI_shot_noise_photoelectron
    SDP_I=lp.current(SDP_probe,-1,Rsolo,UV,[]);

    SDP_shot_noise_photoelectron_bias=[f 2*Units.e*abs(SDP_I.photo)*(SDP_R_plasma_bias^2./(1+(2*pi*f).^2*SDP_R_plasma_bias^2*SDP_C_antenna^2))/SDP_antenna_eff_length^2];
    SDP_shot_noise_photoelectron_bias_low=[f 2*Units.e*abs(SDP_I.photo)*(SDP_R_plasma_bias_low^2./(1+(2*pi*f).^2*SDP_R_plasma_bias_low^2*SDP_C_antenna^2))/SDP_antenna_eff_length^2];
  end
  if 1 % total noise = shot noise photo + thermal
    SDP_total_noise_bias=irf_add(1,SDP_shot_noise_photoelectron_bias,1,SDP_thermal_noise_bias);
    SDP_total_noise_bias_low=irf_add(1,SDP_shot_noise_photoelectron_bias_low,1,SDP_thermal_noise_bias_low);
  end
  if 1 % bit noise
    f_sampling=25e3; % DC up to kHz
    tmunit=15.2e-6; % in V/m for gain=1
    tmrange=0.5;    % in V/m
    SOFI_bit_noise_gain1=[f_range(1) f_range(1)*10]';
    SOFI_bit_noise_gain1(:,2)=(tmunit/2)^2/(f_sampling/2)/3;
  end
end

if 1 % calculate spectra 1AU
  %%
  Vf=400;T=30;B=5;nn=3; % R=1AU
  VA=20*B/sqrt(nn);
  clear SP
  SP.power(1)=1.6e2;% power of signal at first frequency
  SP.Epower(1)=SP.power(1)*1e-12*VA^2;
  SP.EpowerVf(1)=SP.power(1)*1e-12*Vf^2;
  SP.R_AU=1; % distance in AU

  f_dop_ion=Vf/(2*pi*(100*sqrt(2*T)/B));
  f_dop_electron=Vf/(2*pi*(sqrt(10*T)/B));
  SP.R_RS=1/0.00465; % distance in AU
  SP.f=[f_range(1) f_dop_ion f_dop_electron f_dop_electron*10];
  SP.slopes=[-1.6 -2.5 -2.5];
  for i=2:length(SP.f)
    SP.power(i)=SP.power(i-1)*10^(log10(SP.f(i)/SP.f(i-1))*SP.slopes(i-1));
    if i==2
      SP.Epower(i)=SP.power(i)*1e-12*VA^2;   % (V/m)^2
    else
      SP.Epower(i)=SP.Epower(i-1)*10^(log10(SP.f(i)/SP.f(i-1))*SP.slopes(i-1)/2);
    end
    SP.EpowerVf(i)=SP.power(i)*1e-12*Vf^2; % (V/m)^2
  end
  SP.f=[SP.f 1e4 1e5];
  SP.Epower(end+1) = 1e-16; % Expected QTN intensity
  SP.Epower(end+1) = 1e-16;
  SP1AU=SP;
  %figure(2),loglog(SP.f(1:4),SP.power),grid on
end

if 1 % EMC requirements
  HFA_EMC.f =      [10    1e2   1e3   1e4   1e5   2e5];
  HFA_EMC.Epower = [2e-11 5e-13 5e-13 3e-14 3e-14 5e-13]; % (V/m)^2/Hz

  SDP_EMC.f =      [0.1 30    1e2   1e3   1e4   1e5   2e5];
  SDP_EMC.Epower =     [1e-10 1e-15 1e-15 1e-16 1e-17 1e-17 1e-15];
  SDP_EMC_OLD.f =      [0.1 1 10    1e2   1e3   1e4   1e5   2e5];
  SDP_EMC_OLD.Epower = [5e-8 1e-10 2e-13 1e-15 1e-16 1e-17 1e-17 1e-15];
end

if 1 % initialize figure - HFA
  figure(12);clf
  h=irf_plot(1);
  set(h,'position',[0.15 0.1 0.75 0.75]);
  set(gcf,'defaultLineLineWidth',2);
  set(gcf,'PaperUnits','centimeters')
  xSize = 12; ySize = 13;
  xLeft = (21-xSize)/2; yTop = (30-ySize)/2;
  set(gcf,'PaperPosition',[xLeft yTop xSize ySize])
  set(gcf,'Position',[10 10 xSize*50 ySize*50])
end
if 1 % electric field example spectra
  % MMS MP IAWs
  %loglog(MP_example_Espectra_MMS1brstMPaiw(:,1), MP_example_Espectra_MMS1brstMPaiw(:,2),'color',[0.1 0.5 0.0],'linewidth',1);
  %text(129,8e-9,'MMS-burst-iaw-MP','fontsize',10,'color',[0.1 0.5 0.0],'units','data','horizontalalignment','left','verticalalignment','bottom');
  %hca=h(1);
  %hold(hca,'on');

  % MMS shock
  loglog(SW_example_Espectra_MMS4brstShock(:,1), SW_example_Espectra_MMS4brstShock(:,2),'color',[0.8 0.5 0.0],'linewidth',1);
  text(129,9e-9,'MMS-shocks-IAWs','fontsize',10,'color',[0.8 0.5 0.0],'units','data','horizontalalignment','left','verticalalignment','bottom');
  hca=h(1);
  hold(hca,'on');

  % MMS MP Whistlers
  loglog(MP_example_Espectra_MMS1brstMPwhi(:,1), MP_example_Espectra_MMS1brstMPwhi(:,2),'color',[0.1 0.5 0.0],'linewidth',1);
  text(1290,7e-10,'MMS-MP-Whistlers','fontsize',10,'color',[0.1 0.5 0.0],'units','data','horizontalalignment','left','verticalalignment','bottom');
  hca=h(1);
  hold(hca,'on');



  % MMS hmfe SW Langmuir waves
  loglog(SW_example_Espectra_MMS4hmfe(2048:end,1), SW_example_Espectra_MMS4hmfe(2048:end,2),'color',[0.5 0.5 0],'linewidth',1);
  text(1e4,2e-10,'MMS-SW-LWs','fontsize',10,'color',[0.5 0.5 0],'units','data','horizontalalignment','left','verticalalignment','bottom');

  % Cluster example solar wind spectrum from CalibrationReport
  %loglog(SW_example_Espectra_Cluster4(:,1), SW_example_Espectra_Cluster4(:,2),'color',[.49 .18 .56],'linewidth',1);
  %text(3,3e-10,'Cluster-sw','fontsize',10,'color',[.49 .18 .56],'units','data','horizontalalignment','left','verticalalignment','bottom');
end
if 1 % electric field plot
  hca=h(1); hold(hca,'on');set(hca,'XScale','log','YScale','log')
  set(hca,'xlim',[10 1e5])
  set(hca,'ylim',[2e-16 9e-7])
  set(hca,'xtick',10.^[log10(f_range(1)):1:log10(f_range(2))]),
  set(hca,'ytick',10.^[log10(PE_range(1)):2:log10(PE_range(2))]),
  grid(hca,'on');
  set(hca,'xminorgrid','off');
  set(hca,'yminorgrid','off');

  ylabel(hca,'S_E [(V/m)^2/Hz]');
  xlabel(hca,'frequency [Hz]');

  %title(hca,'Expected electric field spectra and HFA noise levels')
  title(hca,'Expected HFA noise levels')
end

if 1 % plot electric field noises

  if 0 % HFA
    HFA_color = [0.8 0.5 0.5];
    loglog(HFA_total_noise_bias,HFA_total_noise_bias(:,2),'color',HFA_color);hold(hca,'on');
    loglog(HFA_total_noise_nobias(:,1),HFA_total_noise_nobias(:,2),'color',HFA_color,'linestyle',':');
    text(HFA_total_noise_nobias(1,1)*200,HFA_total_noise_nobias(1,2),'HFA nobias noise','fontsize',10,'color',HFA_color,'units','data','horizontalalignment','left','verticalalignment','bottom');
    text(HFA_total_noise_bias(1,1)*5000,HFA_total_noise_bias(1,2)*1.2,' plasma+photoelectron\newline fluctuation level','fontsize',10,'color',HFA_color,'units','data','horizontalalignment','left','verticalalignment','bottom');
    loglog(1*[1 5],5e-17*10*[1 1],'color',HFA_color,'linestyle',':');
    loglog(1*[1 5],2e-17*10*[1 1],'color',HFA_color,'linestyle','-');
    text(8,5e-16,['unbiased HFA, R=' num2str(HFA_R_plasma_nobias/1e6,3) 'M\Omega'],'horizontalalignment','left','verticalalignment','middle','color','k');
    text(8,2e-16,['biased HFA, R=' num2str(HFA_R_plasma_bias/1e6,3) 'M\Omega'],'horizontalalignment','left','verticalalignment','middle','color','k');
    loglog(HFA_instr_noise(:,1), HFA_instr_noise(:,2),'color',[0.3 0.3 0.3]); hold(hca,'on');
    text(HFA_instr_noise(2,1),HFA_instr_noise(end,2)*0.9,'preamp noise','fontsize',10,'color','k','units','data','horizontalalignment','left','verticalalignment','top');
    loglog(hca,HFA_EMC.f,HFA_EMC.Epower,'o-','markersize',10,'color',[0 0.8 0]);
    text(0.98,0.55,'EMC req','fontsize',12,'fontweight','demi','color',[0 0.8 0],'units','normalized','horizontalalignment','right','parent',hca);
  end

  if 1 % HFA-D
    HFD_color = [0.0 0.0 0.0];
    if 0 % nobias
      loglog(HFD_total_noise_nobias(:,1),HFD_total_noise_nobias(:,2),'color',HFD_color,'linestyle',':');
      text(HFD_total_noise_nobias(1,1)*200,HFD_total_noise_nobias(1,2),'HFA-D nobias noise','fontsize',10,'color',HFD_color,'units','data','horizontalalignment','left','verticalalignment','bottom');
      loglog(1*[1 5],5e-16*10*[1 1],'color',HFD_color,'linestyle',':');
      text(8,5e-15,['unbiased HFA-D, R=' num2str(HFD_R_plasma_nobias/1e6,3) 'M\Omega'],'horizontalalignment','left','verticalalignment','middle','color','k');
    end
    if 1 % bias
      loglog(HFD_total_noise_bias,HFD_total_noise_bias(:,2),'color',HFD_color);hold(hca,'on');
      text(HFD_total_noise_bias(1,1)*5000,HFD_total_noise_bias(1,2)*1.2,' HFA-D bias noise','fontsize',10,'color',HFD_color,'units','data','horizontalalignment','left','verticalalignment','bottom');
      loglog(1*[1 5],2e-16*10*[1 1],'color',HFD_color,'linestyle','-');
      text(12,2e-15,['biased HFA-D, R=' num2str(HFD_R_plasma_bias/1e6,3) 'M\Omega'],'horizontalalignment','left','verticalalignment','middle','color','k');
    end
    if 0 % preamp
      loglog(HFD_instr_noise(:,1), HFD_instr_noise(:,2),'color',[0.5 0.5 0.5]); hold(hca,'on');
      text(HFD_instr_noise(2,1),HFD_instr_noise(end,2)*0.9,'preamp noise D','fontsize',10,'color','k','units','data','horizontalalignment','left','verticalalignment','top');
    end
    %loglog(hca,HFA_EMC.f,HFA_EMC.Epower,'o-','markersize',10,'color',[0 0.8 0]);
    %text(0.98,0.55,'EMC req','fontsize',12,'fontweight','demi','color',[0 0.8 0],'units','normalized','horizontalalignment','right','parent',hca);
  end

  irf_legend(hca,['ne=' num2str(n(1)/1e6,3) 'cc, Te=' num2str(T_plasma_eV(1),3) 'eV'],[0.98 0.02])
end
irf_legend(0,['PO EFI-HFA noise ' char(datetime("now","Format","uuuu-MM-dd HH:mm:ss"))],[0,0.001],'interpreter','none','color',[0.5 0.5 0.5])
irf_print_fig(['HFA_noise_' char(datetime("now","Format","uuuuMMdd"))],'png')

if 1 % initialize figure - SDP
  figure(13);clf
  h=irf_plot(1);
  set(h,'position',[0.15 0.1 0.75 0.75]);
  set(gcf,'defaultLineLineWidth',2);
  set(gcf,'PaperUnits','centimeters')
  xSize = 12; ySize = 13;
  xLeft = (21-xSize)/2; yTop = (30-ySize)/2;
  set(gcf,'PaperPosition',[xLeft yTop xSize ySize])
  set(gcf,'Position',[10 10 xSize*50 ySize*50])
end
if 1 % electric field plot
  hca=h(1);
  loglog(hca,SP1AU.f,SP1AU.Epower,'b--','markersize',20);
  hold(hca,'on');
  %loglog(hca,SP1AU.f(1:4),SP1AU.EpowerVf,'r.-','markersize',20);

  %set(hca,'xlim',f_range)
  set(hca,'xlim',[9e-2 5e5])
  set(hca,'ylim',PE_range)
  set(hca,'xtick',10.^[log10(f_range(1)):1:log10(f_range(2))]),
  set(hca,'ytick',10.^[log10(PE_range(1)):2:log10(PE_range(2))]),
  grid(hca,'on');
  set(hca,'xminorgrid','off');
  set(hca,'yminorgrid','off');

  ylabel(hca,'S_E [(V/m)^2/Hz]');
  xlabel(hca,'frequency [Hz]');

  text(0.97,0.85,'S_{E} - Expected SW-frame spectrum','fontsize',12,'fontweight','demi','color','b','units','normalized','horizontalalignment','right','parent',hca);
  %text(0.97,0.8,'S_{ExB} - Expected spectrum of V_{SW}xB','fontsize',12,'fontweight','demi','color','r','units','normalized','horizontalalignment','right','parent',hca);

  %text(0.97,0.85,'spectra at R=0.3 AU','fontsize',12,'fontweight','demi','color','r','units','normalized','horizontalalignment','right','parent',hca);
  title(hca,['Predicted E-field spectra and SDP noise levels in SW\newline' ...
    'S_{E} =V_A^2*S_B (inertial rng), V_A=' num2str(VA,'%.1f') ' km/s, S_{E} \sim k^2 S_B (kinetic rng)  \newline'...
    'S_B based on Sahraoui et al, 2009'])
  %'S_{VxB} =V_{SW}^2*S_B, V_{SW}=' num2str(Vf,'%.1f') ' km/s (S_B based on Sahraoui et al, 2009)'])
end
if 1 % electric field example spectra
  %loglog(SW_example_Espectra_THEMIS(:,1), SW_example_Espectra_THEMIS(:,2),'color',[0.5 0.5 0.5],'linewidth',1);
  %text(29,1.2e-13,'THEMIS','fontsize',10,'color',[0.5 0.5 0.5],'units','data','horizontalalignment','left','verticalalignment','bottom');

  % Cluster example solar wind spectrum from CalibrationReport
  loglog(SW_example_Espectra_Cluster4(:,1), SW_example_Espectra_Cluster4(:,2),'color',[.49 .18 .56],'linewidth',1);
  text(3,3e-10,'Cluster','fontsize',10,'color',[.49 .18 .56],'units','data','horizontalalignment','left','verticalalignment','bottom');


  % Cluster example solar wind spectrum from CalibrationReport
  %loglog(SW_example_Espectra_MMS4fast(:,1), SW_example_Espectra_MMS4fast(:,2),'color',[0.5 0.5 0.5],'linewidth',1);
  %text(.3,3e-9,'MMS-fast','fontsize',10,'color',[0.5 0.5 0.5],'units','data','horizontalalignment','left','verticalalignment','bottom');

  % MMS Solar wind example fast data
  loglog(SWsc_example_Espectra_MMS4fast(:,1), SWsc_example_Espectra_MMS4fast(:,2),'color',[0.5 0.5 0.5],'linewidth',1);
  text(20,1e-12,'MMS-fast','fontsize',10,'color',[0.5 0.5 0.5],'units','data','horizontalalignment','left','verticalalignment','bottom');

  % MMS example solar wind burst data
  loglog(SW_example_Espectra_MMS4brst(:,1), SW_example_Espectra_MMS4brst(:,2),'color',[0 0 0],'linewidth',1);
  text(129,1.2e-13,'MMS-burst','fontsize',10,'color',[0 0 0],'units','data','horizontalalignment','left','verticalalignment','bottom');

  % Cluster solar wind frame
  %loglog(SW_example_Espectra_Cluster4SW(:,1), SW_example_Espectra_Cluster4SW(:,2),'color',[0.5 0.5 0.5],'linewidth',1);
  %text(3,3e-10,'Cluster-plasma frame','fontsize',10,'color',[0.5 0.5 0.5],'units','data','horizontalalignment','left','verticalalignment','bottom');


  % the power is based on electric field Cluster spectra from Sahraoui et al 2009.
  %loglog(sahraoui2009_Espectra(:,1),sahraoui2009_Espectra(:,2),'color',[0.5 0.5 0.5],'linewidth',1);
  %text(2,3e-4,'Cluster','fontsize',10,'color',[0.5 0.5 0.5],'units','data','horizontalalignment','left','verticalalignment','bottom');

  % electric field Cluster spectra for interval Sahraoui et al 2009.
  %loglog(Cluster_Espectra(:,1),Cluster_Espectra(:,2),'color',[0.5 0.5 0.5],'linewidth',1);
  %text(3,3e-10,'Cluster','fontsize',10,'color',[0.5 0.5 0.5],'units','data','horizontalalignment','left','verticalalignment','bottom');
end
if 1 % plot electric field noises
  if 1
    SDP_color = [0.8 0.5 0.0];
    loglog(SDP_total_noise_bias,SDP_total_noise_bias(:,2),'color',SDP_color);
    loglog(SDP_total_noise_bias_low(:,1),SDP_total_noise_bias_low(:,2),'color',SDP_color,'linestyle',':');
    text(SDP_total_noise_bias(1,1)*20,SDP_total_noise_bias(1,2)*.001,' plasma+photoelectron\newline fluctuation level','fontsize',10,'color',[0.8 0.5 0.0],'units','data','horizontalalignment','left','verticalalignment','bottom');
    loglog(100*f_range(1)*[1 5],PE_range(1)*3*[1 1],'color',SDP_color,'linestyle','-');
    text(100*f_range(1)*6,PE_range(1)*3,['biased SDP, R=' num2str(SDP_R_plasma_bias/1e6,3) 'M\Omega'],'horizontalalignment','left','verticalalignment','middle','color','k');
    loglog(100*f_range(1)*[1 5],PE_range(1)*15*[1 1],'color',SDP_color,'linestyle',':');
    text(100*f_range(1)*6,PE_range(1)*15,['biased SDP, R=' num2str(SDP_R_plasma_bias_low/1e6,3) 'M\Omega'],'horizontalalignment','left','verticalalignment','middle','color','k');

    loglog(SDP_instr_noise(:,1), SDP_instr_noise(:,2),'k');
    text(SDP_instr_noise(2,1)*2,SDP_instr_noise(end,2)*0.9,'preamp noise','fontsize',10,'color','k','units','data','horizontalalignment','left','verticalalignment','top');
    %loglog(hca,SDP_EMC_OLD.f,SDP_EMC_OLD.Epower,'o-','markersize',10,'color',[0 0.7 0.7]);
    %text(0.97,0.65,'OLD SDP EMC req','fontsize',12,'fontweight','demi','color',[0 0.7 0.7],'units','normalized','horizontalalignment','right','parent',hca);
    loglog(hca,SDP_EMC.f(1:3),SDP_EMC.Epower(1:3),'--','color',[0 0.5 0]);
    %loglog(hca,SDP_EMC.f(2:end),SDP_EMC.Epower(2:end),'o-','markersize',10,'color',[0 0.5 0]);
    loglog(hca,SDP_EMC_OLD.f,SDP_EMC_OLD.Epower,'o-','markersize',10,'color',[0 0.8 0]);
    text(0.97,0.45,'EMC req','fontsize',12,'fontweight','demi','color',[0 0.8 0],'units','normalized','horizontalalignment','right','parent',hca);
    text(0.97,0.40,'goal - -','fontsize',12,'fontweight','demi','color',[0 0.5 0],'units','normalized','horizontalalignment','right','parent',hca);

  end
  %irf_legend(hca,['ne=' num2str(n(1)/1e6,3) 'cc, Te=' num2str(T_plasma_eV(1),3) 'eV'],[0.98 0.02])
end
irf_legend(0,['THOR EFI-SDP noise ' char(datetime("now","Format","uuuu-MM-dd HH:mm:ss"))],[0,0.001],'interpreter','none','color',[0.5 0.5 0.5])
irf_print_fig(['SDP_noise_' char(datetime("now","Format","uuuuMMdd"))],'png')


%% HFA in dB uV/m
if 1
  EMCR = HFA_EMC; ant = 'HFA';
  EMCR.ff =      [EMCR.f, EMCR.f];
  [EMCR.ff,idx] = sort(EMCR.ff); EMCR.ff([1 end]) = [];
  EMCR.EEpower = [EMCR.Epower, EMCR.Epower];
  EMCR.EEpower = EMCR.EEpower(idx); EMCR.EEpower([1 end]) = [];
  EMCR.bw =         [.1 1 10 100 1e3];
  EMCR.bw = sort([EMCR.bw, EMCR.bw]);
  EMCR.Efield = sqrt(EMCR.EEpower.*EMCR.bw);
  HFA_EMC = EMCR;

  if 1 % initialize figure - HFA
    figure(14);clf
    h=irf_plot(1);
    set(h,'position',[0.15 0.1 0.75 0.75]);
    set(gcf,'defaultLineLineWidth',2);
    set(gcf,'PaperUnits','centimeters')
    xSize = 16; ySize = 13;
    xLeft = (21-xSize)/2; yTop = (30-ySize)/2;
    set(gcf,'PaperPosition',[xLeft yTop xSize ySize])
    set(gcf,'Position',[10 10 xSize*50 ySize*50])
  end

  semilogx(EMCR.ff,10*log10(1e6*EMCR.Efield),'r'), hold on
  semilogx(EMCR.ff,30+10*log10(1e6*EMCR.Efield),'b')
  set(gca,'xLim',EMCR.ff([1 end]))
  grid on, set(gca,'xminorgrid','off'), set(gca,'yminorgrid','off')
  title(['EFI-' ant ' requirement'])
  ylabel('E-field [ dB ({\mu}V/m) ]'), xlabel('frequency [Hz]')
  legend('Broadband noise limit','Limit for stable narrowband emissions')
  irf_legend(0,['THOR EFI-' ant ' EMC req ' char(datetime("now","Format","uuuu-MM-dd HH:mm:ss"))],[0,0.001],'interpreter','none','color',[0.5 0.5 0.5])
  irf_print_fig([ant '_req_dBuVm_' char(datetime("now","Format","uuuuMMdd"))],'png')
end

%% SDP in dB uV/m
if 1
  EMCR = SDP_EMC_OLD; ant = 'SDP';
  EMCR.ff =      [EMCR.f, EMCR.f];
  [EMCR.ff,idx] = sort(EMCR.ff); EMCR.ff([1 end]) = [];
  EMCR.EEpower = [EMCR.Epower, EMCR.Epower];
  EMCR.EEpower = EMCR.EEpower(idx); EMCR.EEpower([1 end]) = [];
  EMCR.bw =         [.001 .01 .1 1 10 100 1e3];
  EMCR.bw = sort([EMCR.bw, EMCR.bw]);
  EMCR.Efield = sqrt(EMCR.EEpower.*EMCR.bw);
  SDP_EMC_OLD = EMCR;

  if 1 % initialize figure - SDP
    figure(15);clf
    h=irf_plot(1);
    set(h,'position',[0.15 0.1 0.75 0.75]);
    set(gcf,'defaultLineLineWidth',2);
    set(gcf,'PaperUnits','centimeters')
    xSize = 16; ySize = 13;
    xLeft = (21-xSize)/2; yTop = (30-ySize)/2;
    set(gcf,'PaperPosition',[xLeft yTop xSize ySize])
    set(gcf,'Position',[10 10 xSize*50 ySize*50])
  end

  semilogx(EMCR.ff,10*log10(1e6*EMCR.Efield),'r'), hold on
  semilogx(EMCR.ff,30+10*log10(1e6*EMCR.Efield),'b')
  set(gca,'xLim',EMCR.ff([1 end]))

  if 1 % GOAL
    SDP_EMC.f =      [0.1 1 10 30    1e2   1e3   1e4   1e5   2e5];
    SDP_EMC.Epower =     [1e-10 1e-12 1e-14 1e-15 1e-15 1e-16 1e-17 1e-17 1e-15];
    EMCR = SDP_EMC;
    EMCR.ff =      [EMCR.f, EMCR.f];
    [EMCR.ff,idx] = sort(EMCR.ff); EMCR.ff([1 end]) = [];
    EMCR.EEpower = [EMCR.Epower, EMCR.Epower];
    EMCR.EEpower = EMCR.EEpower(idx); EMCR.EEpower([1 end]) = [];
    EMCR.bw =         [.001 .01 .1 .1 1 10 100 1e3];
    EMCR.bw = sort([EMCR.bw, EMCR.bw]);
    EMCR.Efield = sqrt(EMCR.EEpower.*EMCR.bw);
    SDP_EMC = EMCR;
    semilogx(SDP_EMC.ff,10*log10(1e6*SDP_EMC.Efield),'r--')
    semilogx(SDP_EMC.ff,30+10*log10(1e6*SDP_EMC.Efield),'b--')
  end


  grid on, set(gca,'xminorgrid','off'), set(gca,'yminorgrid','off')
  title(['EFI-' ant ' requirement'])
  ylabel('E-field [ dB ({\mu}V/m) ]'), xlabel('frequency [Hz]')
  legend('Broadband noise limit (goal --)','Limit for stable narrowband emissions  (goal --)')
  irf_legend(0,['THOR EFI-' ant ' EMC req ' char(datetime("now","Format","uuuu-MM-dd HH:mm:ss"))],[0,0.001],'interpreter','none','color',[0.5 0.5 0.5])
  irf_print_fig([ant '_req_dBuVm_' char(datetime("now","Format","uuuuMMdd"))],'png')
end

%% compare to SolO, JUICE

% We scale to HFA location,
% e.g. 5 m distance, similar to SolO and JUICE/PWI sensors

if 1 % initialize figure - comparison
  figure(16);clf
  h=irf_plot(1);
  set(h,'position',[0.15 0.1 0.75 0.75]);
  set(gcf,'defaultLineLineWidth',2);
  set(gcf,'PaperUnits','centimeters')
  xSize = 22; ySize = 13;
  xLeft = (21-xSize)/2; yTop = (30-ySize)/2;
  set(gcf,'PaperPosition',[xLeft yTop xSize ySize])
  set(gcf,'Position',[10 10 xSize*50 ySize*50])
end

SDP_Scale = (50/5)^3;

semilogx(HFA_EMC.ff,10*log10(1e6*HFA_EMC.Efield),'m'), hold on
semilogx(SDP_EMC_OLD.ff,10*log10(1e6*SDP_Scale*SDP_EMC_OLD.Efield),'b')

SDP_EMC.f =      [0.1 1 10 30    1e2   1e3   1e4   1e5   2e5];
SDP_EMC.Epower =     [1e-10 1e-12 1e-14 1e-15 1e-15 1e-16 1e-17 1e-17 1e-15];
EMCR = SDP_EMC;
EMCR.ff =      [EMCR.f, EMCR.f];
[EMCR.ff,idx] = sort(EMCR.ff); EMCR.ff([1 end]) = [];
EMCR.EEpower = [EMCR.Epower, EMCR.Epower];
EMCR.EEpower = EMCR.EEpower(idx); EMCR.EEpower([1 end]) = [];
EMCR.bw =         [.001 .01 .1 .1 1 10 100 1e3];
EMCR.bw = sort([EMCR.bw, EMCR.bw]);
EMCR.Efield = sqrt(EMCR.EEpower.*EMCR.bw);
SDP_EMC = EMCR;
semilogx(SDP_EMC.ff,10*log10(1e6*SDP_Scale*SDP_EMC.Efield),'b--')

% JUICE RPW EID-B (i3.5, Jan 2017), EIDB-S00310 - boadband noise at sensor position
JUICE_PWI_EMC.ff           = [1  10 10 220 220 1e3 1e3 1e4 1e4 1e5 1e5 16e5];
JUICE_PWI_EMC.Efield_dBuVm = [25 5  15 -10 -10 -10 -5  -5  0   0   20  20];
semilogx(JUICE_PWI_EMC.ff,JUICE_PWI_EMC.Efield_dBuVm,'k')

% JUICE EID-A (i2r7, July 2016), EIDA-R003706 - radiated emissions for
% space-exposed and transition equipment at 1m
%JUICE_EMC.ff           = [30 1e7];
%JUICE_EMC.Efield_dBuVm = [-10 -10];
%JUICE_Scale = 10*log10((5/1)^3); % scale to HFA location
%emilogx(JUICE_EMC.ff,JUICE_EMC.Efield_dBuVm+JUICE_Scale,'k--')

% SolO EID-A (i2,r8, Oct 2011), EIDA R-707 - boadband noise at preamp
% location
SOLO_Scale = 10*log10((5/1)^3); % scale to HFA location
SOLO_EMC.ff           = [1    10  100  1e3 1e4 1e5  1e6];
SOLO_EMC.Efield_dBnVm = [15.8 3.8 -8.2 1.8 1.8 -2.9 27.6];
semilogx(SOLO_EMC.ff,SOLO_EMC.Efield_dBnVm-30+SOLO_Scale,'r') % nV->uV

ylabel('E-field [ dB ({\mu}V/m) ]'), xlabel('frequency [Hz]')
legend('HFA','SDP','SDP goal','JUICE RPW EID-B','SOLO EID-A')

set(gca,'XLim',[.1 2e5],'YLim',[-20 38])
grid on, set(gca,'xminorgrid','off'), set(gca,'yminorgrid','off')
title(['EMC reqs on broadband noise (scaled to HFA location)'])

irf_legend(0,['THOR/JUICE/SOLO EMC req ' char(datetime("now","Format","uuuu-MM-dd HH:mm:ss"))],[0,0.001],'interpreter','none','color',[0.5 0.5 0.5])
irf_print_fig(['THOR_EMC_vs_JUICE_SOLO_' char(datetime("now","Format","uuuuMMdd"))],'png')

%% compare to SolO, JUICE in dB (uV/m)^2 /Hz

% We scale to HFA location,
% e.g. 5 m distance, similar to SolO and JUICE/PWI sensors

if 1 % initialize figure - comparison
  figure(17);clf
  h=irf_plot(1);
  set(h,'position',[0.15 0.1 0.75 0.75]);
  set(gcf,'defaultLineLineWidth',2);
  set(gcf,'PaperUnits','centimeters')
  xSize = 22; ySize = 13;
  xLeft = (21-xSize)/2; yTop = (30-ySize)/2;
  set(gcf,'PaperPosition',[xLeft yTop xSize ySize])
  set(gcf,'Position',[10 10 xSize*50 ySize*50])
end

SDP_Scale = (50/1)^3; HFA_Scale = (5/1)^3;

loglog(HFA_EMC.ff,sqrt(HFA_Scale^2*HFA_EMC.EEpower),'m'), hold on
loglog(SDP_EMC_OLD.ff,sqrt(SDP_Scale^2*SDP_EMC_OLD.EEpower),'b')
loglog(SDP_EMC.ff,sqrt(SDP_Scale^2*SDP_EMC.EEpower),'b--')

% JUICE RPWI EID-B (i3.5, Jan 2017), EIDB-S00310 - boadband noise at sensor position
JUICE_PWI_Scale = (3/1)^3;
JUICE_PWI_EMC.ff           = [1  10 10 220 220 1e3 1e3 1e4 1e4 1e5 1e5 16e5];
JUICE_PWI_EMC.bw           = [1   1 10  10  10  10  30  30 100 100 10e3  10e3]; % New corrected numbers from Lennart, 2017-01-19
JUICE_PWI_EMC.Efield_dBuVm = [25 5  15 -10 -10 -10 -5  -5  0   0   20  20]; % New corrected numbers from Lennart, 2017-01-19
JUICE_PWI_EMC.EEpower = 1e-12*10.^(2*(JUICE_PWI_EMC.Efield_dBuVm + 20*log10(JUICE_PWI_Scale))/20)./JUICE_PWI_EMC.bw; % (V/m)^2/Hz
loglog(JUICE_PWI_EMC.ff,sqrt(JUICE_PWI_EMC.EEpower),'k')

% JUICE EID-A (i2r7, July 2016), EIDA-R003706 - radiated emissions for
% space-exposed and transition equipment at 1m
%JUICE_EMC.ff           = [30 1e7];
%JUICE_EMC.Efield_dBuVm = [-10 -10];
%JUICE_Scale = 10*log10((5/1)^3); % scale to HFA location
%emilogx(JUICE_EMC.ff,JUICE_EMC.Efield_dBuVm+JUICE_Scale,'k--')

% SolO EID-A (i5,r0, Mar 2015), EIDA-5145/EIDA R-707 - boadband noise at preamp
% location
%SOLO_Scale = 10*log10((5/1)^3); % scale to HFA location
SOLO_Scale = 0; % It is not appropriate to scale SOLO numbers, as they are give at the premp location
SOLO_EMC.ff           = [1    10  100  1e3 1e4 1e5  1e6];
SOLO_EMC.Efield_dBnVm = [15.8 3.8 -8.2 1.8 1.8 -2.9 27.6];
SOLO_EMC.EEpower = 1e-18*10.^(2*(SOLO_EMC.Efield_dBnVm+SOLO_Scale)/20); % (V/m)^2/Hz
loglog(SOLO_EMC.ff,sqrt(SOLO_EMC.EEpower),'r')

ylabel('noise [ (V/m)/Hz^{1/2} ]'), xlabel('frequency [Hz]')
legend('THOR-HFA','THOR-SDP','THOR-SDP goal','JUICE RPW EID-B','SOLO EID-A')

set(gca,'XLim',[.1 2e5],'YLim',[-20 38])
grid on, set(gca,'xminorgrid','off'), set(gca,'yminorgrid','off')
title(['EMC reqs on broadband noise (scaled to 1m from the S/C)'])

irf_legend(0,['THOR/JUICE/SOLO EMC req ' char(datetime("now","Format","uuuu-MM-dd HH:mm:ss"))],[0,0.001],'interpreter','none','color',[0.5 0.5 0.5])
irf_print_fig(['THOR_EMC_vs_JUICE_SOLO_' char(datetime("now","Format","uuuuMMdd"))],'png')
