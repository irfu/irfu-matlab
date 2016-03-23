%% Turbulence model spectra and noise       // SolO noise turb
% 20100622 version (Milans antenna parameters)
Units=irf_units;
if 1, % parameters
    f_range=[1e-3 9.9e5];
    f_noise_range=10.^(log10(f_range(1)*50):.1:log10(f_range(2)));
    f=f_noise_range';
    PE_range=[1e-19 1e-6]; % (V^/m^)/Hz
    title_text='';
    %Rsolo=0.28;UV=1;R_plasma_nobias=3e6;R_plasma_bias=0.1e6;plasma=getfield(solar_orbiter('plasma'),'perihelion');
    Rsolo=1.0;UV=1;
    T_plasma_eV=10; % plasma temperature in eV
    n = 2*1e6;      % plasma density m^-3
    title_text=[title_text 'T_p=' num2str(T_plasma_eV) ' eV, '];
    % HFA
    HFA_antenna_eff_length=2.5/2; % efficient distance between antennas
    HFA_R_plasma_nobias=467e6;HFA_R_plasma_bias=3e6;
    HFA_probe=lp.default_lprobe('THOR_HFA');
    HFA_C_antenna=HFA_probe.capacitance; % antenna capacitance in F
    title_text=[title_text 'HFAC_{ant}=' num2str(HFA_probe.capacitance/1e-12) ' pF, '];
    A_antenna=HFA_probe.Area.total; % antenna area in m2
    title_text=[title_text 'HFAA_{ant}=' num2str(HFA_probe.Area.total,'%5.2f') ' m^2, '];
    % SDP
    SDP_antenna_eff_length=100; % efficient distance between antennas
    SDP_R_plasma_nobias=916e6;SDP_R_plasma_bias=21e6;
    SDP_probe=lp.default_lprobe('THOR_SDP');
    SDP_C_antenna=SDP_probe.capacitance; % antenna capacitance in F
    title_text=[title_text 'SDPC_{ant}=' num2str(SDP_probe.capacitance/1e-12) ' pF, '];
    A_antenna=HFA_probe.Area.total; % antenna area in m2
    title_text=[title_text 'SDPA_{ant}=' num2str(SDP_probe.Area.total,'%5.2f') ' m^2, '];
    
end
if 1, % read Alexandrova et al. 2009 spectra
    %file reading
    [xx,yy]=textread('alexandrova2009_fig2.dat','%s%s','headerlines',3);
    for j=1:size(xx,1),
        alexandrova2009_Bspectra(j,1)=10^str2num(xx{j});
        alexandrova2009_Bspectra(j,2)=10^str2num(yy{j})*1e-6; % (V/m)^2/Hz
    end
    clear xx yy;
end
if 1, % read Sahraoui et al. 2009 spectra
    %file reading
    [xx,yy]=textread('sahraoui2009_fig4.dat','%s%s','headerlines',3);
    for j=1:size(xx,1),
        sahraoui2009_Espectra(j,1)=10^str2num(xx{j}); 
        sahraoui2009_Espectra(j,2)=10^str2num(yy{j})*1e-6; % (V/m)^2/Hz
    end
    clear xx yy;
end
if 1, % read Yuris spectra for Sahraoui et al. 2009 event
    %file reading
    [xx,yy]=textread('Cluster_yuri.dat','%s%s','headerlines',3);
    for j=1:size(xx,1),
        Cluster_Espectra(j,1)=10^str2num(xx{j}); 
        Cluster_Espectra(j,2)=10^str2num(yy{j})*1e-6; % (V/m)^2/Hz
    end
    clear xx yy;
end
if 1, % read Helios2 noise spectra for Sahraoui et al. 2009 event
    %file reading
    [xx,yy]=textread('Helios2_noise.dat','%s%s','headerlines',3);
    for j=1:size(xx,1),
        Helios2_Espectra(j,1)=10^str2num(xx{j}); 
        Helios2_Espectra(j,2)=10^str2num(yy{j}); % (V/m)^2/Hz
    end
    clear xx yy;
end
if 1, % read example solar wind spectra
    %file reading
    [xx,yy]=textread('THEMIS_solar_wind_example_spectra.txt','%s%s','headerlines',4);
    for j=1:size(xx,1),
        SW_example_Espectra_THEMIS(j,1)=str2num(xx{j}); % Hz
        SW_example_Espectra_THEMIS(j,2)=str2num(yy{j})*1e-6; % (V/m)^2/Hz
    end
    clear xx yy;
end
if 1, % instrument noise calculations HFA
    if 1, % preamp parameters
        preamp_noise=10e-9; % preamplifier noise 4nV/Hz1/2
        preamp_noise_level=(preamp_noise/HFA_antenna_eff_length)^2;
        f_break=400; % transition frequency at which 1/f noise is starting
    end
    if 1, % preamplifier noise
        HFA_instr_noise=[f(1) f_break f(end)]';
        HFA_instr_noise(:,2)=2*preamp_noise_level*[f_break/f(1); 1;  1];
    end
    if 1, % photoelectron thermal noise  S=4kTZ
        T_eV=1; % photoelectron temperature
        HFA_thermal_noise_bias=[f 4*Units.e*T_eV*sqrt(HFA_R_plasma_bias^2./(1+(2*pi*f).^2*HFA_R_plasma_bias^2*HFA_C_antenna^2))/HFA_antenna_eff_length^2];
        HFA_thermal_noise_nobias=[f 4*Units.e*T_eV*sqrt(HFA_R_plasma_nobias^2./(1+(2*pi*f).^2*HFA_R_plasma_nobias^2*HFA_C_antenna^2))/HFA_antenna_eff_length^2];
    end
    if 1, % shot noise plasma and photoelectrons
        % SOFI_shot_noise_plasma XXX WE IGNORE THIS
        %nu=n/2*sqrt(8*Units.e*T_plasma_eV/pi/Units.me)*A_antenna;
        %SOFI_shot_noise_bias=[f 2*Units.e^2*nu*(R_plasma_bias^2./(1+(2*pi*f).^2*R_plasma_bias^2*C_antenna^2))/antenna_eff_length^2];
        %SOFI_shot_noise_nobias=[f 2*Units.e^2*nu*(R_plasma_nobias^2./(1+(2*pi*f).^2*R_plasma_nobias^2*C_antenna^2))/antenna_eff_length^2];
        % SOFI_shot_noise_photoelectron
        HFA_I=lp.current(HFA_probe,-1,Rsolo,UV,[]);
        
        HFA_shot_noise_photoelectron_bias=[f 2*Units.e*abs(HFA_I.photo)*(HFA_R_plasma_bias^2./(1+(2*pi*f).^2*HFA_R_plasma_bias^2*HFA_C_antenna^2))/HFA_antenna_eff_length^2];
        HFA_shot_noise_photoelectron_nobias=[f 2*Units.e*abs(HFA_I.photo)*(HFA_R_plasma_nobias^2./(1+(2*pi*f).^2*HFA_R_plasma_nobias^2*HFA_C_antenna^2))/HFA_antenna_eff_length^2];
    end
    if 1, % total noise = shot noise photo + thermal
        HFA_total_noise_bias=irf_add(1,HFA_shot_noise_photoelectron_bias,1,HFA_thermal_noise_bias);
        HFA_total_noise_nobias=irf_add(1,HFA_shot_noise_photoelectron_nobias,1,HFA_thermal_noise_nobias);
    end
    if 1, % bit noise
        f_sampling=25e3; % DC up to kHz
        tmunit=15.2e-6; % in V/m for gain=1
        tmrange=0.5;    % in V/m
        SOFI_bit_noise_gain1=[f_range(1) f_range(1)*10]';
        SOFI_bit_noise_gain1(:,2)=(tmunit/2)^2/(f_sampling/2)/3;
    end
end

if 1, % instrument noise calculations SDP
    if 1, % preamp parameters
        preamp_noise=10e-9; % preamplifier noise 4nV/Hz1/2
        preamp_noise_level=(preamp_noise/HFA_antenna_eff_length)^2;
        f_break=400; % transition frequency at which 1/f noise is starting
    end
    if 1, % preamplifier noise
        SDP_instr_noise=[f(1) f_break f(end)]';
        SDP_instr_noise(:,2)=2*preamp_noise_level*[f_break/f(1); 1;  1];
    end
    if 1, % photoelectron thermal noise  S=4kTZ
        T_eV=1; % photoelectron temperature
        SDP_thermal_noise_bias=[f 4*Units.e*T_eV*sqrt(SDP_R_plasma_bias^2./(1+(2*pi*f).^2*SDP_R_plasma_bias^2*SDP_C_antenna^2))/SDP_antenna_eff_length^2];
        SDP_thermal_noise_nobias=[f 4*Units.e*T_eV*sqrt(SDP_R_plasma_nobias^2./(1+(2*pi*f).^2*SDP_R_plasma_nobias^2*SDP_C_antenna^2))/SDP_antenna_eff_length^2];
    end
    if 1, % shot noise plasma and photoelectrons
        % SOFI_shot_noise_plasma XXX WE IGNORE THIS
        %nu=n/2*sqrt(8*Units.e*T_plasma_eV/pi/Units.me)*A_antenna;
        %SOFI_shot_noise_bias=[f 2*Units.e^2*nu*(R_plasma_bias^2./(1+(2*pi*f).^2*R_plasma_bias^2*C_antenna^2))/antenna_eff_length^2];
        %SOFI_shot_noise_nobias=[f 2*Units.e^2*nu*(R_plasma_nobias^2./(1+(2*pi*f).^2*R_plasma_nobias^2*C_antenna^2))/antenna_eff_length^2];
        % SOFI_shot_noise_photoelectron
        SDP_I=lp.current(SDP_probe,-1,Rsolo,UV,[]);
        
        SDP_shot_noise_photoelectron_bias=[f 2*Units.e*abs(SDP_I.photo)*(SDP_R_plasma_bias^2./(1+(2*pi*f).^2*SDP_R_plasma_bias^2*SDP_C_antenna^2))/SDP_antenna_eff_length^2];
        SDP_shot_noise_photoelectron_nobias=[f 2*Units.e*abs(SDP_I.photo)*(SDP_R_plasma_nobias^2./(1+(2*pi*f).^2*SDP_R_plasma_nobias^2*SDP_C_antenna^2))/SDP_antenna_eff_length^2];
    end
    if 1, % total noise = shot noise photo + thermal
        SDP_total_noise_bias=irf_add(1,SDP_shot_noise_photoelectron_bias,1,SDP_thermal_noise_bias);
        SDP_total_noise_nobias=irf_add(1,SDP_shot_noise_photoelectron_nobias,1,SDP_thermal_noise_nobias);
    end
    if 1, % bit noise
        f_sampling=25e3; % DC up to kHz
        tmunit=15.2e-6; % in V/m for gain=1
        tmrange=0.5;    % in V/m
        SOFI_bit_noise_gain1=[f_range(1) f_range(1)*10]';
        SOFI_bit_noise_gain1(:,2)=(tmunit/2)^2/(f_sampling/2)/3;
    end
end

if 1, % calculate spectra 1AU
    Vf=630;T=60;B=15;nn=20; % R=1AU
    VA=20*B/sqrt(nn);
    SP.power(1)=2e5*1e-6;% power of signal at first frequency
    SP.Epower(1)=SP.power(1)*1e-6*VA^2;
    SP.EpowerVf(1)=SP.power(1)*1e-6*Vf^2;
    SP.R_AU=1; % distance in AU
    
    f_dop_ion=Vf/(2*pi*(100*sqrt(2*T)/B))*0.5;
    f_dop_electron=Vf/(2*pi*(sqrt(10*T)/B))*0.5;
    SP.R_RS=1/0.00465; % distance in AU
    SP.f=[f_range(1) 1e-3 f_dop_ion f_dop_electron f_dop_electron*10];
    SP.slopes=[-1 -1.6 -2.8 -4];
    for i=2:length(SP.f),
        SP.power(i)=SP.power(i-1)*10^(log10(SP.f(i)/SP.f(i-1))*SP.slopes(i-1));
        SP.Epower(i)=SP.power(i)*1e-6*VA^2;
        SP.EpowerVf(i)=SP.power(i)*1e-6*Vf^2;
    end
    SP1AU=SP;
end
if 1, % calculate spectra R=.3AU
    Vf=385;T=4;B=40;nn=25; % R=.3AU standard conditions
    VA=20*B/sqrt(nn);
    SP.power(1)=3e6*1e-6; % power of signal at first frequency
    SP.Epower(1)=SP.power(1)*1e-6*Vf^2;
    SP.R_AU=.3; % distance in AU
    f_dop_ion=Vf/(2*pi*(100*sqrt(2*T)/B))*0.5;
    f_dop_electron=Vf/(2*pi*(sqrt(10*T)/B))*0.5;
    SP.R_RS=1/0.00465; % distance in AU
    SP.f=[f_range(1) 6e-3 f_dop_ion f_dop_electron f_dop_electron*10];
    SP.slopes=[-1 -1.6 -2.8 -4];
    for i=2:length(SP.f),
        SP.power(i)=SP.power(i-1)*10^(log10(SP.f(i)/SP.f(i-1))*SP.slopes(i-1));
        SP.Epower(i)=SP.power(i)*1e-6*VA^2;
    end
    SP03AU=SP;
end

if 1, % initialize figure
    figure(11);clf
    h=irf_plot(1);
    set(h,'position',[0.15 0.1 0.75 0.75]);
    set(gcf,'defaultLineLineWidth',2);
    set(gcf,'PaperUnits','centimeters')
    xSize = 12; ySize = 13;
    xLeft = (21-xSize)/2; yTop = (30-ySize)/2;
    set(gcf,'PaperPosition',[xLeft yTop xSize ySize])
    set(gcf,'Position',[10 10 xSize*50 ySize*50])
end
if 1, % electric field plot
    hca=h(1);
    loglog(hca,SP1AU.f,SP1AU.Epower,'b.-','markersize',20);
    hold(hca,'on');
    loglog(hca,SP03AU.f,SP03AU.Epower,'r.-','markersize',20);
    
    set(hca,'xlim',f_range)
    set(hca,'ylim',PE_range)
    set(hca,'xtick',10.^[log10(f_range(1)):1:log10(f_range(2))]),
    set(hca,'ytick',10.^[log10(PE_range(1)):2:log10(PE_range(2))]),
    grid(hca,'on');
    set(hca,'xminorgrid','off');
    set(hca,'yminorgrid','off');
    
    ylabel(hca,'S_E [(V/m)^2/Hz]');
    xlabel(hca,'frequency [Hz]');
    
    text(0.97,0.8,'spectra at R=1 AU','fontsize',12,'fontweight','demi','color','b','units','normalized','horizontalalignment','right','parent',hca);
    text(0.97,0.85,'spectra at R=0.3 AU','fontsize',12,'fontweight','demi','color','r','units','normalized','horizontalalignment','right','parent',hca);
    title(hca,['Predicted electric field spectra and noise levels in solar wind \newline' ...
        'Noise levels calculated for distance ' num2str(Rsolo,3) 'AU \newline' ...
        'L_{eff}= ' num2str(HFA_antenna_eff_length) 'm. S_E =V_A^2*S_B (S_B empirical), V_A=' num2str(VA) 'km/s. '])
end
if 1, % electric field example spectra
    loglog(SW_example_Espectra_THEMIS(:,1), SW_example_Espectra_THEMIS(:,2),'color',[0.5 0.5 0.5],'linewidth',1);
    text(10,3e-12,'THEMIS','fontsize',10,'color',[0.5 0.5 0.5],'units','data','horizontalalignment','left','verticalalignment','bottom');
    
    % the power is based on electric field Cluster spectra from Sahraoui et al 2009.
    %loglog(sahraoui2009_Espectra(:,1),sahraoui2009_Espectra(:,2),'color',[0.5 0.5 0.5],'linewidth',1);
    %text(2,3e-4,'Cluster','fontsize',10,'color',[0.5 0.5 0.5],'units','data','horizontalalignment','left','verticalalignment','bottom');
    
    % electric field Cluster spectra for interval Sahraoui et al 2009.
    loglog(Cluster_Espectra(:,1),Cluster_Espectra(:,2),'color',[0.5 0.5 0.5],'linewidth',1);
    text(3,3e-10,'Cluster','fontsize',10,'color',[0.5 0.5 0.5],'units','data','horizontalalignment','left','verticalalignment','bottom');
    
    loglog(Helios2_Espectra(:,1),Helios2_Espectra(:,2),'color',[0.5 0.5 0.5],'linewidth',1);
    text(80,3e-9,'Helios2','fontsize',10,'color',[0.5 0.5 0.5],'units','data','horizontalalignment','left','verticalalignment','bottom');
    
end
if 1, % plot electric field noises
    
  if 1, % HFA
    loglog(HFA_total_noise_bias,HFA_total_noise_bias(:,2),'color',[0.8 0.0 0.0]);
    loglog(HFA_total_noise_nobias(:,1),HFA_total_noise_nobias(:,2),'color',[0.8 0.0 0.0],'linestyle',':');
    text(HFA_total_noise_bias(1,1)*1.5,HFA_total_noise_bias(1,2),'total noise','fontsize',10,'color',[0.8 0.0 0.0],'units','data','horizontalalignment','left','verticalalignment','bottom');
  end
  if 1
    loglog(SDP_total_noise_bias,SDP_total_noise_bias(:,2),'color',[0.8 0.5 0.0]);
    loglog(SDP_total_noise_nobias(:,1),SDP_total_noise_nobias(:,2),'color',[0.8 0.5 0.0],'linestyle',':');
    text(SDP_total_noise_bias(1,1)*1.5,SDP_total_noise_bias(1,2),'total noise','fontsize',10,'color',[0.8 0.5 0.0],'units','data','horizontalalignment','left','verticalalignment','bottom');
  end
    if 0
    loglog(thermal_noise_bias(:,1), thermal_noise_bias(:,2),'color',[0.5 0.5 0]);
    loglog(thermal_noise_nobias(:,1), thermal_noise_nobias(:,2),'color',[0.5 0.5 0],'linestyle',':');
    text(thermal_noise_bias(1,1)*200,thermal_noise_bias(1,2),'thermal noise 1eV','fontsize',10,'color',[0.5 0.5 0],'units','data','verticalalignment','bottom');
    irf_legend(['A_{antenna}=' num2str(A_antenna,3) 'm^2, C_{antenna}=' num2str(C_antenna*1e12,3) 'pF'],[0.98 0.98]);
    
    loglog(SOFI_instr_noise(:,1), SOFI_instr_noise(:,2),'k');
    text(SOFI_instr_noise(2,1),SOFI_instr_noise(end,2),'preamp-noise','fontsize',10,'color','k','units','data','horizontalalignment','left','verticalalignment','top');
    
    loglog(SOFI_bit_noise_gain1(:,1), SOFI_bit_noise_gain1(:,2),'color',[0 0.5 0]);
    text(SOFI_bit_noise_gain1(1,1)*1.5,SOFI_bit_noise_gain1(1,2),...
        ['bit noise gain=1\newline 1tm=' num2str(tmunit*1e6,3) '\muV/m\newlinerange \pm' num2str(tmrange,3) 'V/m\newline '],'fontsize',9,'color',[0 0.5 0],'units','data','verticalalignment','middle');
    gain=5;
    loglog(SOFI_bit_noise_gain1(:,1), SOFI_bit_noise_gain1(:,2)/gain^2,'color',[0 0.5 0]);
    text(SOFI_bit_noise_gain1(1,1)*1.5,SOFI_bit_noise_gain1(1,2)/gain^2,...
        ['bit noise gain=5\newline 1tm=' num2str(tmunit/gain*1e6,3) '\muV/m\newlinerange \pm' num2str(tmrange/gain,3) 'V/m\newline '],'fontsize',9,'color',[0 0.5 0],'units','data','verticalalignment','middle');
    gain=1/15;
    loglog(SOFI_bit_noise_gain1(:,1), SOFI_bit_noise_gain1(:,2)/gain^2,'color',[0 0.5 0]);
    text(SOFI_bit_noise_gain1(1,1)*1.5,SOFI_bit_noise_gain1(1,2)/gain^2,...
        ['bit noise gain=1/15\newline 1tm=' num2str(tmunit/gain*1e6,3) '\muV/m\newlinerange \pm' num2str(tmrange/gain,3) 'V/m\newline '],'fontsize',9,'color',[0 0.5 0],'units','data','verticalalignment','middle');
    
    loglog(SOFI_shot_noise_bias(:,1), SOFI_shot_noise_bias(:,2),'color',[0 0.3 0.3]);
    loglog(SOFI_shot_noise_nobias(:,1), SOFI_shot_noise_nobias(:,2),'color',[0 0.3 0.3],'linestyle',':');
    text(SOFI_shot_noise_bias(1,1)*1.5,SOFI_shot_noise_bias(1,2),'shot noise plasma','fontsize',10,'color',[0 0.3 0.3],'units','data','horizontalalignment','left','verticalalignment','bottom');
    
    loglog(SOFI_shot_noise_photoelectron_bias(:,1), SOFI_shot_noise_photoelectron_bias(:,2),'color',[0 0.5 0.5]);
    loglog(SOFI_shot_noise_photoelectron_nobias(:,1), SOFI_shot_noise_photoelectron_nobias(:,2),'color',[0 0.5 0.5],'linestyle',':');
    text(SOFI_shot_noise_photoelectron_bias(1,1)*1.5,SOFI_shot_noise_photoelectron_bias(1,2),'shot noise photo','fontsize',10,'color',[0 0.5 0.5],'units','data','horizontalalignment','left','verticalalignment','bottom');
    end
    loglog(f_range(1)*[1 10],PE_range(1)*3*[1 1],':');
    loglog(f_range(1)*[1 10],PE_range(1)*10*[1 1],'-');
    text(f_range(1)*10,PE_range(1)*3,['unbiased probe, R=' num2str(HFA_R_plasma_nobias/1e6,3) 'M\Omega'],'horizontalalignment','left','verticalalignment','middle','color','k');
    text(f_range(1)*10,PE_range(1)*10,['biased probe, R=' num2str(HFA_R_plasma_bias/1e6,3) 'M\Omega'],'horizontalalignment','left','verticalalignment','middle','color','k');
    irf_legend(hca,['ne=' num2str(n(1)/1e6,3) 'cc, Te=' num2str(T_plasma_eV(1),3) 'eV'],[0.98 0.02])
end
irf_legend(0,['SolO noise turb ' datestr(now,31)],[0,0.001],'interpreter','none','color',[0.5 0.5 0.5])
