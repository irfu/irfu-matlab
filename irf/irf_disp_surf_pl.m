function irf_disp_surf_pl(kc_x_max,kc_z_max,wfinal,extraparam,surfchoice,colorchoice,axisnorm,w_pe,m_i)
%IRF_DISP_SURF_PL   Plot the cold plasma dispersion surfaces
%
%  IRF_DISP_SURF_PL(K_PERP_MAX,K_PAR_MAX,W,EXTRAPAR,SURFCHOICE,COLCHOICE,AXISNORM)
%  plots the dispersion surfaces calculated by IRF_DISP_SURFCALC.
%
%  This function is essential for IRF_DISP_SURF.
%
%  K_PERP_MAX = max value of k_perpendicular*c/w_c
%  K_PAR_MAX = max value of k_perpendicular*c/w_c
%  W = the dispersion surfaces from IRF_DISP_SURFCALC.
%  EXTRAPAR = the extra parameters from IRF_DISP_SURFCALC
%  SURFCHOICE = vector explaining which of the surfaces to plot
%  COLCHOICE = number telling how the surface will be painted
%  AXISNORM = how axes in the plots are normalized, 1 is electron constants
%  and 2 is ion constants
%

% By Anders Tjulin, last update 4/2-2003.

% First clear the plot

cla

% Make vectors of the wave numbers

kc_z=linspace(0,kc_z_max,35);
kc_x=linspace(0,kc_x_max,35);

% Remember the view settings

[az,el]=view;

% for backwards compatability
if nargin<7; axisnorm = 1; w_pe = nan; m_i = nan; end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Set some variables according to the colorchoice
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
colorMap = 'parula'; % default colorbar
if colorchoice==2
  colorstring='<---- Electrostatic  log10(cB/E)   Electromagnetic ---->';
  colorlimits=[-2,2];
  colorMap = 'bluered';
elseif colorchoice==3
  colorstring='<---- Transversal                 Longitudinal ---->';
  colorlimits=[0,1];
elseif colorchoice==4
  colorstring='<---- Perpendicular                   Parallel ---->';
  colorlimits=[0,1];
elseif colorchoice==5
  colorstring='<---- Perpendicular                   Parallel ---->';
  colorlimits=[0,1];
elseif colorchoice==6
  colorstring='v_g/c';
  colorlimits=[0,1];
elseif colorchoice==7
  colorstring='<---- Left handed              Right handed ---->';
  colorlimits=[-1,1];
  colorMap = 'bluered';
elseif colorchoice==8
  colorstring='<---- Most energy in B    log10(eps0 E^2/(B^2/mu0))  Most energy in E ---->';
  colorlimits=[-3,3];
  colorMap = 'bluered';
elseif colorchoice==9
  colorstring='<---- Most energy in ions     log10(We/Wi)      Most energy in electrons ---->';
  colorlimits=[-3,3];
  colorMap = 'bluered';
elseif colorchoice==10
  colorstring='<---- Fast ions     log10(Ve/Vi)      Faster electrons ---->';
  colorlimits=[-1,1];
  colorMap = 'bluered';
elseif colorchoice==11
  colorstring='<---- Most energy in fields     log10(Wp/Wf)      Most energy in particles ---->';
  colorlimits=[-3,3];
  colorMap = 'bluered';
elseif colorchoice==12
  colorstring='<---- Left handed              Right handed ---->';
  colorlimits=[-1,1];
  colorMap = 'bluered';
elseif colorchoice==13
  colorstring='log10(v_{phase}/v_A)';
  colorlimits=[-2,2];
  colorMap = 'bluered';
elseif colorchoice==14
  colorstring='log10(v_{e,par}/v_{e,perp})';
  colorlimits=[-2,2];
  colorMap = 'bluered';
elseif colorchoice==15
  colorstring='log10(v_{i,par}/v_{i,perp})';
  colorlimits=[-2,2];
  colorMap = 'bluered';
elseif colorchoice==16
  colorstring='<---- Most energy in fields     log10(W_e/W_f)      Most energy in electrons ---->';
  colorlimits=[-2,2];
  colorMap = 'bluered';
elseif colorchoice==17
  colorstring='log10(dn_e/dn_i)';
  colorlimits=[-2,2];
  colorMap = 'bluered';
elseif colorchoice==18
  colorstring='log10[(dn_e/n)/(dB/B)]';
  colorlimits=[-2,2];
  colorMap = 'bluered';
elseif colorchoice==19
  colorstring='log10[(dn_i/n)/(dB/B)]';
  colorlimits=[-2,2];
  colorMap = 'bluered';
elseif colorchoice==20
  colorstring='log10[(dn_e/n)/(dBpar/B)]';
  colorlimits=[-2,2];
  colorMap = 'bluered';
elseif colorchoice==21
  colorstring='log10[(dn_i/n)/(dBpar/B)]';
  colorlimits=[-2,2];
  colorMap = 'bluered';
elseif colorchoice==22
  colorstring='log10[dn_e/(k.E eps0/e)]';
  colorlimits=[-2,2];
  colorMap = 'bluered';
elseif colorchoice==23
  colorstring='S_{||}/|S|';
  colorlimits=[0,1];
else
  colorchoice=1;
  colorlimits=[0,1];
  colorstring=('This colorscale is not implemented yet');
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Create matrices that the normal surf command can handle
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

w1(:,:)=wfinal(10,:,:);color1(:,:)=extraparam(colorchoice,10,:,:);
w2(:,:)=wfinal(9,:,:);color2(:,:)=extraparam(colorchoice,9,:,:);
w3(:,:)=wfinal(8,:,:);color3(:,:)=extraparam(colorchoice,8,:,:);
w4(:,:)=wfinal(7,:,:);color4(:,:)=extraparam(colorchoice,7,:,:);
w5(:,:)=wfinal(6,:,:);color5(:,:)=extraparam(colorchoice,6,:,:);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Convert kc_x, kc_z, w1,...w5 to ion units if requested
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if axisnorm == 2
  % convert k to k*di = k*Va/wci
  wfac = m_i;
  kFac = sqrt(m_i)/w_pe;

  w1 = wfac*w1;
  w2 = wfac*w2;
  w3 = wfac*w3;
  w4 = wfac*w4;
  w5 = wfac*w5;
  kc_x = kFac*kc_x;
  kc_z = kFac*kc_z;
  kc_x_max = kFac*kc_x_max;
  kc_z_max = kFac*kc_z_max;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot the surfaces we want
%%%%%%%%%%%%%%%%%%%%%%%%%%%

if surfchoice(1)
  surf(kc_x,kc_z,w1,color1)
  hold on
end
if surfchoice(2)
  surf(kc_x,kc_z,w2,color2)
  hold on
end
if surfchoice(3)
  surf(kc_x,kc_z,w3,color3)
  hold on
end
if surfchoice(4)
  surf(kc_x,kc_z,w4,color4)
  hold on
end
if surfchoice(5)
  surf(kc_x,kc_z,w5,color5)
  hold on
end

if axisnorm==2
  xlabel('k_\perp V_A/\omega_{ci}')
  ylabel('k_{||} V_A/\omega_{ci}')
  zlabel('\omega/\omega_{ci}')
else
  xlabel('k_\perp c/\omega_{ce}')
  ylabel('k_{||} c/\omega_{ce}')
  zlabel('\omega/\omega_{ce}')
end

set(gca,'xlim',[0 kc_x_max]);
set(gca,'ylim',[0 kc_z_max]);

hold off
caxis(colorlimits);
h_cbar=colorbar('hor');
colormap(irf_colormap(colorMap));

% Make sure that the frequency scale starts at zero

temp=axis;
axis([temp(1:4),0,temp(6)]);

% Print the color explanations
xlabel(h_cbar,colorstring);

% Set the view as it was

view(az,el);

% Make itposible to rotate the view

rotate3d on
