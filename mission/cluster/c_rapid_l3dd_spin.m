function c=c_rapid_l3dd_spin(varargin)
% C_RAPID_L3DD_SPIN display the electron flux as function of azimuthal angle during every spin period
%
%Input parameter (6 in total)
%1. Teve_sta='2007-10-01T07:16:00Z'
%2. C?_CP_RAP_L3DD          (Rapid 3D data from CAA)
%3. C?_CP_FGM_SPIN          (FGM spin resolution from CAA)
%4. C?_CP_RAP_EFLOW_GSE     (Rotation matrix from GSE to SC, from CAA)
%5. Cl_id=?;                 cluster ID, ?=1,2,3 or 4
%6. E_channel=1;             channel 1 is 40.7 keV; channel 2 is 68.1 keV
%
%---------------------------------------
%Example:
%c_rapid_l3dd_spin('2007-10-01T07:16:00Z', C2_CP_RAP_L3DD, C2_CP_FGM_SPIN, C2_CP_RAP_EFLOW_GSE, 2, 1);
%---------------------------------------


[ax,args,nargs] = axescheck(varargin{:});
Teve_sta=args{1};
L3DD_data=args{2};
FGM_data=args{3};
EFLOW_data=args{4};
cl_id=args{5};
E_channel=args{6};


%% Init figure
nrow=4;
ncol=4;

set(0,'DefaultAxesFontSize',8);
set(0,'DefaultLineLineWidth', 0.5);
fn=figure(61);clf;
set(gcf,'PaperUnits','centimeters')
xSize = 25; ySize = 19; coef=floor(min(800/xSize,800/ySize));
xLeft = (21-xSize)/2; yTop = (30-ySize)/2;
set(gcf,'PaperPosition',[xLeft yTop xSize ySize])
set(gcf,'Position',[10 10 xSize*coef ySize*coef])


%% Rapid pitch angle dirstribution

ic=cl_id;
dobjname='L3DD_data';
varname=irf_ssub('Electron_L_Dif_flux_3D__C?_CP_RAP_L3DD',ic);
E_lab=irf_ssub('Dimension_E__C?_CP_RAP_L3DD',ic);
AZ_lab=irf_ssub('Dimension_Az__C?_CP_RAP_L3DD',ic);
Pol_lab=irf_ssub('Dimension_Pol__C?_CP_RAP_L3DD',ic);

var=eval(['getmat(' dobjname ',''' varname ''')']);
E_lab=eval(['getmat(' dobjname ',''' E_lab ''')']);
AZ_lab=eval(['getmat(' dobjname ',''' AZ_lab ''')']);
Pol_lab=eval(['getmat(' dobjname ',''' Pol_lab ''')']);

Time_lab=var.t;
Flux_l3dd=var.data;     %--1D: time, 2D: energy, 3D: azimuzal, 4D: polar
Xaxis=AZ_lab(:,1);
Yaxis=Pol_lab(:,1);

dobjname='FGM_data';
varname=irf_ssub('B_vec_xyz_gse__C?_CP_FGM_SPIN',ic);
B_gse=eval(['getmat(' dobjname ',''' varname ''')']);
Tfgm=B_gse(:,1,1,1);

dobjname='EFLOW_data';
varname=irf_ssub('SC2GSE__C?_CP_RAP_EFLOW_GSE',ic);
EFLOW_struct=eval(['getv(' dobjname ',''' varname ''')']);
Rot_sc2gse=EFLOW_struct.data;
Rot_sc2gse=squeeze(Rot_sc2gse(1,:,:));


Teve_sta=iso2epoch(Teve_sta);
ind_flux=find(Time_lab>Teve_sta);
i_flux=ind_flux(1);


minval=zeros(1,nrow*ncol); maxval=zeros(1,nrow*ncol);

for i_subplot=1:(nrow*ncol)
  h(i_subplot)=irf_subplot(nrow,ncol,-i_subplot);

  T_eve=Time_lab(i_flux);
  Flux=Flux_l3dd(i_flux,E_channel,:,:);
  Flux=double(squeeze(Flux));
  Flux=transpose(log10(Flux));
  Flux(Flux==-Inf)=NaN;

  minval(i_subplot)=min(Flux(:));
  maxval(i_subplot)=max(Flux(:));

  specFlux=pcolor(Flux);
  set(specFlux,'EdgeColor','none');
  hold on;

  %---rotation matrix from SC to GSE
  Rot=Rot_sc2gse;
  %---from GSE to SC
  Rot=transpose(Rot);

  ind_fgm=find(Tfgm>T_eve);
  i_fgm=ind_fgm(1);
  B_vec=B_gse(i_fgm,:,:,:);
  Bgse=B_vec(2:4);

  Bsc(1)=Rot(1,1)*Bgse(1)+Rot(1,2)*Bgse(2)+Rot(1,3)*Bgse(3);
  Bsc(2)=Rot(2,1)*Bgse(1)+Rot(2,2)*Bgse(2)+Rot(2,3)*Bgse(3);
  Bsc(3)=Rot(3,1)*Bgse(1)+Rot(3,2)*Bgse(2)+Rot(3,3)*Bgse(3);

  pit_ang=zeros(length(Yaxis),length(Xaxis));
  for ix=1:length(Xaxis)
    for iy=1:length(Yaxis)
      AZ=Xaxis(ix)*(pi/180);
      Pol=Yaxis(iy)*(pi/180);
      [xfsc, yfsc, zfsc]=sph2cart(AZ-pi,Pol-(pi/2),1.0);
      v1=[xfsc yfsc zfsc];
      v2=Bsc;
      pit_ang(iy,ix)=acosd(dot(v1,v2)/(norm(v1)*norm(v2)));
    end
  end

  [C, PAcontour]=contour(pit_ang,[7 20 60 90 120 160 173]);
  clabel(C, PAcontour, [7 60 90 120 173], 'Rotation',0);
  hold off;

  set(gca,'yscale','lin');
  set(gca,'xtick', 1:2:16, 'ytick', 1:8);
  set(gca,'YDir','reverse');
  ylabel(gca,'Polar ang');
  Time=epoch2iso(T_eve); tlab=Time(12:23);
  text(4.0,8.6, [tlab ' UT']);%title(epoch2iso(T_eve));


  i_flux=i_flux+1;
end

hcol=colorbar('peer',h(2),'North', 'XDir','reverse', 'TickDir','out', 'XAxisLocation','top');
colormap(jet);
lhp=get(h(1),'Position'); rhp=get(h(ncol),'Position');
left=lhp(1); low=lhp(2)+lhp(4)+0.01;
width=rhp(1)+rhp(3)-lhp(1); height=0.04;
set(hcol,'Position',[left low width height]);

for i_subplot=1:(nrow*ncol)
  caxis(h(i_subplot),[min(minval) max(maxval)]);
end



end
%% save figure
% figname=['Figure_Rapid_l3dd_C' num2str(cl_id)];
% print(gcf, '-dpdf', [figname '.pdf']);


